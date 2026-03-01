#include "Model.hpp"
#include<chrono>
/* start and end addresses of the embedded meta.bin file. */
extern "C" {
    extern const unsigned char _binary_meta_bin_start[];
    extern const unsigned char _binary_meta_bin_end[];
}
// Minimum set size for SVM training.
static const int MIN_SET_SIZE = 5;

svm_parameter param = {
    C_SVC,    /* svm_type     */
    RBF,      /* kernel_type  */
    3,        /* degree       */
    1.0/DIM_S,/* gamma        */
    0.0,      /* coef0        */
    200,      /* cache_size   */
    1e-3,     /* eps          */
    100.0,    /* C            */
    0,        /* nr_weight    */
    nullptr,  /* weight_label */
    nullptr,  /* weight       */
    0.0,      /* nu           */
    0.0,      /* p            */
    true,     /* shrinking    */
    false     /* probability  */
};

float W = 1.0F;
/* Model parameters. */
const float *MODELS[N_MODELS];
/** 
 * @brief           Calculate dot product of two vectors using AVX instructions.
 * @param x         The first vector.
 * @param y         The second vector.
 * @param dim       The dimension of the vectors.
 * @return          The dot product of the two vectors.
 */
inline float dot_ff_avx(const float *x, const float *y, int dim) {
#ifdef __AVX__
    __m256 sum = _mm256_setzero_ps();
    int i = 0;

#ifdef __FMA__
    for (; i + 7 < dim; i += 8) {
        __m256 vx = _mm256_loadu_ps(x + i);
        __m256 vy = _mm256_loadu_ps(y + i);
        sum = _mm256_fmadd_ps(vx, vy, sum);
    }
#else
    for (; i + 7 < dim; i += 8) {
        __m256 vx = _mm256_loadu_ps(x + i);
        __m256 vy = _mm256_loadu_ps(y + i);
        sum = _mm256_add_ps(sum, _mm256_mul_ps(vx, vy));
    }
#endif
    __m128 lo = _mm256_castps256_ps128(sum);
    __m128 hi = _mm256_extractf128_ps(sum, 1);
    __m128 s2 = _mm_add_ps(lo, hi);
    s2 = _mm_add_ps(s2, _mm_movehl_ps(s2, s2));
    s2 = _mm_add_ss(s2, _mm_shuffle_ps(s2, s2, 1));
    float result = _mm_cvtss_f32(s2);
    for (; i < dim; ++i)
        result += x[i] * y[i];

    return result;
#else
    float result = 0.0f;
    for (int i = 0; i < dim; ++i)
        result += x[i] * y[i];
    return result;
#endif
}

bool model::init_models() {
    const size_t expected_size = N_MODELS * N_PARAMS * sizeof(float);
    const size_t actual_size = _binary_meta_bin_end - _binary_meta_bin_start;

    if (actual_size != expected_size || 
        reinterpret_cast<uintptr_t>(_binary_meta_bin_start) %
        alignof(float) != 0) 
    {
        std::cerr << "Error: failed to load embedded model data\n";
        return false;
    }

    const float* base = reinterpret_cast<const float*>(_binary_meta_bin_start);

    for (int i = 0; i < N_MODELS; ++i) MODELS[i] = base + i * N_PARAMS;

    return true;
}

void model::mlp_predict(int index, float *data, int size, float *probas) {
    static float tmp = 2.0;
    if (!size) return;
    /* scale data */
    const float *means = MODELS[index], *stds = means + DIM_S;
    float *scaled = encoding::std_scale(data, size, DIM_S, means, stds);
    if (!scaled) return;
    /* mlp process */
    const float *model = stds + DIM_S; // first hidden w
    const float *hid_bs = model + N_HIDDEN*DIM_S;  // first hidden b
    const float *out_ws = hid_bs + N_HIDDEN;  // output w
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < size; i ++) {
        float *x = scaled + i*DIM_S;
        // Hidden Layer (765, ) -> (200, )
        float hid_out[N_HIDDEN] = {0.0};
        for (int j = 0; j < N_HIDDEN; j ++) {
            const float *hid_w = model + j*DIM_S;
            // hid_out = hid_w * x + hid_b
            hid_out[j] = dot_ff_avx(x, hid_w, DIM_S);
            // ReLU activation function
            hid_out[j] = std::max(0.0F, hid_out[j]+(float)hid_bs[j]);
        }
        // Output Layer (200, ) -> (1, )
        float out_b = (float)model[N_PARAMS-1];
        // out = out_w * hid_out + out_b
        probas[i] = (dot_ff_avx(hid_out, out_ws, N_HIDDEN) + out_b) / tmp;
        // Sigmoid activation function
        probas[i] = 1.0/(1.0 + std::exp(-(probas[i])));
    }
    delete[] scaled;
}

svm_model* model::rbf_train(float *params, int size, int dim, 
                  float *i_scores, float *mins, float *maxs) {
    svm_problem prob;
    for (int j = 0; j < dim; ++j) {
        mins[j] =  std::numeric_limits<float>::infinity();
        maxs[j] = -std::numeric_limits<float>::infinity();
    }

    // split data into training set
    int n_pos = 0, n_neg = 0;
    std::vector<int> train_indices;
    prob.y = NEW float[size];
    if (prob.y == nullptr) return nullptr;
    prob.l = 0;
    for (int i = 0; i < size; ++i) {
        float s = i_scores[i];
        if (s > UP_PROBA || (s < 0.05 && s > 1E-6)) {
            float* row = params + i*dim;
            train_indices.push_back(i);

            for (int j = 0; j < dim; ++j) {
                float v = row[j];
                if (v < mins[j]) mins[j] = v;
                if (v > maxs[j]) maxs[j] = v;
            }
            if (s > UP_PROBA) { prob.y[prob.l] = 1.0; n_pos ++; } 
            else { prob.y[prob.l] = -1.0; n_neg ++; }
            ++prob.l;
        }
    }
    if ((float)n_neg/n_pos>=10) W = 1.25-0.108*std::log((float)n_neg/n_pos);

    if (n_pos < MIN_SET_SIZE || n_neg < MIN_SET_SIZE || n_neg / n_pos > 75) {
        return nullptr; 
    }

    // preprocess data (scale data to [0, 1])
    for (int j = 0; j < dim; ++j) {
        float intv = maxs[j] - mins[j];
        if (intv < EPSILON) {
            for (int i = 0; i < size; ++i)
                params[i*dim + j] = 0.0;
        } else {
            float inv_intv = 1.0 / intv;
            float minv = mins[j];

            for (int i = 0; i < size; ++i) {
                float* row = params + i*dim;
                row[j] = (row[j] - minv) * inv_intv;
            }
        }
    }

    // calculate gamma
    float sum = 0.0;
    float sqsum = 0.0;
    const int total = prob.l * dim;
    for (int i = 0; i < prob.l; ++i) {
        float* row = params + train_indices[i]*dim;
        for (int j = 0; j < dim; ++j) {
            float v = row[j];
            sum += v;
            sqsum += v * v;
        }
    }
    float mean = sum / total;
    float mean_sq = mean * mean;
    float var = (sqsum / total) - mean_sq;
    if (var <= 0) var = EPSILON;
    param.gamma = 0.5 / (dim * var);
    // train svm model
    prob.x = NEW float *[prob.l];
    if (prob.x == nullptr) return nullptr;
    for (int i = 0; i < prob.l; i ++)
        prob.x[i] = params + train_indices[i]*dim;
    svm_model *model = svm_train(&prob, &param, dim);

    delete[] prob.x;
    delete[] prob.y;
    return model;
}
/**
 * @brief   calculate the offset of a counter in the profile matrix.
 * 
 * @param order     The order of the counter.
 * @param ps        The pointer to the sequence.
 * @return int      The offset of the counter.
 */
inline int counter_offset(int order, char *ps) {
    int ctr_off = 0;
    for (int p = 0; p <= order; p ++) {
        if (!base2off.count(ps[p])) {
            ctr_off = -1;
            break;
        }
        ctr_off += (1<<(2*(order-p))) * base2off[ps[p]];
    }
    return ctr_off;
}
/**
 * @brief   update the counter in the profile matrix.
 * 
 * @param order     The order of the counter.
 * @param ps        The pointer to the sequence.
 * @param pm        The pointer to the profile matrix.
 * @param bkg       The pointer to the background counter.
 */
static void update_counter(int order, char *ps, float *pm, float *bkg) {
    for (int j = 0; j < (65 - order); j ++) {
        int ctr_off = counter_offset(order, ps+j);
        if (ctr_off >= 0) pm[j*64+ctr_off] ++;
        if (bkg && base2off.count(ps[j])) 
            bkg[base2off[ps[j]]] ++;
    }
}
/**
 * @brief   normalize and log the profile matrix.
 * 
 * @param order     The order of the counter.
 * @param pm        The pointer to the profile matrix.
 * @param nrows     The number of rows in the profile matrix.
 */
static void norm_log(int order, float *pm, int nrows) {
    int ncols = 1 << (2*(order+1));
    for (int i = 0; i < nrows; i ++) {
        float *prow = pm + i*64;
        for (int j = 0; j < ncols; j += 4) {
            float sum = prow[j+0] + prow[j+1] + prow[j+2] + prow[j+3];
            if (sum > 0.0F) {
                prow[j+0] = std::log(prow[j+0]/sum);
                prow[j+1] = std::log(prow[j+1]/sum);
                prow[j+2] = std::log(prow[j+2]/sum);
                prow[j+3] = std::log(prow[j+3]/sum);
            }
        }
    }
}

bool model::mm_train(
    bio::orf_array &orfs, int order, float *params, str_array &starts, 
    int table, float &pFU, float &pFD, int &max_alter
) {
    float bkg[4] = { 0.0F };
    std::fill_n(params, TIS_S, 1.0F);
    for (int i = 0; i < orfs.size(); i ++) {
        // only sampling from long orfs
        if (orfs[i].len < 300) continue;
        int tsp = orfs[i].t_start, h_len = orfs[i].host_len;
        char *pstr = orfs[i].pstr;
        for (int is = 0; is < std::min((int)orfs[i].starts.size(), 3); is ++) {
            // sampling putative true and downstream false
            float *pm = params + TIS_S/3;
            int sp = orfs[i].starts[is];
            if (sp < tsp || sp+15 > h_len || sp-50 < 0) continue;
            else if (sp == tsp) pm = params;
            update_counter(order, pstr+(sp-tsp)-50, pm, nullptr);
        }
        // sampling upstream false, regardless of frame
        float *pm = params + TIS_S/3*2;
        int is = std::max(50,tsp-90), end = std::min(tsp-3, h_len-15);
        while (is < end) {
            char *ps = pstr + (is - tsp);
            if (bio_util::match_codon(ps, starts))
                update_counter(order, ps-50, pm, bkg);
            is += 1;
        }
    }
    // normalize PWM
    for (int i = 0; i < 3; i ++)
        norm_log(order, params+i*TIS_S/3, 65);
    float C = bkg[0] + bkg[1] + bkg[2] + bkg[3];
    if (C > 0.0F) for (int i = 0; i < 4; i ++) bkg[i] /= C;

    float a = bkg[0], c = bkg[1], g = bkg[2], t = bkg[3];
    float s = a*t*g;
    // support for standard
    if (table != 1) s += g*t*g+t*t*g;
    float e = t*a*g+t*a*a;
    // support for mycoplasma
    if (table != 4) e += t*g*a;
    float p = s / (s + e), qpi = e / (s + e), P = 0.0F;

    pFU = 0, max_alter = 0;
    int i = 0;
    while (true) {
        P += qpi;
        if (P > 0.99) {
			max_alter = i + 1;
			break;
		}
        pFU += i * qpi;
        qpi *= p;
        i++;
    }
    pFD = max_alter - 1 - pFU;
    return true;
}

float model::mm_revise(
    bio::orf_array &orfs, int order, float *params,
    float &pFU, float &pFD, int &max_alter
) {
    float *tT = params, *dF = params + TIS_S/3, *uF = params + TIS_S/3*2;
    int n_orfs = orfs.size(), unc = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:unc) schedule(static)
#endif
    for (int i = 0; i < n_orfs; i ++) {
        int tsp = orfs[i].t_start, n_alter = orfs[i].starts.size();
        float max_score = 0.0F;
        char *pstr = orfs[i].pstr;
        for (int is = 0; is < std::min(n_alter, max_alter); is ++) {
            int sp = orfs[i].starts[is], h_len = orfs[i].host_len;
            if (sp+15 > h_len || sp-50 < 0) continue;
            char *ps = pstr + (sp - tsp) - 50;

            float t = 0.0F, d = 0.0F, u = 0.0F;
            for(int j = 0; j < 65-order; ++j){
                int ctr_off = counter_offset(order, ps+j);
                if (ctr_off >= 0) {
                    t += tT[j*64+ctr_off];
                    d += dF[j*64+ctr_off];
                    u += uF[j*64+ctr_off];
                }
            }

            float score = 1.0F / (1 + std::exp(d-t) * pFD + std::exp(u-t) * pFU);
            if (score > max_score) {
                orfs[i].t_start = orfs[i].starts[is];
                max_score = score;
            }
        }
        orfs[i].pstr = pstr + orfs[i].t_start - tsp;
        orfs[i].seq += orfs[i].t_start - tsp;
        orfs[i].len -= orfs[i].t_start - tsp;
        if (orfs[i].pstr == pstr) unc ++;
    }
    return (float) unc / n_orfs;
}