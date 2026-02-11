#include "Model.hpp"
#include<chrono>
/* start and end addresses of the embedded meta.bin file. */
extern "C" {
    extern const unsigned char _binary_meta_bin_start[];
    extern const unsigned char _binary_meta_bin_end[];
}
/* epsilon value for floating point comparisons. */
#define EPSILON 1E-12
// Minimum set size for SVM training.
static const int MIN_SET_SIZE = 5;
/* SVM parameters. */
svm_parameter param = {
    C_SVC,   /* svm_type     */
    RBF   ,  /* kernel_type  */
    3,       /* degree       */
    1/DIM_S, /* gamma        */
    0.0,     /* coef0        */
    200,     /* cache_size   */
    1e-3,    /* eps          */
    100.0,   /* C            */
    0,       /* nr_weight    */
    nullptr, /* weight_label */
    nullptr, /* weight       */
    0.0,     /* nu           */
    0.0,     /* p            */
    true,    /* shrinking    */
    false    /* probability  */
};
/* Model parameters. */
const float *MODELS[N_MODELS];
/**
 * @brief           Calculate dot product of two vectors using AVX instructions.
 * @param x         The first vector.
 * @param y         The second vector.
 * @param dim       The dimension of the vectors.
 * @return          The dot product of the two vectors.
 */
inline double dot_df_avx(const double* x, const float* y, const int dim) {
#ifdef __AVX__
    // ---- AVX version ----
    __m256d sum_vec = _mm256_setzero_pd();
    int d = 0;
    for (; d + 4 <= dim; d += 4) {
        __m256d xv = _mm256_loadu_pd(x + d);
        __m128  yf = _mm_loadu_ps(y + d);
        __m256d yv = _mm256_cvtps_pd(yf);
        __m256d prod = _mm256_mul_pd(xv, yv);
        sum_vec = _mm256_add_pd(sum_vec, prod);
    }
    double buf[4];
    _mm256_storeu_pd(buf, sum_vec);
    double sum = buf[0] + buf[1] + buf[2] + buf[3];
    for (; d < dim; ++d)
        sum += x[d] * y[d];
    return sum;
#else
    // ---- Non-AVX fallback ----
    double sum = 0.0;
    for (int d = 0; d < dim; ++d)
        sum += x[d] * y[d];
    return sum;
#endif
}
/** 
 * @brief           Calculate dot product of two vectors using AVX instructions.
 * @param x         The first vector.
 * @param y         The second vector.
 * @param dim       The dimension of the vectors.
 * @return          The dot product of the two vectors.
 */
inline double dot_dd_avx(const double *x, const double *y, int dim)
{
#ifdef __AVX__
    __m256d sum = _mm256_setzero_pd();
    int i = 0;

    #ifdef __FMA__
    for (; i + 3 < dim; i += 4) {
        __m256d vx = _mm256_loadu_pd(x + i);
        __m256d vy = _mm256_loadu_pd(y + i);
        sum = _mm256_fmadd_pd(vx, vy, sum);
    }
    #else

    for (; i + 3 < dim; i += 4) {
        __m256d vx = _mm256_loadu_pd(x + i);
        __m256d vy = _mm256_loadu_pd(y + i);
        sum = _mm256_add_pd(sum, _mm256_mul_pd(vx, vy));
    }
    #endif

    __m128d lo = _mm256_castpd256_pd128(sum);
    __m128d hi = _mm256_extractf128_pd(sum, 1);
    __m128d s2 = _mm_add_pd(lo, hi);
    double result = _mm_cvtsd_f64(_mm_add_pd(s2, _mm_unpackhi_pd(s2, s2)));
    for (; i < dim; ++i)
        result += x[i] * y[i];

    return result;

#else
    double result = 0.0;
    for (int i = 0; i < dim; ++i)
        result += x[i] * y[i];
    return result;

#endif
}
/** 
 * @brief           Normalize a vector to unit length.
 * @param v         The vector to be normalized.
 * @param n         The dimension of the vector.
 */
static void normalize(double *v, int n) {
    double norm = std::sqrt(dot_dd_avx(v, v, n));
    if (norm < EPSILON) return;
    for (int i = 0; i < n; ++i) v[i] /= norm;
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

double *model::build_rbs_data(
    pch_array &flankings,
    svm_model* cds_model,
    double *params,
    bio::region *rbs_region,
    int plot_range,
    flt_array &exp_lens
) {
    static int dim = 7;
    int rbs_len = rbs_region ? (rbs_region->end-rbs_region->start) : 0;
    int n_alters = (int) flankings.size();
    double *rbs_data = NEW double[n_alters*dim]();
    int start_idx = plot_range / 2;
    if (rbs_data == nullptr) return nullptr;
# ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
# endif
    for (int i = 0; i < n_alters; i ++) {
        if (cds_model) {
            double *p_params = params + i*DIM_S;
            double *q_params = params + (i + n_alters)*DIM_S;
            rbs_data[i*dim+0] = svm_predict_score(cds_model, p_params, DIM_S);
            rbs_data[i*dim+1] = svm_predict_score(cds_model, q_params, DIM_S);
        }
        double *curve = NEW double[plot_range*3]();
        if (curve != nullptr) {
            encoding::z_curve(flankings[i], plot_range, curve);
            rbs_data[i*dim+2] = encoding::get_slope(curve+0*plot_range+rbs_region->start, rbs_len);
            rbs_data[i*dim+3] = encoding::get_slope(curve+1*plot_range+rbs_region->start, rbs_len);
            rbs_data[i*dim+4] = encoding::get_slope(curve+2*plot_range+rbs_region->start, rbs_len);
            delete[] curve;
        }
        char start_codon = flankings[i][start_idx];
        if (start_codon == 'A') {
            rbs_data[i*dim+5] = 0.6;
        } else if (start_codon == 'G') {
            rbs_data[i*dim+5] = 0.3;
        } else {
            rbs_data[i*dim+5] = 0.1;
        }
        rbs_data[i*dim+6] = exp_lens[i];
    }
    return rbs_data;
}

void model::mlp_predict(int index, double *data, int size, double *probas) {
    static double tmp = 2.0;
    if (!size) return;
    /* scale data */
    const float *means = MODELS[index], *stds = means + DIM_S;
    double *scaled = encoding::std_scale(data, size, DIM_S, means, stds);
    if (!scaled) return;
    /* mlp process */
    const float *model = stds + DIM_S; // first hidden w
    const float *hid_bs = model + N_HIDDEN*DIM_S;  // first hidden b
    const float *out_ws = hid_bs + N_HIDDEN;  // output w
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < size; i ++) {
        double *x = scaled + i*DIM_S;
        // Hidden Layer (765, ) -> (200, )
        double hid_out[N_HIDDEN] = {0.0};
        for (int j = 0; j < N_HIDDEN; j ++) {
            const float *hid_w = model + j*DIM_S;
            // hid_out = hid_w * x + hid_b
            hid_out[j] = dot_df_avx(x, hid_w, DIM_S);
            // ReLU activation function
            hid_out[j] = std::max(0.0, hid_out[j]+(double)hid_bs[j]);
        }
        // Output Layer (200, ) -> (1, )
        double out_b = (double)model[N_PARAMS-1];
        // out = out_w * hid_out + out_b
        probas[i] = (dot_df_avx(hid_out, out_ws, N_HIDDEN) + out_b) / tmp;
        // Sigmoid activation function
        probas[i] = 1.0/(1.0 + std::exp(-(probas[i])));
    }
    delete[] scaled;
}


bool model::fisher_train(double **data, int size, int dim, int *labels, double *coef) {
    const double lambda = 1e-8;
    if (!data || !labels || !coef || size <= 1 || dim <= 0)
        return false;

    // ---------- 1. compute class means ----------
    double *mu_pos = new double[dim]();
    double *mu_neg = new double[dim]();
    int n_pos = 0, n_neg = 0;

    for (int i = 0; i < size; ++i) {
        const double *x = data[i];
        if (labels[i] > 0) {
            for (int j = 0; j < dim; ++j)
                mu_pos[j] += x[j];
            ++n_pos;
        } else if (labels[i] < 0) {
            for (int j = 0; j < dim; ++j)
                mu_neg[j] += x[j];
            ++n_neg;
        }
    }

    if (n_pos == 0 || n_neg == 0)
        return false;

    for (int j = 0; j < dim; ++j) {
        mu_pos[j] /= n_pos;
        mu_neg[j] /= n_neg;
    }

    // ---------- 2. compute Sw ----------
    double *Sw = new double[dim * dim]();
    double *diff = new double[dim];

    for (int i = 0; i < size; ++i) {
        const double *x = data[i];
        const double *mu = (labels[i] > 0) ? mu_pos : mu_neg;

        for (int j = 0; j < dim; ++j)
            diff[j] = x[j] - mu[j];

        // Sw += diff * diff^T
        for (int r = 0; r < dim; ++r) {
            double dr = diff[r];
            double *row = Sw + r * dim;
            for (int c = 0; c < dim; ++c)
                row[c] += dr * diff[c];
        }
    }

    for (int i = 0; i < dim; ++i)
        Sw[i * dim + i] += lambda;

    // ---------- 3. mean difference ----------
    double *bvec = new double[dim];
    for (int i = 0; i < dim; ++i)
        bvec[i] = mu_pos[i] - mu_neg[i];

    // ---------- 4. solve Sw * w = bvec (Gaussian elimination) ----------
    // augmented matrix
    double *A = new double[dim * (dim + 1)];
    for (int i = 0; i < dim; ++i) {
        std::memcpy(A + i * (dim + 1), Sw + i * dim, dim * sizeof(double));
        A[i * (dim + 1) + dim] = bvec[i];
    }

    // Gaussian elimination
    for (int i = 0; i < dim; ++i) {
        double pivot = A[i * (dim + 1) + i];
        if (std::fabs(pivot) < EPSILON) pivot = (pivot >= 0 ? 1 : -1) * EPSILON;

        for (int j = i; j < dim + 1; ++j)
            A[i * (dim + 1) + j] /= pivot;

        for (int k = 0; k < dim; ++k) {
            if (k == i) continue;
            double factor = A[k * (dim + 1) + i];
            for (int j = i; j < dim + 1; ++j)
                A[k * (dim + 1) + j] -= factor * A[i * (dim + 1) + j];
        }
    }

    // solution w
    double *w = new double[dim];
    for (int i = 0; i < dim; ++i)
        w[i] = A[i * (dim + 1) + dim];

    // ---------- 5. bias ----------
    // b = -0.5 * w · (mu_pos + mu_neg)
    double *mu_sum = new double[dim];
    for (int i = 0; i < dim; ++i)
        mu_sum[i] = mu_pos[i] + mu_neg[i];

    double bias = -0.5 * dot_dd_avx(w, mu_sum, dim);

    // ---------- 6. normalization ----------
    // ||(w,b)|| = 1
    double norm2 = dot_dd_avx(w, w, dim) + bias * bias;
    double norm = std::sqrt(norm2);

    if (norm <= EPSILON)
        return false;

    for (int i = 0; i < dim; ++i)
        coef[i] = w[i] / norm;
    coef[dim] = bias / norm;

    // ---------- cleanup ----------
    delete[] mu_pos;
    delete[] mu_neg;
    delete[] Sw;
    delete[] diff;
    delete[] bvec;
    delete[] A;
    delete[] w;
    delete[] mu_sum;

    return true;
}

void model::fisher_predict(double *data, int size, int dim, double *coef, double *scores) {
    if (!scores) return;

    for (int i = 0; i < size; ++i) {
        const double *x = data + i * dim;
        double score = dot_dd_avx(x, coef, dim) + coef[dim];
        scores[i] = score;
    }
}

svm_model* model::rbf_train(double *params, int size, int dim, 
                 double *i_scores, double *mins, double *maxs) {
    svm_problem prob;
    for (int j = 0; j < dim; ++j) {
        mins[j] =  std::numeric_limits<double>::infinity();
        maxs[j] = -std::numeric_limits<double>::infinity();
    }

    // split data into training set
    int n_pos = 0, n_neg = 0;
    std::vector<int> train_indices;
    prob.y = NEW double[size];
    if (prob.y == nullptr) return nullptr;
    prob.l = 0;
    for (int i = 0; i < size; ++i) {
        double s = i_scores[i];
        if (s > UP_PROBA || (s < 0.05 && s > 1E-6)) {
            double* row = params + i*dim;
            train_indices.push_back(i);

            for (int j = 0; j < dim; ++j) {
                double v = row[j];
                if (v < mins[j]) mins[j] = v;
                if (v > maxs[j]) maxs[j] = v;
            }
            if (s > UP_PROBA) { prob.y[prob.l] = 1.0; n_pos ++; } 
            else { prob.y[prob.l] = -1.0; n_neg ++; }
            ++prob.l;
        }
    }

    if (n_pos < MIN_SET_SIZE || n_neg < MIN_SET_SIZE) return nullptr; 

    // preprocess data (scale data to [0, 1])
    for (int j = 0; j < dim; ++j) {
        double intv = maxs[j] - mins[j];
        if (intv < EPSILON) {
            for (int i = 0; i < size; ++i)
                params[i*dim + j] = 0.0;
        } else {
            double inv_intv = 1.0 / intv;
            double minv = mins[j];

            for (int i = 0; i < size; ++i) {
                double* row = params + i*dim;
                row[j] = (row[j] - minv) * inv_intv;
            }
        }
    }

    // calculate gamma
    double sum = 0.0;
    double sqsum = 0.0;
    const int total = prob.l * dim;
    for (int i = 0; i < prob.l; ++i) {
        double* row = params + train_indices[i]*dim;
        for (int j = 0; j < dim; ++j) {
            double v = row[j];
            sum += v;
            sqsum += v * v;
        }
    }
    double mean = sum / total;
    double mean_sq = mean * mean;
    double var = (sqsum / total) - mean_sq;
    if (var <= 0) var = EPSILON;
    param.gamma = 0.5 / (dim * var);
    // train svm model
    prob.x = NEW double *[prob.l];
    if (prob.x == nullptr) return nullptr;
    for (int i = 0; i < prob.l; i ++)
        prob.x[i] = params + train_indices[i]*dim;
    svm_model *model = svm_train(&prob, &param, dim);

    delete[] prob.x;
    delete[] prob.y;
    return model;
}