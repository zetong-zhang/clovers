/**
 * @brief       Model functions for Z-curve.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     1.0.3
 * @date        2025-11-30
 * @modified    2026-03-18
 * @license     GNU GPLv3
 * @contact     ylin@tju.edu.cn | fgao@tju.edu.cn
 */
#ifndef ZCURVE_MODEL
#define ZCURVE_MODEL

// The number of mlp models in meta.bin
#define N_MODELS  61
// The number of neural neurons
#define N_HIDDEN  100
// The number of total params in a mlp model
#define N_PARAMS  19479
// The upstream sampling region of TIS
#define U_REGION  50
// The downstream sampling region of TIS
#define D_REGION  15
// The upper limit for threshold of seed ORFs
#define UP_PROBA  0.6

#include <iostream>
#include <fstream>
#include <random>
#include <cstring>
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include "BioUtil.hpp"
#include "svm.h"

/* number of TIS model params */
const int TIS_S = ((U_REGION+D_REGION)*64+16)*3;
/* mininum set size for Markov training */
const int MIN_MARKOV_SET = 3;
/* SVM hyper-params. */
extern svm_parameter param;
/* pretrained TIS params */
extern float tis_params[TIS_S];

namespace model {
    /**
     * @brief               Initialize the mlp models from embedded source.
     * @return              True if the models are successfully initialized.
     */
    bool       init_models() noexcept;
    /**
     * @brief               Predict the output of a mlp model.
     * @param model_id      The index of the mlp model.
     * @param data          The input data.
     * @param size          The size of the input data.
     * @param probas        The output probabilities.
     */
    void       mlp_predict(int index, float *data, int size, float *probas) noexcept;
    /**
     * @brief               Train a SVM model for CDSs.
     * @param params        The SVM parameters.
     * @param size          The size of the input data.
     * @param dim           The dimension of the input data.
     * @param offset        The offset of the input data.
     * @param init_score    The initial scores.
     * @return              The model
     */
    svm_model* rbf_train(float *params, int size, int dim, float *init_score, 
                         float *mins, float *maxs) noexcept;
    /**
     * @brief               Train a Markov model for TISs.
     * @param orfs          The ORFs.
     * @param order         The order of the Markov model.
     * @param params        The Markov model parameters.
     * @param starts        The alternative start codons.
     * @param table         The codon table code.
     * @return              True if the model is successfully trained.
     */
    bool       mm_train (bio::orf_array &orfs, int order, float *params, str_array&starts,
                         int table, float &pFU, float &pFD, int &max_alter) noexcept;
    /**
     * @brief               Revise the TISs using the Markov model.
     * @param orfs          The ORFs.
     * @param order         The order of the Markov model.
     * @param params        The Markov model parameters.
     * @param minlen       The mininum length of an ORF.
     * @return              The rate of unchanged TISs.
     */
    float      mm_revise(bio::orf_array &orfs, int order, float *params,float pFU, 
                         float pFD, int max_alter, int minlen) noexcept;
}

extern float W;

#endif