/**
 * @brief       Model functions for Z-curve.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     0.1.0
 * @date        2025-11-30
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
// The upper limit for threshold of seed ORFs
#define UP_PROBA  0.6

#include <iostream>
#include <fstream>
#include <random>
#include <cstring>
#include <iomanip>
#include <algorithm>
#include "Encoding.hpp"
#include "svm.h"

extern svm_parameter param;

namespace model {
    /**
     * @brief               Initialize the mlp models from embedded source.
     * @return              True if the models are successfully initialized.
     */
    bool       init_models();
    /**
     * @brief               Build RBS data for training.
     * @param flankings     The flankings of the ORFs.
     * @param cds_model     The CDS model.
     * @param params        The Z-curve params.
     * @param rbs_region    The RBS region.
     * @param plot_range    The plot range.
     * @param exp_lens      The expected lengths.
     * @return              The RBS data.
     */
    double *   build_rbs_data(pch_array &flankings,svm_model* cds_model, double *params,
                              bio::region *rbs_region, int plot_range, flt_array &exp_lens);
    /**
     * @brief               Predict the output of a mlp model.
     * @param model_id      The index of the mlp model.
     * @param data          The input data.
     * @param size          The size of the input data.
     * @param probas        The output probabilities.
     */
    void       mlp_predict(int index, double *data, int size, double *probas);
    /**
     * @brief               Train a Fisher linear discriminant.
     * @param data          The input data.
     * @param size          The size of the input data.
     * @param dim           The dimension of the input data.
     * @param labels        The class labels of the input data.
     * @param coef          The output coefficients.
     * @return              True if the training is successful.
     */
    bool       fisher_train(double **data, int size, int dim, int *labels, double *coef);
    /**
     * @brief               Predict the output of a Fisher linear discriminant.
     * @param data          The input data.
     * @param size          The size of the input data.
     * @param dim           The dimension of the input data.
     * @param coef          The coefficients of the Fisher linear discriminant.
     * @return              The output scores.
     */
    void      fisher_predict(double *data, int size, int dim, double *coef, double *scores);
    /**
     * @brief               Train a SVM model.
     * @param params        The SVM parameters.
     * @param size          The size of the input data.
     * @param dim           The dimension of the input data.
     * @param offset        The offset of the input data.
     * @param init_score    The initial scores.
     * @return              The model
     */
    svm_model* rbf_train(double *params, int size, int dim, double *init_score, 
                         double *mins, double *maxs);
}

#endif