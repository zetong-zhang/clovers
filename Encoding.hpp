/**
 * @brief       Encoding functions for Z-curve.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     1.0.0
 * @date        2025-11-30
 * @modified    2026-02-28
 * @license     GNU GPLv3
 * @contact     ylin@tju.edu.cn | fgao@tju.edu.cn
 */
#ifndef ENCODING
#define ENCODING

#include <omp.h>
#include <cfloat>
#include <cmath>
#include <assert.h>
#include <iostream>
#include "BioStruct.hpp"
/** @deprecated
 * dimension of Z-curve params (all-in version).     */
#define DIM_A 765
/* dimension of Z-curve params. (simplified version) */
#define DIM_S 189
/* 
 * Map for converting ASCII chars into one-hot vectors
 *
 * A = [1, 0, 0, 0] G = [0, 1, 0, 0]
 * C = [0, 0, 1, 0] T = [0, 0, 0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
extern float ONE_HOT[][4];

namespace encoding {
    /**
     * @brief           Calculate 1-mer Z-curve params for a given sequence.
     * @param seq       The input sequence.
     * @param len       The length of the input sequence.
     * @param params    The output Z-curve parameters.
     */
    void     mono_trans(const char *seq, int len, float *params);
    /**
     * @brief           Calculate 2-mer Z-curve params for a given sequence.
     * @param seq       The input sequence.
     * @param len       The length of the input sequence.
     * @param params    The output Z-curve parameters.
     */
    void     di_trans  (const char *seq, int len, float *params);
    /**
     * @brief           Calculate 3-mer Z-curve params for a given sequence.
     * @param seq       The input sequence.
     * @param len       The length of the input sequence.
     * @param params    The output Z-curve parameters.
     */
    void     tri_trans (const char *seq, int len, float *params);
    /**
     * @deprecated
     * @brief           Calculate 4-mer Z-curve params for a given sequence.
     * @param seq       The input sequence.
     * @param len       The length of the input sequence.
     * @param params    The output Z-curve parameters.
     */
    void     quart_trans(const char *seq, int len, float *params);
    /**
     * @brief       Encode a given sequence into Z-curve params.
     * @param seq   The input sequence.
     * @param len   The length of the input sequence.
     * @param data  The output Z-curve parameters.
     */
    void     encode     (const char *seq, int len, float *data, int n_trans);
    /**
     * @brief       Encode ORFs into Z-curve params.
     * @param orfs  The input ORF array.
     * @param data  The output Z-curve params.
     */
    void     encode_orfs(bio::orf_array &orfs, float *data, int n_trans);
    /**
     * @brief           Transform Z-curve params to standardized data.
     * @param data      The input Z-curve params.
     * @param n         The number of samples.
     * @param dim       The dimension of Z-curve params.
     * @param means     The means of Z-curve params.
     * @param stds      The standard deviations of Z-curve params.
     * @return          The transformed Z-curve params.
     */
    float * std_scale(float *data, int n, int dim, const float *means, const float *stds);
    /**
     * @brief           Transform Z-curve params to min-max scaled data.
     * @param data      The input Z-curve params.
     * @param n         The number of samples.
     * @param dim       The dimension of Z-curve params.
     * @param mins      The minimums of Z-curve params.
     * @param maxs      The maximums of Z-curve params.
     * @return          The transformed Z-curve params.
     */
    float * minmax_scale(float *data, int n, int dim, float *mins, float *maxs);
}

#endif