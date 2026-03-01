#ifndef _LIBSVM_H
#define _LIBSVM_H

#include <immintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

struct svm_problem
{
	int l;
	float *y;
	float **x;
};

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED }; /* kernel_type */

struct svm_parameter
{
	int svm_type;
	int kernel_type;
	int degree;	/* for poly */
	float gamma;	/* for poly/rbf/sigmoid */
	float coef0;	/* for poly/sigmoid */

	/* these are for training only */
	float cache_size; /* in MB */
	float eps;	/* stopping criteria */
	float C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	int nr_weight;		/* for C_SVC */
	int *weight_label;	/* for C_SVC */
	float* weight;		/* for C_SVC */
	float nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
	float p;	/* for EPSILON_SVR */
	int shrinking;	/* use the shrinking heuristics */
	int probability; /* do probability estimates */
};

//
// svm_model
//
struct svm_model
{
	struct svm_parameter param;	/* parameter */
	int nr_class;		/* number of classes, = 2 in regression/one class svm */
	int l;			/* total #SV */
	float **SV;		/* SVs (SV[l]) */
	float **sv_coef;	/* coefficients for SVs in decision functions (sv_coef[k-1][l]) */
	float *rho;		/* constants in decision functions (rho[k*(k-1)/2]) */
	float *probA;		/* pariwise probability information */
	float *probB;
	float *prob_density_marks;	/* probability information for ONE_CLASS */
	int *sv_indices;        /* sv_indices[0,...,nSV-1] are values in [1,...,num_traning_data] to indicate SVs in the training set */

	/* for classification only */

	int *label;		/* label of each class (label[k]) */
	int *nSV;		/* number of SVs for each class (nSV[k]) */
				/* nSV[0] + nSV[1] + ... + nSV[k-1] = l */
	/* XXX */
	int free_sv;		/* 1 if svm_model is created by svm_load_model*/
				/* 0 if svm_model is created by svm_train */
};

struct svm_model *svm_train(const struct svm_problem *prob, const struct svm_parameter *param, int dim);

float svm_predict_score(const struct svm_model *model, const float *x, int dim);

void svm_free_model_content(struct svm_model *model_ptr);

#ifdef __cplusplus
}
#endif

#endif /* _LIBSVM_H */
