/***************************************************************
*  Copyright (C) 2016 Yang Lu <ylu465@usc.edu>
*  Computational and Molecular Biology, Department of Biological Science
*  University of Southern California, LA, CA 90089, USA
*
*  Related publication:
*  TBA
***************************************************************/
#ifndef _DIST_MODEL_H
#define	_DIST_MODEL_H

#include "kmer.h"

enum dist{
	D2_LS, D2_NGS, D2STAR_LS, D2STAR_NGS, D2SHEPP_LS, D2SHEPP_NGS,
	Eu_LS, Eu_NGS, Ma_LS, Ma_NGS, Ch_LS, Ch_NGS, HAO_LS, HAO_NGS, JS, CHISQ_NGS, CHISQ_LS
};

class L1FreqStrategy : public AbsTupleDistStrategy
{
public:
	L1FreqStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { double d_tmp_diff = src_X_w - trgt_X_w; if (d_tmp_diff < 0) d_tmp_diff = (0 - d_tmp_diff); d_result += d_tmp_diff; }
	double getDist(){ return d_result; }

private:
	double d_result = 0;
};

class L2FreqStrategy : public AbsTupleDistStrategy
{
public:
	L2FreqStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { double d_tmp_diff = src_X_w - trgt_X_w; d_tmp_diff = d_tmp_diff*d_tmp_diff; d_result += d_tmp_diff; }
	double getDist(){ return sqrt(d_result); }

private:
	double d_result = 0;
};

class LInfFreqStrategy : public AbsTupleDistStrategy
{
public:
	LInfFreqStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ d_result = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { double d_tmp_diff = src_X_w - trgt_X_w; if (d_tmp_diff < 0) d_tmp_diff = (0 - d_tmp_diff); d_result = std::max(d_result, d_tmp_diff); }
	double getDist(){ return d_result; }

private:
	double d_result = 0;
};

class D2Strategy : public AbsTupleDistStrategy
{
public:
	D2Strategy(int i_arg_k, bool b_arg_singleStrain) : AbsTupleDistStrategy(i_arg_k, b_arg_singleStrain){ sum_src_X_w_trgt_X_w = 0; sum_src_X_w_sq = 0; sum_trgt_X_w_sq = 0; }
	void dealWithTuple(double src_X_w, double trgt_X_w) { sum_src_X_w_trgt_X_w += (src_X_w * trgt_X_w); sum_src_X_w_sq += (src_X_w * src_X_w); sum_trgt_X_w_sq += (trgt_X_w * trgt_X_w); }
	double getDist(){ return 1.0 - sum_src_X_w_trgt_X_w / (sqrt(sum_src_X_w_sq)*sqrt(sum_trgt_X_w_sq)); }

private:
	double sum_src_X_w_trgt_X_w = 0, sum_src_X_w_sq = 0, sum_trgt_X_w_sq = 0;
};

class D2starStrategy : public AbsQuadStrategy
{
public:
	D2starStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsQuadStrategy(i_arg_k, b_arg_singleStrain){ sum_numerator = 0; sum_src_X_w_tilde_sq_div_EX_w = 0; sum_trgt_X_w_tilde_sq_div_EX_w = 0; }
	void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w);
	double getDist(){ return 0.5*(1.0 - sum_numerator / (sqrt(sum_src_X_w_tilde_sq_div_EX_w)*sqrt(sum_trgt_X_w_tilde_sq_div_EX_w))); }

private:
	double sum_numerator = 0, sum_src_X_w_tilde_sq_div_EX_w = 0, sum_trgt_X_w_tilde_sq_div_EX_w = 0;
};

class D2sheppStrategy : public AbsQuadStrategy
{
public:
	D2sheppStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsQuadStrategy(i_arg_k, b_arg_singleStrain){ sum_numerator = 0; sum_src_X_w_tilde_sq_div_sqr_sum = 0; sum_trgt_X_w_tilde_sq_div_sqr_sum = 0; }
	void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w);
	double getDist(){ return 0.5*(1.0 - sum_numerator / (sqrt(sum_src_X_w_tilde_sq_div_sqr_sum)*sqrt(sum_trgt_X_w_tilde_sq_div_sqr_sum))); }

private:
	double sum_numerator = 0, sum_src_X_w_tilde_sq_div_sqr_sum = 0, sum_trgt_X_w_tilde_sq_div_sqr_sum = 0;
};

class HaoStrategy : public AbsQuadStrategy
{
public:
	HaoStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsQuadStrategy(i_arg_k, b_arg_singleStrain){ sum_numerator = 0; sum_sq_src_X_w_tilde_div_EX_w = 0; sum_sq_trgt_X_w_tilde_div_EX_w = 0; }
	void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w);
	double getDist(){ return 0.5*(1.0 - sum_numerator / (sqrt(sum_sq_src_X_w_tilde_div_EX_w)*sqrt(sum_sq_trgt_X_w_tilde_div_EX_w))); }

private:
	double sum_numerator = 0, sum_sq_src_X_w_tilde_div_EX_w = 0, sum_sq_trgt_X_w_tilde_div_EX_w = 0;
};

class ChiSqStrategy : public AbsQuadStrategy
{
public:
	ChiSqStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsQuadStrategy(i_arg_k, b_arg_singleStrain){ sum = 0; }
	void dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w);
	double getDist(){ return sum; }
private:
	double sum = 0;
};

class JensenShannonStrategy : public AbsMrkvStrategy
{
public:
	JensenShannonStrategy(int i_arg_k, bool b_arg_singleStrain) : AbsMrkvStrategy(i_arg_k, b_arg_singleStrain){ sum_entropy = 0; src_entropy = 0; trgt_entropy = 0; }
	void dealWithMrkv(MarkovModel* src_mrkvModel, MarkovModel* trgt_mrkvModel);
	double getDist(){ return sqrt( - sum_entropy + (src_entropy + trgt_entropy) / 2); }

private:
	double sum_entropy = 0, src_entropy = 0, trgt_entropy = 0;
};


class DistFactory
{
public:
	static DistFactory *getInstance();

	double getL1dist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getL2dist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);
	double getLInfdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getD2dist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getChiSqdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getD2stardist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, 
		int i_src_order, int i_trgt_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
		KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getD2sheppdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
		int i_src_order, int i_trgt_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
		KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getHaodist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
		std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
		KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);

	double getJensenShannondist(int i_arg_k, bool b_arg_singleStrain, 
		int i_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
		KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel);


private:
	DistFactory(){}
	static DistFactory* instance;
};




//class AbsDistModel
//{
//public:
//	AbsDistModel(int i_arg_k, int i_arg_lowerCnt) { i_k = i_arg_k; i_lowerCnt = i_arg_lowerCnt; }
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel) = 0;
//	virtual ~AbsDistModel() { };
//
//public:
//	int i_k, i_lowerCnt;
//	bool b_singleStrain = true;
//};
//
//class AbsMrkvDistModel : public AbsDistModel
//{
//public:
//	AbsMrkvDistModel(int i_arg_k, int i_arg_lowerCnt) : AbsDistModel(i_arg_k, i_arg_lowerCnt) {}
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel) = 0;
//
//public:
//	int i_src_order, i_trgt_order;
//	std::string str_src_saveURLPrefix, str_trgt_saveURLPrefix;
//};
//
//class ChiSqDistModel : public AbsDistModel
//{
//public:
//	ChiSqDistModel(int i_arg_k) : AbsDistModel(i_arg_k, 0){}
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel);
//};
//
//class ChiSqLongSeqDistModel : public ChiSqDistModel
//{
//public:
//	ChiSqLongSeqDistModel(int i_arg_k) : ChiSqDistModel(i_arg_k) { b_singleStrain = true; }
//};
//
//class ChiSqNGSDistModel : public ChiSqDistModel
//{
//public:
//	ChiSqNGSDistModel(int i_arg_k) : ChiSqDistModel(i_arg_k) { b_singleStrain = false; }
//};
//
//class D2DistModel : public AbsDistModel
//{
//public:
//	D2DistModel(int i_arg_k, int i_arg_lowerCnt) : AbsDistModel(i_arg_k, i_arg_lowerCnt){}
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel);
//};
//
//class D2LongSeqDistModel : public D2DistModel
//{
//public:
//	D2LongSeqDistModel(int i_arg_k, int i_arg_lowerCnt) : D2DistModel(i_arg_k, i_arg_lowerCnt) { b_singleStrain = true; }
//};
//
//class D2NGSDistModel : public D2DistModel
//{
//public:
//	D2NGSDistModel(int i_arg_k, int i_arg_lowerCnt) : D2DistModel(i_arg_k, i_arg_lowerCnt) { b_singleStrain = false; }
//};
//
//class D2starDistModel : public AbsMrkvDistModel
//{
//public:
//	D2starDistModel(int i_arg_k, int i_arg_lowerCnt): AbsMrkvDistModel(i_arg_k, i_arg_lowerCnt){}
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel);
//};
//
//class D2starLongSeqDistModel : public D2starDistModel
//{
//public:
//	D2starLongSeqDistModel(int i_arg_k, int i_arg_lowerCnt): D2starDistModel(i_arg_k, i_arg_lowerCnt) { b_singleStrain = true; }
//};
//
//class D2starNGSDistModel : public D2starDistModel
//{
//public:
//	D2starNGSDistModel(int i_arg_k, int i_arg_lowerCnt): D2starDistModel(i_arg_k, i_arg_lowerCnt) { b_singleStrain = false; }
//};
//
//class D2sheppDistModel : public AbsMrkvDistModel
//{
//public:
//	D2sheppDistModel(int i_arg_k, int i_arg_lowerCnt): AbsMrkvDistModel(i_arg_k, i_arg_lowerCnt){}
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel);
//};
//
//class D2sheppLongSeqDistModel : public D2sheppDistModel
//{
//public:
//	D2sheppLongSeqDistModel(int i_arg_k, int i_arg_lowerCnt): D2sheppDistModel(i_arg_k, i_arg_lowerCnt) { b_singleStrain = true; }
//};
//
//class D2sheppNGSDistModel : public D2sheppDistModel
//{
//public:
//	D2sheppNGSDistModel(int i_arg_k, int i_arg_lowerCnt) : D2sheppDistModel(i_arg_k, i_arg_lowerCnt) { b_singleStrain = false; }
//};
//
//class HaoDistModel : public AbsMrkvDistModel
//{
//public:
//	HaoDistModel(int i_arg_k, int i_arg_lowerCnt)
//		: AbsMrkvDistModel(i_arg_k, i_arg_lowerCnt)
//	{
//		//b_singleStrain = true;
//		i_src_order = i_arg_k - 2; i_trgt_order = i_arg_k - 2;
//	}
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel);
//};
//
//class HaoLongSeqDistModel : public HaoDistModel
//{
//public:
//	HaoLongSeqDistModel(int i_arg_k, int i_arg_lowerCnt) : HaoDistModel(i_arg_k, i_arg_lowerCnt) { b_singleStrain = true; }
//};
//
//class HaoNGSDistModel : public HaoDistModel
//{
//public:
//	HaoNGSDistModel(int i_arg_k, int i_arg_lowerCnt) : HaoDistModel(i_arg_k, i_arg_lowerCnt) { b_singleStrain = false; }
//};
//
//class FreqNormDistModel : public AbsDistModel
//{
//public:
//	FreqNormDistModel(int i_arg_k, NORM arg_norm) : AbsDistModel(i_arg_k, 0){ /*b_singleStrain = true;*/ norm = arg_norm; }
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel);
////private: double dealWithL2(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel);
//public:
//	NORM norm;
//};
//
//class FreqNormLongSeqDistModel : public FreqNormDistModel
//{
//public:
//	FreqNormLongSeqDistModel(int i_arg_k, NORM arg_norm) : FreqNormDistModel(i_arg_k, arg_norm) { b_singleStrain = true; }
//};
//
//class FreqNormNGSDistModel : public FreqNormDistModel
//{
//public:
//	FreqNormNGSDistModel(int i_arg_k, NORM arg_norm) : FreqNormDistModel(i_arg_k, arg_norm) { b_singleStrain = false; }
//};
//
//
//class JensenShannonDistModel : public AbsMrkvDistModel
//{
//public:
//	JensenShannonDistModel(int i_arg_k): AbsMrkvDistModel(i_arg_k, 0){ b_singleStrain = true; }
//	virtual double dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel);
//};
//



#endif