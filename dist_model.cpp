#include "dist_model.h"

DistFactory* DistFactory::instance = 0;

DistFactory* DistFactory::getInstance()
{
	if (!instance) instance = new DistFactory();
	return instance;
}

void D2starStrategy::dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w)
{
	if (0 == src_EX_w || 0 == trgt_EX_w) return;

	double src_X_w_tilde = src_X_w - src_EX_w;
	double trgt_X_w_tilde = trgt_X_w - trgt_EX_w;
	double numerator = src_X_w_tilde * trgt_X_w_tilde / sqrt(src_EX_w*trgt_EX_w);
	double src_X_w_tilde_sq_div_EX_w = src_X_w_tilde*src_X_w_tilde / src_EX_w;
	double trgt_X_w_tilde_sq_div_EX_w = trgt_X_w_tilde*trgt_X_w_tilde / trgt_EX_w;

	sum_numerator += numerator;
	sum_src_X_w_tilde_sq_div_EX_w += src_X_w_tilde_sq_div_EX_w;
	sum_trgt_X_w_tilde_sq_div_EX_w += trgt_X_w_tilde_sq_div_EX_w;
}

void D2sheppStrategy::dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w)
{
	double src_X_w_tilde = src_X_w - src_EX_w;
	double trgt_X_w_tilde = trgt_X_w - trgt_EX_w;
	double src_X_w_tilde_sq = src_X_w_tilde*src_X_w_tilde;
	double trgt_X_w_tilde_sq = trgt_X_w_tilde*trgt_X_w_tilde;
	double denominator = sqrt(src_X_w_tilde_sq + trgt_X_w_tilde_sq);

	if (0 == denominator) return;

	sum_numerator += src_X_w_tilde * trgt_X_w_tilde / denominator;
	sum_src_X_w_tilde_sq_div_sqr_sum += src_X_w_tilde_sq / denominator;
	sum_trgt_X_w_tilde_sq_div_sqr_sum += trgt_X_w_tilde_sq / denominator;
}

void HaoStrategy::dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w)
{
	if (0 == src_EX_w || 0 == trgt_EX_w) return;

	double src_X_w_tilde_div_EX_w = ((double)src_X_w - src_EX_w) / src_EX_w;
	double trgt_X_w_tilde_div_EX_w = ((double)trgt_X_w - trgt_EX_w) / trgt_EX_w;

	sum_numerator += src_X_w_tilde_div_EX_w*trgt_X_w_tilde_div_EX_w;
	sum_sq_src_X_w_tilde_div_EX_w += src_X_w_tilde_div_EX_w*src_X_w_tilde_div_EX_w;
	sum_sq_trgt_X_w_tilde_div_EX_w += trgt_X_w_tilde_div_EX_w*trgt_X_w_tilde_div_EX_w;
}

void ChiSqStrategy::dealWithQuad(double src_X_w, double src_EX_w, double trgt_X_w, double trgt_EX_w)
{
	//double tmp1 = src_X_w_1 - src_X_w*all_X_w_1 / all_X_w; sum += (tmp1*tmp1*all_X_w / (src_X_w*all_X_w_1));
	double tmp1 = src_EX_w - src_X_w*trgt_EX_w / trgt_X_w; sum += (tmp1*tmp1*trgt_X_w / (src_X_w*trgt_EX_w));
}

void JensenShannonStrategy::dealWithMrkv(MarkovModel* src_mrkvModel, MarkovModel* trgt_mrkvModel)
{
	int i_src_order = src_mrkvModel->getOrder();

	for (unsigned long long i = 0; i < (unsigned long long)pow(BASE, i_src_order); ++i)
	{
		double sum_entropyOverCol = 0, src_entropyOverCol = 0, trgt_entropyOverCol = 0;

		for (unsigned int j = 0; j < BASE; ++j)
		{
			if (0 != src_mrkvModel->getTransProb(i, j)) src_entropyOverCol += exp(src_mrkvModel->getTransProb(i, j)) * src_mrkvModel->getTransProb(i, j) / LOG2;
			if (0 != trgt_mrkvModel->getTransProb(i, j)) trgt_entropyOverCol += exp(trgt_mrkvModel->getTransProb(i, j)) * trgt_mrkvModel->getTransProb(i, j) / LOG2;
			double d_tmp_sumTransProb = (exp(src_mrkvModel->getTransProb(i, j)) + exp(trgt_mrkvModel->getTransProb(i, j))) / 2;

			if (0 != d_tmp_sumTransProb) sum_entropyOverCol += d_tmp_sumTransProb * log(d_tmp_sumTransProb) / LOG2;
		}
		src_entropy += exp(src_mrkvModel->getMargProb(i)) * src_entropyOverCol;
		trgt_entropy += exp(trgt_mrkvModel->getMargProb(i)) * trgt_entropyOverCol;
		sum_entropy += (exp(src_mrkvModel->getMargProb(i)) + exp(trgt_mrkvModel->getMargProb(i))) / 2 * sum_entropyOverCol;
	}
}


double DistFactory::getL1dist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	L1FreqStrategy* strategy = new L1FreqStrategy(i_arg_k, b_arg_singleStrain);

	double dist = IterFactory::getInstance()->getFreqDist(strategy, 
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getL2dist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	L2FreqStrategy* strategy = new L2FreqStrategy(i_arg_k, b_arg_singleStrain);

	double dist = IterFactory::getInstance()->getFreqDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getLInfdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	LInfFreqStrategy* strategy = new LInfFreqStrategy(i_arg_k, b_arg_singleStrain);

	double dist = IterFactory::getInstance()->getFreqDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, arg_trgtKmerModel->totalKmer());

	delete strategy;
	return dist;
}

double DistFactory::getD2dist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	D2Strategy* strategy = new D2Strategy(i_arg_k, b_arg_singleStrain);

	double dist = IterFactory::getInstance()->getCntDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, 
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getChiSqdist(int i_arg_k, bool b_arg_singleStrain, KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	ChiSqStrategy* strategy = new ChiSqStrategy(i_arg_k, b_arg_singleStrain);

	double dist = IterFactory::getInstance()->getCntDist(strategy,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec,
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec);

	delete strategy;
	return dist;
}

double DistFactory::getD2stardist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
	int i_src_order, int i_trgt_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
	KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	D2starStrategy* strategy = new D2starStrategy(i_arg_k, b_arg_singleStrain);

	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_src_order, str_src_saveURLPrefix);
	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_trgt_order, str_trgt_saveURLPrefix);

	double dist = IterFactory::getInstance()->getCntExpDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, src_markovModel, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, trgt_markovModel, arg_trgtKmerModel->totalKmer());

	delete strategy; delete src_markovModel; delete trgt_markovModel;
	return dist;
}

double DistFactory::getD2sheppdist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
	int i_src_order, int i_trgt_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
	KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	D2sheppStrategy* strategy = new D2sheppStrategy(i_arg_k, b_arg_singleStrain);

	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_src_order, str_src_saveURLPrefix);
	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_trgt_order, str_trgt_saveURLPrefix);

	double dist = IterFactory::getInstance()->getCntExpDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, src_markovModel, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, trgt_markovModel, arg_trgtKmerModel->totalKmer());

	delete strategy; delete src_markovModel; delete trgt_markovModel;
	return dist;
}

double DistFactory::getHaodist(int i_arg_k, bool b_arg_singleStrain, int i_arg_lowerCnt,
	std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
	KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	HaoStrategy* strategy = new HaoStrategy(i_arg_k, b_arg_singleStrain);

	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_arg_k - 2, str_src_saveURLPrefix);
	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_arg_k - 2, str_trgt_saveURLPrefix);

	double dist = IterFactory::getInstance()->getCntExpDist(strategy, i_arg_lowerCnt,
		arg_srcKmerModel->kmerCntUnorderMap, arg_srcKmerModel->kmerVec, src_markovModel, arg_srcKmerModel->totalKmer(),
		arg_trgtKmerModel->kmerCntUnorderMap, arg_trgtKmerModel->kmerVec, trgt_markovModel, arg_trgtKmerModel->totalKmer());

	delete strategy; delete src_markovModel; delete trgt_markovModel;
	return dist;
}

double DistFactory::getJensenShannondist(int i_arg_k, bool b_arg_singleStrain,
	int i_order, std::string str_src_saveURLPrefix, std::string str_trgt_saveURLPrefix,
	KmerModel* arg_srcKmerModel, KmerModel* arg_trgtKmerModel)
{
	JensenShannonStrategy* strategy = new JensenShannonStrategy(i_arg_k, b_arg_singleStrain);

	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_order, str_src_saveURLPrefix);
	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_order, str_trgt_saveURLPrefix);

	double dist = IterFactory::getInstance()->getMrkvDist(strategy, src_markovModel, trgt_markovModel);

	delete strategy; delete src_markovModel; delete trgt_markovModel;
	return dist;
}


//double ChiSqDistModel::dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel)
//{
//	AbsIterator<double>* src_kmerCntIter = arg_srcKmerModel->getKmerCntIterator();
//	AbsIterator<double>* trgt_kmerCntIter = arg_trgtKmerModel->getKmerCntIterator();
//
//	double src_X_w_1, src_X_w_2, src_X_w_3, src_X_w_4, src_X_w;
//	double trgt_X_w_1, trgt_X_w_2, trgt_X_w_3, trgt_X_w_4, trgt_X_w;
//	double all_X_w_1, all_X_w_2, all_X_w_3, all_X_w_4, all_X_w;
//
//	double sum = 0;
//	while (true)
//	{
//		if (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext()) { src_X_w_1 = *(*src_kmerCntIter); (*src_kmerCntIter)++; trgt_X_w_1 = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++; } else break;
//		if (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext()) { src_X_w_2 = *(*src_kmerCntIter); (*src_kmerCntIter)++; trgt_X_w_2 = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++; } else break;
//		if (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext()) { src_X_w_3 = *(*src_kmerCntIter); (*src_kmerCntIter)++; trgt_X_w_3 = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++; } else break;
//		if (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext()) { src_X_w_4 = *(*src_kmerCntIter); (*src_kmerCntIter)++; trgt_X_w_4 = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++; } else break;
//
//		src_X_w = src_X_w_1 + src_X_w_2 + src_X_w_3 + src_X_w_4;
//		trgt_X_w = trgt_X_w_1 + trgt_X_w_2 + trgt_X_w_3 + trgt_X_w_4;
//		all_X_w_1 = src_X_w_1 + trgt_X_w_1;
//		all_X_w_2 = src_X_w_2 + trgt_X_w_2;
//		all_X_w_3 = src_X_w_3 + trgt_X_w_3;
//		all_X_w_4 = src_X_w_4 + trgt_X_w_4;
//		all_X_w = src_X_w + trgt_X_w;
//
//		if (0 == all_X_w) continue;
//
//		if (all_X_w_1 > 0)
//		{
//			if (src_X_w > 0) { double tmp1 = src_X_w_1 - src_X_w*all_X_w_1 / all_X_w; sum += (tmp1*tmp1*all_X_w / (src_X_w*all_X_w_1)); }
//			if (trgt_X_w > 0) { double tmp2 = trgt_X_w_1 - trgt_X_w*all_X_w_1 / all_X_w; sum += (tmp2*tmp2*all_X_w / (trgt_X_w*all_X_w_1)); }
//		}
//		if (all_X_w_2 > 0)
//		{
//			if (src_X_w > 0) { double tmp1 = src_X_w_2 - src_X_w*all_X_w_2 / all_X_w; sum += (tmp1*tmp1*all_X_w / (src_X_w*all_X_w_2)); }
//			if (trgt_X_w > 0) { double tmp2 = trgt_X_w_2 - trgt_X_w*all_X_w_2 / all_X_w; sum += (tmp2*tmp2*all_X_w / (trgt_X_w*all_X_w_2)); }
//		}
//		if (all_X_w_3 > 0)
//		{
//			if (src_X_w > 0) { double tmp1 = src_X_w_3 - src_X_w*all_X_w_3 / all_X_w; sum += (tmp1*tmp1*all_X_w / (src_X_w*all_X_w_3)); }
//			if (trgt_X_w > 0) { double tmp2 = trgt_X_w_3 - trgt_X_w*all_X_w_3 / all_X_w; sum += (tmp2*tmp2*all_X_w / (trgt_X_w*all_X_w_3)); }
//		}
//		if (all_X_w_4 > 0)
//		{
//			if (src_X_w > 0) { double tmp1 = src_X_w_4 - src_X_w*all_X_w_4 / all_X_w; sum += (tmp1*tmp1*all_X_w / (src_X_w*all_X_w_4)); }
//			if (trgt_X_w > 0) { double tmp2 = trgt_X_w_4 - trgt_X_w*all_X_w_4 / all_X_w; sum += (tmp2*tmp2*all_X_w / (trgt_X_w*all_X_w_4)); }
//		}
//	}
//	delete src_kmerCntIter; delete trgt_kmerCntIter;
//	return sum;
//}
//
//double D2DistModel::dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel)
//{
//	AbsIterator<double>* src_kmerCntIter = arg_srcKmerModel->getKmerCntIterator();
//	AbsIterator<double>* trgt_kmerCntIter = arg_trgtKmerModel->getKmerCntIterator();
//
//	double sum_src_X_w_trgt_X_w = 0, sum_src_X_w_sq = 0, sum_trgt_X_w_sq = 0;
//	while (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext())
//	{
//		double src_X_w = *(*src_kmerCntIter); (*src_kmerCntIter)++;
//		double trgt_X_w = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++;
//
//		if (src_X_w < i_lowerCnt || trgt_X_w < i_lowerCnt) continue;
//
//		double tmp1 = src_X_w * trgt_X_w;
//		double tmp2 = src_X_w * src_X_w;
//		double tmp3 = trgt_X_w * trgt_X_w;
//
//		sum_src_X_w_trgt_X_w += tmp1;
//		sum_src_X_w_sq += tmp2;
//		sum_trgt_X_w_sq += tmp3;
//	}
//
//	delete src_kmerCntIter; delete trgt_kmerCntIter;
//	return 1.0 - sum_src_X_w_trgt_X_w / (sqrt(sum_src_X_w_sq)*sqrt(sum_trgt_X_w_sq));
//}
//
//double D2starDistModel::dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel)
//{
//	AbsIterator<double>* src_kmerCntIter = arg_srcKmerModel->getKmerCntIterator();
//	AbsIterator<double>* trgt_kmerCntIter = arg_trgtKmerModel->getKmerCntIterator();
//
//	KmerProbEnsembDelegate* src_kmerProbDelegate = arg_srcKmerModel->getKmerProDelegate(i_src_order, str_src_saveURLPrefix);
//	KmerProbEnsembDelegate* trgt_kmerProbDelegate = arg_trgtKmerModel->getKmerProDelegate(i_trgt_order, str_trgt_saveURLPrefix);
//
//	double sum_numerator = 0, sum_src_X_w_tilde_sq_div_EX_w = 0, sum_trgt_X_w_tilde_sq_div_EX_w = 0;
//
//	double log_src_totalKmerLen = log(arg_srcKmerModel->totalKmer());
//	double log_trgt_totalKmerLen = log(arg_trgtKmerModel->totalKmer());
//
//	unsigned long long kmerIdx = 0;
//	while (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext())
//	{
//		kmerIdx++;
//		double src_X_w = *(*src_kmerCntIter); (*src_kmerCntIter)++;
//		double trgt_X_w = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++;
//
//		if (src_X_w < i_lowerCnt || trgt_X_w < i_lowerCnt) continue;
//
//		//std::cout << "kmer = " << (kmerIdx - 1) << std::endl;
//		double src_prob_X_w = src_kmerProbDelegate->getKmerlogProb(kmerIdx - 1);
//		//std::cout << "src_prob_X_w = " << src_prob_X_w << std::endl;
//		double trgt_prob_X_w = trgt_kmerProbDelegate->getKmerlogProb(kmerIdx - 1);
//		//std::cout << "trgt_prob_X_w = " << trgt_prob_X_w << std::endl;
//
//		if (0 == src_prob_X_w || 0 == trgt_prob_X_w) continue;
//
//		double src_EX_w = exp(log_src_totalKmerLen + src_prob_X_w);
//		double trgt_EX_w = exp(log_trgt_totalKmerLen + trgt_prob_X_w);
//
//		double src_X_w_tilde = (double)src_X_w - src_EX_w;
//		double trgt_X_w_tilde = (double)trgt_X_w - trgt_EX_w;
//		double numerator = src_X_w_tilde * trgt_X_w_tilde / sqrt(src_EX_w*trgt_EX_w);
//		double src_X_w_tilde_sq_div_EX_w = src_X_w_tilde*src_X_w_tilde / src_EX_w;
//		double trgt_X_w_tilde_sq_div_EX_w = trgt_X_w_tilde*trgt_X_w_tilde / trgt_EX_w;
//
//		if (almostEquals(0, src_EX_w) || almostEquals(0, trgt_EX_w)) continue;
//
//		sum_numerator += numerator;
//		sum_src_X_w_tilde_sq_div_EX_w += src_X_w_tilde_sq_div_EX_w;
//		sum_trgt_X_w_tilde_sq_div_EX_w += trgt_X_w_tilde_sq_div_EX_w;
//	}
//
//	delete src_kmerCntIter; delete trgt_kmerCntIter;
//	delete src_kmerProbDelegate; delete trgt_kmerProbDelegate;
//	return 0.5*(1.0 - sum_numerator / (sqrt(sum_src_X_w_tilde_sq_div_EX_w)*sqrt(sum_trgt_X_w_tilde_sq_div_EX_w)));
//}
//
//double D2sheppDistModel::dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel)
//{
//	AbsIterator<double>* src_kmerCntIter = arg_srcKmerModel->getKmerCntIterator();
//	AbsIterator<double>* trgt_kmerCntIter = arg_trgtKmerModel->getKmerCntIterator();
//
//	KmerProbEnsembDelegate* src_kmerProbDelegate = arg_srcKmerModel->getKmerProDelegate(i_src_order, str_src_saveURLPrefix);
//	KmerProbEnsembDelegate* trgt_kmerProbDelegate = arg_trgtKmerModel->getKmerProDelegate(i_trgt_order, str_trgt_saveURLPrefix);
//
//	double sum_numerator = 0, sum_src_X_w_tilde_sq_div_sqr_sum = 0, sum_trgt_X_w_tilde_sq_div_sqr_sum = 0;
//
//	double log_src_totalKmerLen = log(arg_srcKmerModel->totalKmer());
//	double log_trgt_totalKmerLen = log(arg_trgtKmerModel->totalKmer());
//
//	unsigned long long kmerIdx = 0;
//	while (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext())
//	{
//		kmerIdx++;
//		double src_X_w = *(*src_kmerCntIter); (*src_kmerCntIter)++;
//		double trgt_X_w = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++;
//
//		if (src_X_w < i_lowerCnt || trgt_X_w < i_lowerCnt) continue;
//
//		double src_prob_X_w = src_kmerProbDelegate->getKmerlogProb(kmerIdx - 1);
//		double trgt_prob_X_w = trgt_kmerProbDelegate->getKmerlogProb(kmerIdx - 1);
//
//		double src_EX_w = 0, trgt_EX_w = 0;
//		if (0 != src_prob_X_w) src_EX_w = exp(log_src_totalKmerLen + src_prob_X_w);
//		if (0 != trgt_prob_X_w) trgt_EX_w = exp(log_trgt_totalKmerLen + trgt_prob_X_w);
//
//		double src_X_w_tilde = (double)src_X_w - src_EX_w;
//		double trgt_X_w_tilde = (double)trgt_X_w - trgt_EX_w;
//		double src_X_w_tilde_sq = src_X_w_tilde*src_X_w_tilde;
//		double trgt_X_w_tilde_sq = trgt_X_w_tilde*trgt_X_w_tilde;
//		double denominator = sqrt(src_X_w_tilde_sq + trgt_X_w_tilde_sq);
//
//		if (almostEquals(0, denominator)) continue;
//
//		sum_numerator += src_X_w_tilde * trgt_X_w_tilde / denominator;
//		sum_src_X_w_tilde_sq_div_sqr_sum += src_X_w_tilde_sq / denominator;
//		sum_trgt_X_w_tilde_sq_div_sqr_sum += trgt_X_w_tilde_sq / denominator;
//	}
//
//	delete src_kmerCntIter; delete trgt_kmerCntIter;
//	delete src_kmerProbDelegate; delete trgt_kmerProbDelegate;
//	return 0.5*(1.0 - sum_numerator / (sqrt(sum_src_X_w_tilde_sq_div_sqr_sum)*sqrt(sum_trgt_X_w_tilde_sq_div_sqr_sum)));
//}
//
//double HaoDistModel::dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel)
//{
//	AbsIterator<double>* src_kmerCntIter = arg_srcKmerModel->getKmerCntIterator();
//	AbsIterator<double>* trgt_kmerCntIter = arg_trgtKmerModel->getKmerCntIterator();
//
//	KmerProbEnsembDelegate* src_kmerProbDelegate = arg_srcKmerModel->getKmerProDelegate(i_src_order, str_src_saveURLPrefix);
//	KmerProbEnsembDelegate* trgt_kmerProbDelegate = arg_trgtKmerModel->getKmerProDelegate(i_trgt_order, str_trgt_saveURLPrefix);
//
//	double sum_numerator = 0, sum_sq_src_X_w_tilde_div_EX_w = 0, sum_sq_trgt_X_w_tilde_div_EX_w = 0;
//
//	double log_src_totalKmerLen = log(arg_srcKmerModel->totalKmer());
//	double log_trgt_totalKmerLen = log(arg_trgtKmerModel->totalKmer());
//
//	unsigned long long kmerIdx = 0;
//	while (src_kmerCntIter->hasNext() && trgt_kmerCntIter->hasNext())
//	{
//		kmerIdx++;
//		double src_X_w = *(*src_kmerCntIter); (*src_kmerCntIter)++;
//		double trgt_X_w = *(*trgt_kmerCntIter); (*trgt_kmerCntIter)++;
//
//		if (src_X_w < i_lowerCnt || trgt_X_w < i_lowerCnt) continue;
//
//		double src_prob_X_w = src_kmerProbDelegate->getKmerlogProb(kmerIdx - 1);
//		double trgt_prob_X_w = trgt_kmerProbDelegate->getKmerlogProb(kmerIdx - 1);
//
//		if (0 == src_prob_X_w || 0 == trgt_prob_X_w) continue;
//
//		double src_EX_w = exp(log_src_totalKmerLen + src_prob_X_w);
//		double trgt_EX_w = exp(log_trgt_totalKmerLen + trgt_prob_X_w);
//
//		double src_X_w_tilde_div_EX_w = ((double)src_X_w - src_EX_w) / src_EX_w;
//		double trgt_X_w_tilde_div_EX_w = ((double)trgt_X_w - trgt_EX_w) / trgt_EX_w;
//
//		if (almostEquals(0, src_EX_w) || almostEquals(0, trgt_EX_w)) continue;
//
//		sum_numerator += src_X_w_tilde_div_EX_w*trgt_X_w_tilde_div_EX_w;
//		sum_sq_src_X_w_tilde_div_EX_w += src_X_w_tilde_div_EX_w*src_X_w_tilde_div_EX_w;
//		sum_sq_trgt_X_w_tilde_div_EX_w += trgt_X_w_tilde_div_EX_w*trgt_X_w_tilde_div_EX_w;
//	}
//
//	delete src_kmerCntIter; delete trgt_kmerCntIter;
//	delete src_kmerProbDelegate; delete trgt_kmerProbDelegate;
//	return 0.5*(1.0 - sum_numerator / (sqrt(sum_sq_src_X_w_tilde_div_EX_w)*sqrt(sum_sq_trgt_X_w_tilde_div_EX_w)));
//}
//
//double FreqNormDistModel::dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel)
//{
//	//if (L2 == norm) return dealWithL2(arg_srcKmerModel, arg_trgtKmerModel);
//
//	AbsIterator<double>* src_kmerFreqIter = arg_srcKmerModel->getKmerFreqIterator();
//	AbsIterator<double>* trgt_kmerFreqIter = arg_trgtKmerModel->getKmerFreqIterator();
//
//	double d_tmp_result = 0;
//	while (src_kmerFreqIter->hasNext() && trgt_kmerFreqIter->hasNext())
//	{
//		double src_X_freq = *(*src_kmerFreqIter); (*src_kmerFreqIter)++;
//		double trgt_X_freq = *(*trgt_kmerFreqIter); (*trgt_kmerFreqIter)++;
//
//		double d_tmp_diff = src_X_freq - trgt_X_freq; if (d_tmp_diff < 0) d_tmp_diff = (0 - d_tmp_diff);
//		if (L2 == norm) d_tmp_diff = d_tmp_diff*d_tmp_diff;
//
//		if (LINF == norm) d_tmp_result = std::max(d_tmp_result, d_tmp_diff);
//		else d_tmp_result += d_tmp_diff;
//	}
//
//	if (L2 == norm) d_tmp_result = sqrt(d_tmp_result);
//	delete src_kmerFreqIter; delete trgt_kmerFreqIter;
//	return d_tmp_result;
//}
//
///*double FreqNormDistModel::dealWithL2(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel)
//{
//	AbsIterator<double>* src_kmerFreqIter = arg_srcKmerModel->getKmerFreqIterator();
//	AbsIterator<double>* trgt_kmerFreqIter = arg_trgtKmerModel->getKmerFreqIterator();
//
//	double d_tmp_result_log = 0; bool firstTime = true;
//
//	while (src_kmerFreqIter->hasNext() && trgt_kmerFreqIter->hasNext())
//	{
//		double src_X_freq = *(*src_kmerFreqIter); (*src_kmerFreqIter)++;
//		double trgt_X_freq = *(*trgt_kmerFreqIter); (*trgt_kmerFreqIter)++;
//
//		double d_tmp_diff = src_X_freq - trgt_X_freq; if (d_tmp_diff < 0) d_tmp_diff = (0 - d_tmp_diff);
//
//		if (0 == d_tmp_diff) continue;
//
//		if (firstTime && (0 == d_tmp_result_log)) { d_tmp_result_log = 2 * log(d_tmp_diff); firstTime = false; }
//		else d_tmp_result_log = log_sum(d_tmp_result_log, 2 * log(d_tmp_diff));
//	}
//	delete src_kmerFreqIter; delete trgt_kmerFreqIter;
//	return sqrt(exp(d_tmp_result_log));
//}*/
//
//double JensenShannonDistModel::dist(AbsKmerModel* arg_srcKmerModel, AbsKmerModel* arg_trgtKmerModel)
//{
//	MarkovModel* src_markovModel = arg_srcKmerModel->getMarkovModel(i_src_order, str_src_saveURLPrefix);
//	MarkovModel* trgt_markovModel = arg_trgtKmerModel->getMarkovModel(i_trgt_order, str_trgt_saveURLPrefix);
//
//	double sum_entropy = 0, src_entropy = 0, trgt_entropy = 0;
//
//	for (unsigned long long i = 0; i < (unsigned long long)pow(BASE, i_src_order); ++i)
//	{
//		double sum_entropyOverCol = 0, src_entropyOverCol = 0, trgt_entropyOverCol = 0;
//
//		for (unsigned int j = 0; j < BASE; ++j)
//		{
//			//if (0 != src_markovModel->getTransProb(i, j)) src_entropyOverCol += src_markovModel->getTransProb(i, j) * log(src_markovModel->getTransProb(i, j)) / log(2);
//			//if (0 != trgt_markovModel->getTransProb(i, j)) trgt_entropyOverCol += trgt_markovModel->getTransProb(i, j) * log(trgt_markovModel->getTransProb(i, j)) / log(2);
//			//double d_tmp_sumTransProb = (src_markovModel->getTransProb(i, j) + trgt_markovModel->getTransProb(i, j)) / 2;
//
//			if (0 != src_markovModel->getTransProb(i, j)) src_entropyOverCol += exp(src_markovModel->getTransProb(i, j)) * src_markovModel->getTransProb(i, j) / LOG2;
//			if (0 != trgt_markovModel->getTransProb(i, j)) trgt_entropyOverCol += exp(trgt_markovModel->getTransProb(i, j)) * trgt_markovModel->getTransProb(i, j) / LOG2;
//			double d_tmp_sumTransProb = (exp(src_markovModel->getTransProb(i, j)) + exp(trgt_markovModel->getTransProb(i, j))) / 2;
//
//			if (0 != d_tmp_sumTransProb) sum_entropyOverCol += d_tmp_sumTransProb * log(d_tmp_sumTransProb) / LOG2;
//		}
//
//		//src_entropy += src_markovModel->getMargProb(i) * src_entropyOverCol;
//		//trgt_entropy += trgt_markovModel->getMargProb(i) * trgt_entropyOverCol;
//		//sum_entropy += (src_markovModel->getMargProb(i) + trgt_markovModel->getMargProb(i)) / 2 * sum_entropyOverCol;
//
//		src_entropy += exp(src_markovModel->getMargProb(i)) * src_entropyOverCol;
//		trgt_entropy += exp(trgt_markovModel->getMargProb(i)) * trgt_entropyOverCol;
//		sum_entropy += (exp(src_markovModel->getMargProb(i)) + exp(trgt_markovModel->getMargProb(i))) / 2 * sum_entropyOverCol;
//	}
//
//	src_entropy *= -1;
//	trgt_entropy *= -1;
//	sum_entropy *= -1;
//
//	delete src_markovModel; delete trgt_markovModel;
//	return sqrt(sum_entropy - (src_entropy + trgt_entropy) / 2);
//}
