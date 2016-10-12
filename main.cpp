/***************************************************************
*  Copyright (C) 2016 Yang Lu <ylu465@usc.edu>
*  Computational and Molecular Biology, Department of Biological Science
*  University of Southern California, LA, CA 90089, USA
*
*  Related publication:
*  TBA
***************************************************************/
#include "utils.h"
#include "seq_model.h"
#include "kmer.h"
#include "dist_model.h"
#include "output.h"

void print_usage_and_exit()
{
	printf("CAFE:\t aCcelerated Alignment-FrEe sequence analysis\n");
	printf("Description:\t The program provides multiple alignment-free sequence distance measures, including d2, d2star, d2shepp, Manhattan (Ma), Euclidean (Eu), Chebyshev (Ch), CVtree (Hao) and Jensen-Shannon divergence (JS).\n");
	printf("Authors:\t Yang Lu and Prof. Fengzhu Sun, Computational and Molecular Biology, University of Southern California.\n");
	printf("\nusage:\n");
	printf("./cafe [options]* -D <dist> -I <fa_files> -J <jfexe_path> -K <intK>\n");
	printf("\nMain arguments\n");
	printf("\t-D <dist>\tComma-separated list of distance measurements, the options include: D2-LongSeq, D2-NGS, D2star-LongSeq, D2star-NGS, D2shepp-LongSeq, D2shepp-NGS, Eu, Ma, Ch, Hao, CHISQ, and JS. E.g. -D D2star-LongSeq,Ma,Hao. \n");
	printf("\t-I <fa_files>\tComma-separated list of sequence fasta files, e.g. -I speciesA.fa,speciesB.fa,speciesC.fa. Pairwise similarity is calculated based upon the sequences specified with this option. \n");
	printf("\t-J <jfexe_path>\tUse jellyfish to accelerate kmer counting. <jfexe_path> denotes the file path of jellyfish executable file, e.g. jellyfish-2.2.4/bin/./jellyfish \n");
	printf("\t-K <intK>\tKmer Length\n");
	printf("\nOptions\n");
	printf("\t-L <lower-count>\tOnly consider k-mer with occurrence >= <lower-count>. The default value is 0. \n");
	printf("\t-M <order>\tMarkov Order involved in D2star-LongSeq, D2star-NGS, D2shepp-LongSeq, D2shepp-NGS. There are two possible options. The first option is one single value indicating that all the sequences use the same order. The second option is comma-separated list of orders. Notice that the length of the list should match the number of fasta files. The order value could be non-negative integer but less than Kmer length or \"-1\" with the special intention to automatically infer the suitable order (not suitable for JS). The default Markov Order is 0 (i.i.d. model).\n");
	printf("\t-S <dir>\tSave/Load calculated k-mer count binary files to the folder <dir>. Each input fasta file corresponds to particular model. \n");
	printf("\t-O <path>\tOutput results to file at <path> \n");
	printf("\t-T <type>\tThe output type as the input to downstream analysis, including: plain, phylip (as hierarchical clustering), cytoscape (as network analysis) and mds (Multidimensional Scaling as 2D plotting). E.g. -T mds. The default type is plain. \n");
	printf("\t-V <dir>\tSave visualization result to the folder <dir>. \n");
	printf("\nExamples:\n");
	printf("\t./cafe -M 0 -O output_path -S model_dir -T plain -I speciesA.fa,speciesB.fa -J jellyfish-2.2.4/bin/./jellyfish -K 10 -D D2star-LongSeq,Ma\n");
	printf("\t./cafe -M 1,A -S model_dir -I speciesA.fa,speciesB.fa -J jellyfish-2.2.4/bin/./jellyfish -K 10 -D D2star-LongSeq,Ma\n");
	printf("\t./cafe -M 0 -L 2 -I speciesA.fa,speciesB.fa -J jellyfish-2.2.4/bin/./jellyfish -K 10 -D D2star-LongSeq,Ma\n");
	printf("\n");
	exit(0);
}

int main(int argc, char* argv[])
{
	clock_t startTime, endTime;
	startTime = clock();

	std::vector<std::string> vec_distStr;
	std::vector<std::string> vec_orderStr;
	std::vector<std::string> vec_fastaFiles;

	std::vector<int> vec_order;
	std::vector<dist> vec_dist;
	std::vector<std::string> vec_namelist;
	std::vector<std::string> vec_saveURLlist;

	int i_k = 0, i_lowerCnt = 0;
	OUTPUT_TYPE outputType = PLAIN;
	bool singleStrain = false, containChiSq = false, jellyfishValid = true;
	std::string str_save_modelDir = "", str_save_vizDir = "", str_outputFileURL = "", str_jellyfishExeURL = "";
	
	if (argc < 2 || !strcmp(argv[1], "-help") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--usage")) print_usage_and_exit();
	for (int i = 1; i < argc; ++i)
	{
		if (!strcmp(argv[i], "-I") || !strcmp(argv[i], "-i")) split(std::string(argv[++i]), ",", vec_fastaFiles);
		else if (!strcmp(argv[i], "-K") || !strcmp(argv[i], "-k")) i_k = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-M") || !strcmp(argv[i], "-m")) split(std::string(argv[++i]), ",", vec_orderStr); 
		else if (!strcmp(argv[i], "-L") || !strcmp(argv[i], "-l")) i_lowerCnt = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "-s")) str_save_modelDir = std::string(argv[++i]);
		else if (!strcmp(argv[i], "-O") || !strcmp(argv[i], "-o")) str_outputFileURL = std::string(argv[++i]);
		else if (!strcmp(argv[i], "-V") || !strcmp(argv[i], "-v")) str_save_vizDir = std::string(argv[++i]);
		else if (!strcmp(argv[i], "-J") || !strcmp(argv[i], "-j"))
		{
			str_jellyfishExeURL = std::string(argv[++i]);
			if (!file_exists(str_jellyfishExeURL)) { jellyfishValid = false; std::cout << "[error]: Jellyfish executable file not exist at " << str_jellyfishExeURL << std::endl; }
		}
		else if (!strcmp(argv[i], "-T") || !strcmp(argv[i], "-t"))
		{
			++i;
			if (!strcmp(toLowerCase(std::string(argv[i])).c_str(), "plain")) outputType = PLAIN;
			else if (!strcmp(toLowerCase(std::string(argv[i])).c_str(), "phylip")) outputType = PHYLIP;
			else if (!strcmp(toLowerCase(std::string(argv[i])).c_str(), "cytoscape")) outputType = CYTOSCAPE;
			else if (!strcmp(toLowerCase(std::string(argv[i])).c_str(), "mds")) outputType = MDS;
			else printf("[warning]: The output type is unrecognized! ");
		}
		else if (!strcmp(argv[i], "-D") || !strcmp(argv[i], "-d"))
		{
			split(std::string(argv[++i]), ",", vec_distStr);

			for (int j = 0; j < vec_distStr.size(); ++j)
			{
				if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2-longseq")) vec_dist.push_back(D2_LS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2-ngs")) vec_dist.push_back(D2_NGS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2star-longseq")) vec_dist.push_back(D2STAR_LS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2star-ngs")) vec_dist.push_back(D2STAR_NGS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2shepp-longseq")) vec_dist.push_back(D2SHEPP_LS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "d2shepp-ngs")) vec_dist.push_back(D2SHEPP_NGS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "ma")) vec_dist.push_back(Ma_LS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "eu")) vec_dist.push_back(Eu_LS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "ch")) vec_dist.push_back(Ch_LS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "hao")) vec_dist.push_back(HAO_LS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "js")) vec_dist.push_back(JS);
				else if (!strcmp(toLowerCase(vec_distStr[j]).c_str(), "chisq")) 
				{ 
					vec_dist.push_back(CHISQ_LS);
					containChiSq = true;
				}
				else
				{
					printf("[warning]: The distance measurement %s is unrecognized!\n ", vec_distStr[j].c_str());
					vec_distStr.erase(vec_distStr.begin() + j);
				}
			}
		}
	}

	if (!str_save_modelDir.empty() && !dir_exists(str_save_modelDir)) { std::string cmd = "mkdir " + str_save_modelDir; system(cmd.c_str()); }

	if (!str_save_vizDir.empty() && !dir_exists(str_save_vizDir))
	{
		std::string cmd = "mkdir " + str_save_vizDir; system(cmd.c_str());
		std::string cpCssCMD = "cp -r res/css " + str_save_vizDir; system(cpCssCMD.c_str());
		std::string cpLibCMD = "cp -r res/lib " + str_save_vizDir; system(cpLibCMD.c_str());
		std::string cpLogoCMD = "cp -r res/logo.jpg " + str_save_vizDir; system(cpLogoCMD.c_str());
	}

	// validity check for Kmer length
	if (i_k <= 0)
	{
		printf("[error]: Kmer length should be positive! ");
		print_usage_and_exit();
	}

	// validity check for possible markov order
	if (vec_orderStr.empty())  vec_orderStr.push_back("0");
	if (vec_orderStr.size() != vec_fastaFiles.size())
	{
		if (vec_orderStr.size() > 1)
		{
			printf("[error]: The length of the order list should match the number of fasta files! ");
			print_usage_and_exit();
		}
		else if (vec_orderStr.size() == 1)
		{
			while (vec_orderStr.size() < vec_fastaFiles.size()) vec_orderStr.push_back(vec_orderStr.at(0));
		}
	}
	for (int i = 0; i < vec_orderStr.size(); ++i)
	{
		int order = atoi(vec_orderStr[i].c_str());
		if (order >= i_k || order < -1)
		{
			printf("[error]: Markov Order should be non-negative and less than k or '-1' indicating auto inference! ");
			std::cout << "k = " << i_k << "  order = " << order << std::endl;
			print_usage_and_exit();
		}
	}

	// validity check for input fasta files
	for (int i = 0; i < vec_fastaFiles.size(); ++i)
	{
		if (!file_exists(vec_fastaFiles[i]))
		{
			std::cout << "Input file not exist! Skip: " << vec_fastaFiles[i] << std::endl;
			vec_fastaFiles.erase(vec_fastaFiles.begin() + i);
			vec_orderStr.erase(vec_orderStr.begin() + i);
		}
	}


	// validity check and pre-processing for hashs of input fasta files
	for (int i = 0; i < vec_fastaFiles.size(); ++i)
	{
		int order = atoi(vec_orderStr[i].c_str()); 
		std::vector<int> cachedOrderVec;
		if (containChiSq) cachedOrderVec.push_back(i_k + 1);
		cachedOrderVec.push_back(i_k);
		if (0 == order) cachedOrderVec.push_back(1);
		else if (order > 0) { cachedOrderVec.push_back(order + 1); cachedOrderVec.push_back(order); }
		else { for (int minOrder = i_k-1; minOrder > 0 ; --minOrder) cachedOrderVec.push_back(minOrder); }

		std::string str_seqName = getFileName(vec_fastaFiles[i]);
		std::string str_saveDir = str_save_modelDir;
		if (!str_saveDir.empty() && !endsWith(str_saveDir, "/")) str_saveDir.append("/");
		str_saveDir.append("hash_").append(str_seqName).append("_L_").append(std::to_string(i_lowerCnt)).append("_k_");

		for (int orderIdx = 0; orderIdx < cachedOrderVec.size(); ++orderIdx)
		{
			int currK = cachedOrderVec[orderIdx];
			std::string str_saveURL = str_saveDir + std::to_string(currK); if (file_exists(str_saveURL)) continue;
			
			KmerModel* kmerModel = new KmerModel(currK, true);
			if (0 == orderIdx)
			{
				if (jellyfishValid)
				{
					std::string str_jfBinURL = str_saveURL; str_jfBinURL.append(".jf");
					std::string str_jfTabTxtURL = str_saveURL; str_jfTabTxtURL.append(".cnt");

					std::string lowerCntStr = " -L " + std::to_string(i_lowerCnt);
					if (i_lowerCnt < 2) lowerCntStr = "";

					std::string cmd1 = str_jellyfishExeURL + " count -m " + std::to_string(currK) + " -s 100M -t 20" + lowerCntStr + " -o " + str_jfBinURL + " " + vec_fastaFiles[i];
					system(cmd1.c_str());  std::cout << "Execute Command: " << cmd1 << std::endl;

					std::string cmd2 = str_jellyfishExeURL + " dump -t " + str_jfBinURL + lowerCntStr + " > " + str_jfTabTxtURL;
					system(cmd2.c_str());  std::cout << "Execute Command: " << cmd2 << std::endl;

					kmerModel->saveFromJellyFish(str_jfTabTxtURL, str_saveURL);
				}
				else
				{
					kmerModel->saveFromFasta(currK, vec_fastaFiles[i], str_saveURL);
				}
			}
			else
			{
				int prevK = cachedOrderVec[orderIdx - 1];
				std::string str_prevOrderURL = str_saveDir + std::to_string(prevK);
				kmerModel->saveFromLargerK(currK, prevK, str_prevOrderURL, str_saveURL);
			}
			delete kmerModel; std::cout << "Now save model to " << str_saveURL << std::endl;
		}
	}

	for (int i = 0; i < vec_fastaFiles.size(); ++i)
	{
		std::string str_seqName = getFileName(vec_fastaFiles[i]);

		std::string str_saveDir = str_save_modelDir;
		if (!str_saveDir.empty() && !endsWith(str_saveDir, "/"))
			str_saveDir.append("/");
		str_saveDir.append("hash_").append(str_seqName).append("_L_").append(std::to_string(i_lowerCnt)).append("_k_");

		int order = atoi(vec_orderStr[i].c_str());
		if (order < 0) order = getEstMarkovOrder(i_k, str_saveDir, str_seqName); // need auto infer order

		vec_order.push_back(order);
		vec_saveURLlist.push_back(str_saveDir);
		vec_namelist.push_back(str_seqName);
	}

	endTime = clock();
	std::cout << "Time Elapsed: " << ((float)endTime - (float)startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
	startTime = clock();

	for (int distIdx = 0; distIdx < vec_dist.size(); ++distIdx)
	{
		dist currDist = vec_dist[distIdx];
		int hashK = i_k; if (CHISQ_NGS == currDist || CHISQ_LS == currDist) hashK = i_k + 1;
		
		bool b_singleStrain = true;
		if (D2_NGS == currDist || D2STAR_NGS == currDist || D2SHEPP_NGS == currDist || CHISQ_NGS == currDist ||
			HAO_NGS == currDist || Ma_NGS == currDist || Eu_NGS == currDist || Ch_NGS == currDist) b_singleStrain = false;

		smat::Matrix<double>* simMat = new smat::Matrix<double>(vec_fastaFiles.size(), vec_fastaFiles.size(), 0);
		
		KmerModel* src_kmerModel, *trgt_kmerModel;

		for (int i = 0; i < vec_fastaFiles.size() - 1; ++i)
		{
			src_kmerModel = new KmerModel(hashK, b_singleStrain);
			src_kmerModel->load(hashK, vec_saveURLlist[i] + std::to_string(hashK));

			for (int j = i + 1; j < vec_fastaFiles.size(); ++j)
			{
				trgt_kmerModel = new KmerModel(hashK, b_singleStrain);
				trgt_kmerModel->load(hashK, vec_saveURLlist[j] + std::to_string(hashK));

				double distVal = 0;

				if (D2_LS == currDist || D2_NGS == currDist) 
					distVal = DistFactory::getInstance()->getD2dist(hashK, b_singleStrain, i_lowerCnt, src_kmerModel, trgt_kmerModel);
				else if (D2STAR_LS == currDist || D2STAR_NGS == currDist)
					distVal = DistFactory::getInstance()->getD2stardist(hashK, b_singleStrain, i_lowerCnt, vec_order.at(i), vec_order.at(j), vec_saveURLlist[i], vec_saveURLlist[j], src_kmerModel, trgt_kmerModel);
				else if (D2SHEPP_LS == currDist || D2SHEPP_NGS == currDist)
					distVal = DistFactory::getInstance()->getD2sheppdist(hashK, b_singleStrain, i_lowerCnt, vec_order.at(i), vec_order.at(j), vec_saveURLlist[i], vec_saveURLlist[j], src_kmerModel, trgt_kmerModel);
				else if (HAO_LS == currDist || HAO_NGS == currDist)
					distVal = DistFactory::getInstance()->getHaodist(hashK, b_singleStrain, i_lowerCnt, vec_saveURLlist[i], vec_saveURLlist[j], src_kmerModel, trgt_kmerModel);
				else if (Ma_LS == currDist || Ma_NGS == currDist)
					distVal = DistFactory::getInstance()->getL1dist(hashK, b_singleStrain, src_kmerModel, trgt_kmerModel);
				else if (Eu_LS == currDist || Eu_NGS == currDist)
					distVal = DistFactory::getInstance()->getL2dist(hashK, b_singleStrain, src_kmerModel, trgt_kmerModel);
				else if (Ch_LS == currDist || Ch_NGS == currDist)
					distVal = DistFactory::getInstance()->getLInfdist(hashK, b_singleStrain, src_kmerModel, trgt_kmerModel);
				else if (CHISQ_NGS == currDist || CHISQ_NGS == currDist)
					distVal = DistFactory::getInstance()->getChiSqdist(hashK, b_singleStrain, src_kmerModel, trgt_kmerModel);
				else if (JS == currDist)
				{
					if (vec_order.at(i) != vec_order.at(j)) throw std::runtime_error("JS expect same order! ");
					distVal = DistFactory::getInstance()->getJensenShannondist(hashK, b_singleStrain, vec_order.at(i), vec_saveURLlist[i], vec_saveURLlist[j], src_kmerModel, trgt_kmerModel);
				}
				simMat->set(i, j, distVal); simMat->set(j, i, distVal);
				delete trgt_kmerModel;
			}

			delete src_kmerModel;
		}

		std::cout << "-------------------------------------------------" << std::endl;
		std::cout << "Dist: \t" << vec_distStr[distIdx] << std::endl;

		if (!str_outputFileURL.empty())
		{
			std::string postfixedOutputURL = str_outputFileURL;
			postfixedOutputURL.append(".").append(vec_distStr[distIdx]);
			OutputWriter::getInstance()->writeToFile(outputType, simMat, &vec_namelist, postfixedOutputURL);
		}
		else
			OutputWriter::getInstance()->writeToConsole(outputType, simMat, &vec_namelist);

		///start visualization
		std::string str_saveVizTemp = str_save_vizDir;
		if (!str_saveVizTemp.empty() && !endsWith(str_saveVizTemp, "/")) str_saveVizTemp.append("/");
		str_saveVizTemp.append(vec_distStr[distIdx]).append("_k_").append(std::to_string(i_k)).append("_main.html_part2");
		std::ofstream html2Out(str_saveVizTemp.c_str());

		std::string str_outputVizURL = str_save_vizDir;
		if (!str_outputVizURL.empty() && !endsWith(str_outputVizURL, "/")) str_outputVizURL.append("/");
		str_outputVizURL.append(vec_distStr[distIdx]).append("_k_").append(std::to_string(i_k)).append("_main.html");
		
		html2Out << "var genomeNameArr = ["; html2Out << "\"" << vec_namelist.at(0) << "\"";
		for (int col = 1; col < vec_namelist.size(); ++col) html2Out << ",\"" << vec_namelist.at(col) << "\"";
		html2Out << "];" << std::endl;
	
		html2Out << "var genomeMat = [";
		for (int row = 0; row < vec_namelist.size(); ++row)
		{
			html2Out << "[" << simMat->get(row, 0);
			for (int col = 1; col < vec_namelist.size(); ++col) html2Out << "," << simMat->get(row, col);
			html2Out << "]";

			if (row != vec_namelist.size() - 1) html2Out << "," << std::endl;
		}
		html2Out << "];" << std::endl;
		html2Out.flush(); html2Out.close();

		std::string catCMD = "cat res/main.html_part1 " + str_saveVizTemp + " res/main.html_part3 > " + str_outputVizURL;
		system(catCMD.c_str());
		///end visualization

		endTime = clock();
		std::cout << "Time Elapsed: " << ((float)endTime - (float)startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
		startTime = clock();
	}
	
	//system("pause");
	return 0;
}