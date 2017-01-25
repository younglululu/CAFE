#include "output.h"

OutputWriter* OutputWriter::instance = 0;

OutputWriter* OutputWriter::getInstance()
{
	if (!instance) instance = new OutputWriter();
	return instance;
}

void OutputWriter::writeToFile(OUTPUT_TYPE arg_output_type, smat::Matrix<double>* arg_distMat, std::vector<std::string>* arg_nameVec, std::string str_arg_outputFileURL)
{
	std::ofstream tmp_ofsPipe(str_arg_outputFileURL.c_str(), std::ofstream::out);
	if (PLAIN == arg_output_type || CYTOSCAPE == arg_output_type)
	{
		for (int i = 0; i<arg_nameVec->size(); ++i)
			for (int j = i + 1; j<arg_nameVec->size(); ++j)
				tmp_ofsPipe << arg_nameVec->at(i) << "\t" << arg_nameVec->at(j) << "\t" << arg_distMat->get(i, j) << std::endl;
	}
	else if (PHYLIP == arg_output_type)
	{
		tmp_ofsPipe << arg_nameVec->size() << std::endl;
		for (int i = 0; i<arg_nameVec->size(); ++i)
		{
			tmp_ofsPipe << arg_nameVec->at(i);
			for (int j = 0; j<arg_nameVec->size(); ++j)
				tmp_ofsPipe << "\t" << arg_distMat->get(i, j);
			tmp_ofsPipe << std::endl;
		}
	}
	else if (MDS == arg_output_type)
	{
		smat::Matrix<double> * mat = arg_distMat->MDS_UCF(2, 30);
		for (int i = 0; i<arg_nameVec->size(); ++i)
		{
			tmp_ofsPipe << arg_nameVec->at(i);
			for (int j = 0; j<mat->columns(); ++j)
				tmp_ofsPipe << "\t" << mat->get(i, j);
			tmp_ofsPipe << std::endl;
		}
	}

	tmp_ofsPipe.close();
}


void OutputWriter::writeToConsole(OUTPUT_TYPE arg_output_type, smat::Matrix<double>* arg_distMat, std::vector<std::string>* arg_nameVec)
{
	if (PLAIN == arg_output_type || CYTOSCAPE == arg_output_type)
	{
		for (int i = 0; i<arg_nameVec->size(); ++i)
			for (int j = i + 1; j<arg_nameVec->size(); ++j)
				std::cout << arg_nameVec->at(i) << "\t" << arg_nameVec->at(j) << "\t" << arg_distMat->get(i, j) << std::endl;
	}
	else if (PHYLIP == arg_output_type)
	{
		std::cout << arg_nameVec->size() << std::endl;
		for (int i = 0; i<arg_nameVec->size(); ++i)
		{
			std::cout << arg_nameVec->at(i);
			for (int j = 0; j<arg_nameVec->size(); ++j)
				std::cout << "\t" << arg_distMat->get(i, j);
			std::cout << std::endl;
		}
	}
	else if (MDS == arg_output_type)
	{
		smat::Matrix<double> * mat = arg_distMat->MDS_UCF(2, 30);
		for (int i = 0; i<arg_nameVec->size(); ++i)
		{
			std::cout << arg_nameVec->at(i);
			for (int j = 0; j<mat->columns(); ++j)
				std::cout << "\t" << mat->get(i, j);
			std::cout << std::endl;
		}
	}
}
