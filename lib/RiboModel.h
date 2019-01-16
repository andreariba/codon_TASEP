#ifndef RIBOMODEL_H
#define RIBOMODEL_H

#include "Ribosome.h"
#include <fstream>
using std::ofstream;
using std::string;

#include <vector>

using std::vector;

class RiboModel {
public:
	RiboModel();
	~RiboModel();
	void setTranscript(vector<double*> t, int nr);
	void setInitiationRate(double i);
	void simulate(double initialtime, double finaltime);
	vector<double> getTimePoints() {return vectimes;}
	vector<int> getProteinSynthesis() {return vecproteins;}
	vector<double> getRibosomeDensity() {return vecdensity;}
	vector<double> getQueuingEvents() {return vecqueues;}
	vector<double> averageProteinsSynthesisRate(int w);
	vector<double> averageDensities(int w);
	vector<double> averageQueuing(int w);
	vector<double> getProfile();
	void saveResult(string f);
private:
	vector<double*> transcript;
	double* profile;
	int length;
	double initialtime;
	double finaltime;
	double initiationrate;
	int number_of_reactions;
	vector<Ribosome*> ribosome;
	vector<double> vectimes;
	vector<int> vecproteins;
	vector<double> vecdensity;
	vector<double> vecqueues;
	vector<double> vector_averages;
	int synthesizedproteins;
	double ribosomedensity;
	int queuingevents;
};

#endif
