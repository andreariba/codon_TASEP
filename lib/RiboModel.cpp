#include "RiboModel.h"
#include <iostream>
using namespace std;

RiboModel::RiboModel() {
	profile=NULL;
}


RiboModel::~RiboModel() {
	for(int i=0;i<transcript.size();i++) {
		if(transcript[i]!=NULL) {
			delete [] transcript[i];
		}
	}
	if( profile!=NULL ) {
		delete [] profile;
		profile=NULL;
	}
}


void RiboModel::setTranscript(vector<double*> t,int nr) {
	for(int i=0;i<transcript.size();i++) {
		if(transcript[i]!=NULL) {
			delete [] transcript[i];
		}
	}
	length = t.size();
	number_of_reactions = nr;
	for(int i=0;i<length;i++) {
		transcript.push_back(new double[number_of_reactions]);
		for(int r=0;r<number_of_reactions;r++) {
			transcript[i][r] = t[i][r];
			//cout << i << " " << t[i][r];
		}
	}
	profile = new double[length];
	for(int i=0;i<length;i++) {
		profile[i] = 0.0;
	}
}

void RiboModel::setInitiationRate(double i) {
	initiationrate = i;
}

void RiboModel::simulate(double it, double ft) {
	initialtime = it;
	finaltime = ft;
	double t = 0.0;
	int number_of_ribosomes = 0;
	vector<double> reaction_times;
	vector<double> queuing_times;
	vector<int> available_ribosomes;
	//std::cout << initiationrate << std::endl;
	ribosome.clear();
	reaction_times.clear();
	queuing_times.clear();
	available_ribosomes.clear();
	synthesizedproteins = 0;
	ribosomedensity = 0.0;
	queuingevents = 0;
	vecproteins.clear();
	vecdensity.clear();
	vecqueues.clear();
	while( t<=finaltime ) {
		if( !ribosome.empty() ) {
			// add standard reactions for every ribosome
			// ribosome.front() first element
			number_of_ribosomes = ribosome.size();
			for(int i=0; i<number_of_ribosomes; i++) {
				if(i==0) {
					reaction_times.push_back(ribosome[i]->nextReaction(transcript[ribosome[i]->getPosition()]));
					available_ribosomes.push_back(i);
				} else {
					if( ( ribosome[i-1]->getPosition()-ribosome[i]->getPosition() )>ribosome[i]->getSize() || (ribosome[i]->getReaction()+1)<ribosome[i]->getNumberOfReactions() ) {
						reaction_times.push_back(ribosome[i]->nextReaction(transcript[ribosome[i]->getPosition()]));
						available_ribosomes.push_back(i);
					} else {
						if(t>=initialtime) {
							queuing_times.push_back(ribosome[i]->nextReaction(transcript[ribosome[i]->getPosition()]));
						}
					}
				}
			}
		}
		if( ribosome.empty() || (ribosome.back()->getPosition()-ribosome.back()->getSize())>0 ) {
			//add initiation reaction
			double r = ((double)rand() / (double)(RAND_MAX));
			reaction_times.push_back(-1.0/initiationrate*log(r));
			available_ribosomes.push_back(-1);
		}
		//go on
		int event = 0;
		for(int i=0;i<reaction_times.size();i++) {
			if( reaction_times[i]<reaction_times[event] ) {
				event = i;
			}
		}
		t += reaction_times[event];
		if(t>=initialtime) {
			for(int i=0;i<queuing_times.size();i++) {
				if(queuing_times[i]<reaction_times[event]) {
					queuingevents++;
					break;
				}
			}
		}
		event = available_ribosomes[event];
		if(event>=0) {
			ribosome[event]->updatePosition();
			if(ribosome[event]->getPosition()==length) {
				delete ribosome[0];
				ribosome.erase(ribosome.begin());
				if(t>=initialtime) synthesizedproteins++;
			}
		} else if(event==-1) {
			//create new ribosome with reactions
			Ribosome* newribosome = new Ribosome();
			ribosome.push_back(newribosome);
			ribosome.back()->setNumberOfReactions(number_of_reactions);
		} else {
			//you shouldn't be here
			cout << "You should not be here" << endl;
			exit(-1);
		}

		#ifdef DEBUG
		cout << "time:" << t << "\t";
		for(int i=0; i<ribosome.size(); i++) {
			cout << i << ":" << ribosome[i]->getPosition() << "\t";
		}
		cout << endl;
		#endif	

		//Store results
		if(t>initialtime) {
			double tprime = t-initialtime;
			int l = 0;
			for(int i=ribosome.size()-1; i>=0; i--) {
				int temp = ribosome[i]->getPosition();
				while( l<temp ) {
					profile[l] = profile[l]*(tprime-reaction_times[event])/tprime;
					l++;
				}
				if(l==temp) {
					profile[l] = (profile[l]*(tprime-reaction_times[event])+reaction_times[event])/tprime;
					l++;
				} else {
					cout << "You should not be here" << endl;
					exit(-1);
				}
			}
			while( l<length ) {
				profile[l] = profile[l]*(tprime-reaction_times[event])/tprime;
				l++;
			}

			vectimes.push_back(t);
			vecproteins.push_back(synthesizedproteins);
			vecdensity.push_back(static_cast<double>(ribosome.size()));
			vecqueues.push_back(static_cast<double>(queuingevents)/vectimes.size());
		}

		//Clear vector to be used again
		available_ribosomes.clear();
		reaction_times.clear();
		queuing_times.clear();
	}
	for(int i=0; i<ribosome.size(); i++) {
		delete ribosome[i];
	}
}

vector<double> RiboModel::averageProteinsSynthesisRate(int w) {
	//vector<double> vector_averages;
	vector_averages.clear();
	double deltat = static_cast<int>((finaltime-initialtime)/w);
	int i=0;
	int j=0;
	for(double t=initialtime;t<finaltime;t+=deltat) {
		while( vectimes[i]<t ) {
			i++;
		}
		j=i;
		while( vectimes[i]<(t+deltat) ) {
			i++;
		}
		//std::cout << t << "\t" << j << "\t" << i  << "\t" << vecproteins[i] << "-" << vecproteins[j] << "/" << vectimes[i]-vectimes[j] << std::endl;
		vector_averages.push_back(static_cast<double>( (vecproteins[i]-vecproteins[j])/(vectimes[i]-vectimes[j]) ) );
	}
	return vector_averages;
}

vector<double> RiboModel::averageDensities(int w) {
	//vector<double> vector_averages;
	vector_averages.clear();
	double deltat = static_cast<int>((finaltime-initialtime)/w);
	int i=0;
	int j=0;	
	for(double t=initialtime;t<finaltime;t+=deltat) {
		double sum=0.0;
		while( vectimes[i]<t ) {
			i++;
		}
		j=i;
		while( vectimes[i]<(t+deltat) ) {
			i++;
		}
		for(int k=j;k<i;k++) {
			sum += static_cast<double>(vecdensity[k]*(vectimes[k+1]-vectimes[k]));
		}
		//std::cout << t << "\t" << j << "\t" << i  << "\t" << sum/length/(vectimes[i]-vectimes[j]) << std::endl;
		vector_averages.push_back( sum/length/(vectimes[i]-vectimes[j]) );
	}
	return vector_averages;
}

vector<double> RiboModel::averageQueuing(int w) {
	//vector<double> vector_averages;
	vector_averages.clear();
	double deltat = static_cast<int>((finaltime-initialtime)/w);
	int i=0;
	int j=0;
	for(double t=initialtime;t<finaltime;t+=deltat) {
		double sum=0.0;
		while( vectimes[i]<t ) {
			i++;
		}
		j=i;
		while( vectimes[i]<(t+deltat) ) {
			i++;
		}
		for(int k=j;k<i;k++) {
			sum += static_cast<double>(vecqueues[k]*(vectimes[k+1]-vectimes[k]));
		}
		vector_averages.push_back(sum/(vectimes[i]-vectimes[j]));
	}
	return vector_averages;
}

vector<double> RiboModel::getProfile() {
	vector<double> vecprofile;
	vecprofile.clear();
	for(int i=0;i<length;i++) {
		vecprofile.push_back(profile[i]);
	}
	return vecprofile;
}

void RiboModel::saveResult(string f) {
	ofstream outputfile(f.c_str());
	for(int i=0;i<vectimes.size();i++) {
		outputfile << vectimes[i] << "\t" << vecproteins[i] << "\t" << vecdensity[i] << "\t" << vecqueues[i] << "\n";
	}
	outputfile.close();
}

