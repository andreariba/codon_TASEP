#include "CodonTranslator.h"

#include <iostream>
using namespace std;


CodonTranslator::~CodonTranslator() {
	for(int i=0;i<rate.size();i++) {
		if(rate[i]!=NULL) {
			delete [] rate[i];
		}
	}
	for(int i=0;i<transcript.size();i++) {
		if(transcript[i]!=NULL) {
			delete [] transcript[i];
		}
	}
}

void CodonTranslator::setRatePerCodon(vector<string> c, vector<double*>& r, int nr) {
	this->~CodonTranslator();
	codon = c;
	number_of_reactions = nr;
	for(int i=0;i<r.size();i++) {
		rate.push_back(new double[number_of_reactions]);
		cout << i << ": " << codon[i] << flush;
		for(int j=0;j<number_of_reactions;j++) {
			rate[i][j] = r[i][j];	
			//cout << "\t" << rate[i][j] << flush;
		}
		cout << endl;
	}
}

void CodonTranslator::updateRatePerCodon( string co, vector<double> ra ) {
	for(int c=0;c<codon.size();c++) {
		if( codon[c].compare( co )==0 ) {
			for(int r=0;r<number_of_reactions;r++) {
				rate[c][r] = ra[r];
			}
			break;
		}
	}
	this->inputSequence(sequence);
}

void CodonTranslator::inputSequence(string seq) {
	sequence = seq;
	//cout << "Sequence: " << sequence << endl;

	for(int i=0;i<transcript.size();i++) {
		if(transcript[i]!=NULL) {
			delete [] transcript[i];
		}
	}
	transcript.clear();

	for(int i=0;i<sequence.length();i+=3) {
		for(int c=0;c<codon.size();c++) {
			if( codon[c].compare( sequence.substr(i,3) )==0 ) {	
				//cout << i << "\t" << c << ":" << codon[c] << " - " << sequence.substr(i,3) << flush;
				transcript.push_back(new double[number_of_reactions]);
				for(int r=0;r<number_of_reactions;r++) {
					//cout << "\t" << i/3 << " - " << r << "\t" << rate[c][r] << flush;
					transcript[i/3][r] = rate[c][r];
					//cout << " - " << transcript[i/3][r]<< flush;
				}
				//cout << endl;
				break;
			}
		}
	}
}

vector<double*>& CodonTranslator::getTranscript() {
	return transcript;
}

