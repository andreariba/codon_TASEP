#include <iostream>
#include <sstream>
#include <ctime>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <ctype.h>
#include <unistd.h>

#include "../lib/RiboModel.h"
#include "../lib/CodonTranslator.h"

using namespace std;


int main(int argc,char* argv[]) {

	int c;
	string name;
	string sequence;
	string codonfilename;
	string outfilename;
	double factor = 0.0;

	opterr = 1;

	int number_of_options = 0;

	while ((c = getopt (argc, argv, "n:s:c:f:")) != -1) {
		switch(c) {
			case 'n':
				name = optarg;
				outfilename = name;
				cout << "ORF name: " << name << endl;
				number_of_options++;
				break;
			case 's':
				sequence = optarg;
				//cout << "Sequence: " << sequence << endl;
				number_of_options++;
				break;
			case 'f':
				factor = atof(optarg);
				cout << "Factor: " << factor << endl;
				number_of_options++;
				break;
			case 'c':
				codonfilename = optarg;
				cout << "Codon speed file: " << codonfilename << endl;
				number_of_options++;
				break;
			case '?':
				if( optopt=='c' || optopt =='n' || optopt =='s' || optopt =='f' )
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if( isprint(optopt) )
					fprintf(stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
				return 1;
			default:
				abort();
		}
	}

	if(number_of_options!=4) {
		cout << "ERROR! Missing parameters.\nUsage: simulator.exe -c <codonspeedfile> -f <elongationspeed> -n <ORFname> -s <ORFsequence>" << endl;
		exit(-1);
	}

	//Read codon speed from file and convert the transcript
	cout << "Converting sequence to rate: " << flush;
	CodonTranslator translator;
	vector<string> codonstring;
	vector<double*> codonspeed;
	ifstream inputcodon(codonfilename.c_str());
	string line;
	int l = 0;

	//Number of reactions inferred from the file
	int nr = 0;
	while(getline(inputcodon, line)) {
		int c=0;
		string var;
		string codon = "NULL";
		vector<double> speed;
		istringstream ss(line);
		//cout << line << endl;
		while( ss >> var ) {
			if(c==0) {
				codon = var;
				//cout << c << ":" << var;
			} else {
				speed.push_back( factor*atof(var.c_str()) );
				//cout << "\t" << c << ":" << var;
			}
			c++;
		}	
		//cout << endl;
		codonstring.push_back(codon);
		cout << l << ": " << codon;
		codonspeed.push_back( new double[speed.size()] );
		for( int j=0; j<speed.size(); j++) {
			codonspeed[l][j] = speed[j];
			//cout << "\t" << codonspeed[l][j];
		}
		cout << endl;
		nr = speed.size();
		l++;		
	}
	if( l==0 ) { cout << "ERROR with codon speed file!!!" << endl; exit(-1);}
	translator.setRatePerCodon(codonstring, codonspeed, nr);
	translator.inputSequence(sequence);
	cout << "OK" << endl;

	for( int j=0; j<codonspeed.size(); j++) {
		delete [] codonspeed[j];
	}

	//Simulation settings 
	srand(time(NULL));
	vector<double> times;
	
	double initialtime = 1800.0;
	double finaltime = 200*60.0*30.0+initialtime;

	//Output file
	outfilename.append(".dat");
	ofstream outputfile(outfilename.c_str());

	for(double initiationrate=0.01; initiationrate<=factor; initiationrate*=2 ) {

		cout << "(*) " << initiationrate << " ... " << flush;

		RiboModel* TASEPsimulator=new RiboModel();
		//cout << "here" << endl;

		TASEPsimulator->setTranscript( translator.getTranscript(), nr );

		TASEPsimulator->setInitiationRate(initiationrate);

		TASEPsimulator->simulate(initialtime, finaltime);

		cout << "done\nCalculations ... " << flush;

		vector<double> proteinsynthesis = TASEPsimulator->averageProteinsSynthesisRate(100);
		vector<double> densities = TASEPsimulator->averageDensities(100);
		vector<double> queues = TASEPsimulator->averageQueuing(100);
		vector<double> profile = TASEPsimulator->getProfile();

		double mean_protein_synthesis = accumulate(proteinsynthesis.begin(), proteinsynthesis.end(), 0.0)/proteinsynthesis.size();
		double mean_densities = accumulate(densities.begin(), densities.end(), 0.0)/densities.size();
		double mean_queues = accumulate(queues.begin(), queues.end(), 0.0)/queues.size();

		cout << "done\nSynthesized proteins: " << TASEPsimulator->getProteinSynthesis().back() << "\nSynthesized proteins per s: " << TASEPsimulator->getProteinSynthesis().back()/(finaltime-initialtime) << "\nQueuing events: " << TASEPsimulator->getQueuingEvents().back() << endl;

		outputfile << factor << "\t" << initiationrate << "\t" << mean_protein_synthesis << "\t" << mean_densities << "\t" << mean_queues << "\t";
		outputfile << "[";
		for( int i=0; i<profile.size();i++) {
			outputfile << profile[i] ;
			if(i<(profile.size()-1)) outputfile << ",";
		}
		outputfile << "]" << endl;

		cout << factor << "\t" << initiationrate << "\t" << mean_protein_synthesis << "\t" << mean_densities << "\t" << mean_queues << endl;

		delete TASEPsimulator;

	}
	outputfile.close();

	return 0;
}
