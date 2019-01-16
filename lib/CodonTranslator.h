#ifndef CODONTRANSLATOR
#define CODONTRANSLATOR

#include <vector>
using std::vector;
#include <string>
using std::string;

class CodonTranslator {
public:
	~CodonTranslator();
	void setRatePerCodon(vector<string> c, vector<double*>& r, int nr);
	void updateRatePerCodon( string c, vector<double> r );
	void inputSequence(string seq);
	vector<double*>& getTranscript();
private:
	string sequence;
	vector<double*> transcript;
	vector<string> codon;
	vector<double*> rate;
	int number_of_reactions;
};

#endif
