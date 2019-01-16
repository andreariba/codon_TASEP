#ifndef RIBOSOME_H
#define RIBOSOME_H

#include <cstdlib>
#include <cmath>

class Ribosome {
public:
	Ribosome();
	void setNumberOfReactions(int nr);
	int getNumberOfReactions();
	int getReaction();
	int getSize();
	int getPosition();
	double nextReaction(double* rate);
	void updatePosition();

private:
	int position;
	int length;
	int number_of_reactions;
	int reaction_count;
	double reaction_time;
};


#endif
