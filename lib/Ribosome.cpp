#include "Ribosome.h"


Ribosome::Ribosome() {
	position = 0;
	length = 10;
	number_of_reactions = 1;
	reaction_count = 0;
	reaction_time = 0.0;
}

int Ribosome::getNumberOfReactions() {
	return number_of_reactions;
}

void Ribosome::setNumberOfReactions(int nr) {
	 number_of_reactions = nr;	
}

int Ribosome::getReaction() {
	return reaction_count;	
}

int Ribosome::getSize() {
	return length;
}

int Ribosome::getPosition() {
	return position;
}

double Ribosome::nextReaction(double* rate) {
	double r = ((double)rand() / (double)(RAND_MAX));
	reaction_time = -1.0*log(r)/rate[reaction_count];
	return reaction_time;
}

void Ribosome::updatePosition() {
	reaction_count++;
	if(reaction_count==number_of_reactions) {
		reaction_count = 0;
		position++;
	}
}



