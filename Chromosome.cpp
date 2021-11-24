#include "Chromosome.h"

Chromosome::Chromosome() {
	nMut = 0;
	Nho = 0;
	linkNeut[0] = 0.0;
	linkNeut[1] = 0.0;
}

Chromosome::~Chromosome() {
}

void Chromosome::deleteChromo(void) {
	if (!mutations.empty()) mutations.clear();
}

double Chromosome::addDelMutation(int hom, double pos, double ss, double hh) {
	double v;
	mutation mut;
	mut.homol = hom;
	mut.s = ss;
	mut.h = hh;
	mut.pos = pos;

	if (mutations.find(pos) == mutations.end()) {
		nMut++;
		mutations[pos] = mut;
		v = (hh * ss);
	}
	else v = 1.0;

	return v;
}

void Chromosome::add_neutral_mut(int hom, double pos, double ss, double hh) {
	mutation mut;
	mut.homol = hom;
	mut.s = ss;
	mut.h = hh;
	mut.pos = pos;

	if (mutations.find(pos) == mutations.end()) {
		nMut++;
		mutations[pos] = mut;
	}
	else {
		cout << "mutation already present" << endl;
	}
}

//-------------------------------
//GB_07/04/20: add neutral locus linked to deleterious mutations (Keightley & Otto 2006)
void Chromosome::InitNeutral(double s1, double s2)
{
	linkNeut[0] = s1;
	linkNeut[1] = s2;
}
