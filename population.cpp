#include "population.h"
#include "individuals.h"

Population::Population(int xx, int yy) {
	N = 0;
	Noffs = 0;
	DelMut = 0;
	x = xx;
	y = yy;
	meanDelMut = 0.0;
	DelMutSd = 0.0;
	W = 0.0;
	Wsd = 0.0;
	Homoz = 0.0; Hsd = 0.0;
	Vneutral = 0.0; VneutralSD = 0.0;
	Wmax = 0.0;
	Wmin = 0.0;

	if (!popMuts.empty()) popMuts.clear();
}

Population::~Population() {
}

void Population::initialise_pop(double k, double kk, Parameters para, std::gamma_distribution<> finds, std::uniform_real_distribution<> findpos, std::normal_distribution<> neut) {

	N = (int)k;

	std::uniform_real_distribution<> unireal(0.0, 1.0);

	//initialise individuals
	for (int i = 0; i < (int)(N); i++) {
		inds.push_back(Individuals(para));
		inds[i].initialise(kk, para, finds, findpos, neut);
		W += 1.0; //at initialisation every individual has viability = 0 (relative)
	}
}
void Population::computeSums(Parameters para, Individuals ind) {
	DelMut += ind.chromo.nMut;
	DelMutSd += (double)ind.chromo.nMut * (double)ind.chromo.nMut; //sum of squares
	W += ind.w;
	Wsd += ind.w * ind.w; //sum of squares
	if (ind.w > Wmax)(Wmax = ind.w); //save the highest fitness in current generation for scaled sampling of reproduction probability in loadEffect = 1.
	if (ind.w < Wmin)(Wmin = ind.w);
	Vneutral += ind.chromo.linkNeut[0] + ind.chromo.linkNeut[1];
	VneutralSD += (ind.chromo.linkNeut[0] * ind.chromo.linkNeut[0]) + (ind.chromo.linkNeut[1] * ind.chromo.linkNeut[1]); //sum of squares
}
void Population::computeStats(Parameters para) {
	DelMutSd = std::sqrt((DelMutSd - ((double)DelMut * (double)DelMut) / (double)N) / (double)N);
	meanDelMut = (double)DelMut / (double)N;
	Wsd = std::sqrt((Wsd - (W * W) / (double)N) / (double)N);
	W /= (double)N;
	Hsd = std::sqrt((Hsd - (Homoz * Homoz) / (double)N) / (double)N);
	Homoz /= (double)N;
	Vneutral = (VneutralSD - (Vneutral * Vneutral) / (2.0 * (double)para.K)) / (2.0 * (double)para.K);
}
void Population::set2zero(void) {
	//cout << "Max fitness type: " << Wmax << endl;
	DelMut = 0;
	DelMutSd = 0.0;
	W = 0.0;
	Wsd = 0.0;
	Homoz = 0.0; Hsd = 0.0;
	Vneutral = 0.0; VneutralSD = 0.0;
	Wmax = -100000000000000000000000.0;
	Wmin = 100000000000000000000000.0;

}

void Population::deleteAdults(void) {
	vector<Individuals>::iterator iter;

	if (!inds.empty()) {
		for (iter = inds.begin(); iter != inds.end(); iter++) iter->deleteInd();
		inds.clear();
	}
}

void Population::outPop(int r, int g, std::ofstream* out) {
	*out << r << "\t" << g << "\t" << W << "\t" << Vneutral << "\t";
	*out << endl;
}

void Population::addMutation(double pos, double ss, double hh)
{
	pop_muts muts;

	muts.h = hh;
	muts.s = ss;

	if (popMuts.find(pos) == popMuts.end()) {
		muts.count = 1;
		popMuts[pos] = muts;
	}
	else popMuts[pos].count++;
}
void Population::outMutations(int n, int r, int g, std::ofstream* out, Parameters para)
{
	map<double, pop_muts>::iterator iter;

	double count; 
	count = 0.0;

	for (iter = popMuts.begin(); iter != popMuts.end(); iter++) {
		*out << r << "\t" << g << "\t";
		*out << iter->second.s << "\t" << iter->second.h << "\t";
		*out << (double)iter->second.count / (2.0 * (double)para.K) << endl;
		count += 1.0;
	}
	//cout << "total number of mutations: " << count << endl;
}

