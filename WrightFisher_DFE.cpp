#include "WrightFisher_DFE.h"

#if CLUSTER
int main(int argc, char* argv[])
{
	// Get the current directory.
	char* buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; //Current directory path
	dirOut = dir + "Outputs/"; //Outpus folder path


	para.SimNr = std::atoi(argv[1]);
	para.Ud = std::atof(argv[2]);
	para.sh_dist = std::atoi(argv[3]);

	RunModel();

	cout << "Simulation completed" << endl;

	return 0;
}
#else
int _tmain(int argc, _TCHAR* argv[])
{
	// Get the current directory.
	char* buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; //Current directory path
	dirOut = dir + "Outputs\\"; //Outpus folder path

	extime = clock();

	RunModel();

	std::cout << "Simulation completed" << endl;

	return 0;
}
#endif
const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}

void RunModel(void) {
	std::cout << "Simulation nr. " << para.SimNr << endl;
	string name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Para.txt";
	para.outPara(name);

	if (para.out_int > 0) outPop_header();
	if (para.PopMut_interval > 0) {
		outPopMut_header();
		outSFS_sample_header();
	}
	for (r = 0; r < para.rep; r++) {
		std::cout << "rep = " << r << "==================" << endl;

		//Initialisation
		initialisation();

		//Loop through GENERATIONS ------------------------------------
		for (g = 0; g < para.gen; g++) {

			if (g % 50 == 0) {
				std::cout << "gen = " << g << endl;
				extime = clock() - extime;
				std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
				extime = clock();
			}

			//Reproduction
			if (para.loadEffect == 0) reproduction_0();
			else reproduction_1();

			housekeeping();
		}

		//Delete populations and the landscape
		for (int x = 0; x < para.x_max; x++) {
			for (int y = 0; y < para.y_max; y++) {
				if (pop[x][y] != NULL) {
					pop[x][y]->deleteAdults();
					delete pop[x][y]; pop[x][y] = NULL;
				}
			}
			delete[] pop[x]; pop[x] = NULL;
		}
		delete[] pop; pop = NULL;
	}

	if (pops.is_open()) pops.close();
	if (popmut.is_open()) popmut.close();
	if (SFSsample.is_open()) SFSsample.close();
}

void initialisation() {

	pop = new Population * *[1];
	for (int j = 0; j < 1; j++) {
		pop[j] = new Population * [1];
		//do not intialise the single pop. - do it only when there is a population
		for (int jj = 0; jj < 1; jj++) pop[j][jj] = NULL;
	}
	pop[0][0] = new Population(0,0);
	pop[0][0]->initialise_pop(para.K, k, para, s_mild, position, linked_neutral);
}

void reproduction_1(void) { //fertility selection

	int n_mut;
	double p_surv; //survival probability 
	Individuals* ind;
	int mate1, mate2;
	double pr1, pr2;
	int gg; // count number of sampled individuals for SFS
	gg = 0;
	std::poisson_distribution<> n_mildmut(para.Ud);
	std::poisson_distribution<> n_lethal(para.Ul);


	//Loop through the landscape
	for (int x = 0; x < para.x_max; x++) {
		for (int y = 0; y < para.y_max; y++) {
			if (pop[x][y] != NULL) {

				//sample index for individual
				std::uniform_int_distribution<> fe_WF(0, (pop[x][y]->N - 1));
				
				//clear Map of population's deleterious mutations
				if (!pop[x][y]->popMuts.empty()) pop[x][y]->popMuts.clear();


				//MATING
				for (int k = 0; k < (pop[x][y]->N); k++) {
					//pick two random individuals as mates (possibly the same individual)
					mate1 = fe_WF(rdgen);
					mate2 = fe_WF(rdgen);

					//Let the probability of reproducing after being sampled be equal to fitness

					std::uniform_real_distribution<> scaled(0, pop[x][y]->Wmax);

					if (scaled(rdgen) < pop[x][y]->inds[mate1].w)(pop[x][y]->inds[mate1].reproduce = true);
					else (pop[x][y]->inds[mate1].reproduce = false);
					if (scaled(rdgen) < pop[x][y]->inds[mate2].w)(pop[x][y]->inds[mate2].reproduce = true);
					else (pop[x][y]->inds[mate2].reproduce = false);

					if (pop[x][y]->inds[mate1].reproduce && pop[x][y]->inds[mate2].reproduce) {
						ind = new Individuals(para);
						gg += 1;
						if (para.neut_pos == 0) inheritance_0(ind, pop[x][y]->inds[mate1], pop[x][y]->inds[mate2]);
						else inheritance_1(ind, pop[x][y]->inds[mate1], pop[x][y]->inds[mate2]);

						//Mutation in neutral linked locus
						int hom = Bern(rdgen);
						ind->chromo.linkNeut[hom] += linked_neutral(rdgen);

						//mildly deleterius mutations
						n_mut = n_mildmut(rdgen);
						if (n_mut > 0) ind->delet_mutation(n_mut, k, position, s_mild, uniform, para);

						//lethal mutations
						//n_mut = n_lethal(rdgen);
						//if (n_mut > 0) ind->delet_mutation(n_mut, position, para.sl, para.hl);

						if (gg <= para.SFS_sample && g >= para.PopMut_interval) { //sample genomes of individuals for contructions of SFS
							ind->SFS_sample(g, gg, &SFSsample);
						}

						if (ind->w > 0.0) {

							if (ind->alive && k < para.K) {
								pop[x][y]->tmp_inds.push_back(*ind);
							}
							else {
								ind->deleteInd();
							}
						}
						else {
							ind->deleteInd();
						}
						delete ind;
					}
					else k--;
				}
				pop[x][y]->N = 0;
				pop[x][y]->deleteAdults();
			}
		}
	}
}

void reproduction_0(void) {
	int n_mut;
	double p_surv; //survival probability 
	Individuals* ind;
	int mate1, mate2;
	double pr1, pr2;
	int gg; // count number of sampled individuals for SFS
	gg = 0;

	std::poisson_distribution<> n_mildmut(para.Ud);
	std::poisson_distribution<> n_lethal(para.Ul);

	//Loop through the landscape
	for (int x = 0; x < para.x_max; x++) {
		for (int y = 0; y < para.y_max; y++) {
			if (pop[x][y] != NULL) {

				//sample index for individual
				std::uniform_int_distribution<> fe_WF(0, (pop[x][y]->N - 1));

				//clear Map of population's deleterious mutations
				if (!pop[x][y]->popMuts.empty()) pop[x][y]->popMuts.clear();

				//MATING
				for (int k = 0; k < (pop[x][y]->N); k++) {
					//pick two random individuals as mates (possibly the same individual)
					mate1 = fe_WF(rdgen);
					mate2 = fe_WF(rdgen);

					ind = new Individuals(para);
					gg += 1;
					if (para.neut_pos == 0) inheritance_0(ind, pop[x][y]->inds[mate1], pop[x][y]->inds[mate2]);
					else inheritance_1(ind, pop[x][y]->inds[mate1], pop[x][y]->inds[mate2]);

					//Mutation in neutral linked locus
					int hom = Bern(rdgen);
					ind->chromo.linkNeut[hom] += linked_neutral(rdgen);

					//mildly deleterius mutations
					n_mut = n_mildmut(rdgen);
					if (n_mut > 0) ind->delet_mutation(n_mut, k, position, s_mild, uniform, para);

					//lethal mutations
					//n_mut = n_lethal(rdgen);
					//if (n_mut > 0) ind->delet_mutation(n_mut, position, para.sl, para.hl);

					if (gg <= para.SFS_sample && g >= para.PopMut_interval) { //sample genomes of individuals for contructions of SFS
						ind->SFS_sample(g, gg, &SFSsample);
					}

					std::uniform_real_distribution<> scaled(0, pop[x][y]->Wmax);

					if (ind->w > scaled(rdgen) && k < para.K) {
						pop[x][y]->tmp_inds.push_back(*ind);
					}
					else {
						ind->deleteInd();
						k--;
					}
					delete ind;
				}
				pop[x][y]->N = 0;
				pop[x][y]->deleteAdults();
			}
		}
	}
}


void inheritance_0(Individuals* pup, Individuals mom, Individuals dad) {
	int rdn, rdn2, pos, pos2;
	int hom;
	int n_crossovers;
	double cross;
	std::map<double, mutation>::iterator iter, iter2;
	std::set <double> recomSites;
	std::set<double>::iterator itercross;

	bool check_neutral;

	//Recombination of deleterious mutations (see Roze & Rousset 2009, JEB - Appendix 2)
	//Inherit homologue 1 from mate1 (mom)
	check_neutral = true;

	if (mom.chromo.nMut > 0) {
		//sample no. of crossovers
		n_crossovers = crossn(rdgen);
		//sample crossover positions
		for (int i = 0; i < n_crossovers; i++) {
			cross = position(rdgen);
			recomSites.insert(cross);
		}
		itercross = recomSites.begin(); //iterator through crossover positions

		//sample starting homologue
		hom = Bern(rdgen);
		iter = mom.chromo.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {
			if (*itercross > (para.R / 2.0)) { 
				pup->chromo.linkNeut[0] = mom.chromo.linkNeut[hom];
				check_neutral = false;
			}
			if (hom == 0) hom++;
			else hom--;
			itercross++;
			n_crossovers--;
		}
		for (iter = mom.chromo.mutations.begin(); iter != mom.chromo.mutations.end(); iter++) {
			//cross-overs
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}
			if (check_neutral && *itercross > (para.R / 2.0)) {
				pup->chromo.linkNeut[0] = mom.chromo.linkNeut[hom];
				check_neutral = false;
			}

			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {
				pup->chromo.mutations[iter->first] = iter->second;
				pup->chromo.mutations[iter->first].homol = 0; //inherit first homologue from mon

				pup->chromo.nMut++;
				//calculate fitness considering the mutation as heterozygote
				pup->w *= (1.0 - iter->second.h * iter->second.s);		

			}
		}
		if (!recomSites.empty()) recomSites.clear();
		if (check_neutral) {
			pup->chromo.linkNeut[0] = mom.chromo.linkNeut[hom];
			check_neutral = false;
		}
	}
	else {
		hom = Bern(rdgen);
		pup->chromo.linkNeut[0] = mom.chromo.linkNeut[hom];
		check_neutral = false;
	}

	check_neutral = true;

	if (dad.chromo.nMut > 0) {
		//sample no. of crossovers
		n_crossovers = crossn(rdgen);
		//sample crossover positions
		for (int i = 0; i < n_crossovers; i++) {
			cross = position(rdgen);
			recomSites.insert(cross);
		}
		itercross = recomSites.begin(); //iterator through crossover positions
										//sample starting homologue
		hom = Bern(rdgen);
		iter = dad.chromo.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {
			if (*itercross > (para.R / 2.0)) {
				pup->chromo.linkNeut[1] = dad.chromo.linkNeut[hom];
				check_neutral = false;
			}
			if (hom == 0) hom++;
			else hom--;
			itercross++;
			n_crossovers--;
		}
		for (iter = dad.chromo.mutations.begin(); iter != dad.chromo.mutations.end(); iter++) {
			//crossovers
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}
			if (check_neutral && *itercross > (para.R / 2.0)) {
				pup->chromo.linkNeut[1] = dad.chromo.linkNeut[hom];
				check_neutral = false;
			}
			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {
				iter2 = pup->chromo.mutations.find(iter->first);
				//if mutation is already present --> it is homozygous
				if (iter2 != pup->chromo.mutations.end()) {
					iter2->second.homol = 2; //mutation is homozygote
					pup->chromo.Nho++;
					//change fitness effect
					pup->w /= (1.0 - iter2->second.h * iter2->second.s);
					pup->w *= (1.0 - iter2->second.s);
				}
				else { //mutation is heterozygote
					pup->chromo.mutations[iter->first] = iter->second;
					pup->chromo.mutations[iter->first].homol = 1;
					pup->chromo.nMut++;
					//fitness effect
					pup->w *= (1.0 - iter->second.h * iter->second.s);
				}
			}
		}
		if (!recomSites.empty()) recomSites.clear();
		if (check_neutral) {
			pup->chromo.linkNeut[1] = dad.chromo.linkNeut[hom];
			check_neutral = false;
		}
	}
	else {
		hom = Bern(rdgen);
		pup->chromo.linkNeut[1] = dad.chromo.linkNeut[hom];
		check_neutral = false;
	}
}


void inheritance_1(Individuals* pup, Individuals mom, Individuals dad) {
	int rdn, rdn2, pos, pos2;
	int hom;
	int n_crossovers;
	double cross;
	std::map<double, mutation>::iterator iter, iter2;
	std::set <double> recomSites;
	std::set<double>::iterator itercross;

	bool check_neutral;

	//Recombination of deleterious mutations (see Roze & Rousset 2009, JEB - Appendix 2)
	//Inherit homologue 1 from mate1 (mom)
	check_neutral = true;

	if (mom.chromo.nMut > 0) {
		//sample no. of crossovers
		n_crossovers = crossn(rdgen);
		//sample crossover positions
		for (int i = 0; i < n_crossovers; i++) {
			cross = position(rdgen);
			recomSites.insert(cross);
		}
		itercross = recomSites.begin(); //iterator through crossover positions

		//sample starting homologue
		hom = Bern(rdgen);
		iter = mom.chromo.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {
			if (itercross == recomSites.end()) { 
				pup->chromo.linkNeut[0] = mom.chromo.linkNeut[hom];
				check_neutral = false;
			}

			if (hom == 0) hom++;
			else hom--;
			itercross++;
			n_crossovers--;
		}
		for (iter = mom.chromo.mutations.begin(); iter != mom.chromo.mutations.end(); iter++) {
			//cross-overs
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}
			if (check_neutral && itercross == recomSites.end()) {
				pup->chromo.linkNeut[0] = mom.chromo.linkNeut[hom];
				check_neutral = false;
			}
			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {
				pup->chromo.mutations[iter->first] = iter->second;
				pup->chromo.mutations[iter->first].homol = 0; //inherit first homologue from mon

				pup->chromo.nMut++;
				//calculate fitness considering the mutation as heterozygote
				pup->w *= (1.0 - iter->second.h * iter->second.s);
			}
		}
		if (!recomSites.empty()) recomSites.clear();
		if (check_neutral) {
			pup->chromo.linkNeut[0] = mom.chromo.linkNeut[hom];
			check_neutral = false;
		}
	}
	else {
		hom = Bern(rdgen);
		pup->chromo.linkNeut[0] = mom.chromo.linkNeut[hom];
		check_neutral = false;
	}

	check_neutral = true;

	if (dad.chromo.nMut > 0) {
		//sample no. of crossovers
		n_crossovers = crossn(rdgen);
		//sample crossover positions
		for (int i = 0; i < n_crossovers; i++) {
			cross = position(rdgen);
			recomSites.insert(cross);
		}
		itercross = recomSites.begin(); //iterator through crossover positions
	    //sample starting homologue
		hom = Bern(rdgen);
		iter = dad.chromo.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {
			if (itercross == recomSites.end()) {
				pup->chromo.linkNeut[1] = dad.chromo.linkNeut[hom];
				check_neutral = false;
			}
			if (hom == 0) hom++;
			else hom--;
			itercross++;
			n_crossovers--;
		}
		for (iter = dad.chromo.mutations.begin(); iter != dad.chromo.mutations.end(); iter++) {
			//crossovers
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}
			if (check_neutral && itercross == recomSites.end()) {
				pup->chromo.linkNeut[1] = dad.chromo.linkNeut[hom];
				check_neutral = false;
			}
			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {
				iter2 = pup->chromo.mutations.find(iter->first);
				//if mutation is already present --> it is homozygous
				if (iter2 != pup->chromo.mutations.end()) {
					iter2->second.homol = 2; //mutation is homozygote
					pup->chromo.Nho++;
					//change fitness effect
					pup->w /= (1.0 - iter2->second.h * iter2->second.s);
					pup->w *= (1.0 - iter2->second.s);
				}
				else { //mutation is heterozygote
					pup->chromo.mutations[iter->first] = iter->second;
					pup->chromo.mutations[iter->first].homol = 1;
					pup->chromo.nMut++;
					//fitness effect
					pup->w *= (1.0 - iter->second.h * iter->second.s);
				}
			}
		}
		if (!recomSites.empty()) recomSites.clear();
		if (check_neutral) {
			pup->chromo.linkNeut[1] = dad.chromo.linkNeut[hom];
			check_neutral = false;
		}
	}
	else {
		hom = Bern(rdgen);
		pup->chromo.linkNeut[1] = dad.chromo.linkNeut[hom];
		check_neutral = false;
	}
}

void housekeeping() {
	vector<Individuals>::iterator iter;
	std::map<double, mutation>::iterator iter2;

	for (int x = 0; x < para.x_max; x++) {
		for (int y = 0; y < para.y_max; y++) {
			pop[x][y]->set2zero();

			for (iter = pop[x][y]->tmp_inds.begin(); iter != pop[x][y]->tmp_inds.end(); iter++) {
				pop[x][y]->inds.push_back(*iter);
				pop[x][y]->N++;
				pop[x][y]->computeSums(para, *iter);

				//sum deleterious mutation to calculate population frequencies
				if (g % para.PopMut_interval == 0) {
					for (iter2 = iter->chromo.mutations.begin(); iter2 != iter->chromo.mutations.end(); iter2++) {
						pop[x][y]->addMutation(iter2->first, iter2->second.s, iter2->second.h);
						if (iter2->second.homol == 2) pop[x][y]->addMutation(iter2->first, iter2->second.s, iter2->second.h);
					}
				}
			}
			pop[x][y]->computeStats(para);

			//output populations
			if (g > para.out_start - 1 && (g % para.out_int == 0)) {
				pop[x][y]->outPop(r, g, &pops);
			}

			//output frequency of deleterious mutations
			if (g > para.out_start - 1 && g > 0 && g % para.PopMut_interval == 0) pop[x][y]->outMutations(para.K, r, g, &popmut, para);
			
			pop[x][y]->tmp_inds.clear();
		}
	}
}

void outPop_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Pops.txt";
	pops.open(name.c_str());
	pops << "rep\tgen\t\tW\tVN";
	pops << endl;
}

void outPopMut_header(void){
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_PopMut.txt";
	popmut.open(name.c_str());
	popmut << "rep\tgen\ts\th\tfreq" << endl;
}

void outSFS_sample_header(void){
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_SFSsample.txt";
	SFSsample.open(name.c_str());
	SFSsample << "gen\ts\th\tpos\tind\thom" << endl;
}