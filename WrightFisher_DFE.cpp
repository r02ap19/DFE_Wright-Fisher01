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
	para.mean = std::atof(argv[2]);
	para.shape = std::atof(argv[3]);
	para.assumed = std::atoi(argv[4]);

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
		outNeutral_SFS_sample_header();
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
	if (neutralSFSsample.is_open()) neutralSFSsample.close();
}

void initialisation() {

	pop = new Population * *[1];
	for (int j = 0; j < 1; j++) {
		pop[j] = new Population * [1];
		//do not intialise the single pop. - do it only when there is a population
		for (int jj = 0; jj < 1; jj++) pop[j][jj] = NULL;
	}
	pop[0][0] = new Population(0, 0);
	pop[0][0]->initialise_pop(para.K, k, para, s_mild, position, linked_neutral);
}

void reproduction_1(void) { //fertility selection

	int n_mut;
	double p_surv; //survival probability 
	Individuals* ind;
	int mate1, mate2, counter;
	double pr1, pr2;
	int gg; // count number of sampled individuals for SFS
	gg = 0;
	std::poisson_distribution<> n_mildmut(para.Ud);
	std::poisson_distribution<> n_lethal(para.Ul);

	double f_Wmax, m_Wmax, f_Wmin, m_Wmin, count_f, count_m;

	f_Wmax = 0.0;
	m_Wmax = 0.0;
	f_Wmin = 1.0;
	m_Wmin = 1.0;

	count_f = 0.0;
	count_m = 0.0;
	counter = 0;

	vector<Individuals>::iterator it;

	//Loop through the landscape
	for (int x = 0; x < para.x_max; x++) {
		for (int y = 0; y < para.y_max; y++) {
			if (pop[x][y] != NULL) {

				//summing male and female fitness to calculate mean 
				for (int i = 0; i < (para.K -1); i++) {
					count_f += pop[x][y]->inds[i].w;
				}

				count_f = (count_f / para.K);
				
				//Re-calculating max female and male fitness
				for (int i = 0; i < ((para.K / 2) - 1); i++) {
					if (pop[x][y]->inds[i].w > f_Wmax) (f_Wmax = pop[x][y]->inds[i].w);
					if (pop[x][y]->inds[i].w < f_Wmin) (f_Wmin = pop[x][y]->inds[i].w);
				}
				for (int i = ((para.K) / 2); i < (para.K - 1); i++) {
					if (pop[x][y]->inds[i].w > m_Wmax) (m_Wmax = pop[x][y]->inds[i].w);
					if (pop[x][y]->inds[i].w < m_Wmin) (m_Wmin = pop[x][y]->inds[i].w);
				}

				//cout << "max/min female fitness: " << f_Wmax << "\t"<< f_Wmin << endl;
				//cout << "max/min male fitness: " << m_Wmax << "\t" << m_Wmin << endl;


				//scaling after population mean fitness
				for (int i = 0; i < (para.K); i++) {
					if (i < (para.K / 2)) {
						pop[x][y]->inds[i].w = (pop[x][y]->inds[i].w/ (count_f));
					}
					else {
						pop[x][y]->inds[i].w = (pop[x][y]->inds[i].w / (count_f));
					}
				}

				//Re-calculating max female and male fitness
				for (int i = 0; i < ((para.K / 2)-1); i++) {
					if (pop[x][y]->inds[i].w > f_Wmax) (f_Wmax = pop[x][y]->inds[i].w);
					if (pop[x][y]->inds[i].w < f_Wmin) (f_Wmin = pop[x][y]->inds[i].w);
				}
				for (int i = ((para.K) / 2); i < (para.K - 1); i++) {
					if (pop[x][y]->inds[i].w > m_Wmax) (m_Wmax = pop[x][y]->inds[i].w);
					if (pop[x][y]->inds[i].w < m_Wmin) (m_Wmin = pop[x][y]->inds[i].w);
				}

				//cout << "After: max/min female fitness: " << f_Wmax << "\t" << f_Wmin << endl;
				//cout << "After: max/min male fitness: " << m_Wmax << "\t" << m_Wmin << endl;

				//sample index for male and female parents (creates seperate sexes, removes self-fertilization as in SFS_CODE) //31/07/31
				std::uniform_int_distribution<> fe_WF(0, ((pop[x][y]->N)/2)-1);
				std::uniform_int_distribution<> ma_WF(((pop[x][y]->N) / 2), (pop[x][y]->N - 1));

				//clear Map of population's deleterious mutations
				if (!pop[x][y]->popMuts.empty()) pop[x][y]->popMuts.clear();

				//MATING
				for (int k = 0; k < (pop[x][y]->N); k++) {
					//pick two random individuals as mates (possibly the same individual)
					mate1 = fe_WF(rdgen);
					mate2 = ma_WF(rdgen);

					//Let the probability of reproducing after being sampled be relative to fitness

					std::uniform_real_distribution<> scaled_f(0, f_Wmax);
					std::uniform_real_distribution<> scaled_m(0, m_Wmax);

					if (scaled_f(rdgen) <= pop[x][y]->inds[mate1].w)(pop[x][y]->inds[mate1].reproduce = true);
					else (pop[x][y]->inds[mate1].reproduce = false);
					if (scaled_m(rdgen) <= pop[x][y]->inds[mate2].w)(pop[x][y]->inds[mate2].reproduce = true);
					else (pop[x][y]->inds[mate2].reproduce = false);

					if (pop[x][y]->inds[mate1].reproduce && pop[x][y]->inds[mate2].reproduce) {
						ind = new Individuals(para);
						gg += 1;
						if (para.neut_pos == 0) inheritance_0(ind, pop[x][y]->inds[mate1], pop[x][y]->inds[mate2]);
						else inheritance_1(ind, pop[x][y]->inds[mate1], pop[x][y]->inds[mate2]);

						if (para.neutral_genome == true) {
							neutral_inheritance(ind, pop[x][y]->inds[mate1], pop[x][y]->inds[mate2]);
							n_mut = n_mildmut(rdgen);
							if (n_mut > 0) {
								ind->neutral_genome_mut(n_mut, position, uniform, para);
							}
						}

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
							ind->neutralSFS_sample(g, gg, &neutralSFSsample);
						}
						if (ind->w > 0) {
							if (ind->w == 1)( counter += 1 );
							if (ind->alive && k < para.K) {
								pop[x][y]->tmp_inds.push_back(*ind);
							}
							else {
								ind->deleteInd();
							}
							delete ind;
							
						}
						else {
							ind->deleteInd();
							k--;
						}

					}
					else k--;
				}
				//cout << "number of perfect inds generation:" << endl;
				//cout << counter << endl;
				cout << "varience at the neutral linked locus: " << pop[x][y]->Vneutral << endl;
				//cout << "max female fitness: " << f_Wmax << endl;
				//cout << "max male fitness: " << m_Wmax << endl;
				//cout << "pop size: " << pop[x][y]->N << endl;
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
						ind->neutralSFS_sample(g, gg, &neutralSFSsample);
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

void neutral_inheritance(Individuals* pup, Individuals mom, Individuals dad) {
	int rdn, rdn2, pos, pos2;
	int hom;
	int n_crossovers;
	double cross;
	std::map<double, mutation>::iterator iter, iter2;
	std::set <double> recomSites;
	std::set<double>::iterator itercross;

	if (mom.neutral_genome.nMut > 0) {
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
		iter = mom.neutral_genome.mutations.begin();

		//no mutations before cross-overs positions
		while (n_crossovers > 0 && *itercross < iter->first) {
			if (hom == 0) hom++;
			else hom--;
			itercross++;
			n_crossovers--;
		}
		for (iter = mom.neutral_genome.mutations.begin(); iter != mom.neutral_genome.mutations.end(); iter++) {
			//cross-overs
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}

			//if mutation is on the right homologue inherit it, otherwise ignore it
			if (iter->second.homol == hom || iter->second.homol == 2) {
				pup->neutral_genome.mutations[iter->first] = iter->second;
				pup->neutral_genome.mutations[iter->first].homol = 0; //inherit first homologue from mon

				pup->neutral_genome.nMut++;

			}
		}
		if (!recomSites.empty()) recomSites.clear();
	}

	if (dad.neutral_genome.nMut > 0) {
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
		iter = dad.neutral_genome.mutations.begin();

		while (n_crossovers > 0 && *itercross < iter->first) {
			if (hom == 0) hom++;
			else hom--;
			itercross++;
			n_crossovers--;
		}
		for (iter = dad.neutral_genome.mutations.begin(); iter != dad.neutral_genome.mutations.end(); iter++) {
			//crossovers
			while (n_crossovers > 0 && *itercross < iter->first) {
				if (hom == 0) hom++;
				else hom--;
				itercross++;
				n_crossovers--;
			}
			if (iter->second.homol == hom || iter->second.homol == 2) {
				iter2 = pup->neutral_genome.mutations.find(iter->first);
				//if mutation is already present --> it is homozygous
				if (iter2 != pup->neutral_genome.mutations.end()) {
					iter2->second.homol = 2; //mutation is homozygote
					pup->neutral_genome.Nho++;
				}
				else { //mutation is heterozygote
					pup->neutral_genome.mutations[iter->first] = iter->second;
					pup->neutral_genome.mutations[iter->first].homol = 1;
					pup->neutral_genome.nMut++;
				}
			}

		}
		if (!recomSites.empty()) recomSites.clear();
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
				pup->w -= (iter->second.h * iter->second.s);

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
					pup->w += (iter2->second.h * iter2->second.s);
					pup->w -= (2*(iter2->second.s) + pow(iter2->second.s, 2));
					//pup->w -= (2 * (iter2->second.s));
				}
				else { //mutation is heterozygote
					pup->chromo.mutations[iter->first] = iter->second;
					pup->chromo.mutations[iter->first].homol = 1;
					pup->chromo.nMut++;
					//fitness effect
					pup->w -= (iter->second.h * iter->second.s);
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
				pup->w -= (iter->second.h * iter->second.s);
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
					pup->w += (iter2->second.h * iter2->second.s);
					pup->w -= (2*(iter2->second.s) + pow(iter2->second.s,2));
					//pup->w -= (2 * (iter2->second.s));
				}
				else { //mutation is heterozygote
					pup->chromo.mutations[iter->first] = iter->second;
					pup->chromo.mutations[iter->first].homol = 1;
					pup->chromo.nMut++;
					//fitness effect
					pup->w -= (iter->second.h * iter->second.s);
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
				pop[x][y]->N++;
				pop[x][y]->inds.push_back(*iter);
			}

			for (iter = pop[x][y]->tmp_inds.begin(); iter != pop[x][y]->tmp_inds.end(); iter++) {
				//sum deleterious mutation to calculate population frequencies
				if (g % para.PopMut_interval == 0) {
					for (iter2 = iter->chromo.mutations.begin(); iter2 != iter->chromo.mutations.end(); iter2++) {
						pop[x][y]->addMutation(iter2->first, iter2->second.s, iter2->second.h);
						if (iter2->second.homol == 2) pop[x][y]->addMutation(iter2->first, iter2->second.s, iter2->second.h);
					}
				}
			}

			for (iter = pop[x][y]->inds.begin(); iter != pop[x][y]->inds.end(); iter++) {
				pop[x][y]->computeSums(para, *iter); //compute normalized Wmin & Wmax; boundries for reproduction sampling in next generation.
				//cout << iter->w << endl;
			}

			pop[x][y]->computeStats(para); //compute normalized mean fitness for output

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

void outNeutral_SFS_sample_header(void){
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_neutralSFSsample.txt";
	neutralSFSsample.open(name.c_str());
	neutralSFSsample << "gen\ts\th\tpos\tind\thom" << endl;
}

void outPop_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Pops.txt";
	pops.open(name.c_str());
	pops << "rep\tgen\t\tW\tVN";
	pops << endl;
}

void outPopMut_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_PopMut.txt";
	popmut.open(name.c_str());
	popmut << "rep\tgen\ts\th\tfreq" << endl;
}

void outSFS_sample_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_SFSsample.txt";
	SFSsample.open(name.c_str());
	SFSsample << "gen\ts\th\tpos\tind\thom" << endl;
}

