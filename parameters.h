#pragma once

#include <fstream>

using namespace std;

class Parameters
{
public:
	Parameters();
	~Parameters();
	//Simulation
	int SimNr;
	int rep; //replicates
	int gen; //generations
	int out_int; //generations interval for population and trait outputs 
	int out_start;
	int PopMut_interval;
	int SFS_sample; //Sample size for the site frequency spectrum (numbers of individuals). 
	//Landscape 
	double K; //population census size
	int x_max;
	int y_max;

	//Mutations
	int loadEffect; //0 = affects offsping survival; 1 = affects fecundity and fertilization or mating probability in males
	int initial_nMut; //initial number of mutations
	double R; //genome map length (see Roze & Rousset 2009, JEB) - corresponds to recombination rate
	double Ud; //mutation rate for mildly deleterious mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	double Ul; //mutation rate for lethal mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	double mean_sd; //mean selection coefficient for mildly deleterious mutations
	double mean_hd; //mean dominance coefficient for mildly deleterious mutations
	double sl; //selection coefficient for lethal mutations (fix - not sampled from a distribution)	
	double hl; //dominance coefficient for lethal mutations

	//DFE
	int sh_dist; // 0 = s,h sampled as in Spigler et al 2016; 1 = s sampled from unireal (0.0,1.0), h sampled from unireal (0.0,1.0); 2 = as 0, but with distributions for s and h switched. 3 = s sampled uniform, h = 0.5. 4 = s sampled from uniform, h = 0.2. 5 means s=h=0 for all mutations. 6: h=0.5, s gamma dist (0.25, 0.05).
	//sh_dist = 7: samples s coefficients as scaled values (as 4Nes for Ne = N) using gamma dists from Tataru et al. 2017. (re-scale values by 1/2*4Ne) before adding to genome. sh =8, defines gamma dist from mean and shape parameter, scale calculated (Tataru et al. 2017). 
	int neut_pos; //0: Neutral linked locus is centrometic. 1: Neutral linked locus is telomeric.
	double mean; // mean of gamma distribution of selection coefficients
	double shape; //shape of gamma distribution for custom option, sh_dist =8.

	void outPara(string dir); //Parameter output
};