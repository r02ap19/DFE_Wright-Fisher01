#include "parameters.h"

Parameters::Parameters() {

	//Simulation
	SimNr = 99001;
	rep = 1; //replicates
	gen = 10000; //generations
	out_int = 1; //generation interval for output (population)
	out_start = 1; //output start generation
	PopMut_interval = 10000;
	SFS_sample = 10;

	//Landscape
	K = 500; // Population census size
	x_max = 1;
	y_max = 1;

	//Deleterious mutations
	loadEffect = 1; //0 = affects offsping survival; 1 = affects probability of reproducing;
	initial_nMut = 0;
	R = 10.0;
	Ud = 1.0; //mutation rate for mildly deleterious mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	Ul = 0.0000000000001; //mutation rate for lethal mutations / diploid genome / generation (see Spigler et al. 2016, Evolution)
	mean_sd = 0.05; //mean selection coefficient for mildly deleterious mutations (for some values of sh_dist)
	mean_hd = 0.3; //mean dominance coefficient for mildly deleterious mutations (for some values of sh_dist)
	sl = 1.0; //selection coefficient for lethal mutations (fix - not sampled from a distribution)
	hl = 0.02; //mean dominance coefficient for lethal mutations

	//DFE Simulations
	sh_dist = 8;
	neut_pos = 0;
	mean = 50000000000000000000.0;
	shape = 0.4;
	neutral_genome = true;
	assumed_Ne = 500;
}

void Parameters::outPara(string name) {

	ofstream out;

	out.open(name.c_str());

	//Simulation
	out << "SimNr\t" << SimNr << endl;
	out << "rep\t" << rep << endl; //replicates
	out << "gen\t" << gen << endl; //generations
	out << "out_int\t" << out_int << endl;
	out << "out_start\t" << out_start << endl;
	out << "PopMut_interval\t" << PopMut_interval << endl;

	//Landscape	
	out << "K\t" << K << endl; //N (census population size)
	out << "x_max\t" << x_max << endl;
	out << "y_max\t" << y_max << endl;

	//Mutations
	out << "load_effect\t" << loadEffect << endl;
	out << "initial_Nmut\t" << initial_nMut << endl;
	out << "genome_map_length\t" << R << endl;
	out << "Ud\t" << Ud << endl;
	out << "Ul\t" << Ul << endl;
	out << "mean_selection_mild_mut\t" << mean_sd << endl;
	out << "mean_dominance_mild_mut\t" << mean_hd << endl;
	out << "selection_lethal_mut\t" << sl << endl;
	out << "dominance_lethal_mut\t" << hl << endl;

	//DFE
	out << "sh_dist\t" << sh_dist << endl;
	out << "neut_pos\t" << neut_pos << endl;
	out << "SFS_sample\t" << SFS_sample << endl;
	out << "mean\t" << mean << endl;
	out << "shape\t" << shape << endl;
	out << "neutral_genome\t" << neutral_genome << endl;
	out << "Assumed Ne\t" << assumed_Ne << endl;

	out.close();
}
//------------------------------
Parameters::~Parameters() {
}