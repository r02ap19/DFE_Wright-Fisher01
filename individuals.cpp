#include "individuals.h"

Parameters para2;

std::random_device rd2;
std::mt19937 gen(rd2());

std::bernoulli_distribution findhom(0.5); //for sampling the homologue

Individuals::Individuals(Parameters para) {
	alive = true;
	reproduce = true;
	w = 1.0;
}

Individuals::~Individuals() {
}

void Individuals::deleteInd(void) {
	chromo.deleteChromo();
}

void Individuals::initialise(double kk, Parameters para, std::gamma_distribution<> finds, std::uniform_real_distribution<> findpos, std::normal_distribution<> neut) {

	int hom;
	double pos;
	double ss, ss1;
	double hh;

	ss = neut(gen);
	ss1 = neut(gen);
	chromo.InitNeutral(ss, ss1);

	if (para.initial_nMut > 0) {
		for (int i = 0; i < para.initial_nMut; i++) {
			//sample homologue
			hom = findhom(gen);
			//sample selection coefficient & dominance coefficient
			if (para.sh_dist == 0) {
				ss = finds(gen);
				std::uniform_real_distribution<> findh(0.0, std::exp(-kk * ss));
				hh = findh(gen);
			}
			if (para.sh_dist == 1) {
				std::uniform_real_distribution<> finds2(0.0, 1.0);
				ss = finds2(gen);
				std::uniform_real_distribution<> h_rec(0.0, 1.0);
				hh = h_rec(gen);
			}
			if (para.sh_dist == 2) {
				double mean_sd = para.mean_hd;
				double mean_hd = para.mean_sd;
				double k = -log(2.0 * mean_hd) / mean_sd;
				gamma_distribution<double> h_cof(1.0, mean_hd);
				hh = h_cof(gen);
				uniform_real_distribution<double> s_cof(0.0, exp(-k * hh));
				ss = s_cof(gen);
			}
			if (para.sh_dist == 3) {
				std::uniform_real_distribution<> finds2(0.0, 1.0);
				ss = finds2(gen);
				hh = 0.5;
			}
			if (para.sh_dist == 4) {
				std::uniform_real_distribution<> finds2(0.0, 1.0);
				ss = finds2(gen);
				hh = 0.2;
			}
			if (para.sh_dist == 5) {
				ss = 0.0;
				hh = 0.0;
			}
			if (para.sh_dist == 6) {
				gamma_distribution<double> ss_cof(0.25, 0.05);
				ss = ss_cof(gen);
				hh = 0.5;
			}
			if (para.sh_dist == 7) {
				gamma_distribution<double> ss_cof(0.4, 25000);
				ss = ss_cof(gen);
				ss = ss / double(0.5 * 4 * para.K);
				hh = 0.5;
				while (ss > 1.0)
				{
				ss = ss_cof(gen);
				ss = ss / double(0.5 * 4 * para.K);
				}
			}
			if (para.sh_dist == 8) {
				gamma_distribution<double> ss_cof(para.shape, (para.mean / para.shape));
				ss = ss_cof(gen);
				ss = ss / double(0.5 * 4 * para.K);
				hh = 0.5;
				if (ss > 1.0)(ss = 1.0);
				//while (ss > 1.0)
				//{
				//	ss = ss_cof(gen);
				//	ss = ss / double(0.5 * 4 * para.K);
				//}
			}
			//sample position
			pos = findpos(gen);
			//add mutation
			w *= chromo.addDelMutation(hom, pos, ss, hh);
		}//hello Anders!
	}
}

void Individuals::benef_mutation(int nmut, double s, std::uniform_real_distribution<> findpos) {
	int hom;
	double pos;
	double ss;
	std::bernoulli_distribution findhom(0.5); //for sampling the homologue

	for (int i = 0; i < nmut; i++) {
		//sample homologue
		hom = findhom(gen);
		//sample selection coefficient
		ss = s;
		//sample position
		pos = findpos(gen);
		//add mutation
		w *= chromo.addDelMutation(hom, pos, ss, 1.0); //assume complete dominance of beneficial mutations
	}
}
//--------------------------------------------------------------------------
void Individuals::back_mutation(int nmut)
{
	int rdn;
	std::map<double, mutation, std::less<double>>::iterator it;

	std::bernoulli_distribution hom(0.5);

	for (int i = 0; i < nmut; i++) {
		std::uniform_int_distribution<> samp(0, chromo.mutations.size() - 1);
		rdn = samp(gen);
		it = chromo.mutations.begin();
		std::advance(it, rdn);

		if (it->second.homol == 2) { //mutation is homozygote
			it->second.homol = hom(gen);
			chromo.Nho--;
			//update fitness
			w /= (1.0 - it->second.s);
			w *= (1.0 - it->second.h * it->second.s);
		}
		else { //mutation is heterozygote
			chromo.mutations.erase(it);
			chromo.nMut--;
			//update fitness
			w /= (1.0 - it->second.h * it->second.s);
		}
	}
}
//--------------------------------------------------------------------------
void Individuals::delet_mutation(int nmut, double kk, std::uniform_real_distribution<> findpos,
	std::gamma_distribution<> finds, std::uniform_real_distribution<> uniform, Parameters para) {
	int hom;
	double pos;
	double ss;
	double hh;

	for (int i = 0; i < nmut; i++) {
		//sample homologue
		hom = findhom(gen);
		if (para.sh_dist == 0) {
			ss = finds(gen);
			std::uniform_real_distribution<> findh(0.0, std::exp(-kk * ss));
			hh = findh(gen);
		}
		if (para.sh_dist == 1) {
			ss = uniform(gen);
			hh = uniform(gen);
		}
		if (para.sh_dist == 2) {
			double mean_sd = para.mean_hd;
			double mean_hd = para.mean_sd;
			double k = -log(2.0 * mean_hd) / mean_sd;
			gamma_distribution<double> h_cof(1.0, mean_hd);
			hh = h_cof(gen);
			uniform_real_distribution<double> s_cof(0.0, exp(-k * hh));
			ss = s_cof(gen);
		}
		if (para.sh_dist == 3) {
			ss = uniform(gen);
			hh = 0.5;
		}
		if (para.sh_dist == 4) {
			ss = uniform(gen);;
			hh = 0.2;
		}
		if (para.sh_dist == 5) {
			ss = 0.0;
			hh = 0.0;
		}
		if (para.sh_dist == 6) {
			gamma_distribution<double> ss_cof(0.25, 0.05);
			ss = ss_cof(gen);
			hh = 0.5;
		}
		if (para.sh_dist == 7) {
			gamma_distribution<double> ss_cof(0.4, 25000);
			ss = ss_cof(gen);
			ss = ss / double(0.5 * 4 * para.K);
			hh = 0.5;
			while (ss > 1.0)
			{
				ss = ss_cof(gen);
				ss = ss / double(0.5 * 4 * para.K);
			}
		}
		if (para.sh_dist == 8) {
			gamma_distribution<double> ss_cof(para.shape, (para.mean/para.shape));
			ss = ss_cof(gen);
			ss = ss / double(0.5 * 4 * para.K);
			hh = 0.5;
			if (ss > 1.0)(ss = 1.0);
			//while (ss > 1.0)
			//{
			//	ss = ss_cof(gen);
			//	ss = ss / double(0.5 * 4 * para.K);
			//}
		}
		//sample position
		pos = findpos(gen);
		//add mutation
		w *= chromo.addDelMutation(hom, pos, ss, hh);
	}
}
//--------------------------------------------------------------------------
void Individuals::delet_mutation(int nmut, std::uniform_real_distribution<> findpos,
	double ss, double hh) {
	int hom;
	double pos;

	for (int i = 0; i < nmut; i++) {
		hom = findhom(gen);
		pos = findpos(gen);
		w *= chromo.addDelMutation(hom, pos, ss, hh);
	}
}

void Individuals::SFS_sample(int g, int gg, std::ofstream* out) {

	map<double, mutation>::iterator iter;

	for (iter = chromo.mutations.begin(); iter != chromo.mutations.end(); iter++) {
		*out << g << "\t" << iter->second.s << "\t" << iter->second.h << "\t" << iter->second.pos << "\t" << gg << "\t" << iter->second.homol;
		*out << endl;
	}
}