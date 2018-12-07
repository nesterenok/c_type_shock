#pragma once
#include "spectroscopy.h"
#include "coll_rates.h"

// The routine merges the data of CO-H collision rates from different files;
void merge_co_h_coll_data(const std::string & path);

// data on CO-H2 collisions; 41 rotational levels of ground vibrational state of CO are considered, 2 < T < 3000 K;
// Yang et al., ApJ 718, pp. 1062–1069, 2010;
class co_h2_coll_data : public collision_data
{
public:
	co_h2_coll_data(const std::string path, const energy_diagram *, bool is_ortho_h2, int verbosity=1);
};

// data on CO-He collisions, 15 lowest rotational levels of CO are considered, 5 < T < 500 K;
// Cecchi-Pestellini et al., ApJ 571, pp. 1015–1020, 2002;
class co_he_coll_data : public collision_data
{
public:
	co_he_coll_data(const std::string path, const energy_diagram *, int verbosity=1);
};

// data on CO-H collisions, rotational transitions of the ground state; 2 < T < 3000 K; 
// Walker et al., ApJ 811, p. 27 (2015);
class co_h_coll_data : public collision_data
{
public:
	co_h_coll_data(const std::string path, const energy_diagram *, int verbosity=1);
};

// data on CO-H collisions, ro-vibrational transitions, v->v' = 0->0, 1->1, 2->2, 1->0, 2->1; 10 < T < 3000 K; 
// pure rotational transitions - Walker et al., ApJ 811, p. 27 (2015); rovibrational - Song et al., J. Chem. Phys. 142, 204303 (2015); Song et al., ApJ 813, p. 96, 2015; 
class co_h_vibr_coll_data : public collision_data
{
public:
	co_h_vibr_coll_data(const std::string path, const energy_diagram *, int verbosity=1);
};

// The methods of this class calculate the total collisional rates for transitions between levels of CO;
class co_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
		double el_conc, double *&concentration, int *&indices) const;

	void get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
		double temp_neutrals, const double *concentration, const int *indices) const;

	co_collisions(const std::string &, const energy_diagram *, int verbosity =1);
};
