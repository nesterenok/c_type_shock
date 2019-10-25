#pragma once
#include "spectroscopy.h"
#include "coll_rates.h"

// data on OH-H2 collisions; 20 rotational levels of ground vibrational state of OH are considered, 10 < T < 150 K;
// Klos et al., MNRAS 471, 4249 (2017);
// there are two sets of collisional data: for H2 in J = 0, and for H2 in J >= 1;
class oh_h2_coll_data : public collision_data
{
public:
	oh_h2_coll_data(const std::string path, const energy_diagram *, bool is_h2_j0, int verbosity=1);
};

// data on OH-He collisions, 46 lowest rotational levels of ground vibrational state of OH are considered, 5 < T < 500 K; 
// there is a discrepancy betweeb HITRAN 2016 and adopted level list for these data for 45,46 levels,
// Klos et al., Chemical Physics Letters 445, 12 (2007);
class oh_he_coll_data : public collision_data
{
public:
	oh_he_coll_data(const std::string path, const energy_diagram *, int verbosity=1);
};

// Offer et al., J. of Chemical Physics 100, 362 (1994);
// the data in the files are in simple format (only rate coefficients), original data are in LAMBDA format,
class oh_hf_h2_coll_data : public collision_data
{
public:
    oh_hf_h2_coll_data(const std::string path, const energy_diagram*, bool coll_partner_is_ortho, int verbosity = 1);
};


// Methods of this class calculate the total collisional rates for transitions between levels of OH;
class oh_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double h2j0_conc, double h2j1_conc, double h_conc,
		double el_conc, double *&concentration, int *&indices) const;

	void get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
		double temp_neutrals, const double *concentration, const int *indices) const;

	oh_collisions(const std::string &, const energy_diagram *, int verbosity =1);
};

class oh_hf_collisions : public collisional_transitions
{
public:
    void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
        double el_conc, double*& concentration, int*& indices) const;

    void get_rate_neutrals(const energy_level& up_lev, const energy_level& low_lev, double& down_rate, double& up_rate,
        double temp_neutrals, const double* concentration, const int* indices) const;

    oh_hf_collisions(const std::string&, const energy_diagram*, int verbosity = 1);
};
