#pragma once

#include <string>
#include "spectroscopy.h"
#include "coll_rates.h"

#define H2_COLL_CUBIC_SPLINE 1
#define H2_DISS_CUBIC_SPLINE 0

// 0 - data by Wrathmall et al. (2007), Martin & Mandy (1995), see below for full references;
// 1 - data by Lique (2015) for lowest 54 levels, data by Wrathmall et al. (2007), Martin & Mandy (1995);
// 2 - data by Lique (2015) and by Bossion et al. (2018)
#define H2_H_COLL_DATA 2

// 0 - data by Flower & Roueff (1998, 1999);
// 1 - data by Wan et al. (2018);
#define H2_H2_COLL_DATA 1

// the function reads the data for para-H2 and ortho-H2 rate coefficients and combine them in one table for h2 molecule;
void h2_coll_data_process(const std::string &data_path);
// rewrite data given by Dr. Bossion in more compact form; 
void h2_coll_data_process_bossion(const std::string &data_path);

// H2 - ortho-H2 rate coefficients. The data from Flower & Roueff, J. Phys. B 32, 3399 (1999);
// original data at 100-6000 K.
// Additional point at 0 K was added, the data were set to zero at this temperature (the same is true for other data);
// Note, for the cold cloud simulations the low temperature data are needed;
class h2_oh2_flower_data
#if H2_COLL_CUBIC_SPLINE
	: public collision_data_cub_spline
#else
	: public collision_data
#endif
{
public:
	h2_oh2_flower_data(const std::string &path, const energy_diagram *, int verbosity=1);
};

// H2 - para-H2 rate coefficients. The data from Flower & Roueff, J. Phys. B 31, pp. 2935-2947 (1998);
class h2_ph2_flower_data
#if H2_COLL_CUBIC_SPLINE
	: public collision_data_cub_spline
#else
	: public collision_data
#endif
{
public:
	h2_ph2_flower_data(const std::string &path, const energy_diagram *, int verbosity=1);
};

// H2 - He rate coefficients. The data from Flower et al., J. Phys. B 31, 1105-1113 (1998);
class h2_he_flower_data
#if H2_COLL_CUBIC_SPLINE
	: public collision_data_cub_spline
#else
	: public collision_data
#endif
{
public:
	h2_he_flower_data(const std::string &path, const energy_diagram *, int verbosity=1);
};

// The data from Wrathmall et al., MNRAS 382, pp. 133–138 (2007); Wrathmall & Flower, J. Phys. B 40, pp. 3221-3230 (2007);
// http://ccp7.dur.ac.uk/pubs.html
// the contribution of reactive channels to the collisions are added (including those that change ortho-para state), 
// see Le Bourlot et al., MNRAS 305, pp. 802-810 (1999)
class h2_h_wrathmall_data
#if H2_COLL_CUBIC_SPLINE
	: public collision_data_cub_spline
#else
	: public collision_data
#endif
{
public:
	h2_h_wrathmall_data(const std::string &path, const energy_diagram *, bool reactive_channels, int verbosity=1);
};

// The data by Wan et al., ApJ 862, p.132 (2018)
class h2_h2_wan_data
#if H2_COLL_CUBIC_SPLINE
	: public collision_data_cub_spline
#else
	: public collision_data
#endif
{
public:
	h2_h2_wan_data(const std::string &path, const energy_diagram *, bool coll_partner_is_ortho, int verbosity=1);
};


// The data by Lique F., MNRAS 453, 810-818 (2015); 
// 54 levels of ortho- and para-H2, 100 <= T <= 5000 K;
class h2_h_lique_data
#if H2_COLL_CUBIC_SPLINE
	: public collision_data_cub_spline
#else
	: public collision_data
#endif
{
public:
	h2_h_lique_data(const std::string &path, const energy_diagram *, int verbosity=1);
};

// The data by Bossion et al. MNRAS 480, p.3718, 2018;
// 1000 < T <= 5000 K; all levels, but not ordered;
class h2_h_bossion_data
#if H2_COLL_CUBIC_SPLINE
	: public collision_data_cub_spline
#else
	: public collision_data
#endif
{
public:
	h2_h_bossion_data(const std::string &path, const energy_diagram *, int verbosity=1);
};

// The data by Martin & Mandy, ApJ 455, pp. L89-L92 (1995); provided by e-mail by David Flower;
// data are based on quasiclassical calculations, is valid only at high temperatures >> 1000 K;
// linear interpolation is used for the rate coefficients; available only for the full H2 level set (ortho- and para-);
// Note: the levels v=11 j=13 e=36310.61 cm-1; v=11 j=14 e=36593.4 cm-1 are absent in the list of Martin & Mandy (1995) data;
class h2_h_martin_data : public collision_data
{
public:
	h2_h_martin_data(const std::string &path, const energy_diagram *, int verbosity=1);
};

// The rate coefficients are compilated based on various sources:
// Gerjuoy E., Stein S., Physical Review 97, p. 1671 (1955); 
// England et al., Aust. J. Phys. 41, pp. 573-86 (1988);
// Ehrhardt et al., Physical Review 173, p. 222 (1968); 
// Yoon et al., J. Phys. Chem. Ref. Data 37, p. 913 (2008);
class h2_e_data: public collision_data
{
public:
	h2_e_data(const std::string &path, const energy_diagram *, int verbosity=1);
};

// Gonzalez-Lezana & Honvault, MNRAS 467, 1294 (2017),
// there is problem with the data from fig.7, (v=1,j=0)->(v=2,j=1,2,3), the high values at low temperatures
// the spline is not used here, the temperature grid is dense,
class h2_hp_gonzalez_lezana_data : public collision_data
{
public:
    h2_hp_gonzalez_lezana_data(const std::string &path, const energy_diagram *, int verbosity = 1);
};

// The class that calculates collisional rates;
class h2_collisions : public collisional_transitions
{
public:
	void get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
		double temp_neutrals, const double *concentration, const int *indices) const;

    void get_rate_ions(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
        double temp_neutrals, double temp_ions, const double *concentration, const int *indices) const;
	
	void set_gas_param(double temp_neutrals, double el_temp, double he_conc, double ph2_conc, double oh2_conc, double h_conc, 
        double el_conc, double *&concentration, int *&indices) const;
	
    void set_ion_param(double temp_neutrals, double temp_ions, double hp_conc, double h3p_conc, 
        double *&concentration, int *&indices) const;

	void check_spline(int ilev, int flev, const std::string & fname) const;
	h2_collisions(const std::string &, const energy_diagram *, int verbosity =1);
};


// H2 dissociation data
class h2_h_dissociation_data
#if H2_DISS_CUBIC_SPLINE
	: public dissociation_data_cub_spline
#else
	: public dissociation_data
#endif
{
public:
	h2_h_dissociation_data(const std::string & data_path, const energy_diagram *, int verbosity);
};


// H2 excitation processes;
class h2_excitation_process
{
public:
	virtual double get_efficiency(int ilev, int flev) const { return 0.; }
	virtual double get_efficiency(int lev) const { return 0.; }
	virtual double get_efficiency(int lev, double gtemp) const { return 0.; }
};

// H2 excitation due to formation on grains;
class h2_grain_formation : h2_excitation_process
{
protected:
	int nb_lev;
	double average_energy;
	double *popul;

public:
	// the function returns the fraction of H2 molecules at level l, produced on grains:
	double get_efficiency(int l) const { return popul[l]; }

	h2_grain_formation(const energy_diagram *);
	~h2_grain_formation();
};

// H2 excitation due to formation in the gas phase;
class h2_gasphase_formation : h2_excitation_process
{
protected:
	int nb_lev, nb_temp;
	double step, step_inv, max_temp;
	double **popul;

public:
	// the function returns the fraction of H2 molecules at level l, the temperature of the distribution is given:
	double get_efficiency(int l, double gtemp) const;

	h2_gasphase_formation(const energy_diagram *);
	~h2_gasphase_formation();
};

// Is not used in the code.
// H2 excitation due to cosmic rays; the returned values must be multiplied by the cosmic ray or x-ray ionization rates in s-1; 
// Note: the routine does not work yet for high (> 0.01) fractional ionizations;
// at small fractional ionizations < 1.e-6 the data at 1.e-6 are used;
class h2_excit_cosmic_rays : h2_excitation_process
{
protected:
	int nb_ion, nb_fin_lev, nb_init_lev;
	double *ion_grid;
	double **eff, **exit_eff;
	
public:
	void set_ionization(double ionization, int &index, double &param) const;
	int get_nb_lev() const { return nb_init_lev; }
	
	double get_efficiency(int il, int fl, int index, double param) const;
	double get_efficiency(int i, int index, double param) const;

	h2_excit_cosmic_rays(const std::string &, const energy_diagram *, int verbosity = 1);
	~h2_excit_cosmic_rays();
};
