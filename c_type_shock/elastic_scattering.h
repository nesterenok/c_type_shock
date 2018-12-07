#pragma once

class elastic_cross_section
{
public:
	// velocity in cm/s, the product of the cross section and velocity is returned (in cm3 s-1):
	virtual double get(double vel) const = 0;
};

class elastic_cross_section_const : public elastic_cross_section
{
public:
	double cs_osterbrock;
	virtual double get(double vel) const { return cs_osterbrock; }
	// polarizability in [bohr radii ^3], for h2 = 5.315, for h = 4.5, for he = 1.383, for O = 5.326 (Draine, Interstellar medium, 2010);
	elastic_cross_section_const(double polarizability, double neutral_mass, double ion_mass);
};

// The dependence of cross section on energy is taken according to Flower, MNRAS 313, L19, 2000; 
// They performed calculations of the cross-section for momentum transfer in collisions para-H2 - HCO+;
class elastic_cross_section_powerlaw : public elastic_cross_section_const
{
public:
	double vel0, g;
	double get(double vel) const;
	// At collision energies < E_lim (or vel < vel0), Osterbrock cross section is used, 
	// collision energies > E_lim, cross section is proportional to v^g; 
	// limiting collision energy in K;
	elastic_cross_section_powerlaw(double polarizability, double neutral_mass, double ion_mass, 
		double coll_energy_lim = 25., double g = -0.62);  // 24 K - more exact;
};

class elastic_cross_section_table : public elastic_cross_section
{
public:
	int nb_cs;
	double *vel_arr, *cs_arr;
	
	double get(double vel) const;
	// path to the data folder, and file name with the path within the data folder must be given:
	elastic_cross_section_table(const std::string &data_path, const std::string &name);
	~elastic_cross_section_table();
};

// Methods of this class calculate the momentum and energy exchange between two fluids;
// here the constant cross section of elastic scattering is used:
class elastic_scattering
{
public:
	int verbosity;
	double csv, mass1, mass2, reduced_mass, mass_sum;
	
	// calculation of necessary data tables based on cross section data:
	virtual void calc_data(elastic_cross_section *cross_section);
	
	// there is drift between fluids.
	// the values of momentum and energy gain are calculated and are added to the values given with the function call:
	virtual void calc_source_terms(double & mom_gain1, double & mom_gain2, double & energy_gain1, double & energy_gain2, double conc1, 
		double conc2, double vel_12_diff, double temp1, double temp2);
	
	// there is no drift between fluids (energy exchange in static gas): 
	virtual void calc_source_terms(double & energy_gain1, double & energy_gain2, double conc1, double conc2, double temp1, double temp2);
	virtual void save_data(const std::string & fname) {;}

	elastic_scattering(double mass1, double mass2, int verbosity=1);
};

// The calculation of momentum and energy exchange between neutrals and charged particles using specific cross section data;
class elastic_scatt_neutral_charged : public elastic_scattering
{
public:
	int nb_velth, nb_s;
	double *velth_arr, *s_arr, *c4_arr;
	double **c1_arr, **c2_arr, **c3_arr;

	void calc_data(elastic_cross_section *cross_section);
	
	// vel_ni_diff = vel_n - vel_ch
	void calc_source_terms(double & mom_gain_n, double & mom_gain_ch, double & energy_gain_n, double & energy_gain_ch, 
		double conc_n, double conc_ch, double vel_ni_diff, double temp_n, double temp_ch);

	void calc_source_terms(double & energy_gain_n, double & energy_gain_ch, double conc_n, double conc_ch, double temp_n, 
		double temp_ch);
	
	// parameter s is the ratio of drift velocity and "reduced" thermal speed:
	void get_integrals(double &c1, double &c2, double &c3, double vel_th, double s) const;
	void get_integrals(double &c1, double vel_th) const;
	
	void save_data(const std::string & fname);

	// the maximal and minimal "reduced" thermal speeds must be given;
	// for ion-neutrals 9.1e+3 < v < 2.9e+6 cm/s, that corresponds to reduced temperatures 1 < T < 1e+5 K (for H-ion system)
	// for electron-neutrals 3.9e+5 < v < 1.2e+8 cm/s, 1 < T < 1e+5 K
	elastic_scatt_neutral_charged(double neutral_mass, double charged_mass, double velth_min, double velth_max, int verbosity=1);
	~elastic_scatt_neutral_charged();
};

// Auxulary classes for integral calculations:
struct el_func1 
{
	double p1, p2; // p1 is the reduced thermal speed, p2 is the ratio of drift velocity and thermal speed;
	const elastic_cross_section *cs;
	// the function get() returns the product of the cross section and velocity (in cm3 s-1), velocity = x* thermal velocity;
	// the factor exp(-p2*p2) is taken into account here;
	double operator () (double x) const {		
		return x*exp(-(x-p2)*(x-p2)) *cs->get(x*p1);
	}
};

struct el_func2
{
	double p1, p2;
	const elastic_cross_section *cs;

	double operator () (double x) const {
		return x*x *exp(-(x-p2)*(x-p2)) *cs->get(x*p1);
	}
};

struct el_func3
{
	double p1, p2;
	const elastic_cross_section *cs;

	double operator () (double x) const {
		return x*x*x *exp(-(x-p2)*(x-p2)) *cs->get(x*p1);
	}
};

struct el_func4
{
	double p; // p is the reduced thermal speed;
	const elastic_cross_section *cs;

	// the function get() returns the product of the cross section and velocity (in cm3 s-1);
	double operator () (double x) const {
		return x*x *exp(-x*x)*cs->get(x*p);
	}
};
struct el_func5
{
	double p;
	const elastic_cross_section *cs;

	double operator () (double x) const {
		return x*x*x*x *exp(-x*x)*cs->get(x*p);
	}
};

struct el_func6
{
	double p;
	const elastic_cross_section *cs;

	double operator () (double x) const {
		return x*x*x*x*x*x *exp(-x*x)*cs->get(x*p);
	}
};

struct el_func7
{
	double p;
	const elastic_cross_section *cs;

	double operator () (double x) const {
		return x*x*x*x*x*x*x*x *exp(-x*x)*cs->get(x*p);
	}
};
