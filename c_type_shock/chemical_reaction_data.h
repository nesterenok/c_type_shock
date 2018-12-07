#pragma once
#include <string>
#include <cstring>
#include <cmath>

// Auxulary routine to construct data files with gas-grain reactions: adsorption, desorption, sputtering;
// file with the list of adsorbed species is required, in UMIST format, file name with the full path:
void construct_gas_grain_reactions(std::string input_fname, std::string output_path);

// Auxulary routine to construct data files with ion recombination reactions on grains; 
void construct_ion_recomb_grains(std::string output_path);

// H2 and H ionization by ions. Is not used;
void construct_hydrogen_ioniz_reactions(std::string output_path);

// Pure virtual class:
class process_cross_section
{
public:
	// for the neutral-ion reactions, the first mass is the neutral mass (g), energy is in erg;
	double mass1, mass2, reduced_mass, mass_sum, en_min;
	
	// the function returns the product of the cross section and collisional energy (cm2 erg):
	virtual double get(double energy) const = 0;
	
	process_cross_section(double mass); // the mass of second partner is assumed to be extremely large in this case;
	process_cross_section(double mass1, double mass2);
	virtual ~process_cross_section() {;}
};

// Cross section data is stored in the table; see code for the file data format;
// is not used yet;
class chem_reaction_cross_section_table : public process_cross_section
{
public:
	int nb_cs;
	double *en_arr, *cs_arr;
	
	// energy is in erg, returned value in cm2*erg;
	double get(double energy) const;
	chem_reaction_cross_section_table(const std::string &data_path, const std::string &name, double m1, double m2);
	~chem_reaction_cross_section_table();
};

// Reconstruction of the cross section by the method by Draine, Katz, ApJ 306, p.655, 1986;
// is not used yet;
class chem_reaction_cross_section_1 : public process_cross_section
{
private:
	double p, sigma;

public:
	// energy is in erg, returned value in cm2*erg;
	double get(double energy) const;
	// parameters a, b, c are reaction rate parameters, k = a*T^b*exp(-c/T)
	chem_reaction_cross_section_1(double mass1, double mass2, double a, double b, double c);
};

// Draine, Roberge, Dalgarno, ApJ 264, p. 485, 1983;
// cross section approximation for ion-neutral collisions, A + I+ -> A+ + e- + I+
// is not used yet;
class chem_reaction_cross_section_2 : public process_cross_section
{
private:
	double sigma;

public:
	// energy is in erg, returned value in cm2*erg;
	double get(double energy) const;
	// ionization potential (erg) and oscillator strength of the neutral specimen;
	// oscillator strength 0.665 for H; 1 for H2 and He;
	chem_reaction_cross_section_2(double mass1, double mass2, double ioniz_potential, double oscill_strength);
};

// Approximation of sputtering yield by Draine, Salpeter, ApJ 231, 77, 1979;
class sputtering_yield : public process_cross_section
{
protected:
	double norm, bind_en, p; // binding energy is stored in erg;

public:
	// energy is in erg, returned value atom/molecule per projectile multiplied by energy (in erg), 
	// here, "cross section" is sputtering yield (atom per atom);
	double get(double energy) const;
	
	// given mass must be in g, binding energy given must be in K:
	sputtering_yield(double m_proj, double m_target, double binding_energy);
};

class reaction_rate_data
{
private:
	int nb_v1, nb_v2;
	double int_err;
	double *var1, *var2; // 1 - reduced temperature [K], 2 - drift velocity divided by "reduced" thermal velocity;
	double **f_arr;

public:
	// returned value is <sigma*velocity>, for sputtering - "cross section" is sputtering yield (atom per atom);
	// the first parameter is "reduced" temperature (K), second - drift velocity divided by "reduced" thermal velocity:
	double get(double reduced_temp, double norm_drift_vel) const;
	// the data table is evaluated:
	void calc_data(process_cross_section *cross_section);
	void print(std::string rname);
	void delete_data();

	reaction_rate_data();
	reaction_rate_data(const reaction_rate_data &);
	~reaction_rate_data() { delete_data(); }
};

struct ch_func1
{
	// p1 is the "reduced" temperature in erg; p2 is the ratio of drift velocity and "reduced" thermal speed;
	double p1, p2;
	const process_cross_section *cs;
	// the function cs->get() returns the product of the cross section and collisional energy (cm2*erg for chemical reactions);
	// x*x = E/kT;
	double operator () (double x) const {
		return exp(-(x-p2)*(x-p2))*cs->get(x*x*p1);
	}
};

struct ch_func2
{
	// p1 is the "reduced" temperature in erg;
	double p1;
	const process_cross_section *cs;

	double operator () (double x) const {
		return exp(-x)*cs->get(x*p1);
	}
};

// The class contains the data arrays with accretion rate and energy removal rate values. 
// the main equations are given by Draine, Sutin, ApJ 320, p. 803 (1987);
// it is assumed that ion charge is equal to -+1; reduced temperature must be in cm*K;
class accretion_rate_functions
{
protected:
	int nb_rtemp, nb_chr;
	double constant;
	double *rtemp_arr, *chr_arr; // charge ratio = grain charge /ion charge 
	double **accr_rate_func, **energy_removal_func;	// first variable - charge ratio, second - reduced temperature;
	
	double calc_accr_rate(double rtemp, double charge_ratio);
	double calc_energy_removal(double rtemp, double charge_ratio);

public:
	// reduced temperature = temperature *grain radius, must be in cm*K;
	// the function returns "reduced" accretion rate coefficient, int_0^\infty dx x exp(-x) sigma(x*t, nu)/(pi*a*a)
	double get_accr_rate(double rtemp, double charge_ratio) const;
	// "reduced" rate of energy removal from gas through particle sticking on grains:
	double get_energy_removal(double rtemp, double charge_ratio) const;

	void save(const std::string &path) const;

	accretion_rate_functions();
	~accretion_rate_functions();
};
