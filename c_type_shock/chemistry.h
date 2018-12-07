#pragma once
#include <string>
#include <vector>
#include <map>
#include "chemical_reaction_data.h"

#define NB_OF_CHEM_ELEMENTS 13 // H He C N O Si S Fe Na Mg Cl P F
#define NB_OF_CHEM_REACTION_TYPES 45

class chem_specimen;
class chem_reaction;
class chem_network;

const std::string chemical_reaction_types[] = {
	"A + CR -> Ion + e-",			// 0	CP		standard (UMIST/KIDA) format
	"A + CR -> B + Ion + e-",		// 1	CP		standard
	"A + CR -> Ion + Ion",			// 2	CP		standard
	"A + CR -> B + C",				// 3	CP		standard
	"A + CRPhoton -> Ion + e-",		// 4	  CR	standard
	"A + CRPhoton -> B + Ion + e-",	// 5	  CR	standard
	"A + CRPhoton -> B + C + D + E",// 6	  CR	standard
	"Ion + CRPhoton -> B + Ion",	// 7	  CR	standard
	"Ion + CRPhoton -> A + e-",		// 8	  CR	standard
	"A + ISPhoton -> Ion + e-",		// 9	PH		standard
	"A + ISPhoton -> B + Ion + e-",	// 10	PH		standard
	"A + ISPhoton -> B + C + D + E",// 11	PH		standard
	"Ion + ISPhoton -> B + Ion",	// 12	PH		standard
	"Ion + ISPhoton -> A + e-",		// 13	PH		standard
	"A + B -> C + D + E",			// 14			standard
	"A + B -> C + Photon",			// 15			standard
	"A + B -> Ion + e-",			// 16			standard
	"A + e- -> Ion",				// 17			standard
	"A + e- -> B + C + Ion",		// 18			standard
	"A + e- -> B + C + e-",			// 19			standard
	"A + Ion -> Ion + Photon",		// 20			standard
	"A + Ion -> B + C + D + Ion",	// 21			standard
	"A + Ion -> B + C + e-",		// 22			standard
	"Ion + e- -> A + Photon",		// 23			standard
	"Ion + e- -> A + B + C",		// 24			standard
	"Ion + Ion -> A + B + C + D",	// 25			standard
	"H + H + grain -> H2 + grain",  // 26	HFG, H2 formation on grains
	"A + GRAIN -> *A + GRAIN",		// 27	ADS, adsorption
	"*A + GRAIN -> A + GRAIN",		// 28	THD, thermal desorption
	"*A + CR -> A",					// 29	CPD
	"*A + CRPhoton -> A",			// 30	  CRD
	"*A + CRPhoton -> A + B",		// 31	  CRD
	"*A + ISPhoton -> A",			// 32	PHD
	"*A + ISPhoton -> A + B",		// 33	PHD
	"*A + C -> A + C",				// 34	  SPU, sputtering
	"*A + C -> A + B + C",			// 35	  SPU
	"*A + *B -> *C + *D + *E",		// 36	GS, grain surface reactions
	"*A + *B -> C",					// 37	GS
	"*A + *B -> C + D + E",			// 38	GS
	"*A + *B -> *A + B",			// 39	GS
	"*A + CRPhoton -> *B + *C",		// 40	CRDGS1, CR photon dissociation on grain surface via ionization and recombination
	"*A + CRPhoton -> *B + *C",		// 41	CRDGS2, via dissociation
	"*A + ISPhoton -> *B + *C",		// 42	ISDGS1, IS photon dissociation on grain surface via ionization and recombination
	"*A + ISPhoton -> *B + *C",		// 43	ISDGS2, via dissociation
	"Ion + GRAIN -> A + B + GRAIN", // 44	IAG, ion attachment on grains
};

std::string get_reaction_type(const int i);
std::string get_reaction_name(std::string rcode, const std::vector<chem_specimen> &, const int *reactant, int nb_of_reactants, 
	const int *product, int nb_of_products);

// The function rearrange the specimen order in the array: 1 - adsorbed, 2 - neutrals, 3 - ions, 4 - electrons;
void check_order(const std::vector<chem_specimen> &, int *arr, int nb);

// Determination of the reaction type;
int react_type_determ(std::string rcode, const std::vector<chem_specimen> &, const int *reactant, int nb_of_reactants, 
	const int *product, int nb_of_products, int verbosity = 1);

// Calculation of the reaction rate coefficient;
// 1. temperature is measured in K;
// 2. cr_ioniz_rate is in fractions of STANDARD_CR_IONIZ_RATE;
// 3. visual_extinct is the dust extinction at visible wavelengths, at 5500 A = 0.55 um;
// 4. uv_field_strength is in fractions of the standard UV flux, photon energy > 5 eV (Draine, 1978);
// 5. rates are not extrapolated to the temperatures outside the approximation range;
// 6. ads_dust_area = pi *a^2 *n_g, ads_grain_veln2 = (vel_n - vel_g)^2
// 7. unimolecular - rate in s-1, bimolecular - rate in cm3 s-1
double reaction_rate(const std::vector<chem_specimen> & species, const accretion_rate_functions *accr_func, const chem_reaction & reaction, 
	double cr_ioniz_rate, double uv_field_strength, double visual_extinct, double vel_n, double vel_i, 
	double ads_dust_area, double ads_grain_veln2, double temp_n, double temp_i, double temp_e, double temp_d, 
	double photoreact_factor_cr, double desorption_factor_cr, double photodes_factor_cr, double photodes_factor_is, double conc_h_tot, double nb_ads_sites_grain);

// Temperature is in erg;
void chemistry_source_terms(double & mom_gain_n, double & mom_gain_i, double & mom_gain_e, 
	double & energy_gain_n, double & energy_gain_i, double & energy_gain_e, double & energy_gain_d, 
	const std::vector<chem_specimen> &, const chem_reaction &, double vel_n, double vel_i, 
	double ads_grain_velz, double ads_grain_veln2, double temp_n, double temp_i, double temp_e, double temp_d, double rate);

// Matar et al., J. Chem. Phys. 133, 104507, 2010; temperature must be in K;
double get_h_sticking_coeff(double gas_temp, double dust_temp);
double get_h2_sticking_coeff(double gas_temp, double dust_temp);

// Approximation formula of Matar et al. (2010) is used for all species, 
// the critical temperature is set to be proportional to the mass of the particle, mass in g;
double get_sticking_coeff(double gas_temp, double mass);

// He et al., ApJ 823, p.56 (2016); binding energy, dust temperature in K;
double get_sticking_coeff_he2016(double bind_en, double dust_temp, double beta = 0.11, double gamma = 0.042);

// Initialization of initial chemical abundance, file name must contain the full path to the file;
void init_chem_abund(const std::string file_name, const chem_network *, double *abundances);

class chem_element
{
public:
	double mass;
	std::string name;
};

class chem_specimen
{
public:
	bool is_enth_def; // true - enthalpy is defined;
	int charge, n_nb; // n_nb - neutral nb for ions and adsorbed species;
	int formula[NB_OF_CHEM_ELEMENTS];
	// mass in g, bind_en and diff_barrier in K; 
	// parameters bind_en, surf_freq and diff_barrier are for adsorbed species;
	double enthalpy, mass, bind_en, surf_freq, diff_barrier; // masses are considered to be equal if dm < 0.01 a.m.u.
	// type: neutral, ion, e-, adsorbed 
	std::string name, type;

	// the neutrals are the "smallest", then - positive ions, negative ions (and electron), adsorbed species; 
	// within the group, species are sorted by mass, if masses are equal - by name;
	friend bool operator == (const chem_specimen &, const chem_specimen &);
	friend bool operator != (const chem_specimen &, const chem_specimen &);
	friend bool operator < (const chem_specimen &, const chem_specimen &);
	friend bool operator > (const chem_specimen &, const chem_specimen &);

	chem_specimen();
};

class chem_reaction
{
public:
	int type, nb_of_reactants, nb_of_products, nb_of_fits, nb_of_param, cs_reconstr;
	int	*reactant, *product;
	
	// Note: for unimolecular reactions only mass1 is defined;
	double energy_released, mass1, mass2, reduced_mass, mass_sum, max_rate, min_rate;
	double *temp_min, *temp_max, *parameters;
	
	std::string name, rcode;
	const reaction_rate_data *rate_data;
	
	// function deletes data with exception of rate_data:
	void delete_data();
	
	chem_reaction & operator=(const chem_reaction &);
	friend bool operator == (const chem_reaction &, const chem_reaction &);
	chem_reaction();
	chem_reaction(const chem_reaction &);
	~chem_reaction() { delete_data(); }
};

class chem_network
{
public:
	// max_nb_carbon - maximal value of carbon atoms allowed in the chemical specimen;
	int verbosity, nb_of_species, nb_of_gmantle_species, nb_of_reactions, nb_reactions_ion_grains, max_nb_carbon, 
		h2_nb, ah2_nb, h_nb, he_nb, e_nb, oi_nb, ci_nb, cii_nb, co_nb, h2o_nb, oh_nb, nh3_nb, ch3oh_nb, h2_h_diss_nb;

	std::map<int, chem_element> elements;
	std::vector<chem_specimen>	species;
	std::vector<chem_reaction>	reaction_array;	
	std::vector<reaction_rate_data> rate_data_array;
	
	void add_element(const std::string &);
	void set_max_nb_carbon(int i) { max_nb_carbon = i; } 
	
	// File names with the full path must be given:
	// the function reads file with gas-phase species:
	void init_gas_phase_species(const std::string file_name);

	// the function must be called after gas-phase species being initialized:
	// in the file: specimen name, bunding energy in K, comments;
	void init_gmantle_species(const std::string file_name);
	
	// the function excludes chemical species from the list, must be called before initialization of chemical reactions
	// and before initialization of specimen numbers (next routine):
	void exclude_chem_specimen(std::string sname);
	// the function must be called after all species being initialized (must be called in all cases):
	void init_species_nbs();
	
	// the full file name with the chemistry network must be given, data in the format of UMIST database:
	void init_network_umistf(const std::string fname, bool update = true);
	
	// this routine must be called after gas-phase reactions being initialized,
	void init_photoreact_surface_chemistry();
	// Adding neutral-neutral reactions on grain surfaces,
	// the full file name with data on grain surface reactions must be given (NAUTILUS format):
	void init_grain_surface_chemistry(const std::string fname);

	// H2 formation on grains is treated individually,
	void init_h2_formation();  // for "ad-hoc" modelling and empirical formula, see parameters.h
	
	// the function must be called after all reactions being initialized (must be called in all cases):
	void check_reactions();
	void print_network(const std::string &output_path);
	
	int find_specimen(const std::string &name) const;
	int find_reaction(const std::string &name) const;
	
	chem_network(const std::string &path, int verbosity =1);
	~chem_network();
};

// path to the data
void reformat_chemical_data_Belloche2014(const std::string &path);
