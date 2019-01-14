#pragma once

#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <string>
#include <vector>

#include "utils.h"
#include "chemistry.h"
#include "elastic_scattering.h"
#include "spectroscopy.h"
#include "coll_rates.h"
#include "dust_model.h"
#include "radiation_field.h"
#include "lvg_method_functions.h"

// temperatures: neutral, ion, electron; velocities: neutral and ion;
#define NB_MHD_EQUATIONS 5

// two type of dust: large grains and pah molecules; 
// for pah molecules charge distribution is calculated (not average charge);
class evolution_data
{
protected:
	int verbosity;
	int nb_of_equat, nb_of_species, nb_of_gas_species, nb_of_gmantle_species, nb_of_dust_comp, nb_of_grain_charges, nb_dct, nb_mhd, nb_rtd;
	int	nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_onh3, nb_lev_pnh3, nb_lev_ch3oh, nb_lev_oi, nb_lev_ci, nb_lev_cii;
	
	// CR ionization rate is in units of standard value;
	// UV field strength is in units of the "standard" UV flux given by Draine (1978);
	double cr_ioniz_rate, ir_field_strength, uv_field_strength, visual_extinct;
	// heating and cooling variables have negative values for cooling process by definition;
	double neut_heat_chem, neut_heat_dust_coll, pheff_gas_heat, neut_cr_heat, neut_heat_h2, neut_heat_atoms, neut_heat_ph2o, 
		neut_heat_oh2o, neut_heat_co, neut_heat_oh, neut_heat_pnh3, neut_heat_onh3, neut_heat_ch3oh_a, neut_heat_ch3oh_e,
		neut_heat_scatt_ions, neut_heat_scatt_el, el_heat_atoms, el_heat_h2, el_heat_ph2o, el_heat_oh2o, el_heat_scatt_neut, 
		el_heat_scatt_ions, el_heat_chem, ion_heat_scatt_n, ion_heat_scatt_el, ion_heat_chem, rad_energy_loss_h2;
	double vel_n, vel_i, vel_ni_diff, vel_turb, vel_n_grad, vel_i_grad, vel_grad_min;
	double temp_n, temp_i, temp_e, temp_d, temp_n_erg, temp_i_erg, temp_e_erg, temp_d_erg;
	double conc_n, conc_i, conc_e, conc_he, conc_h2, conc_ph2, conc_oh2, conc_h2j0, conc_h, conc_h_tot, conc_h2o, conc_co, 
		conc_oh, conc_nh3, conc_ch3oh, conc_oi, conc_ci, conc_cii, conc_ice;
	double h2_prod, h2_prod_gr, h2_prod_gas, h2o_prod, co_prod, oh_prod, nh3_prod, ch3oh_prod, ci_prod, cii_prod, oi_prod;
	double nb_gain_n, nb_gain_i, nb_gain_e, mass_gain_n, mass_gain_i, mass_gain_e, mom_gain_n, mom_gain_i, mom_gain_e, 
		energy_gain_n, energy_gain_i, energy_gain_e, energy_gain_d;
	double ads_dust_area, ads_grain_velz, ads_grain_veln2;
	double coverage, ion_mass, magn_field_0, magn_field, shock_vel, photodes_factor_cr, photodes_factor_is, desorption_factor_cr, 
		photoem_factor_is_uv, photoem_factor_is_vis;
	double oh2_form_gaschem, oh2_form_grains, oh2_form_hcoll, h2_h_diss_rate;
	
	bool *is_dust_charge_large, *is_charged_dust_as_ions;
	int *indices, *min_grain_charge, *max_grain_charge, *nb_dch;
	// two arrays to store dust heating by molecular emission, one for H2 emission, second - for other molecules;
	double *dh_isrf_arr, *coll_partn_conc, *alpha, *beta, *wt2_arr, *grain_velz, *grain_veln2, *grain_veli2, 
		*grain_velz_g, *av_grain_velz, *grain_conc, *dust_heat_mline, *dust_heat_h2_line, *dust_heat_coll, *dust_heat_chem, 
		*chem_reaction_rates, *chem_heating_rates_n, *phel_rate_cr, *phel_rate_uv, *phel_rate_vis, *ion_neutr_rate, *el_att_rate;
	
	const energy_diagram *OI_di, *CI_di, *CII_di, *h2_di, *ph2o_di, *oh2o_di, *co_di, *oh_di, *onh3_di, *pnh3_di, 
		*ch3oh_a_di, *ch3oh_e_di;
	const einstein_coeff *OI_einst, *CI_einst, *CII_einst, *h2_einst, *ph2o_einst, *oh2o_einst, *co_einst, *oh_einst, 
		*pnh3_einst, *onh3_einst, *ch3oh_a_einst, *ch3oh_e_einst;
	const collisional_transitions *OI_coll, *CI_coll, *CII_coll, *h2_coll, *ph2o_coll, *oh2o_coll, *co_coll, *oh_coll, 
		*onh3_coll, *pnh3_coll, *ch3oh_a_coll, *ch3oh_e_coll;
	
	const h2_h_dissociation_data *h2_h_diss_data;

	const h2_excit_cosmic_rays *h2_excit_cr;
	const h2_grain_formation *h2_excit_gf;
	const h2_gasphase_formation *h2_excit_gasph;
	const accretion_rate_functions *accr_func;
	
	chem_network *network;
	elastic_scattering *elastic_h2_ions, *elastic_h_ions, *elastic_he_ions, *elastic_h2_el, *elastic_he_el, *elastic_h_el;
	
	const dust_model *dust;
	const lvg_method_data *loss_func_line_phot, *loss_func_cont_phot;
	dust_heating_ISRF *dheat_isrf;	
	
	// arrays to save radiative transfer factors:
	// be carefull about dimension of dheat_efficiency,
	double *gamma_factors, *delta_factors, *dheat_efficiency, *esc_prob_int1, *esc_prob_int2;

public:	
	void set_parameters(double vis_ext, double cr_ion, double uv_field, double ir_field);
	
	void set_veln_grad(double);
	void set_veli_grad(double);
	void set_vel_turb(double vt) { vel_turb = vt; }
	void set_grain_charge_ranges(bool *is_av_ch, int *min_gch, int *max_gch, int nb_of_chs);
	void set_magn_field(double m) { magn_field = magn_field_0 = m; } // in Gauss

    void set_tolerances(N_Vector abs_tol);

	int get_nb_of_species() const { return nb_of_species; };
	void get_nb_of_levels(int & nb_lev_h2, int & nb_lev_h2o, int & nb_lev_co, int & nb_lev_oh, int & nb_lev_pnh3, int & nb_lev_onh3,
		int & nb_lev_ch3oh, int & nb_lev_ci, int & nb_lev_oi, int & nb_lev_cii) const;
	void get_nbs(int & nb_of_grain_charges, int & nb_of_equat, int & nb_dct, int & nb_mhd) const;

	// the data for grains of i component are located in variable list having numbers nb1 <= i <= nb2;
	void get_dust_component_nbs(int i, int & nb1, int & nb2) const;
	int get_dust_zmin(int i) const { return min_grain_charge[i]; }
	int get_dust_zmax(int i) const { return max_grain_charge[i]; }

	double get_veln_grad() const { return vel_n_grad; }
	double get_veli_grad() const { return vel_i_grad; }
	double get_vel_grad_min() const { return vel_grad_min; }
	double get_vel_turb() const { return vel_turb; }
	double get_shock_speed() const { return shock_vel; }

	// the average temperature of grains that are able to adsorb species (averaging with surface area as a weight):
	double get_av_dust_temp() const { return temp_d; }
	
	const chem_network *get_network() const { return network; }
	const dust_model *get_dust() const { return dust; }
	
	void get_neutral_heating(double & atomic_n, double & h2_n, double & h2o_n, double & co_n, double & oh_n, double & nh3_n, 
		double & ch3oh_n, double & coll_h, double & chem_h, double & pheff_h, double & cr, double & scatt_i, double & scatt_e, 
		double & rad_en_loss_h2) const;	
	void get_electron_heating(double & atomic_e, double & h2_e, double & h2o_e, double & scatt_n, double & scatt_i, double & chem) const;
	void get_ion_heating(double & scatt_n, double & scatt_e, double & chem) const;

	// dust heating normalized on one dust grain [erg s-1],
	// by interstellar radiation field, surface chemistry, gas-dust collisions, H2 line emission, by all molecules (except H2);
	void get_dust_heating_rates(int i, double & isrf, double & chem, double & coll, double & h2_mol, double & mol) const;
	// grain charging rates [cm-3 s-1]:
	void get_dust_charging_rates(int i, double & el_att, double & ion_neutr, double & phel_uv, double & phel_vis, double & phel_cr) const;
	
	int get_reaction_nb() const { return network->nb_of_reactions; }
	double get_reaction_rate(int i) const;
	double get_h2_form_grains() const { return h2_prod_gr; }

	// the pointer to the data on heating rates of neutral fluid, number of reactions are returned:
	void get_chem_heating_rates(double *& chem_heat_rates, int & nb) const;

	// parameters of H2 molecule chemistry and ortho-para conversion;
	void get_h2_chem(double & h2_gr, double & h2_gas, double & o_grains, double & o_gchem, double & o_hcoll, double & h2_h_diss) const;
	
	// vector defining the ODE system for the chemistry evolution of the static gas.
	virtual int f(realtype t, N_Vector y, N_Vector ydot);

	// Note: some functions of the classes use internal variables;
	// the heating rate variable > 0 for heating and < 0 for cooling;
	// the dimension of heating efficiency dh_eff[] is cm-1 cm-3 s-1, dimension of grain heating rate heating_d[] is cm-1 s-1
    void specimen_population_derivative(const realtype *y_data, realtype *ydot_data, int nb, const energy_diagram *, const einstein_coeff *,
        const collisional_transitions *, double *coll_partn_conc, int *indices, double vel_grad, double & heating_n,
        double & heating_e, double & rad_energy_loss, double *heating_d);
	
    void calc_radiative_coeff(const realtype *y_data, int nb, const energy_diagram *, const einstein_coeff *, 
        const collisional_transitions *, double *coll_partn_conc, int *indices, double vel_grad, double *g_factors, 
        double *d_factors, double *dh_eff, double *ep_int1, double *ep_int2, int nb2);

	double calc_ice_conc(const N_Vector &y) const; // calculation of adsorbed atom/molecule concentration in cm-3
	double calc_conc_ph2(const N_Vector &y) const; // calculation of p-H2 concentration, cm-3
	double calc_conc_h_tot(const N_Vector &y) const;
	double calc_hydrocarbon_conc(const N_Vector &y) const;

	// ion mass density include contribution from PAHs and small grains here:
	void calc_ion_dens(const N_Vector &y, double & ion_conc, double & ion_pah_conc, double & ion_mass_dens) const;
	void calc_neutral_dens(const N_Vector &y, double & neut_conc, double & neut_mass_dens) const;

	// grain charge (in module of electron charge) and grain area (cm2) must be given, 
	// grain velocity in the laboratory frame is calculated (in z direction), cm/s:
	double calc_grain_velocity(const N_Vector &y, double charge, double area, double & a, double & b, double & wt2, double & velx) const;

	// calculation of the parameter wt2/(1+wt^2), i is the number of dust component:
	double calc_dust_ion_coupling(const N_Vector &y, int i) const;

	// calculation of the average velocity of the grains of dust component i (in z direction), cm/s:
	double calc_av_grain_velocity(const N_Vector &y, int i) const;

	// calculation of the average charge of one grain of dust component i:
	double calc_average_grain_charge(const N_Vector &y, int i) const;
	
	// calculation of grain concentration of dust component i, cm-3:
	double calc_grain_conc(const N_Vector &y, int i) const;

	// calculation of the total charge of the gas (without grains), cm-3:
	double calc_total_gas_charge(const N_Vector &y) const;

	// calculation of the total charge of (all) grains, cm-3:
	double calc_total_grain_charge(const N_Vector &y) const;
	
	// recalculation of ranges for grain charge distribution:
	bool recalc_grain_charge_ranges(N_Vector y, std::vector<double> & new_y);
	void recalc_grain_charge_ranges(double temp_e, double conc_e, double conc_h_tot, std::vector<double> & init_ch);

	void create_file_radiative_transfer(const std::string & output_path) const;
	void save_radiative_transfer_data(const std::string & output_path, double ty) const; // ty is evolution time in yrs;

	// reinitialization of arrays that have dimension equal to the sum of sizes of dust component charge distributions;
	void reinit_arrays(int nb_of_grain_charges);

	// path - path to the data files;
	// maximal values for para-H2O - 413, ortho-H2O - 411, H2 - 298 (restriction by the data on the Einstein coeff.):
	evolution_data(const std::string &path, const std::string &output_path, int nb_lev_h2, int nb_vibr_h2o, int nb_lev_h2o, 
		int nb_vibr_co, int nb_lev_co, int nb_lev_pnh3, int nb_lev_onh3, int nb_lev_oh, int nb_vibr_ch3oh, int nb_lev_ch3oh, double c_abund_pah, int verbosity = 1);
	virtual ~evolution_data();
};

class chemistry_evolution_data : public evolution_data
{
public:
	// vector defining the ODE system for the chemistry evolution of the stationary gas cloud.
	int f(realtype t, N_Vector y, N_Vector ydot);
	chemistry_evolution_data(const std::string &path, const std::string &output_path, int nb_lev_h2, int nb_vibr_h2o, int nb_lev_h2o, 
		int nb_vibr_co, int nb_lev_co, int nb_lev_pnh3, int nb_lev_onh3, int nb_lev_oh, int nb_vibr_ch3oh, int nb_lev_ch3oh, double c_abund_pah, int verbosity = 1);

	~chemistry_evolution_data();
};

class mhd_shock_data : public evolution_data
{
public:
	double magn_field_energy, add_el_source, velg_mhd_n, velg_mhd_i, ion_vg_denominator;

	void set_shock_vel(double v) { shock_vel = v; } // in cm/s
	double get_add_electron_sterm() { return add_el_source; }
	double get_velg_mhd_n() { return velg_mhd_n; }
	double get_velg_mhd_i() { return velg_mhd_i; }
    double get_ionvg_deniminator() { return ion_vg_denominator; }

	// vector defining the ODE system for the mhd shock.
	int f(realtype t, N_Vector y, N_Vector ydot);

	mhd_shock_data(const std::string &path, const std::string &output_path, int nb_lev_h2, int nb_vibr_h2o, int nb_lev_h2o, 
		int nb_vibr_co, int nb_lev_co, int nb_lev_pnh3, int nb_lev_onh3, int nb_lev_oh, int nb_vibr_ch3oh, int nb_lev_ch3oh, double c_abund_pah, int verbosity = 1);
	~mhd_shock_data();
};

// The functions used by ODE system solver:
int f_chem(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int f_mhd(realtype t, N_Vector y, N_Vector ydot, void *user_data);
