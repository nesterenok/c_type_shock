
//
// 22.03.2017. Check for errors.
// 02.10.2017. Check for errors.
// 06.02.2018. Check for errors.
// 23.04.2018. Check for errors.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "constants.h"
#include "parameters.h"
#include "parameters.h"
#include "coll_rates_h2.h"
#include "coll_rates_h2o.h"
#include "coll_rates_ions.h"
#include "coll_rates_co.h"
#include "coll_rates_oh.h"
#include "coll_rates_nh3.h"
#include "coll_rates_ch3oh.h"
#include "differential_equations.h"

#define MIN_LINE_OPACITY 1.e-99  // arbitrary very small value;
#define MIN_CARBON_ATOMS_IN_HYDROCARBONS 2 // minimal nb of carbon atoms in hydrocarbon molecules;
// The heat capacity of dust grain is very low, in order to speed up calculations, this factor may be set > 1;
#define DUST_HEAT_CAPACITY_FACTOR 1
using namespace std;

// Structure of the data array:
//  ______________________________________________________________________________________________________________________________________________________________________________
// |Chem. species|H2 pop.  |p-H2O pop.|o-H2O pop.|CO pop.  |OH pop.  | p-NH3     |o-NH3      |CH3OH -A, -E|CI pop.  | OI pop. | CII pop. |grain  | dust temp.    | phys. param.   |
// |nb_of_species|nb_lev_h2|nb_lev_h2o|nb_lev_h2o|nb_lev_co|nb_lev_oh|nb_lev_pnh3|nb_lev_onh3|nb_lev_ch3oh|nb_lev_ci|nb_lev_oi|nb_lev_cii|charge |nb_of_dust_comp| 3 - chem. evol.|
// |_____________|_________|__________|__________|_________|_________|___________|___________|____________|_________|_________|__________|distr__|_______________|_5 - shock______|
//

// The constant necessary for the calculation of energy transfer between ions and electrons;
// Draine, ApJ 241, 1021, 1980; Spitzer, Physics of fully ionized gases, 1962, p. 80;
const double ion_electron_energy_constant = 4.*pow(2.*ELECTRON_MASS*M_PI/BOLTZMANN_CONSTANT, 0.5) *pow(ELECTRON_CHARGE, 4); // *n_e*n_i*(T_e-T_i)/m_i/T_e^1.5
const double log_coulomb_constant = log(1.5*pow(BOLTZMANN_CONSTANT, 1.5)/(sqrt(M_PI)*ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE)); // *T_e^1.5/n_e^0.5
// The parameter value at electron energy 35 eV are taken, H2/He gas mixture is considered, Dalgarno et al., ApJSS 125, 237 (1999);
// at low ionization fraction 0.055 = 1. + (0.055 - 1.)/(1. + 2.17*pow(conc_e/(conc_h2 + conc_he), 0.366));
const double cr_heating_const = 0.055*STANDARD_CR_IONIZ_RATE *35.*EV_TO_ERGS;

// The constants necessary for the calculation of grain velocity, 
// the main gas species H2 and He; 
const double av_mol_mass = (1. + 4.*HE_TO_H_NB_RATIO)*ATOMIC_MASS_UNIT/(0.5 + HE_TO_H_NB_RATIO); // note: may be more accurate calculations?
const double v_const_1 = 9.*M_PI*av_mol_mass/128.; // (v_i - v_e)^2/(k T_n);
const double v_const_2 = v_const_1 *ELECTRON_CHARGE*ELECTRON_CHARGE
	/((1. + 4.*HE_TO_H_NB_RATIO)*(1. + 4.*HE_TO_H_NB_RATIO)*SPEED_OF_LIGHT*SPEED_OF_LIGHT*ATOMIC_MASS_UNIT*ATOMIC_MASS_UNIT); //

//
// The functions for ODE system solver
//
int f_chem(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
	return static_cast<chemistry_evolution_data *>(user_data)->f(t, y, ydot);
}

int f_mhd(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
	return static_cast<mhd_shock_data *>(user_data)->f(t, y, ydot);
}

//
// Base class for the calculation of evolution of chemistry/shock
//

evolution_data::evolution_data(const string &path, const std::string &output_path, int nb_h2, int nb_vibr_h2o, int nb_h2o, 
	int nb_vibr_co, int nb_co, int nb_pnh3, int nb_onh3, int nb_oh, int nb_vibr_ch3oh, int nb_ch3oh, double c_abund_pah, int verb)
	: nb_lev_h2(nb_h2), nb_lev_h2o(nb_h2o), nb_lev_co(nb_co), nb_lev_pnh3(nb_pnh3), nb_lev_onh3(nb_onh3), nb_lev_oh(nb_oh), nb_lev_ch3oh(nb_ch3oh), verbosity(verb), nb_rtd(0),
	cr_ioniz_rate(1.), ir_field_strength(1.), uv_field_strength(1.), visual_extinct(0.), magn_field_0(0.), magn_field(0.), 
	shock_vel(0.), rad_energy_loss_h2(0.), neut_heat_h2(0.), neut_heat_atoms(0.), neut_heat_ph2o(0.), neut_heat_oh2o(0.), 
	neut_heat_co(0.), neut_heat_oh(0.), neut_heat_pnh3(0.), neut_heat_onh3(0.), neut_heat_ch3oh_a(0.), neut_heat_ch3oh_e(0.),
	neut_heat_dust_coll(0.), neut_heat_chem(0.), pheff_gas_heat(0.), neut_cr_heat(0.), neut_heat_scatt_ions(0.), 
	neut_heat_scatt_el(0.), el_heat_atoms(0.), el_heat_h2(0.), el_heat_ph2o(0.), el_heat_oh2o(0.), el_heat_scatt_neut(0.), 
	el_heat_scatt_ions(0.), el_heat_chem(0.), ion_heat_h2(0.), ion_heat_scatt_n(0.), ion_heat_scatt_el(0.), ion_heat_chem(0.), ads_dust_area(0.),
	ads_grain_velz(0.), ads_grain_veln2(0.), desorption_factor_cr(0.), photodes_factor_cr(0.), photodes_factor_is(0.), 
	oh2_form_gaschem(0.), oh2_form_grains(0.), oh2_form_hcoll(0.), h2_h_diss_rate(0.), h2_h_diss_cooling(0.),
	photoem_factor_is_uv(0), photoem_factor_is_vis(0), coll_partn_conc(0), indices(0), dust_heat_h2_line(0), dust_heat_mline(0), 
	dust_heat_coll(0), dust_heat_chem(0), chem_reaction_rates(0), gamma_factors(0), delta_factors(0), dheat_efficiency(0), 
	esc_prob_int1(0), esc_prob_int2(0), ch3oh_a_di(0), ch3oh_e_di(0), ch3oh_a_einst(0), ch3oh_e_einst(0), ch3oh_a_coll(0), 
	ch3oh_e_coll(0)
{
	int i, isotope, angm_ch3oh_max;
	bool update;
	double mass, vmin, vmax, dg_ratio, spin;
	
	// parameter necessary for radiative transfer calculations, 
	// the discussion on this parameter see in Bergin, Tafalla, Annu. Rev. Astron. Astrophys. 45, p. 339 (2007), section 2.5:
	vel_turb = 2.e+4; // cm/s
	
	// the ratio of velocity dispersion divided by characterictic scale, in cm/s/cm;
	// this parameter can be estimated from virial equilibrium for a uniform density sphere, see Goldsmith, ApJ 557, p. 736 (2001)
	// characteristic value for molecular clouds, dv/dz = 1 km s-1 pc-1 = 3.2e-14 cm s-1 cm-1;  
	vel_grad_min = 3.e-14;
	vel_n_grad = vel_i_grad = - vel_grad_min;

	// Initialization of cross section data for elastic scattering:
	// Please, see the paper by Dalgarno 1999 and reference therein (Shimamura et al. 1989),
	elastic_cross_section_table	*elastic_h2_el_cs
		= new elastic_cross_section_table(path, "elastic_scattering/e-h2_mt_cs.txt");

	elastic_cross_section_table	*elastic_he_el_cs
		= new elastic_cross_section_table(path, "elastic_scattering/e-he_mt_cs.txt");

	elastic_cross_section_table	*elastic_h_el_cs
		= new elastic_cross_section_table(path, "elastic_scattering/e-h_mt_cs.txt");
	
	mass = 29.*ATOMIC_MASS_UNIT; 
	
	// The first argument is polarizability in [bohr radii ^3];
	elastic_cross_section_powerlaw *elast_h2_ion_cs
		= new elastic_cross_section_powerlaw(5.315, 2.*ATOMIC_MASS_UNIT, mass);

	elastic_cross_section_const *elast_h_ion_cs
		= new elastic_cross_section_const(4.5, ATOMIC_MASS_UNIT, mass);

	elastic_cross_section_const *elast_he_ion_cs
		= new elastic_cross_section_const(1.383, 4.*ATOMIC_MASS_UNIT, mass);

	// Initialization of elastic scattering objects. Velocity is given in cm/s;
	elastic_h2_el = new elastic_scatt_neutral_charged(2.*ATOMIC_MASS_UNIT, ELECTRON_MASS, vmin = 5.e+5, vmax = 2.e+8, verbosity);
	elastic_he_el = new elastic_scatt_neutral_charged(4.*ATOMIC_MASS_UNIT, ELECTRON_MASS, vmin = 5.e+5, vmax = 2.e+8, verbosity);
	elastic_h_el = new elastic_scatt_neutral_charged(ATOMIC_MASS_UNIT, ELECTRON_MASS, vmin = 5.e+5, vmax = 2.e+8, verbosity);
	
	elastic_h2_el->calc_data(elastic_h2_el_cs);
	elastic_he_el->calc_data(elastic_he_el_cs);
	elastic_h_el->calc_data(elastic_h_el_cs);
	
	elastic_h_ions = new elastic_scattering(ATOMIC_MASS_UNIT, mass, verbosity);
	elastic_he_ions = new elastic_scattering(4.*ATOMIC_MASS_UNIT, mass, verbosity);
	
	if (FLOWER_ELASTIC_APPROX) {
		elastic_h2_ions = new elastic_scatt_neutral_charged(2.*ATOMIC_MASS_UNIT, mass, vmin = 1.e+4, vmax = 3.e+6, verbosity);
	} 
	else { // the approximation of constant cross section (using result by Osterbrock, 1961);	
		elastic_h2_ions = new elastic_scattering(2.*ATOMIC_MASS_UNIT, mass, verbosity);
	}
	
	elastic_h2_ions->calc_data(elast_h2_ion_cs);
	elastic_h_ions->calc_data(elast_h_ion_cs);
	elastic_he_ions->calc_data(elast_he_ion_cs);
	
	// H2 molecule data:
	molecule h2_mol("H2", isotope = 1, 2.*ATOMIC_MASS_UNIT);
	
	h2_di = new h2_diagram(path, h2_mol, nb_lev_h2, verbosity);
	h2_einst = new h2_einstein_coeff(path, h2_di, verbosity);
	
	h2_collisions *h2_coll_aux 
		= new h2_collisions(path, h2_di, verbosity);

	h2_coll = h2_coll_aux;
	// h2_coll_aux->check_spline(2, 0, output_path + "h2_spline.txt");

	h2_h_diss_data = new h2_h_dissociation_data(path, h2_di, verbosity);

	// Calculation of the data of H2 excitation processes:
	h2_excit_gf = new h2_grain_formation(h2_di);
	h2_excit_gasph = new h2_gasphase_formation(h2_di);
	h2_excit_cr = new h2_excit_cosmic_rays(path, h2_di);

	// Ion data:
	mass = 16.*ATOMIC_MASS_UNIT;
	molecule ion_OI("OI", isotope = 1, mass);

	mass = 12.*ATOMIC_MASS_UNIT;
	molecule ion_CI("CI", isotope = 1, mass);
	molecule ion_CII("CII", isotope = 1, mass);

	nb_lev_oi = 5;
	OI_di = new ion_diagram(path, ion_OI, nb_lev_oi, verbosity);
	OI_einst = new ion_einstein_coeff(path, OI_di, verbosity);
	OI_coll = new OI_collisions(path, OI_di, verbosity);

	nb_lev_ci = 3;
	CI_di = new ion_diagram(path, ion_CI, nb_lev_ci, verbosity);
	CI_einst = new ion_einstein_coeff(path, CI_di, verbosity);
	CI_coll = new CI_collisions(path, CI_di, verbosity);

	nb_lev_cii = 2;
	CII_di = new ion_diagram(path, ion_CII, nb_lev_cii, verbosity);
	CII_einst = new ion_einstein_coeff(path, CII_di, verbosity);
	CII_coll = new CII_collisions(path, CII_di, verbosity);

	// H2O molecule data:
	mass = 18.*ATOMIC_MASS_UNIT;
	molecule ph2o_mol("pH2O", isotope = 1, mass, spin = 0.);
	molecule oh2o_mol("oH2O", isotope = 1, mass, spin = 1.);

	h2o_diagram *di = new h2o_diagram(path, ph2o_mol, nb_lev_h2o, nb_vibr_h2o, verbosity);
	ph2o_di = di;
	
	// the pointer to derivative class must be used here (not to base class):
	ph2o_einst = new h2o_einstein_coeff(path, di, verbosity);
	ph2o_coll =	new h2o_collisions(path, ph2o_di, false, verbosity);

	di = new h2o_diagram(path, oh2o_mol, nb_lev_h2o, nb_vibr_h2o, verbosity);
	oh2o_di = di;
	
	oh2o_einst = new h2o_einstein_coeff(path, di, verbosity);
	oh2o_coll =	new h2o_collisions(path, oh2o_di, false, verbosity);

	// CO molecule data:
	mass = 28.*ATOMIC_MASS_UNIT;
	molecule co_mol("CO", isotope = 1, mass, spin = 0.);

	co_di = new co_diagram(path, co_mol, nb_lev_co, nb_vibr_co, verbosity);
	co_einst = new co_einstein_coeff(path, co_di, verbosity);
	co_coll = new co_collisions(path, co_di, verbosity);

	// OH molecule data:
	if (nb_lev_oh > 20)
        nb_lev_oh = 20;
	
    mass = 17.*ATOMIC_MASS_UNIT;
	molecule oh_mol("OH", isotope = 1, mass, spin = 0.);

#if (CALCULATE_POPUL_NH3_OH)
	oh_di = new oh_diagram(path, oh_mol, nb_lev_oh, verbosity);
	oh_einst = new oh_einstein_coeff(path, oh_di, verbosity);
	oh_coll = new oh_collisions(path, oh_di, verbosity);
#endif
	// NH3 molecule data
	if (nb_lev_onh3 > 17)
        nb_lev_onh3 = 17; // ortho-NH3: He coll data - 22, H2 coll data - 17
	
    if (nb_lev_pnh3 > 34)
        nb_lev_pnh3 = 34; // para-NH3: He coll data - 16, H2 coll data - 34

	mass = 17.*ATOMIC_MASS_UNIT;
	molecule onh3_mol("oNH3", isotope = 1, mass, spin = 1.5);
	molecule pnh3_mol("pNH3", isotope = 1, mass, spin = 0.5);

#if (CALCULATE_POPUL_NH3_OH)
	onh3_di = new nh3_diagram(path, onh3_mol, nb_lev_onh3, verbosity);
	onh3_einst = new nh3_einstein_coeff(path, onh3_di, verbosity);
	onh3_coll = new nh3_collisions(path, onh3_di, verbosity);

	pnh3_di = new nh3_diagram(path, pnh3_mol, nb_lev_pnh3, verbosity);
	pnh3_einst = new nh3_einstein_coeff(path, pnh3_di, verbosity);
	pnh3_coll = new nh3_collisions(path, pnh3_di, verbosity);
#endif
	// CH3OH molecule data
	angm_ch3oh_max = 15; // the available collisional data is restricted to this value;
	mass = 32.*ATOMIC_MASS_UNIT;

	molecule ch3oh_a_mol("ch3oh_a", isotope = 1, mass, spin = 1.5);
	molecule ch3oh_e_mol("ch3oh_e", isotope = 1, mass, spin = 0.5);

#if (CALCULATE_POPUL_METHANOL)
	ch3oh_a_di = new ch3oh_diagram(path, ch3oh_a_mol, nb_lev_ch3oh, nb_vibr_ch3oh, angm_ch3oh_max);
	ch3oh_a_einst = new ch3oh_einstein_coeff(path, ch3oh_a_di); 
	ch3oh_a_coll = new ch3oh_collisions(path, ch3oh_a_di, verbosity);
	
	ch3oh_e_di = new ch3oh_diagram(path, ch3oh_e_mol, nb_lev_ch3oh, nb_vibr_ch3oh, angm_ch3oh_max);
	ch3oh_e_einst = new ch3oh_einstein_coeff(path, ch3oh_e_di);
	ch3oh_e_coll = new ch3oh_collisions(path, ch3oh_e_di, verbosity);
#endif

	// LVG data:
	loss_func_line_phot = new lvg_method_data(path, "lvg_loss_func.txt", verbosity);
	// the data is used for which integration range for mu is [0.01;1],
	loss_func_cont_phot = new lvg_method_data(path, "lvg_loss_func_qt1_mu1e-2.txt", verbosity);

	// 2-component dust model:
	dust = new two_component_dust_model(path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, 
		verbosity);
	nb_of_dust_comp = dust->nb_of_comp;

	// Arrays for dust heating rates by interstellar radiation field:
	dh_isrf_arr = new double [nb_of_dust_comp];
	memset(dh_isrf_arr, 0, nb_of_dust_comp*sizeof(double));

	// by molecular emission:
	dust_heat_mline = new double [nb_of_dust_comp]; // atoms and molecules other than H2
	memset(dust_heat_mline, 0, nb_of_dust_comp*sizeof(double));

	dust_heat_h2_line = new double [nb_of_dust_comp];
	memset(dust_heat_h2_line, 0, nb_of_dust_comp*sizeof(double));
	
	// by collisions with gas:
	dust_heat_coll = new double [nb_of_dust_comp];
	memset(dust_heat_coll, 0, nb_of_dust_comp*sizeof(double));

	// by chemical reactions on grain surface:
	dust_heat_chem = new double [nb_of_dust_comp];
	memset(dust_heat_chem, 0, nb_of_dust_comp*sizeof(double));

	dheat_isrf = new dust_heating_ISRF();
	dheat_isrf->calculate(path, dust, STANDARD_NB_CR_PHOTONS);
  
	// parameters have default values at this point:
	dheat_isrf->get_heating(visual_extinct, ir_field_strength, uv_field_strength, cr_ioniz_rate, dh_isrf_arr, nb_of_dust_comp); 
	
	// Initialization of chemistry:
	network = new chem_network(path, verbosity);
	network->add_element("H");
	network->add_element("He");
	network->add_element("O");
	network->add_element("C");
	network->set_max_nb_carbon(11);
	network->add_element("N");
	network->add_element("Si");
	network->add_element("S");
	network->add_element("Fe");
	network->add_element("Na");
	network->add_element("Mg");
	network->add_element("Cl");
	// P and F may be not to be included in order to reduce computation time (as in Walsh et al., ApJ 722, 1607, 2010)
//	network->add_element("P");
//	network->add_element("F");
	
	network->init_gas_phase_species(path + "chemistry/UMIST_2012/species_UMIST2012.txt");
//	network->init_gmantle_species(path + "chemistry/UMIST_2012/surface_binding_energies_NAUTILUS.txt");
	network->init_gmantle_species(path + "chemistry/UMIST_2012/surface_binding_energies_Penteado2017.txt");
	
	// in order to prevent adsorbtion of large amount of H2 molecules on grains:
	if (GRAIN_SURFACE_CHEMISTRY_ON == 0) {
		network->exclude_chem_specimen("*H2");
	}
	network->init_species_nbs();

	// Main network and updates:
	network->init_network_umistf(path + "chemistry/UMIST_2012/rates_UMIST2012.txt", update = false);
	network->init_network_umistf(path + "chemistry/UMIST_2012/rates_update_CRphoto_Heys2017.txt");
	network->init_network_umistf(path + "chemistry/UMIST_2012/rates_update_ISRFphoto_Heys2017.txt");
	network->init_network_umistf(path + "chemistry/UMIST_2012/rates_update_Chabot2013.txt");
	
	network->init_network_umistf(path + "chemistry/reactions_adsorp_desorption.txt");
	network->init_network_umistf(path + "chemistry/reactions_ion_recomb_grains.txt");
	
	// Palau et al., MNRAS 467, p.2723 (2017), arXiv:1701.04802v1;
	network->init_network_umistf(path + "chemistry/reactions_COM_Palau2017.txt");
	// additional collisional dissociation reactions;
	network->init_network_umistf(path + "chemistry/reactions_colldiss.txt");

	// ion neutral reactions for CH3O and CH2OH and other radicals:
	network->init_network_umistf(path + "chemistry/reactions_ion_radicals.txt");
	
	if (H2_FORMATION_MODE > 0)
		network->init_h2_formation();
	
	if (GRAIN_SURFACE_CHEMISTRY_ON > 0) {
		network->init_photoreact_surface_chemistry();
		network->init_grain_surface_chemistry(path + "chemistry/nautilus_networks/bimolecular_gs_reactions.txt");
	//	network->init_grain_surface_chemistry(path + "chemistry/nautilus_networks/bimolecular_gs_reactions_add.txt");
	//	network->init_grain_surface_chemistry(path + "chemistry/Belloche_Garrod_2014/bimolecular_gs_reactions.txt");
	}

	accr_func = new accretion_rate_functions();
	
	nb_of_species = network->nb_of_species;
	nb_of_gmantle_species = network->nb_of_gmantle_species;
	nb_of_gas_species = nb_of_species - nb_of_gmantle_species;  
	
	// auxulary variables that facilitate access to the data;
	nb_dch = new int [nb_of_dust_comp+1];
	memset(nb_dch, 0, (nb_of_dust_comp+1)*sizeof(int));

	min_grain_charge = new int [nb_of_dust_comp];
	max_grain_charge = new int [nb_of_dust_comp];
	
	nb_dch[0] = nb_of_species + nb_lev_h2 + 2*nb_lev_h2o + nb_lev_co + nb_lev_oh + nb_lev_pnh3 + nb_lev_onh3 
		+ 2*nb_lev_ch3oh + nb_lev_ci + nb_lev_oi + nb_lev_cii;

	// default values for min and max charge numbers:
	for (i = 0; i < nb_of_dust_comp; i++)
	{
		if (dust->components[i]->zmin > -NB_CHARGE_BINS_LARGE_GRAINS/2)
			min_grain_charge[i] = dust->components[i]->zmin;
		else min_grain_charge[i] = -NB_CHARGE_BINS_LARGE_GRAINS/2;

		if (dust->components[i]->zmax < NB_CHARGE_BINS_LARGE_GRAINS/2)
			max_grain_charge[i] = dust->components[i]->zmax;
		else max_grain_charge[i] = NB_CHARGE_BINS_LARGE_GRAINS/2;
		
		nb_dch[i+1] = nb_dch[i] + max_grain_charge[i] - min_grain_charge[i] + 1;
	}
	
	nb_of_grain_charges = nb_dch[nb_of_dust_comp] - nb_dch[0];
	nb_dct = nb_dch[nb_of_dust_comp];
	nb_mhd = nb_dct + nb_of_dust_comp;
	
	is_dust_charge_large = new bool [nb_of_dust_comp];
	memset(is_dust_charge_large, 0, nb_of_dust_comp*sizeof(bool));

	// by definition, only small charged grains are treated as ions:
	is_charged_dust_as_ions = new bool [nb_of_dust_comp];
	for (i = 0; i < nb_of_dust_comp; i++) 
	{
		if (dust->get_grain_radius(i) < MAX_ION_COUPLED_GRAIN_RADIUS)
			is_charged_dust_as_ions[i] = true;
		else is_charged_dust_as_ions[i] = false;
	}

	alpha = new double [nb_of_grain_charges];
	beta = new double [nb_of_grain_charges];
	wt2_arr = new double [nb_of_grain_charges];
	grain_velz = new double [nb_of_grain_charges];
	grain_velz_g = new double [nb_of_grain_charges];
	grain_veln2 = new double [nb_of_grain_charges];
	grain_veli2 = new double [nb_of_grain_charges];

	av_grain_velz = new double [nb_of_dust_comp];
	grain_conc = new double [nb_of_dust_comp];

	phel_rate_cr = new double [nb_of_dust_comp];
	phel_rate_uv = new double [nb_of_dust_comp];
	phel_rate_vis = new double [nb_of_dust_comp];
	ion_neutr_rate = new double [nb_of_dust_comp];
	el_att_rate = new double [nb_of_dust_comp];

#if (SAVE_RADIATIVE_TRANSFER_FACTORS)
	nb_rtd = nb_lev_h2*(nb_lev_h2-1)/2 + nb_lev_h2o*(nb_lev_h2o-1) + nb_lev_co*(nb_lev_co-1)/2 + nb_lev_oh*(nb_lev_oh-1)/2 
		+ nb_lev_pnh3*(nb_lev_pnh3-1)/2 + nb_lev_onh3*(nb_lev_onh3-1)/2	+ nb_lev_ch3oh*(nb_lev_ch3oh-1) + 
		nb_lev_ci*(nb_lev_ci-1)/2 + nb_lev_oi*(nb_lev_oi-1)/2 + nb_lev_cii*(nb_lev_cii-1)/2;

	gamma_factors = new double [nb_rtd];
	memset(gamma_factors, 0, nb_rtd*sizeof(double));

	delta_factors = new double [nb_rtd];
	memset(delta_factors, 0, nb_rtd*sizeof(double));

	dheat_efficiency = new double [nb_rtd];
	memset(dheat_efficiency, 0, nb_rtd*sizeof(double));

	esc_prob_int1 = new double [nb_rtd];
	memset(esc_prob_int1, 0, nb_rtd*sizeof(double));

	esc_prob_int2 = new double [nb_rtd];
	memset(esc_prob_int2, 0, nb_rtd*sizeof(double));
#endif

	delete elastic_h2_el_cs;
	delete elastic_he_el_cs;
	delete elastic_h_el_cs;

	delete elast_h2_ion_cs;
	delete elast_h_ion_cs;
	delete elast_he_ion_cs;
}

evolution_data::~evolution_data()
{
	delete elastic_h2_ions;
	delete elastic_h_ions;
	delete elastic_he_ions;
	
	delete elastic_h2_el;
	delete elastic_he_el;
	delete elastic_h_el;

	delete h2_di;
	delete h2_einst;
	delete h2_coll;
	delete h2_excit_gf;
	delete h2_excit_gasph;
	delete h2_excit_cr;
	
	delete OI_di;
	delete OI_einst;
	delete OI_coll;

	delete CI_di;
	delete CI_einst;
	delete CI_coll;

	delete CII_di;
	delete CII_einst;
	delete CII_coll;

	delete ph2o_di;
	delete oh2o_di;

	delete ph2o_einst;
	delete oh2o_einst;

	delete ph2o_coll;
	delete oh2o_coll;

	delete co_di;
	delete co_einst;
	delete co_coll;

	delete oh_di;
	delete oh_einst;
	delete oh_coll;

	delete onh3_di;
	delete onh3_einst;
	delete onh3_coll;

	delete pnh3_di;
	delete pnh3_einst;
	delete pnh3_coll;

	delete ch3oh_a_di;
	delete ch3oh_a_einst;
	delete ch3oh_a_coll;

	delete ch3oh_e_di;
	delete ch3oh_e_einst;
	delete ch3oh_e_coll;

	delete loss_func_line_phot;
	delete loss_func_cont_phot;
	delete dust;
	delete dheat_isrf;
	delete network;
	delete accr_func;

	delete [] dh_isrf_arr;
	delete [] indices;
	delete [] coll_partn_conc;

	delete [] nb_dch;
	delete [] min_grain_charge;
	delete [] max_grain_charge;
	delete [] is_dust_charge_large;
	delete [] is_charged_dust_as_ions;
	delete [] av_grain_velz;
	delete [] grain_conc;

	delete [] alpha;
	delete [] beta;
	delete [] wt2_arr;
	delete [] grain_velz;
	delete [] grain_veln2;
	delete [] grain_veli2;
	delete [] grain_velz_g;
	
	delete [] dust_heat_mline; 
	delete [] dust_heat_h2_line;
	delete [] dust_heat_coll;
	delete [] dust_heat_chem;
	// these arrays are initialized in the daughter classes:
	delete [] chem_reaction_rates;
	delete [] chem_heating_rates_n;

	delete [] phel_rate_cr;
	delete [] phel_rate_uv;
	delete [] phel_rate_vis;
	delete [] ion_neutr_rate;
	delete [] el_att_rate;

	delete [] gamma_factors;
	delete [] delta_factors;
	delete [] dheat_efficiency;
	delete [] esc_prob_int1;
	delete [] esc_prob_int2;
}

void evolution_data::reinit_arrays(int nb_of_grain_charges)
{
	delete [] alpha;
	delete [] beta;
	delete [] wt2_arr;
	delete [] grain_velz;
	delete [] grain_veln2;
	delete [] grain_veli2;
	delete [] grain_velz_g;

	alpha = new double [nb_of_grain_charges];
	beta = new double [nb_of_grain_charges];
	wt2_arr = new double [nb_of_grain_charges];
	grain_velz = new double [nb_of_grain_charges];
	grain_velz_g = new double [nb_of_grain_charges];
	grain_veln2 = new double [nb_of_grain_charges];
	grain_veli2 = new double [nb_of_grain_charges];
}

int evolution_data::f(realtype t, N_Vector y, N_Vector ydot)
{
	int i, j, k, l, nb, nb2;
	double a, b, c, d, rate, mom_n, mom_i, mom_e, en_n, en_i, en_e, en_d, h2_pr1, h2_pr2, h2_d1, el_att, phel_em, nb_ads_sites_grain,
		photoem_factor_cr, photodes_factor_cr2, photoreact_factor_cr, ioniz_fraction;
	realtype *y_data, *ydot_data;
	
	y_data = NV_DATA_S(y);
	ydot_data = NV_DATA_S(ydot);

    // control unphysical negative specimen concentrations (and level population densities)
 /*   for (k = 0; k < nb_of_species; k++) {
        if (y_data[k] < 0.)
            return 1;
    }*/
	for (k = 0; k < nb_of_equat; k++) {
		ydot_data[k] = 0.;
	}

	// velocity values are initialized in the function of daughter class
	temp_n = y_data[nb_mhd];
	temp_i = y_data[nb_mhd + 1];
	temp_e = y_data[nb_mhd + 2];
	
	temp_n_erg = temp_n *BOLTZMANN_CONSTANT;
	temp_i_erg = temp_i *BOLTZMANN_CONSTANT;
	temp_e_erg = temp_e *BOLTZMANN_CONSTANT;
	
	conc_h2 = y_data[network->h2_nb];
	conc_h = y_data[network->h_nb];
	conc_he = y_data[network->he_nb];
	conc_e = y_data[network->e_nb];
	conc_h2o = y_data[network->h2o_nb];
	conc_co	= y_data[network->co_nb];
	conc_oh	= y_data[network->oh_nb];
	conc_nh3 = y_data[network->nh3_nb];
	conc_ch3oh = y_data[network->ch3oh_nb];

	conc_oi = y_data[network->oi_nb]; 
	conc_ci = y_data[network->ci_nb]; 
	conc_cii = y_data[network->cii_nb];

    conc_hp = y_data[network->hp_nb];
    conc_h3p = y_data[network->h3p_nb];
	// only H2, H are taken into account in the calculation of H nuclei concentration:
	conc_h_tot = 2.*conc_h2 + conc_h;
	
	// ion mass and concentration do not include PAH;
	ion_mass = conc_n = conc_i = conc_ice = 0.;
	for (i = 0; i < nb_of_gas_species; i++) {
		if (network->species[i].type == "neutral") {
			conc_n += y_data[i];
		}	
		else if (network->species[i].type == "ion") 
		{ 
			conc_i += y_data[i];
			ion_mass += y_data[i] *network->species[i].mass;
		}
	}
	ion_mass /= conc_i;
	
	for (i = nb_of_gas_species; i < nb_of_gas_species + nb_of_gmantle_species; i++) {
		conc_ice += y_data[i];
	}
	
	// evaluation of the concentrations of ortho- and para-H2:
	conc_ph2 = 0.;
	for (i = nb_of_species; i < nb_of_species + nb_lev_h2; i++) {
		if (rounding(h2_di->lev_array[i-nb_of_species].j)%2 == 0)
			conc_ph2 += y_data[i];
	}
	conc_oh2 = conc_h2 - conc_ph2;
	conc_h2j0 = y_data[nb_of_species]; // concentration of H2 in the level j=0;

	// Calculation of dust component velocities and concentrations:
	for (i = 0; i < nb_of_dust_comp; i++) 
	{
		if (is_charged_dust_as_ions[i]) {
			if (is_dust_charge_large[i])
			{	
				l = nb_dch[i] - nb_dch[0];
				av_grain_velz[i] = grain_velz[l] = vel_i;
				grain_conc[i] = y_data[nb_dch[i]];
				
				grain_veln2[l] = (vel_i - vel_n)*(vel_i - vel_n);
				grain_veli2[l] = 0.;
			}
			else
			{
				grain_conc[i] = av_grain_velz[i] = 0.;
				for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) 
				{
					l = nb_dch[i] - nb_dch[0] + k;
					if (k == -min_grain_charge[i]) {
						grain_velz[l] = vel_n;
						grain_veln2[l] = 0.;
						grain_veli2[l] = (vel_i - vel_n)*(vel_i - vel_n);
					}
					else {
						grain_velz[l] = vel_i;
						grain_veln2[l] = (vel_i - vel_n)*(vel_i - vel_n);
						grain_veli2[l] = 0.;
					}
					grain_conc[i] += y_data[nb_dch[i] + k];
					av_grain_velz[i] += grain_velz[l] *y_data[nb_dch[i] + k];
				}
				av_grain_velz[i] /= grain_conc[i];
			}
		}
		else {
			if (is_dust_charge_large[i])
			{
				l = nb_dch[i]-nb_dch[0];
				grain_conc[i] = y_data[nb_dch[i]];
				
				av_grain_velz[i] = grain_velz[l] = 
					calc_grain_velocity(y, y_data[nb_dch[i] + 1], dust->get_grain_area(i), alpha[l], beta[l], wt2_arr[l], c);
				
				grain_veln2[l] = (grain_velz[l] - vel_n)*(grain_velz[l] - vel_n) + c*c;
				grain_veli2[l] = (grain_velz[l] - vel_i)*(grain_velz[l] - vel_i) + c*c;
			}
			else 
			{
				grain_conc[i] = av_grain_velz[i] = 0.;
				for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) 
				{
					l = nb_dch[i] - nb_dch[0] + k;
		
					grain_velz[l] = 
						calc_grain_velocity(y, min_grain_charge[i] + k, dust->get_grain_area(i), alpha[l], beta[l], wt2_arr[l], c);
				
					grain_conc[i] += y_data[nb_dch[i] + k];
					av_grain_velz[i] += grain_velz[l] *y_data[nb_dch[i] + k];
					
					grain_veln2[l] = (grain_velz[l] - vel_n)*(grain_velz[l] - vel_n) + c*c;
					grain_veli2[l] = (grain_velz[l] - vel_i)*(grain_velz[l] - vel_i) + c*c;
				}
				av_grain_velz[i] /= grain_conc[i];
			}
		}
	}

	// here the adsorption-desorption area and corresponding dust temperature are calculated;
	b = temp_d = ads_dust_area = ads_grain_velz = ads_grain_veln2 = 0.;
	for (i = 0; i < nb_of_dust_comp; i++) {
		if (dust->get_grain_radius(i) > MIN_ADSORPTION_RADIUS) 
		{
			b += grain_conc[i];
			a = dust->get_grain_area(i) *grain_conc[i];
			ads_dust_area += a;
		
			temp_d += y_data[nb_dct + i] *a;
			ads_grain_velz += av_grain_velz[i] *a; //
			
			if (is_dust_charge_large[i]) {
				ads_grain_veln2 += a*grain_veln2[nb_dch[i]-nb_dch[0]];
			}
			else {
				for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) {
					ads_grain_veln2 += dust->get_grain_area(i) *y_data[nb_dch[i] + k] *grain_veln2[nb_dch[i] - nb_dch[0] + k];
				}
			}
		}
	}
	nb_ads_sites_grain = 4.*ads_dust_area*GRAIN_SITES_PER_CM2/b;
	ads_grain_velz /= ads_dust_area;
	ads_grain_veln2 /= ads_dust_area;
	
	temp_d /= ads_dust_area;
	temp_d_erg = temp_d*BOLTZMANN_CONSTANT;

	// only two surface layers of adsorbed species are accounted for;
	coverage = 8.*GRAIN_SITES_PER_CM2*ads_dust_area/(conc_ice + DBL_EPSILON);
	if (coverage > 1.) 
		coverage = 1.;

	// auxulary variables to calculate reaction rates:
	a = 2.*conc_h2/conc_h_tot;
	photodes_factor_cr2 = photodes_factor_cr *a;
	photoreact_factor_cr = cr_ioniz_rate *2.*a; // grain albedo is assumed to equal to 0.5;

	memset(chem_reaction_rates, 0, network->nb_of_reactions*sizeof(double));
	memset(chem_heating_rates_n, 0, network->nb_of_reactions*sizeof(double));
	mom_n = mom_i = mom_e = en_n = en_i = en_e = en_d = h2_pr1 = h2_pr2 = h2_d1 = 0.;

#pragma omp parallel reduction(+: mom_n, mom_i, mom_e, en_n, en_i, en_e, en_d, h2_pr1, h2_pr2, h2_d1) private(i, k, l, rate, a)
	{
		double *arr = new double [nb_of_species];
		memset(arr, 0, nb_of_species*sizeof(double));

		double *arr_chr = new double [network->nb_of_reactions];
		memset(arr_chr, 0, network->nb_of_reactions*sizeof(double));

		double *arr_chhe = new double [network->nb_of_reactions];
		memset(arr_chhe, 0, network->nb_of_reactions*sizeof(double));

        
#pragma omp for schedule(dynamic, 10)
		for (i = 0; i < network->nb_of_reactions - network->nb_reactions_ion_grains; i++)
		{
            const chem_reaction & reaction = network->reaction_array[i];
			rate = reaction_rate(network->species, accr_func, reaction, cr_ioniz_rate, uv_field_strength,
				visual_extinct, vel_n, vel_i, ads_dust_area, ads_grain_veln2, temp_n, temp_i, temp_e, temp_d, photoreact_factor_cr, 
				desorption_factor_cr, photodes_factor_cr2, photodes_factor_is, conc_h_tot, nb_ads_sites_grain);
		
            if (rate < 0.) {
                //cout << left << "warning, rate in reaction is negative: " << setw(13) << rate << network->reaction_array[i].name << endl;
            }

			switch (reaction.type)
			{ 
			// Unimolecular reactions:
			case 0: // "A + CR -> Ion + e-"
			case 1: // "A + CR -> B + Ion + e-"
			case 2: // "A + CR -> Ion + Ion"
			case 3: // "A + CR -> B + C"
			case 4: // "A + CRPhoton -> Ion + e-"
			case 5: // "A + CRPhoton -> B + Ion + e-"
			case 6: // "A + CRPhoton -> B + C + D + E"
			case 7: // "Ion + CRPhoton -> B + Ion"
			case 8: // "Ion + CRPhoton -> A + e-"
			case 9: // "A + ISPhoton -> Ion + e-"
			case 10: // "A + ISPhoton -> B + Ion + e-"
			case 11: // "A + ISPhoton -> B + C + D + E"
			case 12: // "Ion + ISPhoton -> B + Ion"
			case 13: // "Ion + ISPhoton -> A + e-"
			case 27: // "A + grain -> *A + grain"
			case 40: // "*A + CRPhoton -> *B + *C"
			case 41:
			case 42: // "*A + ISPhoton -> *B + *C"
			case 43:
				k = reaction.reactant[0];
				rate *= y_data[k];
				arr[k] -= rate;
				break;

			// Desorption:
			// the desorption is from the top two layers;
			case 28: // "*A + grain -> A + grain"
			case 29: // "*A + CR -> A"
			case 30: // "*A + CRPhoton -> A"
			case 31: // "*A + CRPhoton -> A + B"
			case 32: // "*A + ISPhoton -> A"
			case 33: // "*A + ISPhoton -> A + B"			
				k = reaction.reactant[0];
				rate *= y_data[k]*coverage;
				arr[k] -= rate;
				break;

			// Bimolecular reactions:
			case 14: // "A + B -> C + D + E"
			case 15: // "A + B -> C + Photon"
			case 16: // "A + B -> Ion + e-"
			case 17: // "A + e- -> Ion"
			case 18: // "A + e- -> B + C + Ion"
			case 19: // "A + e- -> B + C + e-"
			case 20: // "A + Ion -> Ion + Photon"
			case 21: // "A + Ion -> B + C + D + Ion"
			case 22: // "A + Ion -> B + C + e-"
			case 23: // "Ion + e- -> A + Photon"
			case 24: // "Ion + e- -> A + B + C"
			case 25: // "Ion + Ion -> A + B + C + D"
			case 36: // "*A + *B -> *C + *D + *E"
			case 37: // "*A + *B -> C"
			case 38: // "*A + *B -> C + D + E"
			case 39: // "*A + *B -> *A + B"
				k = reaction.reactant[0];
				l = reaction.reactant[1];

				rate *= y_data[k] *y_data[l];

                if (H2_H_DISSIOCIATION_DATA) {
                    if (i == network->h2_h_diss_nb)
                        rate = 0.;
                }
                else {
                    if (i == network->h2_h_diss_nb)
                        h2_d1 = rate;
                }

				arr[k] -= rate;
				arr[l] -= rate;
				break;
		
			// Sputtering:
			case 34: // "*A + C -> A + C"
			case 35: // "*A + C -> A + B + C"
				k = reaction.reactant[0];
				l = reaction.reactant[1];

				rate *= y_data[k] *y_data[l] *coverage;

				arr[k] -= rate;
				arr[l] -= rate;
				break;
			
			case 26: // "H + H + grain -> H2 + grain"
				rate *= conc_h;
				arr[network->h_nb] -= 2.*rate;
				break;	
			};

			for (k = 0; k < reaction.nb_of_products; k++) {
				l = reaction.product[k];
				arr[l] += rate;
			}

			// saving rates of chemical reactions (in cm-3 s-1):
			arr_chr[i] = rate;
			
			a = en_n;
			chemistry_source_terms(mom_n, mom_i, mom_e, en_n, en_i, en_e, en_d, network->species, 
				reaction, vel_n, vel_i, ads_grain_velz, ads_grain_veln2, temp_n_erg, temp_i_erg, temp_e_erg, temp_d_erg, rate);

			// saving the rates of chemical heating of neutral fluid (erg cm-3 s-1):
			arr_chhe[i] = en_n - a;
			
			// saving H2 formation rates on grains:
			if (reaction.type >= 36 && reaction.type <= 38) {
				for (k = 0; k < reaction.nb_of_products; k++) {
					if (reaction.product[k] == network->ah2_nb || reaction.product[k] == network->h2_nb)
						h2_pr1 += rate;
				}
			}
			if (reaction.type == 26)
				h2_pr1 += rate;
			
			// saving H2 formation rates in the gas phase (note: the reaction H2+H2->H2+H+H)
			if (reaction.type <= 25) {
				for (k = 0; k < reaction.nb_of_products; k++) {
					if (reaction.product[k] == network->h2_nb)
						h2_pr2 += rate;
				}
			}
		}
		#pragma omp critical
		{
			for (i = 0; i < nb_of_species; i++) {
				ydot_data[i] += arr[i];
			}
			for (i = 0; i < network->nb_of_reactions - network->nb_reactions_ion_grains; i++) {
				chem_reaction_rates[i] += arr_chr[i];
			}
			for (i = 0; i < network->nb_of_reactions - network->nb_reactions_ion_grains; i++) {
				chem_heating_rates_n[i] += arr_chhe[i];
			}
		}
		delete [] arr;
		delete [] arr_chr;
		delete [] arr_chhe;
	}

	mom_gain_n = mom_n;
	mom_gain_i = mom_i;
	mom_gain_e = mom_e;

	neut_heat_chem = energy_gain_n = en_n;
	ion_heat_chem = energy_gain_i = en_i;
	el_heat_chem = energy_gain_e = en_e;
	energy_gain_d = en_d;
	
	h2_prod_gr = h2_pr1;
	h2_prod_gas = h2_pr2;
	h2_h_diss_rate = h2_d1;

	// Chemical reactions changing the charge distribution of grains:
	photoem_factor_cr = cr_ioniz_rate *2.*conc_h2/conc_h_tot; // grain albedo is taken into account yet;
	mom_n = mom_i = mom_e = en_n = en_i = el_att = phel_em = pheff_gas_heat = 0.;
	for (i = 0; i < nb_of_dust_comp; i++) 
	{
		// grain charging rates have dimension cm-3 s-1
		phel_rate_cr[i] = phel_rate_uv[i] = phel_rate_vis[i] = ion_neutr_rate[i] = el_att_rate[i] = 0.;
		d = 0.5*conc_e *sqrt(8.*temp_e_erg/(M_PI*ELECTRON_MASS))*(1. - exp(-1.e+7*dust->get_grain_radius(i)))
			*dust->get_grain_area(i);
		
		// Two quantities are calculated in this case: grain concentration and average charge,
		// Note: there is no restriction on grain charge value increase (zmin, zmax)
		if (is_dust_charge_large[i])
		{	
			// electron attachment:
			rate = d *accr_func->get_accr_rate(dust->get_grain_radius(i)*temp_e, -y_data[nb_dch[i] + 1]);
			ydot_data[nb_dch[i] + 1] -= rate;
			el_att_rate[i] = rate *y_data[nb_dch[i]];
				
			// ion attachment;
			for (l = network->nb_of_reactions - network->nb_reactions_ion_grains; l < network->nb_of_reactions; l++) 
            {
                const chem_reaction & reaction = network->reaction_array[l];

				// sticking coefficient is equal 1, first parameter - branching ratio, second - ion velocity const:
				c = reaction.parameters[0]
					*y_data[reaction.reactant[0]] *dust->get_grain_area(i);

				a = temp_i + PI_DIVBY_EIGHT_BOLTZMANN_CONST *grain_veli2[nb_dch[i]-nb_dch[0]] * reaction.mass1;

				rate = c *sqrt(a)*accr_func->get_accr_rate(dust->get_grain_radius(i)*a, y_data[nb_dch[i] + 1]);
				ydot_data[nb_dch[i] + 1] += rate;
				rate *= y_data[nb_dch[i]];		
					
				// source terms for ions and molecules (momentum, energy)
				b = rate * reaction.mass1;

				mom_i -= b*vel_i;
				mom_n += b*grain_velz[nb_dch[i]-nb_dch[0]];
				
				en_n += 0.5*b *grain_veln2[nb_dch[i]-nb_dch[0]];
				en_i -= 1.5*temp_i_erg*rate;

				ydot_data[reaction.reactant[0]] -= rate;
				for (j = 0; j < reaction.nb_of_products; j++) {
					ydot_data[reaction.product[j]] += rate;
				}
				
				chem_reaction_rates[l] += rate;
				ion_neutr_rate[i] += rate;
			}

			// photoelectron emission:
			a = photoem_factor_is_uv *dust->components[i]->get_isrf_uv_phel_rate(y_data[nb_dch[i] + 1], b);
			ydot_data[nb_dch[i] + 1] += a;

			a *= y_data[nb_dch[i]];
			phel_rate_uv[i] = a;
			
			rate = a;              // total rate of photoemission in cm-3 s-1
			pheff_gas_heat += a*b; // heating rate, erg cm-3 s-1

			a = photoem_factor_is_vis *dust->components[i]->get_isrf_vis_phel_rate(y_data[nb_dch[i] + 1], b);
			ydot_data[nb_dch[i] + 1] += a;

			a *= y_data[nb_dch[i]];
			phel_rate_vis[i] = a;

			rate += a; 
			pheff_gas_heat += a*b; 

			a = photoem_factor_cr *dust->components[i]->get_cr_phel_rate(y_data[nb_dch[i] + 1], b); 
			ydot_data[nb_dch[i] + 1] += a;	
			
			a *= y_data[nb_dch[i]];
			phel_rate_cr[i] = a;
			
			rate += a;
			pheff_gas_heat += a*b;

			// source terms for electron (momentum, energy):
			mom_e += rate*ELECTRON_MASS *grain_velz[nb_dch[i]-nb_dch[0]];	
		}
		else // grain distribution is calculated in this case:
		{ // electron attachment:
			j = nb_dch[i] - nb_dch[0];
			for (k = 0; k < max_grain_charge[i] - min_grain_charge[i]; k++)
			{
				rate = d *y_data[nb_dch[i] + k + 1] 
					*accr_func->get_accr_rate(dust->get_grain_radius(i)*temp_e, -(min_grain_charge[i] + k + 1));

				ydot_data[nb_dch[i] + k] += rate;
				ydot_data[nb_dch[i] + k + 1] -= rate;
				el_att_rate[i] += rate;

				// calculation of source terms due to grain charge fluctuations (PAH momentum and energy change):
				if (is_charged_dust_as_ions[i]) 
				{
					a = rate *dust->get_grain_mass(i);
					if (k + 1 == 1 - min_grain_charge[i]) // charge of the grain that accepts electron is +1
					{		
						mom_n += a*vel_i;
						mom_i -= a*vel_i;
						en_n += 0.5*a*(vel_n - vel_i)*(vel_n - vel_i);
					}
					else if (k + 1 == -min_grain_charge[i]) // grain charge is 0
					{
						mom_n -= a*vel_n;
						mom_i += a*vel_n;
						en_i += 0.5*a*(vel_n - vel_i)*(vel_n - vel_i);
					}
				}
			}
			// ion attachment on grains; 
			for (l = network->nb_of_reactions - network->nb_reactions_ion_grains; l < network->nb_of_reactions; l++)
			{
                const chem_reaction & reaction = network->reaction_array[l];
				c = reaction.parameters[0] *y_data[reaction.reactant[0]] *dust->get_grain_area(i);

				rate = 0.;
				for (k = 0; k < max_grain_charge[i] - min_grain_charge[i]; k++)
				{
					a = temp_i + PI_DIVBY_EIGHT_BOLTZMANN_CONST *grain_veli2[j + k] * reaction.mass1; // temperature in K;

					b = c *sqrt(a)*y_data[nb_dch[i] + k]  
						*accr_func->get_accr_rate(dust->get_grain_radius(i)*a, min_grain_charge[i] + k);

					ydot_data[nb_dch[i] + k + 1] += b;
					ydot_data[nb_dch[i] + k] -= b;
					rate += b;

					// source terms: molecule momentum and energy:
					mom_n += b * reaction.mass1 * grain_velz[j + k]; // velocity of the grain with Z
					en_n += 0.5 * b *grain_veln2[j + k] * reaction.mass1;
						
					// PAH momentum and energy:
					if (is_charged_dust_as_ions[i]) 
					{
						a = b *dust->get_grain_mass(i);
						if (k == -1 - min_grain_charge[i]) // charge of the grain that accepts ion is -1
						{		
							mom_n += a*vel_i;
							mom_i -= a*vel_i;
							en_n += 0.5*a*(vel_n - vel_i)*(vel_n - vel_i);
						}
						else if (k == -min_grain_charge[i]) // grain charge is 0
						{
							mom_n -= a*vel_n;
							mom_i += a*vel_n;
							en_i += 0.5*a*(vel_n - vel_i)*(vel_n - vel_i);
						}
					}
				}
				// source terms: ion momentum and energy;
				mom_i -= rate * reaction.mass1 *vel_i;
				en_i -= 1.5*temp_i_erg*rate;

				ydot_data[reaction.reactant[0]] -= rate;
				for (k = 0; k < reaction.nb_of_products; k++) {
					ydot_data[reaction.product[k]] += rate;
				}
				
				chem_reaction_rates[l] += rate;
				ion_neutr_rate[i] += rate;
			}
			// photoelectron emission:
			for (k = 0; k < max_grain_charge[i] - min_grain_charge[i]; k++)
			{
				a = photoem_factor_is_uv *y_data[nb_dch[i] + k]
					*dust->components[i]->get_isrf_uv_phel_rate(min_grain_charge[i] + k, b); 

				rate = a;
				phel_rate_uv[i] += a;
				pheff_gas_heat += a*b;

				a = photoem_factor_is_vis *y_data[nb_dch[i] + k]
					*dust->components[i]->get_isrf_vis_phel_rate(min_grain_charge[i] + k, b); 

				rate += a;
				phel_rate_vis[i] += a;
				pheff_gas_heat += a*b;

				a = photoem_factor_cr *y_data[nb_dch[i] + k]
					*dust->components[i]->get_cr_phel_rate(min_grain_charge[i] + k, b); 
				
				rate += a;
				phel_rate_cr[i] += a;
				pheff_gas_heat += a*b;

				ydot_data[nb_dch[i] + k] -= rate;
				ydot_data[nb_dch[i] + k + 1] += rate;

				mom_e += rate*ELECTRON_MASS *grain_velz[j + k];

				if (is_charged_dust_as_ions[i]) 
				{
					a = rate *dust->get_grain_mass(i);
					if (k == -1 - min_grain_charge[i]) // grain charge is -1
					{		
						mom_n += a*vel_i;
						mom_i -= a*vel_i;
						en_n += 0.5*a*(vel_n - vel_i)*(vel_n - vel_i);
					}
					else if (k == -min_grain_charge[i]) // grain charge is 0
					{
						mom_n -= a*vel_n;
						mom_i += a*vel_n;
						en_i += 0.5*a*(vel_n - vel_i)*(vel_n - vel_i);
					}
				}
			}
		}
		el_att += el_att_rate[i];
		phel_em += phel_rate_uv[i] + phel_rate_vis[i] + phel_rate_cr[i];
	}
	ydot_data[network->e_nb] += phel_em - el_att;

	mom_gain_n += mom_n;
	mom_gain_i += mom_i;
	mom_gain_e += mom_e - el_att *ELECTRON_MASS *vel_i;
	
	// energy gain due to photoeffect; 
	// at fractional ionization > 1.e-6 the excitation of H2 molecules diminishes and fast electron begin to lose energy on Coulomb scattering;
	// Dalgarno et al. ApJSS 125, p. 237 (1999); critical densities for v = 1 states of H2 are about 10^5 - 10^6 cm-3;
	// For simplicity, all photoeffect electron energy goes to neutral gas here:
	energy_gain_n += en_n + pheff_gas_heat;
	energy_gain_e += 1.5*(phel_em - el_att)*temp_e_erg;	
	energy_gain_i += en_i;
	
	// calculation of number density and mass density gains
	nb_gain_n = nb_gain_i = nb_gain_e = mass_gain_n = mass_gain_i = mass_gain_e = 0.;
	for (i = 0; i < nb_of_gas_species; i++) {
		if (network->species[i].type == "neutral")
		{
			nb_gain_n += ydot_data[i];
			mass_gain_n += ydot_data[i] *network->species[i].mass;
		}
		else if (network->species[i].type == "ion")
		{
			nb_gain_i += ydot_data[i];
			mass_gain_i += ydot_data[i] *network->species[i].mass;
		}
		else {
			nb_gain_e += ydot_data[i];
			mass_gain_e += ydot_data[i] *ELECTRON_MASS;
		}
	}

	// calculation of source terms due to grain charge fluctuations (PAH nb and mass density change):
	for (i = 0; i < nb_of_dust_comp; i++) {
		if (is_charged_dust_as_ions[i] && !is_dust_charge_large[i])
		{
			a = ydot_data[nb_dch[i] - min_grain_charge[i]];
			nb_gain_n += a;
			nb_gain_i -= a;
			
			a *= dust->get_grain_mass(i);
			mass_gain_n += a;
			mass_gain_i -= a;
		}
	}

	// production rates of species, for which level populations are calculated:
    h2_prod = ydot_data[network->h2_nb]; // not accounting for e- excitation and H2-H dissociation
	h2o_prod = ydot_data[network->h2o_nb];
	co_prod = ydot_data[network->co_nb];
	oh_prod = ydot_data[network->oh_nb];
	nh3_prod = ydot_data[network->nh3_nb];
	ch3oh_prod = ydot_data[network->ch3oh_nb];
	ci_prod = ydot_data[network->ci_nb];
	oi_prod = ydot_data[network->oi_nb];
	cii_prod = ydot_data[network->cii_nb];

	memset(dust_heat_h2_line, 0, nb_of_dust_comp*sizeof(double));
	memset(dust_heat_mline, 0, nb_of_dust_comp*sizeof(double));
	
#if (SAVE_RADIATIVE_TRANSFER_FACTORS)
	memset(gamma_factors, 0, nb_rtd*sizeof(double));
	memset(delta_factors, 0, nb_rtd*sizeof(double));
	memset(dheat_efficiency, 0, nb_rtd*sizeof(double));
	memset(esc_prob_int1, 0, nb_rtd*sizeof(double));
	memset(esc_prob_int2, 0, nb_rtd*sizeof(double));
#endif

	// calculation of the heating (cooling) due to H2: 
	h2_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);
    h2_coll->set_ion_param(temp_n, temp_i, conc_hp, conc_h3p, coll_partn_conc, indices);

	specimen_population_derivative(y_data, ydot_data, nb_of_species, h2_di, h2_einst, h2_coll, coll_partn_conc, indices, vel_n_grad,
		neut_heat_h2, el_heat_h2, ion_heat_h2, rad_energy_loss_h2, dust_heat_h2_line);

#if (SAVE_RADIATIVE_TRANSFER_FACTORS)
    calc_radiative_coeff(y_data, nb_of_species, h2_di, h2_einst, h2_coll, coll_partn_conc, indices, vel_n_grad, gamma_factors, 
        delta_factors, dheat_efficiency, esc_prob_int1, esc_prob_int2, 0);
#endif

	// here the excitation of H2 by cosmic rays is taken into acccount;
	if (H2_CR_EXCITATION) {
		ioniz_fraction = conc_e/conc_h_tot;
		// the values of two last parameters are returned (index and derivative for ionization grid)
		h2_excit_cr->set_ionization(ioniz_fraction, l, b);

		for (i = 0; i < h2_excit_cr->get_nb_lev(); i++) 
		{
			a = y_data[nb_of_species + i] *STANDARD_CR_IONIZ_RATE *cr_ioniz_rate;
			ydot_data[nb_of_species + i] -= a *h2_excit_cr->get_efficiency(i, l, b);

			for (k = 0; k < nb_lev_h2; k++) {
				ydot_data[nb_of_species + k] += a *h2_excit_cr->get_efficiency(i, k, l, b);
			}
		}
	}

    if (H2_H_DISSIOCIATION_DATA) {
        h2_h_diss_cooling = h2_h_diss_rate = 0.;
        for (i = nb_of_species; i < nb_of_species + nb_lev_h2; i++)
        {
            c = y_data[i] * conc_h *h2_h_diss_data->get_rate(i - nb_of_species, temp_n);
            h2_h_diss_rate += c;
            ydot_data[i] -= c;

            h2_h_diss_cooling += (network->reaction_array[i].energy_released 
                + h2_di->lev_array[i - nb_of_species].energy * CM_INVERSE_TO_ERG) *c; // < 0 for cooling,
        }
        ydot_data[network->h_nb] += 2.*h2_h_diss_rate;
        ydot_data[network->h2_nb] -= h2_h_diss_rate;

        i = network->h2_h_diss_nb;
        chem_reaction_rates[i] = h2_h_diss_rate; // reaction rate is saved,

        chem_heating_rates_n[i] = h2_h_diss_cooling; // cooling of the gas by this process,
        neut_heat_chem += h2_h_diss_cooling;
        energy_gain_n += h2_h_diss_cooling;
    }

    oh2_form_hcoll = 0.;
    for (i = 0; i < nb_lev_h2; i++)
    {
        if (rounding(h2_di->lev_array[i].spin) == 1) {
            oh2_form_hcoll += ydot_data[nb_of_species + i];
        }
    }

	// other processes of H2 formation are assumed not to change the level population: molecules are formed
	// with level population distrubution proportional to the current level population distribution;
	c = h2_prod/conc_h2;
	for (i = nb_of_species; i < nb_of_species + nb_lev_h2; i++) {
		ydot_data[i] += y_data[i] *c;
	}

	energy_gain_n += neut_heat_h2;
	energy_gain_e += el_heat_h2;
    energy_gain_i += ion_heat_h2;
    nb = nb_of_species + nb_lev_h2;
	
	// calculation of the level population gain for para-H2O molecule (not normalized to H2O concentration);
	ph2o_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);
	
	specimen_population_derivative(y_data, ydot_data, nb, ph2o_di, ph2o_einst, ph2o_coll, coll_partn_conc, indices, vel_n_grad, 
		neut_heat_ph2o, el_heat_ph2o, a, a, dust_heat_mline);

#if (SAVE_RADIATIVE_TRANSFER_FACTORS)
    nb2 = nb_lev_h2*(nb_lev_h2-1)/2;
    calc_radiative_coeff(y_data, nb, ph2o_di, ph2o_einst, ph2o_coll, coll_partn_conc, indices, vel_n_grad, gamma_factors, 
        delta_factors, dheat_efficiency, esc_prob_int1, esc_prob_int2, nb2);
#endif

	// processes of H2O formation are assumed to not change the level population, concentration of H2O must be > 0;
	c = h2o_prod/conc_h2o;
	for (i = 0; i < nb_lev_h2o; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c;
	}

	energy_gain_n += neut_heat_ph2o;
	energy_gain_e += el_heat_ph2o;
	nb += nb_lev_h2o;
	
	// calculation of the level population gain for ortho-H2O molecule
	oh2o_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);
	
	specimen_population_derivative(y_data, ydot_data, nb, oh2o_di, oh2o_einst, oh2o_coll, coll_partn_conc, indices, vel_n_grad, 
		neut_heat_oh2o, el_heat_oh2o, a, a, dust_heat_mline);

#if (SAVE_RADIATIVE_TRANSFER_FACTORS)
    nb2 += nb_lev_h2o * (nb_lev_h2o - 1) / 2;
    calc_radiative_coeff(y_data, nb, oh2o_di, oh2o_einst, oh2o_coll, coll_partn_conc, indices, vel_n_grad, gamma_factors, 
        delta_factors, dheat_efficiency, esc_prob_int1, esc_prob_int2, nb2);
#endif

	for (i = 0; i < nb_lev_h2o; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c; // c is defined higher
	}

	energy_gain_n += neut_heat_oh2o;
	energy_gain_e += el_heat_oh2o;
	nb += nb_lev_h2o;

	// calculation of the level population gain for CO molecule
	co_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);

	specimen_population_derivative(y_data, ydot_data, nb, co_di, co_einst, co_coll, coll_partn_conc, indices, vel_n_grad, 
		neut_heat_co, a, a, a, dust_heat_mline);

#if (SAVE_RADIATIVE_TRANSFER_FACTORS)
    nb2 += nb_lev_h2o*(nb_lev_h2o-1)/2;
    calc_radiative_coeff(y_data, nb, co_di, co_einst, co_coll, coll_partn_conc, indices, vel_n_grad, gamma_factors, 
        delta_factors, dheat_efficiency, esc_prob_int1, esc_prob_int2, nb2);
#endif

	c = co_prod/conc_co;
	for (i = 0; i < nb_lev_co; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c;
	}
	energy_gain_n += neut_heat_co;
	nb += nb_lev_co;
	
	// calculation of the level population gain for OH molecule	
    c = oh_prod/conc_oh;

#if (CALCULATE_POPUL_NH3_OH)
	oh_coll->set_gas_param(temp_n, temp_e, conc_he, conc_h2j0, conc_h2-conc_h2j0, conc_h, conc_e, coll_partn_conc, indices);

	specimen_population_derivative(y_data, ydot_data, nb, oh_di, oh_einst, oh_coll, coll_partn_conc, indices, vel_n_grad, 
		neut_heat_oh, a, a, a, dust_heat_mline);

    energy_gain_n += neut_heat_oh;
#endif

	for (i = 0; i < nb_lev_oh; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c;
	}
	nb += nb_lev_oh;
	
	// calculation of the level population gain for para-NH3 molecule;
    c = nh3_prod/conc_nh3;

#if (CALCULATE_POPUL_NH3_OH)
	pnh3_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);
	
	specimen_population_derivative(y_data, ydot_data, nb, pnh3_di, pnh3_einst, pnh3_coll, coll_partn_conc, indices, vel_n_grad, 
		neut_heat_pnh3, a, a, a, dust_heat_mline);
	
    energy_gain_n += neut_heat_pnh3;
#endif	
	
    for (i = 0; i < nb_lev_pnh3; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c;
	}
	nb += nb_lev_pnh3;
	
	// calculation of the level population gain for ortho-NH3 molecule
#if (CALCULATE_POPUL_NH3_OH)
    onh3_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);
	
	specimen_population_derivative(y_data, ydot_data, nb, onh3_di, onh3_einst, onh3_coll, coll_partn_conc, indices, vel_n_grad, 
		neut_heat_onh3, a, a, a, dust_heat_mline);

    energy_gain_n += neut_heat_onh3;
#endif
	for (i = 0; i < nb_lev_onh3; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c; // c is defined higher
	}
	nb += nb_lev_onh3;
	
	// calculation of the level population gain for CH3OH A and E
	c = ch3oh_prod/conc_ch3oh;

#if (CALCULATE_POPUL_METHANOL)
	ch3oh_a_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);
	
	specimen_population_derivative(y_data, ydot_data, nb, ch3oh_a_di, ch3oh_a_einst, ch3oh_a_coll, coll_partn_conc, indices, vel_n_grad, 
		neut_heat_ch3oh_a, a, a, a, dust_heat_mline);
	
	energy_gain_n += neut_heat_ch3oh_a;
#endif

	for (i = 0; i < nb_lev_ch3oh; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c; // c is defined higher
	}
	nb += nb_lev_ch3oh;
	
#if (CALCULATE_POPUL_METHANOL)
	ch3oh_e_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);
	
	specimen_population_derivative(y_data, ydot_data, nb, ch3oh_e_di, ch3oh_e_einst, ch3oh_e_coll, coll_partn_conc, indices, vel_n_grad, 
		neut_heat_ch3oh_e, a, a, a, dust_heat_mline);
	
	energy_gain_n += neut_heat_ch3oh_e;
#endif

	for (i = 0; i < nb_lev_ch3oh; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c; // c is defined higher
	}
	nb += nb_lev_ch3oh;
	
	// calculation of the level population gain for CI:
	CI_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);
	
	specimen_population_derivative(y_data, ydot_data, nb, CI_di, CI_einst, CI_coll, coll_partn_conc, indices, vel_n_grad, 
		en_n, en_e, a, a, dust_heat_mline);

	c = ci_prod/conc_ci;
	for (i = 0; i < nb_lev_ci; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c;
	}

	neut_heat_atoms = en_n;
	el_heat_atoms = en_e;
	nb += nb_lev_ci;
	
	// calculation of the level population gain for OI:
	OI_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);

	specimen_population_derivative(y_data, ydot_data, nb, OI_di, OI_einst, OI_coll, coll_partn_conc, indices, vel_n_grad, 
		en_n, en_e, a, a, dust_heat_mline);

	c = oi_prod/conc_oi;
	for (i = 0; i < nb_lev_oi; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c;
	}

	neut_heat_atoms += en_n;
	el_heat_atoms += en_e;
	nb += nb_lev_oi;
	
	// calculation of the level population gain for CII:
	CII_coll->set_gas_param(temp_n, temp_e, conc_he, conc_ph2, conc_oh2, conc_h, conc_e, coll_partn_conc, indices);

	specimen_population_derivative(y_data, ydot_data, nb, CII_di, CII_einst, CII_coll, coll_partn_conc, indices, vel_i_grad, 
		en_n, en_e, a, a, dust_heat_mline);

	c = cii_prod/conc_cii;
	for (i = 0; i < nb_lev_cii; i++) {
		ydot_data[nb + i] += y_data[nb + i]*c;
	}

	neut_heat_atoms += en_n;
	el_heat_atoms += en_e;

	energy_gain_n += neut_heat_atoms;
	energy_gain_e += el_heat_atoms;

	// Heating of the gas by cosmic rays - secondary electrons, Dalgarno et al., ApJSS 125, 237 (1999);
	// It is assumed that cosmic ray ionization rate is in s-1 per H_2 (it is accounted by factor 0.5),
	// Production rate of electrons is approximately equal to production rate of H2+ (per H2 molecule), Padovani et al., A&A 501, 619, 2009;
	neut_cr_heat = cr_heating_const*cr_ioniz_rate *0.5*conc_h_tot; 
	energy_gain_n += neut_cr_heat;

	// calculation of source terms due to elastic scattering, neutral-ion scattering:
	if (fabs(vel_ni_diff) > DBL_EPSILON)
	{
		en_n = en_i = 0.;
		elastic_h2_ions->calc_source_terms(mom_gain_n, mom_gain_i, en_n, en_i, conc_h2, conc_i, vel_ni_diff, temp_n, temp_i);
		elastic_h_ions->calc_source_terms(mom_gain_n, mom_gain_i, en_n, en_i, conc_h, conc_i, vel_ni_diff, temp_n, temp_i);
		elastic_he_ions->calc_source_terms(mom_gain_n, mom_gain_i, en_n, en_i, conc_he, conc_i, vel_ni_diff, temp_n, temp_i);

		energy_gain_n += en_n;
		energy_gain_i += en_i;

		neut_heat_scatt_ions = en_n;
		ion_heat_scatt_n = en_i;
		
		// neutral - electron scattering:
		en_n = en_e = 0.;
		elastic_h2_el->calc_source_terms(mom_gain_n, mom_gain_e, en_n, en_e, conc_h2, conc_e, vel_ni_diff, temp_n, temp_e);
		elastic_he_el->calc_source_terms(mom_gain_n, mom_gain_e, en_n, en_e, conc_he, conc_e, vel_ni_diff, temp_n, temp_e);
		elastic_h_el->calc_source_terms(mom_gain_n, mom_gain_e, en_n, en_e, conc_h, conc_e, vel_ni_diff, temp_n, temp_e);

		energy_gain_n += en_n;
		energy_gain_e += en_e;
		
		el_heat_scatt_neut = en_e;
		neut_heat_scatt_el = en_n;
	}
	else
	{
		en_n = en_i = 0.;
		elastic_h2_ions->calc_source_terms(en_n, en_i, conc_h2, conc_i, temp_n, temp_i);
		elastic_h_ions->calc_source_terms(en_n, en_i, conc_h, conc_i, temp_n, temp_i);
		elastic_he_ions->calc_source_terms(en_n, en_i, conc_he, conc_i, temp_n, temp_i);

		energy_gain_n += en_n;
		energy_gain_i += en_i;

		neut_heat_scatt_ions = en_n;
		ion_heat_scatt_n = en_i;

		// neutral - electron scattering:
		en_n = en_e = 0.;
		elastic_h2_el->calc_source_terms(en_n, en_e, conc_h2, conc_e, temp_n, temp_e);
		elastic_he_el->calc_source_terms(en_n, en_e, conc_he, conc_e, temp_n, temp_e);
		elastic_h_el->calc_source_terms(en_n, en_e, conc_h, conc_e, temp_n, temp_e);

		energy_gain_n += en_n;
		energy_gain_e += en_e;

		el_heat_scatt_neut = en_e;
		neut_heat_scatt_el = en_n;
	}

	// energy transfer between ions and electrons:
	en_i = ion_electron_energy_constant *conc_e*conc_i*(temp_e - temp_i) /(ion_mass*pow(temp_e, 1.5)) 
		*(log_coulomb_constant + 0.5*log(temp_e *temp_e *temp_e/conc_e));

	ion_heat_scatt_el = en_i;
	el_heat_scatt_ions = -en_i;

	energy_gain_e -= en_i;
	energy_gain_i += en_i;

	// Equations for dust temperature, heating rates are per one grain, erg s-1,
	neut_heat_dust_coll = 0.;
	c = 2.*conc_n *sqrt(8.*temp_n_erg/(M_PI*av_mol_mass))*BOLTZMANN_CONSTANT;
	energy_gain_d /= ads_dust_area; // area is n_gr*pi*a*a;

	for (i = 0; i < nb_of_dust_comp; i++) 
	{
		// heating by diluted star-light and cooling by self-emission:
		temp_d = y_data[nb_dct + i];
		ydot_data[nb_dct + i] = dh_isrf_arr[i] - dust->components[i]->get_int_emiss(temp_d);
		
		// heating by molecular emission:
		if (GRAIN_HEATING_BY_MOLECULES_ON) {
			ydot_data[nb_dct + i] += (dust_heat_mline[i] + dust_heat_h2_line[i])*CM_INVERSE_TO_ERG;
		}

		// the heat associated with adsorption-desorption (only binding energy), grain surface chemical reactions may be not included (check);
		// normalized on one dust particle:
		if (dust->get_grain_radius(i) > MIN_ADSORPTION_RADIUS ) 
		{
			dust_heat_chem[i] = energy_gain_d *dust->get_grain_area(i);
			if (GRAIN_HEATING_GS_CHEMISTRY)
				ydot_data[nb_dct + i] += dust_heat_chem[i];
		}

		// heat due to gas-dust collisions (eq. from Draine, ApJ 241, p.1021, 1980):
		if (is_dust_charge_large[i])
		{
			b = 0.25*av_mol_mass *grain_veln2[nb_dch[i]-nb_dch[0]]/temp_n_erg;
			// dependence of accomodation factor on gas temperature is taken to be linear on log scale, 
			// and not depending on dust temperature, grain type and size, see Burke & Hollenbach, ApJ 265, p. 223 (1983);  
			a = 0.35 - 0.0833*(log10(temp_n*(1. + M_PI_2*b)) - 1.); // 0.35 at T = 10 K; 0.1 at T = 10^4 K;
			if (a < 0.1)
				a = 0.1;

			dust_heat_coll[i] = a*c*sqrt(1. + M_PI_2*b)*dust->get_grain_area(i)*(temp_n - temp_d + b*temp_n);
		}
		else
		{
			dust_heat_coll[i] = 0.;
			for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) 
			{
				b = 0.25*av_mol_mass *grain_veln2[nb_dch[i] - nb_dch[0] + k]/temp_n_erg;
				a = 0.35 - 0.0833*(log10(temp_n*(1. + M_PI_2*b)) - 1.); 
				if (a < 0.1)
					a = 0.1;

				dust_heat_coll[i] += a*c*sqrt(1. + M_PI_2*b)*dust->get_grain_area(i)*(temp_n - temp_d + b*temp_n) 
					*y_data[nb_dch[i] + k];
			}
			dust_heat_coll[i] /= grain_conc[i];
		}
		ydot_data[nb_dct + i] += dust_heat_coll[i];

		// Debye model of heat capacity of grains,
		// Draine & Li, ApJ 551, p.807, 2001; Cuppen et al., MNRAS 367, p.1757, 2006;
		ydot_data[nb_dct + i] /= DUST_HEAT_CAPACITY_FACTOR*dust->components[i]->heat_capacity_const *temp_d*temp_d*temp_d;

		// Note: the calculation of this parameter must not be included in chemical source terms:
		neut_heat_dust_coll -= dust_heat_coll[i] *grain_conc[i];
	}
	energy_gain_n += neut_heat_dust_coll;
	return 0;
}

// the dimension of heating efficiency dh_eff[] is cm-1 cm-3 s-1, dimension of grain heating rate heating_d[] is cm-1 s-1
void evolution_data::specimen_population_derivative(const realtype *y_data, realtype *ydot_data, int nb, const energy_diagram *mol_di,
    const einstein_coeff *mol_einst, const collisional_transitions *mol_coll, double *coll_partn_conc, int *indices,
    double vel_grad, double & heating_n, double & heating_e, double & heating_i, double & rad_energy_loss, double *heating_d)
{
    int i, l, k, nb_lev;
    double upl, lowl, ep1, ep2, delta, gamma, c, d, vd, down_rate, up_rate, energy, line_em, line_op, dust_op, h_n, h_e, h_i, r_e_l;

    nb_lev = mol_di->nb_lev;
    vd = pow(2.*BOLTZMANN_CONSTANT*temp_n / mol_di->mol.mass + vel_turb * vel_turb, 0.5);
    r_e_l = h_e = h_n = h_i = 0.;

#pragma omp parallel reduction(+: h_n, h_e, h_i, r_e_l) private(i, l, k, upl, lowl, c, d, delta, gamma, ep1, ep2, down_rate, up_rate, energy, line_em, line_op, dust_op)
    {
        double *arr = new double[nb_lev];
        memset(arr, 0, nb_lev * sizeof(double));

        double *h_d = new double[nb_of_dust_comp];
        memset(h_d, 0, nb_of_dust_comp * sizeof(double));

#pragma omp for schedule(dynamic, 1)
        for (l = 0; l < nb_lev - 1; l++)
        {
            lowl = y_data[nb + l];
            for (i = l + 1; i < nb_lev; i++)
            {
                // level energy is in cm-1;
                energy = mol_di->lev_array[i].energy - mol_di->lev_array[l].energy;
                upl = y_data[nb + i];

                mol_coll->get_rate_neutrals(mol_di->lev_array[i], mol_di->lev_array[l], down_rate, up_rate, temp_n,
                    coll_partn_conc, indices);

                down_rate *= upl;
                up_rate *= lowl;

                h_n += (down_rate - up_rate) *energy; // level energy is in cm-1;	
 
                c = down_rate;
                d = up_rate;

                mol_coll->get_rate_electrons(mol_di->lev_array[i], mol_di->lev_array[l], down_rate, up_rate, temp_e,
                    coll_partn_conc, indices);

                down_rate *= upl;
                up_rate *= lowl;

                h_e += (down_rate - up_rate) *energy;

#if (H2_IONS_EXCITATION)
                c += down_rate;
                d += up_rate;
                
                mol_coll->get_rate_ions(mol_di->lev_array[i], mol_di->lev_array[l], down_rate, up_rate, temp_n, temp_i,
                    coll_partn_conc, indices);

                down_rate *= upl;
                up_rate *= lowl;

                h_i += (down_rate - up_rate) *energy;
#endif

                down_rate += c;
                up_rate += d;

                arr[i] += up_rate - down_rate;
                arr[l] += down_rate - up_rate;

                if (mol_einst->arr[i][l] > 1.e-99)
                {
                    c = 1. / (EIGHT_PI *vd*energy*energy*energy);

                    line_em = mol_einst->arr[i][l] * upl *c;
                    line_op = mol_einst->arr[l][i] * lowl *c - line_em + MIN_LINE_OPACITY;
                    dust_op = dust->absorption(energy, grain_conc); // dust absorption in cm-1

                    if (line_op < 0.)
                        line_op *= -0.1;

                    if (line_em < 0.)
                        line_em *= -0.1; // level population may be negative due to simulation errors

                    gamma = fabs(vel_grad) / (vd*line_op);
                    delta = fabs(vel_grad) / (vd*dust_op);

                    // linear interpolation is used in calculating escape probabilities:
                    ep1 = loss_func_line_phot->get_esc_func_2(gamma, delta);

                    c = line_em / line_op * ep1;
                    d = mol_einst->arr[i][l] * upl *(1. + c);

                    arr[i] -= d;
                    arr[l] += d;

                    d = mol_einst->arr[l][i] * lowl *c;

                    arr[i] += d;
                    arr[l] -= d;

                    // dust heating:
                    ep2 = loss_func_cont_phot->get_esc_func_2(gamma, delta);
                    c = mol_einst->arr[i][l] * upl*energy*ep2 / dust_op; // cm-1 s-1 cm-2

                    for (k = 0; k < nb_of_dust_comp; k++) { // absorption in cm2, heating rate of one grain in cm-1 s-1 
                        h_d[k] += c * dust->components[k]->absorption(energy);
                    }

                    r_e_l += mol_einst->arr[i][l] * upl *energy; // radiative energy loss, cm-1 s-1 cm-3
                }
            }
        }
#pragma omp critical
        {
            for (i = 0; i < nb_lev; i++) {
                ydot_data[nb + i] += arr[i];
            }
            for (i = 0; i < nb_of_dust_comp; i++) {
                heating_d[i] += h_d[i];
            }
        }
        delete[] arr;
        delete[] h_d;
    }
    heating_n = h_n * CM_INVERSE_TO_ERG; // heating/cooling units are erg cm-3 s-1;
    heating_e = h_e * CM_INVERSE_TO_ERG;
    heating_i = h_i * CM_INVERSE_TO_ERG;
    rad_energy_loss = r_e_l * CM_INVERSE_TO_ERG;
}


void evolution_data::calc_radiative_coeff(const realtype *y_data, int nb, const energy_diagram *mol_di, const einstein_coeff *mol_einst, 
    const collisional_transitions *mol_coll, double *coll_partn_conc, int *indices, double vel_grad, double *g_factors, 
	double *d_factors, double *dh_eff, double *ep_int1, double *ep_int2, int nb2)
{
	int i, l, nb_lev;
	double upl, lowl, ep1, ep2, delta, gamma, c, vd, energy, line_em, line_op, dust_op;

	nb_lev = mol_di->nb_lev;
	vd = pow(2.*BOLTZMANN_CONSTANT*temp_n /mol_di->mol.mass + vel_turb*vel_turb, 0.5);
	
#pragma omp parallel private(i, l, upl, lowl, c, delta, gamma, ep1, ep2, energy, line_em, line_op, dust_op)
	{
		double **g_arr = alloc_2d_array<double>(nb_lev, nb_lev);
		memset(*g_arr, 0, nb_lev*nb_lev*sizeof(double));

		double **d_arr = alloc_2d_array<double>(nb_lev, nb_lev);
		memset(*d_arr, 0, nb_lev*nb_lev*sizeof(double));

		double **eff_arr = alloc_2d_array<double>(nb_lev, nb_lev);
		memset(*eff_arr, 0, nb_lev*nb_lev*sizeof(double));

		double **e1_arr = alloc_2d_array<double>(nb_lev, nb_lev);
		memset(*e1_arr, 0, nb_lev*nb_lev*sizeof(double));

		double **e2_arr = alloc_2d_array<double>(nb_lev, nb_lev);
		memset(*e2_arr, 0, nb_lev*nb_lev*sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (l = 0; l < nb_lev - 1; l++) {
            lowl = y_data[nb + l];
			for (i = l + 1; i < nb_lev; i++) {					
				if (mol_einst->arr[i][l] > 1.e-99)
				{        
                    upl = y_data[nb + i];                    
                    energy = mol_di->lev_array[i].energy - mol_di->lev_array[l].energy; // level energy is in cm-1;
					c = 1./(EIGHT_PI *vd*energy*energy*energy);
				
					line_em = mol_einst->arr[i][l] *upl *c;
					line_op = mol_einst->arr[l][i] *lowl *c - line_em + MIN_LINE_OPACITY;				
					dust_op = dust->absorption(energy, grain_conc); // dust absorption in cm-1
				
					if (line_op < 0.) 
						line_op *= -0.1;

					if (line_em < 0.)
						line_em *= -0.1; // level population may be negative due to simulation errors
					
					gamma = fabs(vel_grad) /(vd*line_op);
					delta = fabs(vel_grad) /(vd*dust_op);

					// linear interpolation is used in calculating escape probabilities:
					ep1 = loss_func_line_phot->get_esc_func_2(gamma, delta);
					
					// dust heating:
					ep2 = loss_func_cont_phot->get_esc_func_2(gamma, delta);
					c = mol_einst->arr[i][l]*upl*energy*ep2/dust_op; // cm-1 s-1 cm-2
					
		            g_arr[l][i] = gamma;
					d_arr[l][i] = delta;
					eff_arr[l][i] = c*dust_op; // opacity in cm-1, heating rate in cm-1 cm-3 s-1

					e1_arr[l][i] = ep1;
					e2_arr[l][i] = ep2;
				}
			}
		}
#pragma omp critical
		{
			for (l = 0; l < nb_lev-1; l++) {
				for (i = l+1; i < nb_lev; i++) 
				{
					g_factors[nb2 + i*(i-1)/2 + l] += g_arr[l][i];
					d_factors[nb2 + i*(i-1)/2 + l] += d_arr[l][i];
					dh_eff[nb2 + i*(i-1)/2 + l] += eff_arr[l][i];

					ep_int1[nb2 + i*(i-1)/2 + l] += e1_arr[l][i];
					ep_int2[nb2 + i*(i-1)/2 + l] += e2_arr[l][i];
				}
			}
		}
		free_2d_array(g_arr);
		free_2d_array(d_arr);
		free_2d_array(eff_arr);
		free_2d_array(e1_arr);
		free_2d_array(e2_arr);
	}
}

double evolution_data::calc_grain_velocity(const N_Vector &y, double charge, double area, double & a, double & b, double &wt2, double & velx) const
{
	double c(0.);
	if (nb_of_equat == nb_mhd + NB_MHD_EQUATIONS) 
	{// shock simulations;
		c = charge*magn_field/(area *conc_h_tot);
		a = v_const_2*c*c/temp_n_erg;
		b = v_const_1*vel_ni_diff*vel_ni_diff/temp_n_erg + 1.;

		wt2 = 0.5*(a - 1. + sqrt((a-1.)*(a-1.) + 4.*a*b))/b;
		velx = -vel_ni_diff *sqrt(wt2)/(1. + wt2);
		return vel_n - vel_ni_diff *wt2/(1. + wt2);
	}
	else { // static cloud;
		c = charge*magn_field_0/(area *conc_h_tot);	
		a = v_const_2 *c*c/temp_n_erg;
		b = 1.;

		wt2 = 0.5*(a - 1. + sqrt((a-1.)*(a-1.) + 4.*a));
		velx = 0.;
		return 0.;
	}
}

double evolution_data::calc_av_grain_velocity(const N_Vector &y, int i) const
{
	int k;
	double a, b, c, d, conc, vel;

	if (is_dust_charge_large[i]) {
		vel = calc_grain_velocity(y, NV_Ith_S(y, nb_dch[i] + 1), dust->get_grain_area(i), a, b, c, d);
	}
	else {
		conc = vel = 0.;
		for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) 
		{
			vel += NV_Ith_S(y, nb_dch[i] + k) *calc_grain_velocity(y, min_grain_charge[i] + k, dust->get_grain_area(i), a, b, c, d);
			conc += NV_Ith_S(y,nb_dch[i] + k);
		}
		vel /= conc;
	}
	return vel;
}

double evolution_data::calc_dust_kinetic_energy_flux(const N_Vector & y) const
{
    int i, k;
    double a, b, c, d, vel, flux;

    flux = 0.;
    for (i = 0; i < nb_of_dust_comp; i++) {
        if (is_dust_charge_large[i]) {
            vel = calc_grain_velocity(y, NV_Ith_S(y, nb_dch[i] + 1), dust->get_grain_area(i), a, b, c, d);
            flux += 0.5 * vel * vel * vel * NV_Ith_S(y, nb_dch[i]) * dust->components[i]->mass;
        }
        else {
            for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++)
            {
                vel = calc_grain_velocity(y, min_grain_charge[i] + k, dust->get_grain_area(i), a, b, c, d);
                flux += 0.5 * vel * vel * vel * NV_Ith_S(y, nb_dch[i] + k) * dust->components[i]->mass;
            }
        }
    }
    return flux;
}

double evolution_data::calc_conc_ph2(const N_Vector &y) const
{
	const realtype *y_data = NV_DATA_S(y);
	double a(0.);

	for (int i = nb_of_species; i < nb_of_species + nb_lev_h2; i++) {
		if (rounding(h2_di->lev_array[i-nb_of_species].j)%2 == 0)
			a += y_data[i];
	}
	return a;
}

double evolution_data::calc_ice_conc(const N_Vector &y) const
{
	double a(0.);
	const realtype *y_data = NV_DATA_S(y);
	
	for (int i = nb_of_gas_species; i < nb_of_species; i++) {
		a += y_data[i];
	}
	return a;
}

double evolution_data::calc_conc_h_tot(const N_Vector &y) const 
{
	double a = 0.;
	const realtype *y_data = NV_DATA_S(y);

	for (int i = 0; i < nb_of_species; i++) {
		if (network->species[i].formula[0] > 0) {
			a += network->species[i].formula[0] *y_data[i];
		}
	}
	return a;
}

double evolution_data::calc_hydrocarbon_conc(const N_Vector &y) const
{
	bool failed;
	int i, j;
	double a(0.);

	for (i = nb_of_gas_species; i < nb_of_species; i++) {
		failed = false;
		for (j = 0; j < NB_OF_CHEM_ELEMENTS; j++) {
			if (j != 0 && j != 2 && network->species[i].formula[j] != 0) {
				failed = true;
				break;
			}
		}
		if (!failed && network->species[i].formula[2] >= MIN_CARBON_ATOMS_IN_HYDROCARBONS)
			a += NV_Ith_S(y, i);
	}
	return a;
}

void evolution_data::calc_ion_dens(const N_Vector &y, double & ion_conc, double & ion_pah_conc, double & ion_dens, 
    double & ion_pah_dens) const
{
	int i, k;
	const realtype *y_data = NV_DATA_S(y);

	// calculation of ion mass density and concentration due to chemical species:
	ion_conc = ion_dens = 0.;
	for (i = 0; i < nb_of_gas_species; i++) { 
		if (network->species[i].charge != 0 && i != network->e_nb) 
		{
			ion_conc += y_data[i];
			ion_dens += y_data[i] *network->species[i].mass;
		}
	}

	ion_pah_conc = ion_pah_dens = 0.;
	// calculation of charged dust contribution to ion mass:
	for (i = 0; i < dust->nb_of_comp; i++) 
	{ 
		if (is_charged_dust_as_ions[i]) {
			if (is_dust_charge_large[i]) 
			{
				ion_pah_conc += y_data[nb_dch[i]];
				ion_pah_dens += y_data[nb_dch[i]] *dust->get_grain_mass(i);
			}
			else {
				for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) {
					if (min_grain_charge[i] + k != 0) 
					{
						ion_pah_conc += y_data[nb_dch[i] + k];
						ion_pah_dens += y_data[nb_dch[i] + k] *dust->get_grain_mass(i);
					}
				}
			}
		}
	}
}

void evolution_data::calc_neutral_dens(const N_Vector &y, double & neut_conc, double & neut_mass_dens) const
{
	int i;
	const realtype *y_data = NV_DATA_S(y);

	// calculation of neutral mass and concentration due to chemical species:
	neut_conc = neut_mass_dens = 0.;
	for (i = 0; i < nb_of_gas_species; i++) {
		if (network->species[i].charge == 0) 
		{
			neut_conc += y_data[i];
			neut_mass_dens += y_data[i] *network->species[i].mass;
		}
	}

	// calculation of dust contribution to neutral mass:
	for (i = 0; i < dust->nb_of_comp; i++) { 
		if (is_charged_dust_as_ions[i] && !is_dust_charge_large[i]) {
			if (min_grain_charge[i] <= 0 && max_grain_charge[i] >= 0) 
			{
				neut_conc += y_data[nb_dch[i] - min_grain_charge[i]];
				neut_mass_dens += y_data[nb_dch[i] - min_grain_charge[i]] *dust->get_grain_mass(i);
			}
		}
	}
}

double evolution_data::calc_total_gas_charge(const N_Vector &y) const
{
	int i;
	double a;
	const realtype *y_data = NV_DATA_S(y);

	a = -y_data[network->e_nb];
	for (i = 0; i < nb_of_species; i++) 
	{
		if (network->species[i].type == "ion") {
			a += y_data[i] *network->species[i].charge;
		}
	}
	return a;
}

double evolution_data::calc_total_grain_charge(const N_Vector &y) const
{
	int i, j;
	double a(0.);
	const realtype *y_data = NV_DATA_S(y);

	for (i = 0; i < nb_of_dust_comp; i++) 
	{
		if (is_dust_charge_large[i]) {
			a += y_data[nb_dch[i]] *y_data[nb_dch[i] + 1];
		}
		else {
			for (j = 0; j < nb_dch[i+1] - nb_dch[i]; j++) {
				a += (min_grain_charge[i] + j) *y_data[nb_dch[i] + j];
			}
		}
	}
	return a;
}

double evolution_data::calc_average_grain_charge(const N_Vector &y, int i) const
{
	if (is_dust_charge_large[i]) {
		return NV_Ith_S(y, nb_dch[i] + 1);
	}
	else 
	{
		double b(0.), conc(0.);
		for (int k = 0; k < nb_dch[i+1] - nb_dch[i]; k++) 
		{
			conc += NV_Ith_S(y, nb_dch[i] + k);
			b += (min_grain_charge[i] + k) *NV_Ith_S(y, nb_dch[i] + k); 
		}
		return b/conc;
	}
}

double evolution_data::calc_grain_conc(const N_Vector &y, int i) const
{
	if (is_dust_charge_large[i]) {
		return NV_Ith_S(y, nb_dch[i]);
	}
	else 
	{
		double conc(0.);
		for (int k = 0; k < nb_dch[i+1] - nb_dch[i]; k++) {
			conc += NV_Ith_S(y, nb_dch[i] + k);
		}
		return conc;
	}
}

double evolution_data::calc_dust_ion_coupling(const N_Vector &y, int i) const
{
	int k;
	double a, b, c, wt2, r, s;

	if (is_dust_charge_large[i])
	{
		calc_grain_velocity(y, NV_Ith_S(y, nb_dch[i] + 1), dust->get_grain_area(i), a, b, wt2, c);
		r = wt2/(1. + wt2);
	}
	else
	{
		s = r = 0.;
		for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++)
		{
			calc_grain_velocity(y, min_grain_charge[i] + k, dust->get_grain_area(i), a, b, wt2, c);
			r += wt2/(1. + wt2) *NV_Ith_S(y, nb_dch[i] + k);
			s += NV_Ith_S(y, nb_dch[i] + k);
		}
		r /= s;
	}
	return r;
}

void evolution_data::set_grain_charge_ranges(bool *is_av_ch, int *min_gch, int *max_gch, int nb_of_chs)
{
	memcpy(is_dust_charge_large, is_av_ch, nb_of_dust_comp*sizeof(bool));
	memcpy(min_grain_charge, min_gch, nb_of_dust_comp*sizeof(int));
	memcpy(max_grain_charge, max_gch, nb_of_dust_comp*sizeof(int));

	nb_of_equat += nb_of_chs - nb_of_grain_charges;
	nb_of_grain_charges = nb_of_chs;

	for (int i = 0; i < nb_of_dust_comp; i++)
	{
		if (is_av_ch[i]) {
			nb_dch[i+1] = nb_dch[i] + 2;
		}
		else {
			nb_dch[i+1] = nb_dch[i] + max_grain_charge[i] - min_grain_charge[i] + 1;
		}
	}
	nb_dct = nb_dch[nb_of_dust_comp];
	nb_mhd = nb_dct + nb_of_dust_comp;

	reinit_arrays(nb_of_grain_charges);
}

void evolution_data::recalc_grain_charge_ranges(double temp_e, double conc_e, double conc_h_tot, vector<double> & init_ch)
{
	int i, j, k, l, n, zlim, zmin, zmax, *new_nb_dch;
	double a, b, c, z, z1, z2;
	
	init_ch.clear();
	
	new_nb_dch = new int [nb_of_dust_comp+1];
	new_nb_dch[0] = nb_dch[0];
	
	c = 0.5*conc_e *sqrt(8.*temp_e*BOLTZMANN_CONSTANT/(M_PI*ELECTRON_MASS));
	for (i = 0; i < nb_of_dust_comp; i++)
	{
		z1 = zmin = dust->components[i]->zmin;
		z2 = zmax = dust->components[i]->zmax;
		
		while (z2 - z1 > 0.01)
		{
			z = 0.5*(z1 + z2);
			a = -c*(1. - exp(-1.e+7*dust->get_grain_radius(i))) *dust->get_grain_area(i) 
				*accr_func->get_accr_rate(dust->get_grain_radius(i)*temp_e, -z);

			a += photoem_factor_is_uv *dust->components[i]->get_isrf_uv_phel_rate(z, b);
			a += photoem_factor_is_vis *dust->components[i]->get_isrf_vis_phel_rate(z, b);
			// molecularization degree is not taken into account:
			a += cr_ioniz_rate *dust->components[i]->get_cr_phel_rate(z, b); 
				
			if (a > 0.)
				z1 = z;
			else z2 = z;
		}
		
		if (fabs(z) > NB_CHARGE_BINS_LARGE_GRAINS/2)
		{
			new_nb_dch[i+1] = new_nb_dch[i] + 2;
			is_dust_charge_large[i] = true;

			init_ch.push_back(dust->components[i]->conc_perH *conc_h_tot);
			init_ch.push_back(z);
		}
		else
		{
			if (z > 0.) 
				n = (int) (z + 0.5);
			else n = (int) (z - 0.5);

			zlim = abs(n) + 3;
			if (zlim > NB_CHARGE_BINS_LARGE_GRAINS/2)
				zlim = NB_CHARGE_BINS_LARGE_GRAINS/2;

			k = n - zlim;
			if (k < zmin) 
				k = zmin;

			l = n + zlim + 1;
			if (l > zmax)
				l = zmax;

			min_grain_charge[i] = k;
			max_grain_charge[i] = l;

			is_dust_charge_large[i] = false;
			new_nb_dch[i+1] = new_nb_dch[i] + max_grain_charge[i] - min_grain_charge[i] + 1;

			n = (int) z;
			if (z < 0.)
				n--;

			for (j = min_grain_charge[i]; j < n; j++) {
				init_ch.push_back(0.);
			}
			init_ch.push_back(dust->components[i]->conc_perH *conc_h_tot *(n + 1. - z));
			init_ch.push_back(dust->components[i]->conc_perH *conc_h_tot *(z - n));
			
			for (j = n + 2; j <= max_grain_charge[i]; j++) {
				init_ch.push_back(0.);
			}
		}
	}
	nb_of_equat += new_nb_dch[nb_of_dust_comp] - nb_dch[nb_of_dust_comp];
	memcpy(nb_dch, new_nb_dch, (nb_of_dust_comp + 1)*sizeof(int));

	nb_of_grain_charges = nb_dch[nb_of_dust_comp] - nb_dch[0];
	nb_dct = nb_dch[nb_of_dust_comp];
	nb_mhd = nb_dct + nb_of_dust_comp;

	reinit_arrays(nb_of_grain_charges);
	delete [] new_nb_dch;
}

bool evolution_data::recalc_grain_charge_ranges(N_Vector y, vector<double> & new_y)
{
	bool must_be_restarted = false;
	int i, j, k, l, n, zlim, zmin, zmax, *new_nb_dch;
	double a, b, c, z, dz;

    cout << scientific;
    cout.precision(3);

	dz = 0.;
	new_y.clear();

	for (i = 0; i < nb_dch[0]; i++) {
		new_y.push_back(NV_Ith_S(y, i));
	}

	new_nb_dch = new int [nb_of_dust_comp+1];
	new_nb_dch[0] = nb_dch[0];

	for (i = 0; i < nb_of_dust_comp; i++)
	{
		zmin = dust->components[i]->zmin;
		zmax = dust->components[i]->zmax;	
		z = calc_average_grain_charge(y, i);
		
		if (!is_dust_charge_large[i]) {
			if (fabs(z) < NB_CHARGE_BINS_LARGE_GRAINS/2 + 1)
			{
				if (z > 0.) 
					n = (int) (z + 0.5);
				else n = (int) (z - 0.5);
			
				if ( ((max_grain_charge[i] + min_grain_charge[i])/2 <= n - 2 && max_grain_charge[i] < zmax)
					|| ((max_grain_charge[i] + min_grain_charge[i])/2 >= n + 2 && min_grain_charge[i] > zmin) )
				{
					zlim = abs(n) + 3;
					if (zlim > NB_CHARGE_BINS_LARGE_GRAINS/2)
						zlim = NB_CHARGE_BINS_LARGE_GRAINS/2;

					k = n - zlim;
					if (k < zmin) 
						k = zmin;

					l = n + zlim + 1;
					if (l > zmax)
						l = zmax;

					new_nb_dch[i+1] = new_nb_dch[i] + l - k + 1;

					for (j = 0; j < min_grain_charge[i] - k; j++) {
						new_y.push_back(0.);
					}

					a = b = 0.;
					for (j = 0; j < max_grain_charge[i] - min_grain_charge[i] + 1; j++)
					{
						if (min_grain_charge[i] + j >= k && min_grain_charge[i] + j <= l) {
							new_y.push_back(NV_Ith_S(y, nb_dch[i] + j));
						}
						else {
							a += NV_Ith_S(y, nb_dch[i] + j);
							b += (min_grain_charge[i] + j) *NV_Ith_S(y, nb_dch[i] + j);
						}
					}
					for (j = 0; j < l - max_grain_charge[i]; j++) {
						new_y.push_back(0.);
					}
					
					// in order to conserve dust grain mass (but the charge is not conserved):
					a /= l - k + 1;
					for (j = 0; j < l - k + 1; j++) 
					{
						new_y[new_nb_dch[i] + j] += a;
						b -= (k + j)*a;
					}
					dz += b;

					if (verbosity)
					{
						cout << left << "dust component: " << i << endl
							<< "    old distribution; z_min, z_max, conc (cm-3):" << endl
							<< setw(8) << min_grain_charge[i] << setw(8) << max_grain_charge[i];

						for (j = nb_dch[i]; j < nb_dch[i+1]; j++) {
							cout << left << setw(12) << NV_Ith_S(y, j);
						}
						cout << endl;

						cout << "    new distribution; z_min, z_max, conc (cm-3):" << endl
							<< setw(8) << k << setw(8) << l;

						for (j = 0; j < l - k + 1; j++) {
							cout << left << setw(12) << new_y[new_nb_dch[i] + j];
						}
						cout << endl;
						cout << "    grain charge lost dz = " << b << endl;
					}
					
					min_grain_charge[i] = k;
					max_grain_charge[i] = l;		
					must_be_restarted = true;
				}
				else 
				{
					for (j = nb_dch[i]; j < nb_dch[i+1]; j++){
						new_y.push_back(NV_Ith_S(y, j));
					}
					new_nb_dch[i+1] = new_nb_dch[i] + max_grain_charge[i] - min_grain_charge[i] + 1;
				}
			}
			else if ((z > 0. && zmax > NB_CHARGE_BINS_LARGE_GRAINS) || (z < 0. && zmin < -NB_CHARGE_BINS_LARGE_GRAINS))
			{
				a = calc_grain_conc(y, i);			
				new_y.push_back(a);
				new_y.push_back(z);

				new_nb_dch[i+1] = new_nb_dch[i] + 2;

				is_dust_charge_large[i] = true;
				must_be_restarted = true;

				if (verbosity)
				{
					cout << "dust componet: " << i << endl
						<< "    old distribution; min z, max z, conc (cm-3):" << endl;

					cout << left << setw(8) << min_grain_charge[i] << setw(8) << max_grain_charge[i];
					for (j = nb_dch[i]; j < nb_dch[i+1]; j++) {
						cout << left << setw(12) << NV_Ith_S(y, j);
					}
					cout << endl;
					cout << left << "    average parameters; z, conc (cm-3), vel (cm/s): " 
						<< setw(12) << z << setw(12) << a << setw(12) << calc_grain_velocity(y, z, dust->get_grain_area(i), c, c, c, c) << endl;
				}
			}
			else 
			{
				for (j = nb_dch[i]; j < nb_dch[i+1]; j++){
					new_y.push_back(NV_Ith_S(y, j));
				}
				new_nb_dch[i+1] = new_nb_dch[i] + max_grain_charge[i] - min_grain_charge[i] + 1;
			}
		}
		else
		{
			// charge value of 0 must be included in distribution;
			if (fabs(z) < NB_CHARGE_BINS_LARGE_GRAINS/2)
			{
				if (z > 0.) 
					n = (int) (z + 0.5);
				else n = (int) (z - 0.5);
			
				zlim = abs(n) + 3;	
				if (zlim > NB_CHARGE_BINS_LARGE_GRAINS/2)
					zlim = NB_CHARGE_BINS_LARGE_GRAINS/2;

				k = n - zlim;
				if (k < zmin) 
					k = zmin;

				l = n + zlim + 1;
				if (l > zmax)
					l = zmax;
			
				n = (int) z;
				if (z < 0.)
					n--;

				a = NV_Ith_S(y, nb_dch[i]);
				for (j = k; j <= l; j++)
				{
					if (j == n)
						new_y.push_back(a*(n + 1 - z));
					else if (j == n + 1)
						new_y.push_back(a*(z - n));
					else new_y.push_back(0.);
				}
			
				min_grain_charge[i] = k;
				max_grain_charge[i] = l;	
				new_nb_dch[i+1] = new_nb_dch[i] + l - k + 1;
				
				is_dust_charge_large[i] = false;
				must_be_restarted = true;		
		
				if (verbosity) {
					cout << left << "dust componet: " << i << endl
						 << "    average parameters; z, conc (cm-3): " << setw(12) << z << setw(12) << a << endl;
					
					cout << left << "    new distribution; z_min, z_max, conc (cm-3): " << setw(8) << k << setw(8) << l;
					for (j = 0; j < l - k + 1; j++) {
						cout << left << setw(12) << new_y[new_nb_dch[i]+j];
					}
					cout << endl;
				}
			}
			else 
			{
				new_y.push_back(NV_Ith_S(y, nb_dch[i]));
				new_y.push_back(NV_Ith_S(y, nb_dch[i] + 1));
				new_nb_dch[i+1] = new_nb_dch[i] + 2;
			}
		}
	}
	for (i = nb_dch[nb_of_dust_comp]; i < nb_of_equat; i++) {
		new_y.push_back(NV_Ith_S(y, i));
	}

	if (must_be_restarted)
	{
		// a trick for charge conservation: 
		new_y[network->e_nb] -= dz;
		if (verbosity) {
            cout << "new ranges for grain charge distribution were adopted," << endl 
                << "concentration of electrons was changed (to compensate grain charge lost) by (cm-3): " << -dz << endl;
		}

		nb_of_equat += new_nb_dch[nb_of_dust_comp] - nb_dch[nb_of_dust_comp];
		memcpy(nb_dch, new_nb_dch, (nb_of_dust_comp + 1)*sizeof(int));

		nb_of_grain_charges = nb_dch[nb_of_dust_comp] - nb_dch[0];
		nb_dct = nb_dch[nb_of_dust_comp];
		nb_mhd = nb_dct + nb_of_dust_comp;

		reinit_arrays(nb_of_grain_charges);
	}	
	delete [] new_nb_dch;
	return must_be_restarted;
}

void evolution_data::get_nb_of_levels(int & nb1, int & nb2, int & nb3, int & nb4, int & nb5, int & nb6, int & nb7, int & nb8, 
	int & nb9, int & nb10) const
{
	nb1 = nb_lev_h2;
	nb2 = nb_lev_h2o;
	nb3 = nb_lev_co;
	nb4 = nb_lev_oh;
	nb5 = nb_lev_pnh3;
	nb6 = nb_lev_onh3;
	nb7 = nb_lev_ch3oh;
	nb8 = nb_lev_ci; 
	nb9 = nb_lev_oi;
	nb10 = nb_lev_cii;
}

void evolution_data::get_nbs(int & nb1, int & nb2, int & nb3, int & nb4) const
{
	nb1 = nb_of_grain_charges;
	nb2 = nb_of_equat;
	nb3 = nb_dct;
	nb4 = nb_mhd;
}

void evolution_data::get_dust_component_nbs(int i, int & nb1, int & nb2) const
{
	nb1 = nb_dch[i];
	nb2 = nb_dch[i+1];
}

void evolution_data::get_neutral_heating(double & atomic_n, double & h2_n, double & h2o_n, double & co_n, double & oh_n, 
	double &nh3_n, double & ch3oh_n, double & coll_h, double & chem_h, double &pheff_h, double & cr, double & scatt_i, 
	double & scatt_e, double & rel_h2) const
{
	atomic_n = neut_heat_atoms;
	h2_n = neut_heat_h2;
	h2o_n = neut_heat_ph2o + neut_heat_oh2o;
	co_n = neut_heat_co;
	oh_n = neut_heat_oh;
	nh3_n = neut_heat_pnh3 + neut_heat_onh3;
	ch3oh_n = neut_heat_ch3oh_a + neut_heat_ch3oh_e;

	chem_h = neut_heat_chem;
	coll_h = neut_heat_dust_coll;
	pheff_h = pheff_gas_heat;
	cr = neut_cr_heat;

	scatt_i = neut_heat_scatt_ions; 
	scatt_e = neut_heat_scatt_el;

	rel_h2 = rad_energy_loss_h2;
}

void evolution_data::get_electron_heating(double & atomic_e, double & h2_e, double & h2o_e, double & scatt_n, double & scatt_i, double & chem) const
{
	atomic_e = el_heat_atoms;
	h2_e = el_heat_h2;
	h2o_e = el_heat_ph2o + el_heat_oh2o;

	scatt_n = el_heat_scatt_neut; 
	scatt_i = el_heat_scatt_ions;
	chem = el_heat_chem;
}

void evolution_data::get_ion_heating(double & h2_i, double & scatt_n, double & scatt_e, double & chem) const
{
    h2_i = ion_heat_h2;
	scatt_n = ion_heat_scatt_n; 
	scatt_e = ion_heat_scatt_el;
	chem = ion_heat_chem;
}

// dust heating normalized on one dust grain, erg s-1:
void evolution_data::get_dust_heating_rates(int i, double & isrf, double & chem, double &coll, double &h2_mol, double & mol) const
{
	if (i >= 0 && i < nb_of_dust_comp)
	{
		isrf = dh_isrf_arr[i];
		chem = dust_heat_chem[i];
		coll = dust_heat_coll[i];
		h2_mol = dust_heat_h2_line[i] *CM_INVERSE_TO_ERG; // conversion of cm-1 to erg;
		mol = dust_heat_mline[i] *CM_INVERSE_TO_ERG;
	}
}

// grain charging rates have dimension cm-3 s-1
void evolution_data::get_dust_charging_rates(int i, double & el_att, double & ion_neutr, double & phel_uv, 
	double & phel_vis, double & phel_cr) const
{
	if (i >= 0 && i < nb_of_dust_comp)
	{
		el_att = el_att_rate[i];
		ion_neutr = ion_neutr_rate[i];
		phel_uv = phel_rate_uv[i];
		phel_vis = phel_rate_vis[i];
		phel_cr = phel_rate_cr[i];
	}
}

double evolution_data::get_reaction_rate(int i) const
{
	if (i >= 0 && i < network->nb_of_reactions)
		return chem_reaction_rates[i];
	else return 0.;
}

void evolution_data::get_chem_heating_rates(double *& cheat, int & nb) const
{
	cheat = chem_heating_rates_n;
	nb = network->nb_of_reactions;
}

void evolution_data::get_h2_chem(double & h2_gr, double & h2_gas, double & o_grains, double & o_gchem, double & o_hcoll,
	double & h2_h_diss) const 
{ 
	h2_gr = h2_prod_gr;
	h2_gas = h2_prod_gas;
	o_gchem = oh2_form_gaschem; 
	o_grains = oh2_form_grains;
	o_hcoll = oh2_form_hcoll;
	h2_h_diss = h2_h_diss_rate;
}

void evolution_data::set_tolerances(N_Vector abs_tol)
{
    int i, j;
    for (i = 0; i < nb_of_species; i++) {
        NV_Ith_S(abs_tol, i) = ABS_CONCENTRATION_ERROR_SOLVER;
    }

    for (i = nb_of_species; i < nb_dch[0]; i++) {
        NV_Ith_S(abs_tol, i) = ABS_POPULATION_ERROR_SOLVER;
    }
    for (i = 0; i < nb_of_dust_comp; i++) {
        if (is_dust_charge_large[i]) {
            NV_Ith_S(abs_tol, nb_dch[i]) = ABS_CONCENTRATION_ERROR_SOLVER; // grain concentration
            NV_Ith_S(abs_tol, nb_dch[i] + 1) = ABS_PARAMETER_ERROR_SOLVER; // average grain charge
        }
        else {
            for (j = nb_dch[i]; j < nb_dch[i + 1]; j++) {
                NV_Ith_S(abs_tol, j) = ABS_CONCENTRATION_ERROR_SOLVER;
            }
        }
    }
    for (i = nb_dct; i < nb_of_equat; i++) {
        NV_Ith_S(abs_tol, i) = ABS_PARAMETER_ERROR_SOLVER;
    }
}

void evolution_data::set_parameters(double vis_ext, double cr_ion, double uv_field, double ir_field)
{
	visual_extinct = vis_ext;
	cr_ioniz_rate = cr_ion;
	uv_field_strength = uv_field;
	ir_field_strength = ir_field;
	
	// approximate scaling with visual extinction, for dust charging:
	photoem_factor_is_uv = uv_field_strength *exp(-2.*visual_extinct);
	// the rates are calculated for photons < 5 eV, may be more accurate scaling?:
	photoem_factor_is_vis = exp(-visual_extinct);

	photodes_factor_is = uv_field_strength *DRAINE1978_ISRF_FUV *exp(-2.*visual_extinct)/(8.*GRAIN_SITES_PER_CM2);
	photodes_factor_cr = STANDARD_NB_CR_PHOTONS*cr_ioniz_rate/(8.*GRAIN_SITES_PER_CM2); // the dependence on R_V is not taken into account,
	desorption_factor_cr = STANDARD_FLUX_CR_IRON *cr_ioniz_rate/(8.*GRAIN_SITES_PER_CM2);

	dheat_isrf->get_heating(visual_extinct, ir_field_strength, uv_field_strength, cr_ioniz_rate, dh_isrf_arr, nb_of_dust_comp);
}

void evolution_data::set_veln_grad(double vg) 
{
	if (fabs(vg) > vel_grad_min)
		vel_n_grad = vg; // velocity gradient is allowed to be negative;
	else if (vg > 0.) 
		vel_n_grad = vel_grad_min; 
	else vel_n_grad = -vel_grad_min; 
}

void evolution_data::set_veli_grad(double vg) 
{
	if (fabs(vg) > vel_grad_min)
		vel_i_grad = vg;
	else if (vg > 0.) 
		vel_i_grad = vel_grad_min;
	else vel_i_grad = -vel_grad_min;
}

void evolution_data::create_file_radiative_transfer(const string & output_path) const
{
	string fname = output_path + "sim_rad_transf.txt";
	ofstream output;

	output.open(fname.c_str());
	output.close();
}

void evolution_data::save_radiative_transfer_data(const string & output_path, double ty) const
{
	const double min_heff1 = 1.e-4, min_heff2 = 3.e-3; // for H2 and other molecules
	int l, i, k, nb;
	double hr1, hr2;

	string fname = output_path + "sim_rad_transf.txt";
	ofstream output;
	
	hr1 = hr2 = 0.;
	for (i = 0; i < nb_of_dust_comp; i++) {
		hr1 += dust_heat_h2_line[i] *grain_conc[i]; // heating rate of dust by H2 molecules in cm-3 of gas, cm-1 cm-3 s-1
	}
	for (i = 0; i < nb_of_dust_comp; i++) {
		hr2 += dust_heat_mline[i] *grain_conc[i]; // heating by all molecules except H2, cm-1 cm-3 s-1
	}

	output.open(fname.c_str(), ios::app);
	output << scientific;
	output.precision(2);
	
	output << "# Evolution age (years):" << endl;
	output << ty << endl;
	output << "# Dust heating rate by H2 [cm-1 cm-3 s-1]: " << hr1 << endl
		<< "# Dust heating rate by other molecules [cm-1 cm-3 s-1]: " << hr2 << endl;

	// effici - the relative contribution of a line to the heating process;
	output << left << setw(5) << "# nb" << setw(5) << "v_u" << setw(5) << "j_u" << setw(5) << "t_u" 
		<< setw(5) << "v_l" << setw(5) << "j_l" << setw(5) << "t_l" << setw(10) << "en(cm-1)" 
		<< setw(10) << "g" << setw(10) << "d" << setw(10) << "effici" << setw(10) << "ep_int1" << setw(10) << "ep_int2" << endl;

	// H2
	output << "# Radiative transfer data for H2" << endl;
	k = 0;
	for (l = 0; l < h2_di->nb_lev-1; l++) {
		for (i = l+1; i < h2_di->nb_lev; i++) 
		{
			if (h2_einst->arr[i][l] > DBL_EPSILON && dheat_efficiency[i*(i-1)/2+l]/hr1 > min_heff1)
			{
				output << left << setw(5) << k << setw(5) << h2_di->lev_array[i].v << setw(10) << rounding(h2_di->lev_array[i].j) 
					<< setw(5) << h2_di->lev_array[l].v << setw(10) << rounding(h2_di->lev_array[l].j) << setw(10) << h2_di->lev_array[i].energy - h2_di->lev_array[l].energy 
					<< setw(10) << gamma_factors[i*(i-1)/2+l] << setw(10) << delta_factors[i*(i-1)/2+l] << setw(10) << dheat_efficiency[i*(i-1)/2+l]/hr1 
					<< setw(10) << esc_prob_int1[i*(i-1)/2+l] << setw(10) << esc_prob_int2[i*(i-1)/2+l] << endl;
				k++;
			}
		}
	}
	// p-H2O
	output << "# Radiative transfer for p-H2O" << endl;
	k = 0;
	nb = nb_lev_h2*(nb_lev_h2-1)/2;
	
	for (l = 0; l < ph2o_di->nb_lev-1; l++) {
		for (i = l+1; i < ph2o_di->nb_lev; i++) 
		{
			if (ph2o_einst->arr[i][l] > DBL_EPSILON && dheat_efficiency[nb+i*(i-1)/2+l]/hr2 > min_heff2)
			{
				output << left << setw(5) << k << setw(5) << ph2o_di->lev_array[i].v << setw(5) << rounding(ph2o_di->lev_array[i].j) 
					<< setw(5) << rounding(ph2o_di->lev_array[i].k1 - ph2o_di->lev_array[i].k2)
					<< setw(5) << ph2o_di->lev_array[l].v << setw(5) << rounding(ph2o_di->lev_array[l].j) << setw(5) << rounding(ph2o_di->lev_array[l].k1 - ph2o_di->lev_array[l].k2)
					<< setw(10) << ph2o_di->lev_array[i].energy - ph2o_di->lev_array[l].energy 
					<< setw(10) << gamma_factors[nb+i*(i-1)/2+l] << setw(10) << delta_factors[nb+i*(i-1)/2+l] << setw(10) << dheat_efficiency[nb+i*(i-1)/2+l]/hr2 
					<< setw(10) << esc_prob_int1[nb+i*(i-1)/2+l] << setw(10) << esc_prob_int2[nb+i*(i-1)/2+l] << endl;
				k++;
			}
		}
	}

	// o-H2O
	output << "# Radiative transfer for o-H2O" << endl;
	k = 0;
	nb += nb_lev_h2o*(nb_lev_h2o-1)/2;
	
	for (l = 0; l < oh2o_di->nb_lev-1; l++) {
		for (i = l+1; i < oh2o_di->nb_lev; i++) 
		{
			if (oh2o_einst->arr[i][l] > DBL_EPSILON && dheat_efficiency[nb+i*(i-1)/2+l]/hr2 > min_heff2)
			{
				output << left << setw(5) << k << setw(5) << oh2o_di->lev_array[i].v << setw(5) << rounding(oh2o_di->lev_array[i].j) 
					<< setw(5) << rounding(oh2o_di->lev_array[i].k1 - oh2o_di->lev_array[i].k2)
					<< setw(5) << oh2o_di->lev_array[l].v << setw(5) << rounding(oh2o_di->lev_array[l].j) << setw(5) << rounding(oh2o_di->lev_array[l].k1 - oh2o_di->lev_array[l].k2)
					<< setw(10) << oh2o_di->lev_array[i].energy - oh2o_di->lev_array[l].energy 
					<< setw(10) << gamma_factors[nb+i*(i-1)/2+l] << setw(10) << delta_factors[nb+i*(i-1)/2+l] << setw(10) << dheat_efficiency[nb+i*(i-1)/2+l]/hr2 
					<< setw(10) << esc_prob_int1[nb+i*(i-1)/2+l] << setw(10) << esc_prob_int2[nb+i*(i-1)/2+l] << endl;
				k++;
			}
		}
	}
	// CO
	output << "# Radiative transfer for CO" << endl;
	k = 0;
	nb += nb_lev_h2o*(nb_lev_h2o-1)/2;
	
	for (l = 0; l < co_di->nb_lev-1; l++) {
		for (i = l+1; i < co_di->nb_lev; i++) 
		{
			if (co_einst->arr[i][l] > DBL_EPSILON && dheat_efficiency[nb+i*(i-1)/2+l]/hr2 > min_heff2)
			{
				output << left << setw(5) << k << setw(5) << co_di->lev_array[i].v << setw(10) << rounding(co_di->lev_array[i].j) 
					<< setw(5) << co_di->lev_array[l].v << setw(10) << rounding(co_di->lev_array[l].j) << setw(10) << co_di->lev_array[i].energy - co_di->lev_array[l].energy 
					<< setw(10) << gamma_factors[nb+i*(i-1)/2+l] << setw(10) << delta_factors[nb+i*(i-1)/2+l] << setw(10) << dheat_efficiency[nb+i*(i-1)/2+l]/hr2 
					<< setw(10) << esc_prob_int1[nb+i*(i-1)/2+l] << setw(10) << esc_prob_int2[nb+i*(i-1)/2+l] << endl;
				k++;
			}
		}
	}
	output.close();
}

chemistry_evolution_data::chemistry_evolution_data(const std::string &input_path, const std::string &output_path, int nb_h2, 
	int nb_vibr_h2o, int nb_h2o, int nb_vibr_co, int nb_co, int nb_pnh3, int nb_onh3, int nb_oh, int nb_vibr_ch3oh, int nb_lev_ch3oh, double c_abund_pah, int verb)
	: evolution_data(input_path, output_path, nb_h2, nb_vibr_h2o, nb_h2o, nb_vibr_co, nb_co, nb_pnh3, nb_onh3, nb_oh, 
        nb_vibr_ch3oh, nb_lev_ch3oh, c_abund_pah, verb)
{
	// please, check the nb of equations. For chemical evolution calculations:
	nb_of_equat = nb_of_species + nb_lev_h2 + 2*nb_lev_h2o + nb_lev_co + nb_lev_oh + nb_lev_pnh3 + nb_lev_onh3  
		+ 2*nb_lev_ch3oh + nb_lev_ci + nb_lev_oi + nb_lev_cii + nb_of_grain_charges + nb_of_dust_comp + 3;

	network->check_reactions(); // must be called in all cases;
	
	// Saving data, only in dark cloud simulations:
    network->print_network(output_path);
    dust->save_data(output_path);
    dheat_isrf->save_data(output_path);

	// elastic_h2_ions->save_data(output_path + "el_scatt_h2_ions.txt"); // h2-ion scattering rates;
	// accr_func->save(output_path); // charge particle accretion rates on dust grains;

	chem_reaction_rates = new double [network->nb_of_reactions];
	memset(chem_reaction_rates, 0, network->nb_of_reactions*sizeof(double));

	chem_heating_rates_n = new double [network->nb_of_reactions];
	memset(chem_heating_rates_n, 0, network->nb_of_reactions*sizeof(double));
	
	if (verbosity)
		cout << "Total nb of differential equations " << nb_of_equat << endl;
}

chemistry_evolution_data::~chemistry_evolution_data() 
{;}

int chemistry_evolution_data::f(realtype t, N_Vector y, N_Vector ydot)
{
	int i, returned_val;
	// must be defined before parent function call:
	vel_ni_diff = vel_n = vel_i = 0.;
    returned_val = evolution_data::f(t, y, ydot);

	// The equations for neutral gas, ions and electron temperatures:
	if (IS_TEMPERATURE_FIXED) 
	{
		// gas temperature derivative:
		NV_Ith_S(ydot, nb_mhd) = NV_Ith_S(ydot, nb_mhd + 1) = NV_Ith_S(ydot, nb_mhd + 2) = 0.;
		
		// dust temperarure derivative:
		for (i = 0; i < nb_of_dust_comp; i++) {
			NV_Ith_S(ydot, nb_dct + i) = 0.;
		}
	}
	else
	{
		NV_Ith_S(ydot, nb_mhd) = (2.*energy_gain_n/(3.*BOLTZMANN_CONSTANT) - temp_n* nb_gain_n)/conc_n;
		NV_Ith_S(ydot, nb_mhd + 1) = (2.*energy_gain_i/(3.*BOLTZMANN_CONSTANT) - temp_i* nb_gain_i)/conc_i;
		NV_Ith_S(ydot, nb_mhd + 2) = (2.*energy_gain_e/(3.*BOLTZMANN_CONSTANT) - temp_e* nb_gain_e)/conc_e;
	}
	return returned_val;
}

//
// Methods of the class calculate shock evolution
//

mhd_shock_data::mhd_shock_data(const string &input_path, const std::string &output_path, int nb_h2, int nb_vibr_h2o, int nb_h2o, 
	int nb_vibr_co, int nb_co, int nb_pnh3, int nb_onh3, int nb_oh, int nb_vibr_ch3oh, int nb_ch3oh, double c_abund_pah, int verb)
	: evolution_data(input_path, output_path, nb_h2, nb_vibr_h2o, nb_h2o, nb_vibr_co, nb_co, nb_pnh3, nb_onh3, nb_oh, 
        nb_vibr_ch3oh, nb_ch3oh, c_abund_pah, verb),
	magn_field_energy(0.), add_el_source(0.), velg_mhd_n(0.), velg_mhd_i(0.), ion_vg_denominator(0.), neut_vg_denominator(0.)
{
	// calculation of the number of equations:
	nb_of_equat = nb_of_species + nb_lev_h2 + 2*nb_lev_h2o + nb_lev_co + nb_lev_oh + nb_lev_pnh3 + nb_lev_onh3 
		+ 2*nb_lev_ch3oh + nb_lev_ci + nb_lev_oi + nb_lev_cii + nb_of_grain_charges + nb_of_dust_comp + NB_MHD_EQUATIONS;

	// chemical reactions relevant for shock modelling:
	if (GRAIN_MANTLE_SPUTTERING_ON) {
		network->init_network_umistf(input_path + "chemistry/reactions_grain_mantle_sputt.txt");
	}
	
	network->check_reactions(); // must be called in all cases;
//  the size of the file is about 1 Mbyte, in order to save disk place - comment:
//	network->print_network(output_path); 

	chem_reaction_rates = new double [network->nb_of_reactions];
	memset(chem_reaction_rates, 0, network->nb_of_reactions*sizeof(double));

	chem_heating_rates_n = new double [network->nb_of_reactions];
	memset(chem_heating_rates_n, 0, network->nb_of_reactions*sizeof(double));
	
	if (verbosity)
		cout << "Total nb of differential equations " << nb_of_equat << endl;
}

mhd_shock_data::~mhd_shock_data() 
{;}

int mhd_shock_data::f(realtype t, N_Vector y, N_Vector ydot)
{
   	int i, k, l, returned_val;
	double a, b, c, d, e, en_n, nb_e, g_velg, neut_nb_dens, ion_conc, ion_pah_conc, neut_mass_dens, ion_mass_dens, ion_pah_dens,
        ion_vg_terms, neut_vg_terms;

	vel_n = NV_Ith_S(y, nb_mhd + 3);
	vel_i = NV_Ith_S(y, nb_mhd + 4);
	vel_ni_diff = vel_n - vel_i;

	magn_field = magn_field_0 * shock_vel / vel_i;
	magn_field_energy = magn_field * magn_field / (4. * M_PI);

    returned_val = evolution_data::f(t, y, ydot);
	
	realtype *y_data = NV_DATA_S(y);
	realtype *ydot_data = NV_DATA_S(ydot);

	// Gas-dust scattering, vel_ni_diff = vel_n - vel_i;
	en_n = 0.;
	for (i = 0; i < nb_of_dust_comp; i++) {
		if (is_dust_charge_large[i])
		{
			a = (1. + 4.*HE_TO_H_NB_RATIO)*ATOMIC_MASS_UNIT*conc_h_tot *grain_conc[i] *dust->get_grain_area(i) 
				*sqrt(temp_n_erg/v_const_1 + grain_veln2[nb_dch[i]-nb_dch[0]]) *(vel_n - grain_velz[nb_dch[i] - nb_dch[0]]);
		
			mom_gain_n -= a;
			mom_gain_i += a;
			en_n += a *vel_ni_diff;
		}
		else
		{
			for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) 
			{
				l = nb_dch[i] - nb_dch[0] + k;
				a = (1. + 4.*HE_TO_H_NB_RATIO)*ATOMIC_MASS_UNIT*conc_h_tot *dust->get_grain_area(i) 
					*NV_Ith_S(y, nb_dch[i] + k) *sqrt(temp_n_erg/v_const_1 + grain_veln2[l])*(vel_n - grain_velz[l]);
		
				mom_gain_n -= a;
				mom_gain_i += a;
				en_n += a *vel_ni_diff;
			}
		}
	}
	neut_heat_dust_coll += en_n;
	energy_gain_n += en_n;

	// Note:
	// nb   neutral temperature;
	// nb+1 ion temperature;	
	// nb+2 electron temperature;
	// nb+3 neutral velocity;
	// nb+4 ion velocity;

	calc_neutral_dens(y, neut_nb_dens, neut_mass_dens);
	calc_ion_dens(y, ion_conc, ion_pah_conc, ion_mass_dens, ion_pah_dens);
    ion_mass_dens += ion_pah_dens;
    ion_conc += ion_pah_conc;

	// The MHD equations are given by Draine et al., ApJ 264, p.485 (1983); Draine, MNRAS, 220, p.133 (1986);
	// The derivation of derivatives see in Roberge & Draine, ApJ 350, p.700 (1990);
	// 1. The derivative of the neutral gas velocity. 
	// adiabatic index of the neutral and ion gas is taken to be 5/3; 
	// internal energy gain is taken into account implicitly in energy_gain_n;
    neut_vg_terms = -mom_gain_n * vel_n + mass_gain_n * vel_n * vel_n;
    neut_vg_denominator = 1.666666667 * temp_n_erg * neut_nb_dens - neut_mass_dens * vel_n * vel_n;
    velg_mhd_n = ydot_data[nb_mhd + 3] = (0.66666667 *energy_gain_n + neut_vg_terms) / neut_vg_denominator;

	// 2. The derivative of ion velocity, 
	// -mom_gain_n = mom_gain_i + mom_gain_e, 
	// due to adsorption: -mass_gain_n != mass_gain_i + mass_gain_e
    ion_vg_terms = mom_gain_n * vel_i + mass_gain_i * vel_i * vel_i;
    ion_vg_denominator = 1.666666667 * (temp_i_erg * ion_conc + temp_e_erg * conc_e) - ion_mass_dens * vel_i * vel_i
        + magn_field_energy;
	
    velg_mhd_i = ydot_data[nb_mhd + 4] = (0.66666667 * (energy_gain_i + energy_gain_e) + ion_vg_terms)
        / ion_vg_denominator;
	
	// 3. The derivative of the temperature of the neutral gas, the temperature unit is K:
    ydot_data[nb_mhd] = -neut_vg_terms - temp_n_erg * nb_gain_n
        + (temp_n_erg * neut_nb_dens - neut_mass_dens * vel_n * vel_n) * velg_mhd_n;
	
    ydot_data[nb_mhd] /= BOLTZMANN_CONSTANT * neut_nb_dens * vel_n;
		
	// 4. The derivative of the temperature of the ions,
	// there is some discrepancy with Roberge & Draine (1990), may be T_s <-> T_d in their equations?
    ydot_data[nb_mhd + 1] = -ion_vg_terms + 0.66666667 * (energy_gain_i - energy_gain_e) - 2.* temp_i_erg * nb_gain_i
        + (ion_vg_denominator - 1.333333333 * temp_i_erg * ion_conc) * velg_mhd_i;
	
    ydot_data[nb_mhd + 1] /= 2. * BOLTZMANN_CONSTANT * ion_conc * vel_i;
	
	// 5. The derivative of the electron temperature:
    ydot_data[nb_mhd + 2] = -ion_vg_terms + 0.66666667 * (energy_gain_e - energy_gain_i) - 2. * temp_e_erg * nb_gain_e
        + (ion_vg_denominator - 1.333333333 * temp_e_erg * conc_e) * velg_mhd_i;

	ydot_data[nb_mhd + 2] /= 2. * BOLTZMANN_CONSTANT * conc_e * vel_i;
	
	// saving the source term of electrons:
	add_el_source = ydot_data[network->e_nb];

	// the derivative of the concentration of chemical species, only gas species:
	for (i = 0; i < nb_of_gas_species; i++) {
		if (network->species[i].charge == 0) {
			ydot_data[i] = (ydot_data[i] - y_data[i]*velg_mhd_n)/vel_n;
		}
		else {
			ydot_data[i] = (ydot_data[i] - y_data[i]*velg_mhd_i)/vel_i;
		}
	}

	// the derivative of the population densities of atoms and molecules, CII is an ion;
	for (i = nb_of_species; i < nb_dch[0] - nb_lev_cii; i++) {
		ydot_data[i] = (ydot_data[i] - y_data[i] *velg_mhd_n)/vel_n;
	}

	for (i = nb_dch[0] - nb_lev_cii; i < nb_dch[0]; i++) {
		ydot_data[i] = (ydot_data[i] - y_data[i] *velg_mhd_i)/vel_i;
	}

	// the derivative of grain charge distribution/average charge and concentration. Be carefull: vel_ni_diff = vel_n - vel_i;
	// gradients are necessary: neutral (3) and ion (4) velocities, neutral temperature (0);
	nb_e = 0.;
	e = (velg_mhd_n - velg_mhd_i)/vel_ni_diff;

	c = (2.*ydot_data[network->h2_nb] + ydot_data[network->h_nb])/conc_h_tot;
	c = -(2.*c + 2.*velg_mhd_i/vel_i + ydot_data[nb_mhd]/y_data[nb_mhd]);
	d = 2.*e - ydot_data[nb_mhd]/y_data[nb_mhd];
	
	for (i = 0; i < nb_of_dust_comp; i++) {
		if (is_charged_dust_as_ions[i]) {
			if (is_dust_charge_large[i])
			{
				grain_velz_g[nb_dch[i]-nb_dch[0]] = velg_mhd_i;
				
				// d/dz of the average grain charge:
				ydot_data[nb_dch[i] + 1] = ydot_data[nb_dch[i] + 1]/vel_i;

				// d/dz of grain concentration:
				ydot_data[nb_dch[i]] = -y_data[nb_dch[i]] *velg_mhd_i/vel_i;
			}
			else {
				for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) 
				{
					if (k == -min_grain_charge[i]) // charge is 0
						g_velg = velg_mhd_n;
					else g_velg = velg_mhd_i;
				
					ydot_data[nb_dch[i] + k] = (ydot_data[nb_dch[i] + k] - y_data[nb_dch[i] + k] *g_velg)
						/grain_velz[nb_dch[i] - nb_dch[0] + k];

					grain_velz_g[nb_dch[i] - nb_dch[0] + k] = g_velg;
				}
			}
		}
		else {
			if (is_dust_charge_large[i])
			{
				// d/dz of the average grain charge:
				l = nb_dch[i]-nb_dch[0];
				ydot_data[nb_dch[i] + 1] = ydot_data[nb_dch[i] + 1]/grain_velz[l];
			
				if (wt2_arr[l] > DBL_EPSILON)
				{
					a = alpha[l] *(c + 2.*ydot_data[nb_dch[i] + 1]/y_data[nb_dch[i] + 1]);
					b = (beta[l] - 1.)*d;
			
					// note the sign minus before the largest term:
					g_velg = velg_mhd_n + e*(grain_velz[l] - vel_n) 
						- ((a*(wt2_arr[l] + 1.) + alpha[l] *b/beta[l])/(2.*beta[l]*wt2_arr[l] - alpha[l] + 1.) - wt2_arr[l] *b/beta[l])
						*(grain_velz[l] - vel_n)*(grain_velz[l] - vel_n)/(vel_ni_diff*wt2_arr[l]*wt2_arr[l]);
				}
				else g_velg = velg_mhd_n;

				// d/dz of grain concentration:
				ydot_data[nb_dch[i]] = -y_data[nb_dch[i]] *g_velg/grain_velz[l];
			
				nb_e += y_data[nb_dch[i]] *y_data[nb_dch[i] + 1]*(velg_mhd_i - g_velg) 
					+ (y_data[nb_dch[i]] *ydot_data[nb_dch[i] + 1] + ydot_data[nb_dch[i]] *y_data[nb_dch[i] + 1])
					*(vel_i - grain_velz[l]);

				grain_velz_g[nb_dch[i]-nb_dch[0]] = g_velg;
			}
			else {
				for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) 
				{
					l = nb_dch[i] - nb_dch[0] + k;
					if (wt2_arr[l] > DBL_EPSILON)
					{
						a = alpha[l] *c;
						b = (beta[l] - 1.)*d;
			
						g_velg = velg_mhd_n + e*(grain_velz[l] - vel_n) 
							- ((a*(wt2_arr[l] + 1.) + alpha[l] *b/beta[l])/(2.*beta[l]*wt2_arr[l] - alpha[l] + 1.) - wt2_arr[l] *b/beta[l])
							*(grain_velz[l] - vel_n)*(grain_velz[l] - vel_n)/(vel_ni_diff*wt2_arr[l]*wt2_arr[l]);
					}
					else g_velg = velg_mhd_n;

					ydot_data[nb_dch[i] + k] = (ydot_data[nb_dch[i] + k] - y_data[nb_dch[i] + k] *g_velg)
						/grain_velz[l];
				
					nb_e += (min_grain_charge[i] + k) *(y_data[nb_dch[i] + k]*(velg_mhd_i - g_velg) 
						+ ydot_data[nb_dch[i] + k]*(vel_i - grain_velz[l]));

					grain_velz_g[l] = g_velg;
				}
			}
		}
	}
	// For charge neutrality conservation:
	add_el_source = nb_e/fabs(add_el_source); // 
	ydot_data[network->e_nb] += nb_e/vel_i;
	
	// The derivative of the concentration of adsorbed chemical species;
	a = b = 0.;
	for (i = 0; i < nb_of_dust_comp; i++) {
		if (dust->get_grain_radius(i) > MIN_ADSORPTION_RADIUS) 
		{
			if (is_dust_charge_large[i]) 
			{
				l = nb_dch[i]-nb_dch[0];
				c = dust->get_grain_area(i)*y_data[nb_dch[i]] /grain_velz[l];
				a += c;
				b += grain_velz_g[l] *c; 
			}
			else
			{
				for (k = 0; k < max_grain_charge[i] - min_grain_charge[i] + 1; k++) 
				{
					l = nb_dch[i] - nb_dch[0] + k;
					c = dust->get_grain_area(i) *y_data[nb_dch[i]+k] /grain_velz[l];
					a += c;
					b += grain_velz_g[l] *c;
				}
			}
		}
	}
	a /= ads_dust_area;
	b /= ads_dust_area;

	for (i = nb_of_gas_species; i < nb_of_gas_species + nb_of_gmantle_species; i++) {
		ydot_data[i] = ydot_data[i]*a - y_data[i]*b;
	}

	// the derivative of dust temperature:
	for (i = 0; i < nb_of_dust_comp; i++) {
		ydot_data[nb_dct + i] = ydot_data[nb_dct + i]/av_grain_velz[i];
	}
	return returned_val;
}


/*	// here the excitation of H2 in formation on grains is taken into account:
#if (H2_FORMATION_EXCITATION)
	oh2_form_gaschem = oh2_form_grains = 0.; 	
	h2_prod -= h2_prod_gr;

	for (i = 0; i < nb_lev_h2; i++) 
	{
		a = h2_prod_gr *h2_excit_gf->get_efficiency(i);
		ydot_data[nb_of_species + i] += a;
		
		if (rounding(h2_di->lev_array[i].spin) == 1)
			oh2_form_grains += a;
	}
	
	h2_prod -= h2_prod_gas;
	for (i = 0; i < nb_lev_h2; i++) 
	{
		a = h2_prod_gas *h2_excit_gasph->get_efficiency(i, temp_n);
		ydot_data[nb_of_species + i] += a;

		if (rounding(h2_di->lev_array[i].spin) == 1)
			oh2_form_gaschem += a;
	}
#endif */
	
	// H2 dissociation
	
