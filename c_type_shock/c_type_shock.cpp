// 
/* Units:
Temperature is measured in K units; velocity - cm/s;
22.03.2017. Check for errors.
03.10.2017. Check for errors.
12.02.2018. Check for errors.
*/

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#if defined _WIN32 || defined _WIN64 
#include "stdafx.h"
#endif

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <list>
#include <cmath>
#include <ctime>
#include <limits>

#include "utils.h"
#include "integration.h"
#include "constants.h"
#include "parameters.h"
#include "chemistry.h"
#include "photoelectric_emission.h"
#include "chemical_reaction_data.h"
#include "spectroscopy.h"
#include "coll_rates.h"
#include "coll_rates_h2.h"
#include "coll_rates_h2o.h"
#include "coll_rates_co.h"
#include "coll_rates_ions.h"
#include "differential_equations.h"
#include "dust_model.h"
#include "radiation_field.h"
#include "reaction_rate_analyzer.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "c_type_shock.cpp"
using namespace std;

enum SHOCK_STATE_ID {
    SHOCK_STATE_NORMAL = 0,
    SHOCK_STATE_ION_SUPERSONIC,        // ion velocity changed the sign, ion fluid become supersonic
    SHOCK_STATE_NEUTRAL_SUBSONIC,      // neutral gas velocity changed the sign, neutral fluid become subsonic,
    SHOCK_STATE_H2_DISSOCIATION        // 99% of H2 is dissociated,
};

// abc de = abc * 10 ^ (de)
void reformat_floating_value(double x, int & abc, int & de);

// path to the file with cloud parameters must be given,
void assign_cloud_data(const string &path, evolution_data *user_data, vector<double> & iy, double & evol_time, 
	double & visual_extinct, double & cr_ioniz_rate, double & uv_field_strength, double & ir_field_strength, double & conc_h_tot, 
	int verbosity=1);

void calc_chem_evolution(const string &data_path, const string &output_path, double conc_h_tot, double op_ratio_h2,
    double visual_extinct, double cr_ioniz_rate, double uv_field_strength, double ir_field_strength, double c_abund_pah);
// path to data files, path to output for dark cloud simulations, path to output for shock simulations;
SHOCK_STATE_ID calc_shock(const string &data_path, const string &output_path1, const string &output_path2, double shock_vel,
	double magnetic_field, double c_abund_pah, double evol_time);

void calc_cr_dominated_region(const string &data_path, const string &output_path1, const string &output_path2, double c_abund_pah, 
	double evol_time, double cr_ir_factor, double incr_time);

void create_file_cloud_parameters(const string &output_path);
void create_file_specimen_abund(const string & output_path, const chem_network *network);
void create_file_ice_comp(const string & output_path);
void create_file_heating_rates(const string & output_path);
void create_file_energy_fluxes(const string & output_path);
void create_file_chem_hating(const string & output_path);
void create_file_dust_properties(const string & output_path, const dust_model *dust);
void create_file_reaction_rates(const string & output_path, const chem_network *network, bool save_disk_space = false);
void create_file_nautilus_data(const string & output_path);
void create_file_mol_data(const string & output_path, const evolution_data *user_data);
void create_file_h2_chemistry(const string & output_path);

void save_cloud_parameters(const evolution_data *user_data, const string &output_path, double ty, double visual_extinct, 
	double cr_ioniz_rate, double uv_field_strength, double ir_field_strength, double c_abund_pah, const N_Vector &y);

void save_specimen_abund(const string &output_path, int nb_of_species, const N_Vector &y, double conc_h_tot, double var);
void save_ice_comp(const string &output_path, const chem_network *network, const N_Vector &y, double var);
void save_heating_rates(const string &output_path, const evolution_data *user_data, double var);
void save_energy_fluxes(const string &output_path, const evolution_data *user_data, const N_Vector &y, double var, double dvar);
void save_chem_heating_rates(const string &output_path, const evolution_data *user_data, double ty);
void save_dust_properties(const string &output_path, const evolution_data *user_data, const N_Vector &y, double conc_h_tot, double var);
void save_reaction_rates(const string &output_path, const evolution_data *user_data, double var, double temp_n);
void save_nautilus_data(const string & output_path, double t, double av, double n, double gt, double dt);
void save_mol_data(const string & output_path, const evolution_data *user_data, const N_Vector &y, double var);
void save_file_h2_chemistry(const string & output_path, const evolution_data *user_data, const N_Vector &y, double var);

void create_file_mhd_vode(const string & output_path);
void save_mhd_vode(const string & input_path, const string & output_path, const chem_network *network, const N_Vector &y, double conc_h_tot, double var);

void cooling(const string &data_path);
void sputtering();

// function to print final statistics
static void print_stats(void *cvode_mem);

// function to check function return values
static int check_flag(void *flagvalue, char *funcname, int opt);

int main(int argc, char** argv)
{
    SHOCK_STATE_ID shock_state;
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, nb_processors;
	double conc_h_tot, op_ratio_h2, visual_extinct, shock_vel, magnetic_field, cr_ioniz_rate, c_abund_pah, uv_field_strength,
        ir_field_strength, ty, cr_ir_factor, incr_time, max_shock_speed;
	
	string data_path = "C:/input_data/";
	string mode, path, input_path, output_path;
	stringstream ss;
	ifstream input;

/*	ISRF_UV_Mathis1983 *rf1 = new ISRF_UV_Mathis1983();
	rf1->print_parameters(6.*EV_TO_CM_INVERSE, 13.6*EV_TO_CM_INVERSE);

	ISRF_UV_Draine1978 *rf2 = new ISRF_UV_Draine1978();
	rf2->print_parameters(6.*EV_TO_CM_INVERSE, 13.6*EV_TO_CM_INVERSE);

	CR_induced_UV_field *rf3 = new CR_induced_UV_field();
	rf3->print_parameters(7.09*EV_TO_CM_INVERSE, 14.6*EV_TO_CM_INVERSE);*/
	
//	h2_coll_data_process(data_path);
//	h2_coll_data_process_bossion(data_path);
//	reformat_h2o_coll_data(data_path, 0);
//	reformat_h2o_coll_data(data_path, 1);
//	merge_co_h_coll_data(data_path);
	

//	reformat_dust_files(data_path, "dust/draine_data/suvSil_81.txt", "dust/draine_data/dust_silicate_Draine2001.txt");
//	reformat_dust_files(data_path, "dust/draine_data/Gra_81.txt", "dust/draine_data/dust_graphite_Draine1993.txt");
//	reformat_dust_files_PAH(data_path, "dust/draine_data/PAHion_30.txt", "dust/draine_data/Gra_81.txt", "dust/draine_data/dust_PAHion_Draine2001.txt");
//	reformat_dust_files_PAH(data_path, "dust/draine_data/PAHneut_30.txt", "dust/draine_data/Gra_81.txt", "dust/draine_data/dust_PAHneut_Draine2001.txt");
//	reformat_chemical_data_Belloche2014(data_path);

//	path = "";
//	calc_grain_photoelectron_rates(data_path);
//	construct_gas_grain_reactions(data_path + "chemistry/UMIST_2012/surface_binding_energies_Penteado2017.txt", path);
//	construct_ion_recomb_grains(path);

//    path = "C:/Users/Александр/Александр/Данные и графики/paper Chemical evolution in molecular clouds in the vicinity of supernova remnants/";    
//    path += "output_data_2e4/dark_cloud_BEPent_B15A_DB035_QT_CR1-17_mult100/";
    path = "C:/Users/Александр/Documents/Данные и графики/paper C-type shocks - new data on H-H2 collisions/";
//    path += "output_data_2e4/dark_cloud_BEPent_B15A_DB035_QT_CR1-16/";
    path += "output_data_2e5/";
    production_routes(path, path + "shock_cr1-15_15/");

	path = "./output_data_2e4/dark_cloud_BEPent_B15A_DB035_QT_CR3-17/";
//	nautilus_comparison(path);
	
	input.open("input_parameters.txt"); 	
    while (!input.eof())
	{	
		do // comment lines are read:
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> nb_processors;
			
#ifdef _OPENMP
		omp_set_num_threads(nb_processors);
		
#pragma omp parallel 
	{
#pragma omp master 
		{
			cout << "OpenMP is supported" << endl;
			cout << "Nb of threads: " << omp_get_num_threads() << endl;
		}
	}
#endif
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> data_path;

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> mode;

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> c_abund_pah;
		
		// dark cloud parameters	
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> output_path;

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> cr_ioniz_rate;
		cr_ioniz_rate /= STANDARD_CR_IONIZ_RATE;

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> conc_h_tot;

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> visual_extinct;

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> uv_field_strength;

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> ir_field_strength;

        do
            input.getline(text_line, MAX_TEXT_LINE_WIDTH);
        while (text_line[0] == '#');

        ss.clear();
        ss.str(text_line);
        ss >> op_ratio_h2;

		if (mode == "DC") 
		{
			calc_chem_evolution(data_path, output_path, conc_h_tot, op_ratio_h2, visual_extinct, cr_ioniz_rate, uv_field_strength,
				ir_field_strength, c_abund_pah);
            break;
		}
		else 
		{
			do
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			while (text_line[0] == '#');

			ss.clear();
			ss.str(text_line);
			ss >> input_path;

			do
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			while (text_line[0] == '#');

			ss.clear();
			ss.str(text_line);
			ss >> output_path;

			do
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			while (text_line[0] == '#');

			ss.clear();
			ss.str(text_line);
			ss >> ty;

			do
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			while (text_line[0] == '#');

			ss.clear();
			ss.str(text_line);
			ss >> shock_vel;

			do
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			while (text_line[0] == '#');

			ss.clear();
			ss.str(text_line);
			ss >> magnetic_field;

			if (mode == "CS") 
			{ 
                cout << output_path << endl;
                shock_state = calc_shock(data_path, input_path, output_path, shock_vel, magnetic_field, c_abund_pah, ty); 
                cout << "   shock ended with the code " << (int) shock_state << endl;
                break;
            }
            else if (mode == "CS_") {
                // starting shock speed is given in the input data file,
                max_shock_speed = 120.1e+5; // 30.1e+5; 120.1e+5; for test simulations may be lower
                shock_state = SHOCK_STATE_NORMAL;
                for (i = 0; (shock_vel < max_shock_speed) && (shock_state == SHOCK_STATE_NORMAL); i++) {
                    ss.clear();
                    ss.str("");
                    ss << output_path;
                    
                    j = static_cast<int>(1.e-5*shock_vel + 0.1);
                    if (j < 10)
                        ss << "0";
                    ss << j << "/";

                    cout << left << setw(5) << i+1 << ss.str() << endl;
                    shock_state = calc_shock(data_path, input_path, ss.str(), shock_vel, magnetic_field, c_abund_pah, ty);
                    
                    cout << "   shock ended with the code " << (int) shock_state << endl;
                    shock_vel += 5.0e+5;
                }
                break;
            }
			else
			{
				do
					input.getline(text_line, MAX_TEXT_LINE_WIDTH);
				while (text_line[0] == '#');

				ss.clear();
				ss.str(text_line);
				ss >> input_path;

				do
					input.getline(text_line, MAX_TEXT_LINE_WIDTH);
				while (text_line[0] == '#');

				ss.clear();
				ss.str(text_line);
				ss >> output_path;

				do
					input.getline(text_line, MAX_TEXT_LINE_WIDTH);
				while (text_line[0] == '#');

				ss.clear();
				ss.str(text_line);
				ss >> ty;

				do
					input.getline(text_line, MAX_TEXT_LINE_WIDTH);
				while (text_line[0] == '#');

				ss.clear();
				ss.str(text_line);
				ss >> cr_ir_factor;

				do
					input.getline(text_line, MAX_TEXT_LINE_WIDTH);
				while (text_line[0] == '#');

				ss.clear();
				ss.str(text_line);
				ss >> incr_time; // in years;
			
				if (mode == "CR") 
				{
					calc_cr_dominated_region(data_path, input_path, output_path, c_abund_pah, ty, cr_ir_factor, incr_time);
                    break;
				}
				else {
					cout << " undefined mode";
					exit(1);
				}
			}
		}
	}
	return 0;
}

void calc_chem_evolution(const string &data_path, const string &output_path, double conc_h_tot, double op_ratio_h2, 
    double visual_extinct, double cr_ioniz_rate, double uv_field_strength, double ir_field_strength, double c_abund_pah)
{
#ifdef __linux__
	stringstream lin_out;	
	lin_out << output_path;
	lin_out << "out";
	lin_out << "_screen";
	lin_out << ".txt";
	// lin_out.str("/dev/null");

	ofstream outerr(lin_out.str().c_str(), ios::app);
	streambuf *orig_cerr = cerr.rdbuf(outerr.rdbuf());
	
	ofstream out(lin_out.str().c_str(), ios::app);
	streambuf *orig_cout = cout.rdbuf(out.rdbuf());
#endif	
	
	bool is_new_chd;
	int i, flag, nb, nb_of_species, nb_lev_h2o, nb_lev_h2, nb_lev_co, nb_vibr_h2o, nb_vibr_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3,
		nb_vibr_ch3oh, nb_lev_ch3oh, nb_lev_ci, nb_lev_cii, nb_lev_oi, nb_of_dust_comp, nb_of_grain_charges, nb_of_equat, 
		nb_dct, nb_mhd, verbosity;
	long int nb_steps;
	double a, init_temp, conc_e, h2_form_const, t, ty, tfin, tout, rel_tol, tmult, ion_conc, ion_pah_conc, ion_dens, ion_pah_dens, 
        ion_dust_dens, op_h2_ratio;
	
	double *chem_abund(0);
	vector<double> new_y, init_ch;

	string fname;
	ofstream output;
	time_t timer;
	
	cout << scientific;
	cout.precision(3);

	verbosity = 1;
	// initial temperature in K:
	init_temp = 10.;

	// Spectroscopic parameters for H2, H2O and CO molecule:
	// max nb of ortho- and para-H2 levels for which data are given by Wrathmall et al. (2007) - 109;
	nb_lev_h2 = 100; 	
	nb_vibr_h2o = 0;
	nb_lev_h2o = 45; // must be <= 52, otherwise vibrational state levels must be taken into acount;
	nb_vibr_co = 0;
	nb_lev_co = 30; // number of levels lower than first vibrationally excited level is 33;
	
	nb_vibr_ch3oh = 0;
// the number of levels in cold cloud simulations should be no higher than for shock simulations, check it
// may be, it is not necessary to take in to account level nb > 1 for species below 
#if (CALCULATE_POPUL_METHANOL)
	nb_lev_ch3oh = 100;
#else
	nb_lev_ch3oh = 1;
#endif

#if (CALCULATE_POPUL_NH3_OH)
    nb_lev_onh3 = 9; // ortho-NH3: He coll data - 22, H2 coll data - 17
    nb_lev_pnh3 = 16; // para-NH3: He coll data - 16, H2 coll data - 34
    nb_lev_oh = 10; // OH: He coll data - 44, H2 coll data - 20 (without HF splitting)
#else
    nb_lev_onh3 = 1; 
    nb_lev_pnh3 = 1;
    nb_lev_oh = 1;
#endif

	timer = time(NULL);
	cout << ctime(&timer) << "Chemical evolution of static cloud is simulated" << endl;

	chemistry_evolution_data user_data(data_path, output_path, nb_lev_h2, nb_vibr_h2o, nb_lev_h2o, nb_vibr_co, nb_lev_co, 
        nb_lev_pnh3, nb_lev_onh3, nb_lev_oh, nb_vibr_ch3oh, nb_lev_ch3oh, c_abund_pah, verbosity);

	user_data.set_parameters(visual_extinct, cr_ioniz_rate, uv_field_strength, ir_field_strength);
	
	const chem_network *network 
		= user_data.get_network();
	
	const dust_model *dust 
		= user_data.get_dust();

	nb_of_dust_comp = dust->nb_of_comp;
	nb_of_species = user_data.get_nb_of_species();
	// level numbers must be updated:
	user_data.get_nb_of_levels(nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh,
		nb_lev_ci, nb_lev_oi, nb_lev_cii);

	// Initialization of the initial chemical abundances.
	// Note! Concentrations of molecules and atoms, for which the equations for level populations are integrated, must be > 0.;
	chem_abund = new double [nb_of_species];
	fname = data_path + "chemistry/initial_abund_hincelin2011.txt";
	
	init_chem_abund(fname, network, chem_abund);
	conc_e = chem_abund[network->e_nb]*conc_h_tot;

	user_data.recalc_grain_charge_ranges(init_temp, conc_e, conc_h_tot, init_ch);
	user_data.get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
	
	t = 0.;
	tfin = 1.01e+8*YEARS_TO_SECONDS;
	tout = YEARS_TO_SECONDS;
	tmult = pow(10, 0.03125); // 1/16 = 0.0625; 1/32 = 0.03125
	
	N_Vector y, ydot, abs_tol;
	SUNMatrix A(NULL);
	SUNLinearSolver LS(NULL);

	y = N_VNew_Serial(nb_of_equat);
	ydot = N_VNew_Serial(nb_of_equat);
	abs_tol = N_VNew_Serial(nb_of_equat);

	for (i = 0; i < nb_of_equat; i++) {
		NV_Ith_S(y, i) = 0.;
	}
	
	// The definition of initial values of specimen concentrations:
	for (i = 0; i < nb_of_species; i++) 
	{
        NV_Ith_S(y, i) = chem_abund[i]*conc_h_tot;
		if (NV_Ith_S(y, i) < ABS_CONCENTRATION_ERROR_SOLVER)
			NV_Ith_S(y, i) = ABS_CONCENTRATION_ERROR_SOLVER;
	}

	// Initialization of initial abundances of molecules and atoms for which the level populations are calculated, 
	// ortho/para ratio for H2 is arbitrary;
	NV_Ith_S(y, nb_of_species) = 1./(op_ratio_h2 + 1.) *NV_Ith_S(y, network->h2_nb);
	NV_Ith_S(y, nb_of_species + 1) = op_ratio_h2/(op_ratio_h2 + 1.) *NV_Ith_S(y, network->h2_nb);
    nb = nb_of_species + nb_lev_h2; 
	
	// ortho/para ratio for H2O is 3;
	NV_Ith_S(y, nb) = 0.25*NV_Ith_S(y, network->h2o_nb);
	nb += nb_lev_h2o;

	NV_Ith_S(y, nb) = 0.75*NV_Ith_S(y, network->h2o_nb);
	nb += nb_lev_h2o;

	NV_Ith_S(y, nb) = NV_Ith_S(y, network->co_nb);
	nb += nb_lev_co;

	NV_Ith_S(y, nb) = NV_Ith_S(y, network->oh_nb);
	nb += nb_lev_oh;

	// ortho/para statistical weights are 2, but paraNH3 has two times larger number of levels;
	NV_Ith_S(y, nb) = 0.5*NV_Ith_S(y, network->nh3_nb);
	nb += nb_lev_pnh3;

	NV_Ith_S(y, nb) = 0.5*NV_Ith_S(y, network->nh3_nb);
	nb += nb_lev_onh3;

	// A-/E- methanol ratio is assumed 1;
    // see Holdship et al. ApJ 880, p. 138 (2019) 
	NV_Ith_S(y, nb) = 0.5*NV_Ith_S(y, network->ch3oh_nb);
	nb += nb_lev_ch3oh;

	NV_Ith_S(y, nb) = 0.5*NV_Ith_S(y, network->ch3oh_nb);
	nb += nb_lev_ch3oh;

	NV_Ith_S(y, nb) = NV_Ith_S(y, network->ci_nb);
	nb += nb_lev_ci;

	NV_Ith_S(y, nb) = NV_Ith_S(y, network->oi_nb);
	nb += nb_lev_oi;

	NV_Ith_S(y, nb) = NV_Ith_S(y, network->cii_nb);
	nb += nb_lev_cii;
	
	// initialization of dust charges:
	for (i = 0; i < nb_of_grain_charges; i++) {
		NV_Ith_S(y, nb + i) = init_ch[i];
	}

	// initialization of initial dust/gas temperatures:
	for (i = 0; i < nb_of_dust_comp; i++) {
		NV_Ith_S(y, nb_dct + i) = init_temp;
	}
	NV_Ith_S(y, nb_mhd) = NV_Ith_S(y, nb_mhd + 1) = NV_Ith_S(y, nb_mhd + 2) = init_temp;

	a = user_data.calc_total_gas_charge(y) + user_data.calc_total_grain_charge(y);
	NV_Ith_S(y, network->e_nb) += a;

	if (NV_Ith_S(y, network->e_nb) < 0.) {
		cout << "Error during initialization: electron concentration is negative!";
		exit(1);
	}

    // Initialization for tolerances:
    rel_tol = REL_ERROR_SOLVER;
    user_data.set_tolerances(abs_tol);

	// Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration 
	void *cvode_mem = CVodeCreate(CV_BDF);

	// Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(t,y), 
	// the inital time t0, and the initial dependent variable vector y;
	flag = CVodeInit(cvode_mem, f_chem, t, y);

	// Call CVodeSVtolerances to specify the scalar tolerances:
	flag = CVodeSVtolerances(cvode_mem, rel_tol, abs_tol);

	// The maximal number of steps between simulation stops;
	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);

	// specifies the maximum number of error test failures permitted in attempting one step:
	flag = CVodeSetMaxErrTestFails(cvode_mem, MAX_ERR_TEST_FAILS_SOLVER); // default value is 7; 
    flag = CVodeSetMaxConvFails(cvode_mem, MAX_CONV_FAILS_SOLVER); // default value is 10;

	// Create dense SUNMatrix for use in linear solves 
	A = SUNDenseMatrix(nb_of_equat, nb_of_equat);
	
	// Create dense SUNLinearSolver object for use by CVode 
	LS = SUNDenseLinearSolver(y, A);

	// Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode 
	flag = CVDlsSetLinearSolver(cvode_mem, LS, A);

	// Call CVDense to specify the CVDENSE dense linear solver:
	// flag = CVDense(cvode_mem, nb_of_equat);

	// The function attaches the user data block to the solver;
	flag = CVodeSetUserData(cvode_mem, &user_data);
	
	create_file_cloud_parameters(output_path);
	create_file_heating_rates(output_path);
	create_file_chem_hating(output_path);
	create_file_specimen_abund(output_path, network);
	create_file_ice_comp(output_path);
	create_file_dust_properties(output_path, dust);
	create_file_reaction_rates(output_path, network);
	create_file_nautilus_data(output_path);
	create_file_mhd_vode(output_path);

#if (SAVE_RADIATIVE_FACTORS)
	user_data.create_file_radiative_transfer(output_path);
#endif

	fname = output_path + "sim_phys_param.txt";
	output.open(fname.c_str());
	output << "! p1 - total abundance of grain mantle molecules, [molecules/H]" << endl
		<< "! p2 - abundance of hydrocarbon molecules in ice mantles, CnHm, n >=2,  [molecules/H]" << endl
		<< "! p3 - parameter of H2 formation on grains (reaction *H+*H), rate = a*n_H_tot*n_H, a [cm3 s-1]" << endl
		<< "! p4 - electron concentration, [cm-3]" << endl
		<< "! p5 - ion concentration (without PAHs and small grains), [cm-3]" << endl
		<< "! p6 - total electric charge of grains, [cm-3]" << endl;

	output << "!";
	for (i = 0; i < 11; i++) {
		output << left << setw(14) << i + 1;
	}
	output << endl << left << setw(15) << "!time(yrs)" << setw(14) << "T_n" << setw(14) << "T_i" << setw(14) << "T_e" 
		<< setw(14) << "oph2" << setw(14) << "p1" << setw(14) << "p2" << setw(14) << "p3" << setw(14) << "p4" 
        << setw(14) << "p5" << setw(14) << "p6" << endl;
	output.close();
	
	nb = 0;
	timer = time(NULL);
	while (tout < tfin)
	{
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL); // CV_NORMAL or CV_ONE_STEP
		if (flag != CV_SUCCESS) {
			cout << "Error ocurred in solver: " << flag << endl;
			break;
		}
		flag = CVodeGetNumSteps(cvode_mem, &nb_steps);
		check_flag(&flag, "CVodeGetNumSteps", 1);
		
		ty = t/YEARS_TO_SECONDS;
		if (verbosity) {
			cout.precision(3);
			cout << left << "cloud age (yr): " << setw(12) << ty << "calc time (s): " 
				<< setw(8) << (int) (time(NULL) - timer) << "nb of steps: " << nb_steps << endl;
		}

		// the function is called in order to update user data:
		user_data.f(t, y, ydot); 
		
		save_specimen_abund(output_path, nb_of_species, y, conc_h_tot, ty);
		save_ice_comp(output_path, network, y, ty);
		save_heating_rates(output_path, &user_data, ty);
		save_dust_properties(output_path, &user_data, y, conc_h_tot, ty);
		save_nautilus_data(output_path, ty, visual_extinct, conc_h_tot, NV_Ith_S(y, nb_mhd), user_data.get_av_dust_temp());
		
		user_data.calc_ion_dens(y, ion_conc, ion_pah_conc, ion_dens, ion_pah_dens, ion_dust_dens);
		h2_form_const = user_data.get_h2_form_grains()/(conc_h_tot *NV_Ith_S(y, network->h_nb));
        op_h2_ratio = NV_Ith_S(y, network->h2_nb) / user_data.calc_conc_ph2(y) - 1.;

		fname = output_path + "sim_phys_param.txt";
		output.open(fname.c_str(), ios::app);
		output << scientific;

		output.precision(7);
		output << left << setw(15) << ty;

		output.precision(5);		
		output << left << setw(14) << NV_Ith_S(y, nb_mhd) << setw(14) << NV_Ith_S(y, nb_mhd + 1) 
			<< setw(14) << NV_Ith_S(y, nb_mhd + 2) << setw(14) << op_h2_ratio << setw(14) << user_data.calc_ice_conc(y)/conc_h_tot
			<< setw(14) << user_data.calc_hydrocarbon_conc(y)/conc_h_tot << setw(14) << h2_form_const 
			<< setw(14) << NV_Ith_S(y, network->e_nb) << setw(14) << ion_conc
			<< setw(14) << user_data.calc_total_grain_charge(y) << endl;
		output.close();
		
		if (ty > 0.99*1.e+4 && ty < 0.99e+7) 
		{
			save_cloud_parameters(&user_data, output_path, ty, visual_extinct, cr_ioniz_rate, uv_field_strength, 
				ir_field_strength, c_abund_pah, y);
			save_chem_heating_rates(output_path, &user_data, ty);
			save_mhd_vode(data_path, output_path, network, y, conc_h_tot, ty);

#if (SAVE_RADIATIVE_FACTORS)
			user_data.save_radiative_transfer_data(output_path, ty);
#endif
		}
		if (nb%2 == 0) {
			save_reaction_rates(output_path, &user_data, ty, NV_Ith_S(y, nb_mhd));
		}

		is_new_chd = user_data.recalc_grain_charge_ranges(y, new_y);
		if (is_new_chd) 
        {	
			user_data.get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
			
			N_VDestroy_Serial(y);
			N_VDestroy_Serial(ydot);
            N_VDestroy_Serial(abs_tol);

			y = N_VNew_Serial(nb_of_equat);
			ydot = N_VNew_Serial(nb_of_equat);
            abs_tol = N_VNew_Serial(nb_of_equat);

			for (i = 0; i < nb_of_equat; i++) {
				NV_Ith_S(y, i) = new_y[i];
			}

			SUNLinSolFree(LS);
			SUNMatDestroy(A);
			CVodeFree(&cvode_mem);
			
			cvode_mem = CVodeCreate(CV_BDF);
			flag = CVodeInit(cvode_mem, f_chem, t, y);
            
            user_data.set_tolerances(abs_tol);
			flag = CVodeSVtolerances(cvode_mem, rel_tol, abs_tol);
			
			A = SUNDenseMatrix(nb_of_equat, nb_of_equat);
			LS = SUNDenseLinearSolver(y, A);
			flag = CVDlsSetLinearSolver(cvode_mem, LS, A);

			flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
			flag = CVodeSetMaxErrTestFails(cvode_mem, MAX_ERR_TEST_FAILS_SOLVER);
            flag = CVodeSetMaxConvFails(cvode_mem, MAX_CONV_FAILS_SOLVER);
			flag = CVodeSetUserData(cvode_mem, &user_data);
		}
		nb++;
		tout = tout*tmult;
	}
	print_stats(cvode_mem);

	// Free memory;
    CVodeFree(&cvode_mem);
	SUNLinSolFree(LS);
	SUNMatDestroy(A);
	
	N_VDestroy_Serial(y);
	N_VDestroy_Serial(ydot);
    N_VDestroy_Serial(abs_tol);
	delete [] chem_abund;
}

void assign_cloud_data(const string &path, evolution_data *user_data, vector<double> & iy, double & evol_time, double & visual_extinct, 
	double & cr_ioniz_rate, double & uv_field_strength, double & ir_field_strength, double & conc_h_tot, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];	
	int i, j, nb, nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, nb_lev_oi, nb_lev_ci, nb_lev_cii, 
		nb_of_species, nb_of_equat, nb_of_dust_comp, nb_of_grain_charges, nb_dct, nb_mhd;
	double a, time1, time2, g_conc, temp_n, temp_i, temp_e, g_ch, c_abund_pah;
	
	bool *is_av_ch(0);
	int *min_gch(0), *max_gch(0);
	double *temp_d(0);
	
	string fname, sn;
	vector<double> ch_arr, chem_abund, h2_popul_dens, ph2o_popul_dens, oh2o_popul_dens, co_popul_dens, oh_popul_dens, 
		pnh3_popul_dens, onh3_popul_dens, ch3oh_a_popul_dens, ch3oh_e_popul_dens, ci_popul_dens, cii_popul_dens, oi_popul_dens;
	vector<string> sp_names;
	ifstream input;

	fname = path + "sim_cloud_data.txt";
	input.open(fname.c_str(), ios::in);

	time1 = -1;
	while ((time1 < evol_time) && (!input.eof()))
	{		
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> time2; // reading evolution age of the cloud;

		if (time2 > evol_time && time1 > 0.) {
			if ((evol_time - time1) < (time2 - evol_time))
				break;
		}
		time1 = time2;

		ch_arr.clear();
		sp_names.clear();
		chem_abund.clear();
		h2_popul_dens.clear();
		ph2o_popul_dens.clear();
		oh2o_popul_dens.clear();
		co_popul_dens.clear();
		oh_popul_dens.clear();
		pnh3_popul_dens.clear();
		onh3_popul_dens.clear();
		ch3oh_a_popul_dens.clear();
		ch3oh_e_popul_dens.clear();
		ci_popul_dens.clear();
		cii_popul_dens.clear();
		oi_popul_dens.clear();

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> visual_extinct >> cr_ioniz_rate >> uv_field_strength >> ir_field_strength >> conc_h_tot >> c_abund_pah;

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> temp_n >> temp_i >> temp_e;

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb_of_dust_comp;
	
		if (temp_d == 0) 
		{	
			temp_d = new double [nb_of_dust_comp];
			is_av_ch = new bool [nb_of_dust_comp];
	
			min_gch = new int [nb_of_dust_comp];
			max_gch = new int [nb_of_dust_comp];
		}
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		for (i = 0; i < nb_of_dust_comp; i++) {
			input >> nb >> sn >> temp_d[i] >> g_conc >> g_ch >> j;
			if (j == 0) { // average charge;
				ch_arr.push_back(g_conc);
				ch_arr.push_back(g_ch);
				
				is_av_ch[i] = true;
				min_gch[i] = max_gch[i] = 0;
			}
			else {
				input >> min_gch[i] >> max_gch[i];
				for (j = 0; j < max_gch[i] - min_gch[i] + 1; j++) 
				{
					input >> a;
					ch_arr.push_back(a);	
				}
				is_av_ch[i] = false;
			}
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb; // nb of species in the file list;
		
		for (i = 0; i < nb; i++) 
		{
			input >> j >> sn >> a;
			sp_names.push_back(sn);
			chem_abund.push_back(a);
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;
		
		for (i = 0; i < nb; i++) {
			input >> j >> a;
			h2_popul_dens.push_back(a);
		}
		
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;
		
		for (i = 0; i < nb; i++) {
			input >> j >> a;
			ph2o_popul_dens.push_back(a);
		}
		
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;
		
		for (i = 0; i < nb; i++) {
			input >> j >> a;
			oh2o_popul_dens.push_back(a);
		}
		
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;

		for (i = 0; i < nb; i++) {
			input >> j >> a;
			co_popul_dens.push_back(a);
		}
		
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;

		for (i = 0; i < nb; i++) {
			input >> j >> a;
			oh_popul_dens.push_back(a);
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;

		for (i = 0; i < nb; i++) {
			input >> j >> a;
			pnh3_popul_dens.push_back(a);
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;

		for (i = 0; i < nb; i++) {
			input >> j >> a;
			onh3_popul_dens.push_back(a);
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;

		for (i = 0; i < nb; i++) {
			input >> j >> a;
			ch3oh_a_popul_dens.push_back(a);
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;

		for (i = 0; i < nb; i++) {
			input >> j >> a;
			ch3oh_e_popul_dens.push_back(a);
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;
		
		for (i = 0; i < nb; i++) {
			input >> j >> a;
			ci_popul_dens.push_back(a);
		}
		
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;
		
		for (i = 0; i < nb; i++) {
			input >> j >> a;
			oi_popul_dens.push_back(a);
		}
			
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb;
		
		for (i = 0; i < nb; i++) {
			input >> j >> a;
			cii_popul_dens.push_back(a);
		}
		// reading the end of the line:
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	input.close();

	evol_time = time1;
	if (verbosity) {
		cout << "The data on cloud parameters are read, evolution time " << evol_time << endl;
	}

	const dust_model *dust 
		= user_data->get_dust();
	
	const chem_network *network 
		= user_data->get_network();
	
	if (dust->nb_of_comp != nb_of_dust_comp) {
		cout << "Error in data: nbs of dust components differ in file and code;";
		exit(1);
	}
	
	user_data->set_parameters(visual_extinct, cr_ioniz_rate, uv_field_strength, ir_field_strength);
	user_data->set_grain_charge_ranges(is_av_ch, min_gch, max_gch, (int) ch_arr.size());
	
	nb_of_species = user_data->get_nb_of_species();
	user_data->get_nb_of_levels(nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, 
		nb_lev_ci, nb_lev_oi, nb_lev_cii);
	user_data->get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);

	iy.clear();
	for (i = 0; i < nb_of_equat; i++) {
		iy.push_back(0.);
	}
	for (i = 0; i < (int) chem_abund.size(); i++) {
		j = network->find_specimen(sp_names[i]);
		if (j >= 0)
			iy[j] = chem_abund[i];
	}
	// The nb of levels of a specimen in the file data may be != than value given in user_data 
	for (i = 0; i < (int) h2_popul_dens.size() && i < nb_lev_h2; i++) {
		iy[nb_of_species + i] = h2_popul_dens[i];
	}
	nb = nb_of_species + nb_lev_h2;
	
	for (i = 0; i < (int) ph2o_popul_dens.size() && i < nb_lev_h2o; i++) {
		iy[nb + i] = ph2o_popul_dens[i];
	}
	nb += nb_lev_h2o;
	
	for (i = 0; i < (int) oh2o_popul_dens.size() && i < nb_lev_h2o; i++) {
		iy[nb + i] = oh2o_popul_dens[i];
	}
	nb += nb_lev_h2o;
	
	for (i = 0; i < (int) co_popul_dens.size() && i < nb_lev_co; i++) {
		iy[nb + i] = co_popul_dens[i];
	}
	nb += nb_lev_co;
	
	for (i = 0; i < (int) oh_popul_dens.size() && i < nb_lev_oh; i++) {
		iy[nb + i] = oh_popul_dens[i];
	}
	nb += nb_lev_oh;

	for (i = 0; i < (int) pnh3_popul_dens.size() && i < nb_lev_pnh3; i++) {
		iy[nb + i] = pnh3_popul_dens[i];
	}
	nb += nb_lev_pnh3;

	for (i = 0; i < (int) onh3_popul_dens.size() && i < nb_lev_onh3; i++) {
		iy[nb + i] = onh3_popul_dens[i];
	}
	nb += nb_lev_onh3;

	for (i = 0; i < (int) ch3oh_a_popul_dens.size() && i < nb_lev_ch3oh; i++) {
		iy[nb + i] = ch3oh_a_popul_dens[i];
	}
	nb += nb_lev_ch3oh;

	for (i = 0; i < (int) ch3oh_e_popul_dens.size() && i < nb_lev_ch3oh; i++) {
		iy[nb + i] = ch3oh_e_popul_dens[i];
	}
	nb += nb_lev_ch3oh;
	
	for (i = 0; i < (int) ci_popul_dens.size() && i < nb_lev_ci; i++) {
		iy[nb + i] = ci_popul_dens[i];
	}
	nb += nb_lev_ci;
	
	for (i = 0; i < (int) oi_popul_dens.size() && i < nb_lev_oi; i++) {
		iy[nb + i] = oi_popul_dens[i];
	}
	nb += nb_lev_oi;
	
	for (i = 0; i < (int) cii_popul_dens.size() && i < nb_lev_cii; i++) {
		iy[nb + i] = cii_popul_dens[i];
	}
	nb += nb_lev_cii;
	
	for (i = 0; i < nb_of_grain_charges; i++) {
		iy[nb + i] = ch_arr[i];
	}
	
	for (i = 0; i < nb_of_dust_comp; i++) {
		iy[nb_dct + i] = temp_d[i];
	}

	iy[nb_mhd] = temp_n;
	iy[nb_mhd + 1] = temp_i;
	iy[nb_mhd + 2] = temp_e;

	delete [] is_av_ch;
	delete [] min_gch;
	delete [] max_gch;
}

SHOCK_STATE_ID calc_shock(const string &data_path, const string &output_path1, const string &output_path2, double shock_vel,
	double magnetic_field, double c_abund_pah, double evol_time)
{
#ifdef __linux__
	stringstream lin_out;	
	lin_out << output_path2;
	lin_out << "out";
	lin_out << "_screen";
	lin_out << ".txt";
	// lin_out.str("/dev/null");

	ofstream outerr(lin_out.str().c_str(), ios::app);
    streambuf *orig_cerr = cerr.rdbuf(); // original cerr;
    cerr.rdbuf(outerr.rdbuf());
	
	ofstream out(lin_out.str().c_str(), ios::app);
    streambuf *orig_cout = cout.rdbuf();
    cout.rdbuf(out.rdbuf());
#endif	

    SHOCK_STATE_ID shock_state = SHOCK_STATE_NORMAL;
	bool is_post_shock, must_be_stopped, is_new_chd, is_new_vg, save_disk_space;
	int i, nb_saved, nb_not_saved, nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_vibr_h2o, nb_vibr_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_vibr_ch3oh, 
		nb_lev_ch3oh, nb_lev_oi, nb_lev_ci, nb_lev_cii, nb_of_species, nb_of_equat, nb_of_grain_charges, nb_dct, nb_mhd, flag, 
		nb_saved_cloud_param, verbosity;
	long int tot_nb_steps;
	double a, b, visual_extinct, cr_ioniz_rate, uv_field_strength, ir_field_strength, magn_precursor_length, magn_sonic_speed, 
		sound_speed, conc_h_tot, temp_n, temp_i, temp_e, neut_dens, ion_dens, neut_conc, ion_conc, ion_pah_conc, ion_pah_dens, 
        ion_dust_dens, rel_tol, ty, z, zout, zfin, dz, vel_n_grad, vel_i_grad, h2_form_const, dvel_shock_stop, z_saved, dv_to_v_lim;
	double *prev_y(0);
	
	string fname, sn;
	ofstream output;
	time_t timer;
	vector<double> new_y, veln_arr, veli_arr, dz_arr;
	
	SUNMatrix A(NULL);
	SUNLinearSolver LS(NULL);
	N_Vector y, abs_tol;

	verbosity = 1;		
	timer = time(NULL);
	cout << ctime(&timer) << "Shock wave is simulated" << endl;

	nb_lev_h2 = 298; // the maximal number for which Einstein coefficients are provided - 298 levels,
	nb_vibr_h2o = 1;
	nb_lev_h2o = 150;
	nb_vibr_co = 0;
	nb_lev_co = 41; // 41
	
#if (CALCULATE_POPUL_METHANOL)
	nb_vibr_ch3oh = 1;
	nb_lev_ch3oh = 300;
#else
	nb_vibr_ch3oh = 0;
	nb_lev_ch3oh = 1;
#endif

#if (CALCULATE_POPUL_NH3_OH)
	nb_lev_onh3 = 17; // ortho-NH3: He coll data - 22, H2 coll data - 17
	nb_lev_pnh3 = 34; // para-NH3: He coll data - 16, H2 coll data - 34
	nb_lev_oh = 24;  // OH: He coll data - 44, H2 coll data - 20 (without HF splitting); 24 (H2), 56 (He) - with HF splitting
#else
    nb_lev_onh3 = 1; 
    nb_lev_pnh3 = 1;
    nb_lev_oh = 1;
#endif
	// initialization of the data necessary for differential equation integration:
	mhd_shock_data user_data(data_path, output_path2, nb_lev_h2, nb_vibr_h2o, nb_lev_h2o, nb_vibr_co, nb_lev_co, 
        nb_lev_pnh3, nb_lev_onh3, nb_lev_oh, nb_vibr_ch3oh, nb_lev_ch3oh, c_abund_pah, verbosity);

	nb_of_species = user_data.get_nb_of_species();
	user_data.get_nb_of_levels(nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, 
		nb_lev_ci, nb_lev_oi, nb_lev_cii);
	
	// reading initial data for simulations from file, output_path1 - path to the file,
	// data are saved to vector new_y, 
	// physical parameters and parameters of grain charge distribution are assigned in user_data, 
	assign_cloud_data(output_path1, &user_data, new_y, evol_time, visual_extinct, cr_ioniz_rate, uv_field_strength, 
		ir_field_strength, conc_h_tot, verbosity);
	user_data.get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
	
	y = N_VNew_Serial(nb_of_equat);
    abs_tol = N_VNew_Serial(nb_of_equat);
	prev_y = new double [nb_of_equat];

	const dust_model *dust 
		= user_data.get_dust();
	
	const chem_network *network 
		= user_data.get_network();

	user_data.set_magn_field(magnetic_field);
	user_data.set_shock_vel(shock_vel);
	
	for (i = 0; i < nb_of_equat-2; i++) {
		NV_Ith_S(y, i) = new_y[i];
	}
	NV_Ith_S(y, nb_mhd + 3) = shock_vel;
	NV_Ith_S(y, nb_mhd + 4) = 0.9995*shock_vel; // 
	
	temp_n = NV_Ith_S(y, nb_mhd);
	temp_i = NV_Ith_S(y, nb_mhd + 1);
	temp_e = NV_Ith_S(y, nb_mhd + 2);
		
	a = user_data.calc_total_gas_charge(y) + user_data.calc_total_grain_charge(y);
	NV_Ith_S(y, network->e_nb) += a;

	// calculation of the mass and number densities for the neutral and ion fluids;
	user_data.calc_ion_dens(y, ion_conc, ion_pah_conc, ion_dens, ion_pah_dens, ion_dust_dens);
	user_data.calc_neutral_dens(y, neut_conc, neut_dens);

	// Assesment of the magnetic precursor length - the distance that characterizes the decrease of the ion velocity, 
	// Draine, ApJ 241, p. 1021, 1980; Flower & Pineau des Forets, MNRAS 275, p. 1049, 1995; 
	// 'Langevin' cross section is approximately equal to sigma*vel = 2.e-9 cm3 s-1;
	magn_precursor_length = magnetic_field*magnetic_field /(M_PI*neut_dens*shock_vel)
		/(2.e-9*ion_conc + 0.5*dust->area_perH *conc_h_tot *shock_vel);

	// magnetosonic speed, see Draine, ApJ 241, p. 1021, 1980;
	a = magnetic_field*magnetic_field/(4*M_PI);
	magn_sonic_speed = sqrt((a + 5./3.*ion_conc*BOLTZMANN_CONSTANT*(temp_i + temp_e))
			/(ion_dens + ion_pah_dens + a/(SPEED_OF_LIGHT*SPEED_OF_LIGHT)));

	// neutral gas sound speed
	sound_speed = sqrt(5.*neut_conc*BOLTZMANN_CONSTANT*temp_n/(3.*neut_dens));
	 
	if (verbosity) {
        cout << scientific;
        cout.precision(3);
        
        cout << "Characteristic scale of magnetic precursor (cm): " << magn_precursor_length << endl
			<< "Characteristic value of ion velocity gradient (cm/s/cm): " << -shock_vel/magn_precursor_length << endl
			<< "Magnetosonic speed (cm/s): " << magn_sonic_speed << endl
			<< "Neutral gas sound speed (cm/s): " << sound_speed << endl 
			<< "Ion Alfven speed (cm/s): " << magnetic_field/sqrt(4*M_PI*(ion_dens + ion_pah_dens)) << endl
			<< "Alfven speed (cm/s): " << magnetic_field/sqrt(4*M_PI*neut_dens) << endl
			<< "Ionization fraction: " << ion_conc/neut_conc << endl;
	}
	for (i = 0; i < nb_of_equat; i++) {
		prev_y[i] = NV_Ith_S(y, i);
	}

	rel_tol = REL_ERROR_SOLVER;
	user_data.set_tolerances(abs_tol);
	
    ty = z = z_saved = zout = 0.;
    dz = 0.01*magn_precursor_length; // cm
    zfin = 1000.*magn_precursor_length;
    dv_to_v_lim = 0.01; 

    // relative difference between ion and neutral speeds, at which shock stops;
    dvel_shock_stop = 0.03; // for studies of chemical evolution of post-shock gas, set 0.001
    
	// Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration 
	void *cvode_mem = CVodeCreate(CV_BDF);

	// Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(t,y), 
	// the inital time z, and the initial dependent variable vector y:
	flag = CVodeInit(cvode_mem, f_mhd, z, y);

	// Call CVodeSVtolerances to specify the scalar tolerances:
	flag = CVodeSVtolerances(cvode_mem, rel_tol, abs_tol);

	// Create dense SUNMatrix for use in linear solves
	A = SUNDenseMatrix(nb_of_equat, nb_of_equat);
	
	// Create dense SUNLinearSolver object for use by CVode 
	LS = SUNDenseLinearSolver(y, A);

	// Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode 
	flag = CVDlsSetLinearSolver(cvode_mem, LS, A);

	// maximum number of steps between simulation stops:
	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);

	// specifies the maximum number of error test failures permitted in attempting one step:
	flag = CVodeSetMaxErrTestFails(cvode_mem, MAX_ERR_TEST_FAILS_SOLVER); // default value is 7; 
    flag = CVodeSetMaxConvFails(cvode_mem, MAX_CONV_FAILS_SOLVER); // default value is 10;

	// function attaches the user data block to the solver:
	flag = CVodeSetUserData(cvode_mem, &user_data);
		
	fname = output_path2 + "sim_phys_param.txt";
	output.open(fname.c_str());
	
	output << "! nh - concentration of H nuclei, [cm-3]" << endl
		<< "! aic- abundance of molecules adsorbed on grains, mol/H" << endl
		<< "! ae - abundance of electrons, e/H" << endl
		<< "! aio- abundance of ions, ions/H" << endl
		<< "! vgn- velocity gradient of neutral component - averaged based on the fixed dv/v" << endl
        << "!   and instantaneous from MHD equations, [cm/s/cm]" << endl
        << "! vgi- velocity gradient of ion component, [cm/s/cm]" << endl
		<< "! p1 - parameter of H2 formation on grains (reaction *H+*H), rate = a*n_H_tot*n_H, a [cm3 s-1]" << endl
		<< "! p2 - is equal to add_Ne/Ne, add_Ne is the term added to electron production rate Ne to compensate " << endl
		<< "!	the velocity difference between ions and charge grains," << endl
		<< "! p3 - total electric charge of grains, [cm-3]" << endl;

	output << "!";
	for (i = 0; i < 18; i++) {
		output << left << setw(14) << i+1;
	}
	output << endl;

	output << left << setw(18) << "!z(cm)" << setw(14) << "time(yr)" << setw(14) << "T_n" << setw(14) << "T_i" 
		<< setw(14) << "T_e" << setw(14) << "v_n"<< setw(14) << "v_i" << setw(14) << "nh" << setw(14) << "aic" 
		<< setw(14) << "ae" << setw(14) << "aio" << setw(14) << "vgn(av)" << setw(14) << "vgn" 
        << setw(14) << "vgi(av)" << setw(14) << "vgi" << setw(14) << "p1" << setw(14) << "p2" << setw(14) << "p3" << endl;
	output.close();
	
	// create_file_cloud_parameters(output_path2); // is not saved, saving space on the disk,
	create_file_specimen_abund(output_path2, network);
	create_file_heating_rates(output_path2);
	create_file_energy_fluxes(output_path2);
	create_file_chem_hating(output_path2);
	create_file_dust_properties(output_path2, dust);
	create_file_reaction_rates(output_path2, network, save_disk_space = true);
	create_file_mol_data(output_path2, &user_data);
	create_file_h2_chemistry(output_path2);

#if (SAVE_RADIATIVE_FACTORS)
	user_data.create_file_radiative_transfer(output_path2);
#endif

	nb_saved = nb_not_saved = tot_nb_steps = 0;
	nb_saved_cloud_param = -1;
	is_post_shock = must_be_stopped = false;
	
	veln_arr.push_back( NV_Ith_S(y, nb_mhd+3) );
	veli_arr.push_back( NV_Ith_S(y, nb_mhd+4) );

	while (z < zfin && !must_be_stopped) 
	{
		i = 0;
		flag = CV_SUCCESS; 
        zout = z + dz; // updating zout

		while (i < 100 && flag == CV_SUCCESS && z < zout) {
			flag = CVode(cvode_mem, zout, y, &z, CV_ONE_STEP); // CV_NORMAL or CV_ONE_STEP     
     	    i++;
		}

        dz += z - zout;
		if (dz < 10.*DBL_EPSILON*z) // at this moment, dz may be very small or negative due to rounding error;
			dz = 10.*DBL_EPSILON*z;

		tot_nb_steps += i;
		ty += 2.*dz/(YEARS_TO_SECONDS*(prev_y[nb_mhd + 3] + NV_Ith_S(y, nb_mhd + 3)));
  
        if (verbosity) {
            cout << scientific;
            cout.precision(3);

            cout << left << "z (cm): " << setw(12) << z << "dz (cm): " << setw(12) << dz
                << "calc time (s): " << setw(8) << (int)(time(NULL) - timer) << "nb of steps: " << i << endl;
        }

		// Calculation of velocity gradients, some arbitrary parameter for velocity difference,
		// these parameters must be close to the instantaneous velocity gradients
		dz_arr.push_back(dz);
		i = (int) dz_arr.size();
        a = 0.;
		
        do {
			i--;
			a += dz_arr[i];
        } while (i > 0 && fabs(NV_Ith_S(y, nb_mhd + 3) - veln_arr[i]) < dv_to_v_lim *NV_Ith_S(y, nb_mhd + 3));
        vel_n_grad = (NV_Ith_S(y, nb_mhd+3) - veln_arr[i])/a;
        
        a = 0.;
        i = (int)dz_arr.size();
        do {
            i--;
            a += dz_arr[i];
        } while (i > 0 && fabs(NV_Ith_S(y, nb_mhd + 4) - veli_arr[i]) < dv_to_v_lim *NV_Ith_S(y, nb_mhd + 4));
        vel_i_grad = (NV_Ith_S(y, nb_mhd+4) - veli_arr[i])/a;
         
		veln_arr.push_back( NV_Ith_S(y, nb_mhd+3) );
		veli_arr.push_back( NV_Ith_S(y, nb_mhd+4) );

		conc_h_tot = user_data.calc_conc_h_tot(y);
		h2_form_const = user_data.get_h2_form_grains()/(conc_h_tot *NV_Ith_S(y, network->h_nb));

        // Saving data,
        // second condition - in order to look for the cause of crashes
        if (z - z_saved > 0.05 * magn_precursor_length || nb_not_saved > 10) {
            nb_not_saved = 0;
            nb_saved++;
            
            save_specimen_abund(output_path2, nb_of_species, y, conc_h_tot, z);
            save_heating_rates(output_path2, &user_data, z);
            save_energy_fluxes(output_path2, &user_data, y, z, z - z_saved); 
            save_dust_properties(output_path2, &user_data, y, conc_h_tot, z);
            save_mol_data(output_path2, &user_data, y, z);
            save_file_h2_chemistry(output_path2, &user_data, y, z);

            if (nb_saved % 2 == 0 || flag != CV_SUCCESS) {
                save_reaction_rates(output_path2, &user_data, z, NV_Ith_S(y, nb_mhd));
            }

            i = (int)(NV_Ith_S(y, nb_mhd) / 100.); // in this case saving data each 100 K
            if (i != nb_saved_cloud_param)
            {
                nb_saved_cloud_param = i;
                //save_cloud_parameters(&user_data, output_path2, ty, visual_extinct, cr_ioniz_rate, uv_field_strength, ir_field_strength, c_abund_pah, y);
                save_chem_heating_rates(output_path2, &user_data, ty);
#if (SAVE_RADIATIVE_FACTORS)
                user_data.save_radiative_transfer_data(output_path2, ty);
#endif
            }
            user_data.calc_ion_dens(y, ion_conc, ion_pah_conc, ion_dens, ion_pah_dens, ion_dust_dens);

            fname = output_path2 + "sim_phys_param.txt";
            output.open(fname.c_str(), ios::app);
            output << scientific;

            output.precision(10);
            output << left << setw(18) << z;

            output.precision(5);
            output << setw(14) << ty;
            for (i = nb_mhd; i < nb_mhd + NB_MHD_EQUATIONS; i++) {
                output << left << setw(14) << NV_Ith_S(y, i);
            }
            output << left << setw(14) << conc_h_tot
                << setw(14) << user_data.calc_ice_conc(y) / conc_h_tot << setw(14) << NV_Ith_S(y, network->e_nb) / conc_h_tot
                << setw(14) << ion_conc / conc_h_tot << setw(14) << vel_n_grad << setw(14) << user_data.get_velg_mhd_n()
                << setw(14) << vel_i_grad << setw(14) << user_data.get_velg_mhd_i() << setw(14) << h2_form_const
                << setw(14) << user_data.get_add_electron_sterm() << setw(14) << user_data.calc_total_grain_charge(y) << endl;
            output.close();

            z_saved = z;
        }
        else {
            nb_not_saved++;
        }
		// first we save, than we break;
		if (flag != CV_SUCCESS) {
			cout << "Error ocurred in solver" << flag << endl;
			break;
		}
        if (user_data.is_ion_fluid_supersonic()) {
            cout << "Ion fluid is supersonic" << endl;
            shock_state = SHOCK_STATE_ION_SUPERSONIC;
            break;
        }
        if (user_data.is_neut_fluid_subsonic()) {
            cout << "Neutral fluid is subsonic" << endl;
            shock_state = SHOCK_STATE_NEUTRAL_SUBSONIC;
            break;
        }
        if (NV_Ith_S(y, network->h_nb) > 198. * NV_Ith_S(y, network->h2_nb)) {
            cout << "Almost total (99%) H2 dissociation took place" << endl;
            shock_state = SHOCK_STATE_H2_DISSOCIATION;
            break;
        }

        a = 0.;
        for (i = nb_mhd; i < nb_of_equat; i++)
        {
            b = fabs(NV_Ith_S(y, i) - prev_y[i])/prev_y[i];
            if (b > a)
                a = b;
        }

        // updating dz:
        if (a < 0.001) // to verify
            dz *= 2;
        else if (a < 0.02)
            dz *= 1.15;
        else if (a > 0.1)
            dz /= 2.;

		// shock begins:
		if ( fabs(NV_Ith_S(y, nb_mhd + 3) - NV_Ith_S(y, nb_mhd + 4)) > dvel_shock_stop*NV_Ith_S(y, nb_mhd + 3) )
			is_post_shock = true;
		
		// postshock gas comes to the equilibrium state:
		if ( is_post_shock && fabs(NV_Ith_S(y, nb_mhd + 3) - NV_Ith_S(y, nb_mhd + 4)) < dvel_shock_stop * NV_Ith_S(y, nb_mhd + 3) )
			must_be_stopped = true;

        is_new_vg = false;
		is_new_chd = user_data.recalc_grain_charge_ranges(y, new_y);
		if (is_new_chd)
		{
			user_data.get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
			
			N_VDestroy_Serial(y);
			y = N_VNew_Serial(nb_of_equat);

            N_VDestroy_Serial(abs_tol);
            abs_tol = N_VNew_Serial(nb_of_equat);
            user_data.set_tolerances(abs_tol);

			delete [] prev_y;
			prev_y = new double [nb_of_equat];

			for (i = 0; i < nb_of_equat; i++) {
                // check for negative values of abundance and level populations of species: 
                if (i < nb_dct - nb_of_grain_charges && new_y[i] < 0.) { 
                    new_y[i] = 0.;
                }
				NV_Ith_S(y, i) = new_y[i];
			}
			
            is_new_vg = true;
			user_data.set_veln_grad(vel_n_grad);
			user_data.set_veli_grad(vel_i_grad);
            
            SUNLinSolFree(LS);
			SUNMatDestroy(A);
            
            A = SUNDenseMatrix(nb_of_equat, nb_of_equat);
			LS = SUNDenseLinearSolver(y, A);
            
            // number of equations has been changed, 
            CVodeFree(&cvode_mem);
			cvode_mem = CVodeCreate(CV_BDF);

			flag = CVodeInit(cvode_mem, f_mhd, z, y);
            flag = CVodeSVtolerances(cvode_mem, rel_tol, abs_tol);
            flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
			flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
			flag = CVodeSetMaxErrTestFails(cvode_mem, MAX_ERR_TEST_FAILS_SOLVER); // default value is 7;
            flag = CVodeSetMaxConvFails(cvode_mem, MAX_CONV_FAILS_SOLVER); // default value is 10;
			flag = CVodeSetUserData(cvode_mem, &user_data);
		}

		if (!is_new_chd) { 
            a = fabs(vel_n_grad / user_data.get_veln_grad());
            b = fabs(vel_i_grad / user_data.get_veli_grad());

			// this parameter is important for radiative transport modelling (fluctuations of the level populations)
            if (((fabs(vel_n_grad) > user_data.get_vel_grad_min()) && (a > 1.1 || a < 0.9)) ||
                ((fabs(vel_i_grad) > user_data.get_vel_grad_min()) && (b > 1.1 || b < 0.9)))
            {
                is_new_vg = true;
                user_data.set_veln_grad(vel_n_grad);
                user_data.set_veli_grad(vel_i_grad);
                
                // restart of the solver with new values of velocity gradients
                flag = CVodeReInit(cvode_mem, z, y);
            }
		}

        if (verbosity && is_new_vg) {
            cout << scientific;
            cout.precision(3);

            cout << "new velocity gradients are assigned (cm/s/cm)," << endl
                << "    neutrals: " << user_data.get_veln_grad()
                << "    ions: " << user_data.get_veli_grad() << endl;
        }
		
        for (i = 0; i < nb_of_equat; i++) {
			prev_y[i] = NV_Ith_S(y, i);
		}
	}
    if (verbosity) {
        cout << "Total nb of steps: " << tot_nb_steps << endl;
    }
	// Free memory:
    CVodeFree(&cvode_mem);
	SUNLinSolFree(LS);
	SUNMatDestroy(A);
	
	N_VDestroy_Serial(y);
    N_VDestroy_Serial(abs_tol);
	delete [] prev_y;

#ifdef __linux__
    cout.rdbuf(orig_cout);
    cerr.rdbuf(orig_cerr);
#endif
    return shock_state;
}

void calc_cr_dominated_region(const string &data_path, const string &output_path1, const string &output_path2, double c_abund_pah, 
	double evol_time, double cr_ir_factor, double incr_time)
{
#ifdef __linux__
	stringstream lin_out;	
	lin_out << output_path2;
	lin_out << "out";
	lin_out << "_screen";
	lin_out << ".txt";
	// lin_out.str("/dev/null");

	ofstream outerr(lin_out.str().c_str(), ios::app);
	streambuf *orig_cerr = cerr.rdbuf(outerr.rdbuf());
	
	ofstream out(lin_out.str().c_str(), ios::app);
	streambuf *orig_cout = cout.rdbuf(out.rdbuf());
#endif

	bool is_new_chd, linear_time_flow;
	int i, flag, nb, nb_of_species, nb_lev_h2o, nb_lev_h2, nb_lev_co, nb_vibr_h2o, nb_vibr_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3,
		nb_vibr_ch3oh, nb_lev_ch3oh, nb_lev_ci, nb_lev_cii, nb_lev_oi, nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd, verbosity;
	double h2_form_const, t, time_step(0.), ty, tfin, tout, rel_tol, tmult, visual_extinct, cr_ioniz_rate, cr_ioniz_rate0, 
		uv_field_strength, ir_field_strength, conc_h_tot, ion_conc, ion_pah_conc, ion_dens, ion_pah_dens, ion_dust_dens, op_h2_ratio;
	long int nb_steps;
	
	time_t timer;
	string fname;
	ofstream output;
	
	SUNMatrix A(NULL);
	SUNLinearSolver LS(NULL);
	N_Vector y(0), ydot(0), abs_tol(0);
	vector<double> new_y;

	verbosity = 1;
	cout << scientific;
	cout.precision(4);
		
	timer = time(NULL);
	cout << ctime(&timer) << "Study of the chemistry response on the increase of CR flux" << endl;

	// Spectroscopic parameters for H2, H2O and CO molecule - the same as for static cloud
	nb_lev_h2 = 100; 	
	nb_vibr_h2o = 0;
	nb_lev_h2o = 45;
	nb_vibr_co = 0;
	nb_lev_co = 30; 

	nb_vibr_ch3oh = 0;
#if (CALCULATE_POPUL_METHANOL)
	nb_lev_ch3oh = 100;
#else
	nb_lev_ch3oh = 1;
#endif

#if (CALCULATE_POPUL_NH3_OH)
    nb_lev_onh3 = 17; // ortho-NH3: He coll data - 22, H2 coll data - 17
    nb_lev_pnh3 = 34;
    nb_lev_oh = 20;
#else
    nb_lev_onh3 = 1; // ortho-NH3: He coll data - 22, H2 coll data - 17
    nb_lev_pnh3 = 1;
    nb_lev_oh = 1;
#endif
	// initialization of the data necessary for differential equation integration:
	chemistry_evolution_data user_data(data_path, output_path2, nb_lev_h2, nb_vibr_h2o, nb_lev_h2o, nb_vibr_co, nb_lev_co, 
        nb_lev_pnh3, nb_lev_onh3, nb_lev_oh, nb_vibr_ch3oh, nb_lev_ch3oh, c_abund_pah, verbosity);

	nb_of_species = user_data.get_nb_of_species();
	user_data.get_nb_of_levels(nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, 
		nb_lev_ci, nb_lev_oi, nb_lev_cii);
	
	// reading initial data for simulations from file, output_path1 - path to the file,
	// data are saved to vector y, parameters of grain charge distribution are assigned in user_data, 
	assign_cloud_data(output_path1, &user_data, new_y, evol_time, visual_extinct, cr_ioniz_rate, uv_field_strength, 
		ir_field_strength, conc_h_tot, verbosity);
	user_data.get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
	
	y = N_VNew_Serial(nb_of_equat);
	ydot = N_VNew_Serial(nb_of_equat);
    abs_tol = N_VNew_Serial(nb_of_equat);
	
	for (i = 0; i < nb_of_equat; i++) {
		NV_Ith_S(y, i) = new_y[i];
	}

	const dust_model *dust 
		= user_data.get_dust();
	
	const chem_network *network 
		= user_data.get_network();

	cr_ioniz_rate0 = cr_ioniz_rate;
    incr_time *= YEARS_TO_SECONDS;
    evol_time *= YEARS_TO_SECONDS;
	
    t = evol_time; //
	tfin = 1.01e+8 *YEARS_TO_SECONDS;
	tmult = pow(10, 0.03125); // 1/16 = 0.0625; 1/32 = 0.03125
	
    if (incr_time > 1.) // more than 1 sec 
    {
        linear_time_flow = true;     
        time_step = 0.01*incr_time;
    }
    else 
    {
        linear_time_flow = false;
        cr_ioniz_rate = cr_ioniz_rate0 * cr_ir_factor;      
        user_data.set_parameters(visual_extinct, cr_ioniz_rate, uv_field_strength, ir_field_strength);
        time_step = 0.01*YEARS_TO_SECONDS;
    }
    tout = t + time_step;

	rel_tol = REL_ERROR_SOLVER;
    user_data.set_tolerances(abs_tol);
	
	// Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration 
	void *cvode_mem = CVodeCreate(CV_BDF);

	// Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(t,y), 
	// the inital time t0, and the initial dependent variable vector y;
	flag = CVodeInit(cvode_mem, f_chem, t, y);

	// Call CVodeSVtolerances to specify the scalar tolerances:
	flag = CVodeSVtolerances(cvode_mem, rel_tol, abs_tol);

	// The maximal number of steps between simulation stops;
	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);

	// specifies the maximum number of error test failures permitted in attempting one step:
	flag = CVodeSetMaxErrTestFails(cvode_mem, MAX_ERR_TEST_FAILS_SOLVER); // default value is 7;
    flag = CVodeSetMaxConvFails(cvode_mem, MAX_CONV_FAILS_SOLVER); // default value is 10;

	// Create dense SUNMatrix for use in linear solves
	A = SUNDenseMatrix(nb_of_equat, nb_of_equat);
	
	// Create dense SUNLinearSolver object for use by CVode 
	LS = SUNDenseLinearSolver(y, A);

	// Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode 
	flag = CVDlsSetLinearSolver(cvode_mem, LS, A);

	// The function attaches the user data block to the solver;
	flag = CVodeSetUserData(cvode_mem, &user_data);
	
	create_file_cloud_parameters(output_path2);
	create_file_heating_rates(output_path2);
	create_file_chem_hating(output_path2);
	create_file_specimen_abund(output_path2, network);
	create_file_ice_comp(output_path2);
	create_file_dust_properties(output_path2, dust);
	create_file_reaction_rates(output_path2, network);
	create_file_nautilus_data(output_path2);
	create_file_mhd_vode(output_path2);

#if (SAVE_RADIATIVE_FACTORS)
	user_data.create_file_radiative_transfer(output_path2);
#endif

	fname = output_path2 + "sim_phys_param.txt";
	output.open(fname.c_str());
	output << "! p1 - total abundance of grain mantle molecules, [molecules/H]" << endl
		<< "! p2 - abundance of hydrocarbon molecules in ice mantles, CnHm, n >=2,  [molecules/H]" << endl
		<< "! p3 - parameter of H2 formation on grains (reaction *H+*H), rate = a*n_H_tot*n_H, a [cm3 s-1]" << endl
		<< "! p4 - electron concentration, [cm-3]" << endl
		<< "! p5 - ion concentration (without PAHs and small grains), [cm-3]" << endl
		<< "! p6 - total electric charge of grains, [cm-3]" << endl;

	output << "!";
	for (i = 0; i < 11; i++) {
		output << left << setw(14) << i + 1;
	}
	output << endl << left << setw(15) << "!time(yrs)" << setw(14) << "T_n" << setw(14) << "T_i" << setw(14) << "T_e" << setw(14) << "oph2"
		<< setw(14) << "p1" << setw(14) << "p2" << setw(14) << "p3" << setw(14) << "p4" << setw(14) << "p5" << setw(14) << "p6" << endl;
	output.close();
	
	nb = 0;
	timer = time(NULL);
	while (tout < tfin)
	{
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL); // CV_NORMAL or CV_ONE_STEP
		if (flag != CV_SUCCESS) {
			cout << "Error ocurred in solver: " << flag << endl;
			break;
		}
		flag = CVodeGetNumSteps(cvode_mem, &nb_steps);
		check_flag(&flag, "CVodeGetNumSteps", 1);
		
		ty = t/YEARS_TO_SECONDS;
		if (verbosity) {
			cout.precision(3);
			cout << left << "cloud age (yr): " << setw(12) << ty << "calc time (s): " 
				<< setw(8) << (int) (time(NULL) - timer) << "nb of steps: " << nb_steps << endl;
		}

		// the function is called in order to update user data:
		user_data.f(t, y, ydot); 
		
		save_specimen_abund(output_path2, nb_of_species, y, conc_h_tot, ty);
		save_ice_comp(output_path2, network, y, ty);
		save_heating_rates(output_path2, &user_data, ty);
		save_dust_properties(output_path2, &user_data, y, conc_h_tot, ty);
		save_nautilus_data(output_path2, ty, visual_extinct, conc_h_tot, NV_Ith_S(y, nb_mhd), user_data.get_av_dust_temp());
		
		user_data.calc_ion_dens(y, ion_conc, ion_pah_conc, ion_dens, ion_pah_dens, ion_dust_dens);
		h2_form_const = user_data.get_h2_form_grains()/(conc_h_tot *NV_Ith_S(y, network->h_nb));
        op_h2_ratio = NV_Ith_S(y, network->h2_nb) / user_data.calc_conc_ph2(y) - 1.;

		fname = output_path2 + "sim_phys_param.txt";
		output.open(fname.c_str(), ios::app);
		output << scientific;
		
		output.precision(7);
		output << left << setw(15) << ty;

		output.precision(5);		
		output << left << setw(14) << NV_Ith_S(y, nb_mhd) << setw(14) << NV_Ith_S(y, nb_mhd + 1) 
			<< setw(14) << NV_Ith_S(y, nb_mhd + 2) << setw(14) << op_h2_ratio << setw(14) << user_data.calc_ice_conc(y)/conc_h_tot
			<< setw(14) << user_data.calc_hydrocarbon_conc(y)/conc_h_tot << setw(14) << h2_form_const 
			<< setw(14) << NV_Ith_S(y, network->e_nb) << setw(14) << ion_conc
			<< setw(14) << user_data.calc_total_grain_charge(y) << endl;
		output.close();
		
		if (!linear_time_flow && ty < 0.99e+7)
		{
			save_cloud_parameters(&user_data, output_path2, ty, visual_extinct, cr_ioniz_rate, uv_field_strength, ir_field_strength, c_abund_pah, y);
			save_chem_heating_rates(output_path2, &user_data, ty);
			save_mhd_vode(data_path, output_path2, network, y, conc_h_tot, ty);

#if (SAVE_RADIATIVE_FACTORS)
			user_data.save_radiative_transfer_data(output_path2, ty);
#endif
		}
		if (nb%2 == 0) {
			save_reaction_rates(output_path2, &user_data, ty, NV_Ith_S(y, nb_mhd));
		}

        if (t > evol_time + incr_time)
            linear_time_flow = false;

        if (linear_time_flow) {
            cr_ioniz_rate = cr_ioniz_rate0 * (1. + (cr_ir_factor - 1)* (t - evol_time) / incr_time);
            user_data.set_parameters(visual_extinct, cr_ioniz_rate, uv_field_strength, ir_field_strength);
        }

		is_new_chd = user_data.recalc_grain_charge_ranges(y, new_y);
		if (is_new_chd) 
        {	
			user_data.get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
			
			N_VDestroy_Serial(y);
			N_VDestroy_Serial(ydot);
            N_VDestroy_Serial(abs_tol);

			y = N_VNew_Serial(nb_of_equat);
			ydot = N_VNew_Serial(nb_of_equat);
            abs_tol = N_VNew_Serial(nb_of_equat);

            user_data.set_tolerances(abs_tol);
			for (i = 0; i < nb_of_equat; i++) {
				NV_Ith_S(y, i) = new_y[i];
			}

            SUNLinSolFree(LS);
            SUNMatDestroy(A);

            A = SUNDenseMatrix(nb_of_equat, nb_of_equat);
            LS = SUNDenseLinearSolver(y, A);
		}
		
		if (is_new_chd || linear_time_flow) 
        {
			CVodeFree(&cvode_mem);
			cvode_mem = CVodeCreate(CV_BDF);
		
			flag = CVodeInit(cvode_mem, f_chem, t, y);
			flag = CVodeSVtolerances(cvode_mem, rel_tol, abs_tol);
			flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
			flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
			flag = CVodeSetMaxErrTestFails(cvode_mem, MAX_ERR_TEST_FAILS_SOLVER);
            flag = CVodeSetMaxConvFails(cvode_mem, MAX_CONV_FAILS_SOLVER); // default value is 10;
			flag = CVodeSetUserData(cvode_mem, &user_data);
		}

		nb++;
        if (!linear_time_flow)
            time_step *= tmult;

		tout += time_step;
	}
	
	// Free memory;
	SUNLinSolFree(LS);
	SUNMatDestroy(A);
			
	CVodeFree(&cvode_mem);
	N_VDestroy_Serial(y);
	N_VDestroy_Serial(ydot);
    N_VDestroy_Serial(abs_tol);
}

void create_file_cloud_parameters(const string &output_path)
{
	string fname;
	ofstream output;

	fname = output_path + "sim_cloud_data.txt";
	output.open(fname.c_str());
	output.close();
}

void save_cloud_parameters(const evolution_data *user_data, const string &output_path, double ty, double visual_extinct, 
	double cr_ioniz_rate, double uv_field_strength, double ir_field_strength, double c_abund_pah, const N_Vector &y)
{
	int i, k, j, l, nb, nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd, nb_of_dust_comp, nb_of_species, nb_lev_h2, nb_lev_h2o, 
		nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, nb_lev_ci, nb_lev_cii, nb_lev_oi;
	double conc_h_tot, conc_g;

	string fname;
	ofstream output;
	
	const dust_model *dust 
		= user_data->get_dust();

	const chem_network *network 
		= user_data->get_network();
	
	nb_of_species = user_data->get_nb_of_species();
	nb_of_dust_comp = dust->nb_of_comp;

	user_data->get_nb_of_levels(nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, 
		nb_lev_ci, nb_lev_oi, nb_lev_cii);
	user_data->get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
	conc_h_tot = user_data->calc_conc_h_tot(y);

	fname = output_path + "sim_cloud_data.txt";
	output.open(fname.c_str(), ios::app);
	
    output << scientific;
	output.precision(9); 
			
	output << "# Evolution age (years):" << endl;
	output << ty << endl; 

    output.precision(5);
	output << "# Visual extinction, CR ionization rate, UV and IR field strength, H nuclei concentration, carbon abundance in PAH:" << endl;
	output << left << setw(14) << visual_extinct << setw(14) << cr_ioniz_rate << setw(14) << uv_field_strength << setw(14) << ir_field_strength
		<< setw(14) << conc_h_tot << setw(14) << c_abund_pah << endl;

	output << "# Neutral gas temperature, ion and electron temperatures (in K):" << endl;
	output << left << setw(14) << NV_Ith_S(y, nb_mhd) << setw(14) << NV_Ith_S(y, nb_mhd + 1) << 
		setw(14) << NV_Ith_S(y, nb_mhd + 2) << endl;

	output << left << "# Dust temperatures and charges:" << endl << nb_of_dust_comp << endl
		<< setw(3) << "Nb" << setw(10) << "name" << setw(14) << "temp(K)" << setw(14) << "conc(cm-3)" << setw(14) << "av.charge"
		<< setw(11) << "distr(1/0)" << setw(6) << "zmin" << setw(6) << "zmax" << "distribution on charge" << endl;
				
	for (i = 0; i < nb_of_dust_comp; i++) 
	{
		user_data->get_dust_component_nbs(i, j, k);
		conc_g = user_data->calc_grain_conc(y, i);

		output << left << setw(3) << i << setw(10) << dust->components[i]->name << setw(14) << NV_Ith_S(y, nb_dct + i) 
			<< setw(14) << conc_g << setw(14) << user_data->calc_average_grain_charge(y, i) 
			<< setw(11) << ((k - j > 2) ? "1" : "0");

		if (k - j > 2) {
			output << left << setw(6) << user_data->get_dust_zmin(i) << setw(6) << user_data->get_dust_zmax(i);
			for (l = 0; l < k - j; l++) 
			{
				if (fabs(NV_Ith_S(y, j + l)) > MINIMAL_ABUNDANCE*conc_g)
					output << left << setw(14) << NV_Ith_S(y, j + l);
				else output << left << setw(14) << "0.";
			}
		}
		output << endl;
	}
			
	output << "# Concentration of chemical species:" << endl << nb_of_species << endl;
	for (i = 0; i < nb_of_species; i++) 
	{
		output << left << setw(5) << i + 1 << setw(15) << network->species[i].name << " ";
		if (NV_Ith_S(y, i) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i) << endl;
		else output << "0." << endl;
	}

	output << "# The level population densities of H2 molecule:" << endl << nb_lev_h2 << endl;
	for (i = 0; i < nb_lev_h2; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb_of_species) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb_of_species) << endl;
		else output << "0." << endl;
	}
	nb = nb_of_species + nb_lev_h2;

	output << "# The level population densities of para-H2O molecule:" << endl << nb_lev_h2o << endl;
	for (i = 0; i < nb_lev_h2o; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_h2o;
	
	output << "# The level population densities of ortho-H2O molecule:" << endl << nb_lev_h2o << endl;
	for (i = 0; i < nb_lev_h2o; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_h2o;
	
	output << "# The level population densities of CO molecule:" << endl << nb_lev_co << endl;
	for (i = 0; i < nb_lev_co; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_co;
	
	output << "# The level population densities of OH molecule:" << endl << nb_lev_oh << endl;
	for (i = 0; i < nb_lev_oh; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_oh;
	
	output << "# The level population densities of para-NH3 molecule:" << endl << nb_lev_pnh3 << endl;
	for (i = 0; i < nb_lev_pnh3; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_pnh3;

	output << "# The level population densities of ortho-NH3 molecule:" << endl << nb_lev_onh3 << endl;
	for (i = 0; i < nb_lev_onh3; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_onh3;

	output << "# The level population densities of CH3OH-A molecule:" << endl << nb_lev_ch3oh << endl;
	for (i = 0; i < nb_lev_ch3oh; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_ch3oh;

	output << "# The level population densities of CH3OH-E molecule:" << endl << nb_lev_ch3oh << endl;
	for (i = 0; i < nb_lev_ch3oh; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, i + nb) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, i + nb) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_ch3oh;

	output << "# The level population densities of CI:" << endl << nb_lev_ci << endl;
	for (i = 0; i < nb_lev_ci; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, nb + i) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, nb + i) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_ci;
	
	output << "# The level population densities of OI:" << endl << nb_lev_oi << endl;
	for (i = 0; i < nb_lev_oi; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, nb + i) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, nb + i) << endl;
		else output << "0." << endl;
	}
	nb += nb_lev_oi;
	
	output << "# The level population densities of CII:" << endl << nb_lev_cii << endl;
	for (i = 0; i < nb_lev_cii; i++) 
	{
		output << left << setw(4) << i + 1;
		if (NV_Ith_S(y, nb + i) > MINIMAL_ABUNDANCE *conc_h_tot)
			output << NV_Ith_S(y, nb + i) << endl;
		else output << "0." << endl;
	}
	output.close();
}

void create_file_chem_hating(const string &output_path)
{
	string fname = output_path + "sim_chem_heating.txt";
	ofstream output;

	output.open(fname.c_str());
	output << "# Heating and cooling rates of neutral gas component due to chemical reactions, [erg cm-3 s-1]" << endl;
	output.close();
}

void save_chem_heating_rates(const string &output_path, const evolution_data *user_data, double ty)
{
	int i, nb;
	double chem_heat, chem_cool;
	double *chem_heating_rates_n(0);
	
	string fname = output_path + "sim_chem_heating.txt";
	ofstream output;

	const chem_network *network 
		= user_data->get_network();

	chem_heat = chem_cool = 0.;
	user_data->get_chem_heating_rates(chem_heating_rates_n, nb);

	for (i = 0; i < nb; i++) 
	{
		if (chem_heating_rates_n[i] > 0.) 
			chem_heat += chem_heating_rates_n[i];
		else chem_cool += chem_heating_rates_n[i];
	}

	output.open(fname.c_str(), ios::app);
	output << scientific;
	output.precision(3);
			
	output << "# Evolution age (years):" << endl;
	output << ty << endl;

	output << chem_heat << endl;
	for (i = 0; i < nb; i++) {
		if (chem_heating_rates_n[i] >= 0.01*chem_heat)
			output << left << setw(12) << chem_heating_rates_n[i] << network->reaction_array[i].name << endl;
	}

	output << chem_cool << endl;
	for (i = 0; i < nb; i++) {
		if (chem_heating_rates_n[i] <= 0.01*chem_cool)
			output << left << setw(12) << chem_heating_rates_n[i] << network->reaction_array[i].name << endl;
	}
	output.close();
}

void create_file_specimen_abund(const string &output_path, const chem_network *network)
{
	int i, l;
	string fname;
	ofstream output;

	fname = output_path + "sim_specimen_abund.txt";
	output.open(fname.c_str());

	output << "!    ";
	for (i = 0; i < network->nb_of_species + 1; i++) {
		output << left << setw(11) << i + 1;
	}
	output << endl;
	output << left << setw(16) << "!t(yr)/z(cm)"; 
    
    // the length of the name of chemical specimen may be large
	for (i = 0; i < network->nb_of_species; i++) {
        l = (int) network->species[i].name.length();
        l = (l >= 11) ? l + 1 : 11;
		output << left << setw(l) << network->species[i].name;
	}
	output << endl;
	output.close();
}

void save_specimen_abund(const string &output_path, int nb_of_species, const N_Vector &y, double conc_h_tot, double var)
{
	int i;
	string fname;
	ofstream output;

	fname = output_path + "sim_specimen_abund.txt";
	output.open(fname.c_str(), ios::app);
	output << scientific;

	output.precision(8);
	output << left << setw(16) << var;
	output.precision(3);

	for (i = 0; i < nb_of_species; i++) { // abundance values are supposed to have a form (-)x.xxxe-xx 
		output << left << setw(11) << NV_Ith_S(y, i)/conc_h_tot;
	}
	output << endl;
	output.close();
}

void create_file_heating_rates(const string &output_path)
{
	int i;
	string fname;
	ofstream output;

	fname = output_path + "sim_heating_cooling.txt";
	output.open(fname.c_str());

	output << "! na - neutral gas cooling by atomic emission, cooling < 0., [erg cm-3 s-1]" << endl
		<< "! nh2 - neutral gas cooling by H2 emission" << endl
		<< "! nh2o - neutral gas cooling by H2O emission" << endl
		<< "! nco - neutral gas cooling by CO emission" << endl
		<< "! noh - neutral gas cooling by OH emission" << endl
		<< "! nnh3 - neutral gas cooling by p- and o-NH3 emission" << endl
		<< "! nch3oh - neutral gas cooling by CH3OH (A-,E-) emission" << endl
		<< "! ni - neutral gas heating by neutral-ion collisions" << endl
		<< "! ne - neutral gas heating by neutral-electron collisions" << endl
		<< "! nd - neutral gas heating by neutral-dust collisions" << endl
		<< "! nch - neutral gas heating by chemical reactions" << endl
        << "! nh2-h - neutral gas cooling by H2-H dissociation (is included in previous parameter nch)" << endl
		<< "! ph - heating by photoeffect on dust grains" << endl
		<< "! cr - heating by cosmic rays" << endl
		<< "! rlh2 - radiative energy loss by H2 molecule (total and only within state v = 0)" << endl
        << "! ih2 - ion fluid cooling via excitation of H2 via ion-H2 collisions" << endl
		<< "! in - ion gas component heating by neutral-ion collisions" << endl
		<< "! ie - ion gas component heating by ion-electron collisions" << endl
		<< "! ich - ion gas component heating by chemical reactions" << endl
		<< "! ea - electron component cooling by atomic emission" << endl
		<< "! eh2 - electron component cooling by H2 emission" << endl
		<< "! eh2o - electron component cooling by H2O emission" << endl
		<< "! en - electron component heating by neutral-electron collisions" << endl
		<< "! ei - electron component heating by ion-electron collisions" << endl
		<< "! ech - electron component heating by chemical reactions" << endl;

	output << left << "!";
	for (i = 0; i < 27; i++) {
		output << left << setw(12) << i + 1;
	}
	output << endl;

	output << left << setw(13) << "!t(yr)/z(cm)" << setw(12) << "na" << setw(12) << "nh2" << setw(12) << "nh2o" 
		<< setw(12) << "nco" << setw(12) << "noh" << setw(12) << "nnh3" << setw(12) << "nch3oh" << setw(12) << "ni" 
		<< setw(12) << "ne" << setw(12) << "nd" << setw(12) << "nch" << setw(12) << "nh2-h" << setw(12) << "ph" 
        << setw(12) << "cr" << setw(12) << "rlh2" << setw(12) << "rlh2v0" << setw(12) << "ih2" << setw(12) << "in" 
        << setw(12) << "ie" << setw(12) << "ich" << setw(12) << "ea" << setw(12) << "eh2"<< setw(12) << "eh2o" 
        << setw(12) << "en" << setw(12) << "ei" << setw(12) << "ech" << endl;
	output.close();
}

void save_heating_rates(const string &output_path, const evolution_data *user_data, double var)
{
	double neut_heat_atoms, neut_heat_h2, neut_heat_h2o, neut_heat_co, neut_heat_oh, neut_heat_nh3, neut_heat_ch3oh, neut_heat_dust_coll, 
		neut_heat_chem, pheff_gas_heat, neut_cr_heat, neut_heat_scatt_ions, neut_heat_scatt_el, el_heat_atoms, el_heat_h2, 
		el_heat_h2o, el_heat_scatt_neut, el_heat_scatt_ions, el_heat_chem, ion_heat_h2, ion_heat_scatt_n, ion_heat_scatt_el, 
        ion_heat_chem, rad_energy_loss_h2, rad_energy_loss_h2v0, h2_h_diss_cooling;

	string fname;
	ofstream output;
	
	user_data->get_neutral_heating(neut_heat_atoms, neut_heat_h2, neut_heat_h2o, neut_heat_co, neut_heat_oh, neut_heat_nh3, neut_heat_ch3oh, 
		neut_heat_dust_coll, neut_heat_chem, pheff_gas_heat, neut_cr_heat, neut_heat_scatt_ions, neut_heat_scatt_el, rad_energy_loss_h2, 
        rad_energy_loss_h2v0, h2_h_diss_cooling);
	user_data->get_electron_heating(el_heat_atoms, el_heat_h2, el_heat_h2o, el_heat_scatt_neut, el_heat_scatt_ions, el_heat_chem);
	user_data->get_ion_heating(ion_heat_h2, ion_heat_scatt_n, ion_heat_scatt_el, ion_heat_chem);

	fname = output_path + "sim_heating_cooling.txt";
	output.open(fname.c_str(), ios::app);
	output << scientific;

	output.precision(5);
	output << left << setw(13) << var;
	output.precision(3);
	
	output << left << setw(12) << neut_heat_atoms 
		<< setw(12) << neut_heat_h2 
		<< setw(12) << neut_heat_h2o 
		<< setw(12) << neut_heat_co 
		<< setw(12) << neut_heat_oh 
		<< setw(12) << neut_heat_nh3 
		<< setw(12) << neut_heat_ch3oh
		<< setw(12) << neut_heat_scatt_ions 
		<< setw(12) << neut_heat_scatt_el 
		<< setw(12) << neut_heat_dust_coll 
		<< setw(12) << neut_heat_chem 
        << setw(12) << h2_h_diss_cooling
		<< setw(12) << pheff_gas_heat 
		<< setw(12) << neut_cr_heat 
		<< setw(12) << rad_energy_loss_h2 
        << setw(12) << rad_energy_loss_h2v0
        << setw(12) << ion_heat_h2
		<< setw(12) << ion_heat_scatt_n 
		<< setw(12) << ion_heat_scatt_el 
		<< setw(12) << ion_heat_chem 
		<< setw(12) << el_heat_atoms 
		<< setw(12) << el_heat_h2 
		<< setw(12) << el_heat_h2o 
		<< setw(12) << el_heat_scatt_neut 
		<< setw(12) << el_heat_scatt_ions 
		<< setw(12) << el_heat_chem << endl;
	output.close();
}

void create_file_energy_fluxes(const string &output_path)
{
	int i;
	string fname;
	ofstream output;

	fname = output_path + "sim_energy_fluxes.txt";
	output.open(fname.c_str());

	output << "! kin - kinetic energy flux (neutral and ions, dust grains), [erg cm-2 s-1]" << endl
        << "! thm - thermal energy flux (neutrals, ions, electrons)" << endl
        << "! mgn - magnetic energy flux" << endl
        << "! ncool - total integrated gas cooling by line emission (energy transfered via collisions to specimen excitation), cooling < 0., [erg cm-2 s-1]" << endl
        << "! na  - integrated neutral gas cooling by atomic emission" << endl
		<< "! nh2 - neutral gas cooling by H2 emission" << endl
		<< "! nh2o - neutral gas cooling by H2O emission" << endl
		<< "! nco - neutral gas cooling by CO emission" << endl
		<< "! noh - neutral gas cooling by OH emission" << endl
		<< "! nnh3 - neutral gas cooling by p- and o-NH3 emission" << endl
		<< "! nch3oh - neutral gas cooling by CH3OH (A-,E-) emission" << endl
        << "! kin_n - kinetic energy flux of neutrals" << endl
        << "! kin_i - kinetic energy flux of ions (no PAH and small grains)" << endl
        << "! kin_d - kinetic energy flux of dust (including PAH and small grains)" << endl
        << "! thm_n - thermal energy flux of neutrals" << endl
        << "! thm_i - thermal energy flux of ions, electrons, small negative grains (PAH)" << endl
        << "! totfl - total flux, [erg cm-2 s-1]" << endl;

	output << left << "!";
	for (i = 0; i < 18; i++) {
		output << left << setw(12) << i + 1;
	}
	output << endl;

	output << left << setw(13) << "!z(cm)" << setw(12) << "kin" << setw(12) << "thm" << setw(12) << "mgn"
        << setw(12) << "ncool" << setw(12) << "na" << setw(12) << "nh2" << setw(12) << "nh2o" 
        << setw(12) << "nco" << setw(12) << "noh" << setw(12) << "nnh3" << setw(12) << "nch3oh" 
        << setw(12) << "kin_n" << setw(12) << "kin_i" << setw(12) << "kin_d" << setw(12) << "thm_n" << setw(12) << "thm_i"
        << setw(12) << "totfl" << endl;
	output.close();
}

void save_energy_fluxes(const string &output_path, const evolution_data *user_data, const N_Vector &y, double var, double dvar)
{
    int nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd;
	static double neut_heat_atoms(0.), neut_heat_h2(0.), neut_heat_h2o(0.), neut_heat_co(0.), neut_heat_oh(0.), neut_heat_nh3(0.), 
		neut_heat_ch3oh(0.), int_neut_heat_atoms(0.), int_neut_heat_h2(0.), int_neut_heat_h2o(0.), int_neut_heat_co(0.), int_neut_heat_oh(0.), 
		int_neut_heat_nh3(0.), int_neut_heat_ch3oh(0.); 
	double a, x1, x2, x3, x4, x5, x6, x7, neut_heat_dust_coll, neut_heat_chem, pheff_gas_heat, neut_cr_heat, neut_heat_scatt_ions, 
		neut_heat_scatt_el, h2_h_diss_cooling, total_mol_cooling, neut_conc, neut_mass_dens, ion_conc, ion_pah_conc,
        ion_dens, ion_pah_dens, ion_dust_dens, v_n, v_i, kin_energy_flux, kin_energy_flux_ions, kin_energy_flux_neut, 
        kin_energy_flux_dust, thermal_energy_flux, thermal_energy_flux_neut, thermal_energy_flux_ions, magnetic_energy_flux;

	string fname;
	ofstream output;
    
	x1 = neut_heat_atoms;
	x2 = neut_heat_h2;
	x3 = neut_heat_h2o;
	x4 = neut_heat_co;
	x5 = neut_heat_oh;
	x6 = neut_heat_nh3;
	x7 = neut_heat_ch3oh;

	user_data->get_neutral_heating(neut_heat_atoms, neut_heat_h2, neut_heat_h2o, neut_heat_co, neut_heat_oh, neut_heat_nh3, neut_heat_ch3oh, 
		neut_heat_dust_coll, neut_heat_chem, pheff_gas_heat, neut_cr_heat, neut_heat_scatt_ions, neut_heat_scatt_el, a, a, 
        h2_h_diss_cooling);
	
    if (fabs(var - dvar) <= numeric_limits<double>::epsilon() * fabs(var + dvar))
    {
        neut_heat_atoms = neut_heat_h2 = neut_heat_h2o = neut_heat_co = neut_heat_oh = neut_heat_nh3 = neut_heat_ch3oh = 0.;
        int_neut_heat_atoms = int_neut_heat_h2 = int_neut_heat_h2o = int_neut_heat_co = int_neut_heat_oh = int_neut_heat_nh3 =
            int_neut_heat_ch3oh = 0.;
    }

	int_neut_heat_atoms += (x1 + neut_heat_atoms) *0.5*dvar;
	int_neut_heat_h2 += (x2 + neut_heat_h2) *0.5*dvar;
	int_neut_heat_h2o += (x3 + neut_heat_h2o) *0.5*dvar;
	int_neut_heat_co += (x4 + neut_heat_co) *0.5*dvar;
	int_neut_heat_oh += (x5 + neut_heat_oh) *0.5*dvar;
	int_neut_heat_nh3 += (x6 + neut_heat_nh3) *0.5*dvar;
	int_neut_heat_ch3oh += (x7 + neut_heat_ch3oh) *0.5*dvar;

    total_mol_cooling = int_neut_heat_atoms + int_neut_heat_h2 + int_neut_heat_h2o + int_neut_heat_co +
        int_neut_heat_oh + int_neut_heat_nh3 + int_neut_heat_ch3oh;

    user_data->get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
    user_data->calc_neutral_dens(y, neut_conc, neut_mass_dens);
    user_data->calc_ion_dens(y, ion_conc, ion_pah_conc, ion_dens, ion_pah_dens, ion_dust_dens);

	fname = output_path + "sim_energy_fluxes.txt";
	output.open(fname.c_str(), ios::app);
	output << scientific;

	output.precision(5);
	output << left << setw(13) << var;
	output.precision(3);
	
    v_n = NV_Ith_S(y, nb_mhd + 3);
    v_i = NV_Ith_S(y, nb_mhd + 4); 
    
    kin_energy_flux_neut = 0.5 * v_n * v_n * v_n * neut_mass_dens;
    kin_energy_flux_ions = 0.5 * v_i * v_i * v_i * ion_dens;
    kin_energy_flux_dust = user_data->calc_dust_kinetic_energy_flux(y);
    kin_energy_flux = kin_energy_flux_neut + kin_energy_flux_ions + kin_energy_flux_dust;

    thermal_energy_flux_neut = 2.5 * BOLTZMANN_CONSTANT * v_n * neut_conc * NV_Ith_S(y, nb_mhd);
    thermal_energy_flux_ions = 2.5 * BOLTZMANN_CONSTANT * v_i *
        ((ion_conc + ion_pah_conc) * NV_Ith_S(y, nb_mhd + 1) + NV_Ith_S(y, user_data->get_network()->e_nb) * NV_Ith_S(y, nb_mhd + 2));
    thermal_energy_flux = thermal_energy_flux_neut + thermal_energy_flux_ions;
    
    magnetic_energy_flux = user_data->get_magnetic_field() * user_data->get_magnetic_field() * v_i / (4.*M_PI);
    
	output << left << setw(12) << kin_energy_flux
        << setw(12) << thermal_energy_flux 
        << setw(12) << magnetic_energy_flux
        << setw(12) << total_mol_cooling
        << setw(12) << int_neut_heat_atoms 
		<< setw(12) << int_neut_heat_h2 
		<< setw(12) << int_neut_heat_h2o 
		<< setw(12) << int_neut_heat_co 
		<< setw(12) << int_neut_heat_oh 
		<< setw(12) << int_neut_heat_nh3 
		<< setw(12) << int_neut_heat_ch3oh 
        << setw(12) << kin_energy_flux_neut
        << setw(12) << kin_energy_flux_ions
        << setw(12) << kin_energy_flux_dust
        << setw(12) << thermal_energy_flux_neut
        << setw(12) << thermal_energy_flux_ions
        << setw(12) << fabs(total_mol_cooling) + kin_energy_flux + thermal_energy_flux + magnetic_energy_flux << endl;
	output.close();
}

void create_file_dust_properties(const string &output_path, const dust_model *dust)
{
	int i;
	string fname;
	stringstream ss;
	ofstream output;

	fname = output_path + "sim_dust_data.txt";
	output.open(fname.c_str());
	
	output << "! the number in the parameter name is the grain component number" << endl
		<< "! T - temperature of dust component, [K]" << endl
		<< "! ab - grain abundance, [grains/H]" << endl
		<< "! Z - average grain charge" << endl
		<< "! abZ - average grain charge multiplied by grain abundance, [charge/H]" << endl
		<< "! ic - parameter characterizes grain coupling to ions, wt2/(1 + wt2)" << endl
		<< "! gf - flux of dust grains, [cm-2 s-1]" << endl
		<< "! hir - grain heating by IS radiation, [erg s-1]" << endl
		<< "! hch - grain heating by chemical reactions on the surface, [erg s-1]" << endl
		<< "! hc - grain heating by collisions with neutral gas, [erg s-1]" << endl
		<< "! hh2 - grain heating by absorption of H2 emission, [erg s-1]" << endl
		<< "! hm - grain heating by absorption of molecule emission (except H2), [erg s-1]" << endl
		<< "! ae - rate of electron attachment, [cm-3 s-1]" << endl
		<< "! ai - rate of ion attachment, [cm-3 s-1]" << endl
		<< "! puv - rate of photoelectron emission by IS UV radiation, [cm-3 s-1]" << endl
		<< "! pvis - rate of photoelectron emission by IS VIS radiation, [cm-3 s-1]" << endl
		<< "! pcr - rate of photoelectron emission by CR induced radiation, [cm-3 s-1]" << endl
		<< "! at the end, T - average dust temperature, grains with radius > MIN_ADSORPTION_RADIUS are taken into acount" << endl;

	output << left << "!";
	for (i = 0; i < 16*dust->nb_of_comp + 3; i++) {
		output << left << setw(12) << i + 1;
	}
	output << endl;

	output << left << setw(13) << "!t(yr)/z(cm)";
	for (i = 0; i < dust->nb_of_comp; i++) 
	{
		ss.clear();
		ss << i;
		output << left << setw(12) << "T" + ss.str() << setw(12) << "ab" + ss.str() << setw(12) << "Z" + ss.str() << setw(12) << "abZ" + ss.str() 
			<< setw(12) << "ic" + ss.str() << setw(12) << "gf" + ss.str() << setw(12) << "hir" + ss.str() << setw(12) << "hch" + ss.str() 
			<< setw(12) << "hc" + ss.str() << setw(12) << "hh2" + ss.str() << setw(12) << "hm" + ss.str() 
			<< setw(12) << "ae" + ss.str() << setw(12) << "ai" + ss.str() << setw(12) << "puv" + ss.str() << setw(12) << "pvis" + ss.str()
			<< setw(12) << "pcr" + ss.str();
	}
	output << left << setw(12) << "T" << setw(12) << "abZ" << endl;
	output.close();
}

void save_dust_properties(const string &output_path, const evolution_data *user_data, const N_Vector &y, double conc_h_tot, double var)
{
	int i, nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd;
	double abund, charge, h_isrf, h_chem, h_coll, h_h2_mol, h_mol, el_att, ion_neutr, phem_uv, phem_vis, phem_cr, vel, conc;
	
	string fname;
	ofstream output;

	const dust_model *dust 
		= user_data->get_dust();
	user_data->get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);

	fname = output_path + "sim_dust_data.txt";
	output.open(fname.c_str(), ios::app);
	output << scientific;
	
	output.precision(5);
	output << left << setw(13) << var;
	output.precision(3);

	for (i = 0; i < dust->nb_of_comp; i++) 
	{
		charge = user_data->calc_average_grain_charge(y, i);
		abund = user_data->calc_grain_conc(y, i)/conc_h_tot;
		vel = user_data->calc_av_grain_velocity(y, i);
		conc = user_data->calc_grain_conc(y, i);

		user_data->get_dust_heating_rates(i, h_isrf, h_chem, h_coll, h_h2_mol, h_mol);
		user_data->get_dust_charging_rates(i, el_att, ion_neutr, phem_uv, phem_vis, phem_cr);

		output << left << setw(12) << NV_Ith_S(y, nb_dct + i) << setw(12) << abund << setw(12) << charge
			<< setw(12) << abund *charge << setw(12) << user_data->calc_dust_ion_coupling(y, i) << setw(12) << vel*conc
			<< setw(12) << h_isrf << setw(12) << h_chem << setw(12) << h_coll << setw(12) << h_h2_mol << setw(12) << h_mol
			<< setw(12) << el_att << setw(12) << ion_neutr << setw(12) << phem_uv << setw(12) << phem_vis << setw(12) << phem_cr;
	}
	output << left << setw(12) << user_data->get_av_dust_temp() << setw(12) << user_data->calc_total_grain_charge(y)/conc_h_tot << endl;
	output.close();
}

void create_file_reaction_rates(const string & output_path, const chem_network *network, bool save_disk_space)
{
	int i, j, n;
	string fname, path;
	ofstream output;

    j = 0;
    n = (int) output_path.size();
    for (i = n - 1; i >= 0, j < 2; i--) {
        if (output_path[i] == '/')
            j++;
    }

    if (save_disk_space && j == 2)
        path = output_path.substr(0, i + 2);
    else path = output_path;

	fname = path + "sim_species.txt";
	output.open(fname.c_str());
	
	output << left << setw(8) << network->nb_of_species << setw(8) << NB_OF_CHEM_ELEMENTS << endl;
	for (i = 0; i < network->nb_of_species; i++) 
	{
		output << left << setw(5) << i+1 << setw(15) << network->species[i].name;
		for (j = 0; j < NB_OF_CHEM_ELEMENTS; j++) {
			output << left << setw(5) << network->species[i].formula[j];
		}
		output << endl;
	}
	output.close();

	fname = path + "sim_reactions.txt";
	output.open(fname.c_str());
	
	output << network->nb_of_reactions << endl;
	for (i = 0; i < network->nb_of_reactions; i++) 
	{
        const chem_reaction & reaction = network->reaction_array[i];
		
        output << left << setw(8) << i+1 << setw(5) << reaction.reactant[0];
		if (reaction.nb_of_reactants > 1)
			output << left << setw(5) << reaction.reactant[1];
		else output << left << setw(5) << "-1";

		// no more than four products can be in the reaction;
		for (j = 0; j < reaction.nb_of_products; j++) {
			output << left << setw(5) << reaction.product[j];
		}
		for (j = reaction.nb_of_products; j <= 3; j++) {
			output << left << setw(5) << "-1";
		}
		output << reaction.name << endl;
	}
	output.close();
	
    // binary format of the file
	fname = output_path + "sim_reaction_rates.bin";
	output.open(fname.c_str(), ios::binary);
    output.close();
}

void save_reaction_rates(const string &output_path, const evolution_data *user_data, double var, double temp_n)
{
    int i, k, l;
    int16_t abc; 
    int8_t de;

	double a;
	string fname;
	ofstream output;

	fname = output_path + "sim_reaction_rates.bin";
	output.open(fname.c_str(), ios::binary | ios::app);
	
    reformat_floating_value(var, k, l);
    abc = (int16_t) k;
    de = (int8_t)l;

    output.write((char*) & abc, sizeof(abc));
    output.write((char*) & de, sizeof(de));

    reformat_floating_value(temp_n, k, l);
    abc = (int16_t)k;
    de = (int8_t)l;

    output.write((char*)& abc, sizeof(abc));
    output.write((char*)& de, sizeof(de));

    for (i = 0; i < user_data->get_reaction_nb(); i++) {
        a = user_data->get_reaction_rate(i);
        reformat_floating_value(a, k, l);
        abc = (int16_t) k;
        de = (int8_t) l;

        output.write((char*)& abc, sizeof(abc));
        output.write((char*)& de, sizeof(de));
    }
    output.close();
}

// abc de = abc * 10 ^ (de), 
// abc must be fit to int16_t, de must be fit to int8_t, -128 =< de =< 127
void reformat_floating_value(double x, int & abc, int & de)
{
    const int precision = 3; // =< 3
    if (x > 1.e-99)
    {
        if (x < 1.)
            de = ((int)log10(x)) - 1 - precision;
        else de = (int)log10(x) - precision;

        abc = (int) (x *pow(10., -de));
    }
    else {
        abc = de = 0;
    }
}

void create_file_ice_comp(const string & output_path)
{
	string fname;
	ofstream output;

	fname = output_path + "sim_ice_comp.txt";
	output.open(fname.c_str());
	
	output << "! Specimen abundances are normalized on H2O ice abundance, H2O ice concentration is given" << endl;
	output << left << setw(13) << "!t(yr)/z(cm)" << setw(12) << "s-CO" << setw(12) << "s-CO2" << setw(12) << "s-H2CO"
		<< setw(12) << "s-CH3OH" << setw(12) << "s-CH4" << setw(12) << "s-NH3" << setw(12) << "s-N2" << setw(12) << "s-H2O_conc" 
		<< setw(30) << "! main nitrogen species" << endl;
	output.close();
}

void save_ice_comp(const string &output_path, const chem_network *network, const N_Vector &y, double var)
{
	int i, i1, i2, i3;
	double xmax, h2o_conc;
	string fname;
	ofstream output;

	fname = output_path + "sim_ice_comp.txt";
	output.open(fname.c_str(), ios::app);
	
	output << scientific;
	output.precision(5);
	output << left << setw(13) << var;
	
	h2o_conc = NV_Ith_S(y, network->find_specimen("*H2O")) + DBL_EPSILON;

	output.precision(3);
	output << left << setw(12) << NV_Ith_S(y, network->find_specimen("*CO"))/h2o_conc
		<< setw(12) << NV_Ith_S(y, network->find_specimen("*CO2"))/h2o_conc
		<< setw(12) << NV_Ith_S(y, network->find_specimen("*H2CO"))/h2o_conc
		<< setw(12) << NV_Ith_S(y, network->find_specimen("*CH3OH"))/h2o_conc
		<< setw(12) << NV_Ith_S(y, network->find_specimen("*CH4"))/h2o_conc
		<< setw(12) << NV_Ith_S(y, network->find_specimen("*NH3"))/h2o_conc
		<< setw(12) << NV_Ith_S(y, network->find_specimen("*N2"))/h2o_conc
		<< setw(12) << h2o_conc;
	
	// Main mitrogen-bearing species:
	xmax = 0.;
	for (i = network->nb_of_species - network->nb_of_gmantle_species; i < network->nb_of_species; i++) {
		if (NV_Ith_S(y, i) > xmax && network->species[i].formula[3] > 0) {
			i1 = i;
			xmax = NV_Ith_S(y, i);
		}
	}
	
	xmax = 0.;
	for (i = network->nb_of_species - network->nb_of_gmantle_species; i < network->nb_of_species; i++) {
		if (NV_Ith_S(y, i) > xmax && network->species[i].formula[3] > 0 && i != i1) {
			i2 = i;
			xmax = NV_Ith_S(y, i);
		}
	}
	
	xmax = 0.;
	for (i = network->nb_of_species - network->nb_of_gmantle_species; i < network->nb_of_species; i++) {
		if (NV_Ith_S(y, i) > xmax && network->species[i].formula[3] > 0 && i != i1 && i != i2) {
			i3 = i;
			xmax = NV_Ith_S(y, i);
		}
	}
	
	output << left << setw(10) << "! " + network->species[i1].name << setw(10) << network->species[i2].name 
		<< setw(10) << network->species[i3].name << endl;
	output.close();
}

void create_file_nautilus_data(const string & output_path)
{
	string fname;
	ofstream output;

	fname = output_path + "structure_evolution.dat";
	output.open(fname.c_str());
	output << left << setw(14) << "!time" << setw(14) << "log(Av)" << setw(14) << "log(n)" << setw(14) << "log(Tg)" << setw(14) << "log(Td)" << endl
		<< setw(14) << "!(yr)" << setw(14) << "log(mag)" << setw(14) << "log(cm-3)" << setw(14) << "log(K)" << setw(14) << "log(K)" << endl;
	output.close();
}

void save_nautilus_data(const string & output_path, double t, double av, double n, double gt, double dt)
{
	string fname;
	ofstream output;

	fname = output_path + "structure_evolution.dat";
	output.open(fname.c_str(), ios::app);
	output << scientific;
	output.precision(5);

	output << left << setw(14) << t << setw(14) << log10(av) << setw(14) << log10(n) << setw(14) << log10(gt) << setw(14) 
		<< log10(dt) << endl;
	output.close();
}

void create_file_mhd_vode(const string &output_path)
{
	string fname;
	ofstream output;

	fname = output_path + "sim_mhd_vode_data.txt";
	output.open(fname.c_str());
	output.close();
}

void save_mhd_vode(const string & input_path, const string & output_path, const chem_network *network, const N_Vector &y, 
	double conc_h_tot, double var)
{
	const int line_width = 240;
	// abundances in dust grains:
	const double mg = 3.5e-5, fe = 3.0e-5, si = 3.5e-5, ox = 4.e-4, ca = 1.5e-4; 
	char text_line[line_width];

	int i, j;
	double abund;
	// elemental abundances in full list of chemical species, and in the list of the code by Flower & Pineau des Forets (2015);
	double elabund_tot[NB_OF_CHEM_ELEMENTS], elabund_fl[NB_OF_CHEM_ELEMENTS];

	string f1, f2, sname, str, code;
	stringstream ss;

	ifstream input; 
	ofstream output;

	f1 = input_path + "chemistry/species_mhd_vode_template.in";
	input.open(f1.c_str());

	f2 = output_path + "sim_mhd_vode_data.txt";
	output.open(f2.c_str(), ios::app);
	output << scientific;
	output.precision(3);

	memset(elabund_tot, 0, NB_OF_CHEM_ELEMENTS*sizeof(double));
	memset(elabund_fl, 0, NB_OF_CHEM_ELEMENTS*sizeof(double));

	// abundance of elements in the full list of chemical species:
	for (j = 0; j < network->nb_of_species; j++) {
		for (i = 0; i < NB_OF_CHEM_ELEMENTS; i++) {
			elabund_tot[i] += NV_Ith_S(y, j)*network->species[j].formula[i];
		}
	}
	for (i = 0; i < NB_OF_CHEM_ELEMENTS; i++) {
		elabund_tot[i] /= conc_h_tot;
	}

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open template file " << f1 << endl;
	}
	else
	{
		output << var << endl;
		while (!input.eof())
		{
			do {
				input.getline(text_line, line_width);
				if (text_line[0] == '!')
					output << text_line << endl;
			}
			while (text_line[0] == '!');

			if (text_line[0] == '\0')
				break;

			ss.clear();
			ss.str(text_line);

			ss >> i >> sname >> code >> abund;
			output << right << setw(3) << i << "  ";
			output << left << setw(9) << sname << setw(15) << code;

			str = sname;
			if (sname.size() >= 2) {
				if (sname[sname.size()-1] == '*' && sname[sname.size()-2] != '*') 
				{
					str = "*";
					for (j = 0; j < (int) sname.size()-1; j++) {
						str += sname[j];
					}
				}
			}		
			if (sname == "O**")
				abund = ox;
			else if (sname == "Si**")
				abund = si;
			else if (sname == "Mg**")
				abund = mg;
			else if (sname == "C**")
				abund = ca;
			else if (sname == "Fe**")
				abund = fe;
			else abund = 1.e-99;

			// some species presented in the mhd_vode list may be absent in the current specimen list;
			j = network->find_specimen(str);
			if (j >= 0) 
			{
				abund = NV_Ith_S(y, j)/conc_h_tot;
				if (abund < 1.e-99)
					abund = 1.e-99;
				
				for (i = 0; i < NB_OF_CHEM_ELEMENTS; i++) {
					elabund_fl[i] += abund*network->species[j].formula[i];
				}
			}

			// be carefull, in windows the float format is x.xe-xxx, in linux - x.xe-xx
			output << left << setw(11) << abund;
			if (sname[sname.size() - 1] != '*') {
				ss >> code;
				output << right << setw(8) << code;
			}
			else output << "        ";

			for (j = 48; text_line[j] != '\0'; j++){
				output << text_line[j];
			}
			output << endl;
		}
		output << left << "Missing (gas-phase, ice mantles):               O = " << elabund_tot[4] - elabund_fl[4] 
			<< " C = " << elabund_tot[2] - elabund_fl[2]
			<< " N = " << elabund_tot[3] - elabund_fl[3] << endl;
		
		elabund_fl[2] += ca;
		elabund_fl[4] += ox;
		elabund_fl[5] += si;
		elabund_fl[7] += fe;
		elabund_fl[9] += mg;
	
		output << left << "Elemental abundances (for given above species): O = " << elabund_fl[4] << " C = " << elabund_fl[2] 
			<< " N = " << elabund_fl[3] << " Na = " << elabund_fl[8] << " Mg = " << elabund_fl[9] << " S = " << elabund_fl[6] 
			<< " Si = " << elabund_fl[5] << " Fe = " << elabund_fl[7] << endl;
	}
	input.close();
	output.close();
}

void create_file_mol_data(const string & output_path, const evolution_data *user_data)
{
	int i, nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, nb_lev_ci, nb_lev_oi, nb_lev_cii;
	double shock_vel, vturb;
	string fname;
	stringstream comm1, comm2;
	ofstream output;

	user_data->get_nb_of_levels(nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, nb_lev_ci, 
		nb_lev_oi, nb_lev_cii);
	
	shock_vel = user_data->get_shock_speed();
	vturb = user_data->get_vel_turb();

	comm1.clear();
	comm2.clear();
	
	comm1 << scientific;
	comm1.precision(4);
	comm1 << left << "! three parameters are given: shock speed (cm/s), turbulent velocity (cm/s), number of specimen levels: \n"
		<< setw(13) << shock_vel << setw(13) << vturb;
	
	comm2 << left << setw(18) << "!depth(cm)" << setw(13) << "gas_temp(K)" << setw(13) << "el_temp" << setw(13) << "dust_temp(K)" 
		<< setw(13) << "gasvel(cm/s)" << setw(13) << "concHe(cm-3)" << setw(13) << "conc_pH2" << setw(13) << "conc_oH2" 
		<< setw(13) << "conc_H" << setw(13) << "conc_e" << setw(13) << "conc_H_nucl" << setw(13) << "conc_mol";
	
	fname = output_path + "sim_data_h2.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_h2 << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_h2; i++) {
		output << left << setw(13) << i;
	} // no '\n' symbol at the end 
	output.close();

	fname = output_path + "sim_data_ph2o.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_h2o << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_h2o; i++) {
		output << left << setw(13) << i;
	}
	output.close();

	fname = output_path + "sim_data_oh2o.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_h2o << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_h2o; i++) {
		output << left << setw(13) << i;
	}
	output.close();

	fname = output_path + "sim_data_co.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_co << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_co; i++) {
		output << left << setw(13) << i;
	}
	output.close();

#if (CALCULATE_POPUL_NH3_OH)
	fname = output_path + "sim_data_oh.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_oh << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_oh; i++) {
		output << left << setw(13) << i;
	}
	output.close();

	fname = output_path + "sim_data_pnh3.txt";
	output.open(fname.c_str());

	output << comm1.str() << nb_lev_pnh3 << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_pnh3; i++) {
		output << left << setw(13) << i;
	}
	output.close();

	fname = output_path + "sim_data_onh3.txt";
	output.open(fname.c_str());

	output << comm1.str() << nb_lev_onh3 << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_onh3; i++) {
		output << left << setw(13) << i;
	}
	output.close();
#endif

#if (CALCULATE_POPUL_METHANOL)
	fname = output_path + "sim_data_ch3oh_a.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_ch3oh << endl;
	output << comm2.str();
	
	for (i = 1; i <= nb_lev_ch3oh; i++) {
		output << left << setw(13) << i;
	}
	output.close();

	fname = output_path + "sim_data_ch3oh_e.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_ch3oh << endl;
	output << comm2.str();
	
	for (i = 1; i <= nb_lev_ch3oh; i++) {
		output << left << setw(13) << i;
	}
	output.close();
#endif

	fname = output_path + "sim_data_ci.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_ci << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_ci; i++) {
		output << left << setw(13) << i;
	}
	output.close();

	fname = output_path + "sim_data_oi.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_oi << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_oi; i++) {
		output << left << setw(13) << i;
	}
	output.close();

	fname = output_path + "sim_data_cii.txt";
	output.open(fname.c_str());
	
	output << comm1.str() << nb_lev_cii << endl;
	output << comm2.str();

	for (i = 1; i <= nb_lev_cii; i++) {
		output << left << setw(13) << i;
	}
	output.close();
}

void save_mol_data(const string & output_path, const evolution_data *user_data, const N_Vector &y, double var)
{
	int i, nb, nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd, nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, 
		nb_lev_onh3, nb_lev_ch3oh, nb_lev_ci, nb_lev_oi, nb_lev_cii, nb_of_species;
	double a, conc_mol, conc_ph2;

	string fname;
	stringstream ss;
	ofstream output;

	const chem_network *network 
		= user_data->get_network();
	
	user_data->get_nb_of_levels(nb_lev_h2, nb_lev_h2o, nb_lev_co, nb_lev_oh, nb_lev_pnh3, nb_lev_onh3, nb_lev_ch3oh, nb_lev_ci, 
		nb_lev_oi, nb_lev_cii);
	user_data->get_nbs(nb_of_grain_charges, nb_of_equat, nb_dct, nb_mhd);
	
	nb_of_species = user_data->get_nb_of_species();
	conc_ph2 = user_data->calc_conc_ph2(y);

	ss.clear();
	ss << scientific;
	ss.precision(10);
	ss << left << setw(18) << var;
	
	ss.precision(4);
	ss << left << setw(13) << NV_Ith_S(y, nb_mhd) 
		<< setw(13) << NV_Ith_S(y, nb_mhd+2)
		<< setw(13) << user_data->get_av_dust_temp() 	// average temperature of large grains (r > MIN_ADSORPTION_RADIUS),
		<< setw(13) << NV_Ith_S(y, nb_mhd+3)			// neutral speed,
		<< setw(13) << NV_Ith_S(y, network->find_specimen("He"))
		<< setw(13) << conc_ph2
		<< setw(13) << NV_Ith_S(y, network->find_specimen("H2")) - conc_ph2
		<< setw(13) << NV_Ith_S(y, network->find_specimen("H"))
		<< setw(13) << NV_Ith_S(y, network->find_specimen("e-"))
		<< setw(13) << user_data->calc_conc_h_tot(y);

	// H2 molecule
	nb = nb_of_species;
	conc_mol = NV_Ith_S(y, network->find_specimen("H2"));

	fname = output_path + "sim_data_h2.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;
	
	for (i = 0; i < nb_lev_h2; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a; // level population densities are saved, [cm-3]
	}
	output.close();
	nb += nb_lev_h2;

	// para-H2O molecule
	conc_mol = 0.;
	for (i = 0; i < nb_lev_h2o; i++) {
		conc_mol += NV_Ith_S(y, nb + i);
	}

	fname = output_path + "sim_data_ph2o.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;
	
	for (i = 0; i < nb_lev_h2o; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
	nb += nb_lev_h2o;

	// ortho-H2O molecule
	conc_mol = 0.;
	for (i = 0; i < nb_lev_h2o; i++) {
		conc_mol += NV_Ith_S(y, nb + i);
	}

	fname = output_path + "sim_data_oh2o.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;
	
	for (i = 0; i < nb_lev_h2o; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
	nb += nb_lev_h2o;

	// CO molecule
	conc_mol = NV_Ith_S(y, network->find_specimen("CO"));

	fname = output_path + "sim_data_co.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;
	
	for (i = 0; i < nb_lev_co; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
	nb += nb_lev_co;
	
	// OH molecule
#if (CALCULATE_POPUL_NH3_OH)
	conc_mol = NV_Ith_S(y, network->find_specimen("OH"));

	fname = output_path + "sim_data_oh.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;

	for (i = 0; i < nb_lev_oh; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
#endif
	nb += nb_lev_oh;

#if (CALCULATE_POPUL_NH3_OH)
	conc_mol = 0.;
	for (i = 0; i < nb_lev_pnh3; i++) {
		conc_mol += NV_Ith_S(y, nb + i);
	}

	fname = output_path + "sim_data_pnh3.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;

	for (i = 0; i < nb_lev_pnh3; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE * conc_mol)
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
#endif
	nb += nb_lev_pnh3;
	
#if (CALCULATE_POPUL_NH3_OH)
	conc_mol = 0.;
	for (i = 0; i < nb_lev_onh3; i++) {
		conc_mol += NV_Ith_S(y, nb + i);
	}

	fname = output_path + "sim_data_onh3.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;

	for (i = 0; i < nb_lev_onh3; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE * conc_mol)
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
#endif	
	nb += nb_lev_onh3;

	// CH3OH molecule
#if (CALCULATE_POPUL_METHANOL)	
	conc_mol = 0.;
	for (i = 0; i < nb_lev_ch3oh; i++) {
		conc_mol += NV_Ith_S(y, nb + i);
	}
	fname = output_path + "sim_data_ch3oh_a.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;

	for (i = 0; i < nb_lev_ch3oh; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
#endif
	nb += nb_lev_ch3oh;

#if (CALCULATE_POPUL_METHANOL)
	conc_mol = 0.;
	for (i = 0; i < nb_lev_ch3oh; i++) {
		conc_mol += NV_Ith_S(y, nb + i);
	}
	fname = output_path + "sim_data_ch3oh_e.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;
	
	for (i = 0; i < nb_lev_ch3oh; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
#endif
	nb += nb_lev_ch3oh;
	
	// CI ion
	conc_mol = NV_Ith_S(y, network->find_specimen("C"));

	fname = output_path + "sim_data_ci.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;

	for (i = 0; i < nb_lev_ci; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
	nb += nb_lev_ci;
	
	// OI ion
	conc_mol = NV_Ith_S(y, network->find_specimen("O"));

	fname = output_path + "sim_data_oi.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;

	for (i = 0; i < nb_lev_oi; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
	nb += nb_lev_oi;
	
	// CII ion
	conc_mol = NV_Ith_S(y, network->find_specimen("C+"));

	fname = output_path + "sim_data_cii.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(4);
	output << left << endl << ss.str() << setw(13) << conc_mol;

	for (i = 0; i < nb_lev_cii; i++) {
		a = NV_Ith_S(y, nb + i);
		if (a < MINIMAL_ABUNDANCE *conc_mol) 
			a = 0.;
		output << left << setw(13) << a;
	}
	output.close();
}

void create_file_h2_chemistry(const string & output_path)
{
	int i;
	string fname;
	ofstream output;

	fname = output_path + "sim_data_h2_chemistry.txt";
	output.open(fname.c_str());

	output << left << "! o/p-H2 - ortho/para ratio of molecular hydrogen;" << endl
		<< "! o_hcoll - ortho-H2 formation rate due to collisions with H atoms, para-H2 rate = -ortho-H2 rate, [cm-3 s-1]" << endl
		<< "! h2_gr - H2 formation rate on grains, [cm-3 s-1]" << endl
		<< "! h2_gasf - H2 formation rate due to gas-phase chemistry, [cm-3 s-1]" << endl
        << "! h2_gasd - H2 destruction rate due to gas-phase chemistry, [cm-3 s-1]" << endl
		<< "! h2_h_diss - H2 dissociation rate in H2-H collisions, [cm-3 s-1]" << endl
        << "! h2_h2_diss - H2 dissociation rate in H2-H2 collisions, [cm-3 s-1]" << endl
        << "! h2_e_diss - H2 dissociation rate in H2-e collisions, [cm-3 s-1]" << endl
        << "! h2_i_diss - H2 dissociation rate by ions using the method by Wilgenbus et al. (2000)" << endl;
	
	output << left << setw(5) << "!";
	for (i = 0; i < 11; i++) {
		output << left << setw(13) << i;
	}
	output << endl << left << setw(18) << "! depth(cm)" << setw(13) << "o/p-H2" << setw(13) << "o_hcoll" << setw(13) << "h2_gr" 
		<< setw(13) << "h2_gasf" << setw(13) << "h2_gasd" << setw(13) << "h2_h_diss" << setw(13) << "h2_h2_diss" 
        << setw(13) << "vh2_vh2_diss" << setw(13) << "h2_e_diss" << setw(13) << "h2_i_diss" << endl;
	output.close();
}

void save_file_h2_chemistry(const string & output_path, const evolution_data *user_data, const N_Vector &y, double var)
{
	double op_h2_ratio, h2_form_gr, h2_form_gas, h2_destr_gas, oh2_form_hcoll, h2_h_diss, h2_ion_diss, h2_h2_diss, 
        vh2_vh2_diss, h2_e_diss;
	string fname;
	ofstream output;

	const chem_network *network 
		= user_data->get_network();

    // ortho-para ratio of H2
	op_h2_ratio = NV_Ith_S(y, network->h2_nb)/user_data->calc_conc_ph2(y) - 1.;
	user_data->get_h2_chem(h2_form_gr, h2_form_gas, h2_destr_gas, oh2_form_hcoll, h2_h_diss, h2_h2_diss, vh2_vh2_diss, 
        h2_e_diss, h2_ion_diss, y);

	fname = output_path + "sim_data_h2_chemistry.txt";
	output.open(fname.c_str(), ios::app);

	output << scientific;
	output.precision(10);
	output << left << setw(18) << var;
	
	output.precision(4);
	output << left << setw(13) << op_h2_ratio << setw(13) << oh2_form_hcoll << setw(13) << h2_form_gr 
		<< setw(13) << h2_form_gas << setw(13) << h2_destr_gas  << setw(13) << h2_h_diss << setw(13) << h2_h2_diss 
        << setw(13) << vh2_vh2_diss << setw(13) << h2_e_diss << setw(13) << h2_ion_diss << endl;
	output.close();
}

// Auxularu functions to print CVODE statistics;
static void print_stats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  cout << "CVODE Statistics:" << endl
	 << "	Nb of internal steps = " << nst << endl
	 << "	Nb of r. h. function eval. = " << nfe << endl
	 << "	Nb of calls to lin. solver setup function = " << nsetups << endl
	 << "	Nb of function eval. for jacobian calc. = " << nfeLS << endl
	 << "	Nb of jacobian eval. = " << nje << endl
	 << "	Nb of nonlinear solver iter. = " << nni << endl
	 << "	Nb of nonlinear convergence failures = " << ncfn << endl
	 << "	Nb of local error test failures = " << netf << endl
	 << "	Nb of calls to root function = " << nge << endl;
}

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

void cooling(const string &input_data)
{
	int isotope, nb_lev_h2o, nb_lev_oi, nb_lev_ci, nb_lev_cii;
	int *indices(0);
	double mass, t1, t2, cn1, cn2, cn3, cn4, ce1, ce2, ce3, ce4;
	double *oi_pop(0), *ci_pop(0), *cii_pop(0), *h2o_pop(0), *coll_partn_conc(0);

	mass = 18.*ATOMIC_MASS_UNIT;
	molecule h2o_mol("h2o", 1, mass, 0);

	nb_lev_h2o = 413; // para - 413, ortho - 411
	h2o_diagram *h2o_di =
		new h2o_diagram(input_data, h2o_mol, nb_lev_h2o);

	h2o_einstein_coeff *h2o_einst
		= new h2o_einstein_coeff(input_data, h2o_di);

	h2o_collisions *h2o_coll =
		new h2o_collisions(input_data, h2o_di);

	mass = 16.*ATOMIC_MASS_UNIT;
	molecule ion_OI("OI", isotope = 1, mass);

	mass = 12.*ATOMIC_MASS_UNIT;
	molecule ion_CI("CI", isotope = 1, mass);
	molecule ion_CII("CII", isotope = 1, mass);

	nb_lev_oi = 10;
	ion_diagram *OI_di
		= new ion_diagram(input_data, ion_OI, nb_lev_oi);

	nb_lev_ci = 3;
	ion_diagram *CI_di
		= new ion_diagram(input_data, ion_CI, nb_lev_ci);

	nb_lev_cii = 2;
	ion_diagram *CII_di
		= new ion_diagram(input_data, ion_CII, nb_lev_cii);

	ion_einstein_coeff *OI_einst
		= new ion_einstein_coeff(input_data, OI_di);

	ion_einstein_coeff *CI_einst
		= new ion_einstein_coeff(input_data, CI_di);

	ion_einstein_coeff *CII_einst
		= new ion_einstein_coeff(input_data, CII_di);

	OI_collisions *OI_coll
		= new OI_collisions(input_data, OI_di);

	CI_collisions *CI_coll
		= new CI_collisions(input_data, CI_di);

	CII_collisions *CII_coll
		= new CII_collisions(input_data, CII_di);

	oi_pop = new double [nb_lev_oi];
	ci_pop = new double [nb_lev_ci];
	cii_pop = new double [nb_lev_cii];
	h2o_pop = new double [nb_lev_h2o];

	cout << scientific;
	cout.precision(2);

	t1 = 200;
	for (t2 = 10.; t2 < 10000.; t2 *= 1.1)
	{
		OI_coll->set_gas_param(t1, t2, 25., 30., 90., 1., 0.0015, coll_partn_conc, indices);
		opt_thin_pop(oi_pop, OI_di, OI_einst, OI_coll, t1, t2, coll_partn_conc, indices);

		cn1 = 0.08*heating_of_neutral_gas(oi_pop, OI_di, OI_coll, t1, coll_partn_conc, indices);
		ce1 = 0.08*heating_of_electron_gas(oi_pop, OI_di, OI_coll, t2, coll_partn_conc, indices);

		CI_coll->set_gas_param(t1, t2, 25., 30., 90., 1., 0.0015, coll_partn_conc, indices);
		opt_thin_pop(ci_pop, CI_di, CI_einst, CI_coll, t1, t2, coll_partn_conc, indices);
		
		cn2 = 0.03*heating_of_neutral_gas(ci_pop, CI_di, CI_coll, t1, coll_partn_conc, indices);
		ce2 = 0.03*heating_of_electron_gas(ci_pop, CI_di, CI_coll, t2, coll_partn_conc, indices);

		CII_coll->set_gas_param(t1, t2, 25., 30., 90., 1., 0.0015, coll_partn_conc, indices);
		opt_thin_pop(cii_pop, CII_di, CII_einst, CII_coll, t1, t2, coll_partn_conc, indices);

		cn3 = 1e-3*heating_of_neutral_gas(cii_pop, CII_di, CII_coll, t1, coll_partn_conc, indices);
		ce3 = 1e-3*heating_of_electron_gas(cii_pop, CII_di, CII_coll, t2, coll_partn_conc, indices);

		h2o_coll-> set_gas_param(t1, t2, 25., 30., 90., 1., 0.0015, coll_partn_conc, indices);
		opt_thin_pop(h2o_pop, h2o_di, h2o_einst, h2o_coll, t1, t2, coll_partn_conc, indices);

		cn4 = 1e-2*heating_of_neutral_gas(h2o_pop, h2o_di, h2o_coll, t1, coll_partn_conc, indices);
		ce4 = 1e-2*heating_of_electron_gas(h2o_pop, h2o_di, h2o_coll, t2, coll_partn_conc, indices);

	//	CII_coll->get_rate_electrons(CII_di->lev_array[1], CII_di->lev_array[0], a, b, t, coll_partn_conc, indices);
	//	CII_coll->get_rate_neutrals(CII_di->lev_array[1], CII_di->lev_array[0], c, b, t, coll_partn_conc, indices);
	//	cout << left << setw(13) << t << setw(13) << a << setw(13) << c << endl;
		
		cout << left << setw(11) << t2 << setw(11) << cn1 << setw(13) << ce1 << setw(11) << cn2 << setw(13) << ce2 
			<< setw(11) << cn3 << setw(13) << ce3 << setw(13) << cn4 << setw(13) << ce4 << endl;	
	}
}

void sputtering()
{
	int i;
	double rt, s, v;

	// projectile has the first place:
	sputtering_yield *sp_yield 
		= new sputtering_yield(2.*ATOMIC_MASS_UNIT, 18.*ATOMIC_MASS_UNIT, 4200.);
	
	reaction_rate_data *rate_data 
		= new reaction_rate_data();
					
	// rate_data->calc_data(sp_yield);
	rate_data->calc_data(sp_yield);

	rt = 1000.; // in K;
	v = sqrt(2.*rt*BOLTZMANN_CONSTANT/(2.*ATOMIC_MASS_UNIT));
	
	for (i = 0; i < 10; i++)
	{
		s = i;
		cout << left << setw(13) << s << setw(13) << v*s << setw(13) << rate_data->get(rt, s) << endl;
	}
}

/*output << "! abcde = abc*10^(-de) or -abcde = abc*10^(de)" << endl;
    output << left << setw(13) << "!t(yr)/z(cm)";
    for (i = 0; i < network->nb_of_reactions; i++) {
        output << left << setw(6) << i + 1;
    }
    output << endl;
    output.close();*/
    /*	output.setf(ios::scientific, ios::floatfield);
        output.precision(5);
        output << left << setw(13) << var;

        output.setf(ios::fixed, ios::floatfield);
        output.precision(2);

        for (i = 0; i < user_data->get_reaction_nb(); i++) {
            a = user_data->get_reaction_rate(i);
            if (a > 1000*MINIMAL_REACTION_RATE)
            {
                if (a < 1.)
                    j = ((int) log10(a)) - 1 - 2;
                else j = (int) log10(a) - 2;

                k = rounding(a*pow(10., -j));
                if (k == 1000) {
                    k = 100;
                    j++;
                }
                k *= ((j < 0) ? 1 : -1);
                j = abs(j);

                output << left << setw(3) << k;

                if (j >= 10)
                    output << setw(3) << abs(j);
                else if (j > 0)
                    output << "0" << j << " ";
                else output << "00 ";
            }
            else output << left << setw(6) << "0     ";
        }
        output << endl;
        output.close();*/

