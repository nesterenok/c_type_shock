
// 03.10.2017. Check for errors.
// 12.02.2018. Check for errors.
#include <stdio.h>
#include <stdlib.h>

#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <float.h>
#include <algorithm>
#include <sstream>

#include "utils.h"
#include "reaction_rate_analyzer.h"

#define MAX_TEXT_LINE_WIDTH 320
#define RRA_NB_REACTANTS 2
#define RRA_NB_PRODUCTS 4
#define RRA_NB_RATE_VALUES 1500
#define CHEM_ANALYSIS_FOLDER "chem_analysis/"
using namespace std;

const double min_fraction = 1.e-2; // in looking for the important reactions
const double low_err = 0.3; // in comparison between simulated abundances and NAUTILUS data, err = fabs(log10(x1/x2));

int reaction_data::nb_of_reactants = RRA_NB_REACTANTS; 
int reaction_data::nb_of_products = RRA_NB_PRODUCTS; 
int reaction_data::max_nb_of_rate_values = RRA_NB_RATE_VALUES; 

reaction_data::reaction_data() : name("")
{
	reactants = new int [nb_of_reactants];
	products = new int [nb_of_products];
	rates = new double [max_nb_of_rate_values];
	memset(rates, 0, max_nb_of_rate_values*sizeof(double));
}

reaction_data::reaction_data(const reaction_data & obj)
{
	name = obj.name;
	reactants = new int [nb_of_reactants];
	memcpy(reactants, obj.reactants, nb_of_reactants*sizeof(int));

	products = new int [nb_of_products];
	memcpy(products, obj.products, nb_of_products*sizeof(int));

	rates = new double [max_nb_of_rate_values];
	memcpy(rates, obj.rates, max_nb_of_rate_values*sizeof(double));
}

reaction_data & reaction_data::operator=(const reaction_data &obj)
{
	if (this == &obj)
		return *this;
	name = obj.name;
	
	delete [] reactants;
	reactants = new int [nb_of_reactants];
	memcpy(reactants, obj.reactants, nb_of_reactants*sizeof(int));

	delete [] products;
	products = new int [nb_of_products];
	memcpy(products, obj.products, nb_of_products*sizeof(int));

	delete [] rates;
	rates = new double [max_nb_of_rate_values];
	memcpy(rates, obj.rates, max_nb_of_rate_values*sizeof(double));
	
	return *this;
}

reaction_data::~reaction_data()
{
	delete [] reactants;
	delete [] products;
	delete [] rates;
}

//
// Friend functions of class reaction_data
//
bool operator == (const reaction_data & obj1, const reaction_data & obj2) {
	return (obj1.name == obj2.name);
}

bool operator != (const reaction_data & obj1, const reaction_data & obj2)
{ return !(obj1 == obj2); }

bool operator > (const reaction_data & obj1, const reaction_data & obj2)
{
	int i;
	double maxr(0.);

	if (obj1 == obj2)
		return false;

	for (i = 0; i < obj1.max_nb_of_rate_values; i++) {
		if (obj1.rates[i] > maxr)
			maxr = obj1.rates[i];
	}
	for (i = 0; i < obj2.max_nb_of_rate_values; i++) {
		if (obj2.rates[i] > maxr)
			return false;
	}
	return true;
}

bool operator < (const reaction_data & obj1, const reaction_data & obj2)
{
	if (obj1 == obj2) 
		return false;
	return !(obj1 > obj2);
}

/*    specimen_names.push_back("HOCO");
    specimen_names.push_back("HNCO");
    specimen_names.push_back("OCN");

    specimen_names.push_back("CH3OCH2");
    specimen_names.push_back("CH3OCH3");
    specimen_names.push_back("HCOOH");
    specimen_names.push_back("C2H5OH");
    specimen_names.push_back("CH3CHO");
    specimen_names.push_back("C3N");
    specimen_names.push_back("CN");
    specimen_names.push_back("C2");
    specimen_names.push_back("C6H");

    specimen_names.push_back("HCOOH2+");
    specimen_names.push_back("CH3OCH4+");
    specimen_names.push_back("C3N-");
    specimen_names.push_back("C6H-");

    specimen_names.push_back("CH3OH2+");
    specimen_names.push_back("C2O");
    specimen_names.push_back("HC2O");
    specimen_names.push_back("CH2CO");
    specimen_names.push_back("CH3CO");
    specimen_names.push_back("C2H4");
    specimen_names.push_back("C2H5");
    specimen_names.push_back("HCOOCH3");
    specimen_names.push_back("H5C2O2+");
    specimen_names.push_back("C2H5OH2+");*/

void add_cosmicray_chemistry(std::vector<string>& specimen_names) {
    specimen_names.push_back("e-");
    specimen_names.push_back("H2+");
    specimen_names.push_back("H3+");
    specimen_names.push_back("H2O+");
    specimen_names.push_back("OH+");
    specimen_names.push_back("H3O+");
    specimen_names.push_back("C+");
    specimen_names.push_back("HCO+");
    specimen_names.push_back("NH4+");
    specimen_names.push_back("N2H+");
}

void add_oxygen_chemistry(std::vector<string>& specimen_names) {
    specimen_names.push_back("H");
    specimen_names.push_back("O");
    specimen_names.push_back("H2");
    specimen_names.push_back("OH");
    specimen_names.push_back("H2O");
    specimen_names.push_back("O2");
}

void add_carbon_chemistry(std::vector<string>& specimen_names) {
    specimen_names.push_back("C");
    specimen_names.push_back("CO");
    specimen_names.push_back("CO2");
    specimen_names.push_back("HCO");
    specimen_names.push_back("H2CO");
    specimen_names.push_back("CH3O");
    specimen_names.push_back("CH2OH");
    specimen_names.push_back("CH3OH");
    specimen_names.push_back("CH4");
    specimen_names.push_back("CH");
}

void add_nitrogen_chemistry(std::vector<string>& specimen_names) {
    specimen_names.push_back("N2");
    specimen_names.push_back("NO");
    specimen_names.push_back("NH2");
    specimen_names.push_back("NH3");
    specimen_names.push_back("HNC");
    specimen_names.push_back("HCN");
}

void production_routes(string path1, string path2)
{
	bool bo;
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, k, l, nb, nb_of_rate_values, nb_of_species, nb_of_reactions, specimen_nb, nb_of_elem;
    int16_t abc;
    int8_t de;
	double z_arr[RRA_NB_RATE_VALUES], t_arr[RRA_NB_RATE_VALUES], destr_rate[RRA_NB_RATE_VALUES], prod_rate[RRA_NB_RATE_VALUES];

	string str, sn, fn;
	reaction_data rd;

	vector<string> sn_v;
	vector<reaction_data> rd_v;
	vector<reaction_data> destr_rd_v;
	vector<reaction_data> prod_rd_v;
	
	ifstream input;
	ofstream output;
    vector<string> specimen_names;
 
    add_cosmicray_chemistry(specimen_names);
    add_oxygen_chemistry(specimen_names);
    add_carbon_chemistry(specimen_names);
    add_nitrogen_chemistry(specimen_names);

    vector<string>::iterator it;
    it = unique(specimen_names.begin(), specimen_names.end());
    specimen_names.resize(std::distance(specimen_names.begin(), it));

    nb = (int) specimen_names.size();

    fn = path1 + "sim_species.txt";
	input.open(fn.c_str(), ios::in);
	input >> nb_of_species >> nb_of_elem;

	for (i = 0; i < nb_of_species; i++) {
		input >> j >> sn;
		for (j = 0; j < nb_of_elem; j++) {
			input >> k;
		}
		sn_v.push_back(sn);
	}
	input.close();
	
	fn = path1 + "sim_reactions.txt";
	input.open(fn.c_str(), ios::in);
	input >> nb_of_reactions;

	for (i = 0; i < nb_of_reactions; i++) 
	{
		input >> j;
		for (j = 0; j < RRA_NB_REACTANTS; j++) {
			input >> rd.reactants[j];
		}
		for (j = 0; j < RRA_NB_PRODUCTS; j++) {
			input >> rd.products[j];
		}
		getline(input, rd.name);
		rd_v.push_back(rd);
	}
	input.close();

    // reading binary file
	fn = path2 + "sim_reaction_rates.bin";
	input.open(fn.c_str(), ios::binary | ios::in);
    
    i = 0;
    while (!input.eof() && i < RRA_NB_RATE_VALUES)
    {
        input.read((char*)& abc, sizeof(abc));
        if (input.eof()) {
            break;
        }
        input.read((char*)& de, sizeof(de));
        
        z_arr[i] = abc *pow(10., (int) de);

        input.read((char*)& abc, sizeof(abc));
        input.read((char*)& de, sizeof(de));

        t_arr[i] = abc * pow(10., (int) de);

        for (j = 0; j < nb_of_reactions; j++) {
            input.read((char*)& abc, sizeof(abc));
            input.read((char*)& de, sizeof(de));

            rd_v[j].rates[i] = abc * pow(10., (int) de);
        }
        i++;
    }
    input.close();
	nb_of_rate_values = i;

	for (k = 0; k < 2*nb; k++)
	{
		sn = specimen_names[k/2];
		if (k%2 == 0)
			sn = '*' + sn;
		
		destr_rd_v.clear();
		prod_rd_v.clear();

		specimen_nb = -1;
		for (i = 0; i < nb_of_species; i++) {
			if (sn_v[i] == sn) {
				specimen_nb = i;
				break;
			}
		}
		if (specimen_nb == -1)
			cout << "There is not this specimen in the list: " << sn << endl;
		else 
		{
            // formation of the list of production and destruction reactions
			for (i = 0; i < nb_of_reactions; i++) {
                l = 0;
				for (j = 0; j < RRA_NB_REACTANTS; j++) 
				{
					if (rd_v[i].reactants[j] == specimen_nb) {
                        l--;
					}
				}
				for (j = 0; j < RRA_NB_PRODUCTS; j++)
				{
					if (rd_v[i].products[j] == specimen_nb) {
                        l++;
					}
				}
                if (l < 0)
                    destr_rd_v.push_back(rd_v[i]);
                else if (l > 0) prod_rd_v.push_back(rd_v[i]);
			}
	
			for (i = 0; i < nb_of_rate_values; i++) 
			{
				destr_rate[i] = prod_rate[i] = 0.;
				for (j = 0; j < (int) destr_rd_v.size(); j++) {
					destr_rate[i] += destr_rd_v[j].rates[i];
				}
				for (j = 0; j < (int) prod_rd_v.size(); j++) {
					prod_rate[i] += prod_rd_v[j].rates[i];
				}
			}

			l = 0;
			sort(destr_rd_v.begin(), destr_rd_v.end());
			reverse(destr_rd_v.begin(), destr_rd_v.end());
			
			for (j = 0; j < (int) destr_rd_v.size(); j++) 
			{
				bo = true;
				for (i = 0; i < nb_of_rate_values; i++) 
				{
					destr_rd_v[j].rates[i] /= destr_rate[i];
					if (destr_rd_v[j].rates[i] > min_fraction || l < 10) {
						bo = false; // do not use break here;
					}
				}
				if (bo) {
					destr_rd_v.erase(destr_rd_v.begin() + j);
					j--;
				}
				else l++;
			}

			l = 0;
			sort(prod_rd_v.begin(), prod_rd_v.end());
			reverse(prod_rd_v.begin(), prod_rd_v.end());

			for (j = 0; j < (int) prod_rd_v.size(); j++) 
			{
				bo = true;
				for (i = 0; i < nb_of_rate_values; i++) 
				{
					prod_rd_v[j].rates[i] /= prod_rate[i]; 
					if (prod_rd_v[j].rates[i] > min_fraction || l < 10) {
						bo = false;
					}
				}
				if (bo) {
					prod_rd_v.erase(prod_rd_v.begin() + j);
					j--;
				}
				else l++;
			}

			// file name has not to contain '*'
			if (sn[0] == '*')
				sn[0] = 'J';

            fn = path2;
            fn += CHEM_ANALYSIS_FOLDER;
            fn += "destr_";
			fn += sn;
			fn += ".txt";

			output.open(fn.c_str(), ios::out);
			output << scientific;
	
			for (j = 0; j < (int) destr_rd_v.size(); j++) {
				output << left << "! " << setw(5) << j << destr_rd_v[j].name << endl;
			}
	
			output << left << setw(13) << "! z" << setw(10) << "gas_temp" << setw(10) << "tot_rate";
			for (j = 0; j < (int) destr_rd_v.size(); j++) {
				output << left << setw(10) << j;
			}
			output << endl;

			for (i = 0; i < nb_of_rate_values; i++) 
			{
				output.precision(5);
				output << left << setw(13) << z_arr[i];
				output.precision(2); // <=3
				output << left << setw(10) << t_arr[i] << setw(10) << destr_rate[i];
				for (j = 0; j < (int) destr_rd_v.size(); j++) {
					output << left << setw(10) << destr_rd_v[j].rates[i];
				}
				output << endl;
			}
			output.close();

            fn = path2;
            fn += CHEM_ANALYSIS_FOLDER;
            fn += "prod_";
			fn += sn;
			fn += ".txt";

			output.open(fn.c_str(), ios::out);
			output << scientific;
			
			for (j = 0; j < (int) prod_rd_v.size(); j++) {
				output << left << "! " << setw(5) << j << prod_rd_v[j].name << endl;
			}
	
			output << left << setw(13) << "! z" << setw(10) << "gas_temp" << setw(10) << "tot_rate";
			for (j = 0; j < (int) prod_rd_v.size(); j++) {
				output << left << setw(10) << j;
			}
			output << endl;

			for (i = 0; i < nb_of_rate_values; i++) 
			{
				output.precision(5);
				output << left << setw(13) << z_arr[i];
				output.precision(2); // <=3
				output << left << setw(10) << t_arr[i] << setw(10) << prod_rate[i];
				for (j = 0; j < (int) prod_rd_v.size(); j++) {
					output << left << setw(10) << prod_rd_v[j].rates[i];
				}
				output << endl;
			}
			output.close();
		}
	}
}

void nautilus_comparison(string path)
{
	const int line_width = 10000; // must be long
	char text_line[line_width];

	int i, j, k, nb_times, nb_of_species, nb_of_elem, *sp_nb_low, *sp_nb_med, *sp_nb_high;
	double vn, *ty, **rel_err, **sim_arr;

	string fn, sn;
	stringstream ss;
	vector<string> sn_v;
	
	ifstream input, input_n;
	ofstream output;

	fn = path + "sim_species.txt";
	input.open(fn.c_str(), ios::in);
	input >> nb_of_species >> nb_of_elem;

	for (i = 0; i < nb_of_species; i++) {
		input >> j >> sn;
		for (j = 0; j < nb_of_elem; j++) {
			input >> k;
		}
		if (sn[0] == '*')
			sn[0] = 'J';
		sn_v.push_back(sn);
	}
	input.close();

	sim_arr = alloc_2d_array<double>(nb_of_species, RRA_NB_RATE_VALUES);
	memset(*sim_arr, 0, nb_of_species*RRA_NB_RATE_VALUES*sizeof(double));

	ty = new double [RRA_NB_RATE_VALUES];
	memset(ty, 0, RRA_NB_RATE_VALUES*sizeof(double));

	fn = path + "sim_specimen_abund.txt";
	input.open(fn.c_str(), ios::in);
	input.getline(text_line, line_width);
	input.getline(text_line, line_width);

	nb_times = 0;
	while (!input.eof()) 
	{
		input.getline(text_line, line_width);
		if (text_line[0] == '\0')
			break;
		
		ss.clear();
		ss.str(text_line);
		ss >> ty[nb_times];
		
		for (i = 0; i < nb_of_species; i++) {		
			ss >> sim_arr[i][nb_times];
		}
		nb_times++;
	}
	input.close();

	rel_err = alloc_2d_array<double>(nb_of_species, nb_times);
	memset(*rel_err, 0, nb_of_species*nb_times*sizeof(double));

	sp_nb_low = new int [nb_times];
	memset(sp_nb_low, 0, nb_times*sizeof(int));

	sp_nb_med = new int [nb_times];
	memset(sp_nb_med, 0, nb_times*sizeof(int));

	sp_nb_high = new int [nb_times];
	memset(sp_nb_high, 0, nb_times*sizeof(int));

	for (i = 0; i < nb_of_species; i++)
	{
		fn = path + "ab/";
		fn += sn_v[i];
		fn += ".ab";
		input_n.open(fn.c_str(), ios::in);

		j = 0;
		if (input_n.is_open()) 
		{
			// comment lines are read:
			input_n.getline(text_line, line_width);		
			while (!input_n.eof() && j < nb_times) 
			{
				// Reading data from NAUTILUS data file:
				input_n.getline(text_line, line_width);
				if (text_line[0] == '\0')
					break;
				
				ss.clear();
				ss.str(text_line);
				ss >> ty[j] >> vn;
				
				rel_err[i][j] = fabs(sim_arr[i][j]/vn);

				if (fabs(log10(rel_err[i][j])) < low_err)
					sp_nb_low[j]++;
				else if (fabs(log10(rel_err[i][j])) < 1.)
					sp_nb_med[j]++;
				else 
					sp_nb_high[j]++;	
				j++;
			}
		}
		input_n.close();
	}
	
	fn = path + "sim_nautilus_comp.txt";
	output.open(fn.c_str(), ios::out);
	output << scientific;
	output.precision(2);

	output << left << setw(11) << "!time(yr)";
	for (i = 0; i < nb_of_species; i++) {
		output << left << setw(11) << sn_v[i];
	}
	output << endl;

	for (k = 0; k < nb_times; k++)
	{
		output << left << setw(11) << ty[k];
		for (i = 0; i < nb_of_species; i++) {
			output << left << setw(11) << rel_err[i][k];
		}
		output << endl;
	}

	for (j = 0; j < nb_times; j++) 
	{
		output << left << setw(12) << ty[j] << setw(8) << sp_nb_low[j];
		for (i = 0; i < nb_of_species; i++) {
			if (fabs(log10(rel_err[i][j])) < low_err)
				output << left << setw(10) << sn_v[i];
		}
		output << endl;
		
		output << left << setw(12) << "" << setw(8) << sp_nb_med[j]; 
		for (i = 0; i < nb_of_species; i++) {
			if (fabs(log10(rel_err[i][j])) >= low_err && fabs(log10(rel_err[i][j])) < 1.)
				output << left << setw(10) << sn_v[i];
		}
		output << endl;

		output << left << setw(12) << "" << setw(8) << sp_nb_high[j];
		for (i = 0; i < nb_of_species; i++) {
			if (fabs(log10(rel_err[i][j])) >= 1.)
				output << left << setw(10) << sn_v[i];
		}
		output << endl;
	}
	output.close();

	delete [] ty; 
	delete [] sp_nb_low; 
	delete [] sp_nb_med; 
	delete [] sp_nb_high; 
	free_2d_array(rel_err);
	free_2d_array(sim_arr);
}

/*
double depth_temperature_dependence::get_gas_temperature(double d)
{ 
    low = lower_bound(depth.begin() + prev_index, depth.end(), d);   
    prev_index = (low - depth.begin());
    
    if (prev_index < 1)
        prev_index = 1;
    else if (prev_index > (int) depth.size() - 1)
        prev_index = (int)depth.size() - 1;

    return gas_temperature[prev_index -1] + (d - depth[prev_index -1]) 
        * (gas_temperature[prev_index] - gas_temperature[prev_index -1]) / (depth[prev_index] - depth[prev_index -1]);
}

depth_temperature_dependence::depth_temperature_dependence(std::string path) : prev_index(0)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    double z, time, temp;

    string fn;
    stringstream ss;
    ifstream input;

    fn = path + "sim_phys_param.txt";
    input.open(fn.c_str(), ios::in);

    while (!input.eof())
    {
        // comment lines are read:
        do
            input.getline(text_line, MAX_TEXT_LINE_WIDTH);
        while (text_line[0] == '!');

        if (text_line[0] == '\0')
            break;

        ss.clear();
        ss.str(text_line);
        ss >> z >> time >> temp;

        depth.push_back(z);
        gas_temperature.push_back(temp);
    }
    input.close();
}*/
