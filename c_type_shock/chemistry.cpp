//
// 17.03.2017. Check for errors.
// 20.03.2017. Bug was found in calculating the reaction rates for A+Ion->Ion+Ion+e- and A+Ion->B+C+Ion:
//		max rate was not calculated, A coefficient was incorrect for second reaction;
// 28.09.2017. Please, check multichannel reactions with a barrier;
// 31.01.2018. Check for errors.
// 22.03.2018. Check for errors

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include "constants.h"
#include "parameters.h"
#include "utils.h"
#include "integration.h"
#include "special_functions.h"
#include "chemistry.h"
#include "chemical_reaction_data.h"

#define MAX_TEXT_LINE_WIDTH 320 // long lines
#define SOURCE_NAME "chemistry.cpp"
using namespace std;

const double onedivby_8grain_sites_per_cm2 = 0.125/GRAIN_SITES_PER_CM2;

//
// Friend functions of class chem_specimen
//
bool operator == (const chem_specimen & obj1, const chem_specimen & obj2) {
	return (obj1.name == obj2.name);
}

bool operator != (const chem_specimen & obj1, const chem_specimen & obj2)
{ return !(obj1 == obj2); }

bool operator < (const chem_specimen & obj1, const chem_specimen & obj2)
{
	if (obj1 == obj2) 
		return false;

	if (obj1.type == "adsorbed" && obj2.type != "adsorbed")
		return false;
	else if (obj1.type != "adsorbed" && obj2.type == "adsorbed")
		return true;

	// species have the same sign:
	if (obj1.charge == obj2.charge ) 
	{
		if (fabs(obj1.mass - obj2.mass) > 0.01*ATOMIC_MASS_UNIT) {
			return (obj1.mass < obj2.mass); // sort by mass
		}
		else 
		{ // sort by name:
			if (obj1.name.size() != obj2.name.size())
				return (obj1.name.size() < obj2.name.size());
			else {
				for (int i = 0; i < (int) obj1.name.size(); i++) 
				{
					if (obj1.name[i] != obj2.name[i])
						return (obj1.name[i] < obj2.name[i]);
				}
				return false;
			}
		}
	}
	if (obj1.charge == 0) 
		return true;
	
	if (obj2.charge == 0) 
		return false;
	else return (obj1.charge > obj2.charge);
}

bool operator > (const chem_specimen & obj1, const chem_specimen & obj2)
{
	if (obj1 == obj2) 
		return false;
	return !(obj1 < obj2);
}

chem_specimen::chem_specimen() : charge(0), n_nb(0), enthalpy(0.), mass(0.), bind_en(0.), surf_freq(0.), diff_barrier(0.), 
	name(""), type(""), is_enth_def(false) {
	memset(formula, 0, NB_OF_CHEM_ELEMENTS*sizeof(int));
}

//
// Chemical reaction class
//
chem_reaction::chem_reaction() : rate_data(0), reactant(0), product(0), temp_min(0), temp_max(0), parameters(0), 
	type(0), nb_of_reactants(0), nb_of_products(0), nb_of_fits(0), nb_of_param(0), cs_reconstr(0), energy_released(0.), 
	mass1(0.), mass2(0.), reduced_mass(0.), mass_sum(0.), min_rate(0.), max_rate(0.), name(""), rcode("")
{;}

chem_reaction::chem_reaction(const chem_reaction & obj)
{
	type = obj.type;
	nb_of_reactants = obj.nb_of_reactants;
	nb_of_products = obj.nb_of_products;
	nb_of_fits = obj.nb_of_fits;
	nb_of_param = obj.nb_of_param;
	cs_reconstr = obj.cs_reconstr;
	
	energy_released = obj.energy_released;
	mass1 = obj.mass1;
	mass2 = obj.mass2;
	reduced_mass = obj.reduced_mass;
	mass_sum = obj.mass_sum;
	max_rate = obj.max_rate;
	min_rate = obj.min_rate;
	
	name = obj.name;
	rcode = obj.rcode;

	reactant = new int [nb_of_reactants];
	memcpy(reactant, obj.reactant, nb_of_reactants*sizeof(int));

	product = new int [nb_of_products];
	memcpy(product, obj.product, nb_of_products*sizeof(int));

	temp_min = new double [nb_of_fits];
	temp_max = new double [nb_of_fits];

	memcpy(temp_min, obj.temp_min, nb_of_fits*sizeof(double));
	memcpy(temp_max, obj.temp_max, nb_of_fits*sizeof(double));

	parameters = new double [nb_of_param*nb_of_fits];
	memcpy(parameters, obj.parameters, nb_of_param*nb_of_fits*sizeof(double));

	rate_data = obj.rate_data;
}

chem_reaction & chem_reaction::operator=(const chem_reaction &obj)
{
	if (this == &obj)
		return *this;

	delete_data();

	type = obj.type;
	nb_of_reactants = obj.nb_of_reactants;
	nb_of_products = obj.nb_of_products;
	nb_of_fits = obj.nb_of_fits;
	nb_of_param = obj.nb_of_param;
	cs_reconstr = obj.cs_reconstr;
	
	energy_released = obj.energy_released;
	mass1 = obj.mass1;
	mass2 = obj.mass2;
	reduced_mass = obj.reduced_mass;
	mass_sum = obj.mass_sum;
	max_rate = obj.max_rate;
	min_rate = obj.min_rate;

	name = obj.name;
	rcode = obj.rcode;

	reactant = new int [nb_of_reactants];
	memcpy(reactant, obj.reactant, nb_of_reactants*sizeof(int));

	product = new int [nb_of_products];
	memcpy(product, obj.product, nb_of_products*sizeof(int));

	temp_min = new double [nb_of_fits];
	temp_max = new double [nb_of_fits];

	memcpy(temp_min, obj.temp_min, nb_of_fits*sizeof(double));
	memcpy(temp_max, obj.temp_max, nb_of_fits*sizeof(double));

	parameters = new double [nb_of_param*nb_of_fits];
	memcpy(parameters, obj.parameters, nb_of_param*nb_of_fits*sizeof(double));

	rate_data = obj.rate_data;
	return *this;
}

void chem_reaction::delete_data()
{
	type = nb_of_reactants = nb_of_products = nb_of_fits = nb_of_param = cs_reconstr = 0;
	energy_released = mass1 = mass2 = reduced_mass = mass_sum = min_rate = max_rate = 0.;
	name = rcode = "";

	// the object rate_data is not deleted, it is stored and deleted outside this class:
	rate_data = 0;
	
	delete [] reactant;
	delete [] product;

	delete [] temp_min; 
	delete [] temp_max;
	delete [] parameters;

	reactant = product = 0;
	temp_min = temp_max = parameters = 0;
}

bool operator == (const chem_reaction & r1, const chem_reaction & r2)
{
	const int max_nb_products = 5;
	bool b_arr[max_nb_products];
	int k, m;
	
	if (r1.type == r2.type && r1.nb_of_products == r2.nb_of_products && r1.nb_of_reactants == r2.nb_of_reactants)
	{
		m = 0;
		if (r1.nb_of_reactants == 1) { // unimolecular reactions:
			if (r1.reactant[0] == r2.reactant[0])
				m = 1;
		}
		else { // bimolecular reactions:
			if ((r1.reactant[0] == r2.reactant[0] && r1.reactant[1] == r2.reactant[1])
				|| (r1.reactant[1] == r2.reactant[0] && r1.reactant[0] == r2.reactant[1])) {
				m = 1;
			}
		}
		if (m == 1) // reactants coincide;
		{ 
			memset(b_arr, 0, max_nb_products*sizeof(bool));
			for (k = 0; k < r1.nb_of_products; k++) {
				for (m = 0; m < r2.nb_of_products; m++)
				{
					if (r1.product[k] == r2.product[m] && !b_arr[m]) {
						b_arr[m] = true;
						break;
					}
				}
			}
			for (k = 0; k < r1.nb_of_products; k++) {
				if (!b_arr[k])
					break;
			}
			if (k == r1.nb_of_products) {
				return true;
			}
		}
	}
	return false;
}

//
// Class chemical network
//
chem_network::chem_network(const string &path, int verb) : verbosity(verb), max_nb_carbon(11), nb_of_species(0), nb_of_gmantle_species(0), 
	nb_of_reactions(0), nb_reactions_ion_grains(0), h2_nb(0), ah2_nb(0), h_nb(0), he_nb(0), e_nb(0), oi_nb(0), ci_nb(0), cii_nb(0), 
	co_nb(0), h2o_nb(0), oh_nb(0), nh3_nb(0), ch3oh_nb(0), hp_nb(0), h3p_nb(0), h2_h_diss_nb(0), h2_h2_diss_nb(0), h2_e_diss_nb(0)
{;}

chem_network::~chem_network()
{
	elements.clear();
	species.clear();
	reaction_array.clear();
	rate_data_array.clear();
}

// The masses of elements are in atomic mass units, the values are rounded;
void chem_network::add_element(const std::string &name)
{
	chem_element elem;
	elem.name = name;

	if (name == "H") {
		elem.mass = 1.;
		elements.insert(pair<int, chem_element>(0, elem));
	}
	else if (name == "He") {
		elem.mass = 4.;
		elements.insert(pair<int, chem_element>(1, elem));
	}
	else if (name == "C") {
		elem.mass = 12.;
		elements.insert(pair<int, chem_element>(2, elem));
	}
	else if (name == "N") {
		elem.mass = 14.;
		elements.insert(pair<int, chem_element>(3, elem));
	}
	else if (name == "O") {
		elem.mass = 16.;
		elements.insert(pair<int, chem_element>(4, elem));
	}
	else if (name == "Si") {
		elem.mass = 28.;
		elements.insert(pair<int, chem_element>(5, elem));
	}
	else if (name == "S") {
		elem.mass = 32.;
		elements.insert(pair<int, chem_element>(6, elem));
	}
	else if (name == "Fe") {
		elem.mass = 56.;
		elements.insert(pair<int, chem_element>(7, elem));
	}
	else if (name == "Na") {
		elem.mass = 23.;
		elements.insert(pair<int, chem_element>(8, elem));
	}
	else if (name == "Mg") {
		elem.mass = 24.;
		elements.insert(pair<int, chem_element>(9, elem));
	}
	else if (name == "Cl") {
		elem.mass = 35.;
		elements.insert(pair<int, chem_element>(10, elem));
	}
	else if (name == "P") {
		elem.mass = 31.;
		elements.insert(pair<int, chem_element>(11, elem));
	}
	else if (name == "F") {
		elem.mass = 19.;
		elements.insert(pair<int, chem_element>(12, elem));
	}
	else if (verbosity) {
		cout << "The element can not be included in the list: " << name << endl;
		return;
	}

	if (verbosity) 
		cout << "The element was added to the network list: " << name << endl;
}

void chem_network::init_gas_phase_species(const std::string file_name)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	bool failed;
	int i, nb;
	double mass;

	ifstream input;
	stringstream ss;
	chem_specimen specimen;

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with chemical species " << file_name << endl;
		exit(1);
	}
	if (verbosity) {
		cout << "The data is read from file: " << file_name << endl;
	}

	while (!input.eof())
	{
		// comment lines are read:
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		if (text_line[0] == '\0') // check for empty line at the file end;
			break;

		ss.clear();
		ss.str(text_line);
		
		mass = 0.;
		failed = false;
		
		ss >> nb >> specimen.name >> specimen.mass >> specimen.charge;
		for (i = 0; i < NB_OF_CHEM_ELEMENTS; i++) 
		{
			ss >> specimen.formula[i];
			if (specimen.formula[i] > 0) 
			{
				// if the nb of carbon atoms in specimen is > max value, the specimen is not taken into account;
				if ((i == 2 && specimen.formula[i] > max_nb_carbon) || elements.find(i) == elements.end())
					failed = true;
				else mass += specimen.formula[i] *elements[i].mass;
			}
		}
		ss >> specimen.enthalpy;
		if (specimen.name == "H2" || specimen.name == "e-" || specimen.name == "He" || specimen.name == "N2" 
			|| specimen.name == "O2" || fabs(specimen.enthalpy) > DBL_EPSILON) {
				specimen.is_enth_def = true;
		}
		else specimen.is_enth_def = false;

		specimen.enthalpy *= 1.e+10/AVOGADRO_CONSTANT; // kJoule/mol conversion to erg per molecule;
		if (find_specimen(specimen.name) >= 0)
			failed = true;

		if (!failed) 
		{
			if (specimen.mass < 0.99) // it is electron
				specimen.type = specimen.name;
			else
			{ 
				if ( fabs(mass - specimen.mass) > 0.99 && verbosity) 
					cout << "    Please, check the formula of chemical specimen " << specimen.name << endl;
			
				if (specimen.charge == 0) 
					specimen.type = "neutral";
				else specimen.type = "ion";
			} 
			specimen.mass *= ATOMIC_MASS_UNIT; // mass in g;
			species.push_back(specimen);
		}
	}
	input.close();
	nb_of_species = (int) species.size();
}

void chem_network::init_gmantle_species(const std::string file_name)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;
	double b_energy;

	string name;
	stringstream ss, ss2;
	ifstream input;
	
	chem_specimen specimen;
	vector<chem_specimen> gmantle_species;
	
	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with binding energies " << file_name << endl;
		exit(1);
	}
	if (verbosity) {
		cout << "The data is read from file: " << file_name << endl;
	}

	ss2.clear();
	while (!input.eof())
	{
		// comment lines are read:
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#' || text_line[0] == '!');

		if (text_line[0] == '\0') // check for empty line at the file end;
			break;

		ss.str(text_line);
		ss >> name >> b_energy; // in the file, binding energy is in K;
		
		// search for the specimen in the list of gas-phase species:
		i = find_specimen(name);
		if (i >= 0)
		{
			// only neutral species, adsorbed on grains, are considered;
			if (species[i].charge == 0) 
			{
				specimen = species[i];
				specimen.type = "adsorbed";
				specimen.name = '*' + specimen.name; // name of the adsorbed specimen has the symbol '*';
			
				// specimen enthalpy must be in erg, but binding energy is in K;
				specimen.enthalpy -= b_energy*BOLTZMANN_CONSTANT;

				// there is no upper restriction on binding energy here:
				specimen.bind_en = b_energy; // in K;
				// binding energy in K, mass in a.m.u.:
				specimen.surf_freq = SURF_VIBR_FREQ*sqrt(b_energy*ATOMIC_MASS_UNIT/specimen.mass);
				specimen.diff_barrier = DIFFUSION_TO_BINDING_ENERGY_RATIO *b_energy; // must be in K;

#ifdef H_ATOM_DIFF_BARRIER
				if (specimen.name == "*H")
					specimen.diff_barrier = H_ATOM_DIFF_BARRIER;
#endif
				gmantle_species.push_back(specimen);
			}
		}
		else {
			ss2 << " " << name;
		}
	}
	nb_of_gmantle_species = (int) gmantle_species.size();
	nb_of_species += nb_of_gmantle_species;

	for (i = 0; i < nb_of_gmantle_species; i++) {
		species.push_back(gmantle_species[i]);			
	}
	if (verbosity) {
		cout << "  following grain mantle species are not identified: " << endl 
			<< "  " << ss2.str() << endl << endl;
	}
}

void chem_network::exclude_chem_specimen(string sname)
{
	vector<chem_specimen>::iterator it = species.begin();
	while (it != species.end())
	{
		if (it->name == sname) 
		{
			nb_of_species--;
			if (it->type == "adsorbed")
				nb_of_gmantle_species--;
				
			it = species.erase(it);
			if (verbosity) {
				cout << "Specimen was removed from the list: " << sname << endl;
			}
			break;
		}
		it++;
	}
}

void chem_network::init_species_nbs()
{
	int i;
	string n;
	vector<string> nlist;

	sort(species.begin(), species.end());

	// Note! All these species must be included in the specimen list:
	h2_nb = find_specimen("H2");
	h_nb = find_specimen("H");
	he_nb = find_specimen("He");
	e_nb = find_specimen("e-");
	co_nb = find_specimen("CO");
	h2o_nb = find_specimen("H2O");
	oh_nb = find_specimen("OH");
	nh3_nb = find_specimen("NH3");
	ch3oh_nb = find_specimen("CH3OH");
	oi_nb = find_specimen("O");
	ci_nb = find_specimen("C");
	cii_nb = find_specimen("C+");
	ah2_nb = find_specimen("*H2");
    hp_nb = find_specimen("H+");
    h3p_nb = find_specimen("H3+");

	if (h2_nb < 0 || ah2_nb < 0 || h_nb < 0 || he_nb < 0 || e_nb < 0 || h2o_nb < 0 || co_nb < 0 || oh_nb < 0 || nh3_nb < 0 || ch3oh_nb < 0
		|| oi_nb < 0 || ci_nb < 0 || cii_nb < 0 || hp_nb < 0 || h3p_nb < 0) {
		cout << "Error in " << SOURCE_NAME << ": can't find necessary chemical specimen in the list";
		exit(1);
	}

	for (i = 0; i < nb_of_species; i++) 
	{
		if (species[i].type == "neutral") {
			species[i].n_nb = i;
		}
		else if (species[i].type == "ion" )
		{
			n.assign(species[i].name, 0, species[i].name.length()-1); // removing the charge sign;
			species[i].n_nb = find_specimen(n);
		}
		else if (species[i].type == "e-") {
			species[i].n_nb = -1;
		}
		else if (species[i].type == "adsorbed") // it is assumed that adsorbed specimen is neutral;
		{
			n.assign(species[i].name, 1, species[i].name.length()-1); // removing the '*' sign;
			species[i].n_nb = find_specimen(n);
		}
		
		if (species[i].type == "ion" && species[i].n_nb < 0) {
			nlist.push_back(species[i].name);	
		}
	}

	if (verbosity > 0) {	
		if (nlist.size() > 0)
		{
			cout << "Can't find neutral for the ion specimen: ";
			for (i = 0; i < (int) nlist.size(); i++) {
				cout << nlist[i] << ", ";
			}
			cout << endl;
		}
	
		nlist.clear();
		for (i = 0; i < nb_of_species - nb_of_gmantle_species; i++) {
			if (species[i].charge == 0) 
			{
				n = '*' + species[i].name;
				if (find_specimen(n) < 0)
					nlist.push_back(species[i].name);
			}
		}
		if (nlist.size() > 0)
		{
			cout << "Warning: there is no adsorbed specimen for the neutral: ";
			for (i = 0; i < (int) nlist.size(); i++) {
				cout << nlist[i] << ", ";
			}
			cout << endl;
		}
	}
}

void chem_network::init_network_umistf(const std::string file_name, bool update)
{
	// max number of reaction fits:
	const int fnb_max = 3;
	// parameter values are specific to this format:
	const int max_nb_of_reactants = 2, max_nb_of_products = 4, nb_of_param = 3;

	bool failed;
	int i, j, l, nb, nb_of_fits, nb_of_reactants, nb_of_products, type, u_nb(1);
	int reactant[max_nb_of_reactants], product[max_nb_of_products];
	
	char ch, text_line[MAX_TEXT_LINE_WIDTH], str[MAX_TEXT_LINE_WIDTH];
	double b_en, parameters[nb_of_param*fnb_max], temp_min[fnb_max], temp_max[fnb_max];
	
	string name, rcode;
	stringstream ss; 
	ifstream input;

	vector<string> undef_species;
	chem_reaction reaction;
	
	sputtering_yield *sp_yield = 0;
	reaction_rate_data *rate_data = 0;

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with chemical reactions " << file_name << endl;
		exit(1);
	}
	if (verbosity) {
		cout << "Network " << file_name << endl;
	}

	l = 1;
	undef_species.clear();
	while (!input.eof())
	{
		name = "";
		type = -1;
		failed = false;
		
		memset(reactant, 0, max_nb_of_reactants*sizeof(int));
		memset(product, 0, max_nb_of_products*sizeof(int));
		memset(parameters, 0, fnb_max*nb_of_param*sizeof(double));

		// default values (special names "PHOTON", "GRAIN" and etc. are not reactants or products):
		nb_of_reactants = max_nb_of_reactants;
		nb_of_products = max_nb_of_products;
		
		// comment lines are read, comment lines can be located throught the file:
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		if (text_line[0] == '\0') // check for empty line at the file end;
			break;

		ss.clear();
		ss.str(text_line);
		// reading the number of reaction:
		ss >> nb >> ch;
		// reading the reaction type:
		ss.getline(str, MAX_TEXT_LINE_WIDTH, ':');
		rcode = str; // the reaction string codes are specific to the database;

		// reading reactants:
		for (i = 0; i < max_nb_of_reactants; i++)
		{
			ss.getline(str, MAX_TEXT_LINE_WIDTH, ':');
			if (strcmp(str, "PHOTON") != 0 && strcmp(str, "CRPHOT") != 0 && strcmp(str, "CRP") != 0 && strcmp(str, "GRAIN") != 0) 
			{
				reactant[i] = find_specimen(str);
				if (reactant[i] < 0) 
				{
					failed = true;
					for (j = 0; j < (int) undef_species.size(); j++) {
						if (undef_species[j] == str)
							break;
					}
					if (j == (int) undef_species.size())
						undef_species.push_back(str);
				}
			}
			else { // if one of the special words is present in the reactant list:
				nb_of_reactants = 1;
			}
		}
		// reading products:
		for (i = 0; i < max_nb_of_products; i++)
		{	
			ss.getline(str, MAX_TEXT_LINE_WIDTH, ':');
			if (strcmp(str, "") != 0 && strcmp(str, "PHOTON") != 0 && strcmp(str, "GRAIN") != 0) 
			{
				product[i] = find_specimen(str);
				if (product[i] < 0)
				{
					failed = true;
					for (j = 0; j < (int) undef_species.size(); j++) {
						if (undef_species[j] == str)
							break;
					}
					if (j == (int) undef_species.size())
						undef_species.push_back(str);
				}
			}
			else if (nb_of_products == max_nb_of_products) {
				nb_of_products = i;
			}
		}
		// Note: nb of fits must be less or equal to the maximal value fnb_max, defined above:
		ss >> nb_of_fits >> ch;
		for (i = 0; i < nb_of_fits; i++) 
		{
			for (j = 0; j < nb_of_param; j++) {
				ss >> parameters[i*nb_of_param + j] >> ch;
			}	
			ss >> temp_min[i] >> ch >> temp_max[i] >> ch;

			// reading the source type:
			ss.getline(str, MAX_TEXT_LINE_WIDTH, ':');
		
			// reading the accuracy type:
			ss.getline(str, MAX_TEXT_LINE_WIDTH, ':');

			// reading the reference:
			ss >> ch; // reading symbol '"' or ":"
			if (ch == '"') {
				ss.getline(str, MAX_TEXT_LINE_WIDTH, '"');
				ss >> ch;
			}
			ss >> ch;
			if (ch == '"') {
				ss.getline(str, MAX_TEXT_LINE_WIDTH, '"');
				ss >> ch;
			}
			if (ss.eof() || ss.fail() || ss.bad()) {
				cout << SOURCE_NAME << ": error ocurred in string stream in file " << file_name << endl
					<< ss.str();
				exit(1);
			}
		}
		
		if (!failed)
		{
			check_order(species, reactant, nb_of_reactants);
			check_order(species, product, nb_of_products);
						
			// only reactions with defined type are taken into account;
			type = react_type_determ(rcode, species, reactant, nb_of_reactants, product, nb_of_products, verbosity);
			
			// the reaction name is determined by the identificator rcode;
			name = get_reaction_name(rcode, species, reactant, nb_of_reactants, product, nb_of_products);

			if (type == -1 && verbosity) {
				cout << "Can't determine type of the reaction: " << name << endl;
			}
			else
			{
				reaction.delete_data();

				reaction.rcode = rcode;
				reaction.type = type;
				reaction.name = name;
				reaction.nb_of_reactants = nb_of_reactants;
				reaction.nb_of_products = nb_of_products;
				
				reaction.nb_of_fits = nb_of_fits;
				reaction.nb_of_param = nb_of_param;

				reaction.reactant = new int [nb_of_reactants];
				memcpy(reaction.reactant, reactant, nb_of_reactants*sizeof(int));

				reaction.product = new int [nb_of_products];
				memcpy(reaction.product, product, nb_of_products*sizeof(int));

				reaction.temp_min = new double [nb_of_fits];
				reaction.temp_max = new double [nb_of_fits];
				reaction.parameters = new double [nb_of_fits*nb_of_param];

				for (i = 0; i < nb_of_fits; i++) 
				{
					reaction.temp_min[i] = temp_min[i]; // temperature must be in K;
					reaction.temp_max[i] = temp_max[i];
				
					// redefinition of rate constants: temperature is measured in K;
					// for UV induced reactions the parameter[1] > 0 only for CO dissociation;
					if ((reaction.type >= 4 && reaction.type <= 8) || (reaction.type >= 14 && reaction.type <= 25) 
						|| reaction.type == 40 || reaction.type == 41) {
						parameters[i*nb_of_param] = parameters[i*nb_of_param] /pow(300., parameters[i*nb_of_param+1]);
					}
					for (j = 0; j < nb_of_param; j++) {
						reaction.parameters[i*nb_of_param+j] = parameters[i*nb_of_param+j];
					}
				}
				
				// Rates at minimal and maximal temperatures:
				i = nb_of_fits-1; 
				if ((reaction.type >= 4 && reaction.type <= 8) || reaction.type == 40 || reaction.type == 41) 
				{
					reaction.min_rate = reaction.parameters[0] *pow(reaction.temp_min[0], reaction.parameters[1]) *reaction.parameters[2];
					reaction.max_rate = reaction.parameters[i*nb_of_param] *pow(reaction.temp_max[i], reaction.parameters[i*nb_of_param+1])
						*reaction.parameters[i*nb_of_param+2];
				}
				else if (reaction.type >= 14 && reaction.type <= 25) 
				{
					// Redifinition of minimal temperature (if the reaction rate decreases with decreasing temperature):
					if (reaction.temp_min[0] > 10. && reaction.parameters[2] > -DBL_EPSILON 
						&& (reaction.parameters[1] + reaction.parameters[2]/reaction.temp_min[0] > 0.)) 
					{
					/*	if (verbosity) {
							cout << left << setw(4) << l << "Min. temp. is changed to 10 K: " << setw(25) << reaction.name << endl
								<< "   Tmin= " << setw(13) << reaction.temp_min[0] << " b= " << setw(13) << reaction.parameters[1] 
								<< " g= " << reaction.parameters[2] << endl;
						}*/
						reaction.temp_min[0] = 10.;
						l++;
					}
					reaction.min_rate = reaction.parameters[0] *pow(reaction.temp_min[0], reaction.parameters[1]) 
						*exp(-reaction.parameters[2]/reaction.temp_min[0]);

					reaction.max_rate = reaction.parameters[i*nb_of_param] *pow(reaction.temp_max[i], reaction.parameters[i*nb_of_param+1]) 
						*exp(-reaction.parameters[i*nb_of_param+2]/reaction.temp_max[i]);
				}
				else {
					// undefined minimal and maximal rates:
					reaction.min_rate = 0.;
					reaction.max_rate = 1.+99;
				}

				// it is assumed that the neutral is the first reactant in the list for neutral-ion reactions;
				reaction.mass1 = species[ reaction.reactant[0] ].mass;
				if (nb_of_reactants == 2) 
				{
					reaction.mass2 = species[ reaction.reactant[1] ].mass;
					reaction.mass_sum = reaction.mass1 + reaction.mass2;
					reaction.reduced_mass = reaction.mass1 *reaction.mass2 /reaction.mass_sum;
				}

				// the energy that is released in reaction, E_reactants - E_products, enthalpy must be in erg:
				failed = false;
				for (i = 0; i < nb_of_reactants; i++) {
					reaction.energy_released += species[ reactant[i] ].enthalpy;
					if (!species[ reactant[i] ].is_enth_def)
						failed = true;
				}
				for (i = 0; i < nb_of_products; i++) {
					reaction.energy_released -= species[ product[i] ].enthalpy;
					if (!species[ product[i] ].is_enth_def)
						failed = true;
				}
				if (failed)
					reaction.energy_released = 0.;

				// neutral-ion reactions:
				if (reaction.type == 20 || reaction.type == 21 || reaction.type == 22)
				{
					// Reactions with temperature dependence, effective temperature method is employed:
					if (fabs(reaction.parameters[2]) > DBL_EPSILON || fabs(reaction.parameters[1]) > DBL_EPSILON 
						|| nb_of_fits > 1) {
						reaction.cs_reconstr = 1;
					}
					// reactions with constant rate (exothermic reactions):
					else reaction.cs_reconstr = 0;
				}
				// adsorption:
				else if (reaction.type == 27) {
					// rate coefficient is calculated, mass in g:
					reaction.parameters[0] = sqrt( 8.*BOLTZMANN_CONSTANT/(M_PI*species[reactant[0]].mass) ); 
				}
				// thermal desorption:
				else if (reaction.type == 28) 
				{
					b_en = species[reactant[0]].bind_en;
#ifdef SET_MAX_THERMAL_EVAPOR_BARRIER
					if (b_en > SET_MAX_THERMAL_EVAPOR_BARRIER)
						b_en = SET_MAX_THERMAL_EVAPOR_BARRIER;
#endif
					// vibrational frequency of the species in the surface potential well (Hasegawa & Herbst, ApJSS 82, p. 167-195, 1992):
					reaction.parameters[0] = SURF_VIBR_FREQ*sqrt(b_en *ATOMIC_MASS_UNIT/species[reactant[0]].mass);
					reaction.parameters[1] = b_en; // modified binding enrgy;
				}
				// cosmic-ray induced desorption:
				else if (reaction.type == 29) 
				{
					b_en = species[reactant[0]].bind_en;
#ifdef SET_MAX_THERMAL_EVAPOR_BARRIER
					if (b_en > SET_MAX_THERMAL_EVAPOR_BARRIER)
						b_en = SET_MAX_THERMAL_EVAPOR_BARRIER;
#endif
					// Hasegawa, Herbst, MNRAS 261, 83 (1993); see also Bringa, Johnson, ApJ 603, 159 (2004); the constant value by Hasegawa, Herbst (1993) is used;
					reaction.parameters[0] = SURF_VIBR_FREQ*sqrt(b_en*ATOMIC_MASS_UNIT/species[reactant[0]].mass)
						*exp(-b_en/70.) *3.2e-19;
				}
				// sputtering:
				else if (reaction.type == 34 || reaction.type == 35)
				{
					if (verbosity)
						cout << left << "Evaluating sputtering yield: " << reaction.name << endl;
					
					b_en = species[reactant[0]].bind_en;
#ifdef SET_MAX_THERMAL_EVAPOR_BARRIER
					if (b_en > SET_MAX_THERMAL_EVAPOR_BARRIER)
						b_en = SET_MAX_THERMAL_EVAPOR_BARRIER;
#endif
					// projectile has the first place:
					sp_yield = new sputtering_yield(species[reactant[1]].mass, species[reactant[0]].mass, b_en);
					
					rate_data = new reaction_rate_data();
					rate_data->calc_data(sp_yield);
					
					rate_data_array.push_back(*rate_data);
					reaction.rate_data = rate_data;
					
					delete sp_yield;
				}
				// ion attachment on grains
				else if (reaction.type == 44)
					nb_reactions_ion_grains++;
				
				if (update) {
					for (i = 0; i < (int) reaction_array.size(); i++) {
						if (reaction_array[i] == reaction) {
							if (verbosity) {
								cout << u_nb++ << " reaction data were updated for " << reaction_array[i].name << endl;
							}
							reaction_array.erase(reaction_array.begin() + i);
							i--;
						}
					}
				}
				reaction_array.push_back(reaction);	
			}
		}
	}
	nb_of_reactions = (int) reaction_array.size();

	if (verbosity) {
		cout << "    undefined species:" << endl;
		for (i = 0; i < (int) undef_species.size(); i++) {
			cout << "  " << undef_species[i];
		}
		cout << endl;
	}
}

void chem_network::init_h2_formation()
{	
	chem_reaction reaction;

	reaction.nb_of_fits = 1;
	reaction.nb_of_param = 3;
	reaction.rcode = "HFG";

	reaction.nb_of_reactants = 2;
	reaction.nb_of_products = 1;

	reaction.reactant = new int [reaction.nb_of_reactants];
	reaction.product = new int [reaction.nb_of_products];

	reaction.reactant[0] = reaction.reactant[1] = h_nb;
	reaction.product[0] = h2_nb;

	reaction.type = react_type_determ(reaction.rcode, species, reaction.reactant, reaction.nb_of_reactants, 
		reaction.product, reaction.nb_of_products, verbosity);
	
	reaction.name = get_reaction_name(reaction.rcode, species, reaction.reactant, reaction.nb_of_reactants, 
		reaction.product, reaction.nb_of_products);
					
	reaction.temp_min = new double [reaction.nb_of_fits];
	reaction.temp_max = new double [reaction.nb_of_fits];

	reaction.temp_min[0] = 3.;
	reaction.temp_max[0] = 44000.;

	reaction.mass1 = reaction.mass2 = ATOMIC_MASS_UNIT;
	reaction.mass_sum = 2.*ATOMIC_MASS_UNIT;
	reaction.reduced_mass = 0.5*ATOMIC_MASS_UNIT;

	reaction.energy_released = 2.*species[h_nb].enthalpy - species[h2_nb].enthalpy;

	reaction.parameters = new double [reaction.nb_of_param];
	memset(reaction.parameters, 0, reaction.nb_of_param*sizeof(double));

	// if H2_FORMATION_MODE = 0, H2 formation is modelled, reaction *H + *H -> *H2 is taken into account;
	if (H2_FORMATION_MODE == 1) {
		// here, it is suggested that only half of the accreted atoms form H2 molecule,
		// another factor 0.5 due to the fact that two H atoms give one H2:
		reaction.parameters[0] = 0.25*sqrt(8.*BOLTZMANN_CONSTANT/(M_PI*ATOMIC_MASS_UNIT));
	}
	else if (H2_FORMATION_MODE > 1) {
		reaction.parameters[0] = STANDARD_H2_FORMATION_RATE;
	}
	reaction_array.push_back(reaction);
	nb_of_reactions = (int) reaction_array.size();
}

void chem_network::init_photoreact_surface_chemistry()
{
	const int nb_of_fits = 1;
	bool failed;
	int i, j, k, n, l, m, nb_of_param;
	double a;
	
	chem_reaction reaction;
	vector<chem_reaction> new_reactions;
	
	if (verbosity) {
		cout << "Photo-reactions on grain surface (based on gas-phase chemistry)... " << endl;
	}

	nb_of_param = 3;
	reaction.nb_of_param = nb_of_param;
	reaction.parameters = new double [nb_of_param];

	reaction.nb_of_fits = nb_of_fits;
	reaction.temp_min = new double [nb_of_fits];
	reaction.temp_max = new double [nb_of_fits];

	// Gas phase ionization reactions are considered,
	// it is assumed that number of fits for these reactions = 1,
	for (i = 0; i < nb_of_reactions; i++)
	{
		// "A + CRPhoton -> Ion + e-"
		// "A + CRPhoton -> B + Ion + e-"
		// "A + ISPhoton -> Ion + e-"
		// "A + ISPhoton -> B + Ion + e-"
		if (reaction_array[i].type == 4 || reaction_array[i].type == 5 
			|| reaction_array[i].type == 9 || reaction_array[i].type == 10)
		{	
			if (reaction_array[i].type == 4 || reaction_array[i].type == 5) 
				reaction.rcode = "CRDGS1";
			else reaction.rcode = "ISDGS1";

			// Search for dissociative recombination reactions:
			n = 0;
			l = reaction_array[i].nb_of_products - 2;
			a = 0.;
			for (m = 0; m < nb_of_reactions; m++)
			{
				if (reaction_array[m].type == 23 || reaction_array[m].type == 24) {
					if (reaction_array[i].product[l] == reaction_array[m].reactant[0] &&
						reaction_array[i].product[l+1] == reaction_array[m].reactant[1]) 
					{
						failed = false;			
						delete [] reaction.product;
						delete [] reaction.reactant;
						
						reaction.nb_of_products = reaction_array[i].nb_of_products + reaction_array[m].nb_of_products - 2;
						reaction.nb_of_reactants = reaction_array[i].nb_of_reactants; // = 1
						
						reaction.product = new int [reaction.nb_of_products];
						reaction.reactant = new int [reaction.nb_of_reactants]; // nb of reactants is equal to 1;
						
						if (reaction.nb_of_products == 1)
							failed = true;

						k = find_specimen( '*' + species[ reaction_array[i].reactant[0] ].name );
						reaction.reactant[0] = k;
				
						if (k == -1)
							failed = true;

						for (j = 0; j < reaction_array[i].nb_of_products - 2; j++)
						{
							k = find_specimen( '*' + species[ reaction_array[i].product[j] ].name );
							reaction.product[j] = k;
							
							if (k == -1)
								failed = true;
						}
						for (j = 0; j < reaction_array[m].nb_of_products; j++)
						{
							k = find_specimen( '*' + species[ reaction_array[m].product[j] ].name );
							reaction.product[reaction_array[i].nb_of_products - 2 + j] = k;
							
							if (k == -1)
								failed = true;
						}

						if (!failed)
						{
							// only reactions with defined type are taken into account;
							reaction.type = react_type_determ(reaction.rcode, species, 
								reaction.reactant, reaction.nb_of_reactants, reaction.product, reaction.nb_of_products, verbosity);
			
							if (reaction.type >= 0)
							{	
								reaction.name = get_reaction_name(reaction.rcode, species, 
									reaction.reactant, reaction.nb_of_reactants, reaction.product, reaction.nb_of_products);
								
								// only photodissociation reaction "CO + CRPhoton -> O + C" has temperature dependence;
								reaction.temp_min[0] = reaction_array[i].temp_min[0];
								reaction.temp_max[0] = reaction_array[i].temp_max[0];

								reaction.min_rate = reaction_array[i].min_rate;
								reaction.max_rate = reaction_array[i].max_rate;
								
								reaction.mass1 = reaction_array[i].mass1;
								memcpy(reaction.parameters, reaction_array[i].parameters, nb_of_param*sizeof(double));
								
								// for multiple channels:
								n++;
								a += reaction_array[m].parameters[0];
								reaction.parameters[0] *= reaction_array[m].parameters[0];

								new_reactions.push_back(reaction);
							}
						}
					}
				}
			}
			for (m = 0; m < n; m++) {
				new_reactions[(int) new_reactions.size() - m - 1].parameters[0] /= a;
			}
		}
	}

	// Dissociation reactions are considered:
	for (i = 0; i < nb_of_reactions; i++)
	{
		// "A + CRPhoton -> B + C + D + E"
		// "A + ISPhoton -> B + C + D + E"
		if (reaction_array[i].type == 6 || reaction_array[i].type == 11)
		{
			failed = false;
			if (reaction_array[i].type == 6) 
				reaction.rcode = "CRDGS2";
			else reaction.rcode = "ISDGS2";
			
			delete [] reaction.reactant;
			delete [] reaction.product;
			
			reaction.nb_of_products = reaction_array[i].nb_of_products;
			reaction.nb_of_reactants = reaction_array[i].nb_of_reactants; // = 1

			reaction.reactant = new int [reaction.nb_of_reactants];
			reaction.product = new int [reaction.nb_of_products];
			
			k = find_specimen( '*' + species[ reaction_array[i].reactant[0] ].name );
			reaction.reactant[0] = k;
				
			if (k == -1)
				failed = true;

			for (j = 0; j < reaction.nb_of_products; j++)
			{
				k = find_specimen( '*' + species[ reaction_array[i].product[j] ].name );
				reaction.product[j] = k;
				
				if (k == -1)
					failed = true;
			}

			if (!failed)
			{				
				// only reactions with defined type are taken into account;
				reaction.type = react_type_determ(reaction.rcode, species, 
					reaction.reactant, reaction.nb_of_reactants, reaction.product, reaction.nb_of_products, verbosity);
				
				if (reaction.type >= 0)
				{
					reaction.name = get_reaction_name(reaction.rcode, species, 
						reaction.reactant, reaction.nb_of_reactants, reaction.product, reaction.nb_of_products);

					// only photodissociation reaction "CO + CRPhoton -> O + C" has temperature dependence;
					reaction.temp_min[0] = reaction_array[i].temp_min[0];
					reaction.temp_max[0] = reaction_array[i].temp_max[0]; // nb of fits is equal to 1

					reaction.min_rate = reaction_array[i].min_rate;
					reaction.max_rate = reaction_array[i].max_rate;

					reaction.mass1 = reaction_array[i].mass1;

					memcpy(reaction.parameters, reaction_array[i].parameters, nb_of_param*sizeof(double));
					new_reactions.push_back(reaction);
				}
			}
		}
	}

	// the energy that is released in reaction:
	for (i = 0; i < (int) new_reactions.size(); i++)
	{
		failed = false;
		new_reactions[i].energy_released = 0.;
		for (j = 0; j < new_reactions[i].nb_of_reactants; j++) 
		{
			new_reactions[i].energy_released += species[ new_reactions[i].reactant[j] ].enthalpy;
			if (!species[ new_reactions[i].reactant[j] ].is_enth_def)
				failed = true;
		}
		for (j = 0; j < new_reactions[i].nb_of_products; j++) {
			new_reactions[i].energy_released -= species[ new_reactions[i].product[j] ].enthalpy;
			if (!species[ new_reactions[i].product[j] ].is_enth_def)
				failed = true;
		}
		if (failed)
			new_reactions[i].energy_released = 0.;
	}
	
	for (i = 0; i < (int) new_reactions.size(); i++) { 
		// scaling factor is applied for all photodissociation reactions:
		new_reactions[i].parameters[0] *= ICE_PHOTODISS_EFFICIENCY;
		reaction_array.push_back(new_reactions[i]);
	}
	nb_of_reactions = (int) reaction_array.size();
}

void chem_network::init_grain_surface_chemistry(const string fname)
{
	// maximal nb of reactions with the same reactants but different products, and with activation barriers is 3;
	const int nb_of_fits = 1, nb_of_param = 10;
	const int max_nb_reactants = 2, max_nb_products = 4, field_length = 15; // the length of the text field alotted to specimen name;
	char text_line[MAX_TEXT_LINE_WIDTH];
	
	bool failed;
	int i, j, k, l, m, n, nb_of_reactants, nb_of_products;
	double chem_barrier_factor, diff_barrier_factor, barrier;
	int *reactant, *product;
	
	string str;
	ifstream input;

	chem_reaction reaction;
	vector<double> act_barr_arr;
	vector<chem_reaction> new_reactions;
	vector<string> undef_species;

	if (verbosity) {
		cout << "Network " << fname << endl;
	}
		
	// Formula are given by, e.g., Reboussin et al., MNRAS 440, p. 3557 (2014);
	chem_barrier_factor = 4.*CHEM_BARRIER_THICKNESS *M_PI/PLANCK_CONSTANT;
	diff_barrier_factor = 4.*DIFF_BARRIER_THICKNESS *M_PI/PLANCK_CONSTANT;

	reactant = new int [max_nb_reactants];
	product = new int [max_nb_products];

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with chemical reactions " << fname << endl;
		exit(1);
	}
	
	reaction.nb_of_param = nb_of_param;
	reaction.nb_of_fits = nb_of_fits;
	reaction.rcode = "GS";
	reaction.parameters = new double [nb_of_param];

	reaction.temp_min = new double [nb_of_fits];
	reaction.temp_max = new double [nb_of_fits];
	reaction.temp_min[0] = 3.;
	reaction.temp_max[0] = 44000.;
	
	while (!input.eof())
	{	
		failed = false;
		nb_of_reactants = max_nb_reactants;
		nb_of_products = max_nb_products;

		// comment lines are read:
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#' || text_line[0] == '!');
		
		if (text_line[0] == '\0') // empty line at the file end;
			break;
	
		for (j = 0; j < max_nb_reactants; j++)
		{
			i = k = j*field_length;
			str = "";
			while (text_line[i] != ' ' && i - k < field_length) {
				str.push_back(text_line[i]);
				i++;
			}
			
			if (i > k) 
			{
				str[0] = '*';
				reactant[j] = find_specimen(str);
				
				if (reactant[j] < 0) 
				{
					failed = true;
					for (k = 0; k < (int) undef_species.size(); k++) {
						if (undef_species[k] == str)
							break;
					}
					if (k == (int) undef_species.size())
						undef_species.push_back(str);
				}
			}
			else if (nb_of_reactants == max_nb_reactants) {
				nb_of_reactants = j;
			}
		}
		
		for (j = 0; j < max_nb_products; j++) 
		{
			i = k = (max_nb_reactants + j)*field_length;
			str = "";
			while (text_line[i] != ' ') {
				str.push_back(text_line[i]);
				i++;
			}
			
			if (i > k)
			{
				if (str[0] == 'J')
					str[0] = '*';
				product[j] = find_specimen(str);
				
				if (product[j] < 0) 
				{
					failed = true;
					for (k = 0; k < (int) undef_species.size(); k++) {
						if (undef_species[k] == str)
							break;
					}
					if (k == (int) undef_species.size())
						undef_species.push_back(str);
				}
			}
			else if (nb_of_products == max_nb_products) {
				nb_of_products = j;
				break;
			}
		}
		
		k = (max_nb_reactants + max_nb_products)*field_length + 22;
		str = "";
		for (i = k; i < k + 10; i++) {
			str.push_back(text_line[i]);
		}
		barrier = atof(str.c_str());
		
		if (!failed)
		{	
			reaction.nb_of_reactants = nb_of_reactants;
			reaction.nb_of_products = nb_of_products;

			delete [] reaction.reactant;
			delete [] reaction.product;
		
			reaction.reactant = new int [nb_of_reactants];
			memcpy(reaction.reactant, reactant, nb_of_reactants*sizeof(int));

			reaction.product = new int [nb_of_products];
			memcpy(reaction.product, product, nb_of_products*sizeof(int));

			reaction.type = react_type_determ(reaction.rcode, species, reaction.reactant, nb_of_reactants, 
				reaction.product, nb_of_products, verbosity);
			
			if (reaction.type >= 0)
			{
				reaction.name = get_reaction_name(reaction.rcode, species, reaction.reactant, nb_of_reactants, 
					reaction.product, nb_of_products);

				reaction.mass1 = species[ reaction.reactant[0] ].mass;
				reaction.mass2 = species[ reaction.reactant[1] ].mass;
			
				reaction.mass_sum = reaction.mass1 + reaction.mass2;
				reaction.reduced_mass = reaction.mass1 *reaction.mass2 /reaction.mass_sum;
			
				memset(reaction.parameters, 0, nb_of_param*sizeof(double));
				j = reaction.reactant[0];
				k = reaction.reactant[1];

				// Be careful with treatments of homogeneous reactions (involving the same species, H + H -> H2), 
				// their rates by statistical arguments are only half of those for heterogeneous reactions (Semenov et al., A&A 522, A42, 2010);
				if (j == k)
					reaction.parameters[0] = 0.5;
				else reaction.parameters[0] = 1.;
			
				reaction.parameters[1] = barrier; // in K;
				reaction.parameters[2] = 
					diff_barrier_factor*sqrt(2.*species[j].mass *BOLTZMANN_CONSTANT *species[j].diff_barrier);
				
				reaction.parameters[3] = 
					diff_barrier_factor*sqrt(2.*species[k].mass *BOLTZMANN_CONSTANT *species[k].diff_barrier);
				
				reaction.cs_reconstr = 0;
				if ((species[j].name == "*H" || species[j].name == "*H2") && H_H2_QUANTUM_TUNNELING_DIFF_ON)
					reaction.cs_reconstr = 1;
				
				if ((species[k].name == "*H" || species[k].name == "*H2") && H_H2_QUANTUM_TUNNELING_DIFF_ON)
					reaction.cs_reconstr += 2;
						
				if (reaction.parameters[1] > DBL_EPSILON)  // reactions with activation barrier:
				{
					// the larger value of the characteristic frequencies of the two reactants:
					reaction.parameters[4] = (species[j].surf_freq > species[k].surf_freq) ? species[j].surf_freq : species[k].surf_freq;
					// quantum mechanical probability for tunnelling:
					// Note the comment by Cuppen et al. (Space Sci. Rev. 212, p.1, 2017) about tunneling,
					// in reaction CH4 + C2H -> C2H2 + CH3, tunneling species is H atom, and it's mass must be taken into account;
					reaction.parameters[5] = chem_barrier_factor *sqrt(2.*reaction.reduced_mass *BOLTZMANN_CONSTANT *reaction.parameters[1]);
				}
				
				if (reaction.type == 39) 
				{ // desorption encounter mechanism for H2 molecule, Hincelin et al., A&A 2014, arXiv:1410.7375v2;
					if (reaction.name == "*H2 + *H2 -> *H2 + H2") {
						reaction.parameters[4] = H2_ON_H2_BINDING_ENERGY;
						reaction.parameters[5] = H2_ON_H2_BINDING_ENERGY *DIFFUSION_TO_BINDING_ENERGY_RATIO;
						reaction.parameters[6] = 
							diff_barrier_factor*sqrt(2.*species[k].mass *BOLTZMANN_CONSTANT *reaction.parameters[5]);
					}
				}
				else 
				{ // initialisation of chemical desorption (Garrod et al. 2007):
					if (reaction.nb_of_products == 1) 
					{
						if (species[ reaction.product[0] ].type == "neutral")
							reaction.parameters[0] *= CHEM_DESORPTION_FACTOR;
						else reaction.parameters[0] *= 1. - CHEM_DESORPTION_FACTOR;
					}
					else 
					{
						if (species[ reaction.product[0] ].type == "neutral")
							reaction.parameters[0] = 0.; // reactive desorption for reactions with multiple products is not allowed;
					}
				}
				
                // perhaps, chemical desorption factors must be redefined here:
				if (reaction.name == "*H + *H -> *H2" || reaction.name == "*H + *H -> H2") {
					if (H2_FORMATION_MODE == 0)
						new_reactions.push_back(reaction);
				}
				else if (reaction.parameters[0] > DBL_EPSILON) // reactions with zero rate are removed;
					new_reactions.push_back(reaction);
			}
		}
	}

	// the energy that is released in reaction:
	for (i = 0; i < (int) new_reactions.size(); i++) 
	{
		failed = false;
		new_reactions[i].energy_released = 0.;
		for (j = 0; j < new_reactions[i].nb_of_reactants; j++) 
		{
			new_reactions[i].energy_released += species[ new_reactions[i].reactant[j] ].enthalpy;
			if (!species[ new_reactions[i].reactant[j] ].is_enth_def)
				failed = true;
		}
		for (j = 0; j < new_reactions[i].nb_of_products; j++) {
			new_reactions[i].energy_released -= species[ new_reactions[i].product[j] ].enthalpy;
			if (!species[ new_reactions[i].product[j] ].is_enth_def)
				failed = true;
		}
		if (failed)
			new_reactions[i].energy_released = 0.;
	}

	// Calculation of branching ratio:
	// It is postulated, that all reactions with the same reactants have activation barriers or all do not have,
	for (i = 0; i < (int) new_reactions.size(); i++) {
		if (species[ new_reactions[i].product[0] ].type == "adsorbed") 
		{
			m = 1;
			failed = false;
			act_barr_arr.clear();
		
			j = new_reactions[i].reactant[0];
			k = new_reactions[i].reactant[1];
			
			if (new_reactions[i].parameters[1] > DBL_EPSILON) {
				failed = true;
				act_barr_arr.push_back(new_reactions[i].parameters[1]);
			}
			// looking for reactions with the same reactants:
			for (l = 0; l < (int) new_reactions.size(); l++) {
				if (species[ new_reactions[l].product[0] ].type == "adsorbed") 
				{
					if (((new_reactions[l].reactant[0] == j && new_reactions[l].reactant[1] == k) 
						|| (new_reactions[l].reactant[0] == k && new_reactions[l].reactant[1] == j)) && l != i) 
					{
						m++;
						if (new_reactions[l].parameters[1] > DBL_EPSILON) {
							failed = true;
							act_barr_arr.push_back(new_reactions[l].parameters[1]);
						}
					}
				}
			}
			
			if (m > 1) 
			{
				// for reactions without activation barrier;
				if (!failed) {
					new_reactions[i].parameters[0] /= m;
						
					// chemical desorption is applied only to reactions with one product;
					// the reaction *A + *B -> C follows reaction *A + *B -> *C, there is only one product,
					if (new_reactions[i].nb_of_products == 1)
						new_reactions[i+1].parameters[0] /= m;

					if (verbosity) {
						cout << left << setw(5) << i << setw(30) << new_reactions[i].name << " " 
							<< "branching ratio " << 1./m << endl;
					}
				}
				// for reactions with activation barrier;
				// activations barriers of all reactions in the group are saved to each reaction data in the group:
				else { 	
					n = 6;
					for (l = 1; l < (int) act_barr_arr.size() && n < nb_of_param; l++, n += 2) 
					{
						new_reactions[i].parameters[n] = act_barr_arr[l];
						new_reactions[i].parameters[n+1] = chem_barrier_factor 
							*sqrt(2.*new_reactions[i].reduced_mass *BOLTZMANN_CONSTANT *act_barr_arr[l]);
					}

					if (new_reactions[i].nb_of_products == 1) 
						memcpy(new_reactions[i+1].parameters+1, new_reactions[i].parameters+1, (nb_of_param-1)*sizeof(double));

					if (verbosity) {
						cout << left << setw(5) << i << setw(30) << new_reactions[i].name << " activ. barrier ";
						for (l = 0; l < m; l++) {
							cout << left << setw(12) << act_barr_arr[l];
						}
						cout << endl;
					}
				}
			}
		}
	}

	if (verbosity) {
		cout << "    undefined species:" << endl;
		for (i = 0; i < (int) undef_species.size(); i++) {
			cout << undef_species[i] << endl;
		}
		cout << endl;
	}

	for (i = 0; i < (int) new_reactions.size(); i++) {
		reaction_array.push_back(new_reactions[i]);
	}
	nb_of_reactions = (int) reaction_array.size();

	delete [] reactant;
	delete [] product;
}

void chem_network::check_reactions()
{
	bool is_failed;
	int i, j, k, l;

	chem_reaction reaction;
	vector<chem_reaction>::iterator it = reaction_array.begin();
	vector<string> undef_species;

	// reactions of ion neutralization on grains must be at the end of the list;
	for (i = 0; i < nb_of_reactions; i++) {
		if (it->rcode == "IAG") 
		{
			reaction = *it;
			it = reaction_array.erase(it);
			reaction_array.push_back(reaction);
		}
		else it++;
	}
	
	// derivation of reaction numbers must be here:
	h2_h_diss_nb = find_reaction("H + H2 -> H + H + H");
    h2_h2_diss_nb = find_reaction("H2 + H2 -> H2 + H + H");
    h2_e_diss_nb = find_reaction("H2 + e- -> H + H + e-");

	if (h2_h_diss_nb == -1 || h2_h2_diss_nb == -1 || h2_e_diss_nb == -1) {
		cout << "Error in " << SOURCE_NAME << ": can not find reaction..." << endl;
		is_failed = true;
	}

	// check that all reactions have defined type:
	is_failed = false;
	for (i = 0; i < nb_of_reactions; i++) {
		if (reaction_array[i].type < 0) 
		{
			is_failed = true;
			cout << "Error: reaction type is undefined " << reaction_array[i].name << endl;
		}
	}

	// check for charge conservation for all reactions:
	for (i = 0; i < nb_of_reactions; i++) 
	{
		k = 0;
		for (j = 0; j < reaction_array[i].nb_of_reactants; j++) {
			k += species[ reaction_array[i].reactant[j] ].charge;
		}

		for (j = 0; j < reaction_array[i].nb_of_products; j++) {
			k -= species[ reaction_array[i].product[j] ].charge;
		}
		if ((k != 0 && i < nb_of_reactions - nb_reactions_ion_grains) || (k != 1 && i >= nb_of_reactions - nb_reactions_ion_grains)) {
			cout << "Error: no charge conservation in the reaction " << reaction_array[i].name << endl;
			is_failed = true;
		}
	}

	// check for element conservation in each reaction:
	for (i = 0; i < nb_of_reactions; i++) 
	{
		for (l = 0; l < NB_OF_CHEM_ELEMENTS; l++)
		{
			k = 0;
			for (j = 0; j < reaction_array[i].nb_of_reactants; j++) {
				k += species[ reaction_array[i].reactant[j] ].formula[l];
			}

			for (j = 0; j < reaction_array[i].nb_of_products; j++) {
				k -= species[ reaction_array[i].product[j] ].formula[l];
			}
			if (k != 0) {
				cout << "Error: no element " << elements[l].name << " conservation in the reaction " << reaction_array[i].name << endl;
				is_failed = true;
			}
		}
	}

	// check that all species have production AND destruction reactions:
	for (i = 0; i < nb_of_species; i++)
	{
		k = 0;
		for (l = 0; l < nb_of_reactions; l++)
		{
			for (j = 0; j < reaction_array[l].nb_of_reactants; j++) {
				if (reaction_array[l].reactant[j] == i)
					k = 1;
			}
			if (k > 0)
				break;
		}
		if (k == 0) {
			cout << "Error: specimen has not destruction reaction: " << species[i].name << endl;
			is_failed = true;
		}
	}

	for (i = 0; i < nb_of_species; i++)
	{
		k = 0;
		for (l = 0; l < nb_of_reactions; l++)
		{
			for (j = 0; j < reaction_array[l].nb_of_products; j++) {
				if (reaction_array[l].product[j] == i)
					k = 1;
			}
			if (k > 0)
				break;
		}
		if (k == 0) {
			cout << "Error: specimen has not production reaction: " << species[i].name << endl;
			is_failed = true;
		}
	}

	// Check if each adsorbed specimen has destruction reaction other than desorption
	if (GRAIN_SURFACE_CHEMISTRY_ON > 0) {
		for (i = nb_of_species - nb_of_gmantle_species; i < nb_of_species; i++)
		{
			k = 0;
			for (l = 0; l < nb_of_reactions; l++)
			{
				if (reaction_array[l].rcode != "THD" && reaction_array[l].rcode != "CPD" && reaction_array[l].rcode != "CRD"
					&& reaction_array[l].rcode != "PHD" && reaction_array[l].rcode != "SPU")
				{
					for (j = 0; j < reaction_array[l].nb_of_reactants; j++) {
						if (reaction_array[l].reactant[j] == i)
							k = 1;
					}
				}
				if (k > 0)
					break;
			}
			if (k == 0) {
				cout << "Warning: adsorbed specimen has not destruction other than desorption: " << species[i].name << endl;
			}
		}
	}

/*	j = 1;
	cout << "Warning: following reactions have high minimal temperature: " << endl;
	for (i = 0; i < nb_of_reactions; i++) 
	{
		if (reaction_array[i].type >= 14 && reaction_array[i].type <= 25) {
			if (reaction_array[i].temp_min[0] > 10. && (fabs(reaction_array[i].parameters[1]) > DBL_EPSILON || fabs(reaction_array[i].parameters[2]) > DBL_EPSILON)) 
			{
				cout << left << setw(4) << j << setw(25) << reaction_array[i].name << " Tmin= " << setw(12) << reaction_array[i].temp_min[0] 
					<< " b= " << setw(12) << reaction_array[i].parameters[1] << " g= " << reaction_array[i].parameters[2] << endl;
				j++;
			}
		}
	}*/
	
	// Check if species have CR induced photodestruction reaction
	undef_species.clear();
	for (i = 0; i < nb_of_species; i++)
	{
		k = 0;
		for (l = 0; l < nb_of_reactions; l++) {
			if (reaction_array[l].reactant[0] == i) {
				if (species[i].type == "neutral") {
					if (reaction_array[l].type >= 4 && reaction_array[l].type <= 6)
						k++;
				}
			}
		}
		if (k == 0 && species[i].type == "neutral") { // only neutral are saved
			undef_species.push_back(species[i].name);
		}
	}	
	cout << "CR induced photoreactions are absent for neutral species: " << endl << "	";
	for (i = 0; i < (int) undef_species.size(); i++) {
		cout << undef_species[i] << " ";
	}
	cout << endl;

	// Check if species have interstellar UV induced photodestruction reaction
	undef_species.clear();
	for (i = 0; i < nb_of_species; i++)
	{
		k = 0;
		for (l = 0; l < nb_of_reactions; l++) {
			if (reaction_array[l].reactant[0] == i) {
				if (species[i].type == "neutral") {
					if (reaction_array[l].type >= 9 && reaction_array[l].type <= 11)
						k++;
				}
			}
		}
		if (k == 0 && species[i].type == "neutral") { // only neutral are saved
			undef_species.push_back(species[i].name);
		}
	}	
	cout << "Interstellar UV induced photoreactions are absent for neutral species: " << endl << "	";
	for (i = 0; i < (int) undef_species.size(); i++) {
		cout << undef_species[i] << " ";
	}
	cout << endl;

	// Check for identical reactions:
	l = 1;
	cout << "Identical reactions:" << endl;
	for (i = 0; i < nb_of_reactions; i++) {
		for (j = i+1; j < nb_of_reactions; j++) 
		{
			if (reaction_array[i] == reaction_array[j]) {
				cout << left << setw(4) << l << setw(7) << reaction_array[i].rcode 
					<< setw(30) << reaction_array[i].name << " " << reaction_array[j].name << endl;
				l++;
			}
		}
	}
	if (is_failed)
		exit(1);
}

void chem_network::print_network(const string &path)
{
	int i, j;
	string fname;
	ofstream output;

	fname = path + "chemistry.txt";
	output.open(fname.c_str());

	output << "The list of species taken into acount." << endl
		<< "index, name, mass, charge, formula (H He C N O Si S Fe Na Mg Cl P F), entalpy (erg);" << endl;

	output.setf(ios_base::scientific, ios_base::floatfield);
	output.precision(2);
	
	for (i = 0; i < (int) species.size(); i++)
	{
		output << left << setw(5) << i+1 << setw(15) << species[i].name << setw(8) << rounding(species[i].mass/ATOMIC_MASS_UNIT) 
			<< setw(3) << species[i].charge;
		
		for (j = 0; j < NB_OF_CHEM_ELEMENTS; j++) {
			output << left << setw(3) << species[i].formula[j];
		}	
		output << left << species[i].enthalpy << endl;
	}
	output << "Chemical reactions taken into account, nb = " << (int) reaction_array.size() << endl 
		<< "name, energy released (erg per reaction), type, min and max rates:" << endl;

	for (i = 0; i < (int) reaction_array.size(); i++) {
		output << left << setw(5) << i+1 << setw(45) << reaction_array[i].name << setw(11) << reaction_array[i].energy_released 
			<< setw(30) << get_reaction_type(reaction_array[i].type) << setw(10) << reaction_array[i].min_rate 
			<< setw(10) << reaction_array[i].max_rate << endl;
	}
}

int chem_network::find_specimen(const string &name) const
{
	for (int i = 0; i < (int) species.size(); i++) {
		if (species[i].name == name)
			return i;
	}
	return -1;
}

int chem_network::find_reaction(const string &name) const
{
	for (int i = 0; i < (int) reaction_array.size(); i++) {
		if (reaction_array[i].name == name)
			return i;
	}
	if (verbosity) 
		cout << "The reaction is not found: " << name << endl;
	return -1;
}

//
// Functions
//

string get_reaction_type(const int i)
{
	if (i >= 0 && i < NB_OF_CHEM_REACTION_TYPES)
		return chemical_reaction_types[i];
	return "undefined";
}

string get_reaction_name(string rcode, const vector<chem_specimen> &species, const int *reactant, int nb_of_reactants, 
	const int *product, int nb_of_products)
{
	int j;
	stringstream ss("");

	for (j = 0; j < nb_of_reactants; j++) 
	{
		ss << species[ reactant[j] ].name;
		if (j == nb_of_reactants-1) 
		{
			// special words are added (they are not included in the specimen list):
			if (rcode == "CP" || rcode == "CPD") 
				ss << " + CRP";
			else if (rcode == "CR" || rcode == "CRD" || rcode == "CRDGS1" || rcode == "CRDGS2")
				ss << " + CRPhoton";
			else if (rcode == "PH" || rcode == "PHD" || rcode == "ISDGS1" || rcode == "ISDGS2")
				ss << " + ISPhoton";	
			else if (rcode == "ADS" || rcode == "THD" || rcode == "IAG")
				ss << " + GRAIN";
			ss << " -> ";
		}
		else ss << " + ";
	}
	
	for (j = 0; j < nb_of_products; j++) 
	{
		ss << species[ product[j] ].name;
		if (j < nb_of_products-1)
			ss << " + ";
		else if (rcode == "ADS" || rcode == "THD" || rcode == "IAG")
			ss << " + GRAIN";
	}
	return ss.str();
}

// Rearrange the species order: adsorbed - 1, neutrals - 2, ions - 3, electrons - 4;
void check_order(const std::vector<chem_specimen> &species, int *arr, int nb)
{
	int i, j;
	int *aux = new int [nb];
	
	j = 0;
	for (i = 0; i < nb; i++) {
		if (species[ arr[i] ].type == "adsorbed") {
			aux[j] = arr[i];
			j++;
		}
	}
	for (i = 0; i < nb; i++) {
		if (species[ arr[i] ].type == "neutral") {
			aux[j] = arr[i];
			j++;
		}
	}
	for (i = 0; i < nb; i++) {
		if (species[ arr[i] ].type == "ion") {
			aux[j] = arr[i];
			j++;
		}
	}
	for (i = 0; i < nb; i++) {
		if (species[ arr[i] ].type == "e-") {
			aux[j] = arr[i];
			j++;
		}
	}
	memcpy(arr, aux, nb*sizeof(int));
	delete [] aux;
}

// In reactant and product lists adsorbed - 1, neutrals - 2, ions - 3, electrons - 4;
// photons are not included in the product order;
int react_type_determ(std::string rcode, const vector<chem_specimen> &species, const int *reactant, int nb_of_reactants, 
	const int *product, int nb_of_products, int verbosity)
{
	int i, type = -1;
	if (nb_of_reactants > 2) {
		if (verbosity)
			cout << "    Nb of reactants is > 2" << endl;
	}
	else if (rcode == "CP") {
		if (nb_of_products == 2) 
		{
			if (species[ product[0] ].type == "ion" && species[ product[1] ].type == "e-")
				type = 0;
			else if (species[ product[0] ].type == "ion" && species[ product[1] ].type == "ion")
				type = 2;
			else if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "neutral")
				type = 3;
		}
		else if (nb_of_products == 3) 
		{
			if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "ion" 
				&& species[ product[2] ].type == "e-")
				type = 1;
		}
	}
	// Photoreactions induced by photons produced by CR;
	else if (rcode == "CR") {
		if (species[ reactant[0] ].type == "neutral")
		{
			if (species[ product[0] ].type == "ion" && species[ product[1] ].type == "e-")
				type = 4;
			else if (nb_of_products == 3) 
			{
				if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "ion" 
					&& species[ product[2] ].type == "e-")
					type = 5;
			}

			for (i = 0; i < nb_of_products; i++) {
				if (species[ product[i] ].type != "neutral") break;
			}
			if (i == nb_of_products) 
				type = 6;
		}
		else if (species[ reactant[0] ].type == "ion") 
		{
			if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "ion")
				type = 7;
			else if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "e-")
				type = 8;
		}
	}
	// UV cosmic background induced photoreactions;
	else if (rcode == "PH") {
		if (species[ reactant[0] ].type == "neutral")
		{
			if (species[ product[0] ].type == "ion" && species[ product[1] ].type == "e-")
				type = 9;
			else if (nb_of_products == 3) 
			{
				if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "ion" 
					&& species[ product[2] ].type == "e-")
					type = 10;
			}

			for (i = 0; i < nb_of_products; i++) {
				if (species[ product[i] ].type != "neutral") break;
			}
			if (i == nb_of_products) 
				type = 11;
		}
		else if (species[ reactant[0] ].type == "ion") 
		{
			if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "ion")
				type = 12;
			else if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "e-")
				type = 13;
		}
	}
	// bimolecular reactions, 14 =< type <= 25;
	else if ((rcode == "AD" || rcode == "CD" || rcode == "CE" || rcode == "DR" || rcode == "IN" || rcode == "MN"
			|| rcode == "NN" || rcode == "RA" || rcode == "REA" || rcode == "RR") && (nb_of_reactants == 2)) 
	{
		if (species[ reactant[0] ].type == "neutral" && species[ reactant[1] ].type == "neutral") 
		{
			for (i = 0; i < nb_of_products; i++) {
				if (species[ product[i] ].type != "neutral") break;
			}
			if (i == nb_of_products && i >= 2) {
				type = 14;
			}
			else {
				if (nb_of_products == 1)
					type = 15;
				else if (nb_of_products == 2 && species[ product[0] ].type == "ion" && species[ product[1] ].type == "e-")
					type = 16;
			}
		}
		else if (species[ reactant[0] ].type == "neutral" && species[ reactant[1] ].name == "e-") 
		{
			if (nb_of_products == 1) {
				if (species[ product[0] ].type == "ion")
					type = 17;
			}
			else {
				for (i = 0; i < nb_of_products-1; i++) {
					if (species[ product[i] ].type != "neutral") break;
				}
				if (i == nb_of_products-1) {
					if (species[ product[i] ].type == "ion") 
						type = 18;
					else if (species[ product[i] ].type == "e-")
						type = 19;
				}
			}
		}
		else if (species[ reactant[0] ].type == "neutral" && species[ reactant[1] ].type == "ion")
		{
			if (nb_of_products == 1) 
				type = 20;

			if (nb_of_products >= 2) {
				for (i = 0; i < nb_of_products-1; i++) {
					if (species[ product[i] ].type != "neutral") 
						break;
				}
				if (i == nb_of_products-1) { 
					if (species[ product[i] ].type == "ion")
						type = 21;
					else if (species[ product[i] ].type == "e-")
						type = 22;
				}
			}
		}
		else if (species[ reactant[0] ].type == "ion" && species [ reactant[1] ].type == "e-") 
		{
			if (nb_of_products == 1) {
				if (species[ product[0] ].type == "neutral") {
					type = 23;
				}
			}
			else {
				for (i = 0; i < nb_of_products; i++) {
					if (species[ product[i] ].type != "neutral") break;
				}
				if (i == nb_of_products)
					type = 24;
			}
		}
		else if (species[ reactant[0] ].type == "ion" && species[ reactant[1] ].type == "ion") 
		{
			for (i = 0; i < nb_of_products; i++) {
				if (species[ product[i] ].type != "neutral") break;
			}
			if (i == nb_of_products)
				type = 25;
		}
	}
	// H2 formation on grains:
	else if (rcode == "HFG") {
		type = 26;
	}
	// Adsorption on grains
	else if (rcode == "ADS")
	{
		if (species[ reactant[0] ].type == "neutral")
			type = 27;
	}
	// Thermal desorption
	else if (rcode == "THD") {
		type = 28;
	}
	// CR desorption:
	else if (rcode == "CPD") {
		type = 29;
	}
	// Photodesorption by CR induced UV field:
	else if (rcode == "CRD") {
		if (nb_of_products == 1) {
			type = 30;
		}
		else if (nb_of_products == 2) {
			if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "neutral")
				type = 31;
		}
	}
	// Photodesorption by stellar UV field:
	else if (rcode == "PHD") {
		if (nb_of_products == 1) {	
			type = 32;
		}
		else if (nb_of_products == 2) {
			if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "neutral")
				type = 33;
		}
	}
	// sputtering of adsorbed species;
	else if (rcode == "SPU") {
		if (nb_of_products == 2) {	
			type = 34;
		}
		else if (nb_of_products == 3) {
			if (species[ product[0] ].type == "neutral" && species[ product[1] ].type == "neutral")
				type = 35;
		}
	}
	else if (rcode == "GS") {
		for (i = 0; i < nb_of_products; i++) {
			if (species[ product[i] ].type != "adsorbed") 
				break;
		}
		if (i == nb_of_products)
			type = 36;
		
		for (i = 0; i < nb_of_products; i++) {
			if (species[ product[i] ].type != "neutral") 
				break;
		}
		if (i == nb_of_products) {
			if (nb_of_products == 1)
				type = 37;
			else type = 38;			
		}

		if (nb_of_products == 2) {
			if (species[ product[0] ].type == "adsorbed" && species[ product[1] ].type == "neutral")
				type = 39;
		}
	}
	else if (rcode == "CRDGS1" || rcode == "CRDGS2" || rcode == "ISDGS1" || rcode == "ISDGS2") {
		for (i = 0; i < nb_of_products; i++) {
			if (species[ product[i] ].type != "adsorbed") 
				break;
		}
		if (i == nb_of_products) {
			if (rcode == "CRDGS1")
				type = 40;
			else if (rcode == "CRDGS2")
				type = 41;
			else if (rcode == "ISDGS1")
				type = 42;
			else type = 43;
		}
	}
	else if (rcode == "IAG") {
		type = 44;
	}
	else if (verbosity) {
		cout << "Unknown code of the reaction: " << rcode << endl;
	}
	
	if (type == -1 && verbosity) 
	{
		string name = get_reaction_name(rcode, species, reactant, nb_of_reactants, product, nb_of_products);
		cout << "Warning: unknown type of the reaction: " << name << endl;
	}
	return type;
}

double get_h_sticking_coeff(double gas_temp, double dust_temp)
{
	// Hollenbach & McKee, ApJS, 41, 555 (1979); see also Burke & Hollenbach, ApJ 265, p. 223 (1983);
	// return 1./(1. + 0.04*sqrt(gas_temp + dust_temp) + 2.e-3*gas_temp + 8.e-6*gas_temp*gas_temp);
	
	// Matar et al., J. Chem. Phys. 133, 104507, 2010;
	// gas temperature dependent sticking of hydrogen on cold amorphous water ice surfaces, S0 = 1., T0 = 52 K;
	return (1. + 0.0481*gas_temp)/pow(1. + 0.01923*gas_temp, 2.5);
}

double get_h2_sticking_coeff(double gas_temp)
{
	// Matar et al., J. Chem. Phys. 133, 104507, 2010; S0 = 0.76, T0 = 87 K;
	return (0.76 + 0.02184*gas_temp)/pow(1. + 0.01149*gas_temp, 2.5);
}

// mass is given in g;
double get_sticking_coeff(double gas_temp, double mass)
{
	// Matar et al., J. Chem. Phys. 133, 104507, 2010;
	// The critical temperature is set to be proportional to the mass of the particle;
	// T_0 = 0.5*m[a.m.u.]*T_0,H2; T_0,H2 = 87 K; S0 = 1;
	mass = 3.817e-26*gas_temp/mass; // = T_g/T_0 = T_g/(0.5*m[a.m.u.]*T_0,H2) = T_g *ATOMIC_MASS_UNIT/(0.5*m[g]*T_0,H2)
	return (1. + 2.5*mass)/pow(1. + mass, 2.5);
}

double get_sticking_coeff_he2016(double bind_en, double dust_temp, double beta, double gamma)
{
	// He et al., ApJ 823, p.56 (2016);
	// binding energies at low coverage must be used, but in the UMIST database binding energies 
	// at monolayer coverage are given (He et al., ApJ, 825, p.89, 2016), the ratio E_LC/E_ML = 1.5 is adopted;
	double x = exp(2.*beta*(dust_temp - 1.5*gamma*bind_en));
	return 0.5*(1. - (x - 1.)/(x + 1.));
}

double reaction_rate(const vector<chem_specimen> & species, const accretion_rate_functions *accr_func, const chem_reaction &reaction, 
	double cr_ioniz_rate, double uv_field_strength, double visual_extinct, double vel_n, double vel_i, 
	double ads_dust_area, double ads_grain_veln2, double temp_n, double temp_i, double temp_e, double temp_d, 
	double photoreact_factor_cr, double desorption_factor_cr, double photodes_factor_cr, double photodes_factor_is, 
	double conc_h_tot, double nb_ads_sites_grain)
{
	int i, j;
	double k(0.), tempr, teff, s, r;

	switch (reaction.type) 
	{
	case 0: // "A + CR -> Ion + e-"
	case 1: // "A + CR -> B + Ion + e-"
	case 2: // "A + CR -> Ion + Ion"
	case 3: // "A + CR -> B + C"
		k = reaction.parameters[0] *cr_ioniz_rate;
		break;

		// photodissociation or photoionization that proceed via excitation of molecular hydrogen; 
		// the case of CO must be treated individually because of self-shielding, UV extinction is dominated by CO bands and not by dust (Gredel et al. 1987);
		// see the discussion by Flower et al., A&A 474, 923 (2007);
		// Gredel et al. used dust cross section 2.e-21 cm2 per H (at UV range), standard N_H/A_V = 2.e+21 cm-2 per mag 
	case 4: // "A + CRPhoton -> Ion + e-"
	case 5: // "A + CRPhoton -> B + Ion + e-"
	case 6: // "A + CRPhoton -> B + C + D + E"
	case 7: // "Ion + CRPhoton -> B + Ion"
	case 8: // "Ion + CRPhoton -> A + e-"
	case 40: // "*A + CRPhoton -> *B + *C"
	case 41:
		// nb of fits equals 1 for all reactions of these type, temperature dependence exists only for CO dissociation;	
		if (fabs(reaction.parameters[1]) > DBL_EPSILON) {
			if (temp_n < reaction.temp_min[0])
				k = reaction.min_rate *photoreact_factor_cr;
			else
				k = reaction.parameters[0] *pow(temp_n, reaction.parameters[1]) *reaction.parameters[2] *photoreact_factor_cr;
		}
		else {
			k = reaction.parameters[0] *reaction.parameters[2] *photoreact_factor_cr;
		}
		break;

	case 9: // "A + ISPhoton -> Ion + e-"
	case 10: // "A + ISPhoton -> B + Ion + e-"
	case 11: // "A + ISPhoton -> B + C + D + E"
	case 12: // "Ion + ISPhoton -> B + Ion"
	case 13: // "Ion + ISPhoton -> A + e-"
	case 42: // "*A + ISPhoton -> *B + *C"
	case 43:
		k = uv_field_strength *reaction.parameters[0] *exp(-reaction.parameters[2]*visual_extinct);
		break;

	case 14: // "A + B -> C + D + E"
	case 15: // "A + B -> C + Photon"
	case 16: // "A + B -> Ion + e-"
		// it is assumed that maximal value of nb_of_fits is equal to 2;
		if (temp_n > reaction.temp_max[reaction.nb_of_fits-1]) {
			k = reaction.max_rate;
		}
		else if (temp_n < reaction.temp_min[0]) {
			k = reaction.min_rate;
		}
		else {
			if (reaction.nb_of_fits > 1) {
				if (temp_n > reaction.temp_min[1])
					i = reaction.nb_of_param;
				else i = 0;
			}
			else i = 0;

			k = reaction.parameters[i] *pow(temp_n, reaction.parameters[i+1]) *exp(-reaction.parameters[i+2]/temp_n);	
		}
		break;
	
	case 17: // "A + e- -> Ion"
	case 18: // "A + e- -> B + C + Ion"
	case 19: // "A + e- -> B + C + e-"
	case 23: // "Ion + e- -> A + Photon"
	case 24: // "Ion + e- -> A + B + C"
		if (temp_e > reaction.temp_max[reaction.nb_of_fits-1]) {
			k = reaction.max_rate;
		}
		else if (temp_e < reaction.temp_min[0]) {
			// this case is important for dissociative recombination reactions with beta < 0;
			k = reaction.min_rate;
		}
		else {
			if (reaction.nb_of_fits > 1) {
				if (temp_e > reaction.temp_min[1])
					i = reaction.nb_of_param;
				else i = 0;
			}
			else i = 0;
			k = reaction.parameters[i] *pow(temp_e, reaction.parameters[i+1]) * exp(-reaction.parameters[i+2]/temp_e);
		}
		break;
		
		// for all ion-ion reactions the rate is fitted by one formula, nb_of_fits = 1;
	case 25: // "Ion + Ion -> A + B + C + D"
		if (temp_i > reaction.temp_max[0])
			k = reaction.max_rate;
		else if (temp_i < reaction.temp_min[0])
			k = reaction.min_rate;
		else
			k = reaction.parameters[0] *pow(temp_i, reaction.parameters[1]) *exp(-reaction.parameters[2]/temp_i);	
		
		break;

	case 20: // "A + Ion -> Ion + Photon"
	case 21: // "A + Ion -> B + C + D + Ion"
	case 22: // "A + Ion -> B + C + e-"
		switch (reaction.cs_reconstr) {
		// exothermic neutral-ion reactions, the rate is constant;
		case 0:	
			k = reaction.parameters[0];
			break;
		
		// the method of effective temperature for neutral-ion reactions:
		// check if minimal and maximal rates are defined;
		case 1: 
			tempr = (temp_n*reaction.mass2 + temp_i*reaction.mass1)/reaction.mass_sum;
			s = vel_n - vel_i;

			// defenition of effective temperature by Flower et al., MNRAS 216, 775 (1985): Teff = T + m*v_d*v_d/(3k); 
			// by Draine, ApJ 241, 1021 (1980): Teff = T + m*v_d*v_d/(2k);
			teff = tempr + reaction.reduced_mass*s*s/(3.*BOLTZMANN_CONSTANT);
		
			if (teff > reaction.temp_max[reaction.nb_of_fits-1]) {
				k = reaction.max_rate;
			}
			else if (teff < reaction.temp_min[0]) {
				k = reaction.min_rate;
			}
			else {
				if (reaction.nb_of_fits > 1) {
					if (teff > reaction.temp_min[1])
						i = reaction.nb_of_param;
					else i = 0;
				}
				else i = 0;	
				k = reaction.parameters[i] *pow(teff, reaction.parameters[i+1]) *exp(-reaction.parameters[i+2]/teff);	
			}
			break;
		}
		break;

	case 26: // "H + H + grain -> H2 + grain"
		if (H2_FORMATION_MODE == 0) {
			k = 0.; // H2 formation is modelled as surface reaction;
		}
		else if (H2_FORMATION_MODE == 1) 
		{
			teff = temp_n + PI_DIVBY_EIGHT_BOLTZMANN_CONST*ATOMIC_MASS_UNIT*ads_grain_veln2;
			k = reaction.parameters[0] *sqrt(teff) *get_h_sticking_coeff(teff, temp_d) *ads_dust_area;
		}
		else { // here, semi-empirical approximation is used:
			k = reaction.parameters[0] *conc_h_tot;
		}
		break;
	
	case 27: // "A + grain -> *A + grain"
		// area of the grains that are able to absorb chemical species; area is n_gr*<pi*a*a>; 
		teff = temp_n + PI_DIVBY_EIGHT_BOLTZMANN_CONST*reaction.mass1 *ads_grain_veln2;
		k = reaction.parameters[0]*sqrt(teff);

		if (species[ reaction.reactant[0] ].name == "H") 
		{
			k *= get_h_sticking_coeff(teff, temp_d);
			if (H2_FORMATION_MODE)
				k *= 0.5;	// in case of "ad hoc" H2 formation;
		}
		else if (species[ reaction.reactant[0] ].name == "H2") {
			k *= get_h2_sticking_coeff(teff);
		}
		else {
			// common sticking coefficient for other species;
			k *= STICKING_COEFF_NEUTRALS*get_sticking_coeff(teff, species[ reaction.reactant[0] ].mass);
		}
		// see He et al., ApJ 823, p.56 (2016),
		// at this moment, all molecules have the same parameters beta and gamma; reactants have not binding energy initialized;
		k *= ads_dust_area *get_sticking_coeff_he2016(species[ reaction.product[0] ].bind_en, temp_d);
		break;
	
	// Thermal desorption:
	case 28: // "*A + grain -> A + grain"
		// Note: the temperature of the dust grains that are capable to adsorb must be used, 
		k = reaction.parameters[0]*exp(-reaction.parameters[1]/temp_d);
		break;

	// CR desorption:
	case 29: // "*A + CR -> A"
		// alternative approach; Roberts et al., MNRAS 382, 733–742, 2007;
        // see also Dartois et al., arXiv:2001.06349v1 (2020); Faure et al. MNRAS 487, 3392 (2019), appendix A;
		if (species[ reaction.reactant[0] ].bind_en < CR_DESORPTION_LIM_BENERGY)
			k = CR_DESORPTION_YIELD*desorption_factor_cr; 
		
		// standard approach:
		// k = reaction.parameters[0]*cr_ioniz_rate;
		break;

	// CR induced photodesorption:
	case 30: // "*A + CRPhoton -> A"
#ifdef COMMON_PHOTODES_YIELD
		k = COMMON_PHOTODES_YIELD*photodes_factor_cr;
#else
		k = reaction.parameters[0]*photodes_factor_cr;
#endif
		break;
	case 31: // "*A + CRPhoton -> A + B"
#ifdef COMMON_PHOTODES_YIELD
		k = 0.;
#else
		k = reaction.parameters[0]*photodes_factor_cr;
#endif
		break;

	// Stellar UV photodesorption:
	case 32: // "*A + ISPhoton -> A"
#ifdef COMMON_PHOTODES_YIELD
		k = COMMON_PHOTODES_YIELD*photodes_factor_is;
#else
		k = reaction.parameters[0]*photodes_factor_is;
#endif
		break;
	case 33: // "*A + ISPhoton -> A + B"
#ifdef COMMON_PHOTODES_YIELD
		k = 0.;
#else
		k = reaction.parameters[0]*photodes_factor_is;
#endif
		break;

	case 34: // "*A + C -> A + C"
	case 35: // "*A + C -> A + B + C"
		// Note, data in array reaction.parameters[] is not used here, mass of the projectile must be given:
		s = sqrt(0.5*reaction.mass2 *ads_grain_veln2/(BOLTZMANN_CONSTANT*temp_n));
		k = onedivby_8grain_sites_per_cm2 *reaction.rate_data->get(temp_n, s); // ->get() must return sputtering yield multiplied by thermal velocity; 
		break;

	case 36: // "*A + *B -> *C + *D + *E"
	case 37: // "*A + *B -> C"
	case 38: // "*A + *B -> C + D + E"
	case 39: // "*A + *B -> *A + B"
		i = reaction.reactant[0];
		j = reaction.reactant[1];
		
		switch (reaction.cs_reconstr) {
		case 0:
			k = species[i].surf_freq *exp(-species[i].diff_barrier/temp_d) 
				+ species[j].surf_freq *exp(-species[j].diff_barrier/temp_d);
			break;
		case 1:
			s = species[i].diff_barrier/temp_d;
			k = (reaction.parameters[2] < s) ? species[i].surf_freq*exp(-reaction.parameters[2]) : species[i].surf_freq*exp(-s);

			k += species[j].surf_freq *exp(-species[j].diff_barrier/temp_d);
			break;
		case 2:
			s = species[j].diff_barrier/temp_d;
			k = (reaction.parameters[3] < s) ? species[j].surf_freq*exp(-reaction.parameters[3]) : species[j].surf_freq*exp(-s);

			k += species[i].surf_freq *exp(-species[i].diff_barrier/temp_d);
			break;
		case 3:
			s = species[i].diff_barrier/temp_d;
			k = (reaction.parameters[2] < s) ? species[i].surf_freq*exp(-reaction.parameters[2]) : species[i].surf_freq*exp(-s);

			s = species[j].diff_barrier/temp_d;
			k += (reaction.parameters[3] < s) ? species[j].surf_freq*exp(-reaction.parameters[3]) : species[j].surf_freq*exp(-s);
			break;
		};
		
		// reaction with a barrier:
		if (reaction.parameters[1] > DBL_EPSILON)
		{
			r = species[i].surf_freq *exp(-species[i].bind_en/temp_d) // evaporation
				+ species[j].surf_freq *exp(-species[j].bind_en/temp_d);
			
			// choose fastest:
			s = reaction.parameters[1]/temp_d;
			s = (reaction.parameters[5] < s) ? reaction.parameters[4]*exp(-reaction.parameters[5]) : reaction.parameters[4]*exp(-s); 

			// competition among channels of multiple channel reactions:
			// Note, may be wrong: ask prof. Herbst about this moment;
			if (reaction.parameters[6] > DBL_EPSILON) {
				r += (reaction.parameters[7] < reaction.parameters[6]/temp_d) 
					? reaction.parameters[4]*exp(-reaction.parameters[7]) 
					: reaction.parameters[4]*exp(-reaction.parameters[6]/temp_d);

				if (reaction.parameters[8] > DBL_EPSILON) {
					r += (reaction.parameters[9] < reaction.parameters[8]/temp_d) 
					? reaction.parameters[4]*exp(-reaction.parameters[9]) 
					: reaction.parameters[4]*exp(-reaction.parameters[8]/temp_d);
				}
			}
			k = k*s/(k + r + s); // the competition among reaction, diffusion and evaporation (Ruaud et al., MNRAS 459, p.3756, 2016);
		}
		else if (reaction.type == 39) // desorption encounter mechanism, without activation barrier (!)
		{
			s = species[j].surf_freq *exp(-reaction.parameters[4]/temp_d); // thermal desorption rate
			r = species[j].surf_freq *exp(-reaction.parameters[5]/temp_d); // hoping rate
			
			// Note, in Hincelin et al. (2014), the diffusion rate is used, not hoping rate, 
			// the hoping rate must be divided by number of sites on grain:
			k = k*s/(s + r/nb_ads_sites_grain);
		}
		k *= reaction.parameters[0]/(4.*GRAIN_SITES_PER_CM2 *ads_dust_area); // ads_dust_area = pi*a*a*n_g
		
		break;
	default: break;
	}	
	return k;
}

// Temperature is in erg;
// Introduce the efficiency coefficient of energy defect conversion to kinetic energy
void chemistry_source_terms(double & mom_gain_n, double & mom_gain_i, double & mom_gain_e, double & energy_gain_n, 
	double & energy_gain_i, double & energy_gain_e, double & energy_gain_d, const vector<chem_specimen> & species, 
	const chem_reaction & reaction, double vel_n, double vel_i, double ads_grain_velz, double ads_grain_veln2, 
	double temp_n, double temp_i, double temp_e, double temp_d, double rate)
{
	double a, b, m, ma, mi, mj, t;
	// probably, the energy gain in photoreactions and photon-producing reactions must be taken into account accurately;
	// the energy gain includes the translational motions of the molecule; 
	// heat capacity of the neutral and ion gas is taken to be equal 3/2*k_B, (3 degrees of freedom);
	switch(reaction.type) 
	{
	case 0: // "A + CR -> Ion + e-"
		a = species[ reaction.product[0] ].mass *rate;
		mom_gain_n -= a*vel_n;
		mom_gain_i += a*vel_n;
		mom_gain_e += ELECTRON_MASS *vel_n *rate;
	
		energy_gain_n -= 1.5*temp_n *rate;
		energy_gain_i += 1.5*temp_n *rate + 0.5*(vel_n - vel_i)*(vel_n - vel_i)*a;
		// the kinetic energy of the electron goes to: 
		// a) the heating of electron gas and b) excitation of H2 molecules - are taken into account separately;
		energy_gain_e += 1.5*temp_e *rate; // minimal value of electron energy
		break;

	case 1: // "A + CR -> B + Ion + e-"	
		// It is assumed that only one neutral and one ion form, ion has a second place in the product order;
		ma = species[ reaction.reactant[0] ].mass;
		mi = species[ reaction.product[1] ].mass;
		
		a = mi *rate;
		mom_gain_n -= a *vel_n;
		mom_gain_i += a *vel_n;
		mom_gain_e += ELECTRON_MASS *vel_n *rate;
		
		energy_gain_n += -1.5*temp_n*mi/ma *rate;
		energy_gain_i += (1.5*temp_n/ma + 0.5*(vel_n - vel_i)*(vel_n - vel_i))*mi*rate;
		energy_gain_e += 1.5*temp_e *rate;
		break;
	
	case 2: // "A + CR -> Ion + Ion"
		a = species[ reaction.reactant[0] ].mass;
		b = a*rate*vel_n;
		
		mom_gain_n -= b;
		mom_gain_i += b;
		
		energy_gain_n -= 1.5*temp_n *rate;
		energy_gain_i += (1.5*temp_n + 0.5*a *(vel_n - vel_i)*(vel_n - vel_i))*rate;
		break;
	
	case 3: // "A + CR -> B + C"
		break;
	
	case 4: // "A + CRPhoton -> Ion + e-"
	case 9: // "A + ISPhoton -> Ion + e-"
		a = species[ reaction.product[0] ].mass;
		b = a*rate*vel_n;
		
		mom_gain_n -= b;
		mom_gain_i += b;
		mom_gain_e += ELECTRON_MASS *vel_n *rate;

		energy_gain_n -= 1.5*temp_n *rate;
		energy_gain_i += (1.5*temp_n + 0.5*a *(vel_n - vel_i)*(vel_n - vel_i))*rate;
		// What is about electron energy gain?
		energy_gain_e += 1.5*temp_e *rate; // - minimal value of electron energy
		break;

	case 5: // "A + CRPhoton -> B + Ion + e-"
	case 10: // "A + ISPhoton -> B + Ion + e-"
		m = species[ reaction.reactant[0] ].mass;

		mi = species[ reaction.product[1] ].mass; // ion has a right place in the product order;
		b = mi*rate*vel_n;

		mom_gain_n -= b;
		mom_gain_i += b;
		mom_gain_e += ELECTRON_MASS *vel_n *rate;

		energy_gain_n -= 1.5*temp_n *mi/m *rate;
		energy_gain_i += (1.5*temp_n/m + 0.5*(vel_n - vel_i)*(vel_n - vel_i))*mi *rate;
		energy_gain_e += 1.5*temp_e *rate;
		break;
	
	case 6: // "A + CRPhoton -> B + C + D + E"
	case 11: // "A + ISPhoton -> B + C + D + E"
		// the energy gain of neutrals due to photon absorption is not taken into account;
		break;
	
	case 7: // "Ion + CRPhoton -> B + Ion"
	case 12: // "Ion + ISPhoton -> B + Ion"
		m = species[ reaction.reactant[0] ].mass;
		
		// neutral has a first place in the product order;
		a = species[ reaction.product[0] ].mass;
		b = a*rate*vel_i;

		mom_gain_n += b;
		mom_gain_i -= b;

		// the energy gain due to photon excitation is not taken into account;
		energy_gain_n += (0.5*(vel_n - vel_i)*(vel_n - vel_i) + 1.5*temp_i/m)*a*rate;
		energy_gain_i -= 1.5*temp_i*a/m *rate;
		break;

	case 8: // "Ion + CRPhoton -> A + e-"
	case 13: // "Ion + ISPhoton -> A + e-"
		// Only one neutral is formed here, otherwise expressions must be changed;
		a = species[ reaction.product[0] ].mass;
		b = a*rate*vel_i;

		mom_gain_n += b;
		mom_gain_i -= b;
		mom_gain_e += ELECTRON_MASS *vel_i *rate;

		energy_gain_i -= 1.5*temp_i *rate;
		energy_gain_n += (1.5*temp_i + 0.5*a *(vel_n - vel_i)*(vel_n - vel_i)) *rate;
		energy_gain_e += 1.5*temp_e *rate;
		break;
			
	case 14: // "A + B -> C + D + E"
		energy_gain_n += reaction.energy_released *rate;
		break;
	
	case 15: // "A + B -> C + Photon"
		energy_gain_n -= 1.5*temp_n *rate;
		break;
	
	case 16: // "A + B -> Ion + e-", there is only one reaction of this type, is exothermic, CH + O -> HCO+ + e-
		a = species[ reaction.product[0] ].mass;
		b = a *rate*vel_n;

		mom_gain_n -= b;
		mom_gain_i += b;
		mom_gain_e += ELECTRON_MASS *vel_n *rate;

		energy_gain_n -= 3.*temp_n *rate;
		energy_gain_i += (1.5*temp_n + 0.5*a *(vel_n - vel_i)*(vel_n - vel_i))*rate;
		energy_gain_e += (1.5*temp_n + reaction.energy_released)*rate;		
		break;
	
	case 17: // "A + e- -> Ion"
		a = species[ reaction.product[0] ].mass;
		b = a*rate*vel_n;

		mom_gain_n -= b;
		mom_gain_i += b;
		mom_gain_e -= ELECTRON_MASS *vel_i *rate;

		energy_gain_n -= 1.5*temp_n *rate;
		energy_gain_i += (1.5*temp_n + 0.5*a *(vel_n - vel_i)*(vel_n - vel_i))*rate;
		energy_gain_e -= 1.5*temp_e *rate; // kinetic energy is not conserved here;
		break;
	
	case 18: // "A + e- -> B + C + Ion"
		// The energy partition between ion and neutrals is treated approximately: as B and C form a complex;
		ma = species[ reaction.reactant[0] ].mass;
		mi = species[ reaction.product[reaction.nb_of_products-1] ].mass;

		b = mi*rate*vel_n;
		mom_gain_n -= b;
		mom_gain_i += b;
		mom_gain_e -= ELECTRON_MASS *vel_i *rate;

		energy_gain_n += (-1.5*temp_n + 1.5*temp_e + reaction.energy_released)*mi/ma *rate;
		energy_gain_i += ( (1.5*temp_n*mi + (1.5*temp_e + reaction.energy_released)*(ma - mi))/ma + 0.5*mi*(vel_n - vel_i)*(vel_n - vel_i) )*rate;
		energy_gain_e -= 1.5*temp_e *rate;
		break;
	
	case 19: // "A + e- -> B + C + e-"
		// The energy_released goes to the electron, the minimal value of the parameter is equal to the difference
		// of the enthalpies of the products B and C and reactant A (not including their kinetic energies);
		// there is only one reaction of this type, endothermic, H2 + e- -> H + H + e-
		energy_gain_e += reaction.energy_released *rate;
		break;
	
	case 20: // "A + Ion -> Ion + Photon"
		m = species[ reaction.product[0] ].mass;		
		a = species[ reaction.reactant[0] ].mass;
		b = a *rate*vel_n;

		mom_gain_n -= b;
		mom_gain_i += b;
		
		energy_gain_n -= 1.5*temp_n *rate;
		energy_gain_i += 0.5*a/m *( 3.*(temp_n - temp_i) + a*(vel_n - vel_i)*(vel_n - vel_i) )*rate;
		break;
	
	case 21: // "A + Ion -> B + C + D + Ion"
		// The neutrals are treated as a single particle;
		ma =  species[ reaction.reactant[0] ].mass;
		mi = species[ reaction.reactant[1] ].mass;
		mj = species[ reaction.product[reaction.nb_of_products-1] ].mass;
		m = ma + mi;

		b = (-ma*mj*vel_n + mi*(m - mj) *vel_i)/m *rate;
		mom_gain_n += b;
		mom_gain_i -= b;
		
		t = 1.5*(temp_i - temp_n)*((m - mj)*mi + mj*ma)/(m*m);

		energy_gain_n += ( t + reaction.energy_released*mj/m + 0.5*(ma*mi*mj + mi*mi*(m - mj)) *(vel_n - vel_i)*(vel_n - vel_i)/(m*m) )*rate;
		energy_gain_i += (-t + reaction.energy_released*(m - mj)/m + 0.5*(ma*mi*(m - mj) + ma*ma*mj) *(vel_n - vel_i)*(vel_n - vel_i)/(m*m) )*rate;
		break;
	
	case 22: // "A + Ion -> B + C + e-"
		// The expressions for neutral fluid are derived based on those from previous section, assuming that B and C form a complex;
		ma = species[ reaction.reactant[0] ].mass;
		mi = species[ reaction.reactant[1] ].mass;
		m = ma + mi;

		b = mi *rate*vel_i;
		mom_gain_n += b;
		mom_gain_i -= b;
		mom_gain_e += ELECTRON_MASS*(vel_i*mi + vel_n*ma)/m *rate;

		energy_gain_n += 0.5*mi/m*( 3*(temp_i - temp_n)  + mi*(vel_n - vel_i)*(vel_n - vel_i) )*rate;
		energy_gain_i -= 1.5*temp_i *rate;

		b = ((1.5*(temp_n *mi + temp_i*ma) + 0.5*ma*mi*(vel_n - vel_i)*(vel_n - vel_i))/m + reaction.energy_released)*rate;
		energy_gain_e += (b > 0.) ? b : 0.;
		break;
	
	case 23: // "Ion + e- -> A + Photon"
		a = species[ reaction.product[0] ].mass;
		b = a*rate*vel_i;

		mom_gain_n += b;
		mom_gain_i -= b;
		mom_gain_e -= ELECTRON_MASS*vel_i *rate;

		energy_gain_n += (1.5*temp_i + 0.5*a*(vel_n - vel_i)*(vel_n - vel_i))*rate;
		energy_gain_i -= 1.5*temp_i*rate;
		energy_gain_e -= temp_e*rate; // see Draine, MNRAS 220, 133 (1986); 
		break;
	
	case 24: // "Ion + e- -> A + B + C"
		a = species[ reaction.reactant[0] ].mass;
		b = a*rate*vel_i;

		mom_gain_n += b;
		mom_gain_i -= b;
		mom_gain_e -= ELECTRON_MASS*vel_i*rate;

		// it is considered that about 0.3 energy defect goes to the kinetic energy of products;
		// Hamberg et al., J. Phys. Chem. 118, p. 6034 (2014); Ojekull et al., J. Chem. Phys. 120, p. 22, (2004); 
		energy_gain_n += (1.5*temp_i + temp_e + 0.5*a*(vel_n - vel_i)*(vel_n - vel_i) + 0.3*reaction.energy_released)*rate;
		energy_gain_i -= 1.5*temp_i *rate;
		energy_gain_e -= temp_e *rate;
		break;
	
	case 25: // "Ion + Ion -> A + B + C + D"
		m = species[ reaction.reactant[0] ].mass + species[ reaction.reactant[1] ].mass;
		
		b = m *rate*vel_i;
		mom_gain_n += b;
		mom_gain_i -= b;

		energy_gain_n += (3.*temp_i + 0.5*m*(vel_n - vel_i)*(vel_n - vel_i) + reaction.energy_released) *rate;
		energy_gain_i -= 3.*temp_i *rate;
		break;

	case 26: // "H + H + grain -> H2 + grain"
		// it is assumed that 1/3 of the released energy goes into the kinetic energy of H2 molecule; 
		// 1/3 of energy into heating the dust; 1/3 to H2 molecule excitation (last term is taken into account separately)
		// the kinetic energy of adsorbed H atoms is not taken into account here;
		// source term for momentum is ignored;
		b = species[reaction.product[0]].mass *rate;
		a = 0.333333*reaction.energy_released*rate;
		energy_gain_n += a + 0.5*b *ads_grain_veln2;
		energy_gain_d += a;
		break;
	
	case 27: // "A + grain -> *A + grain"
		// source terms for momentum and energy are ignored (they are taken into account separately, Draine, ApJ 241, p. 1021, 1980);
		// energy_gain_n -= 1.5*temp_n*rate;
		energy_gain_d += reaction.energy_released *rate;
		break;
	
	case 28: // "*A + grain -> A + grain"
		// source terms for momentum and energy are ignored;
		// a = species[ reaction.reactant[0] ].mass *rate;
		// energy_gain_n += 1.5*temp_d*rate + 0.5*a*ads_grain_veln2;
		energy_gain_d += reaction.energy_released*rate;
		break;
	
	// Photodesorption by CR particles, CR induced UV field, photodesorption by stellar UV field:
	// source terms for momentum and energy are ignored;
	case 29: // "*A + CR -> A"
	case 30: // "*A + CRPhoton -> A"
	case 32: // "*A + ISPhoton -> A"
	case 31: // "*A + CRPhoton -> A + B"
	case 33: // "*A + ISPhoton -> A + B"
		// a = species[ reaction.reactant[0] ].mass *rate;
		// energy_gain_n += 0.5*a*ads_grain_veln2 + 1.5*reaction.nb_of_products*temp_d*rate;
		break;
	
	// sputtering
	case 34: // "*A + C -> A + C"
	case 35: // "*A + C -> A + B + C"
		a = species[ reaction.reactant[0] ].mass *rate;
		energy_gain_n += 0.5*a*ads_grain_veln2 + rate *reaction.energy_released; 
		break;
	
	case 36: // "*A + *B -> *C + *D + *E"
		energy_gain_d += rate *reaction.energy_released;
		break;

	case 37: // "*A + *B -> C"
	case 38: // "*A + *B -> C + D + E"
		a = (species[ reaction.reactant[0] ].mass + species[ reaction.reactant[1] ].mass)*rate;
		energy_gain_n += 0.5*a*ads_grain_veln2 + 0.333333*rate*reaction.energy_released; // as for H2 molecule formation:
		energy_gain_d += 0.333333*rate*reaction.energy_released;
		break;

		// reaction products: one adsorbed and one neutral specimen;
	case 39: // "*A + *B -> *A + B"
		// a = rate*species[ reaction.product[1] ].mass;
		// energy_gain_n += 0.5*a*ads_grain_veln2 + 1.5*temp_d*rate;
		energy_gain_d += rate *reaction.energy_released;

	case 40: // "*A + CRPhoton -> *B + *C"
	case 41:
	case 42: // "*A + ISPhoton -> *B + *C"
	case 43:
		break;

	default: break;
	}
}

void init_chem_abund(const std::string fname, const chem_network *network, double *abundances)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int nb;
	double a, charge_density;
	
	string sn;
	stringstream ss;
	ifstream input;

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with initial abundances of chemical species " << fname << endl;
		exit(1);
	}
	charge_density = 0.;
	memset(abundances, 0, network->nb_of_species*sizeof(double));
	
	// Note: the precision of the saved values is crucial for physical parameter balance;
	while (!input.eof())
	{
		// comment lines are read:
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		if (text_line[0] == '\0')
			break;

		ss.clear();
		ss.str(text_line);
		ss >> nb >> sn >> a;
		
		nb = network->find_specimen(sn);
		if (nb >= 0) 
		{
			abundances[nb] = a;
			if (network->species[nb].charge != 0 && nb != network->e_nb)
				charge_density += a*network->species[nb].charge;
		}
		
		if (ss.eof() || ss.fail() || ss.bad()) {
			cout << "Error ocurred while reading the file " << fname << endl;
		}
	}
	// the data must be checked that charge_density > 0;
	abundances[network->e_nb] = charge_density;
	input.close();
}


struct rdata {
	string s[6];
	double act_en, width;
	bool operator==(const rdata & obj)
	{ return (s[0] == obj.s[0] && s[1] == obj.s[1] && s[2] == obj.s[2] 
		&& s[3] == obj.s[3] && s[4] == obj.s[4] && s[5] == obj.s[5]); }
};

void reformat_chemical_data_Belloche2014(const std::string &path)
{
	rdata rd;
	vector<rdata> rd_v;

	const int max_nb_reactants = 2, max_nb_products = 5;
	char text_line[MAX_TEXT_LINE_WIDTH];
	
	int i, j, k;
	double act_en, width;
	
	string str, fname, sM = "CH2", sQ = "NH2", sY = "CH2OH", sX = "CH3O", sZ = "CO";
	ifstream input1, input2;
	ofstream output;

	fname = path + "chemistry/Belloche_Garrod_2014/gs_reactions.txt";
	input1.open(fname.c_str(), ios_base::in);
	if (!input1) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with chemical reactions " << fname << endl;
		exit(1);
	}
	fname = path + "chemistry/Belloche_Garrod_2014/activation_energies.txt";
	input2.open(fname.c_str(), ios_base::in);
	if (!input2) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with chemical reactions " << fname << endl;
		exit(1);
	}
	
	while (!input1.eof())
	{	
		do // comment lines are read:
			input1.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#' || text_line[0] == '!');
		
		if (text_line[0] == '\0') // empty line at the file end;
			break;
	
		rd.s[0] = rd.s[1] = rd.s[2] = rd.s[3] = rd.s[4] = rd.s[5] = "";
		for (j = 0; j < max_nb_reactants; j++)
		{
			i = k = j*8;
			str = "";
			while (text_line[i] != ' ' && i - k < 8) {
				str.push_back(text_line[i]);
				i++;
			}	
			if (i > k) {
				rd.s[j] = str;
			}
		}
		
		for (j = 0; j < max_nb_products; j++) 
		{
			i = k = 24 + j*8;
			str = "";
			while (text_line[i] != ' ' && i - k < 8) {
				str.push_back(text_line[i]);
				i++;
			}	
			if (i > k) {
				rd.s[max_nb_reactants + j] = str;
			}
		}
		rd.act_en = rd.width = 0.;
		rd_v.push_back(rd);
	}
	input1.close();

	while (!input2.eof())
	{	
		do // comment lines are read:
			input2.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#' || text_line[0] == '!');
		
		if (text_line[0] == '\0') // empty line at the file end;
			break;
	
		rd.s[0] = rd.s[1] = rd.s[2] = rd.s[3] = rd.s[4] = rd.s[5] = "";
		for (j = 0; j < max_nb_reactants; j++)
		{
			i = k = j*8;
			str = "";
			while (text_line[i] != ' ' && i - k < 8) {
				str.push_back(text_line[i]);
				i++;
			}	
			if (i > k) {
				rd.s[j] = str;
			}
		}
		
		for (j = 0; j < max_nb_products; j++) 
		{
			i = k = 16 + j*8;
			str = "";
			while (text_line[i] != ' ' && i - k < 8) {
				str.push_back(text_line[i]);
				i++;
			}	
			if (i > k) {
				rd.s[max_nb_reactants + j] = str;
			}
		}
		k = 49;
		str = "";
		for (i = k; i < k + 10; i++) {
			str.push_back(text_line[i]);
		}
		act_en = atof(str.c_str());
		
		k = 59;
		str = "";
		for (i = k; i < k + 10; i++) {
			str.push_back(text_line[i]);
		}
		width = atof(str.c_str());

		for (i = 0; i < (int) rd_v.size(); i++) {
			if (rd_v[i] == rd) {
				rd_v[i].act_en = act_en;
				rd_v[i].width = width;
			}
		}
	}
	input2.close();

	fname = path + "chemistry/Belloche_Garrod_2014/bimolecular_gs_reactions_.txt";
	output.open(fname.c_str(), ios_base::out);
	output << "# Network of grain surface reactions," << endl
		<< "# Belloche & Garrod, Science 345, p.1584, 2014; http://www.astro.cornell.edu/~rgarrod/resources/" << endl
		<< "# changes of specimen names must be done:" << endl
		<< "# C3H3-> ?; C4H2 -> HC4H (diacetylene); C4H4 -> CH2CHCCH; C5H4 -> ?; C2H6 -> CH3CH3 (ethane); COOH -> HOCO" << endl;

	for (i = 0; i < (int) rd_v.size(); i++) {
		output << left << setw(15) << rd_v[i].s[0] << setw(15) << rd_v[i].s[1] << setw(15) << rd_v[i].s[2] <<
			setw(15) << rd_v[i].s[3] << setw(15) << rd_v[i].s[4] << setw(15) << rd_v[i].s[5] <<
			setw(11) << "1." << setw(11) << rd_v[i].width << setw(11) << rd_v[i].act_en;
		
		if (i+1 < (int) rd_v.size())
			output << endl;
	}
	output.close();
}

/*				if (text_line[i] == 'M' && text_line[i+1] != 'g')
					str = str + sM;
				else if (text_line[i] == 'Q')
					str = str + sQ;
				else if (text_line[i] == 'Y')
					str = str + sY;
				else if (text_line[i] == 'X')
					str = str + sX;
				else if (text_line[i] == 'Z')
					str = str + sZ;
				else */
/*
if (reaction.name == "*H + *H -> *H2") // this reaction is ignored: rate is equal 0.
							reaction.parameters[0] = 0.;
						
						else if (reaction.name == "*H + *H -> H2") // H2 molecule formed is desorbed immediately;
							reaction.parameters[0] = 1.;
						
// The cross section is evaluated by the method by Draine, Katz, ApJ 306 (1986) (must be checked carefully before usage),
						// there is no restriction on temperature:
						if (REACTION_RATE_CALC_CS && 
							nb_of_fits == 1 && reaction.parameters[1] > -1.5 && reaction.parameters[2] > DBL_EPSILON)
						{
							if (verbosity)
								cout << left << "Evaluating cross section: " << setw(25) << reaction.name 
									<< " " << reaction.parameters[1] << endl;
							
							reaction.cs_reconstr = 1;
							cross_section = new chem_reaction_cross_section_1(reaction.mass1, reaction.mass2, 
								reaction.parameters[0], reaction.parameters[1], reaction.parameters[2]);

							rate_data = new reaction_rate_data();
							rate_data->calc_data(cross_section);

							rate_data_array.push_back(*rate_data);
							reaction.rate_data = rate_data;
							
							delete cross_section;
						}

						
		case 1: 
			// endothermic neutral-ion reaction; rate reconstruction by the method by Draine & Katz (1986);
			// parameter mass1 corresponds to the neutral particle; for this case nb_of_fits = 1 by definition;
			tempr = (temp_n*reaction.mass2 + temp_i*reaction.mass1) /reaction.mass_sum;
			s = sqrt(0.5*reaction.reduced_mass/(BOLTZMANN_CONSTANT*tempr)) *fabs(vel_n - vel_i);
	
			k = reaction.rate_data->get(tempr, s);
		
			if (k > reaction.max_rate) 
				k = reaction.max_rate;
			break;

case 43: // "A + Ion -> Ion + Ion + e-"
		// the reaction will preferentially deplete particular parts of the reactant velocity distribution - it is neglected. 
		// the reaction is intense due to ion-neutral streaming, not thermal motion, KT << m(v_n-v_i)^2,
		ma = species[ reaction.reactant[0] ].mass;
		mi = species[ reaction.reactant[1] ].mass;
		m = ma + mi;

		b = ma*vel_n*rate;
		mom_gain_n -= b;
		mom_gain_i += b;
		mom_gain_e += ELECTRON_MASS*rate*vel_n;
	
		energy_gain_n -= 1.5*temp_n*rate;
		energy_gain_i += (0.5*ma*(vel_n - vel_i)*(vel_n - vel_i) + 1.5*temp_n + reaction.energy_released)*rate;
		energy_gain_e += 0.5*ELECTRON_MASS*(vel_n - vel_i)*(vel_n - vel_i)*rate;
		break;
*/
/*
	output << endl << "Endothermic neutral-ion reactions." << endl 
		<< "Reconstructed cross section (Draine, Katz, ApJ 306, 1986) vs effective temperature method:" << endl;
	
	for (j = 0; j < snb; j++) {
		output << left << setw(15) << sa[j];
	}
	output << endl;
	
	k = 1;
	for (i = 0; i < (int) reaction_array.size(); i++) 
	{
		if ((reaction_array[i].type == 20 || reaction_array[i].type == 21 || reaction_array[i].type == 22) 
			&& reaction_array[i].cs_reconstr == 2) 
		{
			output << left << setw(5) << k << reaction_array[i].name << endl;
			// here, temperature is in K;
			for (t = 10.; t < 0.5*reaction_array[i].temp_max[0]; t *= 1.079775)
			{
				output << left << setw(15) << t;
				for (j = 0; j < snb; j++)
				{
				// defenition of effective temperature by Flower et al., MNRAS 216, 775 (1985): 
				//		Teff = T + m*v_d*v_d/(3k) = T*(1. + 2*s*s/3), s*s = v_d*v_d/(2kT/m); 
				// by Draine, ApJ 241, 1021 (1980): 
				//		Teff = T + m*v_d*v_d/(2k) = T*(1 + s*s);
					et = t *(1. + 2.*sa[j]*sa[j]/3.);
					
					r1 = reaction_array[i].rate_data->get(t, sa[j]);
					if (r1 > reaction_array[i].max_rate)
						r1 = reaction_array[i].max_rate;

					r2 = reaction_array[i].rate_data->get(et, 0);
					if (r2 > reaction_array[i].max_rate)
						r2 = reaction_array[i].max_rate;

					output << setw(15) << r1 << setw(15) << r2;
				}
				output << endl;
			}
			k++;
		}
	}
	*/

