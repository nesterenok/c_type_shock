#pragma once

#include <string>
#include <cstring>
#include <vector>
#include <algorithm>

// the array of gas-phase chemical species is given in input, adsorbed species are considered automatically in the routine;
// path1 - the location of files with chemical species and reactions, 
// path2 - the location of data file, the output is written in this location
void production_routes(std::string path1, std::string path2);
void add_cosmicray_chemistry(std::vector<std::string>& specimen_names);
void add_oxygen_chemistry(std::vector<std::string>& specimen_names);
void add_carbon_chemistry(std::vector<std::string>& specimen_names);
void add_nitrogen_chemistry(std::vector<std::string>& specimen_names);


// The time grid in NAUTILUS data and in the compared data must be the same:
void nautilus_comparison(std::string path);

class reaction_data 
{
private:
	static int nb_of_reactants, nb_of_products, max_nb_of_rate_values;
	
public:
	
	int *reactants, *products;
	double *rates;
	std::string name;

	reaction_data();
	reaction_data(const reaction_data &);
	reaction_data & operator=(const reaction_data &);
	~reaction_data();

	friend bool operator == (const reaction_data &, const reaction_data &);
	friend bool operator != (const reaction_data &, const reaction_data &);
	friend bool operator < (const reaction_data &, const reaction_data &);
	friend bool operator > (const reaction_data &, const reaction_data &);
};

/*
class depth_temperature_dependence
{
private:
    int prev_index;
    std::vector<double>::iterator low;
    std::vector<double> depth, gas_temperature;

public:

    double get_gas_temperature(double depth);
    depth_temperature_dependence(std::string path);
};*/
