#pragma once

#include <string>
#include <cstring>

// the array of gas-phase chemical species is given in input, adsorbed species are considered automatically in the routine;
void production_routes(std::string path);

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
