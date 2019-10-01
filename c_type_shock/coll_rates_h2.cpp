//
// 04.03.2017. Check for errors. Errors were found in the data file for e-H2 collisions.
// 10.09.2017. Check for errors. Error was found in the H2-e collision data.
//		Note: there is a problem when number of H2 levels is small;

#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cfloat>
#include <sstream>

#include "constants.h"
#include "utils.h"
#include "linear_algebra.h"
#include "coll_rates_h2.h"

#define MAX_TEXT_LINE_WIDTH 3000
#define SOURCE_NAME "coll_rates_h2.cpp"
using namespace std;

//
// Classes with collisional coefficient data
//

h2_oh2_flower_data::h2_oh2_flower_data(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, li, lf, f, nb, nb_lines;
	double **temp;

	string file_name;
	ifstream input;

	if (rounding(h2_di->mol.spin) == 0) {
		file_name = data_path + "coll_h2/coll_ph2_oh2.txt";
	}
	else if (rounding(h2_di->mol.spin) == 1) {
		file_name = data_path + "coll_h2/coll_oh2_oh2.txt";
	}
	else {
		// the file, which contain rate coefficients for both ortho- and para-H2:
		file_name = data_path + "coll_h2/coll_h2_oh2.txt";
	}
	
	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) { 
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev;
	
	jmax = 10; // the T = 0 K point is included;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	temp = alloc_2d_array<double>(nb_lev*nb_lev, jmax);
	memset(*temp, 0, nb_lev*nb_lev *jmax *sizeof(double));

	nb_lines = nb_lev*nb_lev;
	tgrid[0] = 0.; // 0 K point, what one can say about rate coefficients at small T?

	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> li >> lf >> f >> f;
		for (j = 1; j < jmax; j++) {
			input >> temp[i][j];
		}
	}
	input.close();
	
	for (li = 1; li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			i = li *nb_lev + lf; // li -> lf
			f = lf *nb_lev + li; // lf -> li
			nb = li*(li-1)/2 + lf;

			if (li < h2_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *exp( (h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j] ) 
						*h2_di->lev_array[lf].g / ((double) h2_di->lev_array[li].g)); 
				}
			}
		}
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);
	
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_ph2_flower_data::h2_ph2_flower_data(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, li, lf, f, nb, nb_lines;
	double **temp;

	string file_name;
	ifstream input;

	if (rounding(h2_di->mol.spin) == 0) {
		// 27 level data were completed by zero rows; coefficients for 28->i were calculated from i->28 and added;
		file_name = data_path + "coll_h2/coll_ph2_ph2.txt";
	}
	else if (rounding(h2_di->mol.spin) == 1) {
		// the zero row was added for the H2 transition 1 -> 1;
		file_name = data_path + "coll_h2/coll_oh2_ph2.txt";
	}
	else {
		file_name = data_path + "coll_h2/coll_h2_ph2.txt";
	}

	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev;

	jmax = 10;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	temp = alloc_2d_array<double>(nb_lev*nb_lev, jmax);
	memset(*temp, 0, nb_lev*nb_lev *jmax *sizeof(double));

	nb_lines = nb_lev*nb_lev;
	tgrid[0] = 0.; // point at 0 K;

	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> li >> lf >> f >> f;
		for (j = 1; j < jmax; j++) {
			input >> temp[i][j];
		}
	}
	input.close();
		
	for (li = 1; li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			i = li *nb_lev + lf;
			f = lf *nb_lev + li;
			nb = li*(li-1)/2 + lf;

			if (li < h2_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *exp((h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]) 
						*h2_di->lev_array[lf].g / ((double) h2_di->lev_array[li].g)); 
				}
			}
		}
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);
	
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Rate coefficients involving level 24 of ortho-H2 (v=1,j) are set to zero, one can set them equal to those involving level 25 (v=3,j=7). 
// Rate coefficients involving level 27 of para-H2 (v=1,j) are set to zero, one can set them equal to those involving level 28 (v=3, j=8).
h2_he_flower_data::h2_he_flower_data(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, li, lf, f, nb, nb_lines;
	double **temp;

	string file_name;
	ifstream input;

	if (rounding(h2_di->mol.spin) == 0) {
		file_name = data_path + "coll_h2/coll_ph2_he.txt";
	}
	else if (rounding(h2_di->mol.spin) == 1) {
		file_name = data_path + "coll_h2/coll_oh2_he.txt";
	}
	else {
		file_name = data_path + "coll_h2/coll_h2_he.txt";
	}

	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev;
	
	jmax = 10;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	temp = alloc_2d_array<double>(nb_lev*nb_lev, jmax);
	memset(*temp, 0, nb_lev*nb_lev *jmax *sizeof(double));

	nb_lines = nb_lev*nb_lev;
	tgrid[0] = 0.; // additional point at 0 K;

	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> li >> lf >> f >> f;
		for (j = 1; j < jmax; j++) {
			input >> temp[i][j];
		}
	}
	input.close();
		
	for (li = 1; li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			i = li *nb_lev + lf;
			f = lf *nb_lev + li;
			nb = li*(li-1)/2 + lf;

			if (li < h2_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *exp((h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]) 
						*h2_di->lev_array[lf].g /((double) h2_di->lev_array[li].g)); 
				}
			}
		}
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);
	
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_h_wrathmall_data::h2_h_wrathmall_data(const std::string &data_path, const energy_diagram *h2_di, bool reactive_channels, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, k, li, lf, f, nb;
	double de, a, **temp;

	string file_name;
	ifstream input;

	if (rounding(h2_di->mol.spin) == 0) {
		file_name = data_path + "coll_h2/coll_ph2_h_w.txt";
	}
	else if (rounding(h2_di->mol.spin) == 1) {
		file_name = data_path + "coll_h2/coll_oh2_h_w.txt";
	}
	else {
		file_name = data_path + "coll_h2/coll_h2_h_w.txt";
	}
	
	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev; // number of levels for which data are given in the file,

	jmax = 61; // point at 0 K;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	temp = alloc_2d_array<double>(nb_lev*nb_lev, jmax);
	memset(*temp, 0, nb_lev*nb_lev *jmax *sizeof(double));

	tgrid[0] = 0.; // point at 0 K;		
	for (j = 1; j < jmax; j++) 
	{
		input >> tgrid[j];		
		for (lf = 0; lf < nb_lev; lf++) {
			for (li = 0; li < nb_lev; li++)
			{
				// columns of the matrix correspond to the initial level of the transition and 
				// rows correspond to the final level of the transition;
				i = li*nb_lev + lf;
				input >> temp[i][j];
				
				if (temp[i][j] < 0.) // some coefficient values in the data are < 0;
					temp[i][j] = 0.;
			}
		}
	}
	input.close();
		
	for (li = 1; li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			i = li *nb_lev + lf;
			f = lf *nb_lev + li;
			nb = (li*(li-1) >> 1) + lf;

			if (li < h2_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *exp((h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]) 
						*h2_di->lev_array[lf].g/ ((double) h2_di->lev_array[li].g)); 
				}
			}
		}
	}
	
	// the contribution of reactive channels to the collisions are added (Le Bourlot et al. 1999):
    if (reactive_channels) {
        for (li = 1; li < h2_di->nb_lev && li < nb_lev; li++) {
            for (lf = 0; lf < li; lf++)
            {
                nb = (li*(li - 1) >> 1) + lf;
                if (h2_di->lev_array[li].v == h2_di->lev_array[lf].v) {
                    if (h2_di->lev_array[li].j - h2_di->lev_array[lf].j == 2) {
                        for (j = 1; j < jmax; j++) {
                            // total rate coefficient was obtained as the sum of the reactive and non-reactive contributions
                            coeff[nb][j] += 8.e-11*exp(-3900. / tgrid[j]); // Schofield, Planet. Space Sci. 15, p.643 (1967);
                        }
                    }
                    else if (h2_di->lev_array[li].j - h2_di->lev_array[lf].j == 1) {
                        if (rounding(h2_di->lev_array[li].j) % 2 == 0) {
                            for (j = 1; j < jmax; j++) {
                                coeff[nb][j] = 8.e-11*exp(-3900. / tgrid[j]);
                            }
                        }
                        else {
                            for (j = 1; j < jmax; j++) {
                                coeff[nb][j] = 2.667e-11*exp(-3900. / tgrid[j]); // reduced by a factor of 3 when J_up is odd
                            }
                        }
                    }
                }
                else {
                    for (j = 1; j < jmax; j++) {
                        de = (3900. - (h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy)*CM_INVERSE_TO_KELVINS) / tgrid[j];
                        if (de < 0.)
                            de = 0.;

                        if (abs(rounding(h2_di->lev_array[li].j - h2_di->lev_array[lf].j)) % 2 == 0) { // |j - j'| is even
                            coeff[nb][j] += coeff[nb][j] * exp(-de);
                        }
                        else { // | j - j'| is odd
                            f = h2_di->get_nb(h2_di->lev_array[lf].v, h2_di->lev_array[lf].j - 1);
                            k = h2_di->get_nb(h2_di->lev_array[lf].v, h2_di->lev_array[lf].j + 1);

                            a = 0.;
                            if (f >= 0 && f < nb_lev) {
                                f = (li*(li - 1) >> 1) + f;
                                a = coeff[f][j];
                            }
                            if (k >= 0 && k < nb_lev) {
                                k = (li*(li - 1) >> 1) + k;
                                a += coeff[k][j];
                            }
                            if (rounding(h2_di->lev_array[li].j) % 2 == 0) // even
                                coeff[nb][j] = 0.5 *a *exp(-de);
                            else // odd
                                coeff[nb][j] = 0.1667 *a *exp(-de);
                        }
                    }
                }
            }
        }
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);
	
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_h2_wan_data::h2_h2_wan_data(const std::string &path, const energy_diagram *h2_di, bool coll_partner_is_ortho, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, u, n, v, nb_lev_f, nb_lines, nb_coll_partners;
	int *lev_index_arr;
	double a;
	
	string file_name, str;
	stringstream ss;
	ifstream input;

	file_name = path + "coll_h2/coll_h2_h2_wan2018.txt";
	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	// comment lines are read:
	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	// the line with molecule name is ignored
	
	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');
	
	// the line with molecule weight is ignored 

	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	ss.clear();
	ss.str(text_line);
	
	ss >> nb_lev_f; // nb of levels, given in the data file
	lev_index_arr = new int [nb_lev_f];
	
	nb_lev = (nb_lev_f < h2_di->nb_lev) ? nb_lev_f : h2_di->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;

	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');
	
	for (i = 0; i < nb_lev_f; i++) {
		ss.clear();
		ss.str(text_line);
		
		ss >> l >> a >> l >> j >> v;
		lev_index_arr[i] = h2_di->get_nb(v, j);

		input.getline(text_line, MAX_TEXT_LINE_WIDTH); // note, the line is read in the cycle end
	}
	
	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	ss.clear();
	ss.str(text_line);
	ss >> nb_lines; // nb of spectroscopic lines
	
	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	for (i = 0; i < nb_lines-1; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}

	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	ss.clear();
	ss.str(text_line);
	ss >> nb_coll_partners; // nb of collisional partners

	for (n = 0; n < nb_coll_partners; n++) {
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		// the line with data details is ignored

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		ss.clear();
		ss.str(text_line);
		ss >> nb_lines; // nb of collisional transitions

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		ss.clear();
		ss.str(text_line);
		
		ss >> jmax; // nb of temperature points
		jmax++; // one point for zero temperature
		
		if (n == 0) {	
			tgrid = new double [jmax];
			tgrid[0] = 0.;

			coeff = alloc_2d_array<double>(imax, jmax);
			memset(*coeff, 0, jmax *imax *sizeof(double));
	
			coeff_deriv = alloc_2d_array<double>(imax, jmax);
			memset(*coeff_deriv, 0, jmax *imax *sizeof(double));
		}

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		ss.clear();
		ss.str(text_line);
		for (j = 1; j < jmax; j++) {	
			ss >> tgrid[j];
		}

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		if ( ((n > 0) && coll_partner_is_ortho) || (n == 0 && !coll_partner_is_ortho) ) {
			for (i = 0; i < nb_lines; i++) {
				ss.clear();
				ss.str(text_line);
				
				ss >> j >> u >> l;
				u = lev_index_arr[u-1]; // the first level has index 0 in the code
				l = lev_index_arr[l-1];

				if (u >= 0 && l >= 0) {
					v = u*(u-1)/2 + l;
					for (j = 1; j < jmax; j++) {
						ss >> coeff[v][j];
					}
				}
				else {
					for (j = 1; j < jmax; j++) {
						ss >> a;
					}
				}
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			}
		}
		else {
			for (i = 0; i < nb_lines; i++) {
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			} // reading all data lines + one next line with comment
		}
	}
	input.close();

	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_h_lique_data::h2_h_lique_data(const string &path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, v, k, l;
    double a;
	string file_name;
	ifstream input;
	
	file_name = path + "coll_h2/coll_h2_h_l.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_lev;
	jmax = 51; // +1 for zero point;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	for (j = 0; j < jmax; j++) {
		tgrid[j] = 100.*j;
	}

	i = 0;
	while (i < imax)
	{
		input >> v >> j;
		k = h2_di->get_nb(v, j);

		input >> v >> j;
		l = h2_di->get_nb(v, j);

        if (k != -1 && l != -1) {
            l = k * (k - 1) / 2 + l; // it is assumed that k > l;
            for (j = 1; j < jmax; j++) {
                input >> coeff[l][j];
            }
        }
        else {
            for (j = 1; j < jmax; j++) {
                input >> a;
            }
        }
		i++;
	}
	input.close();

	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// perhaps, the data must be stitched with the data by Lique(2015)
h2_h_bossion_data::h2_h_bossion_data(const string &path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, k, l, j, n, vi, ji, vf, jf, nb_lines;
	double a;
	string file_name;
    stringstream ss;
	ifstream input;
	
	file_name = path + "coll_h2/coll_h2_h_b.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
    // comments:
    do
        input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    while (text_line[0] == '#');
	
    ss.clear();
    ss.str(text_line);

	ss >> nb_lines >> jmax; // check in the file the presence of these parameters;
	jmax++; // +1 for zero point;

	nb_lev = h2_di->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) { // temperature data are in new line
		input >> tgrid[j];
	}
	
	n = 0;
	while (!input.eof() && n < nb_lines)
	{
		input >> vi >> ji >> vf >> jf;
		l = h2_di->get_nb(vi, ji); // initial level
		k = h2_di->get_nb(vf, jf); // final

		if (l > -1 && k > -1 && l > k) 
		{
			i = l*(l-1)/2 + k;
			for (j = 1; j < jmax; j++) 
			{
				input >> a;
				coeff[i][j] = a;
			}
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> a;
			}
		}
		n++;
	}
	input.close();
	
	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Maximal temperature for which coefficients are calculated is 20000 K;
h2_h_martin_data::h2_h_martin_data(const string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, li, lf, nb_lines;
	double c1, c2, c3, c4, d, t, g;

	string file_name;
	ifstream input;

	nb_lev = h2_di->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 41;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		tgrid[j] = tgrid[j-1] + 500.; // is valid only at high temperatures >> 1000 K;
	}
	
	file_name = data_path + "coll_h2/coll_h2_h_m.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines;

	for (i = 0; i < nb_lines && i < imax; i++)
	{
		input >> li >> lf;
		input >> c1 >> c2 >> c3 >> c4;
		
		if (abs(rounding(h2_di->lev_array[li-1].j - h2_di->lev_array[lf-1].j))%2 == 0)
			d = 1.;
		else d = 0.3;

		for (j = 1; j < jmax; j++) 
		{	
			if (tgrid[j] > 50000.)
				t = d + 50.;
			else t = d + 0.001*tgrid[j];

			g = c1 + c2/t + c3/(t*t) + c4*t;
			coeff[i][j] = pow(10., g);
		}
	}
	input.close();

	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_e_data::h2_e_data(const string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, f, nb_lev_f;
	double a;

	string file_name;
	ifstream input;

	file_name = data_path + "coll_h2/coll_h2_e.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lev_f >> jmax;
	nb_lev = (nb_lev_f < h2_di->nb_lev) ? nb_lev_f : h2_di->nb_lev;

	jmax++; // temperature value at 0 K is added;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	l = 0;
	for (i = 0; i < nb_lev_f; i++) {
		for (f = 0; f < i; f++)
		{
			input >> j >> j >> j >> j;
			if (i < nb_lev && f < nb_lev)
			{
				for (j = 1; j < jmax; j++) {
					input >> coeff[l][j];
				}
				l++;
			}
			else {
				for (j = 1; j < jmax; j++) {
					input >> a;
				}
			}
		}
	}
	input.close();

	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_hp_gonzalez_lezana_data::h2_hp_gonzalez_lezana_data(const std::string &path, const energy_diagram *h2_di, int verbosity)
{
    const int nb_files = 26;
    char text_line[MAX_TEXT_LINE_WIDTH];
    int i, j, n, li, lf, f, vi, vf, ji, jf, t, nb_lines;
    double temp, rate;
    
    string file_name;
    string file_suff[nb_files] = { "Fig3_rateEQM_H3p_v0j1_vf0jf0", "Fig3_rateEQM_H3p_v0j1_vf1jf0", "Fig3_rateEQM_H3p_v0j1_vf2jf0", "Fig3_rateEQM_H3p_v0j1_vf3jf0",
                            "Fig4_rateEQM_H3p_v1j1_vf1jf0", "Fig4_rateEQM_H3p_v2j1_vf2jf0", "Fig4_rateEQM_H3p_v3j1_vf3jf0",
                            "Fig5_rateEQM_H3p_v1j1_vf0jf0", "Fig5_rateEQM_H3p_v1j1_vf2jf0", "Fig5_rateEQM_H3p_v1j1_vf3jf0",
                            "Fig6_rateSQM_H3p_v0j0_vf0jf2", "Fig6_rateSQM_H3p_v0j0_vf0jf3",
                            "Fig7_rateSQM_H3p_v1j0_vf0jf2", "Fig7_rateSQM_H3p_v1j0_vf0jf3", "Fig7_rateSQM_H3p_v1j0_vf1jf2", "Fig7_rateSQM_H3p_v1j0_vf1jf3",
                            "Fig7_rateSQM_H3p_v1j0_vf2jf1", "Fig7_rateSQM_H3p_v1j0_vf2jf2", "Fig7_rateSQM_H3p_v1j0_vf2jf3",  // !
                            "Fig7_rateSQM_H3p_v1j0_vf3jf1", "Fig7_rateSQM_H3p_v1j0_vf3jf2", "Fig7_rateSQM_H3p_v1j0_vf3jf3", 
                            "Fig8_rateSQM_H3p_v2j0_vf2jf2", "Fig8_rateSQM_H3p_v2j0_vf2jf3", 
                            "Fig9_rate_SQM_H3p_v3j0_vf3jf2", "Fig9_rate_SQM_H3p_v3j0_vf3jf3"};
    ifstream input;

    nb_lev = 42; // maximal H2 level, for which data exist, is (v=3,j=3)
    imax = nb_lev * (nb_lev - 1) / 2;
    jmax = 302; // T = 0, 1, 10, 20,..., 3000 K
    
    tgrid = new double[jmax];  
    tgrid[0] = 0.;
    tgrid[1] = 1.;
    for (j = 1; j < jmax-1; j++) {
        tgrid[j+1] = 10.*j;
    }
    
    coeff = alloc_2d_array<double>(imax, jmax);
    memset(*coeff, 0, jmax *imax * sizeof(double));

    coeff_deriv = alloc_2d_array<double>(imax, jmax);
    memset(*coeff_deriv, 0, jmax *imax * sizeof(double));
 
    for (f = 0; f < nb_files; f++) {
        file_name = path + "coll_h2/gonzalez_lezana2017/";
        file_name += file_suff[f];
        file_name += ".dat";
        
        input.open(file_name.c_str(), std::ios_base::in);
        if (!input) {
            cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
            exit(1);
        }
        input.getline(text_line, MAX_TEXT_LINE_WIDTH);
        input.getline(text_line, MAX_TEXT_LINE_WIDTH);
        
        input >> vi >> ji >> vf >> jf;
        input >> nb_lines;

        li = h2_di->get_nb(vi, ji); // initial level
        lf = h2_di->get_nb(vf, jf); // final

        if (li != -1 && lf != -1) {
            if (li > lf) {
                i = li * (li - 1) / 2 + lf;
                for (j = 0; j < nb_lines; j++) {
                    input >> temp >> rate;
                    t = rounding(temp);

                    if (t == 1 || t % 10 == 0) {
                        coeff[i][t / 10 + 1] = rate;
                    }
                }
            }
            else if (li < lf) {
                n = 0;
                i = lf * (lf - 1) / 2 + li;
                
                for (j = 0; j < nb_lines; j++) {
                    input >> temp >> rate;
                    t = rounding(temp);
                   
                    if (t == 1 || t % 10 == 0) {
                        if (rate < 1.e-30) { // some arbitrary small value
                            n = t / 10 + 1;
                        }
                        rate *= exp((h2_di->lev_array[lf].energy - h2_di->lev_array[li].energy) *CM_INVERSE_TO_KELVINS / temp)
                            * h2_di->lev_array[li].g / ((double) h2_di->lev_array[lf].g);
                        coeff[i][t / 10 + 1] = rate;
                    }
                }
                // data from Fig.7, 
                if (li == h2_di->get_nb(1, 0) && 
                    (lf == h2_di->get_nb(2, 1) || lf == h2_di->get_nb(2, 2) || lf == h2_di->get_nb(2, 3))) {
                    n = 300; // n+1 - the last value of the array,
                }
                // be carefull - the rate value may be very low at low temperatures,
                n++;
                for (j = 1; j < n; j++) {
                    coeff[i][j] = coeff[i][n];
                }
            }
        }
        input.close();
    }

    calc_coeff_deriv();
    if (verbosity) {
        cout << "  data have been read for H2-H+ " << endl
            << "  temperature range " << (int)tgrid[1] << " - " << (int)tgrid[jmax - 1] << endl;
    }
}

//
// The class that calculates collisional rates
//

h2_collisions::h2_collisions(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
	: collisional_transitions()
{
	bool coll_partner_is_ortho, reactive_channels(true), ethermal_data(true);
	if (verbosity) {
		cout << "H2 collisional rate coefficients are being initializing..." << endl;
	}
	nb_lev = h2_di->nb_lev;

	coll_data.push_back( new h2_he_flower_data(data_path, h2_di, verbosity) ); // there is new data in the literature

	coll_data.push_back( new h2_h2_wan_data(data_path, h2_di, coll_partner_is_ortho = false, verbosity) );
	coll_data.push_back( new h2_ph2_flower_data(data_path, h2_di, verbosity) );
	coll_data.push_back( new h2_h2_wan_data(data_path, h2_di, coll_partner_is_ortho = true, verbosity) );
	coll_data.push_back( new h2_oh2_flower_data(data_path, h2_di, verbosity) );

	coll_data.push_back( new h2_h_lique_data(data_path, h2_di, verbosity) );
	coll_data.push_back( new h2_h_bossion_data(data_path, h2_di, verbosity) );	
	coll_data.push_back( new h2_h_wrathmall_data(data_path, h2_di, reactive_channels, verbosity) );
	coll_data.push_back( new h2_h_martin_data(data_path, h2_di, verbosity) );

	nb1 = (int) coll_data.size();
    if (ethermal_data) {
        coll_data.push_back( new h2_e_data(data_path, h2_di, verbosity) );
    }
	nb2 = (int) coll_data.size();

	// data on H+ collisions must be here;
    coll_data.push_back( new h2_hp_gonzalez_lezana_data(data_path, h2_di, verbosity));
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void h2_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc, 
    double el_conc, double *&concentration, int *&indices) const
{
	// must be called first:
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, indices);

	concentration[0] = he_conc;
	concentration[1] = concentration[2] = ph2_conc;
	concentration[3] = concentration[4] = oh2_conc;
	concentration[5] = concentration[6] = concentration[7] = concentration[8] = h_conc;  
}

void h2_collisions::set_ion_param(double temp_neutrals, double temp_ions, double hp_conc, double h3p_conc, double *&concentration, int *&indices) const
{
    // H2-H+ data 
    indices[nb2] = coll_data[nb2]->locate( 0.333333*(temp_neutrals + 2.*temp_ions) );
    concentration[nb2] = hp_conc; // +h3p_conc;
}

// The energy of the first level is higher, up_lev.nb > low_lev.nb;
void h2_collisions::get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, 
	double &up_rate, double temp_neutrals, const double *concentration, const int *indices) const
{
	down_rate = 0.;
	if (rounding(low_lev.spin - up_lev.spin) == 0) {
		// collisions with He
		if (up_lev.nb < coll_data[0]->nb_lev) {
			down_rate += coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) 
				*concentration[0];	
		}
		// collisions with p-H2
#if (H2_H2_COLL_DATA == 1)
		if (low_lev.v == 0 && up_lev.v == 0 && up_lev.nb < coll_data[1]->nb_lev) { // Wan et al. (2018)
			down_rate += coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) 
				*concentration[1];	
		}
		else
#endif
        if (up_lev.nb < coll_data[2]->nb_lev) { // Flower & Roueff (1998)
			down_rate += coll_data[2]->get_rate(up_lev.nb, low_lev.nb, indices[2], (temp_neutrals < max_temp[2]) ? temp_neutrals : max_temp[2]) 
				*concentration[2];	
		}
		// collisions with o-H2
#if (H2_H2_COLL_DATA == 1)
		if (low_lev.v == 0 && up_lev.v == 0 && up_lev.nb < coll_data[3]->nb_lev) { // Wan et al. (2018)
			down_rate += coll_data[3]->get_rate(up_lev.nb, low_lev.nb, indices[3], (temp_neutrals < max_temp[3]) ? temp_neutrals : max_temp[3]) 
				*concentration[3];	
		}
		else
#endif
        if (up_lev.nb < coll_data[4]->nb_lev) { // Flower & Roueff (1999)
			down_rate += coll_data[4]->get_rate(up_lev.nb, low_lev.nb, indices[4], (temp_neutrals < max_temp[4]) ? temp_neutrals : max_temp[4]) 
				*concentration[4];	
		}
	}

// is not negligible at small (of the order 10 K) temperatures, provides ortho-/para-H2 conversion noticeable at large evolution times
#if (H2_H_COLL_DATA == 0)
	// H2-H data by Wrathmall et al. (2007):
	if (up_lev.nb < coll_data[7]->nb_lev) {
		down_rate += coll_data[7]->get_rate(up_lev.nb, low_lev.nb, indices[7], (temp_neutrals < max_temp[7]) ? temp_neutrals : max_temp[7]) 
			*concentration[7];
	}
	// H2-H data by Mandy & Martin (1995):
	else if (up_lev.nb < coll_data[8]->nb_lev) {
		down_rate += coll_data[8]->get_rate(up_lev.nb, low_lev.nb, indices[8], (temp_neutrals < max_temp[8]) ? temp_neutrals : max_temp[8]) 
			*concentration[8];
	}
#elif (H2_H_COLL_DATA == 1)
	// H2-H data by Lique (2015) is used:
	if (up_lev.nb < coll_data[5]->nb_lev) {
		down_rate += coll_data[5]->get_rate(up_lev.nb, low_lev.nb, indices[5], (temp_neutrals < max_temp[5]) ? temp_neutrals : max_temp[5]) 
			*concentration[5];
	}
	// H2-H data by Wrathmall et al. (2007):
	else if (up_lev.nb < coll_data[7]->nb_lev) {
		down_rate += coll_data[7]->get_rate(up_lev.nb, low_lev.nb, indices[7], (temp_neutrals < max_temp[7]) ? temp_neutrals : max_temp[7]) 
			*concentration[7];
	}
	// H2-H data by Mandy & Martin (1995):
	else if (up_lev.nb < coll_data[8]->nb_lev) {
		down_rate += coll_data[8]->get_rate(up_lev.nb, low_lev.nb, indices[8], (temp_neutrals < max_temp[8]) ? temp_neutrals : max_temp[8]) 
			*concentration[8];
	}
#elif (H2_H_COLL_DATA >= 2)
    // H2-H data by Lique (2015) is used:
    if (up_lev.nb < coll_data[5]->nb_lev) {
		down_rate += coll_data[5]->get_rate(up_lev.nb, low_lev.nb, indices[5], (temp_neutrals < max_temp[5]) ? temp_neutrals : max_temp[5]) 
			*concentration[5];
	}
	// H2-H data by Bossion et al. MNRAS 480, p.3718, 2018;
	else if (up_lev.nb < coll_data[6]->nb_lev) {
		down_rate += coll_data[6]->get_rate(up_lev.nb, low_lev.nb, indices[6], (temp_neutrals < max_temp[6]) ? temp_neutrals : max_temp[6]) 
			*concentration[6];
	}
#endif

	if (down_rate > MIN_COLLISION_RATE) {
		up_rate = down_rate *exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS/temp_neutrals) *up_lev.g / ((double) low_lev.g);
	}
	else down_rate = up_rate = 0.;
}

void h2_collisions::get_rate_ions(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
    double temp_neutrals, double temp_ions, const double *concentration, const int *indices) const
{
    down_rate = up_rate = 0.;
    if (up_lev.nb < coll_data[nb2]->nb_lev) 
    {
        // t = (t_n*m_i + t_i*m_n)/(m_n + m_i), it is possible that there is no need to such calculations, t_i >> t_n in shock
        temp_ions = 0.333333*(temp_neutrals + 2.*temp_ions);
        down_rate = coll_data[nb2]->get_rate(up_lev.nb, low_lev.nb, indices[nb2], (temp_ions < max_temp[nb2]) ? temp_ions : max_temp[nb2])
            *concentration[nb2];

        if (down_rate > MIN_COLLISION_RATE) {
            up_rate = down_rate * exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS / temp_ions) *up_lev.g / ((double)low_lev.g);
        }
        else down_rate = up_rate = 0.;
    }
}

//
// H2 dissociation
//

h2_h_dissociation_bossion2018::h2_h_dissociation_bossion2018(const std::string & path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, n, vi, ji, nb_lines;
	double a;
    stringstream ss;
	string file_name;
	ifstream input;
	
	file_name = path + "coll_h2/diss_h2_h_b.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}

    do
        input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    while (text_line[0] == '#');
	
    ss.clear();
    ss.str(text_line);

	ss >> nb_lines >> jmax; // check in the file the presence of these parameters;
	jmax++; // +1 for zero point;
	imax = h2_di->nb_lev;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	
	n = 0;
	while (!input.eof() && n < nb_lines)
	{
		input >> vi >> ji;
		i = h2_di->get_nb(vi, ji);
	
		if (i > -1) {
			for (j = 1; j < jmax; j++) {
				input >> coeff[i][j];
			}	
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> a;
			}
		}
		n++;
	}
	input.close();
	
	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_h2_dissociation_martin1998::h2_h2_dissociation_martin1998(const std::string & path, const energy_diagram *h2_di, int verbosity)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int i, j;
    double a;

    string file_name;
    stringstream ss;
    ifstream input;

    file_name = path + "coll_h2/diss_h2-h2_martin1998.txt";
    input.open(file_name.c_str(), std::ios_base::in);

    if (!input) {
        cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
        exit(1);
    }

    do
        input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    while (text_line[0] == '#');

    ss.clear();
    ss.str(text_line);

    ss >> jmax; // check in the file the presence of this parameter;
    jmax++; // +1 for zero point;
    imax = h2_di->nb_lev;

    tgrid = new double[jmax];
    coeff = alloc_2d_array<double>(imax, jmax);
    memset(*coeff, 0, imax*jmax * sizeof(double));

    coeff_deriv = alloc_2d_array<double>(imax, jmax);
    memset(*coeff_deriv, 0, imax*jmax * sizeof(double));

    tgrid[0] = 0.;
    for (j = 1; j < jmax; j++) {
        input >> tgrid[j] >> a;
        for (i = 0; i < imax; i++) { // the rate is identical for all levels
            coeff[i][j] = a;
        }

    }
    input.close();

    calc_coeff_deriv();
    if (verbosity) {
        cout << "  data have been read from file " << file_name << endl
            << "  temperature range " << (int)tgrid[1] << " - " << (int)tgrid[jmax - 1] << endl;
    }
}

h2_h2_dissociation_ceballos2002::h2_h2_dissociation_ceballos2002(const std::string & path, int verbosity)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int i, j;
  
    string file_name;
    ifstream input;

    min_vibrq = 5; // min <= v <= max
    max_vibrq = 14;

    file_name = path + "coll_h2/diss_h2-h2_ceballos2002.txt";
    input.open(file_name.c_str(), std::ios_base::in);

    if (!input) {
        cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
        exit(1);
    }
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);

    input >> imax >> jmax; // check in the file the presence of these parameters;
    jmax++; // +1 for zero point;
    
    tgrid = new double[jmax];
    coeff = alloc_2d_array<double>(imax, jmax);
    memset(*coeff, 0, imax*jmax * sizeof(double));

    coeff_deriv = alloc_2d_array<double>(imax, jmax);
    memset(*coeff_deriv, 0, imax*jmax * sizeof(double));

    tgrid[0] = 0.;
    for (j = 1; j < jmax; j++) {
        input >> tgrid[j];
    }

    for (i = 0; i < imax; i++) {
        input >> j >> j;
        for (j = 1; j < jmax; j++) {
            input >> coeff[i][j];
        }
    }
    input.close();

    calc_coeff_deriv();
    if (verbosity) {
        cout << "  data have been read from file " << file_name << endl
            << "  temperature range " << (int)tgrid[1] << " - " << (int)tgrid[jmax - 1] << endl;
    }
}

double h2_h2_dissociation_ceballos2002::get_rate(int v, double temp, double *vibr_h2_conc) const
{
    if (v < min_vibrq || v > max_vibrq) 
        return 0.;

    // upper limit on temperature
    if (temp > 1.5*tgrid[jmax - 1])
        temp = 1.5*tgrid[jmax - 1];

    double rate(0.);
    int j, l = 0, r = jmax - 1;
    while (r - l > 1)
    {
        j = l + ((r - l) >> 1);
        if (tgrid[j] < temp)
            l = j;
        else r = j;
    }
    
    v = ((v - min_vibrq) / 2) * 5; // 
    for (j = min_vibrq, r = 0; j < nb_vibr_states_h2 && j < max_vibrq; j += 2, r++) 
    {
        rate += (coeff[v + r][l] + coeff_deriv[v + r][l] * (temp - tgrid[l])) 
            * (vibr_h2_conc[j] + vibr_h2_conc[j+1]);
    }
    return 0.5*rate; // 0.5 due to statistical considerations
}


h2_e_dissociation_tennyson::h2_e_dissociation_tennyson()
{;}

void h2_e_dissociation_tennyson::get_rate(double temp_e, double *vibr_h2_rates) const
{
    vibr_h2_rates[0] = 4.4886e-9 * pow(temp_e, 0.1091) *exp(-101858.0 / temp_e);
    vibr_h2_rates[1] = 5.6559e-9 * pow(temp_e, 0.0694) *exp(-85371.7 / temp_e);
    vibr_h2_rates[2] = 2.1755e-9 * pow(temp_e, 0.1491) *exp(-71582.2 / temp_e);
    vibr_h2_rates[3] = 0.9581e-9 * pow(temp_e, 0.2171) *exp(-60136.4 / temp_e);
    vibr_h2_rates[4] = 0.4223e-9 * pow(temp_e, 0.2857) *exp(-50186.9 / temp_e);
    // what about v > 4?
}

//
// H2 excitation processes
//

h2_grain_formation::h2_grain_formation(const energy_diagram *h2_di)
{
	int i;
	double sum = 0.;

	nb_lev = h2_di->nb_lev;
	popul = new double [nb_lev];
    // dissociation energy of H2 - 36118.11 cm-1 = 4.4781 eV  Herzberg & Monfils, J. Molecular Spectroscopy 5, no.1–6, p.482-498 (1961) 
    average_energy = 0.3333 *4.4781 *EV_TO_CM_INVERSE; 

	for (i = 0; i < nb_lev; i++) 
	{
		popul[i] = h2_di->lev_array[i].g *exp(-h2_di->lev_array[i].energy/average_energy);
		sum += popul[i];
	}
	for (i = 0; i < nb_lev; i++) {
		popul[i] /= sum;
	}
}

h2_grain_formation::~h2_grain_formation() {
	delete [] popul;
}

h2_gasphase_formation::h2_gasphase_formation(const energy_diagram *h2_di)
{
	int i, j;
	double gtemp, sum = 0.;

	max_temp = 20000.; // in K
	step = 100.; // in K
	step_inv = 1./step;

	nb_lev = h2_di->nb_lev;
	nb_temp = (int) (max_temp/step) + 1;
	
	popul = alloc_2d_array<double>(nb_lev, nb_temp);
	memset(*popul, 0, nb_lev*nb_temp*sizeof(double));
	
	popul[0][0] = 1.; // zero temperature distribution;
	for (j = 1; j < nb_temp; j++) {
		gtemp = step *j/CM_INVERSE_TO_KELVINS; // in cm-1
		for (i = 0; i < nb_lev; i++) 
		{
			popul[i][j] = h2_di->lev_array[i].g *exp(-h2_di->lev_array[i].energy/gtemp);
			sum += popul[i][j];
		}
		for (i = 0; i < nb_lev; i++) {
			popul[i][j] /= sum;
		}
	}
}

double h2_gasphase_formation::get_efficiency(int l, double gtemp) const
{
	int j;
	double x;

	if (gtemp > max_temp) 
		gtemp = max_temp;
	
	x = step_inv*gtemp;
	j = (int) x;
	x = x - j;
	return (1 - x)*popul[l][j] + x*popul[l][j+1];
}

h2_gasphase_formation::~h2_gasphase_formation() {
	free_2d_array(popul);
}

h2_excit_cosmic_rays::h2_excit_cosmic_rays(const std::string &data_path, const energy_diagram *h2_di, int verbosity) 
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, v, li, lf, dj, imax, jmax;
	double a;

	string file_name = data_path + "xray_efficiencies.txt";
	ifstream input(file_name.c_str());
	
	nb_ion = 4;
	ion_grid = new double [nb_ion];
	
	nb_init_lev = h2_di->get_nb(0, 9.) + 1; // be carefull here, the levels (v,j) = (1,0)-(1,3) are below (0,9);
	nb_fin_lev = h2_di->nb_lev;
	
	eff = alloc_2d_array<double>(nb_ion, nb_fin_lev*nb_init_lev);
	memset(*eff, 0, nb_ion*nb_fin_lev*nb_init_lev*sizeof(double));

	exit_eff = alloc_2d_array<double>(nb_ion, nb_init_lev);
	memset(*exit_eff, 0, nb_ion*nb_init_lev*sizeof(double));
	
	if (input) 
	{
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> imax >> jmax;

		for (l = 0; l < nb_ion; l++)
		{
			input >> ion_grid[l];
			for (i = 0; i < imax; i++)
			{
				input >> v >> dj;
				for (j = 0; j < jmax; j++) 
				{
					input >> a;		
					li = h2_di->get_nb(0, (double) j);
					lf = h2_di->get_nb(v, (double) (j+dj));

					if (li != -1 && lf != -1)
					{
						eff[l][li*nb_fin_lev + lf] = a;
						exit_eff[l][li] += a; // re-entry is not excluded
					}
				}
			}
		}
	}
	else {
		cout << "Error in " << SOURCE_NAME << ": can't open file " << file_name << endl;
		exit(1);
	}
	input.close();

	if (verbosity) 
		cout << "  data have been read from file " << file_name << endl;
}

h2_excit_cosmic_rays::~h2_excit_cosmic_rays()
{
	delete [] ion_grid;
	free_2d_array<double>(eff);
	free_2d_array<double>(exit_eff);
}

// The maximal value of the index is nb_ion-2; the minimal value is -1;
// log-normal approximation is used;
void h2_excit_cosmic_rays::set_ionization(double ionization, int & index, double & param) const
{ 
	int i;
	for (i = 0; i < nb_ion; i++) {
		if (ionization < ion_grid[i])
			break;
	}
    if (i == nb_ion) 
        i--;
	i--;

	if (i >= 0)
		param = log(ionization/ion_grid[i]) /log(ion_grid[i+1]/ion_grid[i]);
	else param = 0.;
	index = i;
}

double h2_excit_cosmic_rays::get_efficiency(int il, int fl, int index, double param) const
{
	if (il > nb_init_lev-1) return 0.;
	double answ;

	if (index < 0)
		answ = eff[0][il*nb_fin_lev + fl];
    else if (index < nb_ion - 1) {
        int i = il * nb_fin_lev + fl;
        answ = eff[index][i] + param * (eff[index + 1][i] - eff[index][i]);
    }
	else // the case of large ionization fractions, > 0.01, must be considered accurately;
		answ = 0.;
	
	return answ;
}

double h2_excit_cosmic_rays::get_efficiency(int l, int index, double param) const
{
	if (l > nb_init_lev-1) return 0.;
	
	double answ;
	if (index < 0)
		answ = exit_eff[0][l];
	else if (index < nb_ion-1)
		answ = exit_eff[index][l] + param*(exit_eff[index+1][l] - exit_eff[index][l]);
	else // the case of large ionization fractions, > 0.01, must be considered accurately;
		answ = 0.;

	return answ;
}

//
// Functions
//

void h2_coll_data_process(const std::string &data_path)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, f, k, l, m, jmax, nb_lev_h2, nb_lev_ph2, nb_lev_oh2;
	double *tgrid, **temp;

	string fnp, fno, fn;
	string n1_arr[3] = {"coll_ph2_oh2.txt", "coll_ph2_ph2.txt", "coll_ph2_he.txt"};
	string n2_arr[3] = {"coll_oh2_oh2.txt", "coll_oh2_ph2.txt", "coll_oh2_he.txt"};
	string n3_arr[3] = {"coll_h2_oh2.txt", "coll_h2_ph2.txt", "coll_h2_he.txt"};
	
	ifstream input;
	ofstream out;
	nb_lev_h2 = nb_lev_ph2 = nb_lev_oh2 = 320;

	molecule ph2_mol("ph2", 1, 2.*ATOMIC_MASS_UNIT, 0);
	molecule oh2_mol("oh2", 1, 2.*ATOMIC_MASS_UNIT, 1);
	molecule h2_mol("h2", 1, 2.*ATOMIC_MASS_UNIT);
	 
	h2_diagram *ph2_di
		= new h2_diagram(data_path, ph2_mol, nb_lev_ph2); // nb of levels is redefined here

	h2_diagram *oh2_di 
		= new h2_diagram(data_path, oh2_mol, nb_lev_oh2);
	
	h2_diagram *h2_di 
		= new h2_diagram(data_path, h2_mol, nb_lev_h2);

	// Nb of temperatures for the data in question:
	jmax = 9;
	tgrid = new double [jmax];
	temp = alloc_2d_array<double>(nb_lev_h2*nb_lev_h2, jmax);
	
	for (m = 0; m < 3; m++) 
	{
		memset(*temp, 0, nb_lev_h2*nb_lev_h2*jmax *sizeof(double));

		fnp = data_path + "coll_h2/";
		fnp += n1_arr[m];
	
		fno = data_path + "coll_h2/";
		fno += n2_arr[m];

		input.open(fnp.c_str(), std::ios_base::in);
		input.seekg (0, ios_base::beg);

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb_lev_ph2;

		for (j = 0; j < jmax; j++) {
			input >> tgrid[j];
		}
		for (l = 0; l < nb_lev_ph2*nb_lev_ph2; l++)
		{
			input >> i >> f >> j >> j;
			// Calculation of the level numbers in the total H2 level list:
			i = h2_di->get_nb(ph2_di->lev_array[i-1].v, ph2_di->lev_array[i-1].j);
			f = h2_di->get_nb(ph2_di->lev_array[f-1].v, ph2_di->lev_array[f-1].j);

			for (j = 0; j < jmax; j++) {	
				input >> temp[i*nb_lev_h2+f][j];
			}
		}
		input.close();

		input.open(fno.c_str(), std::ios_base::in);
		input.seekg (0, ios_base::beg);

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb_lev_oh2;

		for (j = 0; j < jmax; j++) {
			input >> tgrid[j];
		}
		for (l = 0; l < nb_lev_oh2*nb_lev_oh2; l++)
		{
			input >> i >> f >> j >> j;
			i = h2_di->get_nb(oh2_di->lev_array[i-1].v, oh2_di->lev_array[i-1].j);
			f = h2_di->get_nb(oh2_di->lev_array[f-1].v, oh2_di->lev_array[f-1].j);

			for (j = 0; j < jmax; j++) {	
				input >> temp[i*nb_lev_h2+f][j];
			}
		}
		input.close();

		// Calculation of the number of the highest level, for which rate coefficients are available:
		i = h2_di->get_nb(ph2_di->lev_array[nb_lev_ph2-1].v, ph2_di->lev_array[nb_lev_ph2-1].j);
		j = h2_di->get_nb(oh2_di->lev_array[nb_lev_oh2-1].v, oh2_di->lev_array[nb_lev_oh2-1].j);
		k = (i > j ) ? i+1 : j+1;
	
		fn = data_path + "coll_h2/";
		fn += n3_arr[m];

		out.open(fn.c_str());
		out << "# joined data for o-H2 and p-H2;" << endl << k << endl;

		for (j = 0; j < jmax; j++) {
			out << left << setw(12) << tgrid[j];
		}
		out << endl;
	
		for (i = 0; i < k; i++) {
			for (f = 0; f < k; f++) 
			{
				out << left << setw(5) << i+1 << setw(5) << f+1 << setw(5) << "1" << setw(5) << "1";
				for (j = 0; j < jmax; j++) {
					 out << left << setw(12) << temp[i*nb_lev_h2+f][j];
				}
				if (i < k-1 || f < k-1) 
					out << endl;
			}
		}
		out.close();
	}

	delete [] tgrid;
	free_2d_array(temp);

	jmax = 60;
	tgrid = new double [jmax];

	temp = alloc_2d_array<double>(nb_lev_h2*nb_lev_h2, jmax);
	memset(*temp, 0, nb_lev_h2*nb_lev_h2*jmax *sizeof(double));

	fnp = data_path + "coll_h2/coll_ph2_h_w.txt";
	fno = data_path + "coll_h2/coll_oh2_h_w.txt";
	
	input.open(fnp.c_str(), std::ios_base::in);
	input.seekg(0, ios_base::beg);

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev_ph2;

	for (j = 0; j < jmax; j++) 
	{
		input >> tgrid[j];
		for (f = 0; f < nb_lev_ph2; f++) {
			for (i = 0; i < nb_lev_ph2; i++)
			{
				l = h2_di->get_nb(ph2_di->lev_array[i].v, ph2_di->lev_array[i].j); // initial level
				m = h2_di->get_nb(ph2_di->lev_array[f].v, ph2_di->lev_array[f].j);
				input >> temp[l*nb_lev_h2+m][j];
			}
		}
	}
	input.close();

	input.open(fno.c_str(), std::ios_base::in);
	input.seekg (0, ios_base::beg);

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev_oh2;

	for (j = 0; j < jmax; j++) 
	{
		input >> tgrid[j];
		for (f = 0; f < nb_lev_oh2; f++) {
			for (i = 0; i < nb_lev_oh2; i++)
			{
				l = h2_di->get_nb(oh2_di->lev_array[i].v, oh2_di->lev_array[i].j);
				m = h2_di->get_nb(oh2_di->lev_array[f].v, oh2_di->lev_array[f].j);
				input >> temp[l*nb_lev_h2+m][j];
			}
		}
	}
	input.close();

	i = h2_di->get_nb(ph2_di->lev_array[nb_lev_ph2-1].v, ph2_di->lev_array[nb_lev_ph2-1].j);
	j = h2_di->get_nb(oh2_di->lev_array[nb_lev_oh2-1].v, oh2_di->lev_array[nb_lev_oh2-1].j);
	k = (i > j ) ? i+1 : j+1; // 
	
	fn = data_path + "coll_h2/coll_h2_h_w.txt";	
	out.open(fn.c_str());
	out.precision(4);

	out << "# S.A. Wrathmall et al., Mon. Not. R. Astron. Soc. 382, 133 (2007); http://ccp7.dur.ac.uk/pubs.html" << endl
		<< "# columns correspond to the initial level, rows - to the final level; joined data for o-H2 and p-H2;" << endl 
        << k << endl;

	for (j = 0; j < jmax; j++) 
	{
		out << tgrid[j] << endl;
		for (f = 0; f < k; f++) {
			for (i = 0; i < k; i++) {
				out << left << setw(12) << temp[i*nb_lev_h2+f][j];
			}
			if (f < k-1 || j < jmax-1) 
				out << endl;
		}
	}
	out.close();

	delete ph2_di;
	delete oh2_di;
	delete h2_di;

	delete [] tgrid;
	free_2d_array(temp);
}

void h2_coll_data_process_bossion(const std::string &data_path)
{ 
    const int max_temp = 20000; // maximal temperature in the output file
	int k, i, li, lf, vi, vf, ji, jf, nbt, nb_lev, nb_lines;
	double r;
    double **rate_coeff(0), **diss_coeff(0);

	string fname;
	ifstream input;
	ofstream output;

	nb_lev = 298; // the highest level for which data exists is v,j=12,10
	molecule h2_mol("H2", 1, 2.*ATOMIC_MASS_UNIT);
	
	h2_diagram *h2_di 
		= new h2_diagram(data_path, h2_mol, nb_lev);
  
    nbt = 200;
    nb_lines = nb_lev*(nb_lev - 1)/2;

    rate_coeff = alloc_2d_array<double>(nb_lines, nbt);
    memset(*rate_coeff, 0, nb_lines*nbt*sizeof(double));

	// please, check the data format, presence of the comments, nb of temperature values;
	fname = data_path + "coll_h2/H_H2_StS_rates_100-10000.txt";
	input.open(fname.c_str());

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
 
	while (!input.eof())
	{
		input >> vi >> ji >> vf >> jf;
		if (input.eof())
			break;
        
        li = h2_di->get_nb(vi, ji);
        lf = h2_di->get_nb(vf, jf);

		if  (li != -1 && lf != -1 && li > lf) {			
            i = (li*(li - 1) >> 1) + lf;

			for (k = 0; k < 100; k++) {
                input >> r; 
				r = (r > 0.) ? r : 0.; 
                rate_coeff[i][k] = r;
			}
		}
		else {
			for (k = 0; k < 100; k++) {
				input >> r;
			}
		}
	}
	input.close();

    fname = data_path + "coll_h2/H_H2_StS_rates_10100-20000.txt";
    input.open(fname.c_str());

    if (!input) {
        cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
        exit(1);
    }

    while (!input.eof())
    {
        input >> vi >> ji >> vf >> jf;
        if (input.eof())
            break;
        
        li = h2_di->get_nb(vi, ji);
        lf = h2_di->get_nb(vf, jf);

        if (li != -1 && lf != -1 && li > lf) {
            i = (li*(li - 1) >> 1) + lf;

            for (k = 100; k < 200; k++) {
                input >> r;
                r = (r > 0.) ? r : 0.;
                rate_coeff[i][k] = r;
            }
        }
        else {
            for (k = 0; k < 100; k++) {
                input >> r;
            }
        }
    }
    input.close();

    // output file
    fname = data_path + "coll_h2/coll_h2_h_b_.txt";
    output.open(fname.c_str());

    output << scientific;
    output.precision(2);
    output << "# The table contains the H2-H collisional rate coefficients. Bossion et al. MNRAS 480, p.3718, 2018;" << endl
        << "# the highest level for which data exists is (v,j) = (12,10), number 298;" << endl
        << "# vi ji -> vf jf, k(cm3 s-1) (T); nb of data lines, nb of temperatures (T in K):" << endl;

    if (max_temp/100 < nbt)
        nbt = max_temp/100;

    output << left << setw(8) << nb_lines << setw(8) << nbt << endl << left << setw(16) << "";
    for (k = 1; k <= nbt; k++) {
        output << left << setw(10) << k*100;
    }
    
    for (li = 1; li < nb_lev; li++) {
        for (lf = 0; lf < li; lf++) {
            output << endl << left << setw(4) << h2_di->lev_array[li].v << setw(4) << rounding(h2_di->lev_array[li].j) 
                << setw(4) << h2_di->lev_array[lf].v << setw(4) << rounding(h2_di->lev_array[lf].j);
            
            i = (li*(li - 1) >> 1) + lf;
            for (k = 0; k < nbt; k++) {
                output << left << setw(10) << rate_coeff[i][k];
            }
        }
    }
	output.close();

	// please, check the data format, presence of the comments, nb of temperature values;
	nbt = 200;
    diss_coeff = alloc_2d_array<double>(nb_lev, nbt);
    memset(*diss_coeff, 0, nb_lev*nbt * sizeof(double));

	fname = data_path + "coll_h2/H_H2_diss_rates_100-10000.txt";
	input.open(fname.c_str());

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with dissociation data " << fname << endl;
		exit(1);
	}

	while (!input.eof())
	{
		input >> vi >> ji;
		if (input.eof())
			break;
        
        li = h2_di->get_nb(vi, ji);
        if (li != -1) {
            for (k = 0; k < 100; k++) {
                input >> r;
                r = (r > 0.) ? r : 0.;
                diss_coeff[li][k] = r;
            }
        }
        else {
            for (k = 0; k < 100; k++) {
                input >> r;
            }
        }
	}
	input.close();
    
    fname = data_path + "coll_h2/H_H2_diss_rates_10100-20000.txt";
    input.open(fname.c_str());

    if (!input) {
        cout << "Error in " << SOURCE_NAME << ": can't open file with dissociation data " << fname << endl;
        exit(1);
    }

    while (!input.eof())
    {
        input >> vi >> ji;
        if (input.eof())
            break;

        li = h2_di->get_nb(vi, ji);
        if (li != -1) {
            for (k = 100; k < 200; k++) {
                input >> r;
                r = (r > 0.) ? r : 0.;
                diss_coeff[li][k] = r;
            }
        }
        else {
            for (k = 0; k < 100; k++) {
                input >> r;
            }
        }
    }
    input.close();

    fname = data_path + "coll_h2/diss_h2_h_b_.txt";
	output.open(fname.c_str());
	
	output << scientific;
	output.precision(3);
	output << "# The table contains the H2+H->H+H+H dissociation rate coefficients. Bossion et al. MNRAS 480, p.3718, 2018;" << endl
        << "# levels are included that present in Dabrowski I., Can. J. Phys. 62, p. 1639, 1984" << endl 
        << "# the highest level for which data exists is (v,j) = (12,10), number 298;" << endl
		<< "# vi ji, k(cm3 s-1) (T), temperature in K; " << endl;

	output << left << setw(11) << nb_lev << setw(11) << nbt << endl << left << setw(8) << " ";
	for (k = 1; k <= nbt; k++) {
		output << left << setw(11) << k*100;
	}

    for (li = 0; li < nb_lev; li++) {
        output << endl << left << setw(4) << h2_di->lev_array[li].v << setw(4) << rounding(h2_di->lev_array[li].j);
        for (k = 0; k < nbt; k++) {
            output << left << setw(11) << diss_coeff[li][k];
        }
    }
	output.close();

    free_2d_array(diss_coeff);
    free_2d_array(rate_coeff);
}
