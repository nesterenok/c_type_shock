#pragma once
#include "spectroscopy.h"
#include "coll_rates.h"

// General-purpose class for ion-electron collisions; the collision strength data are given in the file;
// for OI - e: Bell, Berrington and Thomas, MNRAS 293, p. L83 (1998); 50 < T < 3000 K;
class ion_electron_coll_data : public collision_data
{
public:
	ion_electron_coll_data(const std::string fname, const energy_diagram *, int verbosity=1);
};

// General-purpose class for ion-neutral collisions; the rate coefficient data are given in the file;
// for OI - he: Monteiro & Flower, MNRAS 228, pp. 101-107, 1987; 100 < T < 1000 K;
// for CI - he: Staemmler & Flower, J. Phys. B 24, p. 2343 (1991); 10 < T < 150 K; 
// for CI - h2: Schroder et al., J. Phys. B 14, p. 2487 (1991); 10 < T < 1200 K;
// for CII - h: Barinovs et al., ApJ 620, p. 537 (2005); 20 < T < 2000 K;
class ion_neutral_coll_data : public collision_data
{
public:
	ion_neutral_coll_data(const std::string fname, int verbosity=1);
};

//
// Data on collisional rates for OI atom;
//
// Approximation: Draine, "Interstellar medium" (2011); Pequignot, A&A 231, p. 499 (1990);
// data: Berrington & Burke, Planet Space Science 29, p. 377 (1981); Berrington, J. Phys. B 21, p. 1083 (1988); 50 < T < 1000 K;
class OI_electrons_berrington_data : public collision_data
{
public:
	OI_electrons_berrington_data(int verbosity=1);
};

// Abrahamsson & Krems, ApJ 654, p. 1171, (2007); for three lower levels; 30 < T < 1000 K;
// Krems, Jamieson, Dalgarno, ApJ 647, p. 1531 (2006); for transitions including fourth level;
class OI_h_abrahamsson_data : public collision_data
{
public:
	OI_h_abrahamsson_data(int verbosity=1);
};

// There are several sources of data;
// Data on OI - H2 collisions are given by Jaquet et al., J.of Phys. B 25, p. 285 (1992); 20 < T < 1500 K;
// The approximations are given by Glover & Jappsen, ApJ 666, p. 1 (2007); and by Draine, "Interstellar medium" (2011);
class OI_oh2_data : public collision_data
{
public:
	OI_oh2_data(int verbosity=1);
};

// See the comments above;
class OI_ph2_data : public collision_data
{
public:
	OI_ph2_data(int verbosity=1);
};

// Methods of this class calculate rates of collisional transitions of OI;
class OI_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc, 
		double el_conc, double *&concentration, int *&indices) const;
	
	OI_collisions(const std::string &, const energy_diagram *, int verbosity =1);
};

//
// Data on collisional rates for CI atom;
//
// Johnson, Burke & Kingston, J. Phys. B 20, p. 2553 (1987); 7.5 < T < 10000 K;
class CI_e_johnson_data : public collision_data
{
public:
	CI_e_johnson_data(int verbosity=1);
};

// Abrahamsson & Krems, ApJ 654, p. 1171, (2007); 5 < T < 1000 K;
class CI_h_abrahamsson_data : public collision_data
{
public:
	CI_h_abrahamsson_data(int verbosity=1);
};

// Methods of this class calculate rates of collisional transitions of CI;
class CI_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc, 
		double el_conc, double *&concentration, int *&indices) const;
	
	CI_collisions(const std::string &, const energy_diagram *, int verbosity =1);
};

//
// Data on collisional rates for CII atom;
//
// For lowest two transitions: Tayal, A&A 486, pp. 629–636 (2008); fit by Draine, "Interstellar medium" (2011); 1000 < T < 300000 K;
// Only data at high electron temperature are available >1000 K, we assume that collisional strength is constant at lower energies;
class CII_e_tayal2008_data : public collision_data
{
public:
	CII_e_tayal2008_data(int verbosity=1);
};

// Lique et al., J. of Chem. Phys. 138, 204314 (2013); Wiesenfeld & Goldsmith ApJ 780, p.183 (2014)
class CII_h2_lique2013 : public collision_data
{
public:
	CII_h2_lique2013(bool is_para_h2, int verbosity = 1);
};

// Flower & Launay, J. Phys. B 10, p. 3673 (1977); cross section data and cooling rate coefficients were given;
// approximation formulae of collisional rate coefficient by Draine, "Interstellar medium" (2011); 0 < T < 250 K;
class CII_h2_flower1977 : public collision_data
{
public:
	CII_h2_flower1977(bool is_para_h2, int verbosity = 1);
};

// Methods of this class calculate rates of collisional transitions of CII;
class CII_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
		double el_conc, double *&concentration, int *&indices) const;

	CII_collisions(const std::string &, const energy_diagram *, int verbosity =1);
};
