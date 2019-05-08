#pragma once

//---------------------------------------------
// Parameters related to chemistry
//---------------------------------------------

// Helium abundance, this parameter must be equal to the value given in the file with element abundances,
#define HE_TO_H_NB_RATIO 0.09

// The standard CR ionization rate is linked to the chemistry database UMIST 2012 (McElroy ET AL., A&A 550, A36, 2013, see UMIST2012 data file);
// Draine ("Interstellar medium", 2010, p.136) points out the characteristic value of ionization rate of 1.e-16 s-1;
// see also Dalgarno, PNAS 103, pp. 12269-12273 (5.e-17 - 1.e-16 s-1);
#define STANDARD_CR_IONIZ_RATE 1.3e-17  // s-1 per H2

// Standard flux of CR iron nuclei;
// Roberts et al., MNRAS 382, 733�742 (2007) used value 2.e-3 cm-2 s-1;
// Fe/H = 4.e-4 (> 5 GeV), Grieder, Cosmic Rays at Earth (2005);
// Fe/H = 1.6e-4 (> 0.02 GeV/n), Leger et al. A&A 144, 147, 1985);
// Fe/H = 7.e-4, Meyer et al., Space Sci. Rev. 86, 179 (1998);
// LIS proton flux is about 18 cm-2 s-1 (Webber, ApJ 506, 329, 1998; Webber, Yushak, ApJ 275, 391, 1983) that corresponds to CR ionization of 2.e-17 s-1 (Draine 2011);
// our estimate of iron flux for CR ionization rate 1.3e-17 s-1 and for Fe/H = 1.6e-4:
#define STANDARD_FLUX_CR_IRON 2.e-3 // cm-2 s-1

// UV photons cm-2 s-1, generated by CR particles;
// Cecchi-Pestellini, Aiello, MNRAS 258, p. 125 (1992); see also the work by Gredel et al., ApJ 347, p. 289, 1989;
// given value corresponds to N_H2/A_v = 0.9*10^21 cm-2 mag-1, Rv = 5.5, CR ionization rate 1.36e-17 s-1; dust abedo 0.5 is taken into account;
#define STANDARD_NB_CR_PHOTONS 3000.
// Comments: Flower & Pineau des Forets, MNRAS 343, p. 390 (2003):
// nb of photons generated per unit volume and time = 0.15*n_H*cr_ion_rate_perH [ph cm-3 s-1] - Cecchi-Pestellini, Aiello, MNRAS 258, p. 125 (1992);
// mean grain opacity = (0.5-1)*1.e-21 n_H cm-1, at UV wavelength (10^4-10^5 cm-1), Rv = 5.5
// source function = emissivity/opacity = 2000-4000 cm-2 s-1 (angle- and frequency-integrated)
// Dalgarno et al. (ApJSS 125, p.237, 1999):
// nb of excitations per ion pair at high energy limit in H2-He gas, He/H = 0.1: B1S+u - 0.35, C1Pu - 0.3
// Kalvans, MNRAS 478, p.2753 (2018) used the value 4875 cm-2 s-1 at cosmic ray i.r. 1.3e-17 s-1; 

// see Draine, 2011, p.123;
// Draine, ApJS 36, p.595 (1978); F=1.94e+8 is integrated by angle and energy over 6-13.6 eV (given approximation is over 5-13.6 eV);
#define DRAINE1978_ISRF_FUV 2.e+8  // ph cm-2 s-1, 
// Mathis et al., A&A 128, p. 212 (1983); F=1.3e+7, E=6-13.6 eV
#define MATHIS1983_ISRF_FUV 1.3e+8  // ph cm-2 s-1, 

// alternative approach; Roberts et al., MNRAS 382, 733�742, 2007; 
// it is assumed that only the volatile species (CO, N2, NO, O2, C2 and CH4) would be desorbed during the transient heating, 
// check limiting binding energy:
#define CR_DESORPTION_LIM_BENERGY 1251. // in K;
// Roberts et al. (2007) used 1.e+5 as a standard value; the yield is estimated based on Bringa, Johnson, ApJ 603, 159, 2004; 
#define CR_DESORPTION_YIELD 1.e+5

// Hasegawa & Herbst, ApJSS 82, p. 167-195 (1992) N = 1.5e+15 sites per cm2;
// Cuppen et al., Space Sci. Rev. 212, 1 (2017) N = 1.e+15 sites per cm2;
// In some experiments and calculations - even lower (Congiu et al., MNRAS 397, L96, 2009)
// (depends on surface material)
#define GRAIN_SITES_PER_CM2 1.e+15

// characteristic vibration frequency constant for the adsorbed species, Reboussin et al., MNRAS 440, 3557 (2014);
// constant must be multiplied by sqrt(binding energy in K/ mass in a.m.u.);
#define SURF_VIBR_FREQ 4104.71*sqrt(GRAIN_SITES_PER_CM2)

// Garrod et al. A&A 457, p.927 (2006) introduces value of 0.5; Reboussin et al. MNRAS 440, p.3557 (2014) also uses this value;
// Karssemeijer & Cuppen, A&A 569, p. 107 (2014) obtained values 0.3 for CO and 0.4 for CO2; 
#define DIFFUSION_TO_BINDING_ENERGY_RATIO 0.35
// Garrod et al., A&A 467, p. 1103 (2007), description of chemical desorption mechanism;
#define CHEM_DESORPTION_FACTOR 0.01

// see discussion of this parameter by Garrod, Pauly, ApJ 735, p. 15 (2011), they suggested to use 2 A;
// Penteado et al., ApJ 844, p. 71 (2017) used value 1.5 A for diffusion barrier thickness;
#define CHEM_BARRIER_THICKNESS 1.5e-8 // in cm
#define DIFF_BARRIER_THICKNESS 1.5e-8

// see discussion by Ruaud et al., MNRAS 459, p. 3756 (2016); Bertin et al., ApJ 779, p. 120 (2013) 
// - photodesorption is indirect process, individual values are not appropriate, in molecules per photon:
#define COMMON_PHOTODES_YIELD 1.e-4

// Maximal value of the desorption barrier is set to be equal to that of water - the evaporation of water, the primary constituent 
// of icy grain mantles, should result in the codesorption of other species (Garrod et al., ApJ 682, p. 283, 2008). 
// Diffusion barriers are unaffected by this adjustment;
#define SET_MAX_THERMAL_EVAPOR_BARRIER 4800. // in K;

// 230 K - Ruaud et al., MNRAS 459, p. 3756 (2016),
// 255 K - Wakelam et al., Molecular Astrophysics 9, 1 (2017);
// if one uses data on binding energy 65O K by Penteado et al., ApJ 844, p.71 (2017), difussion to binding ratio 0.4, this patameter is equal to 260 K;
// #define H_ATOM_DIFF_BARRIER 255. // in K 

// Cuppen, Herbst, ApJ 668, p. 294, 2007; Hincelin et al. A&A, 2014, arXiv:1410.7375v2;
#define H2_ON_H2_BINDING_ENERGY 23. // in K

// Kalvans, MNRAS 478, p. 2753 (2018), Oberg, Chemical Reviews 116, p.9631, 2016;
// gas-phase photodissociation rate coefficients have to be reduced by a factor of 0.2-0.4 before applying them to icy species;
// additional factor 0.5 is attributed to grain shielding;
// The experiments are all suggestive of that photodissociation in ice and gas are of the same order of magnitude but that 
// the effective cross sections are much lower due to fast back-reactions in pure ices (Oberg et al. 2016);
#define ICE_PHOTODISS_EFFICIENCY 0.1 // estimate


//-----------------------------------------------
// Dust model parameters
//-----------------------------------------------

// Sticking coefficient for neutral species depends on gas and dust temperature; H atom and H2 molecule are special cases; 
// different workers use different values, usually 0.5 or 1.;
// Burke & Hollenbach, ApJ 265, p. 223 (1983); 
#define STICKING_COEFF_NEUTRALS 1.

// minimal radius of grains that can adsorb chemical species (in cm), 
// see the discussion by Hollenbach et al., ApJ 690, p. 1497 (2009):
#define MIN_ADSORPTION_RADIUS 2.e-7

// the grains with radius less then this values are strongly coupled to ion fluid;
#define MAX_ION_COUPLED_GRAIN_RADIUS 3.e-7  // in cm

// this number must be odd (to include zero charge);
#define NB_CHARGE_BINS_LARGE_GRAINS 31

//-------------------------------------------------
// Numerical parameters
//-------------------------------------------------

#define REL_ERROR_SOLVER 1.e-6
#define ABS_CONCENTRATION_ERROR_SOLVER 1.e-18 // specimen concentration (cm-3)
#define ABS_POPULATION_H2_ERROR_SOLVER 1.e-14 // specimen level population density (cm-3)
#define ABS_POPULATION_ERROR_SOLVER 1.e-11
#define ABS_PARAMETER_ERROR_SOLVER 1.e-11 // temperature (K), gas velocity (cm/s), average grain charge
#define MAX_CONV_FAILS_SOLVER 100 // default value is 10;
#define MAX_ERR_TEST_FAILS_SOLVER 14 // default value is 7;

#define MINIMAL_ABUNDANCE 1.e-99 // for saving in file
#define MINIMAL_REACTION_RATE 1.e-99
//--------------------------------------------------
// Switchers
//--------------------------------------------------

// 0 - results by Osterbrock, ApJ 134, p. 270 (1961), are used in calculating elastic scattering between ions and all neutrals
//    (the rate coefficient is constant);
// 1 - results by Flower, MNRAS 313, L19 (2000) are used (Osterbrock's cross section at low energies) for ion-H2 system;
#define FLOWER_ELASTIC_APPROX 1

// The ad-hoc formation of H2 on grain surface assumes that each accretion events of two H atoms leads to the formation of H2. 
// In this prescription when the ad-hoc H2 formation is activated, 50% of the adsorbed H are available for grain reactions 
// (other than H2 formation) and 50% for the formation of H2 (as in NAUTILUS code).
// 0 - H2 formation is modelled (grain surface chemistry must be on); 
// 1 - "ad hoc" (as described higher); 2 - semiempirical H2 formation rate; 
#define H2_FORMATION_MODE 0

// Hollenbach, Tielens, Reviews of Modern Physics 71, p. 173, 1999; k = c*n_H*n_H_tot, c = (1-3)*1e-17 cm3 s-1;
#define STANDARD_H2_FORMATION_RATE 3.e-17

// 1 - surface chemistry is on ("minimal" model - 0) 
#define GRAIN_SURFACE_CHEMISTRY_ON 1

// 1 - H and H2 migrate via quantum tunneling or thermal hopping depending on which is faster;
#define H_H2_QUANTUM_TUNNELING_DIFF_ON 1

// Grain heating by molecular emission,
// 0 - no heating taken into account ("minimal" model - 0)
#define GRAIN_HEATING_BY_MOLECULES_ON 0

// Dust heating by chemical reactions on the surface of dust grains
// may lead to the explosive release of the chemical energy stored as free radicals (Shen et al., A&A 415, p.203, 2004)
// 0 - is not allowed, 1 - is on
#define GRAIN_HEATING_GS_CHEMISTRY 0 

// 0 - no sputtering of grain mantle, this parameter is only for shock ("minimal" model - 0)
#define GRAIN_MANTLE_SPUTTERING_ON 1

// 1 - gas and dust temperatures are fixed in chemical evolution simulations, for shock this parameter does not work;
#define IS_TEMPERATURE_FIXED 0

// the data on radiative transfer parameters gamma and delta, dust heating efficiencies of molecule lines are saved;
// 1 - radiative transfer factors are saved;
#define SAVE_RADIATIVE_TRANSFER_FACTORS 0

// Calculate methanol levels, 1 - yes
#define CALCULATE_POPUL_METHANOL 0

// Calculate OH and NH3 populations, 1 - yes
#define CALCULATE_POPUL_NH3_OH 0

// 0 - switch off
// 1 - data are used by Tine et al., ApJ 481, p.282 (1997)
#define H2_CR_EXCITATION_ON 1

// 0 - switch off
// 1 - data on H2-H+ collisions are used by Gonzalez-Lezana & Honvault (2017),
// the question - are rates for collisions with H+ and H3+ equal?
// check this parameter in simulations of chemical evolution of dark cloud
#define H2_IONS_EXCITATION_ON 1

// 0 - Le Bourlot et al., MNRAS 332, 985, 2002
// 1 - Bossion et al. MNRAS 480, p.3718, 2018
#define H2_H_DISSOCIATION_MODE 1

// 0 - swith off
#define H2_FORMATION_EXCITATION_ON 1