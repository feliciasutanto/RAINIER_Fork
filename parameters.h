////////////////////////////////////////////////////////////////////////////////
////////////////////// Input Parameters ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#ifndef PARAMETERS_H
#define PARAMETERS_H

////////////////////// Discrete Settings ///////////////////////////////////////
//#define bPrintLvl // print both discrete and constructed lvl schemes
const int g_nZ = 60; // proton number
const int g_nAMass = 144; // proton + neutron number
const int g_nDisLvlMax = 16; // only trust level scheme to here, sets ECrit

///////////////////// Constructed Level Scheme Settings ////////////////////////
///// Bins /////
const int g_nConSpbMax = 21; // constructed # spin bins, small for light ion rxn

#define bForceBinNum // else force bin spacing
#ifdef bForceBinNum
double g_dConESpac; // constructed E bin spacing
const int g_nConEBin = 400; // number of energy bins in constructed scheme
#else
const double g_dConESpac = 0.01; // MeV; wont matter if forcing bin number
int g_nConEBin;
#endif // bForceBinNum

///// Level Density, LD, model /////
// choose one, fill in corresponding parameters
//#define bLD_CTM
#ifdef bLD_CTM
const double g_dTemp =  0.48473; // MeV
const double g_dE0   = -1.31817; // MeV
const double g_dDeuPair = 0.62834; // MeV; 
//const double g_dE0   = -1.004 + 0.5 * g_dDeuPair; // MeV; von egidy09 fit
#endif

#define bLD_BSFG
#ifdef bLD_BSFG
const double g_dE1 = 0.968; // MeV, excitation energy shift
const double g_dDeuPair = 2.698; // MeV; can get from ROBIN: Pa_prime
//const double g_dE1 = g_dDeuPair * 0.5 - 0.381; // von Egidy fit
#endif
// deuteron pairing energy from mass table, related to backshift in BSFG or CTM
// used for effective energy in LD, spincut, GSF models

// Level Density Parameter "a" = pi^2 / (6 * (n + p orbital spacing) )
//#define bLDaEx // a(Ex) = aAsym * (1 + dW * (1 - exp(-Gam * Eff) / dEff) )
#ifdef bLDaEx
const double g_dLDaAsym   = 14.58; // MeV^-1; Asymptotic value, a(Ex->Inf)
const double g_dDampGam   = 0.0; // Damping Parameter
const double g_dShellDelW = 0.0; // MeV; M_exp - M_LDM ~ shell correction
#else // const a
const double g_dLDa = 14.58; // MeV^-1 aka "LD parameter a" 
#endif

///// Spin Cutoff /////
// choose one:
//#define bJCut_VonEgidy05 // low-energy model
//#define bJCut_SingPart // single particle model
#define bJCut_RigidSph // rigid sphere model
//#define bJCut_VonEgidy09 // empirical fit 
//#define bJCut_Other // so many models of spin cut... code your own

///// Pairity Dependence /////
// choose one:
#define bPar_Equipar // equal +/- states, usually a good approx
//#define bPar_Edep // exponential asymptotic dependence
#ifdef bPar_Edep // 0.5 (1 +/- 1 / (1 + exp( g_dParC * (dEx - g_dParD) ) ) )
const double g_dParC = 3.0; // MeV^-1
const double g_dParD = 0.0; // MeV
#endif

///// Level spacing distribution /////
#define bPoisson // good approx to lvl spacing, might have more severe fluct
//#define bWigner // more representative of nuc lvl spacing, but more t to init

/////////////////////// Gamma Strength Function, GSF, Settings /////////////////
///// Width Fluctuations Distribution (WFD)-reason statistical codes exist /////
// choose one:
#define bWFD_PTD // fastest version of the Porter Thomas Distribution, nu = 1
//#define bWFD_nu // set the chi^2 degrees of freedom, nu - slow
//#define bWFD_Off // no fluctuations - nearly TALYS with level spac fluct
#ifdef bWFD_nu
const double g_dNu = 0.5; // See Koehler PRL105,072502(2010): measured nu~0.5
#endif

// Parameters for 56Fe from TALYS defaults usually an acceptable start
///// fE1 /////
// choose one:
#define bE1_GenLor // General Lorentzian
//#define bE1_EGLO // Enhanced Generalized Lorentzian for A>148
//#define b_GenLor_CT // Use a constant temperature in EGLO
//#define bE1_KMF // Kadmenskij Markushev Furman model
//#define bE1_KopChr // Kopecky Chrien model
//#define bE1_StdLor // standard Lorentzian
const double g_adSigE1[] = {317.00, 0.00}; // mb magnitude
const double g_adEneE1[] = { 15.05, 0.01}; // MeV centroid energy, non-zero
const double g_adGamE1[] = {  5.30, 0.00}; // MeV GDR width
//                                  ^^^^ for a 2nd resonance
#ifdef b_GenLor_CT // Use a constant temperature in EGLO
const double g_GenLor_CT = 0.44; // in MeV
#endif

///// fM1 /////
//#define bM1StrUpbend // Oslo observed low energy upbend aka enhancement
#ifdef bM1StrUpbend // soft pole behavior: C * exp(-A * Eg)
const double g_dUpbendM1Const = 5e-8; // C
const double g_dUpbendM1Exp = 1.0; // (positive) A
#endif

///// fE2 /////
//#define bE2_SingPart // single particle
#define bE2_StdLor // standard Lorentzian, parameterized by Prestwich
#ifdef bE2_SingPart
const double g_dSpSigE2 = 4e-11; // MeV^-5
#endif

////////////////////// Internal Conversion Coefficient, ICC, Settings //////////
//#define bUseICC // ICC = 0.0 otherwise
// if issues, turn ICC off and get your briccs to work outside RAINIER first
const char g_sBrIccModel[] = "BrIccFO"; // Conversion data table
const int g_nBinICC = 100; // Energy bins of BrIcc - more takes lot init time
const double g_dICCMin = g_dConESpac / 2.0; // uses 1st Ebin ICC val below this
const double g_dICCMax = 1.0; // MeV; Uses last Ebin ICC value for higher E

////////////////////// Run Settings ////////////////////////////////////////////
const int g_nReal = 1; // number of realizations of nuclear level scheme
const int g_nEvent = 1e4; // number of events per excitation and realization
const int g_nEvUpdate = 1e2; // print progress to screen at this interval

////////////////////// Excitation Settings /////////////////////////////////////
// choose one, fill in corresponding params:
//#define bExSingle  // single population input
//#define bExSelect // like Beta decay
//#define bExSpread  // ejectile detected input
#define bExFullRxn // no ejectile detected input

#ifdef bExSingle // similar to (n,g)
const double g_dExIMax = 7.8174; // MeV, Ei - "capture state energy"
const double g_dSpI    = 3.0; // hbar, Ji - "capture state spin"
const double g_dParI   = 0; // Pi - "capture state parity" 0=(-), 1=(+)
#endif

#ifdef bExSelect // similar to Beta decay
const double g_adExI[] = {0, 0.88489, 1.78662, 2.69324, 2.89503, 3.0383, 3.3418,
  3.47668, 3.59912, 3.78828, 3.8483, 3.9039, 4.0016, 4.06142, 4.2646, 4.30903, 
  4.46481, 4.5144, 4.5583, 4.5889, 4.7101, 4.7917, 4.8492, 5.0613}; // MeV
const double g_adSpI[]  = {0, 2, 4, 4, 6, 5, 3, 5, 6, 6, 5, 5, 6, 5, 6, 5, 7, 6,
  5, 7, 7, 5, 5, 5}; // hbar
const double g_anParI[] = {1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1,
  1, 0, 0, 1, 1, 0}; // Pi
const double g_adBRI[]   = {0, 0, 0, 0, 0.014, 0.58, 0.0094, 0.061, 0.0103, 
  0.0325, 0.049, 0.0057, 0.0074, 0.0141, 0.0216, 0.0236, 0.0415, 0.0135, 0.0298,
  0.014, 0.0104, 0.0057, 0.0519, 0.005}; // Branching Ratio: sums to 1.0
#endif

#ifdef bExSpread // similar to (p,p'), (a,a'), (he3,a'), etc.
const double g_dExIMax = 8.0; // MeV; constructed lvl scheme built up to this
// dont exceed with init excitations - gaus might sample higher than expected
const double g_adExIMean[] = {3.0, 4.0, 5.0, 6.0, 7.0}; // MeV
const double g_dSpIMean = 4.0; // poisson mean for initial spin state spread
const double g_dExISpread = 0.2 / 2.355; // MeV; std dev sigma = FWHM / 2.355 
const double g_dExRes = 0.2; // excitation resolution on g_ah2ExEg for analysis
#endif

#ifdef bExFullRxn // from a TALYS output file if available
const double g_dExIMax = 16.0; // MeV; above max population energy
const char popFile[] = "Nd144Pop.dat"; // made from TALYS "outpopulation y"
// make sure to match # of discrete bins
const double g_dExRes = 0.2 / 2.355; // excitation resolution on g_ah2ExEg
#endif

////////////////////// Analysis Settings ///////////////////////////////////////
const double g_dPlotSpMax = 10.0;
const int nExJBin = 100;  // n Ex Bins for plotting underlying J
///// JPop Analysis /////
//const int g_anPopLvl[] = {4,6,8,10,7,9};// low-ly populated lvls,0 = gs //56Fe
const int g_anPopLvl[] = {13,8,14,10,6,11}; // 144Nd
///// DRTSC Analysis
//const int g_anDRTSC[] = {1,3,5,8,12,13,15}; // 56Fe
const int g_anDRTSC[] = {1,4,6,15}; // 144Nd
const int g_nEgBin = 500;

////////////////////////////////////////////////////////////////////////////////
////////////////////// End Input Parameters ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/////////////////////////////// Program Includes ///////////////////////////////
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <iomanip> 
#include "math.h"
#include "TTimeStamp.h"
using namespace std;
#include "TRandom2.h" 
// 2=Tausworthe is faster and smaller than 3=Mersenne Twister (MT19937)
// search and replace "TRandom2" to "TRandom3" to change PRNG
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include <TROOT.h>
#include "TF1.h"
#include "TF2.h"
// determine OS for briccs
#ifdef __linux__
char cbriccs[] = "briccs";
#endif
#ifdef __MACH__
char cbriccs[] = "briccsMac";
#endif
#ifdef _WIN32
char cbriccs[] = "BrIccS.exe"; // haven't done any windows testing yet
#endif

/////////////////////////////// Parallel Settings //////////////////////////////
// should handle itself, email me if you get it to work on Mac or PC
#ifdef __CLING__
// cling in root6 won't parse omp.h. but man, root5 flies with 24 cores!
#else
#ifdef __linux__ // MacOS wont run omp.h by default, might exist workaround
    //#define bParallel // Parallel Option
    // ROOT hisograms not thread safe, but only miss ~1e-5 events
#endif // linux
#endif // cling

#ifdef bParallel
#ifndef CINT 
#include "omp.h" // for parallel on shared memory machine (not cluster yet)
#endif // cint
#endif // parallel

#ifdef bE2_StdLor
// Prestwich Physics A Atoms and Nuclei 315, 103-111 (1984)
const double g_dEneE2 = 63.0 * pow(g_nAMass,-1/3.0);
const double g_dGamE2 = 6.11 - 0.012 * g_nAMass;
const double g_dSigE2 = 1.5e-4 * pow(g_nZ,2.0) * g_dEneE2 / (pow(g_nAMass, 1/3.0) * g_dGamE2);
#endif

const double g_dKX1 = 8.673592583E-08; // mb^-1 MeV^-2;  = 1/(3*(pi*hbar*c)^2) 
const double g_adEneM1[] = {41.0 * pow(g_nAMass,-1/3.0), 0.01}; // MeV cent E
const double g_adGamM1[] = {4.00                       , 0.00}; // MeV GDR width
const double g_adSigM1[] = {
  // yes it is this complicated and phenomenological: fM1 = 1.58e-9 A^0.47
  1.58e-9 * pow(g_nAMass, 0.47) / (g_dKX1 * 7.0 * pow(g_adGamM1[0],2)) 
  * ( pow( pow(7.0,2) - pow(g_adEneM1[0],2),2) + pow(7.0 * g_adGamM1[0],2) )
                                                       , 0.00}; // mb magnitude
//                                                       ^^^^ for a 2nd res.
// I might put in the RIPL ratio formula as well: 
// just enter your own value for SigM1 if you have something else

#ifdef bExSelect
const int g_nStateI  = sizeof(g_adExI)     / sizeof(double);
const double g_dExIMax = g_adExI[g_nStateI-1] + 0.25; // build above last val
const double g_adExIMean[] = {0.0}; // unused in select
const double g_dExISpread = 0.0; // unused in select
const double g_dExRes = 0.0; // unused in select
#endif
#ifdef bExSingle
const double g_adExIMean[] = {0.0}; // unused in single
const double g_dExISpread = 0.0; // unused in single
const double g_dExRes = 0.0; // unused in single
#endif
#ifdef bExFullRxn
const double g_adExIMean[] = {0.0}; // unsed in full rxn
const double g_dExISpread = 0.0; // unsed in full rxn
#endif
const int g_nExIMean = sizeof(g_adExIMean) / sizeof(double); // self adjusting
const int g_nParE1   = sizeof(g_adSigE1)   / sizeof(double);
const int g_nParM1   = sizeof(g_adSigM1)   / sizeof(double);
const int g_nPopLvl  = sizeof(g_anPopLvl)  / sizeof(int);
const int g_nDRTSC   = sizeof(g_anDRTSC)   / sizeof(int);

const bool g_bIsEvenA = !(g_nAMass % 2);
const int g_nDisLvlGamMax = 15; // max gammas in for a discrete lvl

#endif //PARAMETERS_H