#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <omp.h>

#include "reprimand/con2prim_imhd.h"
#include "reprimand/eos_hybrid.h"
#include "reprimand/eos_barotr_pwpoly.h"

using namespace EOS_Toolkit;


int main(int argc, char **argv)
{
  if(argc != 2) {
    std::cerr << "usage: " << argv[0] << " datafile.txt\n";
    return 1;
  }

  std::string datafilename(argv[1]);

  //Get some EOS
  const real_t max_eps = 11.;
  const real_t max_rho = 1e6; //fixmy
  
  // M1 from Read et al. (based on Mueller and Serot) from Kastaun et al. paper
  // on RePrimand, using just a single crust piece.
  // Refer to Fig2 and the text referring to it in Read et al. for how to
  // compute these, but *note* incorrect value for SLy crust EOS given there,
  // KCrust is 3.594473308408096 in cgs units for the Gamma = 1.35692 crust.
  // See MS1.ipynb.
  //
  // Specifically I use *old* values of G and Msun:
  // G =  6.67300*1.e-8  % cm^3 g^-1 s^-2
  // Msun = 1.98892*1e33 % g
  // (semi-)independent table of coeffs also available from:
  // http://www.computational-relativity.org/assets/eos/EOS_MS1.txt
  // RePrimIand defines a "polytropic density scale rho_p" on
  // https://wokast.github.io/RePrimAnd/eos_barotr_available.html
  //
  // P = rho_P (rho/rho_P)**Gamma
  //
  // and comparing to the common
  //
  // P = K rho**Gamma
  //
  // we have
  //
  // K = rho_p/rho_P**Gamma = 1/(rho_P**(Gamma-1))
  //
  // or alternatively
  //
  // rho_P = 1/K**(1/(Gamma-1))
  const int nsegs = 4;
  const double K0 = 0.08950758861673326;
  std::vector<real_t> rho_bounds(nsegs);
  rho_bounds[0] = 0.; // always zero for first segment
  rho_bounds[1] = 0.00015247493312376816;  // determined by intersection point with crust EOS
  rho_bounds[2] = 0.000811456143270882; // fixed for all EOS
  rho_bounds[3] = 0.00161906786291838;  // fixed for all EOS
  std::vector<real_t> gammas(nsegs);
  gammas[0] = 1.35692;
  gammas[1] = 3.224;
  gammas[2] = 3.033;
  gammas[3] = 1.325;
  const real_t rmdp0 = 1/pow(K0, gammas[0]-1);
  auto eos_c = make_eos_barotr_pwpoly(rmdp0, rho_bounds, gammas, max_rho);

  // thermal bit
  const real_t gamma_th = 1.8; // used by Kastaun et al. no other reason
  auto eos = make_eos_hybrid(eos_c, gamma_th, max_eps, max_rho);
  
  //Set up atmosphere
  // minimal density
  real_t atmo_rho = 1e-20;
  real_t atmo_eps = eos_c.at_rho(atmo_rho).eps();
  real_t atmo_ye = 0.5;
  real_t atmo_cut = atmo_rho * 1.01;
  real_t atmo_p = eos_c.at_rho(atmo_rho).press();

  atmosphere atmo{atmo_rho, atmo_eps, atmo_ye, atmo_p, atmo_cut};

  //Primitive recovery parameters 
  // density below which RePrimAnd's con2prim becomes more lenient
  real_t rho_strict = 1e-11;
  bool  ye_lenient = false;
  int max_iter = 30;
  real_t c2p_acc = 1e-8;
  real_t max_b = 10.;
  real_t max_z = 1e3;
  
  //Get a recovery function
  con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, max_z, max_b, 
                     atmo, c2p_acc, max_iter);


  // conserved vairables to recover
  std::vector<real_t> dens, tau, scon;
  std::ifstream datafile(datafilename);
#if 1
  // one set of data per line:
  // D tau S
  double densval, tauval, sconval;
  while(datafile >> densval >> tauval >> sconval) {
    dens.push_back(densval);
    tau.push_back(tauval);
    scon.push_back(sconval);
  }
#endif
#if 0 // compute our own conserved values
  double rhoval, vval;
  while(datafile >> rhoval >> vval) {
    real_t epsval = eos_c.at_rho(rhoval).eps();
    real_t pval = eos_c.at_rho(rhoval).press();

    // hard-code metric to be Minkowsky
    real_t w_lorentzval = 1./sqrt(1.-vval*vval);

    real_t densval = rhoval * w_lorentzval;
    real_t hval = 1. + epsval + pval/rhoval;
    real_t enval = rhoval*hval*w_lorentzval*w_lorentzval - pval;
    real_t tauval = enval - densval;
    real_t sconval = rhoval*hval*w_lorentzval*w_lorentzval*vval;

    dens.push_back(densval);
    tau.push_back(tauval);
    scon.push_back(sconval);
  }
#endif
  datafile.close();

  // output
  std::vector<real_t> rho(dens.size());
  std::vector<real_t> eps(dens.size());
  std::vector<real_t> press(dens.size());
  std::vector<real_t> vel(dens.size());

  sm_metric3 g;
  g.minkowski();
  double start_time = omp_get_wtime();
  for(size_t i = 0 ; i < dens.size() ; ++i) {
    //collect
    const real_t Y_e = 0.5;
    cons_vars_mhd evolved{dens[i], tau[i], dens[i]*Y_e, 
                          {scon[i],0.,0.}, {0.,0.,0.}};    
    prim_vars_mhd primitives;
    con2prim_mhd::report rep;
    
    //recover
    cv2pv(primitives, evolved, g, rep);
    
    //check
    if (rep.failed()) {
      std::cerr << rep.debug_message(); 
      //abort
    } else {
      // debug output
#if 0
      std::cout << std::setprecision(18);
      std::cout << "rho: " << primitives.rho << "\n";
      std::cout << "eps: " << primitives.eps << "\n";
      std::cout << "press: " << primitives.press << "\n";
      std::cout << "vel: " << primitives.vel(0) << " "
                << primitives.vel(1) << " " << primitives.vel(2) << "\n";
#endif

      //write back primitive vars
      rho[i] = primitives.rho;
      eps[i] = primitives.eps;
      press[i] = primitives.press;
      vel[i] = primitives.vel(0);
    }
  }
  double end_time = omp_get_wtime();
  std::cout << "took: " << (end_time - start_time) << " seconds\n\n";

  std::cerr << "ignore any 'invalid next size' errors due to Python(?)\n\n";

  // write results to disk
  std::ofstream results("output.txt");
  results << std::setprecision(18);
  for(size_t i = 0 ; i < dens.size() ; ++i) {
    results << dens[i] << " " << tau[i] << " " << scon[i]
            << " " << rho[i]
            << " " << eps[i]
            << " " << press[i]
            << " " << vel[i]
            << "\n";
  }
  results.close();
  
  return 0;
}
