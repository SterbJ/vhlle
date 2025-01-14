/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  Comput. Phys. Commun. 185 (2014), 3016   arXiv:1312.4160                   *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it.                   *
*
 
 NEEDS SOME LINEARIZATION IN NB, NS, NQ !!!
 
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <unistd.h>

#include <complex>
#include <vector>
#include <random>

#include "hdo_box.h"
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "eos.h"
#include "cll.h"
#include "trancoeff.h"

#include "icBox.h"

extern "C" {
#include "../kiss_fft130/kiss_fft.h"
#include "../kiss_fft130/tools/kiss_fftnd.h"
}

using namespace std;

// extern bool debugRiemann ;

double sign(double x) {
 if (x > 0)
  return 1.;
 else if (x < 0.)
  return -1.;
 else
  return 0.;
}

// this version contains NO PRE-ADVECTION for the IS solution

// enable this to use formal solution for the relaxation part of
// Israel-Stewart equations (not recommended)
//#define FORMAL_SOLUTION
// else: use 1st order finite difference update

Hydro::Hydro(Fluid *_f, EoS *_eos, TransportCoeff *_trcoeff, double _t0,
             double _dt) {
 eos = _eos;
 trcoeff = _trcoeff;
 f = _f;
 dt = _dt;
 #ifdef CARTESIAN
 t = _t0;
 tau = 1.0;
 #else
 tau = _t0;
 #endif
}

Hydro::~Hydro() {}

void Hydro::setDtau(double deltaTau) {
 dt = deltaTau;
 if (dt > f->getDx() / 2. ||
     dt > f->getDy() / 2. /*|| dt>tau*f->getDz()/2. */) {
  cout << "too big delta_tau " << dt << "  " << f->getDx() << "  " << f->getDy()
       << "  " << tau * f->getDz() << endl;
  exit(1);
 }
}

void Hydro::hlle_flux(Cell *left, Cell *right, int direction, int mode, int ix, int iy, int iz) {
 // for all variables, suffix "l" = left state, "r" = right state
 // with respect to the cell boundary
 // idea is that e.g. el corresponds to the fluctuation
 // the background variable is usually (except for some special cases) denoted by 0 or _0
 double d_el, d_er, d_pl, d_pr, d_nbl, d_nql, d_nsl, d_nbr, d_nqr, d_nsr, d_vxl, d_vxr, d_vyl, d_vyr, d_vzl,
     d_vzr, bl = 0., br = 0., csb, vb, vb_const, vb_delta, El, Er, dx = 0.;
 double el_bck, er_bck, pl_bck, pr_bck, nbl_bck, nql_bck, nsl_bck, nbr_bck, nqr_bck, nsr_bck, vxl_bck, vxr_bck, vyl_bck, vyr_bck, vzl_bck, vzr_bck;
 double d_Ftl = 0., d_Fxl = 0., d_Fyl = 0., d_Fzl = 0., d_Fbl = 0., d_Fql = 0., d_Fsl = 0.,
        d_Ftr = 0., d_Fxr = 0., d_Fyr = 0., d_Fzr = 0., d_Fbr = 0., d_Fqr = 0., d_Fsr = 0.;
 double d_U1l, d_U2l, d_U3l, d_U4l, d_Ubl, d_Uql, d_Usl, d_U1r, d_U2r, d_U3r, d_U4r, d_Ubr, d_Uqr, d_Usr;
 double d_flux[7];
 const double dta = mode == 0 ? dt / 2. : dt;
 double tauFactor = 1.0;  // fluxes are also multiplied by tau, 1 for Cartesian
    

 if (mode == PREDICT) {
  // get primitive quantities from Q_{i+} at previous timestep
  left->getPrimVarRightQbck(eos, tau, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck,
                           direction);
  left->getPrimVarRight(eos, tau, d_el, d_pl, d_nbl, d_nql, d_nsl, d_vxl, d_vyl, d_vzl,
                        direction, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck);
  // ... and Q_{(i+1)-}
  right->getPrimVarLeftQbck(eos, tau, er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck,
                           direction);
  right->getPrimVarLeft(eos, tau, d_er, d_pr, d_nbr, d_nqr, d_nsr, d_vxr, d_vyr, d_vzr,
                        direction, er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck);

  El = (el_bck + pl_bck) / (1 - vxl_bck * vxl_bck - vyl_bck * vyl_bck - vzl_bck * vzl_bck);
  Er = (er_bck + pr_bck) / (1 - vxr_bck * vxr_bck - vyr_bck * vyr_bck - vzr_bck * vzr_bck);
  #ifndef CARTESIAN
  tauFactor = tau + 0.25 * dt;
  #endif

 } else {
  // use half-step updated Q's for corrector step
  left->getPrimVarHRightQbck(eos, tau, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck,
                            direction);
  left->getPrimVarHRight(eos, tau, d_el, d_pl, d_nbl, d_nql, d_nsl, d_vxl, d_vyl, d_vzl,
                         direction, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck);
  right->getPrimVarHLeftQbck(eos, tau, er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck,
                            direction);
  right->getPrimVarHLeft(eos, tau, d_er, d_pr, d_nbr, d_nqr, d_nsr, d_vxr, d_vyr, d_vzr,
                         direction, er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck);

  El = (el_bck + pl_bck) / (1 - vxl_bck * vxl_bck - vyl_bck * vyl_bck - vzl_bck * vzl_bck);
  Er = (er_bck + pr_bck) / (1 - vxr_bck * vxr_bck - vyr_bck * vyr_bck - vzr_bck * vzr_bck);
  #ifndef CARTESIAN
  tauFactor = tau + 0.5 * dt;
  #endif
 }

// if (d_el < -el_bck || el_bck < 0.) {//modified condition - do we allow fluctuations in empty cell?
//     d_el = -el_bck;
//     d_pl = -pl_bck;
// }
// if (d_er < -er_bck || er_bck < 0.) {//modified condition
//     d_er = -er_bck;
//     d_pr = -pr_bck;
// }
//    if (d_el > el_bck) {//modified condition - do we allow fluctuations in empty cell?
//        d_el = el_bck;
//        d_pl = pl_bck;
//    }
//    if (d_er > er_bck) {//modified condition
//        d_er = er_bck;
//        d_pr = pr_bck;
//    }

 if (d_el+el_bck > 1e10 || el_bck > 1e10) {//modified condition
  cout << "e>1e10; debug info below:\n";
  left->Dump(tau);
  // debugRiemann = true ;
     if (mode == PREDICT){
         left->getPrimVarRightQbck(eos, tau, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck,
                               direction);
         left->getPrimVarRight(eos, tau, d_el, d_pl, d_nbl, d_nql, d_nsl, d_vxl, d_vyl, d_vzl,
                               direction, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck);
     }
     else{
         left->getPrimVarHRightQbck(eos, tau, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck,
                                direction);
         left->getPrimVarHRight(eos, tau, d_el, d_pl, d_nbl, d_nql, d_nsl, d_vxl, d_vyl, d_vzl,
                                direction, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck);
     }
  // debugRiemann = false ;
  exit(0);
 }

 // skip the procedure for two empty cells
 //if ((d_el+el_bck < 0. && d_er+er_bck < 0.) || (el_bck == 0. && er_bck == 0.)) return;//modified condition - thershold!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (d_pr + pr_bck < 0. || pr_bck < 0.) {//modified condition
//     return;
  cout << "Negative pressure    " << d_pr << "      " << d_er << "  " << ix << "    " << iy << "    " << iz << endl;
//  left->getPrimVarRightQbck(eos, tau, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck,
//                           direction);
//  left->getPrimVarRight(eos, tau, d_el, d_pl, d_nbl, d_nql, d_nsl, d_vxl, d_vyl, d_vzl,
//                        direction, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck);
//  right->getPrimVarLeftQbck(eos, tau, er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck,
//                           direction);
//  right->getPrimVarLeft(eos, tau, d_er, d_pr, d_nbr, d_nqr, d_nsr, d_vxr, d_vyr, d_vzr,
//                        direction, er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck);
 }

 // skip the procedure for two partially vacuum cells
 if (left->getM(direction) < 1. && right->getM(direction) < 1.) return;
    
    double gammal_bck = 1./ sqrt(1. - vxl_bck * vxl_bck - vyl_bck * vyl_bck - vzl_bck * vzl_bck);
    double gammar_bck = 1./ sqrt(1. - vxr_bck * vxr_bck - vyr_bck * vyr_bck - vzr_bck * vzr_bck);
    
    d_U1l = (el_bck + pl_bck) * (gammal_bck * gammal_bck * d_vxl ) + (d_el + d_pl) * gammal_bck * gammal_bck * vxl_bck;
    d_U2l = (el_bck + pl_bck) * (gammal_bck * gammal_bck * d_vyl ) + (d_el + d_pl) * gammal_bck * gammal_bck * vyl_bck;
    d_U3l = (el_bck + pl_bck) * (gammal_bck * gammal_bck * d_vzl ) + (d_el + d_pl) * gammal_bck * gammal_bck * vzl_bck;
    d_U4l = (d_el + d_pl) * gammal_bck * gammal_bck - d_pl;

    d_U1r = (er_bck + pr_bck) * (gammar_bck * gammar_bck * d_vxr ) + (d_er + d_pr) * gammar_bck * gammar_bck * vxr_bck;
    d_U2r = (er_bck + pr_bck) * (gammar_bck * gammar_bck * d_vyr ) + (d_er + d_pr) * gammar_bck * gammar_bck * vyr_bck;
    d_U3r = (er_bck + pr_bck) * (gammar_bck * gammar_bck * d_vzr ) + (d_er + d_pr) * gammar_bck * gammar_bck * vzr_bck;
    d_U4r = (d_er + d_pr) * gammar_bck * gammar_bck - d_pr;

    d_Ubl = gammal_bck * d_nbl; // not linearized
    d_Uql = gammal_bck * d_nql; // not linearized
    d_Usl = gammal_bck * d_nsl; // not linearized

    d_Ubr = gammar_bck * d_nbr; // not linearized
    d_Uqr = gammar_bck * d_nqr; // not linearized
    d_Usr = gammar_bck * d_nsr; // not linearized

 if (direction == X_) {
     
     d_Ftl = (el_bck + pl_bck) * (gammal_bck * gammal_bck * d_vxl ) + (d_el + d_pl) * gammal_bck * gammal_bck * vxl_bck;
     d_Fxl = (el_bck + pl_bck) * (gammal_bck * vxl_bck * gammal_bck * d_vxl + gammal_bck * vxl_bck * gammal_bck * d_vxl) + (d_el + d_pl) * gammal_bck * vxl_bck * gammal_bck * vxl_bck + d_pl;
     d_Fyl = (el_bck + pl_bck) * (gammal_bck * vyl_bck * gammal_bck * d_vxl + gammal_bck * vxl_bck * gammal_bck * d_vyl) + (d_el + d_pl) * gammal_bck * vyl_bck * gammal_bck * vxl_bck;
     d_Fzl = (el_bck + pl_bck) * (gammal_bck * vzl_bck * gammal_bck * d_vxl + gammal_bck * vxl_bck * gammal_bck * d_vzl) + (d_el + d_pl) * gammal_bck * vzl_bck * gammal_bck * vxl_bck;
     
     d_Ftr = (er_bck + pr_bck) * (gammar_bck * gammar_bck * d_vxr ) + (d_er + d_pr) * gammar_bck * gammar_bck * vxr_bck;
     d_Fxr = (er_bck + pr_bck) * (gammar_bck * vxr_bck * gammar_bck * d_vxr + gammar_bck * vxr_bck * gammar_bck * d_vxr) + (d_er + d_pr) * gammar_bck * vxr_bck * gammar_bck * vxr_bck + d_pr;
     d_Fyr = (er_bck + pr_bck) * (gammar_bck * vyr_bck * gammar_bck * d_vxr + gammar_bck * vxr_bck * gammar_bck * d_vyr) + (d_er + d_pr) * gammar_bck * vyr_bck * gammar_bck * vxr_bck;
     d_Fzr = (er_bck + pr_bck) * (gammar_bck * vzr_bck * gammar_bck * d_vxr + gammar_bck * vxr_bck * gammar_bck * d_vzr) + (d_er + d_pr) * gammar_bck * vzr_bck * gammar_bck * vxr_bck;
     
     d_Fbl = d_Ubl * d_vxl; // not linearized
     d_Fql = d_Uql * d_vxl; // not linearized
     d_Fsl = d_Usl * d_vxl; // not linearized

     d_Fbr = d_Ubr * d_vxr; // not linearized
     d_Fqr = d_Uqr * d_vxr; // not linearized
     d_Fsr = d_Usr * d_vxr; // not linearized

     // for the case of constant c_s only - takes into account only background
     csb = sqrt(eos->cs2() +
                0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vxl_bck - vxr_bck, 2));
     vb = (sqrt(El) * vxl_bck + sqrt(Er) * vxr_bck) / (sqrt(El) + sqrt(Er));
     bl = min(0., min((vb - csb) / (1 - vb * csb),
                      (vxl_bck - eos->cs()) / (1 - vxl_bck * eos->cs())));
     br = max(0., max((vb + csb) / (1 + vb * csb),
                      (vxr_bck + eos->cs()) / (1 + vxr_bck * eos->cs())));

     dx = f->getDx();

     // bl or br in the case of boundary with vacuum
     if (d_el + el_bck < 10e-18 || el_bck == 0.) bl = -1.;//modified condition
     if (d_er + er_bck < 10e-18 || er_bck == 0.) br = 1.;//modified condition
 }
 if (direction == Y_) {
     
     d_Ftl = (el_bck + pl_bck) * (gammal_bck * gammal_bck * d_vyl ) + (d_el + d_pl) * gammal_bck * gammal_bck * vyl_bck;
     d_Fxl = (el_bck + pl_bck) * (gammal_bck * vxl_bck * gammal_bck * d_vyl + gammal_bck * vyl_bck * gammal_bck * d_vxl) + (d_el + d_pl) * gammal_bck * vxl_bck * gammal_bck * vyl_bck;
     d_Fyl = (el_bck + pl_bck) * (gammal_bck * vyl_bck * gammal_bck * d_vyl + gammal_bck * vyl_bck * gammal_bck * d_vyl) + (d_el + d_pl) * gammal_bck * vyl_bck * gammal_bck * vyl_bck + d_pl;
     d_Fzl = (el_bck + pl_bck) * (gammal_bck * vzl_bck * gammal_bck * d_vyl + gammal_bck * vyl_bck * gammal_bck * d_vzl) + (d_el + d_pl) * gammal_bck * vzl_bck * gammal_bck * vyl_bck;
     
     d_Ftr = (er_bck + pr_bck) * (gammar_bck * gammar_bck * d_vyr ) + (d_er + d_pr) * gammar_bck * gammar_bck * vyr_bck;
     d_Fxr = (er_bck + pr_bck) * (gammar_bck * vxr_bck * gammar_bck * d_vyr + gammar_bck * vyr_bck * gammar_bck * d_vxr) + (d_el + d_pr) * gammar_bck * vxr_bck * gammar_bck * vyr_bck;
     d_Fyr = (er_bck + pr_bck) * (gammar_bck * vyr_bck * gammar_bck * d_vyr + gammar_bck * vyr_bck * gammar_bck * d_vyr) + (d_el + d_pr) * gammar_bck * vyr_bck * gammar_bck * vyr_bck + d_pr;
     d_Fzr = (er_bck + pr_bck) * (gammar_bck * vzr_bck * gammar_bck * d_vyr + gammar_bck * vyr_bck * gammar_bck * d_vzr) + (d_el + d_pr) * gammar_bck * vzr_bck * gammar_bck * vyr_bck;
     
     d_Fbl = d_Ubl * d_vyl; // not linearized
     d_Fql = d_Uql * d_vyl; // not linearized
     d_Fsl = d_Usl * d_vyl; // not linearized

     d_Fbr = d_Ubr * d_vyr; // not linearized
     d_Fqr = d_Uqr * d_vyr; // not linearized
     d_Fsr = d_Usr * d_vyr; // not linearized

     // for the case of constant c_s only - takes into account only background
     csb = sqrt(eos->cs2() +
                0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vyl_bck - vyr_bck, 2));
     vb = (sqrt(El) * vyl_bck + sqrt(Er) * vyr_bck) / (sqrt(El) + sqrt(Er));
     bl = min(0., min((vb - csb) / (1 - vb * csb),
                      (vyl_bck - eos->cs()) / (1 - vyl_bck * eos->cs())));
     br = max(0., max((vb + csb) / (1 + vb * csb),
                      (vyr_bck + eos->cs()) / (1 + vyr_bck * eos->cs())));

     dx = f->getDy();

     // bl or br in the case of boundary with vacuum
     if (d_el + el_bck < 10e-18 || el_bck == 0.) bl = -1.;//modified condition
     if (d_er + er_bck < 10e-18 || er_bck == 0.) br = 1.;//modified condition
 }
 if (direction == Z_) {
     double tau1 = tauFactor;
     
     d_Ftl = (el_bck + pl_bck) * (gammal_bck * gammal_bck * d_vzl ) / tau1 + (d_el + d_pl) * gammal_bck * gammal_bck * vzl_bck / tau1;
     d_Fxl = (el_bck + pl_bck) * (gammal_bck * vxl_bck * gammal_bck * d_vzl + gammal_bck * vzl_bck * gammal_bck * d_vxl) / tau1 + (d_el + d_pl) * gammal_bck * vxl_bck * gammal_bck * vzl_bck / tau1;
     d_Fyl = (el_bck + pl_bck) * (gammal_bck * vyl_bck * gammal_bck * d_vzl + gammal_bck * vzl_bck * gammal_bck * d_vyl) / tau1 + (d_el + d_pl) * gammal_bck * vyl_bck * gammal_bck * vzl_bck / tau1;
     d_Fzl = (el_bck + pl_bck) * (gammal_bck * vzl_bck * gammal_bck * d_vzl + gammal_bck * vzl_bck * gammal_bck * d_vzl) / tau1 + (d_el + d_pl) * gammal_bck * vzl_bck * gammal_bck * vzl_bck / tau1 + d_pl / tau1;
     
     d_Ftr = (er_bck + pr_bck) * (gammar_bck * gammar_bck * d_vzr ) / tau1 + (d_er + d_pr) * gammar_bck * gammar_bck * vzr_bck / tau1;
     d_Fxr = (er_bck + pr_bck) * (gammar_bck * vxr_bck * gammar_bck * d_vzr + gammar_bck * vzr_bck * gammar_bck * d_vxr) / tau1 + (d_er + d_pr) * gammar_bck * vxr_bck * gammar_bck * vzr_bck / tau1;
     d_Fyr = (er_bck + pr_bck) * (gammar_bck * vyr_bck * gammar_bck * d_vzr + gammar_bck * vzr_bck * gammar_bck * d_vyr) / tau1 + (d_er + d_pr) * gammar_bck * vyr_bck * gammar_bck * vzr_bck / tau1;
     d_Fzr = (er_bck + pr_bck) * (gammar_bck * vzr_bck * gammar_bck * d_vzr + gammar_bck * vzr_bck * gammar_bck * d_vzr) / tau1 + (d_er + d_pr) * gammar_bck * vzr_bck * gammar_bck * vzr_bck / tau1 + d_pr / tau1;
     
     d_Fbl = d_Ubl * d_vzl / tau1; // not linearized
     d_Fql = d_Uql * d_vzl / tau1; // not linearized
     d_Fsl = d_Usl * d_vzl / tau1; // not linearized

     d_Fbr = d_Ubr * d_vzr / tau1; // not linearized
     d_Fqr = d_Uqr * d_vzr / tau1; // not linearized
     d_Fsr = d_Usr * d_vzr / tau1; // not linearized

     // for the case of constant c_s only
     // factor 1/tau accounts for eta-coordinate

     // different estimate - takes into account only background
     csb = sqrt(eos->cs2() +
                0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vzl_bck - vzr_bck, 2));
     vb = (sqrt(El) * vzl_bck + sqrt(Er) * vzr_bck) / (sqrt(El) + sqrt(Er));
     bl = 1. / tau * min(0., min((vb - csb) / (1 - vb * csb),
                                 (vzl_bck - eos->cs()) / (1 - vzl_bck * eos->cs())));
     br = 1. / tau * max(0., max((vb + csb) / (1 + vb * csb),
                                 (vzr_bck + eos->cs()) / (1 + vzr_bck * eos->cs())));

     dx = f->getDz();

     // bl or br in the case of boundary with vacuum
     if (d_el + el_bck < 10e-18 || el_bck == 0.) bl = -1. / tau;//modified condition
     if (d_er + er_bck < 10e-18 || er_bck == 0.) br = 1. / tau;//modified condition
 }

 if (bl == 0. && br == 0.) return;

//  finally, HLLE formula for the fluxes
 d_flux[T_] = tauFactor * dta / dx *
            (-bl * br * (d_U4l - d_U4r) + br * d_Ftl - bl * d_Ftr) / (-bl + br);
 d_flux[X_] = tauFactor * dta / dx *
            (-bl * br * (d_U1l - d_U1r) + br * d_Fxl - bl * d_Fxr) / (-bl + br);
 d_flux[Y_] = tauFactor * dta / dx *
            (-bl * br * (d_U2l - d_U2r) + br * d_Fyl - bl * d_Fyr) / (-bl + br);
 d_flux[Z_] = tauFactor * dta / dx *
            (-bl * br * (d_U3l - d_U3r) + br * d_Fzl - bl * d_Fzr) / (-bl + br);
 d_flux[NB_] = tauFactor * dta / dx *
             (-bl * br * (d_Ubl - d_Ubr) + br * d_Fbl - bl * d_Fbr) / (-bl + br);
 d_flux[NQ_] = tauFactor * dta / dx *
             (-bl * br * (d_Uql - d_Uqr) + br * d_Fql - bl * d_Fqr) / (-bl + br);
 d_flux[NS_] = tauFactor * dta / dx *
             (-bl * br * (d_Usl - d_Usr) + br * d_Fsl - bl * d_Fsr) / (-bl + br);

//    cout << d_flux[T_] << endl;
 if (d_flux[NB_] != d_flux[NB_]) {  // if things failed
  cout << "---- error in hlle_flux: f_nb undefined!\n";
  cout << setw(12) << d_U4l << setw(12) << d_U1l << setw(12) << d_U2l << setw(12)
       << d_U3l << endl;
  cout << setw(12) << d_U4r << setw(12) << d_U1r << setw(12) << d_U2r << setw(12)
       << d_U3r << endl;
  cout << setw(12) << d_Ubl << setw(12) << d_Uql << setw(12) << d_Usl << endl;
  cout << setw(12) << d_Ubr << setw(12) << d_Uqr << setw(12) << d_Usr << endl;
  cout << setw(12) << d_Ftl << setw(12) << d_Fxl << setw(12) << d_Fyl << setw(12)
       << d_Fzl << endl;
  cout << setw(12) << d_Ftr << setw(12) << d_Fxr << setw(12) << d_Fyr << setw(12)
       << d_Fzr << endl;
  exit(1);
 }

 // update the cumulative fluxes in both neighbouring cells
    double Tl[7];
    double Tr[7];
    left->getQ(Tl);
    right->getQ(Tr);
    if ((abs(Tl[0]/6 - d_flux[0]) > el_bck/6) || (abs(Tr[0]/6 + d_flux[0]) > er_bck/6)) {
        double max = max( abs( (Tl[0]/6 - d_flux[0])/el_bck/6 ), abs( (Tr[0]/6 + d_flux[0])/er_bck/6 ) );
        d_flux[0] = d_flux[0] / max * 0.9;
        N_id++;
        //                cout << "v" << endl;
    }
    for (int i=1; i<4; i++) {
        if ((abs(Tl[i] - d_flux[i]) > abs(Tl[0] - d_flux[0])/2./sqrt(3.)) || (abs(Tr[i] + d_flux[i]) > abs(Tr[0] + d_flux[0])/2./sqrt(3.))) {
            double maximum = std::max( abs( (Tl[i] - d_flux[i])/(Tl[0] - d_flux[0]) ), abs( (Tr[i] + d_flux[i])/(Tr[0] + d_flux[0]) ) );
            d_flux[i] = d_flux[i] / maximum * 0.9;
        }
    }
         left->addFlux(-d_flux[T_], -d_flux[X_], -d_flux[Y_], -d_flux[Z_], 0., 0., 0.);
         right->addFlux(d_flux[T_], d_flux[X_], d_flux[Y_], d_flux[Z_], 0., 0., 0.);
}


void Hydro::source(double tau1, double x, double y, double z, double Q[7],
                   double S[7]) {
 #ifdef CARTESIAN
 // geometrical source term is zero in Cartesian frame
 for (int i = 0; i < 7; i++) S[i] = 0.0;
//    double e_0=1., nb=0., nq=0., ns=0.;
//    double T0, mub,muq,mus, p_0;
//    eos->eos(e_0, nb, nq, ns, T0, mub, muq, mus, p_0);
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    double mean = 0.;
//    double sigma = 10;
//    std::normal_distribution<double> d(mean, sigma);
//    double Xi = d(gen);
//    double dXi = - gaussian_derivative(Xi, mean, sigma); // gaussian_derivative gives absolute value - minus sign has to be added
//    S[X_] = dXi;
 #else
 // not linearized
 double _Q[7], _Q_bck[7], e, p, nb, nq, ns, vx, vy, vz;
    for (int i = 0; i < 7; i++){
        _Q[i] = Q[i] / tau1;  // no tau factor in  _Q
        _Q_bck[i] = Q_bck[i] / tau1;
    }
 transformPV(eos, _Q, _Q_bck, e, p, nb, nq, ns, vx, vy, vz);
 S[T_] = -_Q[T_] * vz * vz - p * (1. + vz * vz);
 S[X_] = 0.;
 S[Y_] = 0.;
 S[Z_] = -_Q[Z_];
 S[NB_] = 0.;
 S[NQ_] = 0.;
 S[NS_] = 0.;
 #endif
}

void Hydro::source_step(int ix, int iy, int iz, int mode) {
 double _dt;
 if (mode == PREDICT)
  _dt = dt / 2.;
 else
  _dt = dt;

 double tau1;
 double Q[7];
 double k[7];

 double x = f->getX(ix), y = f->getY(iy), z = f->getZ(iz);
 Cell *c = f->getCell(ix, iy, iz);

 if (mode == PREDICT) {
  c->getQ(Q);
  #ifdef CARTESIAN
  tau1 = t;
  #else
  tau1 = tau;
  #endif
 } else {
  c->getQh(Q);
  #ifdef CARTESIAN
  tau1 = t + 0.5 * dt;
  #else
  tau1 = tau + 0.5 * dt;
  #endif
 }
 source(tau1, x, y, z, Q, k);
 for (int i = 0; i < 7; i++) k[i] *= _dt;

 if (k[NB_] != k[NB_]) {  // something failed
  cout << "---- error in source_step: k_nb undefined!\n";
  cout << setw(12) << k[0] << setw(12) << k[1] << setw(12) << k[2] << setw(12)
       << k[3] << endl;
  cout << setw(12) << k[4] << setw(12) << k[5] << setw(12) << k[6] << endl;
  exit(1);
 }
 c->addFlux(k[T_], k[X_], k[Y_], k[Z_], k[NB_], k[NQ_], k[NS_]);
}

void Hydro::visc_source_step(int ix, int iy, int iz) {
 double e, p, nb, nq, ns, vx, vy, vz;
    double e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck;
 double uuu[4];
    double uuubck[4];
 double k[7];

 #ifdef CARTESIAN
  // there is no geometric viscous source term in Cartesian frame
 #else
  // not linearized
 Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenterQbck(eos, tau - dt / 2., e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck);
    c->getPrimVarHCenter(eos, tau - dt / 2., e, p, nb, nq, ns, vx, vy, vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck);  // TODO Cartesian
 if (e <= 0.) return;
    uuu[0] = 1. / sqrt(1. - vx * vx - vy * vy - vz * vz);
    uuu[1] = uuu[0] * vx;
    uuu[2] = uuu[0] * vy;
    uuu[3] = uuu[0] * vz;

 k[T_] = -c->getpiH(3, 3) + c->getPiH() * (-1.0 - uuu[3] * uuu[3]);
 k[X_] = 0.;
 k[Y_] = 0.;
 k[Z_] = -(c->getpiH(0, 3) + c->getPiH() * uuu[0] * uuu[3]);
 for (int i = 0; i < 4; i++) k[i] *= dt;
 c->addFlux(k[T_], k[X_], k[Y_], k[Z_], 0., 0., 0.);
 #endif
}

// for the procedure below, the following approximations are used:
// dv/d(tau) = v^{t+dt}_ideal - v^{t}
// dv/dx_i ~ v^{x+dx}-v{x-dx},
// which makes sense after non-viscous step
void Hydro::NSquant(int ix, int iy, int iz, double d_pi[4][4], double &d_Pi, double d_dmu[4][4], double dmu_bck[4][4], double &d_du, double &du_bck) {
    const double VMIN = 1e-2;
    const double UDIFF = 3.0;
    double d_e0, d_e1, d_p, d_nb, d_nq, d_ns, d_vx1, d_vy1, d_vz1, d_vx0, d_vy0, d_vz0, d_vxH, d_vyH, d_vzH;//fluctuations
    double e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck, e1_bck, vx1_bck, vy1_bck, vz1_bck, vxH_bck, vyH_bck, vzH_bck;//background
    double d_ut0, d_ux0, d_uy0, d_uz0, d_ut1, d_ux1, d_uy1, d_uz1;//fluctuations
    double ut0_bck, ux0_bck, uy0_bck, uz0_bck, ut1_bck, ux1_bck, uy1_bck, uz1_bck;//background
    //    double dmu [4][4] ; // \partial_\mu u^\nu matrix
    // coordinates: 0=tau, 1=x, 2=y, 3=eta
    double d_Z[4][4][4][4];  // Z[mu][nu][lambda][rho] - fluctuating part of Z
    double Z_bck[4][4][4][4]; // Z matrix for background - background part of Z
    double d_uuu[4];         // fluctuation of the 4-velocity
    double uuu_bck[4];        // the background 4-velocity
    double gmunu[4][4] = {{1, 0, 0, 0},
        {0, -1, 0, 0},
        {0, 0, -1, 0},
        {0, 0, 0, -1}};  // omit 1/tau^2 in g^{eta,eta}
    Cell *c = f->getCell(ix, iy, iz);
    double dx = f->getDx(), dy = f->getDy(), dz = f->getDz();
    // check if the cell is next to vacuum from +-x, +-y side:
    if (f->getCell(ix + 1, iy, iz)->getMaxM() <= 0.9 ||
        f->getCell(ix, iy + 1, iz)->getMaxM() <= 0.9 ||
        f->getCell(ix - 1, iy, iz)->getMaxM() <= 0.9 ||
        f->getCell(ix, iy - 1, iz)->getMaxM() <= 0.9 ||
        f->getCell(ix + 1, iy + 1, iz)->getMaxM() <= 0.9 ||
        f->getCell(ix + 1, iy - 1, iz)->getMaxM() <= 0.9 ||
        f->getCell(ix - 1, iy + 1, iz)->getMaxM() <= 0.9 ||
        f->getCell(ix - 1, iy - 1, iz)->getMaxM() <= 0.9) {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                d_pi[i][j] = 0.;
                d_dmu[i][j] = 0.;
                dmu_bck[i][j] = 0.;
            }
        d_Pi = d_du = du_bck = 0.;
        return;
    }
    // calculation of \partial_\mu u^\nu matrix
    // mu=first index, nu=second index
    // centered differences with respect to the values at (it+1/2, ix, iy, iz)
    // d_tau u^\mu
    
#ifdef CARTESIAN
    c->getPrimVarPrevQbck(eos, 1.0, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
    c->getPrimVarPrev(eos, 1.0, d_e0, d_p, d_nb, d_nq, d_ns, d_vx0, d_vy0, d_vz0, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
    c->getPrimVarQbck(eos, 1.0, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);
    c->getPrimVar(eos, 1.0, d_e1, d_p, d_nb, d_nq, d_ns, d_vx1, d_vy1, d_vz1, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);
    c->getPrimVarHCenterQbck(eos, 1.0, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vxH_bck, vyH_bck, vzH_bck);
    c->getPrimVarHCenter(eos, 1.0, d_e1, d_p, d_nb, d_nq, d_ns, d_vxH, d_vyH, d_vzH, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vxH_bck, vyH_bck, vzH_bck);
    double tauPlusHalf = 1.0;
#else
    c->getPrimVarPrevQbck(eos, tau - dt, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
    c->getPrimVarPrev(eos, tau - dt, d_e0, d_p, d_nb, d_nq, d_ns, d_vx0, d_vy0, d_vz0, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
    c->getPrimVarQbck(eos, tau, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);
    c->getPrimVar(eos, tau, d_e1, d_p, d_nb, d_nq, d_ns, d_vx1, d_vy1, d_vz1, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);
    c->getPrimVarHCenterQbck(eos, tau - 0.5 * dt, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vxH_bck, vyH_bck, vzH_bck);
    c->getPrimVarHCenter(eos, tau - 0.5 * dt, d_e1, d_p, d_nb, d_nq, d_ns, d_vxH, d_vyH, d_vzH, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vxH_bck, vyH_bck, vzH_bck);
    double tauPlusHalf = tau + 0.5 * dt;
#endif
    
    //############## get transport coefficients
    double T_bck, mub, muq, mus;
    double etaS, zetaS; // constants
    double h = 0.01; // arbitrary shift for derivative
    double d_s; // definition of fluctiation in enthropy
    double s_bck = eos->s(e1_bck, nb_bck, nq_bck, ns_bck);  // background entropy density in the current cell
    double dS = (eos->s(e1_bck + h, nb_bck, nq_bck, ns_bck) - s_bck ) / h; // derivative of enthropy with respect to the energy density at e_0
    eos->eos(e1_bck, nb_bck, nq_bck, ns_bck, T_bck, mub, muq, mus, p0_bck); // this gives (not just) temperature at the background energy density
    d_s = dS * d_e1; // fluctuation in enthropy
    trcoeff->getEta(e1_bck, nb_bck, T_bck, etaS, zetaS); // gets etaS, zetaS at background energy density
    //##############
    // if(e1<0.00004) s=0. ; // negative pressure due to pi^zz for small e
    
    //background
    ut0_bck = 1./ sqrt(1. - vx0_bck * vx0_bck - vy0_bck * vy0_bck - vz0_bck * vz0_bck);
    ux0_bck = ut0_bck * vx0_bck;
    uy0_bck = ut0_bck * vy0_bck;
    uz0_bck = ut0_bck * vz0_bck;
    ut1_bck = 1./ sqrt(1. - vx1_bck * vx1_bck - vy1_bck * vy1_bck - vz1_bck * vz1_bck);
    ux1_bck = ut1_bck * vx1_bck;
    uy1_bck = ut1_bck * vy1_bck;
    uz1_bck = ut1_bck * vz1_bck;
    uuu_bck[0] = 1./ sqrt(1. - vxH_bck * vxH_bck - vyH_bck * vyH_bck - vzH_bck * vzH_bck);
    uuu_bck[1] = uuu_bck[0] * vxH_bck;
    uuu_bck[2] = uuu_bck[0] * vyH_bck;
    uuu_bck[3] = uuu_bck[0] * vzH_bck;
    
    dmu_bck[0][0] = (ut1_bck - ut0_bck) / dt;
    dmu_bck[0][1] = (ux1_bck - ux0_bck) / dt;
    dmu_bck[0][2] = (uy1_bck - uy0_bck) / dt;
    dmu_bck[0][3] = (uz1_bck - uz0_bck) / dt;
    
    //fluctuations
    d_ut0 = 0.;
    d_ux0 = ut0_bck * d_vx0;
    d_uy0 = ut0_bck * d_vy0;
    d_uz0 = ut0_bck * d_vz0;
    d_ut1 = 0.;
    d_ux1 = ut1_bck * d_vx1;
    d_uy1 = ut1_bck * d_vy1;
    d_uz1 = ut1_bck * d_vz1;
    d_uuu[0] = 0.;
    d_uuu[1] = uuu_bck[0] * d_vxH;
    d_uuu[2] = uuu_bck[0] * d_vyH;
    d_uuu[3] = uuu_bck[0] * d_vzH;
    
    d_dmu[0][0] = (d_ut1 - d_ut0) / dt;
    d_dmu[0][1] = (d_ux1 - d_ux0) / dt;
    d_dmu[0][2] = (d_uy1 - d_uy0) / dt;
    d_dmu[0][3] = (d_uz1 - d_uz0) / dt;
    
    if ((d_e1+e1_bck <= 0. || d_e0+e0_bck <= 0.) || (e1_bck <= 0. || e0_bck <= 0.)) {  // modified condition - matter-vacuum
        d_dmu[0][0] = d_dmu[0][1] = d_dmu[0][2] = d_dmu[0][3] = 0.;
        dmu_bck[0][0] = dmu_bck[0][1] = dmu_bck[0][2] = dmu_bck[0][3] = 0.;
    }
    // d_x u^\mu
    f->getCell(ix + 1, iy, iz)->getPrimVarHCenterQbck(eos, tau, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);//0 i 1 maji stejne v_const
    f->getCell(ix + 1, iy, iz)->getPrimVarHCenter(eos, tau, d_e1, d_p, d_nb, d_nq, d_ns, d_vx1, d_vy1, d_vz1, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);//0 i 1 maji stejne v_const
    f->getCell(ix - 1, iy, iz)->getPrimVarHCenterQbck(eos, tau, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
    f->getCell(ix - 1, iy, iz)->getPrimVarHCenter(eos, tau, d_e0, d_p, d_nb, d_nq, d_ns, d_vx0, d_vy0, d_vz0, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
    
    if (d_e1+e1_bck > 0. && d_e0+e0_bck > 0.) {//modified condition
        // background
        ut0_bck = 1./ sqrt(1. - vx0_bck * vx0_bck - vy0_bck * vy0_bck - vz0_bck * vz0_bck);
        ux0_bck = ut0_bck * vx0_bck;
        uy0_bck = ut0_bck * vy0_bck;
        uz0_bck = ut0_bck * vz0_bck;
        ut1_bck = 1./ sqrt(1. - vx1_bck * vx1_bck - vy1_bck * vy1_bck - vz1_bck * vz1_bck);
        ux1_bck = ut1_bck * vx1_bck;
        uy1_bck = ut1_bck * vy1_bck;
        uz1_bck = ut1_bck * vz1_bck;
        
        dmu_bck[1][0] = 0.5 * (ut1_bck - ut0_bck) / dx;
        dmu_bck[1][1] = 0.5 * (ux1_bck - ux0_bck) / dx;
        dmu_bck[1][2] = 0.5 * (uy1_bck - uy0_bck) / dx;
        dmu_bck[1][3] = 0.5 * (uz1_bck - uz0_bck) / dx;
        // fluctuations
        d_ut0 = 0.;
        d_ux0 = ut0_bck * d_vx0;
        d_uy0 = ut0_bck * d_vy0;
        d_uz0 = ut0_bck * d_vz0;
        d_ut1 = 0.;
        d_ux1 = ut1_bck * d_vx1;
        d_uy1 = ut1_bck * d_vy1;
        d_uz1 = ut1_bck * d_vz1;
        
        d_dmu[1][0] = 0.5 * (d_ut1 - d_ut0) / dx;
        d_dmu[1][1] = 0.5 * (d_ux1 - d_ux0) / dx;
        d_dmu[1][2] = 0.5 * (d_uy1 - d_uy0) / dx;
        d_dmu[1][3] = 0.5 * (d_uz1 - d_uz0) / dx;

         } else {  // matter-vacuum
          d_dmu[1][0] = d_dmu[1][1] = d_dmu[1][2] = d_dmu[1][3] = 0.;
          dmu_bck[1][0] = dmu_bck[1][1] = dmu_bck[1][2] = dmu_bck[1][3] = 0.;

        }
        if (fabs(d_dmu[1][3]) > 1e+10)
            cout << "dmu[1][3]:  " << d_uz1 << "  " << d_uz0 << "  " << d_uuu[3] << endl;
        // d_y u^\mu
        f->getCell(ix, iy + 1, iz)->getPrimVarHCenterQbck(eos, tau, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);
        f->getCell(ix, iy + 1, iz)->getPrimVarHCenter(eos, tau, d_e1, d_p, d_nb, d_nq, d_ns, d_vx1, d_vy1, d_vz1, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);
        f->getCell(ix, iy - 1, iz)->getPrimVarHCenterQbck(eos, tau, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
        f->getCell(ix, iy - 1, iz)->getPrimVarHCenter(eos, tau, d_e0, d_p, d_nb, d_nq, d_ns, d_vx0, d_vy0, d_vz0, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
    
        if (d_e1+e1_bck > 0. && d_e0+e0_bck > 0.) {// modified condition
            // background
            ut0_bck = 1./ sqrt(1. - vx0_bck * vx0_bck - vy0_bck * vy0_bck - vz0_bck * vz0_bck);
            ux0_bck = ut0_bck * vx0_bck;
            uy0_bck = ut0_bck * vy0_bck;
            uz0_bck = ut0_bck * vz0_bck;
            ut1_bck = 1./ sqrt(1. - vx1_bck * vx1_bck - vy1_bck * vy1_bck - vz1_bck * vz1_bck);
            ux1_bck = ut1_bck * vx1_bck;
            uy1_bck = ut1_bck * vy1_bck;
            uz1_bck = ut1_bck * vz1_bck;
            
            dmu_bck[2][0] = 0.5 * (ut1_bck - ut0_bck) / dy;
            dmu_bck[2][1] = 0.5 * (ux1_bck - ux0_bck) / dy;
            dmu_bck[2][2] = 0.5 * (uy1_bck - uy0_bck) / dy;
            dmu_bck[2][3] = 0.5 * (uz1_bck - uz0_bck) / dy;
            //fluctuations
            d_ut0 = 0.;
            d_ux0 = ut0_bck * d_vx0;
            d_uy0 = ut0_bck * d_vy0;
            d_uz0 = ut0_bck * d_vz0;
            d_ut1 = 0.;
            d_ux1 = ut1_bck * d_vx1;
            d_uy1 = ut1_bck * d_vy1;
            d_uz1 = ut1_bck * d_vz1;
            
            d_dmu[2][0] = 0.5 * (d_ut1 - d_ut0) / dy;
            d_dmu[2][1] = 0.5 * (d_ux1 - d_ux0) / dy;
            d_dmu[2][2] = 0.5 * (d_uy1 - d_uy0) / dy;
            d_dmu[2][3] = 0.5 * (d_uz1 - d_uz0) / dy;
            
        } else {  // matter-vacuum
            d_dmu[2][0] = d_dmu[2][1] = d_dmu[2][2] = d_dmu[2][3] = 0.;
            dmu_bck[2][0] = dmu_bck[2][1] = dmu_bck[2][2] = dmu_bck[2][3] = 0.;

        }
        // d_z u^\mu
        f->getCell(ix, iy, iz + 1)->getPrimVarHCenterQbck(eos, tau, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);
        f->getCell(ix, iy, iz + 1)->getPrimVarHCenter(eos, tau, d_e1, d_p, d_nb, d_nq, d_ns, d_vx1, d_vy1, d_vz1, e1_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx1_bck, vy1_bck, vz1_bck);
        f->getCell(ix, iy, iz - 1)->getPrimVarHCenterQbck(eos, tau, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
        f->getCell(ix, iy, iz - 1)->getPrimVarHCenter(eos, tau, d_e0, d_p, d_nb, d_nq, d_ns, d_vx0, d_vy0, d_vz0, e0_bck, p0_bck, nb_bck, nq_bck, ns_bck, vx0_bck, vy0_bck, vz0_bck);
    
        if (d_e1+e1_bck > 0. && d_e0+e0_bck > 0.) {// modified condition
            // background
            ut0_bck = 1./ sqrt(1. - vx0_bck * vx0_bck - vy0_bck * vy0_bck - vz0_bck * vz0_bck);
            ux0_bck = ut0_bck * vx0_bck;
            uy0_bck = ut0_bck * vy0_bck;
            uz0_bck = ut0_bck * vz0_bck;
            ut1_bck = 1./ sqrt(1. - vx1_bck * vx1_bck - vy1_bck * vy1_bck - vz1_bck * vz1_bck);
            ux1_bck = ut1_bck * vx1_bck;
            uy1_bck = ut1_bck * vy1_bck;
            uz1_bck = ut1_bck * vz1_bck;
            
            dmu_bck[3][0] = 0.5 * (ut1_bck - ut0_bck) / dz / tauPlusHalf;
            dmu_bck[3][1] = 0.5 * (ux1_bck - ux0_bck) / dz / tauPlusHalf;
            dmu_bck[3][2] = 0.5 * (uy1_bck - uy0_bck) / dz / tauPlusHalf;
            dmu_bck[3][3] = 0.5 * (uz1_bck - uz0_bck) / dz / tauPlusHalf;
            // fluctuations
            d_ut0 = 0.;
            d_ux0 = ut0_bck * d_vx0;
            d_uy0 = ut0_bck * d_vy0;
            d_uz0 = ut0_bck * d_vz0;
            d_ut1 = 0.;
            d_ux1 = ut1_bck * d_vx1;
            d_uy1 = ut1_bck * d_vy1;
            d_uz1 = ut1_bck * d_vz1;
            
            d_dmu[3][0] = 0.5 * (d_ut1 - d_ut0) / dz / tauPlusHalf;
            d_dmu[3][1] = 0.5 * (d_ux1 - d_ux0) / dz / tauPlusHalf;
            d_dmu[3][2] = 0.5 * (d_uy1 - d_uy0) / dz / tauPlusHalf;
            d_dmu[3][3] = 0.5 * (d_uz1 - d_uz0) / dz / tauPlusHalf;
            
        } else {  // matter-vacuum
            d_dmu[3][0] = d_dmu[3][1] = d_dmu[3][2] = d_dmu[3][3] = 0.;
            dmu_bck[3][0] = dmu_bck[3][1] = dmu_bck[3][2] = dmu_bck[3][3] = 0.;
        }
        // additional terms from Christoffel symbols :)
        // are not considered in this code for Cartesian coordinates - they are 0 anyway
#ifndef CARTESIAN
        d_dmu[3][0] += d_uuu[3] / (tau - 0.5 * dt);
        d_dmu[3][3] += d_uuu[0] / (tau - 0.5 * dt);
        dmu_bck[3][0] += uuu_bck[3] / (tau - 0.5 * dt);
        dmu_bck[3][3] += uuu_bck[0] / (tau - 0.5 * dt);
#endif

        // introduction of Z_bck[mu][nu][lambda][rho]
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    for (int l = 0; l < 4; l++) Z_bck[i][j][k][l] = 0.0;
        // filling Z_bck matrix - background
        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
                for (int lam = 0; lam < 4; lam++)
                    for (int rho = 0; rho < 4; rho++) {
                        if (nu == rho)
                            Z_bck[mu][nu][lam][rho] += 0.5 * (gmunu[mu][lam] - uuu_bck[mu] * uuu_bck[lam]);
                        
                        if (mu == rho)
                            Z_bck[mu][nu][lam][rho] += 0.5 * (gmunu[nu][lam] - uuu_bck[nu] * uuu_bck[lam]);
                        
                        if (lam == rho)
                            Z_bck[mu][nu][lam][rho] -= (gmunu[mu][nu] - uuu_bck[mu] * uuu_bck[nu]) / 3.0;
                        
                    }
    
    // introduction of Z[mu][nu][lambda][rho]
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                for (int l = 0; l < 4; l++) d_Z[i][j][k][l] = 0.0;
    // filling Z matrix - perturbations
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            for (int lam = 0; lam < 4; lam++)
                for (int rho = 0; rho < 4; rho++) {
                    if (nu == rho)
                        d_Z[mu][nu][lam][rho] += 0.5 * (d_uuu[mu] * uuu_bck[lam] + uuu_bck[mu] * d_uuu[lam]);
                    if (mu == rho)
                        d_Z[mu][nu][lam][rho] += 0.5 * (d_uuu[nu] * uuu_bck[lam] + uuu_bck[nu] * d_uuu[lam]);
                    if (lam == rho)
                        d_Z[mu][nu][lam][rho] -= (d_uuu[mu] * uuu_bck[nu] + uuu_bck[mu] * d_uuu[nu]) / 3.0;
                }
    
        // calculating sigma[mu][nu] - pi[mu][nu] respectively
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                d_pi[i][j] = 0.0;
                for (int k = 0; k < 4; k++)
                    for (int l = 0; l < 4; l++) {
                        d_pi[i][j] += ( Z_bck[i][j][k][l] * d_dmu[k][l] - d_Z[i][j][k][l] * dmu_bck[k][l] ) * 2.0 * etaS * s_bck / 5.068
                                    + Z_bck[i][j][k][l] * dmu_bck[k][l] * 2.0 * etaS * d_s / 5.068 ;
                    }
            }

        d_Pi = -zetaS * s_bck * (d_dmu[0][0] + d_dmu[1][1] + d_dmu[2][2] + d_dmu[3][3]) / 5.068 - zetaS * d_s * (dmu_bck[0][0] + dmu_bck[1][1] + dmu_bck[2][2] + dmu_bck[3][3]) / 5.068;  // 5.068 -> fm^{-4} --> GeV/fm^3
        d_du = d_dmu[0][0] + d_dmu[1][1] + d_dmu[2][2] + d_dmu[3][3]; // trace of fluctuation
        du_bck = dmu_bck[0][0] + dmu_bck[1][1] + dmu_bck[2][2] + dmu_bck[3][3]; // trace of background
        //--------- debug part: NaN/inf check, trace check, diag check, transversality
        // check
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                if (d_pi[i][j] != 0. && fabs(1.0 - d_pi[j][i] / d_pi[i][j]) > 1e-10)
                    cout << "non-diag: " << d_pi[i][j] << "  " << d_pi[j][i] << endl;
                if (std::isinf(d_pi[i][j]) || std::isnan(d_pi[i][j])) {
                    cout << "hydro:NSquant: inf/nan i " << i << " j " << j << endl;
                    exit(1);
                }
            }
}

//void Hydro::setNSvalues() {
//double e, p, nb, nq, ns, vx, vy, vz, piNS[4][4], PiNS, dmu[4][4], du;
//double e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0;
// for (int ix = 0; ix < f->getNX(); ix++)
//  for (int iy = 0; iy < f->getNY(); iy++)
//   for (int iz = 0; iz < f->getNZ(); iz++) {
//    Cell *c = f->getCell(ix, iy, iz);
//    c->getPrimVarQbck(eos, tau, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
//    c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
//    if (e <= 0.) continue;
//    // NSquant(ix, iy, iz, piNS, PiNS, dmu, du) ;
//    //############## set NS values assuming initial zero flow + Bjorken z
//    // flow
//    double T, mub, muq, mus;
//    double etaS, zetaS;
//    double s = eos->s(e_0, nb, nq, ns);  // entropy density in the current cell//ZKONTROLOVAT_____________________::::::::::::::::::::::::::::::::::::::::::::
//    eos->eos(e_0, nb, nq, ns, T, mub, muq, mus, p_0);//??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
//    trcoeff->getEta(e_0,nb, T, etaS, zetaS);//Zkontrolovat_________________________________________________________::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//    for (int i = 0; i < 4; i++)
//     for (int j = 0; j < 4; j++) piNS[i][j] = 0.0;  // reset piNS
//    piNS[1][1] = piNS[2][2] = 2.0 / 3.0 * etaS * s / tau / 5.068;
//    piNS[3][3] = -2.0 * piNS[1][1];
//    PiNS = 0.0;
//    for (int i = 0; i < 4; i++)
//        for (int j = 0; j <= i; j++){
//            c->setpi(i, j, piNS[i][j]);
//        }
//    c->setPi(PiNS);
//   }
// cout << "setNS done\n";
//}

void Hydro::ISformal() {
 double d_e, d_p, d_nb, d_nq, d_ns, d_vx, d_vy, d_vz, mub, muq, mus; //fluctuations
 double e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck, T_bck; //background
 double d_piNS[4][4], d_sigNS[4][4], d_PiNS, d_dmu[4][4], d_du, d_pi[4][4], d_piH[4][4], d_Pi, d_PiH; // fluctuations
 double piNS_bck[4][4], sigNS_bck[4][4], PiNS_bck, dmu_bck[4][4], du_bck, pi_bck[4][4], piH_bck[4][4], Pi_bck, PiH_bck; // background
 double dpi_bck[4][4][4], dpiH_bck[4][4][4]; // derivatives of background pi^munu
 const double gmumu[4] = {1., -1., -1., -1.};
 #ifdef CARTESIAN
 double tauMinusHalf = 1.0;
 double tauMinusDt = 1.0;
 #else
 double tauMinusHalf = tau - 0.5 * dt;
 double tauMinusDt = tau - dt;
 #endif
 double dx = f->getDx(), dy = f->getDy(), dz = f->getDz();
 double gmunu[4][4] = {{1, 0, 0, 0},
                      {0, -1, 0, 0},
                      {0, 0, -1, 0},// tau factor not needed in Cartesian coordinates
                      {0, 0, 0, -1}};
 // loop #1 (relaxation+source terms)
    double pr[4][4];
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            pr[i][j] = 0.;
        }
    }
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenterQbck(eos, tauMinusHalf, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck);
    c->getPrimVarHCenter(eos, tauMinusHalf, d_e, d_p, d_nb, d_nq, d_ns, d_vx, d_vy,
                         d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck);  // instead of getPrimVar()
       
    if (d_e+e_bck < 0. || e_bck <= 0.) {// modified condition             // empty cell?
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
       c->setpiH0(i, j, 0.0);
       c->setpi0(i, j, 0.0);
      }
     c->setPiH0(0.0);
     c->setPi0(0.0);
    } else {  // non-empty cell
     // 1) relaxation(pi)+source(pi) terms for half-step
        double gamma = 1.0 / sqrt(1.0 - vx_bck * vx_bck - vy_bck * vy_bck - vz_bck * vz_bck);
        double u_bck[4];//background velocity
        u_bck[0] = gamma;
        u_bck[1] = u_bck[0] * vx_bck;
        u_bck[2] = u_bck[0] * vy_bck;
        u_bck[3] = u_bck[0] * vz_bck;
        
        double d_u[4];//fluctuation of velocity
        d_u[0] = 0.;
        d_u[1] = gamma * d_vx;
        d_u[2] = gamma * d_vy;
        d_u[3] = gamma * d_vz;
            
        // source term  + tau*delta_Q_i/delta_tau
        double d_flux[4];
        for (int i = 0; i < 4; i++){
            d_flux[i] = tauMinusDt * (c->getpi(0, i) + c->getPi() * ( d_u[0] * u_bck[i] + u_bck[0] * d_u[i] ) );
        }
        d_flux[0] += -tauMinusDt * c->getPi();
        c->addFlux(d_flux[0], d_flux[1], d_flux[2], d_flux[3], 0., 0., 0.);
        // now calculating viscous terms in NS limit
        NSquant(ix, iy, iz, d_piNS, d_PiNS, d_dmu, dmu_bck, d_du, du_bck);
        PiNS_bck = 0.; // setting the background NS variables 0 here - need to get it from the full hydro code when run together
        for (int i = 0; i < 4; i++)
            for (int j = 0; j <= i; j++) {
                piNS_bck[i][j]=0.;
            }
// derivative of background pi^munu
#ifdef CARTESIAN
    double tauPlusHalf_0 = 1.0;
#else
    double tauPlusHalf_0 = tau + 0.5 * dt;
#endif
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            dpi_bck[i][j][0] = (f->getCell(ix, iy, iz)->getpi_bck(i,j) - f->getCell(ix, iy, iz)->getpi_bck_prev(i,j)) / dt;
            dpi_bck[i][j][1] = 0.5 * (f->getCell(ix + 1, iy, iz)->getpi_bck(i,j) - f->getCell(ix - 1, iy, iz)->getpi_bck(i,j)) / dx;
            dpi_bck[i][j][2] = 0.5 * (f->getCell(ix, iy + 1, iz)->getpi_bck(i,j) - f->getCell(ix, iy - 1, iz)->getpi_bck(i,j)) / dy;
            dpi_bck[i][j][3] = 0.5 * (f->getCell(ix, iy, iz + 1)->getpi_bck(i,j) - f->getCell(ix, iy, iz - 1)->getpi_bck(i,j)) / dz / tauPlusHalf_0;
        }
    }
        
     double T1; // just for the derivative
     double h = 0.01; // shift for derivative
     eos->eos(e_bck, nb_bck, nq_bck, ns_bck, T_bck, mub, muq, mus, p_bck);
     eos->eos(e_bck* (1 + h), nb_bck, nq_bck, ns_bck, T1, mub, muq, mus, p_bck);
     double dT = (T1 - T_bck) / (e_bck *h) ; // derivative of temperature wrt energy density
     double etaS, zetaS; // constant ratio
     trcoeff->getEta(e_bck, nb_bck, T_bck, etaS, zetaS); // obtains eta and zeta
     const double s_bck = eos->s(e_bck, nb_bck, nq_bck, ns_bck); // background enthropy
     double dS = (eos->s(e_bck * (1 + h), nb_bck, nq_bck, ns_bck) - s_bck ) / (e_bck *h); // derivative of enthropy wtr to energy density
     double d_s = dS * d_e; // fluctuation of enthropy
     const double eta_bck = etaS * s_bck; // eta0 - background viscosity
//        cout << eta_bck << endl;
     const double d_eta = etaS * d_s; // fluctuation in viscosity
     double eta = eta_bck + d_eta; // whole shear viscosity
     // auxiliary variable sigmaNS = piNS / (2*eta),
     // mainly to protect against division by zero in the eta=0 case.
     for(int i=0; i<4; i++)
     for(int j=0; j<4; j++) {
      sigNS_bck[i][j] = 0.; // background sigma NS set to zero - need to get from the full code
      d_sigNS[i][j] = 0.5 * (d_piNS[i][j] * 5.068 - 2. * d_eta * sigNS_bck[i][j] ) / eta_bck ;
         if(eta<=0.0){
             d_sigNS[i][j] = 0.0;
             sigNS_bck[i][j] = 0.;
         }
     }
     //############# get relaxation times
     double taupi_bck, tauPi_bck; // the background values of relaxation times
     double d_taupi, d_tauPi; // fluctuations in relaxation times
     trcoeff->getTau(e_bck, d_e, d_nb, T_bck, dT, taupi_bck, d_taupi, tauPi_bck, d_tauPi);
     double deltapipi, taupipi, lambdapiPi, phi7, phi7_bck, d_phi7, delPiPi, lamPipi; // coefficients for source terms in relaxation equations
     trcoeff->getOther(e_bck, d_nb, d_nq, d_ns, deltapipi, taupipi, lambdapiPi, phi7);
     phi7_bck = phi7/taupi_bck;  // dividing by tau_pi here, to avoid NaNs when tau_pi==0
     d_phi7 = phi7 / (taupi_bck * taupi_bck) * d_taupi; // fluctuation in phi7 coeff - it is the only one which does not include taupi
     if(taupi_bck < 0.5 * dt) // modified condition - what about this condition?
      deltapipi = taupipi = lambdapiPi = phi7 = 0.0;
     trcoeff->getOtherBulk(e_bck, d_nb, d_nq, d_ns, delPiPi, lamPipi);
     if(tauPi_bck < 0.5 * dt)// modified condition
      delPiPi = lamPipi = 0.0;
     //#############
        //noise
        double xi[4][4];
        double delta[4][4];
        double Trxi = 0.;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                delta[i][j] = gmunu[i][j] - u_bck[i]*u_bck[j];
                xi[i][j] = 0.;
            }
        }
        double size_x = f->getDx();
        double size_y = f->getDy();
        double size_z = f->getDz();
        double volume = size_x * size_y * size_z;
        for (int i=0; i<4; i++) {
            for (int j=0; j<=i; j++) {
                std::random_device rd;
                std::mt19937 gen(rd());
                double mean = 0.;
                double sigma = sqrt(2 * eta_bck * T_bck * (delta[i][i] * delta[j][j] + delta[i][j] * delta[j][i])) / (sqrt(dt) * pow(volume, 1./2.)) * sqrt(0.197);
                std::normal_distribution<double> d(mean, sigma);
                xi[i][j] = d(gen);
                if (i!=j) {
                    xi[j][i] = xi[i][j];
                }
//                cout << xi[i][j] << "   " << i << "     " << j << endl;
                if (i==j) {
                    Trxi += xi[i][j];
                }
            }
        }
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                if (i==j && i!=0) {
                    xi[i][j] -= 1./3. * Trxi;
                }
            }
        }
        //######################
     double Delta[10]; // corresponds to background Delta
     // relaxation term, piH,PiH-->half-step
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
              Delta[index44(i, j)] = - u_bck[i] * u_bck[j];
       if (i == j) Delta[index44(i, j)] += gmumu[i];
          
#ifdef FORMAL_SOLUTION
       c->setpiH0(i, j, (c->getpi(i, j) - d_piNS[i][j]) *
                                exp(-dt / 2.0 / gamma / taupi_bck) +
                            d_piNS[i][j]);
       c->addpiH0(i, j, (c->getpi_bck(i, j) - piNS_bck[i][j]) * dt / 2.0  / gamma / (taupi_bck*taupi_bck) * d_taupi ); // this should be the same even for formal solution, I think
#else
          if(taupi_bck > 0.5 * dt){// modified condition - what the condition should be?
              c->setpiH0(i, j, c->getpi(i, j) -
                         (c->getpi(i, j) - d_piNS[i][j] - xi[i][j]) * dt / 2.0 / gamma / taupi_bck);
              c->addpiH0(i, j, (c->getpi_bck(i, j) - piNS_bck[i][j]) * dt / 2.0 / gamma / (taupi_bck * taupi_bck) * d_taupi); // source term from delta tau_pi
//              c->addpiH0(i, j, xi[i][j] * dt / 2.0 / gamma / taupi_bck);//adding the noise term
          }
          else
              c->setpiH0(i, j, d_piNS[i][j] + xi[i][j]);//adding the noise term
#endif
      }
#ifdef FORMAL_SOLUTION
     c->setPiH0((c->getPi() - d_PiNS) * exp(-dt / 2.0 / gamma / tauPi_bck) + d_PiNS);
     c->setPiH0( (c->getPi_bck() - PiNS_bck) * dt / 2.0 / gamma / (tauPi_bck * tauPi_bck) * d_tauPi );
#else
        if(tauPi_bck > 0.5 * dt){// modified condition
            c->setPiH0(c->getPi() - (c->getPi() - d_PiNS) * dt / 2.0 / gamma / tauPi_bck);
            c->setPiH0( (c->getPi_bck() - PiNS_bck) * dt / 2.0 / gamma / (tauPi_bck * tauPi_bck) * d_tauPi );
        }
    else
     c->setPiH0(d_PiNS);
#endif
     // sources from Christoffel symbols from \dot pi_munu - only in tau-eta coordinate frame
     #ifndef CARTESIAN
     double tau1 = tau - dt * 0.75;
     c->addpiH0(0, 0,
                -2. * vz * c->getpi(0, 3) / tau1 * dt / 2.);  // *gamma/gamma
     c->addpiH0(3, 3, -(2. * vz / tau1 * c->getpi(0, 3)) * dt / 2.);
     c->addpiH0(
         3, 0,
         -(vz / tau1 * c->getpi(0, 0) + vz / tau1 * c->getpi(3, 3)) * dt / 2.);
     c->addpiH0(1, 0, -vz / tau1 * c->getpi(1, 3) * dt / 2.);
     c->addpiH0(2, 0, -vz / tau1 * c->getpi(2, 3) * dt / 2.);
     c->addpiH0(3, 1, -(vz / tau1 * c->getpi(0, 1)) * dt / 2.);
     c->addpiH0(3, 2, -(vz / tau1 * c->getpi(0, 2)) * dt / 2.);
     #endif
     // source from full IS equations (see  draft for the description)
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
//        now transversality and cross terms
//          c->addpiH0(i, j, -xi[i][j]);
          c->addpiH0(i, j, (- deltapipi * ( c->getpi(i, j) * du_bck + c->getpi_bck(i,j) * d_du ) ) / gamma * 0.5 * dt ); // the 4/3 source term
          c->addpiH0(i, j, lambdapiPi * ( c->getPi() * sigNS_bck[i][j] + c->getPi_bck() * d_sigNS[i][j] ) / gamma * 0.5 * dt); // lambdapiPi
       for (int k = 0; k < 4; k++) {
//         parts of terms with one internal summation index
        c->addpiH0(i, j, phi7_bck * ( c->getpi_bck(i, k) * c->getpi(j, k) + c->getpi(i, k) * c->getpi_bck(j, k) ) * gmumu[k] / gamma * 0.5 * dt ); // first part of phi7
        c->addpiH0(i, j, taupipi * 0.5 * ( c->getpi_bck(i, k) * d_sigNS[j][k] + c->getpi(i, k) * sigNS_bck[j][k] + c->getpi_bck(j, k) * d_sigNS[i][k] + c->getpi(j, k) * sigNS_bck[i][k] ) * gmumu[k] / gamma * 0.5 * dt ); // first part of taupipi
           c->addpiH0(i,j, - d_u[k] * dpi_bck[i][j][k] / gamma * 0.5 * dt );
           c->addpiH0(i, j, d_phi7 * c->getpi_bck(i,k) * c->getpi_bck(j,k) * gmumu[k] / gamma * 0.5 * dt); // delta phi7 part
//         parts of terms with two internal summation indexes
        for (int l = 0; l < 4; l++){
            c->addpiH0(i, j, - (c->getpi(i, k) * u_bck[j] + c->getpi(j, k) * u_bck[i]) * u_bck[l] * dmu_bck[l][k] * gmumu[k] / gamma * 0.5 * dt
                             - (c->getpi_bck(i, k) * u_bck[j] + c->getpi_bck(j, k) * u_bck[i]) * u_bck[l] * d_dmu[l][k] * gmumu[k] / gamma * 0.5 * dt
                             - d_u[k] * (u_bck[j] * u_bck[l] * dpi_bck[k][i][l] + u_bck[i] * u_bck[l] * dpi_bck[k][j][l]) * gmumu[k] / gamma * 0.5 * dt
                             - (u_bck[j] * c->getpi_bck(i, k) + u_bck[i] * c->getpi_bck(j, k)) * d_u[l] * dmu_bck[l][k] * gmumu[k] / gamma * 0.5 * dt );
            c->addpiH0(i, j, phi7_bck * ( ( c->getpi_bck(i, l) * u_bck[j] + c->getpi_bck(j, l) * u_bck[i] ) * u_bck[k] * c->getpi(k, l)
                                      - 2./3. * Delta[index44(i,j)] * (c->getpi_bck(k, l) * c->getpi(k, l) ) ) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt ); // second part of phi7
            c->addpiH0(i, j, phi7_bck * ( - ( u_bck[j] * c->getpi_bck(i, l) + u_bck[i] * c->getpi_bck(j, l) ) * d_u[k] * c->getpi_bck(k, l)
                                      + 1./3 * ( d_u[i] * u_bck[j] + u_bck[i] * d_u[j] ) * c->getpi_bck(k, l) * c->getpi_bck(k, l) ) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt ); // delta part of brackets for phi7

            c->addpiH0(i, j, taupipi * ( 0.5 * ( c->getpi_bck(i, l) * u_bck[j] + c->getpi_bck(j, l) * u_bck[i] ) * u_bck[k] * d_sigNS[k][l]
                                       + 0.5 * ( sigNS_bck[j][l] * u_bck[i] + sigNS_bck[i][l] * u_bck[j] ) * u_bck[k] * c->getpi(k, l)
                                       - 1./3. * Delta[index44(i,j)] * (c->getpi_bck(k, l) * d_sigNS[k][l] + c->getpi(k, l) * sigNS_bck[k][l] ) ) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt ); // second part of taupipi
            c->addpiH0(i, j, taupipi * ( - 0.5 * ( u_bck[j] * c->getpi_bck(i, l) + u_bck[i] * c->getpi_bck(j, l) ) * d_u[k] * sigNS_bck[k][l]
                                         - 0.5 * ( u_bck[i] * sigNS_bck[j][l] + u_bck[j] * sigNS_bck[i][l] ) * d_u[k] * c->getpi_bck(k, l)
                                         + 1./3 * ( d_u[i] * u_bck[j] + u_bck[i] * d_u[j] ) * c->getpi_bck(k, l) * sigNS_bck[k][l] ) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt ); // delta part of brackets for taupipi
            
            c->addPiH0(lamPipi * (c->getpi(k, l) * sigNS_bck[k][l] + c->getpi_bck(k, l) * d_sigNS[k][l] ) / gamma * 0.5 * dt);
            for (int r = 0; r < 4; r++) {//3 summation indeces
                c->addpiH0(i, j, 2./3 * (gmunu[i][j] + 2 * u_bck[i] * u_bck[j]) * u_bck[k] * c->getpi(k, r) * u_bck[l] * dmu_bck[l][r] * gmumu[k] * gmumu[r] / gamma * 0.5 * dt );
                c->addpiH0(i, j, ( 1./2 * ( ( d_u[j] * u_bck[r] + u_bck[j] * d_u[r] ) * ( gmunu[i][k] + u_bck[i] * u_bck[k] ) + ( d_u[i] * u_bck[k] + u_bck[i] * d_u[k] ) * ( gmunu[j][r] + u_bck[j] * u_bck[r] ) )
                                 + 1./2 * ( ( d_u[i] * u_bck[r] + u_bck[i] * d_u[r] ) * ( gmunu[j][k] + u_bck[j] * u_bck[k] ) + ( d_u[j] * u_bck[k] + u_bck[j] * d_u[k] ) * ( gmunu[i][r] + u_bck[i] * u_bck[r] ) )
                                 - 1./3 * ( ( d_u[k] * u_bck[r] + u_bck[k] * d_u[r] ) * ( gmunu[i][j] + u_bck[i] * u_bck[j] ) + ( d_u[i] * u_bck[j] + u_bck[i] * d_u[j] ) * ( gmunu[k][r] + u_bck[k] * u_bck[r] ) )
                                 ) * u_bck[l] * dpi_bck[k][r][l] * gmumu[k] * gmumu[r] / gamma * 0.5 * dt );
            }
        }
       }
      }
     c->addPiH0(-delPiPi * ( c->getPi() * du_bck + c->getPi_bck() * d_du ) / gamma * 0.5 * dt);
        
     for(int i = 0; i < 4; i++){//calculation of derivative of background pi^munu at the halfstep
         for(int j = 0; j < 4; j++){
             dpiH_bck[i][j][0] = (f->getCell(ix, iy, iz)->getpiH0_bck(i,j) - f->getCell(ix, iy, iz)->getpiH0_bck_prev(i,j)) / dt;
             dpiH_bck[i][j][1] = 0.5 * (f->getCell(ix + 1, iy, iz)->getpiH0_bck(i,j) - f->getCell(ix - 1, iy, iz)->getpiH0_bck(i,j)) / dx;
             dpiH_bck[i][j][2] = 0.5 * (f->getCell(ix, iy + 1, iz)->getpiH0_bck(i,j) - f->getCell(ix, iy - 1, iz)->getpiH0_bck(i,j)) / dy;
             dpiH_bck[i][j][3] = 0.5 * (f->getCell(ix, iy, iz + 1)->getpiH0_bck(i,j) - f->getCell(ix, iy, iz - 1)->getpiH0_bck(i,j)) / dz / tauPlusHalf_0;
         }
     }
        
     // 1) relaxation(piH)+source(piH) terms for full-step
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
#ifdef FORMAL_SOLUTION
       c->setpi0(i, j,
                 (c->getpi(i, j) - d_piNS[i][j]) * exp(-dt / gamma / taupi_bck) +
                     d_piNS[i][j]);
       c->addpi0(i, j, (c->getpiH0_bck(i, j) - piNS_bck[i][j]) * dt / gamma / (taupi_bck*taupi_bck) * d_taupi );//H0 or not?

#else
      if(taupi_bck > 0.5 * dt){// modified condition
              c->setpi0(i, j, c->getpi(i, j) -
                        (c->getpiH0(i, j) - d_piNS[i][j] - xi[i][j]) * dt / gamma / taupi_bck);
              c->addpi0(i, j, (c->getpiH0_bck(i, j) - piNS_bck[i][j]) * dt / gamma / (taupi_bck*taupi_bck) * d_taupi );
//              c->addpi0(i, j, xi[i][j] * dt / gamma / taupi_bck);//adding the noise term
          }
      else
          c->setpi0(i, j, d_piNS[i][j] + xi[i][j]);//adding the noise term
#endif
      }
        
#ifdef FORMAL_SOLUTION
     c->setPi0((c->getPi() - d_PiNS) * exp(-dt / gamma / tauPi_bck) + d_PiNS);
     c->setPi0( (c->getPiH0_bck() - PiNS_bck) * dt / gamma / (tauPi_bck * tauPi_bck) * d_tauPi );
#else
        if(tauPi_bck > 0.5 * dt){// modified condition
            c->setPi0(c->getPi() - (c->getPiH0() - d_PiNS) * dt / gamma / tauPi_bck);
            c->setPi0( (c->getPiH0_bck() - PiNS_bck) * dt / gamma / (tauPi_bck * tauPi_bck) * d_tauPi );
        }
        else
            c->setPi0(d_PiNS);
#endif
     #ifndef CARTESIAN
     tau1 = tau - dt * 0.5;
     c->addpi0(0, 0, -2. * vz / tau1 * c->getpiH0(0, 3) * dt);  // *gamma/gamma
     c->addpi0(3, 3, -(2. * vz / tau1 * c->getpiH0(0, 3)) * dt);
     c->addpi0(
         3, 0,
         -(vz / tau1 * c->getpiH0(0, 0) + vz / tau1 * c->getpiH0(3, 3)) * dt);
     c->addpi0(1, 0, -vz / tau1 * c->getpiH0(1, 3) * dt);
     c->addpi0(2, 0, -vz / tau1 * c->getpiH0(2, 3) * dt);
     c->addpi0(3, 1, -(vz / tau1 * c->getpiH0(0, 1)) * dt);
     c->addpi0(3, 2, -(vz / tau1 * c->getpiH0(0, 2)) * dt);
     #endif
     // source from full IS equations (see draft for the description)
        for (int i = 0; i < 4; i++)
         for (int j = 0; j <= i; j++) {
    //        now transversality and cross terms
             c->addpi0(i, j, (- deltapipi * ( c->getpiH0(i, j) * du_bck + c->getpiH0_bck(i,j) * d_du ) ) / gamma * dt );
             c->addpi0(i, j, lambdapiPi * ( c->getPiH0() * sigNS_bck[i][j] + c->getPiH0_bck() * d_sigNS[i][j] ) / gamma * dt);
          for (int k = 0; k < 4; k++) {
    //         parts of terms with one internal summation index
              c->addpi0(i, j, phi7_bck * ( c->getpiH0_bck(i, k) * c->getpiH0(j, k) + c->getpiH0(i, k) * c->getpiH0_bck(j, k) ) * gmumu[k] / gamma * dt ); // first part of phi7
              c->addpi0(i, j, taupipi * 0.5 * ( c->getpiH0_bck(i, k) * d_sigNS[j][k] + c->getpiH0(i, k) * sigNS_bck[j][k] + c->getpiH0_bck(j, k) * d_sigNS[i][k] + c->getpiH0(j, k) * sigNS_bck[i][k] ) * gmumu[k] / gamma * dt ); // first part of taupipi
              c->addpi0(i, j, - d_u[k] * dpiH_bck[i][j][k] / gamma * dt );
              c->addpi0(i, j, d_phi7 * c->getpiH0_bck(i,k) * c->getpiH0_bck(j,k) * gmumu[k] / gamma * dt); // delta part of phi7
              
    //         parts of terms with two internal summation indexes
           for (int l = 0; l < 4; l++){
               c->addpi0(i, j, - (c->getpiH0(i, k) * u_bck[j] + c->getpiH0(j, k) * u_bck[i]) * u_bck[l] * dmu_bck[l][k] * gmumu[k] / gamma * dt
                            - (c->getpiH0_bck(i, k) * u_bck[j] + c->getpiH0_bck(j, k) * u_bck[i]) * u_bck[l] * d_dmu[l][k] * gmumu[k] / gamma * dt
                            - d_u[k] * (u_bck[j] * u_bck[l] * dpiH_bck[k][i][l] + u_bck[i] * u_bck[l] * dpiH_bck[k][j][l]) * gmumu[k] / gamma * dt
                            - (u_bck[j] * c->getpiH0_bck(i, k) + u_bck[i] * c->getpiH0_bck(j, k)) * d_u[l] * dmu_bck[l][k] * gmumu[k] / gamma * dt );
               c->addpi0(i, j, phi7_bck * ( ( c->getpiH0_bck(i, l) * u_bck[j] + c->getpiH0_bck(j, l) * u_bck[i] ) * u_bck[k] * c->getpiH0(k, l)
                                        - 2./3. * Delta[index44(i,j)] * (c->getpiH0_bck(k, l) * c->getpiH0(k, l) ) ) * gmumu[k] * gmumu[l] / gamma * dt ); // second part of phi7
               c->addpi0(i, j, phi7_bck * ( - ( u_bck[j] * c->getpiH0_bck(i, l) + u_bck[i] * c->getpiH0_bck(j, l) ) * d_u[k] * c->getpiH0_bck(k, l)
                                        + 1./3 * ( d_u[i] * u_bck[j] + u_bck[i] * d_u[j] ) * c->getpiH0_bck(k, l) * c->getpiH0_bck(k, l) ) * gmumu[k] * gmumu[l] / gamma * dt ); // delta part of brackets for phi7

               c->addpi0(i, j, taupipi * ( 0.5 * ( c->getpiH0_bck(i, l) * u_bck[j] + c->getpiH0_bck(j, l) * u_bck[i] ) * u_bck[k] * d_sigNS[k][l]
                                        + 0.5 * ( sigNS_bck[j][l] * u_bck[i] + sigNS_bck[i][l] * u_bck[j] ) * u_bck[k] * c->getpiH0(k, l)
                                        - 1./3. * Delta[index44(i,j)] * (c->getpiH0_bck(k, l) * d_sigNS[k][l] + c->getpiH0(k, l) * sigNS_bck[k][l] ) ) * gmumu[k] * gmumu[l] / gamma * dt ); // second part of taupipi
               c->addpi0(i, j, taupipi * ( - 0.5 * ( u_bck[j] * c->getpiH0_bck(i, l) + u_bck[i] * c->getpiH0_bck(j, l) ) * d_u[k] * sigNS_bck[k][l]
                                        - 0.5 * ( u_bck[i] * sigNS_bck[j][l] + u_bck[j] * sigNS_bck[i][l] ) * d_u[k] * c->getpiH0_bck(k, l)
                                        + 1./3 * ( d_u[i] * u_bck[j] + u_bck[i] * d_u[j] ) * c->getpiH0_bck(k, l) * sigNS_bck[k][l] ) * gmumu[k] * gmumu[l] / gamma * dt ); // delta part of brackets for taupipi

               c->addpi0(i, j, - 1. / 3. * Delta[index44(i,j)] * c->getpiH0_bck(k, l) * d_phi7 * c->getpiH0_bck(k, l) / gamma * 0.5 * dt);
               c->addPi0(lamPipi * (c->getpiH0(k, l) * sigNS_bck[k][l] + c->getpiH0_bck(k, l) * d_sigNS[k][l] ) / gamma * dt);
               
               for (int r = 0; r < 4; r++) {//3 internal indices
                   c->addpi0(i, j, 2./3 * (gmunu[i][j] + 2 * u_bck[i] * u_bck[j]) * u_bck[k] * c->getpiH0(k, r) * u_bck[l] * dmu_bck[l][r] * gmumu[k] * gmumu[r] / gamma * dt );//?????????????????????
                   c->addpi0(i, j, ( 1./2 * ( ( d_u[j] * u_bck[r] + u_bck[j] * d_u[r] ) * ( gmunu[i][k] + u_bck[i] * u_bck[k] ) + ( d_u[i] * u_bck[k] + u_bck[i] * d_u[k] ) * ( gmunu[j][r] + u_bck[j] * u_bck[r] ) )
                                + 1./2 * ( ( d_u[i] * u_bck[r] + u_bck[i] * d_u[r] ) * ( gmunu[j][k] + u_bck[j] * u_bck[k] ) + ( d_u[j] * u_bck[k] + u_bck[j] * d_u[k] ) * ( gmunu[i][r] + u_bck[i] * u_bck[r] ) )
                                - 1./3 * ( ( d_u[k] * u_bck[r] + u_bck[k] * d_u[r] ) * ( gmunu[i][j] + u_bck[i] * u_bck[j] ) + ( d_u[i] * u_bck[j] + u_bck[i] * d_u[j] ) * ( gmunu[k][r] + u_bck[k] * u_bck[r] ) )
                                ) * u_bck[l] * dpiH_bck[k][r][l] * gmumu[k] * gmumu[r] / gamma * dt );
               }
           }
          }
         }
        c->addPi0(-delPiPi * ( c->getPiH0() * du_bck + c->getPiH0_bck() * d_du ) / gamma * dt);
    }  // end non-empty cell
   }   // end loop #1
    
 // 3) -- advection ---
 // takes into account only the background velocity - from the look of the term with delta pi^munu time derivative
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenterQbck(eos, tauMinusHalf, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck);
    c->getPrimVarHCenter(eos, tauMinusHalf, d_e, d_p, d_nb, d_nq, d_ns, d_vx, d_vy,
                         d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck);  // getPrimVar() before
       
    if (e_bck + d_e <= 0. || e_bck <= 0.) continue;// modified condition
    double xm = -vx_bck * dt / f->getDx();
    double ym = -vy_bck * dt / f->getDy();
    double zm = -vz_bck * dt / f->getDz() / tauMinusHalf;
    double xmH = -vx_bck * dt / f->getDx() / 2.0;
    double ymH = -vy_bck * dt / f->getDy() / 2.0;
    double zmH = -vz_bck * dt / f->getDz() / tauMinusHalf / 2.0;
    double wx[2] = {(1. - fabs(xm)), fabs(xm)};
    double wy[2] = {(1. - fabs(ym)), fabs(ym)};
    double wz[2] = {(1. - fabs(zm)), fabs(zm)};
    double wxH[2] = {(1. - fabs(xmH)), fabs(xmH)};
    double wyH[2] = {(1. - fabs(ymH)), fabs(ymH)};
    double wzH[2] = {(1. - fabs(zmH)), fabs(zmH)};
       
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) {
      d_pi[i][j] = d_piH[i][j] = 0.0;
      pi_bck[i][j] = piH_bck[i][j] = 0.0;
     }
    d_Pi = d_PiH = Pi_bck = PiH_bck = 0.0;
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       // pi,Pi-->full step, piH,PiH-->half-step
       Cell *c1 = f->getCell(ix + jx * sign(xm), iy + jy * sign(ym),
                             iz + jz * sign(zm));
       for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
         d_pi[i][j] += wx[jx] * wy[jy] * wz[jz] * c1->getpi0(i, j);
         d_piH[i][j] += wxH[jx] * wyH[jy] * wzH[jz] * c1->getpiH0(i, j);
        }
        d_Pi += wx[jx] * wy[jy] * wz[jz] * c1->getPi0();
        d_PiH += wxH[jx] * wyH[jy] * wzH[jz] * c1->getPiH0();
      }
    //--------- debug part: trace check, diag check, transversality check
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) {
         //cout << pi[i][j] << endl;
      if (d_pi[i][j] != 0. && fabs(1.0 - d_pi[j][i] / d_pi[i][j]) > 1e-10)
       cout << "non-diag: " << d_pi[i][j] << "  " << d_pi[j][i] << endl;
     }
    //------ end debug
    //======= hydro applicability check (viscous corrections limiter):
    double maxT0 = max((e_bck + p_bck) / (1. - vx_bck * vx_bck - vy_bck * vy_bck - vz_bck * vz_bck) - p_bck,
                       (e_bck + p_bck) * (vx_bck * vx_bck + vy_bck * vy_bck + vz_bck * vz_bck) /
                               (1. - vx_bck * vx_bck - vy_bck * vy_bck - vz_bck * vz_bck) +
                           p_bck);
//     double maxpi = max(fabs(pi[1][1]),fabs(pi[2][2])) ;
    double maxpi = 0.;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++){
            if (fabs(d_pi[i][j]) > maxpi) maxpi = fabs(d_pi[i][j]); // modified condition
        }
    bool rescaled = false;
    if (maxT0 / maxpi < 1.0) { // I am not sure how this rescaling should work - if the fluctuation is comparable to background, then it should be rescaled?
//     for (int i = 0; i < 4; i++)
//      for (int j = 0; j < 4; j++) {
//       d_pi[i][j] = 0.1 * d_pi[i][j] * maxT0 / maxpi;
//       d_piH[i][j] = 0.1 * d_piH[i][j] * maxT0 / maxpi;
//      }
//     rescaled = true;
    }
    if (fabs(d_Pi) > d_p) {
     if (d_Pi != 0.) d_Pi = 1 * d_Pi / fabs(d_Pi) * d_p;     //modified .1 -> 1
     if (d_PiH != 0.) d_PiH = 1 * d_PiH / fabs(d_PiH) * d_p; //modified .1 -> 1
     rescaled = true;
    }
    if (rescaled)
     c->setViscCorrCutFlag(maxT0 / maxpi);
    else
     c->setViscCorrCutFlag(1.);
    // updating to the new values
    for (int i = 0; i < 4; i++)
     for (int j = 0; j <= i; j++) {
             c->setpi(i, j, d_pi[i][j]);
             c->setpiH(i, j, d_piH[i][j]);
     }
    c->setPi(d_Pi);
    c->setPiH(d_PiH);
    // source term  - (tau+dt)*delta_Q_(i+1)/delta_tau
//       cout << "Trpi:   " << d_pi[0][0] + d_pi[1][1] + d_pi[2][2] + d_pi[3][3] << endl;
       
    double gamma = 1.0 / sqrt(1.0 - vx_bck * vx_bck - vy_bck * vy_bck - vz_bck * vz_bck);
    double u_bck[4]; // background
    double d_u[4];// fluctuations
    u_bck[0] = gamma;
    u_bck[1] = u_bck[0] * vx_bck;
    u_bck[2] = u_bck[0] * vy_bck;
    u_bck[3] = u_bck[0] * vz_bck;
       
    d_u[0] = 0.;
    d_u[1] = gamma * d_vx;
    d_u[2] = gamma * d_vy;
    d_u[3] = gamma * d_vz;
            
    double d_flux[4];
       for (int i = 0; i < 4; i++){
           d_flux[i] = -tau * (c->getpi(0, i) + c->getPi() * ( u_bck[0] * d_u[i] + d_u[0] * u_bck[i] ) );
       }
    d_flux[0] += tau * c->getPi();
    c->addFlux(d_flux[0], d_flux[1], d_flux[2], d_flux[3], 0., 0., 0.);
   }  // advection loop (all cells)
}

// this procedure explicitly uses T_==0, X_==1, Y_==2, Z_==3
void Hydro::visc_flux(Cell *left, Cell *right, int direction, int ix, int iy, int iz) {
 double d_flux[4];
 int ind2 = 0;
 double dxa = 0.;
 #ifdef CARTESIAN
 double tauMinusHalf = 1.0;
 double tauPlusHalf = 1.0;
 #else
 double tauMinusHalf = tau - 0.5 * dt;
 double tauPlusHalf = tau + 0.5 * dt;
 #endif
 // exit if noth cells are not full with matter
 if (left->getM(direction) < 1. && right->getM(direction) < 1.) return;

 if (direction == X_)
  dxa = f->getDx();
 else if (direction == Y_)
  dxa = f->getDy();
 else if (direction == Z_)
  dxa = f->getDz() * tauPlusHalf;
 double d_e, d_p, d_nb, d_nq, d_ns, d_vxl, d_vyl, d_vzl, d_vxr, d_vyr, d_vzr;
 double  el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck;
 double  er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck;
 // we need to know the velocities at both cell centers at (n+1/2) in order to
 // interpolate to
 // get the value at the interface
 left->getPrimVarHCenterQbck(eos, tauMinusHalf, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck);
 left->getPrimVarHCenter(eos, tauMinusHalf, d_e, d_p, d_nb, d_nq, d_ns, d_vxl, d_vyl, d_vzl, el_bck, pl_bck, nbl_bck, nql_bck, nsl_bck, vxl_bck, vyl_bck, vzl_bck);
 right->getPrimVarHCenterQbck(eos, tauMinusHalf, er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck);
 right->getPrimVarHCenter(eos, tauMinusHalf, d_e, d_p, d_nb, d_nq, d_ns, d_vxr, d_vyr, d_vzr, er_bck, pr_bck, nbr_bck, nqr_bck, nsr_bck, vxr_bck, vyr_bck, vzr_bck);
 d_vxl = 0.5 * (d_vxl + d_vxr);
 d_vyl = 0.5 * (d_vyl + d_vyr);
 d_vzl = 0.5 * (d_vzl + d_vzr);
 
 vxl_bck = 0.5 * (vxl_bck + vxr_bck);
 vyl_bck = 0.5 * (vyl_bck + vyr_bck);
 vzl_bck = 0.5 * (vzl_bck + vzr_bck);
    
 double v_bck = sqrt(vxl_bck * vxl_bck + vyl_bck * vyl_bck + vzl_bck * vzl_bck);//background norm
 double d_v = sqrt(vxl_bck * vxl_bck + vyl_bck * vyl_bck + vzl_bck * vzl_bck + 2 * vxl_bck * d_vxl + 2 * vxl_bck * d_vxl + 2 * vxl_bck * d_vxl);//linearized full norm
// if (v > 1.) {
//     vxl = 0.99 * vxl / v; How to define this velocity rescaling???
//     vyl = 0.99 * vyl / v;
//     vzl = 0.99 * vzl / v;
// }
    
 double gamma = 1. / sqrt(1. - v_bck * v_bck);
 double uuu_bck[4] = {gamma, gamma * vxl_bck, gamma * vyl_bck, gamma * vzl_bck}; // background
 double d_uuu[4] = {0., gamma * d_vxl, gamma * d_vyl, gamma * d_vzl}; // fluctuation
    
 double gmumu[4] = {1., -1., -1., -1.};
 if (direction == X_)
  ind2 = 1;
 else if (direction == Y_)
  ind2 = 2;
 else if (direction == Z_)
  ind2 = 3;
 for (int ind1 = 0; ind1 < 4; ind1++) {
  d_flux[ind1] = 0.5 * (left->getpiH(ind1, ind2) + right->getpiH(ind1, ind2));//tady zmena?
 if (ind1 == ind2){
     d_flux[ind1] += -0.5 * (left->getPiH() + right->getPiH()) *
     gmumu[ind1];  // gmunu is diagonal
 }
     d_flux[ind1] += 0.5 * (left->getPiH() + right->getPiH()) * uuu_bck[ind1] * uuu_bck[ind2]
             + 0.5 * (left->getPiH_bck() + right->getPiH_bck()) * ( uuu_bck[ind1] * d_uuu[ind2] + d_uuu[ind1] * uuu_bck[ind2] );
 }
 for (int i = 0; i < 4; i++) d_flux[i] = d_flux[i] * tauMinusHalf * dt / dxa;
            double Tl[7];
            double Tr[7];
            left->getQ(Tl);
            right->getQ(Tr);
            if ((abs(Tl[0]/6 - d_flux[0]) > el_bck/6) || (abs(Tr[0]/6 + d_flux[0]) > er_bck/6)) {
                double max = max( abs( (Tl[0]/6 - d_flux[0])/el_bck/6 ), abs( (Tr[0]/6 + d_flux[0])/er_bck/6 ) );
                d_flux[0] = d_flux[0] / max * 0.9;
                N_visc++;
                //                cout << "v" << endl;
            }
            for (int i=1; i<4; i++) {
                if ((abs(Tl[i] - d_flux[i]) > abs(Tl[0] - d_flux[0])/2./sqrt(3.)) || (abs(Tr[i] + d_flux[i]) > abs(Tr[0] + d_flux[0])/2./sqrt(3.))) {
                    double maximum = std::max( abs( (Tl[i] - d_flux[i])/(Tl[0] - d_flux[0]) ), abs( (Tr[i] + d_flux[i])/(Tr[0] + d_flux[0]) ) );
                    d_flux[i] = d_flux[i] / maximum * 0.9;
                }
            }
//             else{
                 left->addFlux(-d_flux[T_], -d_flux[X_], -d_flux[Y_], -d_flux[Z_], 0., 0., 0.);
                 right->addFlux(d_flux[T_], d_flux[X_], d_flux[Y_], d_flux[Z_], 0., 0., 0.);
//             }
}

void Hydro::performStep(double ctime) {
 // debugRiemann = false ; // turn off debug output
// ofstream myfile;
 f->updateM(tau, dt);

 tau_z = dt / 2. / log(1 + dt / 2. / tau);
    int nx = f->getNX();
    int ny = f->getNY();
    int nz = f->getNZ();
    int total = nx * ny * nz;
    int dims[3] = {nx, ny, nz};
    double T00=0.;
    double T_mean[7]={0.,0.,0.,0.,0.,0.,0.};
    double variance_e = 0.;
    
//  Some stuff to print out the values
    std::vector<double> values; // vector for FFT
    for (int ix = 0; ix < f->getNX(); ix++) {
        for (int iy = 0; iy < f->getNY(); iy++){
            for (int iz = 0; iz < f->getNZ(); iz++){
//                double T00=0.;
                double e, p, nb, nq, ns, vx, vy, vz, cs, T, mub, muq, mus;
                double e__0, p__0, nb__0, nq__0, ns__0, vx__0, vy__0, vz__0, cs__0, T__0, mub__0, muq__0, mus__0;
                Cell *c = f->getCell(ix, iy, iz);
                c -> getPrimVarQbck(eos, tau, e__0, p__0, nb__0, nq__0, ns__0, vx__0, vy__0, vz__0);
                c -> getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz, e__0, p__0, nb__0, nq__0, ns__0, vx__0, vy__0, vz__0);
//                cs = eos->cs();
                eos->eos(e__0, nb__0, nq__0, ns__0, T__0, mub, muq, mus, p__0);
                double Tmunu[7];
                double T_bck[7];
                c->getQ(Tmunu);
                c->getQbck(T_bck);
                double T00id = Tmunu[0];
                double T00visc = c->getpi0(0, 0);// nebo H?
//                cout << T_bck[0] << endl;
//                cout << ctime << "  " << T00id << "     " << c->getpiH(0, 0) << "     " << c->getpi(0, 0) << "     " << xi00_cell[ix][iy][iz] << "    " << ix << "    " << iy << "    " << iz << endl;
                T00 += T00id + T00visc;
                for (int i=0; i<7; i++) {
                    T_mean[i] += Tmunu[i]/total;
                }
//                double etaS, zetaS;
//                trcoeff->getEta(e__0, nb__0, T__0, etaS, zetaS);
//                values.push_back(T00);//toto
                variance_e = e*e/total - e*e/total/total;
                values.push_back(e);//toto
            }
        }
    }
    ofstream myfile3;
    myfile3.open ("variance_e.dat", ios::app);
    myfile3 << ctime << "      " << variance_e << endl;
    myfile3.close();
    
    ofstream myfile2;
    myfile2.open ("ENERGY_CONSERVATION_60_100.dat", ios::app);
    myfile2 << T_mean[0] << "      " << T_mean[1] << "     " << T_mean[2] << "     " << T_mean[3] << endl;
    myfile2.close();

//##################### FFT ################################
    if (ctime>59.98) {
        ofstream myfile;
        kiss_fftnd_cfg cfg = kiss_fftnd_alloc(dims, 3, false, NULL, NULL);
        kiss_fft_cpx *in = new kiss_fft_cpx[total];
        kiss_fft_cpx *out = new kiss_fft_cpx[total];
        for (int i=0; i<total; i++) {
            in[i].r = values.at(i);
            in[i].i = 0.;
//            out[i].r = 0.;
//            out[i].i = 0.;
        }
        kiss_fftnd(cfg, in, out);
        free(cfg);
        
        double xi_FT[nx][ny][nz];
        double phase[nx][ny][nz];
                double S_K[nx];
        for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
                for (int iz=0; iz<nz; iz++) {
                    xi_FT[ix][iy][iz] = 0.;
//                    phase[ix][iy][iz] = 0.;
                }
            }
        }
        myfile.open ("FT_e_60_100.dat", ios::app);
        int i = 0;
        for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
                for (int iz=0; iz<nz; iz++) {
                    xi_FT[ix][iy][iz] = (out[i].r*out[i].r + out[i].i*out[i].i)/total;
                    i++;
                }
            }
        }
        delete[] in;//
        delete[] out;//

        for (int ix=0; ix<nx/2+1; ix++) {
            for (int iy=0; iy<ny/2+1; iy++) {
                for (int iz=0; iz<nz/2+1; iz++) {
                    double absK = sqrt(ix*ix+iy*iy+iz*iz);
                    myfile << ctime << "    " << xi_FT[ix][iy][iz] << "   " << ix << "   " << iy << "    " << iz << "    " << absK << endl;
                }
            }
        }
        myfile.close();
    }
//####################### FFT end ##################################
    
 //-----PREDICTOR-ideal
 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->saveQprev();
    c->clearFlux();
   }
 // X dir
 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix + 1, iy, iz), X_, PREDICT, ix, iy, iz);
   }
 //    cout << "predictor X done\n" ;
 // Y dir
 for (int iz = 0; iz < f->getNZ(); iz++)
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, PREDICT, ix, iy, iz);
   }
 //    cout << "predictor Y done\n" ;
 // Z dir
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz + 1), Z_, PREDICT, ix, iy, iz);
   }
 //    cout << "predictor Z done\n" ;

 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
    Cell *c = f->getCell(ix, iy, iz);
    source_step(ix, iy, iz, PREDICT);
    c->updateQtoQhByFlux();
    c->clearFlux();
   }

 //----CORRECTOR-ideal

 tau_z = dt / log(1 + dt / tau);
 // X dir
 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix + 1, iy, iz), X_, CORRECT, ix, iy, iz);
   }
 //    cout << "corrector X done\n" ;
 // Y dir
 for (int iz = 0; iz < f->getNZ(); iz++)
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, CORRECT, ix, iy, iz);
   }
 //    cout << "corrector Y done\n" ;
 // Z dir
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz + 1), Z_, CORRECT, ix, iy, iz);
   }
 //    cout << "corrector Z done\n" ;

 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
    Cell *c = f->getCell(ix, iy, iz);
    source_step(ix, iy, iz, CORRECT);
    c->updateByFlux();
    c->clearFlux();
   }
 #ifdef CARTESIAN
 t += dt;
 #else
 tau += dt;
 #endif
 //f->correctImagCells();  // disabled in box mode

 //===== viscous part ======
 if (trcoeff->isViscous()) {
  ISformal();  // evolution of viscous quantities according to IS equations
  // X dir
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++)
    for (int ix = 0; ix < f->getNX(); ix++) {
//        cout << ix << "     " << iy << "    " << iz << endl;
     visc_flux(f->getCell(ix, iy, iz), f->getCell(ix + 1, iy, iz), X_, ix, iy, iz);
    }
  //    cout << "visc_flux X done\n" ;
  // Y dir
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++) {
//        cout << ix << "     " << iy << "    " << iz << endl;
     visc_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, ix, iy, iz);
    }
  //    cout << "visc_flux Y done\n" ;
  // Z dir
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++)
    for (int iz = 0; iz < f->getNZ(); iz++) {
//        cout << ix << "     " << iy << "    " << iz << endl;
     visc_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz + 1), Z_, ix, iy, iz);
    }

  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++)
    for (int ix = 0; ix < f->getNX(); ix++) {
        Cell *c = f->getCell(ix, iy, iz);
     visc_source_step(ix, iy, iz);
     f->getCell(ix, iy, iz)->updateByViscFlux();
     f->getCell(ix, iy, iz)->clearFlux();
    }
     cout << N_id << "      " << N_visc << endl;
     N_id = 0;
     N_visc = 0;
 } else {  // end viscous part
 }
 //==== finishing work ====
 //f->correctImagCellsFull();  // disabled in box mode
}
