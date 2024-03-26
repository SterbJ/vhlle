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
 
 NEEDS SOME LINEARIZATION IN NB, NS, NQ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <unistd.h>
#include "hdo_box.h"
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "eos.h"
#include "cll.h"
#include "trancoeff.h"

#include "icBox.h"

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

void Hydro::hlle_flux(Cell *left, Cell *right, int direction, int mode, int ix) {
 // for all variables, suffix "l" = left state, "r" = right state
 // with respect to the cell boundary
 double el, er, pl, pr, nbl, nql, nsl, nbr, nqr, nsr, vxl, vxr, vyl, vyr, vzl,
     vzr, bl = 0., br = 0., csb, vb, vb_const, vb_delta, El, Er, dx = 0.;
 double el_0, er_0, pl_0, pr_0, nbl_0, nql_0, nsl_0, nbr_0, nqr_0, nsr_0, vxl_0, vxr_0, vyl_0, vyr_0, vzl_0, vzr_0;
 double Ftl = 0., Fxl = 0., Fyl = 0., Fzl = 0., Fbl = 0., Fql = 0., Fsl = 0.,
        Ftr = 0., Fxr = 0., Fyr = 0., Fzr = 0., Fbr = 0., Fqr = 0., Fsr = 0.;
 double U1l, U2l, U3l, U4l, Ubl, Uql, Usl, U1r, U2r, U3r, U4r, Ubr, Uqr, Usr;
 double v0_L, v0_R, delta_vL, delta_vR, cs;
 double br_const, bl_const, P, Q, C, B;
 double flux[7];
 const double dta = mode == 0 ? dt / 2. : dt;
 double tauFactor = 1.0;  // fluxes are also multiplied by tau, 1 for Cartesian
    

 if (mode == PREDICT) {
  // get primitive quantities from Q_{i+} at previous timestep
  left->getPrimVarRightQ0(eos, tau, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0,
                           direction);
  left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                        direction, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0);
  // ... and Q_{(i+1)-}
  right->getPrimVarLeftQ0(eos, tau, er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0,
                           direction);
  right->getPrimVarLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                        direction, er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0);
     
//     cout << pr_0 << endl;


     El = (el_0 + pl_0) / (1 - vxl_0 * vxl_0 - vyl_0 * vyl_0 - vzl_0 * vzl_0);
     Er = (er_0 + pr_0) / (1 - vxr_0 * vxr_0 - vyr_0 * vyr_0 - vzr_0 * vzr_0);
  #ifndef CARTESIAN
  tauFactor = tau + 0.25 * dt;
  #endif
     //cout << "nbl=    " << nbl << "   nbr=      " << nbr << endl;

 } else {
  // use half-step updated Q's for corrector step
  left->getPrimVarHRightQ0(eos, tau, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0,
                            direction);
  left->getPrimVarHRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                         direction, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0);
  right->getPrimVarHLeftQ0(eos, tau, er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0,
                            direction);
  right->getPrimVarHLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                         direction, er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0);

     double deltaVxl_ = vxl;
     double deltaVyl_ = vyl;
     double deltaVzl_ = vzl;
     double deltaVxr_ = vxr;
     double deltaVyr_ = vyr;
     double deltaVzr_ = vzr;

  El = (el_0 + pl_0) / (1 - vxl_0 * vxl_0 - vyl_0 * vyl_0 - vzl_0 * vzl_0);
  Er = (er_0 + pr_0) / (1 - vxr_0 * vxr_0 - vyr_0 * vyr_0 - vzr_0 * vzr_0);
  #ifndef CARTESIAN
  tauFactor = tau + 0.5 * dt;
  #endif
 }

 if (el + el_0 < 0. || el_0 < 0.) {//??????????????????????????????????????????????????
  el = 0.;
  el_0 = 0.;
  pl = 0.;
  pl_0 = 0.;
 }
 if (er + er_0 < 0. || er_0 < 0.) {
  er = 0.;
  er_0 = 0.;
  pr = 0.;
  pr_0 = 0.;
 }

 if (el+el_0 > 1e10) {
  cout << "e>1e10; debug info below:\n";
  left->Dump(tau);
  // debugRiemann = true ;
     if (mode == PREDICT){
         left->getPrimVarRightQ0(eos, tau, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0,
                               direction);
         left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                               direction, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0);
     }
     else{
         left->getPrimVarHRightQ0(eos, tau, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0,
                                direction);
         left->getPrimVarHRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                                direction, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0);
     }
  // debugRiemann = false ;
  exit(0);
 }

 // skip the procedure for two empty cells
 if (el+el_0 == 0. && er+er_0 == 0.) return;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!________________________________________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (pr + pr_0 < 0.) {
  cout << "Negative pressure" << endl;
//     cout << "      "<< ix << "     " << pr << "    " << pr_0 << "      " << pr+pr_0 << endl;
  left->getPrimVarRightQ0(eos, tau, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0,
                           direction);
  left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                        direction, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0);
  right->getPrimVarLeftQ0(eos, tau, er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0,
                           direction);
  right->getPrimVarLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                        direction, er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0);
 }

 // skip the procedure for two partially vacuum cells
 if (left->getM(direction) < 1. && right->getM(direction) < 1.) return;


    // idea is that e.g. el = el0 + deltaEl
//    double el0 = 0, er0 = 0, pl0 = 0, pr0 = 0, vxl0 = 0, vxr0 = 0, vyl0 = 0, vyr0 = 0, vzl0 = 0, vzr0 = 0;//hodnoty spocitane jako prumer vsech bunek pri kazdem casovem kroku
//fluktuace - promenne - v podstate by to melo byt e.g. deltaEl = (el-el0)

    double el0, pl0, nbl0, nql0, nsl0, vxl0, vyl0, vzl0;
    double er0, pr0, nbr0, nqr0, nsr0, vxr0, vyr0, vzr0;

    left->getPrimVarRightQ0(eos, tau, el0, pl0, nbl0, nql0, nsl0, vxl0, vyl0, vzl0,
                          direction);
    right->getPrimVarLeftQ0(eos, tau, er0, pr0, nbr0, nqr0, nsr0, vxr0, vyr0, vzr0,
                          direction);
    
    double gammal_const = 1./ sqrt(1. - vxl0 * vxl0 - vyl0 * vyl0 - vzl0 * vzl0);
    double gammar_const = 1./ sqrt(1. - vxr0 * vxr0 - vyr0 * vyr0 - vzr0 * vzr0);
    double gammal_delta = 0;//vxl0 * vxl + vyl0 * vyl + vzl0 * vzl;
    double gammar_delta = 0;//vxr0 * vxr + vyr0 * vyr + vzr0 * vzr;
    
    U1l = (el0 + pl0) * (gammal_const * gammal_const * vxl ) + (el + pl) * gammal_const * gammal_const * vxl0;
    U2l = (el0 + pl0) * (gammal_const * gammal_const * vyl ) + (el + pl) * gammal_const * gammal_const * vyl0;
    U3l = (el0 + pl0) * (gammal_const * gammal_const * vzl ) + (el + pl) * gammal_const * gammal_const * vzl0;
    U4l = (el + pl) * gammal_const * gammal_const - pl;

    U1r = (er0 + pr0) * (gammar_const * gammar_const * vxr ) + (er + pr) * gammar_const * gammar_const * vxr0;
    U2r = (er0 + pr0) * (gammar_const * gammar_const * vyr ) + (er + pr) * gammar_const * gammar_const * vyr0;
    U3r = (er0 + pr0) * (gammar_const * gammar_const * vzr ) + (er + pr) * gammar_const * gammar_const * vzr0;
    U4r = (er + pr) * gammar_const * gammar_const - pr;


 Ubl = gammal_const * nbl;
 Uql = gammal_const * nql;
 Usl = gammal_const * nsl;

 Ubr = gammar_const * nbr;
 Uqr = gammar_const * nqr;
 Usr = gammar_const * nsr;

 if (direction == X_) {
     
     Ftl = (el0 + pl0) * (gammal_const * gammal_const * vxl ) + (el + pl) * gammal_const * gammal_const * vxl0;
     Fxl = (el0 + pl0) * (gammal_const * vxl0 * gammal_const * vxl + gammal_const * vxl0 * gammal_const * vxl) + (el + pl) * gammal_const * vxl0 * gammal_const * vxl0 + pl;
     Fyl = (el0 + pl0) * (gammal_const * vyl0 * gammal_const * vxl + gammal_const * vxl0 * gammal_const * vyl) + (el + pl) * gammal_const * vyl0 * gammal_const * vxl0;
     Fzl = (el0 + pl0) * (gammal_const * vzl0 * gammal_const * vxl + gammal_const * vxl0 * gammal_const * vzl) + (el + pl) * gammal_const * vzl0 * gammal_const * vxl0;
     
     Ftr = (er0 + pr0) * (gammar_const * gammar_const * vxr ) + (er + pr) * gammar_const * gammar_const * vxr0;
     Fxr = (er0 + pr0) * (gammar_const * vxr0 * gammar_const * vxr + gammar_const * vxr0 * gammar_const * vxr) + (er + pr) * gammar_const * vxr0 * gammar_const * vxr0 + pr;
     Fyr = (er0 + pr0) * (gammar_const * vyr0 * gammar_const * vxr + gammar_const * vxr0 * gammar_const * vyr) + (er + pr) * gammar_const * vyr0 * gammar_const * vxr0;
     Fzr = (er0 + pr0) * (gammar_const * vzr0 * gammar_const * vxr + gammar_const * vxr0 * gammar_const * vzr) + (er + pr) * gammar_const * vzr0 * gammar_const * vxr0;
     
  Fbl = Ubl * vxl;
  Fql = Uql * vxl;
  Fsl = Usl * vxl;

  Fbr = Ubr * vxr;
  Fqr = Uqr * vxr;
  Fsr = Usr * vxr;

  // for the case of constant c_s only
  csb = sqrt(eos->cs2() +
             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vxl0 - vxr0, 2));
  vb = (sqrt(El) * vxl0 + sqrt(Er) * vxr0) / (sqrt(El) + sqrt(Er));
  bl = min(0., min((vb - csb) / (1 - vb * csb),
                      (vxl0 - eos->cs()) / (1 - vxl0 * eos->cs())));
  br = max(0., max((vb + csb) / (1 + vb * csb),
                      (vxr0 + eos->cs()) / (1 + vxr0 * eos->cs())));


  dx = f->getDx();

  // bl or br in the case of boundary with vacuum
  if (el + el_0 == 0.) bl = -1.;
  if (er + er_0 == 0.) br = 1.;
 }
 if (direction == Y_) {
     
     Ftl = (el0 + pl0) * (gammal_const * gammal_const * vyl ) + (el + pl) * gammal_const * gammal_const * vyl0;
     Fxl = (el0 + pl0) * (gammal_const * vxl0 * gammal_const * vyl + gammal_const * vyl0 * gammal_const * vxl) + (el + pl) * gammal_const * vxl0 * gammal_const * vyl0;
     Fyl = (el0 + pl0) * (gammal_const * vyl0 * gammal_const * vyl + gammal_const * vyl0 * gammal_const * vyl) + (el + pl) * gammal_const * vyl0 * gammal_const * vyl0 + pl;
     Fzl = (el0 + pl0) * (gammal_const * vzl0 * gammal_const * vyl + gammal_const * vyl0 * gammal_const * vzl) + (el + pl) * gammal_const * vzl0 * gammal_const * vyl0;
     
     Ftr = (er0 + pr0) * (gammar_const * gammar_const * vyr ) + (er + pr) * gammar_const * gammar_const * vyr0;
     Fxr = (er0 + pr0) * (gammar_const * vxr0 * gammar_const * vyr + gammar_const * vyr0 * gammar_const * vxr) + (el + pr) * gammar_const * vxr0 * gammar_const * vyr0;
     Fyr = (er0 + pr0) * (gammar_const * vyr0 * gammar_const * vyr + gammar_const * vyr0 * gammar_const * vyr) + (el + pr) * gammar_const * vyr0 * gammar_const * vyr0 + pr;
     Fzr = (er0 + pr0) * (gammar_const * vzr0 * gammar_const * vyr + gammar_const * vyr0 * gammar_const * vzr) + (el + pr) * gammar_const * vzr0 * gammar_const * vyr0;
     
  Fbl = Ubl * vyl;
  Fql = Uql * vyl;
  Fsl = Usl * vyl;

  Fbr = Ubr * vyr;
  Fqr = Uqr * vyr;
  Fsr = Usr * vyr;

  // for the case of constant c_s only
  csb = sqrt(eos->cs2() +
             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vyl_0 - vyr_0, 2));
  vb = (sqrt(El) * vyl_0 + sqrt(Er) * vyr_0) / (sqrt(El) + sqrt(Er));
  bl = min(0., min((vb - csb) / (1 - vb * csb),
                   (vyl_0 - eos->cs()) / (1 - vyl_0 * eos->cs())));
  br = max(0., max((vb + csb) / (1 + vb * csb),
                   (vyr_0 + eos->cs()) / (1 + vyr_0 * eos->cs())));


  dx = f->getDy();

  // bl or br in the case of boundary with vacuum
  if (el+el_0 == 0.) bl = -1.;
  if (er+er_0 == 0.) br = 1.;
 }
 if (direction == Z_) {
  double tau1 = tauFactor;
     
     Ftl = (el0 + pl0) * (gammal_const * gammal_const * vzl ) / tau1 + (el + pl) * gammal_const * gammal_const * vzl0 / tau1;
     Fxl = (el0 + pl0) * (gammal_const * vxl0 * gammal_const * vzl + gammal_const * vzl0 * gammal_const * vxl) / tau1 + (el + pl) * gammal_const * vxl0 * gammal_const * vzl0 / tau1;
     Fyl = (el0 + pl0) * (gammal_const * vyl0 * gammal_const * vzl + gammal_const * vzl0 * gammal_const * vyl) / tau1 + (el + pl) * gammal_const * vyl0 * gammal_const * vzl0 / tau1;
     Fzl = (el0 + pl0) * (gammal_const * vzl0 * gammal_const * vzl + gammal_const * vzl0 * gammal_const * vzl) / tau1 + (el + pl) * gammal_const * vzl0 * gammal_const * vzl0 / tau1 + pl / tau1;
     
     Ftr = (er0 + pr0) * (gammar_const * gammar_const * vzr ) / tau1 + (er + pr) * gammar_const * gammar_const * vzr0 / tau1;
     Fxr = (er0 + pr0) * (gammar_const * vxr0 * gammar_const * vzr + gammar_const * vzr0 * gammar_const * vxr) / tau1 + (er + pr) * gammar_const * vxr0 * gammar_const * vzr0 / tau1;
     Fyr = (er0 + pr0) * (gammar_const * vyr0 * gammar_const * vzr + gammar_const * vzr0 * gammar_const * vyr) / tau1 + (er + pr) * gammar_const * vyr0 * gammar_const * vzr0 / tau1;
     Fzr = (er0 + pr0) * (gammar_const * vzr0 * gammar_const * vzr + gammar_const * vzr0 * gammar_const * vzr) / tau1 + (er + pr) * gammar_const * vzr0 * gammar_const * vzr0 / tau1 + pr / tau1;
     
  Fbl = Ubl * vzl / tau1;
  Fql = Uql * vzl / tau1;
  Fsl = Usl * vzl / tau1;

  Fbr = Ubr * vzr / tau1;
  Fqr = Uqr * vzr / tau1;
  Fsr = Usr * vzr / tau1;

  // for the case of constant c_s only
  // factor 1/tau accounts for eta-coordinate

  // different estimate
  csb = sqrt(eos->cs2() +
             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vzl_0 - vzr_0, 2));
  vb = (sqrt(El) * vzl_0 + sqrt(Er) * vzr_0) / (sqrt(El) + sqrt(Er));
  bl = 1. / tau * min(0., min((vb - csb) / (1 - vb * csb),
                              (vzl_0 - eos->cs()) / (1 - vzl_0 * eos->cs())));
  br = 1. / tau * max(0., max((vb + csb) / (1 + vb * csb),
                              (vzr_0 + eos->cs()) / (1 + vzr_0 * eos->cs())));

  dx = f->getDz();

  // bl or br in the case of boundary with vacuum
  if (el+el_0 == 0.) bl = -1. / tau;
  if (er+er_0 == 0.) br = 1. / tau;
 }

 if (bl == 0. && br == 0.) return;

//  finally, HLLE formula for the fluxes
 flux[T_] = tauFactor * dta / dx *
            (-bl * br * (U4l - U4r) + br * Ftl - bl * Ftr) / (-bl + br);
 flux[X_] = tauFactor * dta / dx *
            (-bl * br * (U1l - U1r) + br * Fxl - bl * Fxr) / (-bl + br);
 flux[Y_] = tauFactor * dta / dx *
            (-bl * br * (U2l - U2r) + br * Fyl - bl * Fyr) / (-bl + br);
 flux[Z_] = tauFactor * dta / dx *
            (-bl * br * (U3l - U3r) + br * Fzl - bl * Fzr) / (-bl + br);
 flux[NB_] = tauFactor * dta / dx *
             (-bl * br * (Ubl - Ubr) + br * Fbl - bl * Fbr) / (-bl + br);
 flux[NQ_] = tauFactor * dta / dx *
             (-bl * br * (Uql - Uqr) + br * Fql - bl * Fqr) / (-bl + br);
 flux[NS_] = tauFactor * dta / dx *
             (-bl * br * (Usl - Usr) + br * Fsl - bl * Fsr) / (-bl + br);
    
//    if(direction==X_){
//
//        cout << flux[X_] << "       " << flux[Y_] << "      " << flux[Z_] << endl;
//    }
//        if(ix<=40){
//            cout.precision(20);
//        }
//        else{
//            cout.precision(20);
//        }


 if (flux[NB_] != flux[NB_]) {  // if things failed
  cout << "---- error in hlle_flux: f_nb undefined!\n";
  cout << setw(12) << U4l << setw(12) << U1l << setw(12) << U2l << setw(12)
       << U3l << endl;
  cout << setw(12) << U4r << setw(12) << U1r << setw(12) << U2r << setw(12)
       << U3r << endl;
  cout << setw(12) << Ubl << setw(12) << Uql << setw(12) << Usl << endl;
  cout << setw(12) << Ubr << setw(12) << Uqr << setw(12) << Usr << endl;
  cout << setw(12) << Ftl << setw(12) << Fxl << setw(12) << Fyl << setw(12)
       << Fzl << endl;
  cout << setw(12) << Ftr << setw(12) << Fxr << setw(12) << Fyr << setw(12)
       << Fzr << endl;
  exit(1);
 }

 // update the cumulative fluxes in both neighbouring cells
 left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
               -flux[NQ_], -flux[NS_]);
 right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_], flux[NQ_],
                flux[NS_]);
}


void Hydro::source(double tau1, double x, double y, double z, double Q[7],
                   double S[7]) {
 #ifdef CARTESIAN
 // geometrical source term is zero in Cartesian frame
 for (int i = 0; i < 7; i++) S[i] = 0.0;
 #else //linearizovat!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 double _Q[7], _Q0[7], e, p, nb, nq, ns, vx, vy, vz;
    for (int i = 0; i < 7; i++){
        _Q[i] = Q[i] / tau1;  // no tau factor in  _Q
        _Q0[i] = Q0[i] / tau1;
    }
 transformPV(eos, _Q, _Q0, e, p, nb, nq, ns, vx, vy, vz);
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
    double e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0;
 double uuu[4];
    double uuu0[4];
 double k[7];

 #ifdef CARTESIAN
  // there is no geometric viscous source term in Cartesian frame
 #else
 Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenterQ0(eos, tau - dt / 2., e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    c->getPrimVarHCenter(eos, tau - dt / 2., e, p, nb, nq, ns, vx, vy, vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);  // TODO Cartesian
 if (e <= 0.) return;
 uuu0[0] = 1. / sqrt(1. - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0);
 uuu0[1] = uuu[0] * vx_0;
 uuu0[2] = uuu[0] * vy_0;
 uuu0[3] = uuu[0] * vz_0;
    uuu[0] = 0;//vx_0 * vx + vy_0 * vy + vz_0 * vz;
    uuu[1] = uuu0[0] * vx;
    uuu[2] = uuu0[0] * vy;
    uuu[3] = uuu0[0] * vz;

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
void Hydro::NSquant(int ix, int iy, int iz, double pi[4][4], double &Pi, double dmu[4][4], double dmu0[4][4], double &du, double &du0) {
    const double VMIN = 1e-2;
    const double UDIFF = 3.0;
    double e0, e1, p, nb, nq, ns, vx1, vy1, vz1, vx0, vy0, vz0, vxH, vyH, vzH;
    double e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0, e_01, vx_01, vy_01, vz_01, vx_0H, vy_0H, vz_0H;
    double ut0, ux0, uy0, uz0, ut1, ux1, uy1, uz1;//fluctuations
    double ut0_0, ux0_0, uy0_0, uz0_0, ut1_0, ux1_0, uy1_0, uz1_0;//background
    //	double dmu [4][4] ; // \partial_\mu u^\nu matrix
    // coordinates: 0=tau, 1=x, 2=y, 3=eta
    double Z[4][4][4][4];  // Z[mu][nu][lambda][rho]
    double Z0[4][4][4][4]; // Z matrix for background
    double uuu[4];         // the 4-velocity
    double uuu0[4];        // the background 4-velocity
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
                pi[i][j] = 0.;
                dmu[i][j] = 0.;
            }
        Pi = du = 0.;
        return;
    }
    // calculation of \partial_\mu u^\nu matrix
    // mu=first index, nu=second index
    // centered differences with respect to the values at (it+1/2, ix, iy, iz)
    // d_tau u^\mu
    
#ifdef CARTESIAN
    c->getPrimVarPrevQ0(eos, 1.0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    c->getPrimVarPrev(eos, 1.0, e0, p, nb, nq, ns, vx0, vy0, vz0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    c->getPrimVarQ0(eos, 1.0, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);
    c->getPrimVar(eos, 1.0, e1, p, nb, nq, ns, vx1, vy1, vz1, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);
    c->getPrimVarHCenterQ0(eos, 1.0, e_01, p_0, nb_0, nq_0, ns_0, vx_0H, vy_0H, vz_0H);
    c->getPrimVarHCenter(eos, 1.0, e1, p, nb, nq, ns, vxH, vyH, vzH, e_01, p_0, nb_0, nq_0, ns_0, vx_0H, vy_0H, vz_0H);
    double tauPlusHalf = 1.0;
#else
    c->getPrimVarPrevQ0(eos, tau - dt, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    c->getPrimVarPrev(eos, tau - dt, e0, p, nb, nq, ns, vx0, vy0, vz0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    c->getPrimVarQ0(eos, tau, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);
    c->getPrimVar(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);
    c->getPrimVarHCenterQ0(eos, tau - 0.5 * dt, e_01, p_0, nb_0, nq_0, ns_0, vx_0H, vy_0H, vz_0H);
    c->getPrimVarHCenter(eos, tau - 0.5 * dt, e1, p, nb, nq, ns, vxH, vyH, vzH, e_01, p_0, nb_0, nq_0, ns_0, vx_0H, vy_0H, vz_0H);
    double tauPlusHalf = tau + 0.5 * dt;
#endif
    
    //############## get transport coefficients
    double T0, mub, muq, mus;
//    double T1;
    double etaS, zetaS;
    double h = 0.01;
    double delta_s;
    double s0 = eos->s(e_01, nb, nq, ns);  // background entropy density in the current cell
    double dS = (eos->s(e_01 + h, nb, nq, ns) - s0 ) / h;
//    eos->eos(e_01 + h, nb, nq, ns, T1, mub, muq, mus, p_0);
    eos->eos(e_01, nb, nq, ns, T0, mub, muq, mus, p_0);
//    double dT = (T1 - T0) / h;
    delta_s = dS * e1;
    
    trcoeff->getEta(e_01, nb, T0, etaS, zetaS);
    //##############
    // if(e1<0.00004) s=0. ; // negative pressure due to pi^zz for small e
    
    //background
    ut0_0 = 1./ sqrt(1. - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0);
    ux0_0 = ut0_0 * vx_0;
    uy0_0 = ut0_0 * vy_0;
    uz0_0 = ut0_0 * vz_0;
    ut1_0 = 1./ sqrt(1. - vx_01 * vx_01 - vy_01 * vy_01 - vz_01 * vz_01);
    ux1_0 = ut1_0 * vx_01;
    uy1_0 = ut1_0 * vy_01;
    uz1_0 = ut1_0 * vz_01;
    uuu0[0] = 1./ sqrt(1. - vx_0H * vx_0H - vy_0H * vy_0H - vz_0H * vz_0H);
    uuu0[1] = uuu0[0] * vx_0H;
    uuu0[2] = uuu0[0] * vy_0H;
    uuu0[3] = uuu0[0] * vz_0H;
    //fluctuations
    ut0 = 0.;
    ux0 = ut0_0 * vx0;
    uy0 = ut0_0 * vy0;
    uz0 = ut0_0 * vz0;
    ut1 = 0.;
    ux1 = ut1_0 * vx1;
    uy1 = ut1_0 * vy1;
    uz1 = ut1_0 * vz1;
    uuu[0] = 0.;
    uuu[1] = uuu0[0] * vxH;
    uuu[2] = uuu0[0] * vyH;
    uuu[3] = uuu0[0] * vzH;
    
    dmu0[0][0] = (ut1_0 - ut0_0) / dt;
    dmu0[0][1] = (ux1_0 - ux0_0) / dt;
    dmu0[0][2] = (uy1_0 - uy0_0) / dt;
    dmu0[0][3] = (uz1_0 - uz0_0) / dt;
    
    dmu[0][0] = (ut1 - ut0) / dt;
    dmu[0][1] = (ux1 - ux0) / dt;
    dmu[0][2] = (uy1 - uy0) / dt;
    dmu[0][3] = (uz1 - uz0) / dt;
    
    
    // dmu[0][0] = (ut1 * ut1 - ut0 * ut0) / 2. / uuu[0] / dt;
    // dmu[0][1] = (ux1 * ux1 - ux0 * ux0) / 2. / uuu[1] / dt;
    // dmu[0][2] = (uy1 * uy1 - uy0 * uy0) / 2. / uuu[2] / dt;
    // dmu[0][3] = (uz1 * uz1 - uz0 * uz0) / 2. / uuu[3] / dt;
    // if (fabs(0.5 * (ut1 + ut0) / ut1) > UDIFF) dmu[0][0] = (ut1 - ut0) / dt;
    // if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / ux1) > UDIFF)
    //  dmu[0][1] = (ux1 - ux0) / dt;
    // if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uy1) > UDIFF)
    //  dmu[0][2] = (uy1 - uy0) / dt;
    // if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uz1) > UDIFF)
    //  dmu[0][3] = (uz1 - uz0) / dt;
    if (e1+e_01 <= 0. || e0+e_0 <= 0.) {  // matter-vacuum
        dmu[0][0] = dmu[0][1] = dmu[0][2] = dmu[0][3] = 0.;
        dmu0[0][0] = dmu0[0][1] = dmu0[0][2] = dmu0[0][3] = 0.;
    }
    // d_x u^\mu
    f->getCell(ix + 1, iy, iz)->getPrimVarHCenterQ0(eos, tau, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);//0 i 1 maji stejne v_const
    f->getCell(ix + 1, iy, iz)->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);//0 i 1 maji stejne v_const
    f->getCell(ix - 1, iy, iz)->getPrimVarHCenterQ0(eos, tau, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    f->getCell(ix - 1, iy, iz)->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    
    if (e1+e_01 > 0. && e0+e_0 > 0.) {
        
        ut0_0 = 1./ sqrt(1. - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0);
        ux0_0 = ut0_0 * vx_0;
        uy0_0 = ut0_0 * vy_0;
        uz0_0 = ut0_0 * vz_0;
        ut1_0 = 1./ sqrt(1. - vx_01 * vx_01 - vy_01 * vy_01 - vz_01 * vz_01);
        ux1_0 = ut1_0 * vx_01;
        uy1_0 = ut1_0 * vy_01;
        uz1_0 = ut1_0 * vz_01;
        
        ut0 = 0.;
        ux0 = ut0_0 * vx0;
        uy0 = ut0_0 * vy0;
        uz0 = ut0_0 * vz0;
        ut1 = 0.;
        ux1 = ut1_0 * vx1;
        uy1 = ut1_0 * vy1;
        uz1 = ut1_0 * vz1;
        
        dmu0[1][0] = 0.5 * (ut1_0 - ut0_0) / dx;//why 0.5??????????
        dmu0[1][1] = 0.5 * (ux1_0 - ux0_0) / dx;
        dmu0[1][2] = 0.5 * (uy1_0 - uy0_0) / dx;
        dmu0[1][3] = 0.5 * (uz1_0 - uz0_0) / dx;

        
        dmu[1][0] = 0.5 * (ut1 - ut0) / dx;//why 0.5??????????
        dmu[1][1] = 0.5 * (ux1 - ux0) / dx;
        dmu[1][2] = 0.5 * (uy1 - uy0) / dx;
        dmu[1][3] = 0.5 * (uz1 - uz0) / dx;
        
        //  dmu[1][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dx;
        //  dmu[1][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dx;
        //  dmu[1][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dx;
        //  dmu[1][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dx;
        //  if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
        //   dmu[1][0] = 0.5 * (ut1 - ut0) / dx;
        //  if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
        //   dmu[1][1] = 0.5 * (ux1 - ux0) / dx;
        //  if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
        //   dmu[1][2] = 0.5 * (uy1 - uy0) / dx;
        //  if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
        //   dmu[1][3] = 0.5 * (uz1 - uz0) / dx;
         } else {  // matter-vacuum
          dmu[1][0] = dmu[1][1] = dmu[1][2] = dmu[1][3] = 0.;
          dmu0[1][0] = dmu0[1][1] = dmu0[1][2] = dmu0[1][3] = 0.;

        }
        if (fabs(dmu[1][3]) > 1e+10)
            cout << "dmu[1][3]:  " << uz1 << "  " << uz0 << "  " << uuu[3] << endl;
        // d_y u^\mu
        f->getCell(ix, iy + 1, iz)->getPrimVarHCenterQ0(eos, tau, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);
        f->getCell(ix, iy + 1, iz)->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);
        f->getCell(ix, iy - 1, iz)->getPrimVarHCenterQ0(eos, tau, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
        f->getCell(ix, iy - 1, iz)->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    
        if (e1+e_01 > 0. && e0+e_0 > 0.) {
            
            ut0_0 = 1./ sqrt(1. - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0);
            ux0_0 = ut0_0 * vx_0;
            uy0_0 = ut0_0 * vy_0;
            uz0_0 = ut0_0 * vz_0;
            ut1_0 = 1./ sqrt(1. - vx_01 * vx_01 - vy_01 * vy_01 - vz_01 * vz_01);
            ux1_0 = ut1_0 * vx_01;
            uy1_0 = ut1_0 * vy_01;
            uz1_0 = ut1_0 * vz_01;
            
            ut0 = 0.;
            ux0 = ut0_0 * vx0;
            uy0 = ut0_0 * vy0;
            uz0 = ut0_0 * vz0;
            ut1 = 0.;
            ux1 = ut1_0 * vx1;
            uy1 = ut1_0 * vy1;
            uz1 = ut1_0 * vz1;
            
            dmu0[2][0] = 0.5 * (ut1_0 - ut0_0) / dy;
            dmu0[2][1] = 0.5 * (ux1_0 - ux0_0) / dy;
            dmu0[2][2] = 0.5 * (uy1_0 - uy0_0) / dy;
            dmu0[2][3] = 0.5 * (uz1_0 - uz0_0) / dy;
            
            dmu[2][0] = 0.5 * (ut1 - ut0) / dy;
            dmu[2][1] = 0.5 * (ux1 - ux0) / dy;
            dmu[2][2] = 0.5 * (uy1 - uy0) / dy;
            dmu[2][3] = 0.5 * (uz1 - uz0) / dy;
            
            //  dmu[2][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dy;
            //  dmu[2][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dy;
            //  dmu[2][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dy;
            //  dmu[2][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dy;
            //  if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
            //   dmu[2][0] = 0.5 * (ut1 - ut0) / dy;
            //  if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
            //   dmu[2][1] = 0.5 * (ux1 - ux0) / dy;
            //  if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
            //   dmu[2][2] = 0.5 * (uy1 - uy0) / dy;
            //  if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
            //   dmu[2][3] = 0.5 * (uz1 - uz0) / dy;
        } else {  // matter-vacuum
            dmu[2][0] = dmu[2][1] = dmu[2][2] = dmu[2][3] = 0.;
            dmu0[2][0] = dmu0[2][1] = dmu0[2][2] = dmu0[2][3] = 0.;

        }
        // d_z u^\mu
        f->getCell(ix, iy, iz + 1)->getPrimVarHCenterQ0(eos, tau, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);
        f->getCell(ix, iy, iz + 1)->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1, e_01, p_0, nb_0, nq_0, ns_0, vx_01, vy_01, vz_01);
        f->getCell(ix, iy, iz - 1)->getPrimVarHCenterQ0(eos, tau, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
        f->getCell(ix, iy, iz - 1)->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    
        if (e1+e_01 > 0. && e0+e_0 > 0.) {
            
            ut0_0 = 1./ sqrt(1. - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0);
            ux0_0 = ut0_0 * vx_0;
            uy0_0 = ut0_0 * vy_0;
            uz0_0 = ut0_0 * vz_0;
            ut1_0 = 1./ sqrt(1. - vx_01 * vx_01 - vy_01 * vy_01 - vz_01 * vz_01);
            ux1_0 = ut1_0 * vx_01;
            uy1_0 = ut1_0 * vy_01;
            uz1_0 = ut1_0 * vz_01;
            
            ut0 = 0.;
            ux0 = ut0_0 * vx0;
            uy0 = ut0_0 * vy0;
            uz0 = ut0_0 * vz0;
            ut1 = 0.;
            ux1 = ut1_0 * vx1;
            uy1 = ut1_0 * vy1;
            uz1 = ut1_0 * vz1;
            
            dmu0[3][0] = 0.5 * (ut1_0 - ut0_0) / dz / tauPlusHalf;
            dmu0[3][1] = 0.5 * (ux1_0 - ux0_0) / dz / tauPlusHalf;
            dmu0[3][2] = 0.5 * (uy1_0 - uy0_0) / dz / tauPlusHalf;
            dmu0[3][3] = 0.5 * (uz1_0 - uz0_0) / dz / tauPlusHalf;
            
            dmu[3][0] = 0.5 * (ut1 - ut0) / dz / tauPlusHalf;
            dmu[3][1] = 0.5 * (ux1 - ux0) / dz / tauPlusHalf;
            dmu[3][2] = 0.5 * (uy1 - uy0) / dz / tauPlusHalf;
            dmu[3][3] = 0.5 * (uz1 - uz0) / dz / tauPlusHalf;
            
            //  dmu[3][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dz / tauPlusHalf;
            //  dmu[3][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dz / tauPlusHalf;
            //  dmu[3][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dz / tauPlusHalf;
            //  dmu[3][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dz / tauPlusHalf;
            //  if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
            //   dmu[3][0] = 0.5 * (ut1 - ut0) / dz / tauPlusHalf;
            //  if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
            //   dmu[3][1] = 0.5 * (ux1 - ux0) / dz / tauPlusHalf;
            //  if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
            //   dmu[3][2] = 0.5 * (uy1 - uy0) / dz / tauPlusHalf;
            //  if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
            //   dmu[3][3] = 0.5 * (uz1 - uz0) / dz / tauPlusHalf;
        } else {  // matter-vacuum
            dmu[3][0] = dmu[3][1] = dmu[3][2] = dmu[3][3] = 0.;
            dmu0[3][0] = dmu0[3][1] = dmu0[3][2] = dmu0[3][3] = 0.;
        }
        // additional terms from Christoffel symbols :)
#ifndef CARTESIAN
        dmu[3][0] += uuu[3] / (tau - 0.5 * dt);
        dmu[3][3] += uuu[0] / (tau - 0.5 * dt);
    dmu0[3][0] += uuu0[3] / (tau - 0.5 * dt);
    dmu0[3][3] += uuu0[0] / (tau - 0.5 * dt);
#endif

        // calculation of Z0[mu][nu][lambda][rho]
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    for (int l = 0; l < 4; l++) Z0[i][j][k][l] = 0.0;
        // filling Z0 matrix - background
        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
                for (int lam = 0; lam < 4; lam++)
                    for (int rho = 0; rho < 4; rho++) {
                        if (nu == rho)
                            Z0[mu][nu][lam][rho] += 0.5 * (gmunu[mu][lam] - uuu0[mu] * uuu0[lam]);
                        
                        if (mu == rho)
                            Z0[mu][nu][lam][rho] += 0.5 * (gmunu[nu][lam] - uuu0[nu] * uuu0[lam]);
                        
                        if (lam == rho)
                            Z0[mu][nu][lam][rho] -= (gmunu[mu][nu] - uuu0[mu] * uuu0[nu]) / 3.0;
                        
                    }
    
    // calculation of Z[mu][nu][lambda][rho]
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                for (int l = 0; l < 4; l++) Z[i][j][k][l] = 0.0;
    // filling Z matrix - perturbations
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            for (int lam = 0; lam < 4; lam++)
                for (int rho = 0; rho < 4; rho++) {
                    if (nu == rho)
                        Z[mu][nu][lam][rho] += 0.5 * (uuu[mu] * uuu0[lam] + uuu0[mu] * uuu[lam]);
                    if (mu == rho)
                        Z[mu][nu][lam][rho] += 0.5 * (uuu[nu] * uuu0[lam] + uuu0[nu] * uuu[lam]);
                    if (lam == rho)
                        Z[mu][nu][lam][rho] -= (uuu[mu] * uuu0[nu] + uuu0[mu] * uuu[nu]) / 3.0;
                }
    
        // calculating sigma[mu][nu]
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                pi[i][j] = 0.0;
                for (int k = 0; k < 4; k++)
                    for (int l = 0; l < 4; l++) {
                        pi[i][j] += ( Z0[i][j][k][l] * dmu[k][l] - Z[i][j][k][l] * dmu0[k][l] ) * 2.0 * etaS * s0 / 5.068
                                    + Z0[i][j][k][l] * dmu0[k][l] * 2.0 * etaS * delta_s / 5.068 ;
//                        cout << Z0[i][j][k][l] << "     " << dmu[k][l] << "     " << i << "     " << j << "     " << k << "     " << l << endl;
//                        if(i==1&&j==1){
//                            pi[i][j]=0.;
//                        }
                    }
//                if((ix==20 || ix ==19 || ix ==21 || ix ==60 || ix ==61 || ix ==59)&&iy==0&&iz==0&&i==j&&i==1){
//                    cout << ix << "     " << pi[i][j] << "   " << i << "     " << j << "    " << dmu[i][j] << endl;
//                }
            }

        Pi = -zetaS * s0 * (dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3]) / 5.068 - zetaS * delta_s * (dmu0[0][0] + dmu0[1][1] + dmu0[2][2] + dmu0[3][3]) / 5.068;  // 5.068 -> fm^{-4} --> GeV/fm^3
        du = dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3];
        du0 = dmu0[0][0] + dmu0[1][1] + dmu0[2][2] + dmu0[3][3];
        //--------- debug part: NaN/inf check, trace check, diag check, transversality
        // check
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                if (pi[i][j] != 0. && fabs(1.0 - pi[j][i] / pi[i][j]) > 1e-10)
                    cout << "non-diag: " << pi[i][j] << "  " << pi[j][i] << endl;
                if (std::isinf(pi[i][j]) || std::isnan(pi[i][j])) {
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
//    c->getPrimVarQ0(eos, tau, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
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
 double e, p, nb, nq, ns, vx, vy, vz, T0, mub, muq, mus;
 double e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0;
 double piNS[4][4], piNS0[4][4], sigNS[4][4], sigNS0[4][4], PiNS, PiNS0, dmu[4][4], dmu0[4][4], du, du0, pi[4][4], piH[4][4], Pi, PiH, dpi0[4][4][4], dpi0H[4][4][4];
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
                      {0, 0, -1, 0},//????????????????????????????????????????????????????????????,, tau factor????????
                      {0, 0, 0, -1}};
 // loop #1 (relaxation+source terms)
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenterQ0(eos, tauMinusHalf, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    c->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vx, vy,
                         vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);  // instead of getPrimVar()
       
    if (e+e_0 <= 0.) {             // empty cell?
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
       c->setpiH0(i, j, 0.0);
       c->setpi0(i, j, 0.0);
      }
     c->setPiH0(0.0);
     c->setPi0(0.0);
    } else {  // non-empty cell
     // 1) relaxation(pi)+source(pi) terms for half-step
     double gamma = 1.0 / sqrt(1.0 - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0);
     double u0[4];
     u0[0] = gamma;
     u0[1] = u0[0] * vx_0;
     u0[2] = u0[0] * vy_0;
     u0[3] = u0[0] * vz_0;
        
        double u[4];
        u[0] = 0.;
        u[1] = gamma * vx;
        u[2] = gamma * vy;
        u[3] = gamma * vz;
            
        // source term  + tau*delta_Q_i/delta_tau
        double flux[4];
        for (int i = 0; i < 4; i++){
            flux[i] = tauMinusDt * (c->getpi(0, i) + c->getPi() * ( u[0] * u0[i] + u0[0] * u[i] ) ); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        }
        flux[0] += -tauMinusDt * c->getPi();
          // cout << c->getPi() << endl;
        c->addFlux(flux[0], flux[1], flux[2], flux[3], 0., 0., 0.);
        // now calculating viscous terms in NS limit
        NSquant(ix, iy, iz, piNS, PiNS, dmu, dmu0, du, du0);
        PiNS0 = 0.;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j <= i; j++) {
                piNS0[i][j]=0.;
            }
// derivative of background pi^munu
#ifdef CARTESIAN
    double tauPlusHalf_0 = 1.0;
#else
    double tauPlusHalf_0 = tau + 0.5 * dt;
#endif
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            dpi0[i][j][0] = (f->getCell(ix, iy, iz)->getpi_bck(i,j) - f->getCell(ix, iy, iz)->getpi_bck_prev(i,j)) / dt;
            dpi0[i][j][1] = 0.5 * (f->getCell(ix + 1, iy, iz)->getpi_bck(i,j) - f->getCell(ix - 1, iy, iz)->getpi_bck(i,j)) / dx;
            dpi0[i][j][2] = 0.5 * (f->getCell(ix, iy + 1, iz)->getpi_bck(i,j) - f->getCell(ix, iy - 1, iz)->getpi_bck(i,j)) / dy;
            dpi0[i][j][3] = 0.5 * (f->getCell(ix, iy, iz + 1)->getpi_bck(i,j) - f->getCell(ix, iy, iz - 1)->getpi_bck(i,j)) / dz / tauPlusHalf_0;
        }
    }
        
        
     double T1;
     double h = 0.01;
     eos->eos(e_0, nb, nq, ns, T0, mub, muq, mus, p_0);
     eos->eos(e+e_0, nb, nq, ns, T1, mub, muq, mus, p_0);
     double dT = (T1 - T0) / h ;
     double etaS, zetaS;
     trcoeff->getEta(e_0, nb, T0, etaS, zetaS);
     const double s0 = eos->s(e_0, nb, nq, ns);
     double dS = (eos->s(e_0 + h, nb, nq, ns) - s0 ) / h;
     double delta_s = dS * e;
     const double eta0 = etaS * s0;
     const double delta_eta = etaS * delta_s;
     double eta = eta0 + delta_eta;
     // auxiliary variable sigmaNS = piNS / (2*eta), 
     // mainly to protect against division by zero in the eta=0 case.
     for(int i=0; i<4; i++)
     for(int j=0; j<4; j++) {
      sigNS0[i][j] = 0.;
      sigNS[i][j] = 0.5 * (piNS[i][j] * 5.068 - 2. * delta_eta * sigNS0[i][j] ) / eta0 ;
//      sigNS0[i][j] = 0.5 * piNS0[i][j] / eta * 5.068;
         if(eta<=0.0){
             sigNS[i][j] = 0.0;
             sigNS0[i][j] = 0.;
         }
     }
     //############# get relaxation times
     double taupi0, tauPi0;
     double delta_taupi, delta_tauPi;
     trcoeff->getTau(e_0, e, nb, T0, dT, taupi0, delta_taupi, tauPi0, delta_tauPi);
//        cout << taupi << endl;
     double deltapipi, taupipi, lambdapiPi, phi7, phi70, delta_phi7, delPiPi, lamPipi;
     trcoeff->getOther(e_0, nb, nq, ns, deltapipi, taupipi, lambdapiPi, phi7);
     phi70 = phi7/taupi0;  // dividing by tau_pi here, to avoid NaNs when tau_pi==0
     delta_phi7 = phi7 / (taupi0 * taupi0) * delta_taupi;
     if(taupi0<0.5*dt)
      deltapipi = taupipi = lambdapiPi = phi7 = 0.0;
     trcoeff->getOtherBulk(e_0, nb, nq, ns, delPiPi, lamPipi);
     if(tauPi0<0.5*dt)
      delPiPi = lamPipi = 0.0;
     //#############
     double Delta[10];
     // relaxation term, piH,PiH-->half-step
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
              Delta[index44(i, j)] = - u0[i] * u0[j];//???????why minus????????
       if (i == j) Delta[index44(i, j)] += gmumu[i];

#ifdef FORMAL_SOLUTION
       c->setpiH0(i, j, (c->getpi(i, j) - piNS[i][j]) *
                                exp(-dt / 2.0 / gamma / taupi0) +// how to do this?
                            piNS[i][j]);
       c->setpiH0(i, j, (c->getpi_bck(i, j) - piNS0[i][j]) * dt / 2.0  / gamma / (taupi0*taupi0) * delta_taupi );
#else
          if(taupi0>0.5*dt){// or waht the condition should be?
              c->setpiH0(i, j, c->getpi(i, j) -
                         (c->getpi(i, j) - piNS[i][j]) * dt / 2.0 / gamma / taupi0);
              c->setpiH0(i, j, (c->getpi_bck(i, j) - piNS0[i][j]) * dt / 2.0 / gamma / (taupi0 * taupi0) * delta_taupi);
          }
      else
       c->setpiH0(i, j, piNS[i][j]);
#endif
      }
#ifdef FORMAL_SOLUTION
     c->setPiH0((c->getPi() - PiNS) * exp(-dt / 2.0 / gamma / tauPi0) + PiNS);
     c->setPiH0( (c->getPi_bck() - PiNS0) * dt / 2.0 / gamma / (tauPi0 * tauPi0) * delta_tauPi );
#else
        if(tauPi0>0.5*dt){
            c->setPiH0(c->getPi() - (c->getPi() - PiNS) * dt / 2.0 / gamma / tauPi0);
            c->setPiH0( (c->getPi_bck() - PiNS0) * dt / 2.0 / gamma / (tauPi0 * tauPi0) * delta_tauPi );
        }
    else
     c->setPiH0(PiNS);
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
          c->addpiH0(i, j, (- deltapipi * ( c->getpi(i, j) * du0 + c->getpi_bck(i,j) * du ) ) / gamma * 0.5 * dt );
          c->addpiH0(i, j, lambdapiPi * ( c->getPi() * sigNS0[i][j] + c->getPi_bck() * sigNS[i][j] ) / gamma * 0.5 * dt);
       for (int k = 0; k < 4; k++) {
//         parts of terms with one internal summation index
        c->addpiH0(i, j, phi70 * ( c->getpi_bck(i, k) * c->getpi(j, k) + c->getpi(i, k) * c->getpi_bck(j, k) ) * gmumu[k] / gamma * 0.5 * dt ); // first part of phi7
        c->addpiH0(i, j, taupipi * 0.5 * ( c->getpi_bck(i, k) * sigNS[j][k] + c->getpi(i, k) * sigNS0[j][k] + c->getpi_bck(j, k) * sigNS[i][k] + c->getpi(j, k) * sigNS0[i][k] ) * gmumu[k] / gamma * 0.5 * dt ); // first part of taupipi
           c->addpiH0(i,j, - u[k] * dpi0[i][j][k] / gamma * 0.5 * dt );
           c->addpiH0(i, j, delta_phi7 * c->getpi_bck(i,k) * c->getpi_bck(j,k) * gmumu[k] / gamma * 0.5 * dt);
//         parts of terms with two internal summation indexes
        for (int l = 0; l < 4; l++){
            c->addpiH0(i, j, - (c->getpi(i, k) * u0[j] + c->getpi(j, k) * u0[i]) * u0[l] * dmu0[l][k] * gmumu[k] / gamma * 0.5 * dt
                             - (c->getpi_bck(i, k) * u0[j] + c->getpi_bck(j, k) * u0[i]) * u0[l] * dmu[l][k] * gmumu[k] / gamma * 0.5 * dt
                             - u[k] * (u0[j] * u0[l] * dpi0[k][i][l] + u0[i] * u0[l] * dpi0[k][j][l]) * gmumu[k] / gamma * 0.5 * dt
                             - (u0[j] * c->getpi_bck(i, k) + u0[i] * c->getpi_bck(j, k)) * u[l] * dmu0[l][k] * gmumu[k] / gamma * 0.5 * dt );
            c->addpiH0(i, j, phi70 * ( ( c->getpi_bck(i, l) * u0[j] + c->getpi_bck(j, l) * u0[i] ) * u0[k] * c->getpi(k, l)
                                      - 2./3. * Delta[index44(i,j)] * (c->getpi_bck(k, l) * c->getpi(k, l) ) ) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt ); // second part of phi7
            c->addpiH0(i, j, phi70 * ( - ( u0[j] * c->getpi_bck(i, l) + u0[i] * c->getpi_bck(j, l) ) * u[k] * c->getpi_bck(k, l)
                                      + 1./3 * ( u[i] * u0[j] + u0[i] * u[j] ) * c->getpi_bck(k, l) * c->getpi_bck(k, l) ) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt ); // delta part of brackets for phi7

            c->addpiH0(i, j, taupipi * ( 0.5 * ( c->getpi_bck(i, l) * u0[j] + c->getpi_bck(j, l) * u0[i] ) * u0[k] * sigNS[k][l]
                                       + 0.5 * ( sigNS0[j][l] * u0[i] + sigNS0[i][l] * u0[j] ) * u0[k] * c->getpi(k, l)
                                       - 1./3. * Delta[index44(i,j)] * (c->getpi_bck(k, l) * sigNS[k][l] + c->getpi(k, l) * sigNS0[k][l] ) ) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt ); // second part of taupipi
            c->addpiH0(i, j, taupipi * ( - 0.5 * ( u0[j] * c->getpi_bck(i, l) + u0[i] * c->getpi_bck(j, l) ) * u[k] * sigNS0[k][l]
                                         - 0.5 * ( u0[i] * sigNS0[j][l] + u0[j] * sigNS0[i][l] ) * u[k] * c->getpi_bck(k, l)
                                         + 1./3 * ( u[i] * u0[j] + u0[i] * u[j] ) * c->getpi_bck(k, l) * sigNS0[k][l] ) * gmumu[k] * gmumu[l] / gamma * 0.5 * dt ); // delta part of brackets for taupipi
            
            c->addPiH0(lamPipi * (c->getpi(k, l) * sigNS0[k][l] + c->getpi_bck(k, l) * sigNS[k][l] ) / gamma * 0.5 * dt);
            for (int r = 0; r < 4; r++) {//3 summation indeces
                c->addpiH0(i, j, 2./3 * (gmunu[i][j] + 2 * u0[i] * u0[j]) * u0[k] * c->getpi(k, r) * u0[l] * dmu0[l][r] * gmumu[k] * gmumu[r] / gamma * 0.5 * dt );//?????????????????????
                c->addpiH0(i, j, ( 1./2 * ( ( u[j] * u0[r] + u0[j] * u[r] ) * ( gmunu[i][k] + u0[i] * u0[k] ) + ( u[i] * u0[k] + u0[i] * u[k] ) * ( gmunu[j][r] + u0[j] * u0[r] ) )
                                 + 1./2 * ( ( u[i] * u0[r] + u0[i] * u[r] ) * ( gmunu[j][k] + u0[j] * u0[k] ) + ( u[j] * u0[k] + u0[j] * u[k] ) * ( gmunu[i][r] + u0[i] * u0[r] ) )
                                 - 1./3 * ( ( u[k] * u0[r] + u0[k] * u[r] ) * ( gmunu[i][j] + u0[i] * u0[j] ) + ( u[i] * u0[j] + u0[i] * u[j] ) * ( gmunu[k][r] + u0[k] * u0[r] ) )
                                 ) * u0[l] * dpi0[k][r][l] * gmumu[k] * gmumu[r] / gamma * 0.5 * dt );
            }
        }
       }
      }
     c->addPiH0(-delPiPi * ( c->getPi() * du0 + c->getPi_bck() * du ) / gamma * 0.5 * dt);
        
     for(int i = 0; i < 4; i++){
         for(int j = 0; j < 4; j++){
             dpi0H[i][j][0] = (f->getCell(ix, iy, iz)->getpiH0_bck(i,j) - f->getCell(ix, iy, iz)->getpiH0_bck_prev(i,j)) / dt;
             dpi0H[i][j][1] = 0.5 * (f->getCell(ix + 1, iy, iz)->getpiH0_bck(i,j) - f->getCell(ix - 1, iy, iz)->getpiH0_bck(i,j)) / dx;
             dpi0H[i][j][2] = 0.5 * (f->getCell(ix, iy + 1, iz)->getpiH0_bck(i,j) - f->getCell(ix, iy - 1, iz)->getpiH0_bck(i,j)) / dy;
             dpi0H[i][j][3] = 0.5 * (f->getCell(ix, iy, iz + 1)->getpiH0_bck(i,j) - f->getCell(ix, iy, iz - 1)->getpiH0_bck(i,j)) / dz / tauPlusHalf_0;
         }
     }
        
     // 1) relaxation(piH)+source(piH) terms for full-step
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
#ifdef FORMAL_SOLUTION
       c->setpi0(i, j,
                 (c->getpi(i, j) - piNS[i][j]) * exp(-dt / gamma / taupi0) +
                     piNS[i][j]);
       c->setpi0(i, j, (c->getpiH0_bck(i, j) - piNS0[i][j]) * dt / gamma / (taupi0*taupi0) * delta_taupi );//H0????????????????????,

#else
      if(taupi0>0.5*dt){
              c->setpi0(i, j, c->getpi(i, j) -
                        (c->getpiH0(i, j) - piNS[i][j]) * dt / gamma / taupi0);
              c->setpi0(i, j, (c->getpiH0_bck(i, j) - piNS0[i][j]) * dt / gamma / (taupi0*taupi0) * delta_taupi );
          }
      else
       c->setpi0(i, j, piNS[i][j]);
#endif
//          cout << c->getpi(i,j) << endl;
      }
        
#ifdef FORMAL_SOLUTION
     c->setPi0((c->getPi() - PiNS) * exp(-dt / gamma / tauPi0) + PiNS);
     c->setPi0( (c->getPiH0_bck() - PiNS0) * dt / gamma / (tauPi0 * tauPi0) * delta_tauPi );
#else
        if(tauPi0>0.5*dt){
            c->setPi0(c->getPi() - (c->getPiH0() - PiNS) * dt / gamma / tauPi0);
            c->setPi0( (c->getPiH0_bck() - PiNS0) * dt / gamma / (tauPi0 * tauPi0) * delta_tauPi );
        }
    else
     c->setPi0(PiNS);
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
             c->addpi0(i, j, (- deltapipi * ( c->getpiH0(i, j) * du0 + c->getpiH0_bck(i,j) * du ) ) / gamma * dt );
             c->addpi0(i, j, lambdapiPi * ( c->getPiH0() * sigNS0[i][j] + c->getPiH0_bck() * sigNS[i][j] ) / gamma * dt);
          for (int k = 0; k < 4; k++) {
    //         parts of terms with one internal summation index
              c->addpi0(i, j, phi70 * ( c->getpiH0_bck(i, k) * c->getpiH0(j, k) + c->getpiH0(i, k) * c->getpiH0_bck(j, k) ) * gmumu[k] / gamma * dt ); // first part of phi7
              c->addpi0(i, j, taupipi * 0.5 * ( c->getpiH0_bck(i, k) * sigNS[j][k] + c->getpiH0(i, k) * sigNS0[j][k] + c->getpiH0_bck(j, k) * sigNS[i][k] + c->getpiH0(j, k) * sigNS0[i][k] ) * gmumu[k] / gamma * dt ); // first part of taupipi
              c->addpi0(i, j, - u[k] * dpi0H[i][j][k] / gamma * dt );
              
              c->addpi0(i, j, delta_phi7 * c->getpiH0_bck(i,k) * c->getpiH0_bck(j,k) * gmumu[k] / gamma * dt);
              
    //         parts of terms with two internal summation indexes
           for (int l = 0; l < 4; l++){
               c->addpi0(i, j, - (c->getpiH0(i, k) * u0[j] + c->getpiH0(j, k) * u0[i]) * u0[l] * dmu0[l][k] * gmumu[k] / gamma * dt
                            - (c->getpiH0_bck(i, k) * u0[j] + c->getpiH0_bck(j, k) * u0[i]) * u0[l] * dmu[l][k] * gmumu[k] / gamma * dt
                            - u[k] * (u0[j] * u0[l] * dpi0H[k][i][l] + u0[i] * u0[l] * dpi0H[k][j][l]) * gmumu[k] / gamma * dt //?????????????ten faktor?
                            - (u0[j] * c->getpiH0_bck(i, k) + u0[i] * c->getpiH0_bck(j, k)) * u[l] * dmu0[l][k] * gmumu[k] / gamma * dt );
               c->addpi0(i, j, phi70 * ( ( c->getpiH0_bck(i, l) * u0[j] + c->getpiH0_bck(j, l) * u0[i] ) * u0[k] * c->getpiH0(k, l)
                                        - 2./3. * Delta[index44(i,j)] * (c->getpiH0_bck(k, l) * c->getpiH0(k, l) ) ) * gmumu[k] * gmumu[l] / gamma * dt ); // second part of phi7
               c->addpi0(i, j, phi70 * ( - ( u0[j] * c->getpiH0_bck(i, l) + u0[i] * c->getpiH0_bck(j, l) ) * u[k] * c->getpiH0_bck(k, l)
                                        + 1./3 * ( u[i] * u0[j] + u0[i] * u[j] ) * c->getpiH0_bck(k, l) * c->getpiH0_bck(k, l) ) * gmumu[k] * gmumu[l] / gamma * dt ); // delta part of brackets for phi7

               c->addpi0(i, j, taupipi * ( 0.5 * ( c->getpiH0_bck(i, l) * u0[j] + c->getpiH0_bck(j, l) * u0[i] ) * u0[k] * sigNS[k][l]
                                        + 0.5 * ( sigNS0[j][l] * u0[i] + sigNS0[i][l] * u0[j] ) * u0[k] * c->getpiH0(k, l)
                                        - 1./3. * Delta[index44(i,j)] * (c->getpiH0_bck(k, l) * sigNS[k][l] + c->getpiH0(k, l) * sigNS0[k][l] ) ) * gmumu[k] * gmumu[l] / gamma * dt ); // second part of taupipi
               c->addpi0(i, j, taupipi * ( - 0.5 * ( u0[j] * c->getpiH0_bck(i, l) + u0[i] * c->getpiH0_bck(j, l) ) * u[k] * sigNS0[k][l]
                                        - 0.5 * ( u0[i] * sigNS0[j][l] + u0[j] * sigNS0[i][l] ) * u[k] * c->getpiH0_bck(k, l)
                                        + 1./3 * ( u[i] * u0[j] + u0[i] * u[j] ) * c->getpiH0_bck(k, l) * sigNS0[k][l] ) * gmumu[k] * gmumu[l] / gamma * dt ); // delta part of brackets for taupipi

               c->addpi0(i, j, - 1. / 3. * Delta[index44(i,j)] * c->getpiH0_bck(k, l) * delta_phi7 * c->getpiH0_bck(k, l) / gamma * 0.5 * dt);
               
               c->addPi0(lamPipi * (c->getpiH0(k, l) * sigNS0[k][l] + c->getpiH0_bck(k, l) * sigNS[k][l] ) / gamma * dt);
               for (int r = 0; r < 4; r++) {
                   c->addpi0(i, j, 2./3 * (gmunu[i][j] + 2 * u0[i] * u0[j]) * u0[k] * c->getpiH0(k, r) * u0[l] * dmu0[l][r] * gmumu[k] * gmumu[r] / gamma * dt );//?????????????????????
                   c->addpi0(i, j, ( 1./2 * ( ( u[j] * u0[r] + u0[j] * u[r] ) * ( gmunu[i][k] + u0[i] * u0[k] ) + ( u[i] * u0[k] + u0[i] * u[k] ) * ( gmunu[j][r] + u0[j] * u0[r] ) )
                                + 1./2 * ( ( u[i] * u0[r] + u0[i] * u[r] ) * ( gmunu[j][k] + u0[j] * u0[k] ) + ( u[j] * u0[k] + u0[j] * u[k] ) * ( gmunu[i][r] + u0[i] * u0[r] ) )
                                - 1./3 * ( ( u[k] * u0[r] + u0[k] * u[r] ) * ( gmunu[i][j] + u0[i] * u0[j] ) + ( u[i] * u0[j] + u0[i] * u[j] ) * ( gmunu[k][r] + u0[k] * u0[r] ) )
                                ) * u0[l] * dpi0H[k][r][l] * gmumu[k] * gmumu[r] / gamma * dt );
               }
           }
          }
         }
        c->addPi0(-delPiPi * ( c->getPiH0() * du0 + c->getPiH0_bck() * du ) / gamma * dt);
    }  // end non-empty cell
   }   // end loop #1
    
    
 // 3) -- advection ---
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenterQ0(eos, tauMinusHalf, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
    c->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vx, vy,
                         vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);  // getPrimVar() before
       
    if (e_0 <= 0.) continue;
    double xm = -vx_0 * dt / f->getDx();
    double ym = -vy_0 * dt / f->getDy();
    double zm = -vz_0 * dt / f->getDz() / tauMinusHalf;
    double xmH = -vx_0 * dt / f->getDx() / 2.0;
    double ymH = -vy_0 * dt / f->getDy() / 2.0;
    double zmH = -vz_0 * dt / f->getDz() / tauMinusHalf / 2.0;
    double wx[2] = {(1. - fabs(xm)), fabs(xm)};
    double wy[2] = {(1. - fabs(ym)), fabs(ym)};
    double wz[2] = {(1. - fabs(zm)), fabs(zm)};
    double wxH[2] = {(1. - fabs(xmH)), fabs(xmH)};
    double wyH[2] = {(1. - fabs(ymH)), fabs(ymH)};
    double wzH[2] = {(1. - fabs(zmH)), fabs(zmH)};
       
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) {
      pi[i][j] = piH[i][j] = 0.0;
     }
    Pi = PiH = 0.0;
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       // pi,Pi-->full step, piH,PiH-->half-step
       Cell *c1 = f->getCell(ix + jx * sign(xm), iy + jy * sign(ym),
                             iz + jz * sign(zm));
       for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
         pi[i][j] += wx[jx] * wy[jy] * wz[jz] * c1->getpi0(i, j); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!????????????????.....................
         piH[i][j] += wxH[jx] * wyH[jy] * wzH[jz] * c1->getpiH0(i, j);
        }
        Pi += wx[jx] * wy[jy] * wz[jz] * c1->getPi0();
        PiH += wxH[jx] * wyH[jy] * wzH[jz] * c1->getPiH0();
      }
    //--------- debug part: trace check, diag check, transversality check
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) {
         //cout << pi[i][j] << endl;
      if (pi[i][j] != 0. && fabs(1.0 - pi[j][i] / pi[i][j]) > 1e-10)
       cout << "non-diag: " << pi[i][j] << "  " << pi[j][i] << endl;
     }
    //------ end debug
    //======= hydro applicability check (viscous corrections limiter):
    // double maxT0 = max(fabs((e+p)*vx*vx/(1.-vx*vx-vy*vy-vz*vz)+p),
    //   fabs((e+p)*vy*vy/(1.-vx*vx-vy*vy-vz*vz)+p)) ;
    double maxT0 = max((e_0 + p_0) / (1. - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0) - p_0,
                       (e_0 + p_0) * (vx_0 * vx_0 + vy_0 * vy_0 + vz_0 * vz_0) /
                               (1. - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0) +
                           p_0);
    // double maxpi = max(fabs(pi[1][1]),fabs(pi[2][2])) ;
    double maxpi = 0.;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++){
            if (fabs(pi[i][j]) > maxpi) maxpi = fabs(pi[i][j]); // tohle asi taky ne ..... porovna va se pozadi s fluktuaci
        }
    bool rescaled = false;
    if (maxT0 / maxpi < 1.0) {
     for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
       pi[i][j] = 0.1 * pi[i][j] * maxT0 / maxpi;
       piH[i][j] = 0.1 * piH[i][j] * maxT0 / maxpi;
      }
     rescaled = true;
    }
    if (fabs(Pi) > p) {
     if (Pi != 0.) Pi = 1 * Pi / fabs(Pi) * p;     //modified .1 -> 1
     if (PiH != 0.) PiH = 1 * PiH / fabs(PiH) * p; //modified .1 -> 1
     rescaled = true;
    }
    if (rescaled)
     c->setViscCorrCutFlag(maxT0 / maxpi);
    else
     c->setViscCorrCutFlag(1.);
    // updating to the new values
    for (int i = 0; i < 4; i++)
     for (int j = 0; j <= i; j++) {
//         if(i==1&&j==1){
//             c->setpi(i, j, 0.);
//             c->setpiH(i, j, 0.);
//         }
//         else{
             c->setpi(i, j, pi[i][j]);
             c->setpiH(i, j, piH[i][j]);
//         }
//         cout << pi[i][j] << "  " << i << " " << j << endl;
     }
    c->setPi(Pi);
    c->setPiH(PiH);
    // source term  - (tau+dt)*delta_Q_(i+1)/delta_tau
       
    double gamma = 1.0 / sqrt(1.0 - vx_0 * vx_0 - vy_0 * vy_0 - vz_0 * vz_0);
    double u0[4];
    double u[4];
    u0[0] = gamma;
    u0[1] = u0[0] * vx_0;
    u0[2] = u0[0] * vy_0;
    u0[3] = u0[0] * vz_0;
       
       u[0] = 0.;
       u[1] = gamma * vx;
       u[2] = gamma * vy;
       u[3] = gamma * vz;
            
    double flux[4];
       for (int i = 0; i < 4; i++){
           flux[i] = -tau * (c->getpi(0, i) + c->getPi() * ( u0[0] * u[i] + u[0] * u0[i] ) );
       }
    flux[0] += tau * c->getPi();
    c->addFlux(flux[0], flux[1], flux[2], flux[3], 0., 0., 0.);
   }  // advection loop (all cells)
}

// this procedure explicitly uses T_==0, X_==1, Y_==2, Z_==3
void Hydro::visc_flux(Cell *left, Cell *right, int direction, double ix, double iy, double iz) {
 double flux[4];
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
 double e, p, nb, nq, ns, vxl, vyl, vzl, vxr, vyr, vzr;
 double  el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0;
 double  er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0;
 // we need to know the velocities at both cell centers at (n+1/2) in order to
 // interpolate to
 // get the value at the interface
 left->getPrimVarHCenterQ0(eos, tauMinusHalf, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0);
 left->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vxl, vyl, vzl, el_0, pl_0, nbl_0, nql_0, nsl_0, vxl_0, vyl_0, vzl_0);
 right->getPrimVarHCenterQ0(eos, tauMinusHalf, er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0);
 right->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vxr, vyr, vzr, er_0, pr_0, nbr_0, nqr_0, nsr_0, vxr_0, vyr_0, vzr_0);
 vxl = 0.5 * (vxl + vxr);
 vyl = 0.5 * (vyl + vyr);
 vzl = 0.5 * (vzl + vzr);
 
    vxl_0 = 0.5 * (vxl_0 + vxr_0);
    vyl_0 = 0.5 * (vyl_0 + vyr_0);
    vzl_0 = 0.5 * (vzl_0 + vzr_0);
    
 double v0 = sqrt(vxl_0 * vxl_0 + vyl_0 * vyl_0 + vzl_0 * vzl_0);//??????????????????????????????????????????????????????????????????????????????
 double v = sqrt(vxl_0 * vxl_0 + vyl_0 * vyl_0 + vzl_0 * vzl_0 + 2 * vxl_0 * vxl + 2 * vxl_0 * vxl + 2 * vxl_0 * vxl);
// if (v > 1.) {
//     vxl = 0.99 * vxl / v;//or 0???????????????????????????????????????????????????????????????????????????????????????????????????????????????????
//     vyl = 0.99 * vyl / v;
//     vzl = 0.99 * vzl / v;
// }
    
 double gamma = 1. / sqrt(1. - v0 * v0);
 double uuu0[4] = {gamma, gamma * vxl_0, gamma * vyl_0, gamma * vzl_0};
 double gammal_delta = 0.;
 double uuu[4] = {0., gamma * vxl, gamma * vyl, gamma * vzl};
    
 double gmumu[4] = {1., -1., -1., -1.};
 if (direction == X_)
  ind2 = 1;
 else if (direction == Y_)
  ind2 = 2;
 else if (direction == Z_)
  ind2 = 3;
 for (int ind1 = 0; ind1 < 4; ind1++) {
  flux[ind1] = 0.5 * (left->getpiH(ind1, ind2) + right->getpiH(ind1, ind2));
//     cout << left->getpiH(ind1, ind2) << "      " << right->getpiH(ind1, ind2) << endl;
 if (ind1 == ind2){
     flux[ind1] += -0.5 * (left->getPiH() + right->getPiH()) *
     gmumu[ind1];  // gmunu is diagonal
 }
     flux[ind1] += 0.5 * (left->getPiH() + right->getPiH()) * uuu0[ind1] * uuu0[ind2]
             + 0.5 * (left->getPiH_bck() + right->getPiH_bck()) * ( uuu0[ind1] * uuu[ind2] + uuu[ind1] * uuu0[ind2] );
 }
 for (int i = 0; i < 4; i++) flux[i] = flux[i] * tauMinusHalf * dt / dxa;
//    if ((ix==20 || ix ==19 || ix ==21 || ix ==60 || ix ==61 || ix ==59)&&iy==0&&iz==0) {
//        cout << ix << " x: " << flux[X_] << "   y: " << flux[Y_] << "  z:  " << flux[Z_] << endl;
//    }
 left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
 right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
}

void Hydro::performStep(void) {
 // debugRiemann = false ; // turn off debug output

 f->updateM(tau, dt);

 tau_z = dt / 2. / log(1 + dt / 2. / tau);

//    // vypis profilu energie
    for (int ix = 0; ix < f->getNX(); ix++) {
//     for (int iy = 0; iy < f->getNY(); iy++){
        
        double e, p, nb, nq, ns, vx, vy, vz, cs, T, mub, muq, mus;
        double e__0, p__0, nb__0, nq__0, ns__0, vx__0, vy__0, vz__0, cs__0, T__0, mub__0, muq__0, mus__0;
        double e_H, p_H, nb_H, nq_H, ns_H, vx_H, vy_H, vz_H, cs_H, T_H, mub_H, muq_H, mus_H;
        double eH, pH, nbH, nqH, nsH, vxH, vyH, vzH;

         
        Cell *c = f->getCell(ix, 0, 0);
        
        c -> getPrimVarQ0(eos, tau, e__0, p__0, nb__0, nq__0, ns__0, vx__0, vy__0, vz__0);
        c -> getPrimVarHCenterQ0(eos, tau, e_H, p_H, nb_H, nq_H, ns_H, vx_H, vy_H, vz_H);
        c -> getPrimVarHCenter(eos, tau, eH, pH, nbH, nqH, nsH, vxH, vyH, vzH, e_H, p_H, nb_H, nq_H, ns_H, vx_H, vy_H, vz_H);
        c -> getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz, e__0, p__0, nb__0, nq__0, ns__0, vx__0, vy__0, vz__0);

        cs = eos->cs();
//        eos->eos(e, 0., 0., 0., T, mub, mus, muq, p);
        
        if(ix<=40){
            cout.precision(14);
        }
        else{
            cout.precision(15);
        }
        
//                      if(ix==20 || ix ==19 || ix ==21 || ix ==60 || ix ==61 || ix ==59){
//        if(ix==39 || ix== 40 || ix==41){
        if(ix==40){
//                  cout << vx+vx__0 << "   " << f->getX(ix) << "      " << e+e__0 << "       " << p+p__0 << "     " << cs << endl;
//                  cout << vx << "   " << vy << "    " << vz << "   " << f->getX(ix) << "      " << e+e__0 << "       " << p+p__0 << "     " << cs << endl;
                  cout << vx << "   " << f->getX(ix) << "      " << e+e__0 << "       " << p+p__0 << "     " << cs << endl;


              }

          
//              if(ix==0){
//                  cout << vx << "     " << p << "    " << ix << "    " << f->getX(ix) << "      " << e << "       " << cs << "     " << T << endl;
//          if(ix<40){
//              cout << f->getX(ix) << "      " << e - 1. << endl;
//          }
//          else{
//              cout << f->getX(ix) << "      " << e - 1. << endl;
//          }
//              }
//          }
      }
    
    //cout << "E =        " << e_const;//*(f->getNY())*(f->getNX())*(f->getNZ()) << endl;
   //cout << 1. / sqrt(1 - vx_const * vx_const - vy_const * vy_const - vz_const * vz_const) << endl;
   // cout << "TADY:     " << vxH_const << "      " << vxPrev_const << endl;
    
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
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix + 1, iy, iz), X_, PREDICT, ix);
   }
 //	cout << "predictor X done\n" ;
 // Y dir
 for (int iz = 0; iz < f->getNZ(); iz++)
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, PREDICT, ix);
   }
 //	cout << "predictor Y done\n" ;
 // Z dir
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz + 1), Z_, PREDICT, ix);
   }
 //	cout << "predictor Z done\n" ;

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
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix + 1, iy, iz), X_, CORRECT, ix);
   }
 //	cout << "corrector X done\n" ;
 // Y dir
 for (int iz = 0; iz < f->getNZ(); iz++)
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, CORRECT, ix);
   }
 //	cout << "corrector Y done\n" ;
 // Z dir
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz + 1), Z_, CORRECT, ix);
   }
 //	cout << "corrector Z done\n" ;

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
  //	cout << "visc_flux X done\n" ;
  // Y dir
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++) {
//        cout << ix << "     " << iy << "    " << iz << endl;
     visc_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, ix, iy, iz);
    }
  //	cout << "visc_flux Y done\n" ;
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
     visc_source_step(ix, iy, iz);
     f->getCell(ix, iy, iz)->updateByViscFlux();
     f->getCell(ix, iy, iz)->clearFlux();
    }
 } else {  // end viscous part
 }
 //==== finishing work ====
 //f->correctImagCellsFull();  // disabled in box mode
}
