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

void Hydro::hlle_flux(Cell *left, Cell *right, int direction, int mode, double e_const, double p_const, double vx_const, double vy_const, double vz_const, int ix) {
 // for all variables, suffix "l" = left state, "r" = right state
 // with respect to the cell boundary
 double el, er, pl, pr, nbl, nql, nsl, nbr, nqr, nsr, vxl, vxr, vyl, vyr, vzl,
     vzr, bl = 0., br = 0., csb, vb, vb_const, vb_delta, El, Er, dx = 0.;
 double Ftl = 0., Fxl = 0., Fyl = 0., Fzl = 0., Fbl = 0., Fql = 0., Fsl = 0.,
        Ftr = 0., Fxr = 0., Fyr = 0., Fzr = 0., Fbr = 0., Fqr = 0., Fsr = 0.;
 double U1l, U2l, U3l, U4l, Ubl, Uql, Usl, U1r, U2r, U3r, U4r, Ubr, Uqr, Usr;
 double flux_T_UU, flux_X_UU, flux_Y_UU, flux_Z_UU, flux_NB_UU, flux_NQ_UU, flux_NS_UU, flux_T_Fr, flux_X_Fr, flux_Y_Fr, flux_Z_Fr, flux_NB_Fr, flux_NQ_Fr, flux_NS_Fr, flux_T_Fl, flux_X_Fl, flux_Y_Fl, flux_Z_Fl, flux_NB_Fl, flux_NQ_Fl, flux_NS_Fl;
 double v0_L, v0_R, delta_vL, delta_vR, cs;
 double br_const, bl_const, P, Q, C, B;
 double flux[7];
 const double dta = mode == 0 ? dt / 2. : dt;
 double tauFactor = 1.0;  // fluxes are also multiplied by tau, 1 for Cartesian
 if (mode == PREDICT) {
  // get primitive quantities from Q_{i+} at previous timestep
  left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                        direction);
  // ... and Q_{(i+1)-}
  right->getPrimVarLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                        direction);
     double deltaVxl_ = vxl - vx_const;
     double deltaVyl_ = vyl - vy_const;
     double deltaVzl_ = vzl - vz_const;
     double deltaVxr_ = vxr - vx_const;
     double deltaVyr_ = vyr - vy_const;
     double deltaVzr_ = vzr - vz_const;

//  El = (el + pl) / (1 - vxl * vxl - vyl * vyl - vzl * vzl);
//  Er = (er + pr) / (1 - vxr * vxr - vyr * vyr - vzr * vzr);
  El = (el + pl) * (1 + vx_const * vx_const + vy_const * vy_const + vz_const * vz_const) + 2 * ( e_const + p_const ) * ( vx_const * deltaVxl_ + vy_const * deltaVyl_ + vz_const * deltaVzl_ );
  Er = (er + pr) * (1 + vx_const * vx_const + vy_const * vy_const + vz_const * vz_const) + 2 * ( e_const + p_const ) * ( vx_const * deltaVxr_ + vy_const * deltaVyr_ + vz_const * deltaVzr_ );
  #ifndef CARTESIAN
  tauFactor = tau + 0.25 * dt;
  #endif
     //cout << "nbl=    " << nbl << "   nbr=      " << nbr << endl;

 } else {
  // use half-step updated Q's for corrector step
  left->getPrimVarHRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                         direction);
  right->getPrimVarHLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                         direction);

     double deltaVxl_ = vxl - vx_const;
     double deltaVyl_ = vyl - vy_const;
     double deltaVzl_ = vzl - vz_const;
     double deltaVxr_ = vxr - vx_const;
     double deltaVyr_ = vyr - vy_const;
     double deltaVzr_ = vzr - vz_const;

//  El = (el + pl) / (1 - vxl * vxl - vyl * vyl - vzl * vzl);
//  Er = (er + pr) / (1 - vxr * vxr - vyr * vyr - vzr * vzr);
     double E_const = (e_const + p_const) * (1 + vx_const * vx_const + vy_const * vy_const + vz_const * vz_const);
     double El_delta = (el - e_const + pl - p_const) * (1 + vx_const * vx_const + vy_const * vy_const + vz_const * vz_const) + 2 * ( e_const + p_const ) * ( vx_const * deltaVxl_ + vy_const * deltaVyl_ + vz_const * deltaVzl_ );
     double Er_delta = (er - e_const + pr - p_const) * (1 + vx_const * vx_const + vy_const * vy_const + vz_const * vz_const) + 2 * ( e_const + p_const ) * ( vx_const * deltaVxr_ + vy_const * deltaVyr_ + vz_const * deltaVzr_ );
  El = E_const + El_delta;
  Er = E_const + Er_delta;
  #ifndef CARTESIAN
  tauFactor = tau + 0.5 * dt;
  #endif
 }

 if (el < 0.) {
  el = 0.;
  pl = 0.;
 }
 if (er < 0.) {
  er = 0.;
  pr = 0.;
 }

 if (el > 1e10) {
  cout << "e>1e10; debug info below:\n";
  left->Dump(tau);
  // debugRiemann = true ;
  if (mode == PREDICT)
   left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                         direction);
  else
   left->getPrimVarHRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                          direction);
  // debugRiemann = false ;
  exit(0);
 }

 // skip the procedure for two empty cells
 if (el == 0. && er == 0.) return;
 if (pr < 0.) {
  cout << "Negative pressure" << endl;
  left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl,
                        direction);
  right->getPrimVarLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr,
                        direction);
 }

 // skip the procedure for two partially vacuum cells
 if (left->getM(direction) < 1. && right->getM(direction) < 1.) return;

// double gammal = 1. / sqrt(1 - vxl * vxl - vyl * vyl - vzl * vzl);//potreba linearizovat!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// double gammar = 1. / sqrt(1 - vxr * vxr - vyr * vyr - vzr * vzr);

    // idea is that e.g. el = el0 + deltaEl
//    double el0 = 0, er0 = 0, pl0 = 0, pr0 = 0, vxl0 = 0, vxr0 = 0, vyl0 = 0, vyr0 = 0, vzl0 = 0, vzr0 = 0;//hodnoty spocitane jako prumer vsech bunek pri kazdem casovem kroku
//fluktuace - promenne - v podstate by to melo byt e.g. deltaEl = (el-el0)
    double deltaVxl = vxl - vx_const;
    double deltaVyl = vyl - vy_const;
    double deltaVzl = vzl - vz_const;
//    double deltaVtl = 0.01;
    //double deltaVtl = vtl - vtl0;
    double deltaVxr = vxr - vx_const;
    double deltaVyr = vyr - vy_const;
    double deltaVzr = vzr - vz_const;
//    double deltaVtr = 0.01;
    //double deltaVtr = vtr - vtr0;
    double deltaEl = el - e_const;
    double deltaEr = er - e_const;
    double deltaPl = pl - p_const;
    double deltaPr = pr - p_const;

    double gamma_const = 1. + 0.5 * (vx_const * vx_const + vy_const * vy_const + vz_const * vz_const);
    double gammal_delta = vx_const * deltaVxl + vy_const * deltaVyl + vz_const * deltaVzl;
    double gammar_delta = vx_const * deltaVxr + vy_const * deltaVyr + vz_const * deltaVzr;

    double gammal = gamma_const + gammal_delta;
    double gammar = gamma_const + gammar_delta;

//    //?????????
//    double deltaVtl = gammal - 1. / sqrt(1 - vx_const * vx_const - vy_const * vy_const - vz_const * vz_const);
//    double deltaVtr = gammar - 1. / sqrt(1 - vx_const * vx_const - vy_const * vy_const - vz_const * vz_const);// tohle je asi v pohode??..........................................................................
//    //?????????

    double deltaVtl = gammal - 1. - 0.5 * (vx_const * vx_const + vy_const * vy_const + vz_const * vz_const);
    double deltaVtr = gammar - 1. - 0.5 * (vx_const * vx_const + vy_const * vy_const + vz_const * vz_const);

    //cout << deltaVtl << endl;

    double U1l_const_v, U1r_const_v, U2l_const_v, U2r_const_v, U3l_const_v, U3r_const_v, U4l_const_v, U4r_const_v, U1l_flux, U1r_flux, U2l_flux, U2r_flux, U3l_flux, U3r_flux, U4l_flux, U4r_flux;// U's to be multiplied by constant component of velocity
    double U1l_delta_v, U1r_delta_v, U2l_delta_v, U2r_delta_v, U3l_delta_v, U3r_delta_v, U4l_delta_v, U4r_delta_v, U1r_flux_delta, U1l_flux_delta, U2l_flux_delta, U2r_flux_delta, U3l_flux_delta, U3r_flux_delta, U4l_flux_delta, U4r_flux_delta;// U's to be multiplied by the fluctuation term of velocity

    double el0 = el - deltaEl;
    double pl0 = pl - deltaPl;
    double vxl0 = vxl - deltaVxl;
    double vyl0 = vyl - deltaVyl;
    double vzl0 = vzl - deltaVzl;

    double er0 = er - deltaEr;
    double pr0 = pr - deltaPr;
    double vxr0 = vxr - deltaVxr;
    double vyr0 = vyr - deltaVyr;
    double vzr0 = vzr - deltaVzr;

 U1l_const_v = gamma_const * (gammal + gammal_delta) * (el0 + pl0) * vxl0 + gamma_const * gamma_const * ( (el0 + pl0) * deltaVxl + (deltaEl + deltaPl) * vxl0 );
 U2l_const_v = gamma_const * (gammal + gammal_delta) * (el0 + pl0) * vyl0 + gamma_const * gamma_const * ( (el0 + pl0) * deltaVyl + (deltaEl + deltaPl) * vyl0 );
 U3l_const_v = gamma_const * (gammal + gammal_delta) * (el0 + pl0) * vzl0 + gamma_const * gamma_const * ( (el0 + pl0) * deltaVzl + (deltaEl + deltaPl) * vzl0 );
 U4l_const_v = gamma_const * (gammal + gammal_delta) * (el0 + pl0) + gamma_const * gamma_const * ( (el0 + pl0) * deltaVtl + (deltaEl + deltaPl) ); //POZOR NA TLAK! - pl

    U1l_flux = gamma_const * gamma_const * (el0 + pl0) * vxl0;
    U2l_flux = gamma_const * gamma_const * (el0 + pl0) * vyl0;
    U3l_flux = gamma_const * gamma_const * (el0 + pl0) * vzl0;
    U4l_flux = gamma_const * gamma_const * (el0 + pl0);

    U1l_flux_delta = 2 * gamma_const * gammal_delta * (el0 + pl0) * vxl0 + gamma_const * gamma_const * ( (el0 + pl0) * deltaVxl + (deltaEl + deltaPl) * vxl0 );
    U2l_flux_delta = 2 * gamma_const * gammal_delta * (el0 + pl0) * vyl0 + gamma_const * gamma_const * ( (el0 + pl0) * deltaVyl + (deltaEl + deltaPl) * vyl0 );
    U3l_flux_delta = 2 * gamma_const * gammal_delta * (el0 + pl0) * vzl0 + gamma_const * gamma_const * ( (el0 + pl0) * deltaVzl + (deltaEl + deltaPl) * vzl0 );
    U4l_flux_delta = 2 * gamma_const * gammal_delta * (el0 + pl0) + gamma_const * gamma_const * ( (el0 + pl0) * deltaVtl + (deltaEl + deltaPl) );

 U1l_delta_v = gamma_const * gamma_const * (el0 + pl0) * vxl0;
 U2l_delta_v = gamma_const * gamma_const * (el0 + pl0) * vyl0;
 U3l_delta_v = gamma_const * gamma_const * (el0 + pl0) * vzl0;
 U4l_delta_v = gamma_const * gamma_const * (el0 + pl0);

 U1l = U1l_const_v + U1l_delta_v ;
 U2l = U2l_const_v + U2l_delta_v ;
 U3l = U3l_const_v + U3l_delta_v ;
 U4l = U4l_const_v + U4l_delta_v  - pl; // nebo to rozdelit jako (p0+deltaP)???
 Ubl = gammal * nbl;
 Uql = gammal * nql;
 Usl = gammal * nsl;

 U1r_const_v = gamma_const * (gammar + gammar_delta) * (er0 + pr0) * vxr0 + gamma_const * gamma_const * ( (er0 + pr0) * deltaVxr + (deltaEr + deltaPr) * vxr0 );
 U2r_const_v = gamma_const * (gammar + gammar_delta) * (er0 + pr0) * vyr0 + gamma_const * gamma_const * ( (er0 + pr0) * deltaVyr + (deltaEr + deltaPr) * vyr0 );
 U3r_const_v = gamma_const * (gammar + gammar_delta) * (er0 + pr0) * vzr0 + gamma_const * gamma_const * ( (er0 + pr0) * deltaVzr + (deltaEr + deltaPr) * vzr0 );
 U4r_const_v = gamma_const * (gammar + gammar_delta) * (er0 + pr0) + gamma_const * gamma_const * ( (er0 + pr0) * deltaVtr + (deltaEr + deltaPr) ); //POZOR NA TLAK! - pl

    U1r_flux = gamma_const * gamma_const * (er0 + pr0) * vxr0;
    U2r_flux = gamma_const * gamma_const * (er0 + pr0) * vyr0;
    U3r_flux = gamma_const * gamma_const * (er0 + pr0) * vzr0;
    U4r_flux = gamma_const * gamma_const * (er0 + pr0);

    U1r_flux_delta = 2 * gamma_const * gammar_delta * (er0 + pr0) * vxr0 + gamma_const * gamma_const * ( (er0 + pr0) * deltaVxr + (deltaEr + deltaPr) * vxr0 );
    U2r_flux_delta = 2 * gamma_const * gammar_delta * (er0 + pr0) * vyr0 + gamma_const * gamma_const * ( (er0 + pr0) * deltaVyr + (deltaEr + deltaPr) * vyr0 );
    U3r_flux_delta = 2 * gamma_const * gammar_delta * (er0 + pr0) * vzr0 + gamma_const * gamma_const * ( (er0 + pr0) * deltaVzr + (deltaEr + deltaPr) * vzr0 );
    U4r_flux_delta = 2 * gamma_const * gammar_delta * (er0 + pr0) + gamma_const * gamma_const * ( (er0 + pr0) * deltaVtr + (deltaEr + deltaPr) );

 U1r_delta_v = gamma_const * gamma_const * (er0 + pr0) * vxr0;
 U2r_delta_v = gamma_const * gamma_const * (er0 + pr0) * vyr0;
 U3r_delta_v = gamma_const * gamma_const * (er0 + pr0) * vzr0;
 U4r_delta_v = gamma_const * gamma_const * (er0 + pr0);

 U1r = U1r_const_v + U1r_delta_v ;
 U2r = U2r_const_v + U2r_delta_v ;
 U3r = U3r_const_v + U3r_delta_v ;
 U4r = U4r_const_v + U4r_delta_v  - pr; // nebo to rozdelit jako (p0+deltaP)???
 Ubr = gammar * nbr;
 Uqr = gammar * nqr;
 Usr = gammar * nsr;

 if (direction == X_) {
  Ftl = U4l_const_v * vxl0 + U4l_delta_v * deltaVxl;
  Fxl = U1l_const_v * vxl0 + U1l_delta_v * deltaVxl + pl;
  Fyl = U2l_const_v * vxl0 + U2l_delta_v * deltaVxl;
  Fzl = U3l_const_v * vxl0 + U3l_delta_v * deltaVxl;
  Fbl = Ubl * vxl;
  Fql = Uql * vxl;
  Fsl = Usl * vxl;

  Ftr = U4r_const_v * vxr0 + U4r_delta_v * deltaVxr;
  Fxr = U1r_const_v * vxr0 + U1r_delta_v * deltaVxr + pr;
  Fyr = U2r_const_v * vxr0 + U2r_delta_v * deltaVxr;
  Fzr = U3r_const_v * vxr0 + U3r_delta_v * deltaVxr;
  Fbr = Ubr * vxr;
  Fqr = Uqr * vxr;
  Fsr = Usr * vxr;

  // for the case of constant c_s only
//  csb = sqrt(eos->cs2() +
//             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
//                 pow(vxl - vxr, 2));
//  vb = (sqrt(El) * vxl + sqrt(Er) * vxr) / (sqrt(El) + sqrt(Er));
//  bl = min(0., min((vb - csb) / (1 - vb * csb),
//                      (vxl - eos->cs()) / (1 - vxl * eos->cs())));
//  br = max(0., max((vb + csb) / (1 + vb * csb),
//                      (vxr + eos->cs()) / (1 + vxr * eos->cs())));
     csb = sqrt(eos->cs2());
     vb_const = vx_const;
     vb_delta = 0.5 * (deltaVxl + deltaVxr);
     vb = vb_const + vb_delta;
     bl = min(0., min(1 / (1 - vb_const * csb) * (vb_const - csb + (1 + csb * (vb_const - csb) / (1 - vb_const * csb) ) * vb_delta ),
                      1 / (1 - vx_const * eos->cs()) * (vx_const - eos->cs() + (eos->cs() * (vx_const - eos->cs()) / (1 - vx_const * eos->cs()) + 1 ) * deltaVxl ) ) );
     br = max(0., max(1 / (1 + vb_const * csb) * (vb_const + csb - (1 - csb * (vb_const + csb) / (1 + vb_const * csb) ) * vb_delta ),
                      1 / (1 + vx_const * eos->cs()) * (vx_const + eos->cs() - (1 - eos->cs() * (vx_const + eos->cs()) / (1 + vx_const * eos->cs()) ) * deltaVxr ) ) );

//     cout << br << "     " << bl << endl;

     if (bl > 0.) {
         cout << 1 / (1 - vb_const * csb) * (vb_const - csb + (1 + csb * (vb_const - csb) / (1 - vb_const * csb) ) * vb_delta ) << "        " << 1 / (1 - vz_const * eos->cs()) * (vz_const - eos->cs() + (eos->cs() * (vz_const - eos->cs()) / (1 - vz_const * eos->cs()) + 1 ) * deltaVzl ) << endl;
         sleep(5);
     }

    if (bl == 1 / (1 - vb_const * csb) * (vb_const - csb + (1 + csb * (vb_const - csb) / (1 - vb_const * csb) ) * vb_delta ) ){

         v0_L = vb_const;
         delta_vL = vb_delta;
         cs = csb;
         bl_const = (vb_const - csb) / (1 - vb_const * csb);

    }
    if (bl == 1 / (1 - vx_const * eos->cs()) * (vx_const - eos->cs() + (eos->cs() * (vx_const - eos->cs()) / (1 - vx_const * eos->cs()) + 1 ) * deltaVxl ) ){

        v0_L = vx_const;
        delta_vL = deltaVxl;
        cs = eos->cs(); // tohle je spatne definovane, pro jedno to muze byt cs, pro druhe csb - pokud plati csb = sqrt(eos->cs2()), tak je to jedno!!!!!!!!!!!!!!!!!!!!!
        bl_const = (vx_const - eos->cs()) / (1 - vx_const * eos->cs());

    }
    if (br == 1 / (1 + vb_const * csb) * (vb_const + csb - (1 - csb * (vb_const + csb) / (1 + vb_const * csb) ) * vb_delta ) ){

        v0_R = vb_const;
        delta_vR = vb_delta;
        cs = csb;
        br_const = (vb_const + csb) / (1 + vb_const * csb);

    }
    if (br == 1 / (1 + vx_const * eos->cs()) * (vx_const + eos->cs() - (1 - eos->cs() * (vx_const + eos->cs()) / (1 + vx_const * eos->cs()) ) * deltaVxr ) ){

        v0_R = vx_const;
        delta_vR = deltaVxr;
        cs = eos->cs();
        br_const = (vx_const + eos->cs()) / (1 + vx_const * eos->cs());

    }

     P = (v0_L - cs) / ( (1-v0_L * cs) * (1 + v0_R * cs) ) * ( 1 - cs * (v0_R + cs) / (1 + v0_R * cs) );
     Q = (v0_R + cs) / ( (1-v0_L * cs) * (1 + v0_R * cs) ) * ( 1 + cs * (v0_L - cs) / (1 - v0_L * cs) );
     B = 1 / (1 + v0_R * cs) * ( 1 - cs * (v0_R + cs) / (1 + v0_R * cs) );
     C = 1 / (1 - v0_L * cs) * ( 1 + cs * (v0_L - cs) / (1 - v0_L * cs) );

    flux_T_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U4l_flux + U4l_delta_v - U4r_flux - U4r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U4l - U4r);
    flux_X_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U1l_flux + U1l_delta_v - U1r_flux - U1r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U1l - U1r);
    flux_Y_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U2l_flux + U2l_delta_v - U2r_flux - U2r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U2l - U2r);
    flux_Z_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U3l_flux + U3l_delta_v - U3r_flux - U3r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U3l - U3r);

    flux_T_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Ftr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U4r_flux * vx_const;
    flux_X_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fxr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * (U1r_flux * vx_const + p_const);
    flux_Y_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fyr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U2r_flux * vx_const;
    flux_Z_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fzr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U3r_flux * vx_const;

    flux_T_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Ftl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U4l_flux * vx_const;
    flux_X_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fxl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * (U1l_flux * vx_const + p_const);
    flux_Y_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fyl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U2l_flux * vx_const;
    flux_Z_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fzl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U3l_flux * vx_const;

     if (br == 0) {
         flux_T_UU = 0. ;
         flux_X_UU = 0. ;
         flux_Y_UU = 0. ;
         flux_Z_UU = 0. ;

         flux_T_Fl = 0. ;
         flux_X_Fl = 0. ;
         flux_Y_Fl = 0. ;
         flux_Z_Fl = 0. ;
     }

     if (bl == 0) {
         flux_T_UU = 0. ;
         flux_X_UU = 0. ;
         flux_Y_UU = 0. ;
         flux_Z_UU = 0. ;

         flux_T_Fr = 0. ;
         flux_X_Fr = 0. ;
         flux_Y_Fr = 0. ;
         flux_Z_Fr = 0. ;
     }


  dx = f->getDx();

  // bl or br in the case of boundary with vacuum
  if (el == 0.) bl = -1.;
  if (er == 0.) br = 1.;
 }
 if (direction == Y_) {
  Ftl = U4l_const_v * vyl0 + U4l_delta_v * deltaVyl;
  Fxl = U1l_const_v * vyl0 + U1l_delta_v * deltaVyl;
  Fyl = U2l_const_v * vyl0 + U2l_delta_v * deltaVyl + pl;
  Fzl = U3l_const_v * vyl0 + U3l_delta_v * deltaVyl;
 // Fzl = U3l * vyl;
  Fbl = Ubl * vyl;
  Fql = Uql * vyl;
  Fsl = Usl * vyl;

  Ftr = U4r_const_v * vyr0 + U4r_delta_v * deltaVyr;
  Fxr = U1r_const_v * vyr0 + U1r_delta_v * deltaVyr;
  Fyr = U2r_const_v * vyr0 + U2r_delta_v * deltaVyr + pr;
  Fzr = U3r_const_v * vyr0 + U3r_delta_v * deltaVyr;
  //Fzr = U3r * vyr;
  Fbr = Ubr * vyr;
  Fqr = Uqr * vyr;
  Fsr = Usr * vyr;

  // for the case of constant c_s only
//  csb = sqrt(eos->cs2() +
//             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
//                 pow(vyl - vyr, 2));
//  vb = (sqrt(El) * vyl + sqrt(Er) * vyr) / (sqrt(El) + sqrt(Er));
//  bl = min(0., min((vb - csb) / (1 - vb * csb),
//                   (vyl - eos->cs()) / (1 - vyl * eos->cs())));
//  br = max(0., max((vb + csb) / (1 + vb * csb),
//                   (vyr + eos->cs()) / (1 + vyr * eos->cs())));
     csb = sqrt(eos->cs2());
     vb_const = vy_const;
     vb_delta = 0.5 * (deltaVyl + deltaVyr);
     vb = vb_const + vb_delta;
     bl = min(0., min(1 / (1 - vb_const * csb) * (vb_const - csb + (1 + csb * (vb_const - csb) / (1 - vb_const * csb) ) * vb_delta ),
                     1 / (1 - vy_const * eos->cs()) * (vy_const - eos->cs() + (eos->cs() * (vy_const - eos->cs()) / (1 - vy_const * eos->cs()) + 1 ) * deltaVyl ) ) );
     br = max(0., max(1 / (1 + vb_const * csb) * (vb_const + csb - (1 - csb * (vb_const + csb) / (1 + vb_const * csb) ) * vb_delta ),
                     1 / (1 + vy_const * eos->cs()) * (vy_const + eos->cs() - (1 - eos->cs() * (vy_const + eos->cs()) / (1 + vy_const * eos->cs()) ) * deltaVyr ) ) );

     if (bl > 0.) {
         cout << 1 / (1 - vb_const * csb) * (vb_const - csb + (1 + csb * (vb_const - csb) / (1 - vb_const * csb) ) * vb_delta ) << "        " << 1 / (1 - vz_const * eos->cs()) * (vz_const - eos->cs() + (eos->cs() * (vz_const - eos->cs()) / (1 - vz_const * eos->cs()) + 1 ) * deltaVzl ) << endl;
         sleep(5);
     }

     if (bl == 1 / (1 - vb_const * csb) * (vb_const - csb + (1 + csb * (vb_const - csb) / (1 - vb_const * csb) ) * vb_delta ) ){

         v0_L = vb_const;
         delta_vL = vb_delta;
         cs = csb;
         bl_const = (vb_const - csb) / (1 - vb_const * csb);

     }
     if (bl == 1 / (1 - vy_const * eos->cs()) * (vy_const - eos->cs() + (eos->cs() * (vy_const - eos->cs()) / (1 - vy_const * eos->cs()) + 1 ) * deltaVyl ) ){

         v0_L = vy_const;
         delta_vL = deltaVyl;
         cs = eos->cs();
         bl_const = (vy_const - eos->cs()) / (1 - vy_const * eos->cs());

     }
     if (br == 1 / (1 + vb_const * csb) * (vb_const + csb - (1 - csb * (vb_const + csb) / (1 + vb_const * csb) ) * vb_delta ) ){

         v0_R = vb_const;
         delta_vR = vb_delta;
         cs = csb;
         br_const = (vb_const + csb) / (1 + vb_const * csb);

     }
     if (br == 1 / (1 + vy_const * eos->cs()) * (vy_const + eos->cs() - (1 - eos->cs() * (vy_const + eos->cs()) / (1 + vy_const * eos->cs()) ) * deltaVyr ) ){

         v0_R = vy_const;
         delta_vR = deltaVyr;
         cs = eos->cs();
         br_const = (vy_const + eos->cs()) / (1 + vy_const * eos->cs());

     }

     P = (v0_L - cs) / ( (1-v0_L * cs) * (1 + v0_R * cs) ) * ( 1 - cs * (v0_R + cs) / (1 + v0_R * cs) );
     Q = (v0_R + cs) / ( (1-v0_L * cs) * (1 + v0_R * cs) ) * ( 1 + cs * (v0_L - cs) / (1 - v0_L * cs) );
     B = 1 / (1 + v0_R * cs) * ( 1 - cs * (v0_R + cs) / (1 + v0_R * cs) );
     C = 1 / (1 - v0_L * cs) * ( 1 + cs * (v0_L - cs) / (1 - v0_L * cs) );

     flux_T_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U4l_flux + U4l_delta_v - U4r_flux - U4r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U4l - U4r);
     flux_X_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U1l_flux + U1l_delta_v - U1r_flux - U1r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U1l - U1r);
     flux_Y_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U2l_flux + U2l_delta_v - U2r_flux - U2r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U2l - U2r);
     flux_Z_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U3l_flux + U3l_delta_v - U3r_flux - U3r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U3l - U3r);

     flux_T_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Ftr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U4r_flux * vy_const;
     flux_X_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fxr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U1r_flux * vy_const;
     flux_Y_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fyr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * (U2r_flux * vy_const + p_const);
     flux_Z_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fzr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U3r_flux * vy_const;

     flux_T_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Ftl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U4l_flux * vy_const;
     flux_X_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fxl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U1l_flux * vy_const;
     flux_Y_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fyl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * (U2l_flux * vy_const + p_const);
     flux_Z_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fzl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U3l_flux * vy_const;

     if (br == 0) {
         flux_T_UU = 0. ;
         flux_X_UU = 0. ;
         flux_Y_UU = 0. ;
         flux_Z_UU = 0. ;

         flux_T_Fl = 0. ;
         flux_X_Fl = 0. ;
         flux_Y_Fl = 0. ;
         flux_Z_Fl = 0. ;
     }

     if (bl == 0) {
         flux_T_UU = 0. ;
         flux_X_UU = 0. ;
         flux_Y_UU = 0. ;
         flux_Z_UU = 0. ;

         flux_T_Fr = 0. ;
         flux_X_Fr = 0. ;
         flux_Y_Fr = 0. ;
         flux_Z_Fr = 0. ;
     }


  dx = f->getDy();

  // bl or br in the case of boundary with vacuum
  if (el == 0.) bl = -1.;
  if (er == 0.) br = 1.;
 }
 if (direction == Z_) {
  double tau1 = tauFactor;

  Ftl = U4l_const_v * vzl0 / tau1 + U4l_delta_v * deltaVzl / tau1;
  Fxl = U1l_const_v * vzl0 / tau1 + U1l_delta_v * deltaVzl / tau1;
  Fyl = U2l_const_v * vzl0 / tau1 + U2l_delta_v * deltaVzl / tau1;
  Fzl = U3l_const_v * vzl0 / tau1 + U3l_delta_v * deltaVzl / tau1 + pl / tau1;
  Fbl = Ubl * vzl / tau1;
  Fql = Uql * vzl / tau1;
  Fsl = Usl * vzl / tau1;

  Ftr = U4r_const_v * vzr0 / tau1 + U4r_delta_v * deltaVzr / tau1;
  Fxr = U1r_const_v * vzr0 / tau1 + U1r_delta_v * deltaVzr / tau1;
  Fyr = U2r_const_v * vzr0 / tau1 + U2r_delta_v * deltaVzr / tau1;
  Fzr = U3r_const_v * vzr0 / tau1 + U3r_delta_v * deltaVzr / tau1 + pr / tau1;
  Fbr = Ubr * vzr / tau1;
  Fqr = Uqr * vzr / tau1;
  Fsr = Usr * vzr / tau1;

  // for the case of constant c_s only
  // factor 1/tau accounts for eta-coordinate

  // different estimate
//  csb = sqrt(eos->cs2() +
//             0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
//                 pow(vzl - vzr, 2));
//  vb = (sqrt(El) * vzl + sqrt(Er) * vzr) / (sqrt(El) + sqrt(Er));
//  bl = 1. / tau * min(0., min((vb - csb) / (1 - vb * csb),
//                              (vzl - eos->cs()) / (1 - vzl * eos->cs())));
//  br = 1. / tau * max(0., max((vb + csb) / (1 + vb * csb),
//                              (vzr + eos->cs()) / (1 + vzr * eos->cs())));

     csb = sqrt(eos->cs2());
     vb_const = vz_const;
     vb_delta = 0.5 * (deltaVzl + deltaVzr);
     vb = vb_const + vb_delta;
     bl = 1. / tau * min(0., min(1 / (1 - vb_const * csb) * (vb_const - csb + (1 + csb * (vb_const - csb) / (1 - vb_const * csb) ) * vb_delta ),
                     1 / (1 - vz_const * eos->cs()) * (vz_const - eos->cs() + (eos->cs() * (vz_const - eos->cs()) / (1 - vz_const * eos->cs()) + 1 ) * deltaVzl ) ) );
     br = 1. / tau * max(0., max(1 / (1 + vb_const * csb) * (vb_const + csb - (1 - csb * (vb_const + csb) / (1 + vb_const * csb) ) * vb_delta ),
                     1 / (1 + vz_const * eos->cs()) * (vz_const + eos->cs() - (1 - eos->cs() * (vz_const + eos->cs()) / (1 + vz_const * eos->cs()) ) * deltaVzr ) ) );

    if (bl == 1 / (1 - vb_const * csb) * (vb_const - csb + (1 + csb * (vb_const - csb) / (1 - vb_const * csb) ) * vb_delta ) ){

        v0_L = vb_const;
        delta_vL = vb_delta;
        cs = csb;
        bl_const = (vb_const - csb) / (1 - vb_const * csb);

    }
    if (bl == 1 / (1 - vz_const * eos->cs()) * (vz_const - eos->cs() + (eos->cs() * (vz_const - eos->cs()) / (1 - vz_const * eos->cs()) + 1 ) * deltaVzl ) ){

        v0_L = vz_const;
        delta_vL = deltaVzl;
        cs = eos->cs();
        bl_const = (vz_const - eos->cs()) / (1 - vz_const * eos->cs());

    }
    if (br == 1 / (1 + vb_const * csb) * (vb_const + csb - (1 - csb * (vb_const + csb) / (1 + vb_const * csb) ) * vb_delta ) ){

        v0_R = vb_const;
        delta_vR = vb_delta;
        cs = csb;
        br_const = (vb_const + csb) / (1 + vb_const * csb);


    }
    if (br == 1 / (1 + vz_const * eos->cs()) * (vz_const + eos->cs() - (1 - eos->cs() * (vz_const + eos->cs()) / (1 + vz_const * eos->cs()) ) * deltaVzr ) ){

        v0_R = vz_const;
        delta_vR = deltaVzr;
        cs = eos->cs();
        br_const = (vz_const + eos->cs()) / (1 + vz_const * eos->cs());

    }

     P = (v0_L - cs) / ( (1-v0_L * cs) * (1 + v0_R * cs) ) * ( 1 - cs * (v0_R + cs) / (1 + v0_R * cs) );
     Q = (v0_R + cs) / ( (1-v0_L * cs) * (1 + v0_R * cs) ) * ( 1 + cs * (v0_L - cs) / (1 - v0_L * cs) );
     B = 1 / (1 + v0_R * cs) * ( 1 - cs * (v0_R + cs) / (1 + v0_R * cs) );
     C = 1 / (1 - v0_L * cs) * ( 1 + cs * (v0_L - cs) / (1 - v0_L * cs) );

    flux_T_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U4l_flux + U4l_delta_v - U4r_flux - U4r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U4l - U4r);
    flux_X_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U1l_flux + U1l_delta_v - U1r_flux - U1r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U1l - U1r);
    flux_Y_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U2l_flux + U2l_delta_v - U2r_flux - U2r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U2l - U2r);
    flux_Z_UU = ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U3l_flux + U3l_delta_v - U3r_flux - U3r_delta_v) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U3l - U3r);

    flux_T_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Ftr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U4r_flux * vz_const;
    flux_X_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fxr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U1r_flux * vz_const;
    flux_Y_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fyr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * U2r_flux * vz_const;
    flux_Z_Fr = (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) ) * Fzr + ( (v0_L - cs) / ( (1 - v0_L * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) + C * delta_vL / (br_const - bl_const) ) * (U3r_flux * vz_const + p_const);

    flux_T_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Ftl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U4l_flux * vz_const;
    flux_X_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fxl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U1l_flux * vz_const;
    flux_Y_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fyl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * U2l_flux * vz_const;
    flux_Z_Fl = (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) ) * Fzl + ( (v0_R + cs) / ( (1 + v0_R * cs) * (br_const - bl_const) * (br_const - bl_const) ) * (B * delta_vR + C * delta_vL) - B * delta_vR / (br_const - bl_const) ) * (U3l_flux * vz_const + p_const);

    if (br == 0.) {
        flux_T_UU = 0. ;
        flux_X_UU = 0. ;
        flux_Y_UU = 0. ;
        flux_Z_UU = 0. ;

        flux_T_Fl = 0. ;
        flux_X_Fl = 0. ;
        flux_Y_Fl = 0. ;
        flux_Z_Fl = 0. ;
    }

    if (bl == 0.) {
        flux_T_UU = 0. ;
        flux_X_UU = 0. ;
        flux_Y_UU = 0. ;
        flux_Z_UU = 0. ;

        flux_T_Fr = 0. ;
        flux_X_Fr = 0. ;
        flux_Y_Fr = 0. ;
        flux_Z_Fr = 0. ;
    }

  dx = f->getDz();

  // bl or br in the case of boundary with vacuum
  if (el == 0.) bl = -1. / tau;
  if (er == 0.) br = 1. / tau;
 }

 if (bl == 0. && br == 0.) return;

//  finally, HLLE formula for the fluxes
// flux[T_] = tauFactor * dta / dx *
//            (-bl * br * (U4l - U4r) + br * Ftl - bl * Ftr) / (-bl + br);
// flux[X_] = tauFactor * dta / dx *
//            (-bl * br * (U1l - U1r) + br * Fxl - bl * Fxr) / (-bl + br);
// flux[Y_] = tauFactor * dta / dx *
//            (-bl * br * (U2l - U2r) + br * Fyl - bl * Fyr) / (-bl + br);
// flux[Z_] = tauFactor * dta / dx *
//            (-bl * br * (U3l - U3r) + br * Fzl - bl * Fzr) / (-bl + br);
// flux[NB_] = tauFactor * dta / dx *
//             (-bl * br * (Ubl - Ubr) + br * Fbl - bl * Fbr) / (-bl + br);
// flux[NQ_] = tauFactor * dta / dx *
//             (-bl * br * (Uql - Uqr) + br * Fql - bl * Fqr) / (-bl + br);
// flux[NS_] = tauFactor * dta / dx *
//             (-bl * br * (Usl - Usr) + br * Fsl - bl * Fsr) / (-bl + br);

    flux[T_] = tauFactor * dta / dx *
               (flux_T_UU + flux_T_Fl - flux_T_Fr);
    flux[X_] = tauFactor * dta / dx *
               (flux_X_UU + flux_X_Fl - flux_X_Fr);
    flux[Y_] = tauFactor * dta / dx *
               (flux_Y_UU + flux_Y_Fl - flux_Y_Fr);
    flux[Z_] = tauFactor * dta / dx *
               (flux_Z_UU + flux_Z_Fl - flux_Z_Fr);
    flux[NB_] = tauFactor * dta / dx *
                (-bl * br * (Ubl - Ubr) + br * Fbl - bl * Fbr) / (-bl + br);
    flux[NQ_] = tauFactor * dta / dx *
                (-bl * br * (Uql - Uqr) + br * Fql - bl * Fqr) / (-bl + br);
    flux[NS_] = tauFactor * dta / dx *
                (-bl * br * (Usl - Usr) + br * Fsl - bl * Fsr) / (-bl + br);
    
        if(ix<=40){
            cout.precision(20);
        }
        else{
            cout.precision(20);
        }
//    if(direction==1){
//        if(ix==20 || ix ==19 || ix ==21 || ix ==60 || ix ==61 || ix ==59){
            //    if(ix==40 || ix==39 || ix==41){
//            cout << mode << setw(12) << direction << endl << endl;
//            cout << ix << endl;
            //            cout << flux[T_] << setw(12) << flux[X_] << setw(12) << flux[Y_] << setw(12) << flux[Z_] << setw(12) << "Fxr=" << setw(3) << Fxr << setw(12) << "Fxl=" << setw(3) << Fxl << endl << endl;
            ////            cout << "flux_X_UU=" << setw(3) << flux_X_UU << setw(10) << "flux_X_Fl=" << setw(3) << flux_X_Fl << setw(10) << "flux_X_Fr=" << setw(3) << flux_X_Fr << endl;
            //            cout << "br=" << setw(3) << br << setw(12) << "br_const=" << setw(3) << br_const << setw(12) << "bl=" << setw(3) << bl << setw(12) << "bl_const=" << setw(3) << bl_const << setw(12) << "Fxl=" << setw(3) << Fxl << setw(12) << "Fxr=" << setw(3) << Fxr << setw(12) << "U1r_flux=" << setw(3) << U1r_flux << setw(12) << "U1l_flux=" << setw(3) << U1l_flux << setw(12) << endl;
            //            cout << "U1l_const_v=" << setw(3) << U1l_const_v << setw(12) << "U1r_const_v=" << setw(3) << U1r_const_v << setw(12) << "U1l_delta_v=" << setw(3) << U1l_delta_v << setw(12) << "U1r_delta_v=" << setw(3) << U1r_delta_v << setw(12) << "pl=" << setw(3) << pl << setw(12) << "pr=" << setw(3) << pr << endl;
            ////            cout << "U1l_const_v=" << setw(3) << U1l_const_v << setw(3)<< "U1l_delta_v=" << setw(3) << U1l_delta_v << setw(3)<< "U1r_const_v=" << setw(3) << U1r_const_v << setw(3)<< "U1r_delta_v=" << setw(3) << U1r_delta_v << endl;
            //        cout << "el=" << setw(3) << el << setw(10) << "er=" << setw(3) << er << endl;
//            cout << "U4l_flux=" << "    " << U4l_flux << "      " << "U4r_flux=" << "     " << U4r_flux << "    " << "U4l_const_v=" << "     " << U4l_const_v << "        " << "U4r_const_v=" << "      " << U4r_const_v << endl;
//            cout << "U4l_delta_v=" << "     " << U4l_delta_v << "       " << "U4r_delta_v=" << "    " << U4l_delta_v << endl;
            //            cout << "U1l=" << "  " << U1l << "      " << "U1r=" << "    " << U1r << endl;
            //            cout << "(v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U1l - U1r)=" << "     " << (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U1l - U1r) << endl;
//                        cout << "U1l-U1r=" << "     " << U4l - U4r << endl;
            //            cout << U1l_const_v + U1l_delta_v - (U1r_const_v + U1r_delta_v) << endl;
            //            cout << "flux_X_UU=" << "  " << flux_X_UU << endl;
//            cout << U4l_flux + U4l_delta_v - U4r_flux - U4r_delta_v << endl;
//            cout << "vxl0=" << "    " << vxl0 << "      " << "vxr0=" << "   " << vxr0 << endl;
//            cout << "Fxl=" << "     " << Fxl << "     " << "Fxr=" << "  " << Fxr << "       " << "U1l_flux=" << "   " << U1l_flux << "      " << "U1r_flux=" << "   " << U1r_flux << endl;
//            cout << U4l_const_v * vxl0 << endl;
//            cout << "flux_X_UU=" << "   " << flux_X_UU << endl;
//            cout << ( (P * delta_vR - Q * delta_vL) / (br_const - bl_const) - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) * (B * delta_vR + C * delta_vL) / ( (br_const - bl_const) * (br_const - bl_const) ) ) * (U1l_flux + U1l_delta_v - U1r_flux - U1r_delta_v) << endl;
//            cout << - (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) / (br_const - bl_const) * (U1l - U1r) << endl;
//            cout << (v0_L - cs) * (v0_R + cs) / ( (1 - v0_L * cs) * (1 + v0_R * cs) ) << endl;
//            cout << U1l - U1r << endl;
//            cout << U1l_const_v << "    " << U1l_delta_v << endl;
//            cout << U1r_const_v << "    " << U1r_delta_v << endl;
//            cout << gamma_const * (gammar + gammar_delta) * (er0 + pr0) * vxr0  << "    " << gamma_const * gamma_const * ( (er0 + pr0) * deltaVxr + (deltaEr + deltaPr) * vxr0 ) << endl;
//            cout << gamma_const * gamma_const * ( (er0 + pr0) * deltaVxr + gamma_const * gamma_const *(deltaEr + deltaPr) * vxr0 ) endl;

            
//        }
//    }

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
 #else
 double _Q[7], e, p, nb, nq, ns, vx, vy, vz;
 for (int i = 0; i < 7; i++) _Q[i] = Q[i] / tau1;  // no tau factor in  _Q
 transformPV(eos, _Q, e, p, nb, nq, ns, vx, vy, vz);
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
 double uuu[4];
 double k[7];

 #ifdef CARTESIAN
  // there is no geometric viscous source term in Cartesian frame
 #else
 Cell *c = f->getCell(ix, iy, iz);

 c->getPrimVarHCenter(eos, tau - dt / 2., e, p, nb, nq, ns, vx, vy, vz);  // TODO Cartesian
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
void Hydro::NSquant(int ix, int iy, int iz, double pi[4][4], double &Pi,
                    double dmu[4][4], double dmu_const[4][4], double &du, double e_const, double p_const, double vx_const, double vy_const, double vz_const, double eH_const, double pH_const, double vxH_const, double vyH_const, double vzH_const, double ePrev_const, double pPrev_const, double vxPrev_const, double vyPrev_const, double vzPrev_const) {
    const double VMIN = 1e-2;
    const double UDIFF = 3.0;
    double e0, e1, p, nb, nq, ns, vx1, vy1, vz1, vx0, vy0, vz0, vxH, vyH, vzH;
    double ut0, ux0, uy0, uz0, ut1, ux1, uy1, uz1;
    //	double dmu [4][4] ; // \partial_\mu u^\nu matrix
    // coordinates: 0=tau, 1=x, 2=y, 3=eta
    double Z[4][4][4][4];  // Z[mu][nu][lambda][rho]
    double Z_dmu[4][4]; // linearized product of Z and dmu
    double uuu[4];         // the 4-velocity
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
    c->getPrimVarPrev(eos, 1.0, e0, p, nb, nq, ns, vx0, vy0, vz0);
    c->getPrimVar(eos, 1.0, e1, p, nb, nq, ns, vx1, vy1, vz1);
    c->getPrimVarHCenter(eos, 1.0, e1, p, nb, nq, ns, vxH, vyH, vzH);
    double tauPlusHalf = 1.0;
#else
    c->getPrimVarPrev(eos, tau - dt, e0, p, nb, nq, ns, vx0, vy0, vz0);
    c->getPrimVar(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
    c->getPrimVarHCenter(eos, tau - 0.5 * dt, e1, p, nb, nq, ns, vxH, vyH, vzH);
    double tauPlusHalf = tau + 0.5 * dt;
#endif
    
    //double dmu_const[4][4];
    double deltaVx0, deltaVy0, deltaVz0, deltaVx1, deltaVy1, deltaVz1, deltaVxH, deltaVyH, deltaVzH;
    
    deltaVx0 = vx0 - vxPrev_const;
    deltaVy0 = vy0 - vyPrev_const;
    deltaVz0 = vz0 - vzPrev_const;

    deltaVx1 = vx1 - vx_const;
    deltaVy1 = vy1 - vy_const;
    deltaVz1 = vz1 - vz_const;

    deltaVxH = vxH - vxH_const;
    deltaVyH = vyH - vyH_const;
    deltaVzH = vzH - vzH_const;
    
    //############## get transport coefficients
    double T, mub, muq, mus;
    double etaS, zetaS;
    double s = eos->s(e1, nb, nq, ns);  // entropy density in the current cell
    eos->eos(e1, nb, nq, ns, T, mub, muq, mus, p);
    trcoeff->getEta(e1, nb, T, etaS, zetaS);
    //##############
    // if(e1<0.00004) s=0. ; // negative pressure due to pi^zz for small e
//    ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
//    ux0 = ut0 * vx0;// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    uy0 = ut0 * vy0;
//    uz0 = ut0 * vz0;
//    ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
//    ux1 = ut1 * vx1;
//    uy1 = ut1 * vy1;
//    uz1 = ut1 * vz1;
//    uuu[0] = 1.0 / sqrt(1.0 - vxH * vxH - vyH * vyH - vzH * vzH);
//    uuu[1] = uuu[0] * vxH;
//    uuu[2] = uuu[0] * vyH;
//    uuu[3] = uuu[0] * vzH;
    
    double gamma_u0_const;
    double gamma_u0_delta;

    gamma_u0_const = 1. + 0.5 * (vxPrev_const * vxPrev_const + vyPrev_const * vyPrev_const + vzPrev_const * vzPrev_const);
    gamma_u0_delta = vxPrev_const * deltaVx0 + vyPrev_const * deltaVy0 + vzPrev_const * deltaVz0;
    ut0 = gamma_u0_const + gamma_u0_delta;
    ux0 = ut0 * vxPrev_const + deltaVx0 * gamma_u0_const;
    uy0 = ut0 * vyPrev_const + deltaVy0 * gamma_u0_const;
    uz0 = ut0 * vzPrev_const + deltaVz0 * gamma_u0_const;
    
    double gamma_u_const;
    double gamma_u_delta;
    
    gamma_u_const = 1. + 0.5 * (vx_const * vx_const + vy_const * vy_const + vz_const * vz_const);
    gamma_u_delta = vx_const * deltaVx1 + vy_const * deltaVy1 + vz_const * deltaVz1;
    ut1 = gamma_u_const + gamma_u_delta;
    ux1 = ut1 * vx_const + deltaVx1 * gamma_u_const;
    uy1 = ut1 * vy_const + deltaVy1 * gamma_u_const;
    uz1 = ut1 * vz_const + deltaVz1 * gamma_u_const;
    
    double gamma_const = 1. + 0.5 * (vxH_const * vxH_const + vyH_const * vyH_const + vzH_const * vzH_const);
    double gamma_delta = vxH_const * deltaVxH + vyH_const * deltaVyH + vzH_const * deltaVzH;
    uuu[0] = gamma_const + gamma_delta;
    uuu[1] = uuu[0] * vxH_const + deltaVxH * gamma_const;
    uuu[2] = uuu[0] * vyH_const + deltaVyH * gamma_const;
    uuu[3] = uuu[0] * vzH_const + deltaVzH * gamma_const;
    
    dmu[0][0] = (ut1 - ut0) / dt;
    dmu[0][1] = (ux1 - ux0) / dt;
    dmu[0][2] = (uy1 - uy0) / dt;
    dmu[0][3] = (uz1 - uz0) / dt;
    
    dmu_const[0][0] = (gamma_u_const - gamma_u0_const) / dt;
    dmu_const[0][1] = (gamma_u_const * vx_const- gamma_u0_const * vxPrev_const) / dt;
    dmu_const[0][2] = (gamma_u_const * vy_const- gamma_u0_const * vyPrev_const) / dt;
    dmu_const[0][3] = (gamma_u_const * vz_const- gamma_u0_const * vzPrev_const) / dt;
    
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
    if (e1 <= 0. || e0 <= 0.) {  // matter-vacuum
        dmu[0][0] = dmu[0][1] = dmu[0][2] = dmu[0][3] = 0.;
    }
    // d_x u^\mu
    f->getCell(ix + 1, iy, iz)
    ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);//0 i 1 maji stejne v_const
    f->getCell(ix - 1, iy, iz)
    ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
    
    deltaVx0 = vx0 - vxPrev_const;
    deltaVy0 = vy0 - vyPrev_const;
    deltaVz0 = vz0 - vzPrev_const;
    deltaVx1 = vx1 - vx_const;
    deltaVy1 = vy1 - vy_const;
    deltaVz1 = vz1 - vz_const;
    
    if (e1 > 0. && e0 > 0.) {
        
        gamma_u0_const = 1. + 0.5 * (vxPrev_const * vxPrev_const + vyPrev_const * vyPrev_const + vzPrev_const * vzPrev_const);
        gamma_u0_delta = vxPrev_const * deltaVx0 + vyPrev_const * deltaVy0 + vzPrev_const * deltaVz0;
        ut0 = gamma_u0_const + gamma_u0_delta;
        ux0 = ut0 * vxPrev_const + deltaVx0 * gamma_u0_const;
        uy0 = ut0 * vyPrev_const + deltaVy0 * gamma_u0_const;
        uz0 = ut0 * vzPrev_const + deltaVz0 * gamma_u0_const;
        
        gamma_u_const = 1. + 0.5 * (vx_const * vx_const + vy_const * vy_const + vz_const * vz_const);
        gamma_u_delta = vx_const * deltaVx1 + vy_const * deltaVy1 + vz_const * deltaVz1;
        ut1 = gamma_u_const + gamma_u_delta;
        ux1 = ut1 * vx_const + deltaVx1 * gamma_u_const;
        uy1 = ut1 * vy_const + deltaVy1 * gamma_u_const;
        uz1 = ut1 * vz_const + deltaVz1 * gamma_u_const;
                
        dmu[1][0] = 0.5 * (ut1 - ut0) / dx;//why 0.5??????????
        dmu[1][1] = 0.5 * (ux1 - ux0) / dx;
        dmu[1][2] = 0.5 * (uy1 - uy0) / dx;
        dmu[1][3] = 0.5 * (uz1 - uz0) / dx;
        
        dmu_const[1][0] = 0.5 * (gamma_u_const - gamma_u0_const) / dx;
        dmu_const[1][1] = 0.5 * (gamma_u_const * vx_const- gamma_u0_const * vxPrev_const) / dx;
        dmu_const[1][2] = 0.5 * (gamma_u_const * vy_const- gamma_u0_const * vyPrev_const) / dx;
        dmu_const[1][3] = 0.5 * (gamma_u_const * vz_const- gamma_u0_const * vzPrev_const) / dx;
        
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
        }
        if (fabs(dmu[1][3]) > 1e+10)
            cout << "dmu[1][3]:  " << uz1 << "  " << uz0 << "  " << uuu[3] << endl;
        // d_y u^\mu
        f->getCell(ix, iy + 1, iz)
        ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
        f->getCell(ix, iy - 1, iz)
        ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
    
        deltaVx0 = vx0 - vxPrev_const;
        deltaVy0 = vy0 - vyPrev_const;
        deltaVz0 = vz0 - vzPrev_const;
        deltaVx1 = vx1 - vx_const;
        deltaVy1 = vy1 - vy_const;
        deltaVz1 = vz1 - vz_const;
    
        if (e1 > 0. && e0 > 0.) {
            
            gamma_u0_const = 1. + 0.5 * (vxPrev_const * vxPrev_const + vyPrev_const * vyPrev_const + vzPrev_const * vzPrev_const);
            gamma_u0_delta = vxPrev_const * deltaVx0 + vyPrev_const * deltaVy0 + vzPrev_const * deltaVz0;
            ut0 = gamma_u0_const + gamma_u0_delta;
            ux0 = ut0 * vxPrev_const + deltaVx0 * gamma_u0_const;
            uy0 = ut0 * vyPrev_const + deltaVy0 * gamma_u0_const;
            uz0 = ut0 * vzPrev_const + deltaVz0 * gamma_u0_const;
            
            gamma_u_const = 1. + 0.5 * (vx_const * vx_const + vy_const * vy_const + vz_const * vz_const);
            gamma_u_delta = vx_const * deltaVx1 + vy_const * deltaVy1 + vz_const * deltaVz1;
            ut1 = gamma_u_const + gamma_u_delta;
            ux1 = ut1 * vx_const + deltaVx1 * gamma_u_const;
            uy1 = ut1 * vy_const + deltaVy1 * gamma_u_const;
            uz1 = ut1 * vz_const + deltaVz1 * gamma_u_const;
                        
            dmu[2][0] = 0.5 * (ut1 - ut0) / dy;
            dmu[2][1] = 0.5 * (ux1 - ux0) / dy;
            dmu[2][2] = 0.5 * (uy1 - uy0) / dy;
            dmu[2][3] = 0.5 * (uz1 - uz0) / dy;
            
            dmu_const[2][0] = 0.5 * (gamma_u_const - gamma_u0_const) / dy;
            dmu_const[2][1] = 0.5 * (gamma_u_const * vx_const- gamma_u0_const * vxPrev_const) / dy;
            dmu_const[2][2] = 0.5 * (gamma_u_const * vy_const- gamma_u0_const * vyPrev_const) / dy;
            dmu_const[2][3] = 0.5 * (gamma_u_const * vz_const- gamma_u0_const * vzPrev_const) / dy;
            
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
        }
        // d_z u^\mu
        f->getCell(ix, iy, iz + 1)
        ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
        f->getCell(ix, iy, iz - 1)
        ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
    
        deltaVx0 = vx0 - vxPrev_const;
        deltaVy0 = vy0 - vyPrev_const;
        deltaVz0 = vz0 - vzPrev_const;
        deltaVx1 = vx1 - vx_const;
        deltaVy1 = vy1 - vy_const;
        deltaVz1 = vz1 - vz_const;
    
        if (e1 > 0. && e0 > 0.) {
            
            gamma_u0_const = 1. + 0.5 * (vxPrev_const * vxPrev_const + vyPrev_const * vyPrev_const + vzPrev_const * vzPrev_const);
            gamma_u0_delta = vxPrev_const * deltaVx0 + vyPrev_const * deltaVy0 + vzPrev_const * deltaVz0;
            ut0 = gamma_u0_const + gamma_u0_delta;
            ux0 = ut0 * vxPrev_const + deltaVx0 * gamma_u0_const;
            uy0 = ut0 * vyPrev_const + deltaVy0 * gamma_u0_const;
            uz0 = ut0 * vzPrev_const + deltaVz0 * gamma_u0_const;
            
            gamma_u_const = 1. + 0.5 * (vx_const * vx_const + vy_const * vy_const + vz_const * vz_const);
            gamma_u_delta = vx_const * deltaVx1 + vy_const * deltaVy1 + vz_const * deltaVz1;
            ut1 = gamma_u_const + gamma_u_delta;
            ux1 = ut1 * vx_const + deltaVx1 * gamma_u_const;
            uy1 = ut1 * vy_const + deltaVy1 * gamma_u_const;
            uz1 = ut1 * vz_const + deltaVz1 * gamma_u_const;
            
            dmu[3][0] = 0.5 * (ut1 - ut0) / dz / tauPlusHalf;
            dmu[3][1] = 0.5 * (ux1 - ux0) / dz / tauPlusHalf;
            dmu[3][2] = 0.5 * (uy1 - uy0) / dz / tauPlusHalf;
            dmu[3][3] = 0.5 * (uz1 - uz0) / dz / tauPlusHalf;
            
            dmu_const[3][0] = 0.5 * (gamma_u_const - gamma_u0_const) / dz / tauPlusHalf;
            dmu_const[3][1] = 0.5 * (gamma_u_const * vx_const- gamma_u0_const * vxPrev_const) / dz / tauPlusHalf;
            dmu_const[3][2] = 0.5 * (gamma_u_const * vy_const- gamma_u0_const * vyPrev_const) / dz / tauPlusHalf;
            dmu_const[3][3] = 0.5 * (gamma_u_const * vz_const- gamma_u0_const * vzPrev_const) / dz / tauPlusHalf;
            
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
        }
        // additional terms from Christoffel symbols :)
#ifndef CARTESIAN
        dmu[3][0] += uuu[3] / (tau - 0.5 * dt);
        dmu[3][3] += uuu[0] / (tau - 0.5 * dt);
#endif
    
    double Tr_dmu = dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3];
    double Tr_dmu_const = dmu_const[0][0] + dmu_const[1][1] + dmu_const[2][2] + dmu_const[3][3];
    
    double uuu_product[4][4];
    //diagonal
    uuu_product[0][0] = gamma_const * (uuu[0] + gamma_delta);
    uuu_product[1][1] = gamma_const * ( (uuu[0] + gamma_delta) * vxH_const * vxH_const + 2 * gamma_const * vxH_const * deltaVxH );
    uuu_product[2][2] = gamma_const * ( (uuu[0] + gamma_delta) * vyH_const * vyH_const + 2 * gamma_const * vyH_const * deltaVyH );
    uuu_product[3][3] = gamma_const * ( (uuu[0] + gamma_delta) * vzH_const * vzH_const + 2 * gamma_const * vzH_const * deltaVzH );
    //off-diagonal
    uuu_product[0][1] = gamma_const * ( (uuu[0] + gamma_delta) * vxH_const + gamma_const * deltaVxH );
    uuu_product[0][2] = gamma_const * ( (uuu[0] + gamma_delta) * vyH_const + gamma_const * deltaVyH );
    uuu_product[0][3] = gamma_const * ( (uuu[0] + gamma_delta) * vzH_const + gamma_const * deltaVzH );

    uuu_product[1][0] = uuu_product[0][1];
    uuu_product[1][2] = gamma_const * ( (uuu[0] + gamma_delta) * vxH_const * vyH_const + gamma_const * vxH_const * deltaVyH + gamma_const * vyH_const * deltaVxH );
    uuu_product[1][3] = gamma_const * ( (uuu[0] + gamma_delta) * vxH_const * vzH_const + gamma_const * vxH_const * deltaVzH + gamma_const * vzH_const * deltaVxH );

    uuu_product[2][0] = uuu_product[0][2];
    uuu_product[2][1] = uuu_product[1][2];
    uuu_product[2][3] = gamma_const * ( (uuu[0] + gamma_delta) * vyH_const * vzH_const + gamma_const * vyH_const * deltaVzH + gamma_const * vzH_const * deltaVyH );
    
    uuu_product[3][0] = uuu_product[0][3];
    uuu_product[3][1] = uuu_product[1][3];
    uuu_product[3][2] = uuu_product[2][3];
        
    
    double uuu_product_dmu[4][4]; //matrix of linearized product of u*u*dmu
    
    uuu_product_dmu[0][0] =  gamma_const * ( gamma_const * (dmu[0][0] + vxH_const * dmu[1][0] + vyH_const * dmu[2][0] + vzH_const * dmu[3][0]) +
                                2 * gamma_delta * dmu_const[0][0] + (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * dmu_const[1][0] +
                                (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * dmu_const[2][0] + (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * dmu_const[3][0] ) ;
    uuu_product_dmu[0][1] =  gamma_const * ( gamma_const * (dmu[0][1] + vxH_const * dmu[1][1] + vyH_const * dmu[2][1] + vzH_const * dmu[3][1]) +
                                2 * gamma_delta * dmu_const[0][1] + (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * dmu_const[1][1] +
                                (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * dmu_const[2][1] + (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * dmu_const[3][1] ) ;
    uuu_product_dmu[0][2] =  gamma_const * ( gamma_const * (dmu[0][2] + vxH_const * dmu[1][2] + vyH_const * dmu[2][2] + vzH_const * dmu[3][2]) +
                                2 * gamma_delta * dmu_const[0][2] + (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * dmu_const[1][2] +
                                (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * dmu_const[2][2] + (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * dmu_const[3][2] ) ;
    uuu_product_dmu[0][3] =  gamma_const * ( gamma_const * (dmu[0][3] + vxH_const * dmu[1][3] + vyH_const * dmu[2][3] + vzH_const * dmu[3][3]) +
                                2 * gamma_delta * dmu_const[0][3] + (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * dmu_const[1][3] +
                                (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * dmu_const[2][3] + (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * dmu_const[3][3] ) ;
    uuu_product_dmu[1][0] =  gamma_const * ( gamma_const * (vxH_const * dmu[0][0] + vxH_const * vxH_const * dmu[1][0] + vxH_const * vyH_const * dmu[2][0] + vxH_const * vzH_const * dmu[3][0]) +
                                (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * dmu_const[0][0] +
                                2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVxH) * dmu_const[1][0] +
                                (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH)) * dmu_const[2][0] +
                                (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH)) * dmu_const[3][0] ) ;
    uuu_product_dmu[1][1] =  gamma_const * ( gamma_const * (vxH_const * dmu[0][1] + vxH_const * vxH_const * dmu[1][1] + vxH_const * vyH_const * dmu[2][1] + vxH_const * vzH_const * dmu[3][1]) +
                                (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * dmu_const[0][1] +
                                2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVxH) * dmu_const[1][1] +
                                (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH)) * dmu_const[2][1] +
                                (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH)) * dmu_const[3][1] ) ;
    uuu_product_dmu[1][2] =  gamma_const * ( gamma_const * (vxH_const * dmu[0][2] + vxH_const * vxH_const * dmu[1][2] + vxH_const * vyH_const * dmu[2][2] + vxH_const * vzH_const * dmu[3][2]) +
                                (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * dmu_const[0][2] +
                                2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVxH) * dmu_const[1][2] +
                                (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH)) * dmu_const[2][2] +
                                (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH)) * dmu_const[3][2] ) ;
    uuu_product_dmu[1][3] =  gamma_const * ( gamma_const * (vxH_const * dmu[0][3] + vxH_const * vxH_const * dmu[1][3] + vxH_const * vyH_const * dmu[2][3] + vxH_const * vzH_const * dmu[3][3]) +
                                (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * dmu_const[0][3] +
                                2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVxH) * dmu_const[1][3] +
                                (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH)) * dmu_const[2][3] +
                                (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH)) * dmu_const[3][3] ) ;
    uuu_product_dmu[2][0] =  gamma_const * ( gamma_const * (vyH_const * dmu[0][0] + vxH_const * vyH_const * dmu[1][0] + vyH_const * vyH_const * dmu[2][0] + vyH_const * vzH_const * dmu[3][0]) +
                                (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * dmu_const[0][0] +
                                (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH)) * dmu_const[1][0] +
                                2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVyH) * dmu_const[2][0] +
                                (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH)) * dmu_const[3][0] ) ;
    uuu_product_dmu[2][1] =  gamma_const * ( gamma_const * (vyH_const * dmu[0][1] + vxH_const * vyH_const * dmu[1][1] + vyH_const * vyH_const * dmu[2][1] + vyH_const * vzH_const * dmu[3][1]) +
                                (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * dmu_const[0][1] +
                                (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH)) * dmu_const[1][1] +
                                2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVyH) * dmu_const[2][1] +
                                (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH)) * dmu_const[3][1] ) ;
    uuu_product_dmu[2][2] =  gamma_const * ( gamma_const * (vyH_const * dmu[0][2] + vxH_const * vyH_const * dmu[1][2] + vyH_const * vyH_const * dmu[2][2] + vyH_const * vzH_const * dmu[3][2]) +
                                (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * dmu_const[0][2] +
                                (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH)) * dmu_const[1][2] +
                                2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVyH) * dmu_const[2][2] +
                                (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH)) * dmu_const[3][2] ) ;
    uuu_product_dmu[2][3] =  gamma_const * ( gamma_const * (vyH_const * dmu[0][3] + vxH_const * vyH_const * dmu[1][3] + vyH_const * vyH_const * dmu[2][3] + vyH_const * vzH_const * dmu[3][3]) +
                                (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * dmu_const[0][3] +
                                (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH)) * dmu_const[1][3] +
                                2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVyH) * dmu_const[2][3] +
                                (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH)) * dmu_const[3][3] ) ;
    uuu_product_dmu[3][0] =  gamma_const * ( gamma_const * (vzH_const * dmu[0][0] + vxH_const * vzH_const * dmu[1][0] + vyH_const * vzH_const * dmu[2][0] + vzH_const * vzH_const * dmu[3][0]) +
                                (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * dmu_const[0][0] +
                                (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH)) * dmu_const[1][0] +
                                (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH)) * dmu_const[2][0] +
                                2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVzH) * dmu_const[3][0] ) ;
    uuu_product_dmu[3][1] =  gamma_const * ( gamma_const * (vzH_const * dmu[0][1] + vxH_const * vzH_const * dmu[1][1] + vyH_const * vzH_const * dmu[2][1] + vzH_const * vzH_const * dmu[3][1]) +
                                (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * dmu_const[0][1] +
                                (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH)) * dmu_const[1][1] +
                                (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH)) * dmu_const[2][1] +
                                2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVzH) * dmu_const[3][1] ) ;
    uuu_product_dmu[3][2] =  gamma_const * ( gamma_const * (vzH_const * dmu[0][2] + vxH_const * vzH_const * dmu[1][2] + vyH_const * vzH_const * dmu[2][2] + vzH_const * vzH_const * dmu[3][2]) +
                                (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * dmu_const[0][2] +
                                (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH)) * dmu_const[1][2] +
                                (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH)) * dmu_const[2][2] +
                                2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVzH) * dmu_const[3][2] ) ;
    uuu_product_dmu[3][3] =  gamma_const * ( gamma_const * (vzH_const * dmu[0][3] + vxH_const * vzH_const * dmu[1][3] + vyH_const * vzH_const * dmu[2][3] + vzH_const * vzH_const * dmu[3][3]) +
                                (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * dmu_const[0][3] +
                                (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH)) * dmu_const[1][3] +
                                (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH)) * dmu_const[2][3] +
                                2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVzH) * dmu_const[3][3] ) ;
    
    double uuu_Tr_dmu[4][4];
    
    uuu_Tr_dmu[0][0] = gamma_const * (gamma_const * Tr_dmu + 2 * gamma_delta * Tr_dmu_const);
    uuu_Tr_dmu[0][1] = gamma_const * (gamma_const * vxH_const * Tr_dmu + (2 * gamma_delta * vxH_const + gamma_const * deltaVxH) * Tr_dmu_const );
    uuu_Tr_dmu[0][2] = gamma_const * (gamma_const * vyH_const * Tr_dmu + (2 * gamma_delta * vyH_const + gamma_const * deltaVyH) * Tr_dmu_const );
    uuu_Tr_dmu[0][3] = gamma_const * (gamma_const * vzH_const * Tr_dmu + (2 * gamma_delta * vzH_const + gamma_const * deltaVzH) * Tr_dmu_const );
    
    uuu_Tr_dmu[1][0] = uuu_Tr_dmu[0][1];
    uuu_Tr_dmu[1][1] = gamma_const * (gamma_const * vxH_const * vxH_const * Tr_dmu + 2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVxH) * Tr_dmu_const );
    uuu_Tr_dmu[1][2] = gamma_const * (gamma_const * vxH_const * vyH_const * Tr_dmu + (2 * gamma_delta * vxH_const * vyH_const +
                                                                                    gamma_const * (vxH_const * deltaVyH + vyH_const * deltaVxH) ) * Tr_dmu_const );
    uuu_Tr_dmu[1][3] = gamma_const * (gamma_const * vxH_const * vzH_const * Tr_dmu + (2 * gamma_delta * vxH_const * vzH_const +
                                                                                    gamma_const * (vxH_const * deltaVzH + vzH_const * deltaVxH) ) * Tr_dmu_const );
    uuu_Tr_dmu[2][0] = uuu_Tr_dmu[0][2];
    uuu_Tr_dmu[2][1] = uuu_Tr_dmu[1][2];
    uuu_Tr_dmu[2][2] = gamma_const * (gamma_const * vyH_const * vyH_const * Tr_dmu + 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVyH) * Tr_dmu_const );
    uuu_Tr_dmu[2][3] = gamma_const * (gamma_const * vyH_const * vzH_const * Tr_dmu + (2 * gamma_delta * vyH_const * vzH_const +
                                                                                    gamma_const * (vyH_const * deltaVzH + vzH_const * deltaVyH) ) * Tr_dmu_const );
    
    uuu_Tr_dmu[3][0] = uuu_Tr_dmu[0][3];
    uuu_Tr_dmu[3][1] = uuu_Tr_dmu[1][3];
    uuu_Tr_dmu[3][2] = uuu_Tr_dmu[2][3];
    uuu_Tr_dmu[3][3] = gamma_const * (gamma_const * vzH_const * vzH_const * Tr_dmu + 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVzH) * Tr_dmu_const );

        // calculation of Z[mu][nu][lambda][rho]
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    for (int l = 0; l < 4; l++) Z[i][j][k][l] = 0.0;
        // filling Z matrix
        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
                for (int lam = 0; lam < 4; lam++)
                    for (int rho = 0; rho < 4; rho++) {
//                        if (nu == rho)
//                            Z[mu][nu][lam][rho] += 0.5 * (gmunu[mu][lam] - uuu[mu] * uuu[lam]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//                        if (mu == rho)
//                            Z[mu][nu][lam][rho] += 0.5 * (gmunu[nu][lam] - uuu[nu] * uuu[lam]);
//                        if (lam == rho)
//                            Z[mu][nu][lam][rho] -= (gmunu[mu][nu] - uuu[mu] * uuu[nu]) / 3.0;
                        if (nu == rho)
                            Z[mu][nu][lam][rho] += 0.5 * (gmunu[mu][lam] - uuu_product[mu][lam]);
                        if (mu == rho)
                            Z[mu][nu][lam][rho] += 0.5 * (gmunu[nu][lam] - uuu_product[nu][lam]);
                        if (lam == rho)
                            Z[mu][nu][lam][rho] -= (gmunu[mu][nu] - uuu_product[mu][nu]) / 3.0;
                    }

    // calculation of Z_dmu[mu][nu][lambda][rho]
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) Z_dmu[i][j] = 0.0;
//    // filling Z matrix
//    for (int mu = 0; mu < 4; mu++)
//        for (int nu = 0; nu < 4; nu++)
//            for (int lam = 0; lam < 4; lam++)
//                for (int rho = 0; rho < 4; rho++) {
//                    if (nu == rho)
//                        Z_dmu[mu][nu][lam][rho] += 0.5 * (gmunu[mu][lam] * dmu[lam][nu] - uuu_product_dmu[mu][nu]);
//                    if (mu == rho)
//                        Z_dmu[mu][nu][lam][rho] += 0.5 * (gmunu[nu][lam] * dmu[lam][mu] - uuu_product_dmu[nu][mu]);
//                    if (lam == rho)
//                        Z_dmu[mu][nu][lam][rho] -= (gmunu[mu][nu] * dmu[lam][lam] - uuu_Tr_dmu[mu][nu]) / 3.0;
//                }
//
    // filling Z matrix?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    for (int mu = 0; mu < 4; mu++){
        for (int nu = 0; nu < 4; nu++){
            Z_dmu[mu][nu] = - 0.5 * (uuu_product_dmu[mu][nu] + uuu_product_dmu[nu][mu]) + uuu_Tr_dmu[mu][nu]/3.0;// this part is already summed over lambda, since I have computed the whole product
            for (int lam = 0; lam < 4; lam++){
                
                Z_dmu[mu][nu] += 0.5 * (gmunu[mu][lam] * dmu[lam][nu] + gmunu[nu][lam] * dmu[lam][mu]) - gmunu[mu][nu] * dmu[lam][lam] / 3.0; // this part needs to be summed over lambda
            }
        }
    }
    
//        // calculating sigma[mu][nu]
//        for (int i = 0; i < 4; i++)
//            for (int j = 0; j < 4; j++) {
//                pi[i][j] = 0.0;
//                for (int k = 0; k < 4; k++)
//                    for (int l = 0; l < 4; l++) {
//                        pi[i][j] += Z[i][j][k][l] * dmu[k][l] * 2.0 * etaS * s / 5.068; //TADY SE NASOBI DELTA A DEL!!!!!!!!!!!!!!!!!!!!!!!......................
//                    }
//            }
    
//    // calculating sigma[mu][nu]
//    for (int i = 0; i < 4; i++)
//        for (int j = 0; j < 4; j++) {
//            pi[i][j] = 0.0;
//            for (int k = 0; k < 4; k++)
//                for (int l = 0; l < 4; l++) {
//                    pi[i][j] += Z_dmu[i][j][k][l] * 2.0 * etaS * s / 5.068;
//                }
//        }
    
    // calculating sigma[mu][nu]
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            pi[i][j] = 0.0;
        }
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {

            pi[i][j] += Z_dmu[i][j] * 2.0 * etaS * s / 5.068;
        }
        Pi = -zetaS * s * (dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3]) /
        5.068;  // fm^{-4} --> GeV/fm^3
        du = dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3];
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
void Hydro::setNSvalues() {
double e, p, nb, nq, ns, vx, vy, vz, piNS[4][4], PiNS, dmu[4][4], du;
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
    if (e <= 0.) continue;
    // NSquant(ix, iy, iz, piNS, PiNS, dmu, du) ;
    //############## set NS values assuming initial zero flow + Bjorken z
    // flow
    double T, mub, muq, mus;
    double etaS, zetaS;
    double s = eos->s(e, nb, nq, ns);  // entropy density in the current cell
    eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
    trcoeff->getEta(e,nb, T, etaS, zetaS);
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) piNS[i][j] = 0.0;  // reset piNS
    piNS[1][1] = piNS[2][2] = 2.0 / 3.0 * etaS * s / tau / 5.068;
    piNS[3][3] = -2.0 * piNS[1][1];
    PiNS = 0.0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++){
            c->setpi(i, j, piNS[i][j]);
        }
    c->setPi(PiNS);
   }
 cout << "setNS done\n";
}

void Hydro::ISformal(double e_const, double p_const, double vx_const, double vy_const, double vz_const, double eH_const, double pH_const, double vxH_const, double vyH_const, double vzH_const, double ePrev_const, double pPrev_const, double vxPrev_const, double vyPrev_const, double vzPrev_const) {
 double e, p, nb, nq, ns, vx, vy, vz, T, mub, muq, mus;
 double piNS[4][4], sigNS[4][4], PiNS, dmu[4][4], dmu_const[4][4], du, pi[4][4], piH[4][4], Pi, PiH;
 const double gmumu[4] = {1., -1., -1., -1.};
 #ifdef CARTESIAN
 double tauMinusHalf = 1.0;
 double tauMinusDt = 1.0;
 #else
 double tauMinusHalf = tau - 0.5 * dt;
 double tauMinusDt = tau - dt;
 #endif
 // loop #1 (relaxation+source terms)
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vx, vy,
                         vz);  // instead of getPrimVar()
       
    double deltaVx = vx - vxH_const;
    double deltaVy = vy - vyH_const;
    double deltaVz = vz - vzH_const;
       
    if (e <= 0.) {             // empty cell?
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
       c->setpiH0(i, j, 0.0);
       c->setpi0(i, j, 0.0);
      }
     c->setPiH0(0.0);
     c->setPi0(0.0);
    } else {  // non-empty cell
     // 1) relaxation(pi)+source(pi) terms for half-step
//     double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
     double gamma_const = 1. + 0.5 * (vxH_const * vxH_const + vyH_const * vyH_const + vzH_const * vzH_const);
     double gamma_delta = vxH_const * deltaVx + vyH_const * deltaVy + vzH_const * deltaVz;
     double gamma = gamma_const + gamma_delta;
     double inverse_gamma_const = 1. - 0.5 * (vxH_const * vxH_const + vyH_const * vyH_const + vzH_const * vzH_const);
     double inverse_gamma = inverse_gamma_const - gamma_delta;
     double u[4];
//     u[0] = gamma;
//     u[1] = u[0] * vx; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//     u[2] = u[0] * vy;
//     u[3] = u[0] * vz;
     u[0] = gamma;
     u[1] = u[0] * vxH_const + deltaVx * gamma_const;
     u[2] = u[0] * vyH_const + deltaVy * gamma_const;
     u[3] = u[0] * vzH_const + deltaVz * gamma_const;
        
//     double gamma_u_product[4]; //i.e.u[0] * u[i] for i=0, i<4
//     gamma_u_product[0] = gamma_const * (gamma + gamma_delta);
//     gamma_u_product[1] = gamma_const * ( (gamma + gamma_delta) * vxH_const + gamma_const * deltaVx );
//     gamma_u_product[2] = gamma_const * ( (gamma + gamma_delta) * vyH_const + gamma_const * deltaVy );
//     gamma_u_product[3] = gamma_const * ( (gamma + gamma_delta) * vzH_const + gamma_const * deltaVz );
        
     double u_u_product[4][4];
     //diagonal
     u_u_product[0][0] = gamma_const * (u[0] + gamma_delta);
     u_u_product[1][1] = gamma_const * ( (u[0] + gamma_delta) * vxH_const * vxH_const + 2 * gamma_const * vxH_const * deltaVx );
     u_u_product[2][2] = gamma_const * ( (u[0] + gamma_delta) * vyH_const * vyH_const + 2 * gamma_const * vyH_const * deltaVy );
     u_u_product[3][3] = gamma_const * ( (u[0] + gamma_delta) * vzH_const * vzH_const + 2 * gamma_const * vzH_const * deltaVz );
     //off-diagonal
     u_u_product[0][1] = gamma_const * ( (u[0] + gamma_delta) * vxH_const + gamma_const * deltaVx );
     u_u_product[0][2] = gamma_const * ( (u[0] + gamma_delta) * vyH_const + gamma_const * deltaVy );
     u_u_product[0][3] = gamma_const * ( (u[0] + gamma_delta) * vzH_const + gamma_const * deltaVz );

     u_u_product[1][0] = u_u_product[0][1];
     u_u_product[1][2] = gamma_const * ( (u[0] + gamma_delta) * vxH_const * vyH_const + gamma_const * vxH_const * deltaVy + gamma_const * vyH_const * deltaVx );
     u_u_product[1][3] = gamma_const * ( (u[0] + gamma_delta) * vxH_const * vzH_const + gamma_const * vxH_const * deltaVz + gamma_const * vzH_const * deltaVx );

     u_u_product[2][0] = u_u_product[0][2];
     u_u_product[2][1] = u_u_product[1][2];
     u_u_product[2][3] = gamma_const * ( (u[0] + gamma_delta) * vyH_const * vzH_const + gamma_const * vyH_const * deltaVz + gamma_const * vzH_const * deltaVy );
    
     u_u_product[3][0] = u_u_product[0][3];
     u_u_product[3][1] = u_u_product[1][3];
     u_u_product[3][2] = u_u_product[2][3];
    
     double u_u_product_inverse_gamma[4][4];
     u_u_product_inverse_gamma[0][0] = inverse_gamma_const * u_u_product[0][0] - gamma_delta * gamma_const * gamma_const;
     u_u_product_inverse_gamma[1][1] = inverse_gamma_const * u_u_product[1][1] - gamma_delta * gamma_const * gamma_const * vxH_const * vxH_const;
     u_u_product_inverse_gamma[2][2] = inverse_gamma_const * u_u_product[2][2] - gamma_delta * gamma_const * gamma_const * vyH_const * vyH_const;
     u_u_product_inverse_gamma[3][3] = inverse_gamma_const * u_u_product[3][3] - gamma_delta * gamma_const * gamma_const * vzH_const * vzH_const;
     //off-diagonal
     u_u_product_inverse_gamma[0][1] = inverse_gamma_const * u_u_product[0][1] - gamma_delta * gamma_const * gamma_const * vxH_const;
     u_u_product_inverse_gamma[0][2] = inverse_gamma_const * u_u_product[0][2] - gamma_delta * gamma_const * gamma_const * vyH_const;
     u_u_product_inverse_gamma[0][3] = inverse_gamma_const * u_u_product[0][3] - gamma_delta * gamma_const * gamma_const * vzH_const;

     u_u_product_inverse_gamma[1][0] = u_u_product_inverse_gamma[0][1];
     u_u_product_inverse_gamma[1][2] = inverse_gamma_const * u_u_product[1][2] - gamma_delta * gamma_const * gamma_const * vxH_const * vyH_const;
     u_u_product_inverse_gamma[1][3] = inverse_gamma_const * u_u_product[1][3] - gamma_delta * gamma_const * gamma_const * vxH_const * vzH_const;

     u_u_product_inverse_gamma[2][0] = u_u_product_inverse_gamma[0][2];
     u_u_product_inverse_gamma[2][1] = u_u_product_inverse_gamma[1][2];
     u_u_product_inverse_gamma[2][3] = inverse_gamma_const * u_u_product[2][3] - gamma_delta * gamma_const * gamma_const * vyH_const * vzH_const;
       
     u_u_product_inverse_gamma[3][0] = u_u_product_inverse_gamma[0][3];
     u_u_product_inverse_gamma[3][1] = u_u_product_inverse_gamma[1][3];
     u_u_product_inverse_gamma[3][2] = u_u_product_inverse_gamma[2][3];
        
        // source term  + tau*delta_Q_i/delta_tau
        double flux[4];
        for (int i = 0; i < 4; i++)
   //      flux[i] = tauMinusDt * (c->getpi(0, i) + c->getPi() * u[0] * u[i]); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        flux[i] = tauMinusDt * (c->getpi(0, i) + c->getPi() * u_u_product[0][i]);
        flux[0] += -tauMinusDt * c->getPi();
          // cout << c->getPi() << endl;
        c->addFlux(flux[0], flux[1], flux[2], flux[3], 0., 0., 0.);
        // now calculating viscous terms in NS limit
        NSquant(ix, iy, iz, piNS, PiNS, dmu, dmu_const, du, e_const, p_const, vx_const, vy_const, vz_const, eH_const, pH_const, vxH_const, vyH_const, vzH_const, ePrev_const, pPrev_const, vxPrev_const, vyPrev_const, vzPrev_const);
         
     double u_u_product_dmu[4][4];
     u_u_product_dmu[0][0] =  gamma_const * ( gamma_const * (dmu[0][0] + vxH_const * dmu[1][0] + vyH_const * dmu[2][0] + vzH_const * dmu[3][0]) +
                                 2 * gamma_delta * dmu_const[0][0] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[1][0] +
                                 (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[2][0] + (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[3][0] ) ;
     u_u_product_dmu[0][1] =  gamma_const * ( gamma_const * (dmu[0][1] + vxH_const * dmu[1][1] + vyH_const * dmu[2][1] + vzH_const * dmu[3][1]) +
                                 2 * gamma_delta * dmu_const[0][1] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[1][1] +
                                 (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[2][1] + (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[3][1] ) ;
     u_u_product_dmu[0][2] =  gamma_const * ( gamma_const * (dmu[0][2] + vxH_const * dmu[1][2] + vyH_const * dmu[2][2] + vzH_const * dmu[3][2]) +
                                 2 * gamma_delta * dmu_const[0][2] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[1][2] +
                                 (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[2][2] + (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[3][2] ) ;
     u_u_product_dmu[0][3] =  gamma_const * ( gamma_const * (dmu[0][3] + vxH_const * dmu[1][3] + vyH_const * dmu[2][3] + vzH_const * dmu[3][3]) +
                                 2 * gamma_delta * dmu_const[0][3] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[1][3] +
                                 (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[2][3] + (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[3][3] ) ;
     u_u_product_dmu[1][0] =  gamma_const * ( gamma_const * (vxH_const * dmu[0][0] + vxH_const * vxH_const * dmu[1][0] + vxH_const * vyH_const * dmu[2][0] + vxH_const * vzH_const * dmu[3][0]) +
                                 (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][0] +
                                 2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVx) * dmu_const[1][0] +
                                 (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVy + vyH_const * deltaVx)) * dmu_const[2][0] +
                                 (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVz + vzH_const * deltaVx)) * dmu_const[3][0] ) ;
     u_u_product_dmu[1][1] =  gamma_const * ( gamma_const * (vxH_const * dmu[0][1] + vxH_const * vxH_const * dmu[1][1] + vxH_const * vyH_const * dmu[2][1] + vxH_const * vzH_const * dmu[3][1]) +
                                 (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][1] +
                                 2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVx) * dmu_const[1][1] +
                                 (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVy + vyH_const * deltaVx)) * dmu_const[2][1] +
                                 (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVz + vzH_const * deltaVx)) * dmu_const[3][1] ) ;
     u_u_product_dmu[1][2] =  gamma_const * ( gamma_const * (vxH_const * dmu[0][2] + vxH_const * vxH_const * dmu[1][2] + vxH_const * vyH_const * dmu[2][2] + vxH_const * vzH_const * dmu[3][2]) +
                                 (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][2] +
                                 2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVx) * dmu_const[1][2] +
                                 (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVy + vyH_const * deltaVx)) * dmu_const[2][2] +
                                 (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVz + vzH_const * deltaVx)) * dmu_const[3][2] ) ;
     u_u_product_dmu[1][3] =  gamma_const * ( gamma_const * (vxH_const * dmu[0][3] + vxH_const * vxH_const * dmu[1][3] + vxH_const * vyH_const * dmu[2][3] + vxH_const * vzH_const * dmu[3][3]) +
                                 (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][3] +
                                 2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const * deltaVx) * dmu_const[1][3] +
                                 (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVy + vyH_const * deltaVx)) * dmu_const[2][3] +
                                 (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVz + vzH_const * deltaVx)) * dmu_const[3][3] ) ;
     u_u_product_dmu[2][0] =  gamma_const * ( gamma_const * (vyH_const * dmu[0][0] + vxH_const * vyH_const * dmu[1][0] + vyH_const * vyH_const * dmu[2][0] + vyH_const * vzH_const * dmu[3][0]) +
                                 (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[0][0] +
                                 (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVy + vyH_const * deltaVx)) * dmu_const[1][0] +
                                 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVy) * dmu_const[2][0] +
                                 (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVz + vzH_const * deltaVy)) * dmu_const[3][0] ) ;
     u_u_product_dmu[2][1] =  gamma_const * ( gamma_const * (vyH_const * dmu[0][1] + vxH_const * vyH_const * dmu[1][1] + vyH_const * vyH_const * dmu[2][1] + vyH_const * vzH_const * dmu[3][1]) +
                                 (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[0][1] +
                                 (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVy + vyH_const * deltaVx)) * dmu_const[1][1] +
                                 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVy) * dmu_const[2][1] +
                                 (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVz + vzH_const * deltaVy)) * dmu_const[3][1] ) ;
     u_u_product_dmu[2][2] =  gamma_const * ( gamma_const * (vyH_const * dmu[0][2] + vxH_const * vyH_const * dmu[1][2] + vyH_const * vyH_const * dmu[2][2] + vyH_const * vzH_const * dmu[3][2]) +
                                 (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[0][2] +
                                 (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVy + vyH_const * deltaVx)) * dmu_const[1][2] +
                                 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVy) * dmu_const[2][2] +
                                 (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVz + vzH_const * deltaVy)) * dmu_const[3][2] ) ;
     u_u_product_dmu[2][3] =  gamma_const * ( gamma_const * (vyH_const * dmu[0][3] + vxH_const * vyH_const * dmu[1][3] + vyH_const * vyH_const * dmu[2][3] + vyH_const * vzH_const * dmu[3][3]) +
                                 (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[0][3] +
                                 (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const * deltaVy + vyH_const * deltaVx)) * dmu_const[1][3] +
                                 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const * deltaVy) * dmu_const[2][3] +
                                 (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVz + vzH_const * deltaVy)) * dmu_const[3][3] ) ;
     u_u_product_dmu[3][0] =  gamma_const * ( gamma_const * (vzH_const * dmu[0][0] + vxH_const * vzH_const * dmu[1][0] + vyH_const * vzH_const * dmu[2][0] + vzH_const * vzH_const * dmu[3][0]) +
                                 (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[0][0] +
                                 (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVz + vzH_const * deltaVx)) * dmu_const[1][0] +
                                 (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVz + vzH_const * deltaVy)) * dmu_const[2][0] +
                                 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVz) * dmu_const[3][0] ) ;
     u_u_product_dmu[3][1] =  gamma_const * ( gamma_const * (vzH_const * dmu[0][1] + vxH_const * vzH_const * dmu[1][1] + vyH_const * vzH_const * dmu[2][1] + vzH_const * vzH_const * dmu[3][1]) +
                                 (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[0][1] +
                                 (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVz + vzH_const * deltaVx)) * dmu_const[1][1] +
                                 (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVz + vzH_const * deltaVy)) * dmu_const[2][1] +
                                 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVz) * dmu_const[3][1] ) ;
     u_u_product_dmu[3][2] =  gamma_const * ( gamma_const * (vzH_const * dmu[0][2] + vxH_const * vzH_const * dmu[1][2] + vyH_const * vzH_const * dmu[2][2] + vzH_const * vzH_const * dmu[3][2]) +
                                 (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[0][2] +
                                 (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVz + vzH_const * deltaVx)) * dmu_const[1][2] +
                                 (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVz + vzH_const * deltaVy)) * dmu_const[2][2] +
                                 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVz) * dmu_const[3][2] ) ;
     u_u_product_dmu[3][3] =  gamma_const * ( gamma_const * (vzH_const * dmu[0][3] + vxH_const * vzH_const * dmu[1][3] + vyH_const * vzH_const * dmu[2][3] + vzH_const * vzH_const * dmu[3][3]) +
                                 (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[0][3] +
                                 (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const * deltaVz + vzH_const * deltaVx)) * dmu_const[1][3] +
                                 (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const * deltaVz + vzH_const * deltaVy)) * dmu_const[2][3] +
                                 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const * deltaVz) * dmu_const[3][3] ) ;
        
        double u_u_product_dmu_L[4][4][4];
        u_u_product_dmu_L[0][0][0] = gamma_const * ( gamma_const * dmu[0][0] + 2 * gamma_delta * dmu_const[0][0] );
        u_u_product_dmu_L[0][0][1] = gamma_const * ( gamma_const * vxH_const * dmu[1][0] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[1][0] );
        u_u_product_dmu_L[0][0][2] = gamma_const * ( gamma_const * vyH_const * dmu[2][0] + (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[2][0] );
        u_u_product_dmu_L[0][0][3] = gamma_const * ( gamma_const * vzH_const * dmu[3][0] + (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[3][0] );
        
        u_u_product_dmu_L[0][1][0] = gamma_const * ( gamma_const * dmu[0][1] + 2 * gamma_delta * dmu_const[0][1] );
        u_u_product_dmu_L[0][1][1] = gamma_const * ( gamma_const * vxH_const * dmu[1][1] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[1][0] );
        u_u_product_dmu_L[0][1][2] = gamma_const * ( gamma_const * vyH_const * dmu[2][1] + (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[2][0] );
        u_u_product_dmu_L[0][1][3] = gamma_const * ( gamma_const * vzH_const * dmu[3][1] + (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[3][0] );
        
        u_u_product_dmu_L[0][2][0] = gamma_const * ( gamma_const * dmu[0][2] + 2 * gamma_delta * dmu_const[0][2] );
        u_u_product_dmu_L[0][2][1] = gamma_const * ( gamma_const * vxH_const * dmu[1][2] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[1][0] );
        u_u_product_dmu_L[0][2][2] = gamma_const * ( gamma_const * vyH_const * dmu[2][2] + (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[2][0] );
        u_u_product_dmu_L[0][2][3] = gamma_const * ( gamma_const * vzH_const * dmu[3][2] + (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[3][0] );
        
        u_u_product_dmu_L[0][3][0] = gamma_const * ( gamma_const * dmu[0][3] + 2 * gamma_delta * dmu_const[0][3] );
        u_u_product_dmu_L[0][3][1] = gamma_const * ( gamma_const * vxH_const * dmu[1][3] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[1][0] );
        u_u_product_dmu_L[0][3][2] = gamma_const * ( gamma_const * vyH_const * dmu[2][3] + (2 * gamma_delta * vyH_const + gamma_const * deltaVy) * dmu_const[2][0] );
        u_u_product_dmu_L[0][3][3] = gamma_const * ( gamma_const * vzH_const * dmu[3][3] + (2 * gamma_delta * vzH_const + gamma_const * deltaVz) * dmu_const[3][0] );
                                    

        u_u_product_dmu_L[1][0][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][0] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][0]);
                                                                   
        u_u_product_dmu_L[1][0][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][0] + 2 * (gamma_delta * vxH_const * vxH_const + gamma_const *                                            vxH_const * deltaVx) * dmu_const[1][0] );
        u_u_product_dmu_L[1][0][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][0] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                                 * deltaVy + vyH_const * deltaVx)) * dmu_const[2][0] );
        u_u_product_dmu_L[1][0][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][0] + (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const                                 * deltaVz + vzH_const * deltaVx)) * dmu_const[3][0] );
                                    
        u_u_product_dmu_L[1][1][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][1] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][1]);
                                                                                                              
        u_u_product_dmu_L[1][1][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][1] + 2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const                                  * deltaVx) * dmu_const[1][1] );
        u_u_product_dmu_L[1][1][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][1] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                                 * deltaVy + vyH_const * deltaVx)) * dmu_const[2][1] );
        u_u_product_dmu_L[1][1][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][1] + (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const                                 * deltaVz + vzH_const * deltaVx)) * dmu_const[3][1] );
                                                  
        u_u_product_dmu_L[1][2][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][2] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][2]);

        u_u_product_dmu_L[1][2][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][2] + 2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const                         * deltaVx) * dmu_const[1][2] );
        u_u_product_dmu_L[1][2][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][2] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[2][2] );
        u_u_product_dmu_L[1][2][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][2] + (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const                        * deltaVz + vzH_const * deltaVx)) * dmu_const[3][2] );
                                                   
        u_u_product_dmu_L[1][3][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][3] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][3]);

        u_u_product_dmu_L[1][3][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][3] + 2 * (gamma_delta * vxH_const * vxH_const + gamma_const * vxH_const                         * deltaVx) * dmu_const[1][3] );
        u_u_product_dmu_L[1][3][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][3] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[2][3] );
        u_u_product_dmu_L[1][3][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][3] + (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const                        * deltaVz + vzH_const * deltaVx)) * dmu_const[3][3] );
        
        
        u_u_product_dmu_L[2][0][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][0] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][0]);

        u_u_product_dmu_L[2][0][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][0] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[1][0] );
        u_u_product_dmu_L[2][0][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][0] + 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const                         * deltaVy) * dmu_const[2][0] );
        u_u_product_dmu_L[2][0][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][0] + (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const                        * deltaVz + vzH_const * deltaVx)) * dmu_const[3][0] );
        
        u_u_product_dmu_L[2][1][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][1] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][1]);

        u_u_product_dmu_L[2][1][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][1] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[1][1] );
        u_u_product_dmu_L[2][1][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][1] + 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const                         * deltaVy) * dmu_const[2][1] );
        u_u_product_dmu_L[2][1][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][1] + (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const                        * deltaVz + vzH_const * deltaVx)) * dmu_const[3][1] );
        
        u_u_product_dmu_L[2][2][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][2] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][2]);

        u_u_product_dmu_L[2][2][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][2] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[1][2] );
        u_u_product_dmu_L[2][2][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][2] + 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const                         * deltaVy) * dmu_const[2][2] );
        u_u_product_dmu_L[2][2][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][2] + (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const                        * deltaVz + vzH_const * deltaVx)) * dmu_const[3][2] );
        
        u_u_product_dmu_L[2][3][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][3] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][3]);

        u_u_product_dmu_L[2][3][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][3] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[1][3] );
        u_u_product_dmu_L[2][3][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][3] + 2 * (gamma_delta * vyH_const * vyH_const + gamma_const * vyH_const                         * deltaVy) * dmu_const[2][3] );
        u_u_product_dmu_L[2][3][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][3] + (2 * gamma_delta * vxH_const * vzH_const + gamma_const * (vxH_const                        * deltaVz + vzH_const * deltaVx)) * dmu_const[3][3] );
        
        u_u_product_dmu_L[3][0][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][0] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][0]);

        u_u_product_dmu_L[3][0][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][0] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[1][0] );
        u_u_product_dmu_L[3][0][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][0] + (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const                        * deltaVz + vzH_const * deltaVy)) * dmu_const[2][0] );
        u_u_product_dmu_L[3][0][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][0] + 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const                         * deltaVz) * dmu_const[3][0] );
        
        u_u_product_dmu_L[3][1][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][1] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][1]);

        u_u_product_dmu_L[3][1][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][1] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[1][1] );
        u_u_product_dmu_L[3][1][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][1] + (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const                        * deltaVz + vzH_const * deltaVy)) * dmu_const[2][1] );
        u_u_product_dmu_L[3][1][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][1] + 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const                         * deltaVz) * dmu_const[3][1] );
        
        u_u_product_dmu_L[3][2][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][2] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][2]);

        u_u_product_dmu_L[3][2][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][2] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[1][2] );
        u_u_product_dmu_L[3][2][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][2] + (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const                        * deltaVz + vzH_const * deltaVy)) * dmu_const[2][2] );
        u_u_product_dmu_L[3][2][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][2] + 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const                         * deltaVz) * dmu_const[3][2] );
        
        u_u_product_dmu_L[3][3][0] = gamma_const * ( gamma_const * vxH_const * dmu[0][3] + (2 * gamma_delta * vxH_const + gamma_const * deltaVx) * dmu_const[0][3]);

        u_u_product_dmu_L[3][3][1] = gamma_const * ( gamma_const * vxH_const * vxH_const * dmu[1][3] + (2 * gamma_delta * vxH_const * vyH_const + gamma_const * (vxH_const                        * deltaVy + vyH_const * deltaVx)) * dmu_const[1][3] );
        u_u_product_dmu_L[3][3][2] = gamma_const * ( gamma_const * vxH_const * vyH_const * dmu[2][3] + (2 * gamma_delta * vyH_const * vzH_const + gamma_const * (vyH_const                        * deltaVz + vzH_const * deltaVy)) * dmu_const[2][3] );
        u_u_product_dmu_L[3][3][3] = gamma_const * ( gamma_const * vxH_const * vzH_const * dmu[3][3] + 2 * (gamma_delta * vzH_const * vzH_const + gamma_const * vzH_const                         * deltaVz) * dmu_const[3][3] );
        
 
        
     double u_u_product_dmu_inverse_gamma[4][4];
     u_u_product_dmu_inverse_gamma[0][0] = inverse_gamma_const * u_u_product_dmu[0][0] - gamma_delta * gamma_const * ( gamma_const * (dmu[0][0] + vxH_const * dmu[1][0] + vyH_const * dmu[2][0] + vzH_const * dmu[3][0]) ) ;
     u_u_product_dmu_inverse_gamma[0][1] = inverse_gamma_const * u_u_product_dmu[0][1] - gamma_delta * gamma_const * ( gamma_const * (dmu[0][1] + vxH_const * dmu[1][1] + vyH_const * dmu[2][1] + vzH_const * dmu[3][1]) ) ;
     u_u_product_dmu_inverse_gamma[0][2] = inverse_gamma_const * u_u_product_dmu[0][2] - gamma_delta * gamma_const * ( gamma_const * (dmu[0][2] + vxH_const * dmu[1][2] + vyH_const * dmu[2][2] + vzH_const * dmu[3][2]) ) ;
     u_u_product_dmu_inverse_gamma[0][3] = inverse_gamma_const * u_u_product_dmu[0][3] - gamma_delta * gamma_const * ( gamma_const * (dmu[0][3] + vxH_const * dmu[1][3] + vyH_const * dmu[2][3] + vzH_const * dmu[3][3]) ) ;
     u_u_product_dmu_inverse_gamma[1][0] = inverse_gamma_const * u_u_product_dmu[1][0] - gamma_delta * gamma_const * ( gamma_const * (vxH_const * dmu[0][0] + vxH_const * vxH_const * dmu[1][0] + vxH_const * vyH_const *
                                                                                                                         dmu[2][0] + vxH_const * vzH_const * dmu[3][0]) ) ;
     u_u_product_dmu_inverse_gamma[1][1] = inverse_gamma_const * u_u_product_dmu[1][1] - gamma_delta * gamma_const * ( gamma_const * (vxH_const * dmu[0][1] + vxH_const * vxH_const * dmu[1][1] + vxH_const * vyH_const *
                                                                                                                         dmu[2][1] + vxH_const * vzH_const * dmu[3][1]) ) ;
     u_u_product_dmu_inverse_gamma[1][2] = inverse_gamma_const * u_u_product_dmu[1][2] - gamma_delta * gamma_const * ( gamma_const * (vxH_const * dmu[0][2] + vxH_const * vxH_const * dmu[1][2] + vxH_const * vyH_const *
                                                                                                                         dmu[2][2] + vxH_const * vzH_const * dmu[3][2]) ) ;
     u_u_product_dmu_inverse_gamma[1][3] = inverse_gamma_const * u_u_product_dmu[1][3] - gamma_delta * gamma_const * ( gamma_const * (vxH_const * dmu[0][3] + vxH_const * vxH_const * dmu[1][3] + vxH_const * vyH_const *
                                                                                                                         dmu[2][3] + vxH_const * vzH_const * dmu[3][3]) ) ;
     u_u_product_dmu_inverse_gamma[2][0] = inverse_gamma_const * u_u_product_dmu[2][0] - gamma_delta * gamma_const * ( gamma_const * (vyH_const * dmu[0][0] + vxH_const * vyH_const * dmu[1][0] + vyH_const * vyH_const *
                                                                                                                         dmu[2][0] + vyH_const * vzH_const * dmu[3][0]) ) ;
     u_u_product_dmu_inverse_gamma[2][1] = inverse_gamma_const * u_u_product_dmu[2][1] - gamma_delta * gamma_const * ( gamma_const * (vyH_const * dmu[0][1] + vxH_const * vyH_const * dmu[1][1] + vyH_const * vyH_const *
                                                                                                                         dmu[2][1] + vyH_const * vzH_const * dmu[3][1]) ) ;
     u_u_product_dmu_inverse_gamma[2][2] = inverse_gamma_const * u_u_product_dmu[2][2] - gamma_delta * gamma_const * ( gamma_const * (vyH_const * dmu[0][2] + vxH_const * vyH_const * dmu[1][2] + vyH_const * vyH_const *
                                                                                                                         dmu[2][2] + vyH_const * vzH_const * dmu[3][2]) ) ;
     u_u_product_dmu_inverse_gamma[2][3] = inverse_gamma_const * u_u_product_dmu[2][3] - gamma_delta * gamma_const * ( gamma_const * (vyH_const * dmu[0][3] + vxH_const * vyH_const * dmu[1][3] + vyH_const * vyH_const *
                                                                                                                         dmu[2][3] + vyH_const * vzH_const * dmu[3][3]) ) ;
     u_u_product_dmu_inverse_gamma[3][0] = inverse_gamma_const * u_u_product_dmu[3][0] - gamma_delta * gamma_const * ( gamma_const * (vzH_const * dmu[0][0] + vxH_const * vzH_const * dmu[1][0] + vyH_const * vzH_const *
                                                                                                                         dmu[2][0] + vzH_const * vzH_const * dmu[3][0]) ) ;
     u_u_product_dmu_inverse_gamma[3][1] = inverse_gamma_const * u_u_product_dmu[3][1] - gamma_delta * gamma_const * ( gamma_const * (vzH_const * dmu[0][1] + vxH_const * vzH_const * dmu[1][1] + vyH_const * vzH_const *
                                                                                                                         dmu[2][1] + vzH_const * vzH_const * dmu[3][1]) ) ;
     u_u_product_dmu_inverse_gamma[3][2] = inverse_gamma_const * u_u_product_dmu[3][2] - gamma_delta * gamma_const * ( gamma_const * (vzH_const * dmu[0][2] + vxH_const * vzH_const * dmu[1][2] + vyH_const * vzH_const *
                                                                                                                         dmu[2][2] + vzH_const * vzH_const * dmu[3][2]) ) ;
     u_u_product_dmu_inverse_gamma[3][3] = inverse_gamma_const * u_u_product_dmu[3][3] - gamma_delta * gamma_const * ( gamma_const * (vzH_const * dmu[0][3] + vxH_const * vzH_const * dmu[1][3] + vyH_const * vzH_const *
                                                                                                                           dmu[2][3] + vzH_const * vzH_const * dmu[3][3]) ) ;
                                                                                                                      
//     // source term  + tau*delta_Q_i/delta_tau
//     double flux[4];
//     for (int i = 0; i < 4; i++)
////      flux[i] = tauMinusDt * (c->getpi(0, i) + c->getPi() * u[0] * u[i]); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//     flux[i] = tauMinusDt * (c->getpi(0, i) + c->getPi() * u_u_product[0][i]);
//     flux[0] += -tauMinusDt * c->getPi();
//       // cout << c->getPi() << endl;
//     c->addFlux(flux[0], flux[1], flux[2], flux[3], 0., 0., 0.);
//     // now calculating viscous terms in NS limit
//     NSquant(ix, iy, iz, piNS, PiNS, dmu, dmu_const, du, e_const, p_const, vx_const, vy_const, vz_const, eH_const, pH_const, vxH_const, vyH_const, vzH_const, ePrev_const, pPrev_const, vxPrev_const, vyPrev_const, vzPrev_const);
     eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
     double etaS, zetaS;
     trcoeff->getEta(e,nb, T, etaS, zetaS);
     const double s = eos->s(e, nb, nq, ns);
     const double eta = etaS * s;
     // auxiliary variable sigmaNS = piNS / (2*eta), 
     // mainly to protect against division by zero in the eta=0 case.
     for(int i=0; i<4; i++)
     for(int j=0; j<4; j++) {
      sigNS[i][j] = 0.5 * piNS[i][j] / eta * 5.068;
      if(eta<=0.0) sigNS[i][j] = 0.0;
     }
     //############# get relaxation times
     double taupi, tauPi;
     trcoeff->getTau(e, nb, T, taupi, tauPi);
     double deltapipi, taupipi, lambdapiPi, phi7, delPiPi, lamPipi;
     trcoeff->getOther(e, nb, nq, ns, deltapipi, taupipi, lambdapiPi, phi7);
     phi7 = phi7/taupi;  // dividing by tau_pi here, to avoid NaNs when tau_pi==0
     if(taupi<0.5*dt)
      deltapipi = taupipi = lambdapiPi = phi7 = 0.0;
     trcoeff->getOtherBulk(e, nb, nq, ns, delPiPi, lamPipi);
     if(tauPi<0.5*dt)
      delPiPi = lamPipi = 0.0;
     //#############
     double Delta[10];
     double Delta_inverse_gamma[10];
     // relaxation term, piH,PiH-->half-step
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
//       Delta[index44(i, j)] = -u[i] * u[j]; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!?????????????
         Delta[index44(i, j)] = -u_u_product[i][j];
         Delta_inverse_gamma[index44(i, j)] = -u_u_product_inverse_gamma[i][j];
//       if (i == j) Delta[index44(i, j)] += gmumu[i];
         if (i == j){
            Delta[index44(i, j)] += gmumu[i];
            Delta_inverse_gamma[index44(i, j)] += gmumu[i];
         }
#ifdef FORMAL_SOLUTION
       c->setpiH0(i, j, (c->getpi(i, j) - piNS[i][j]) *
                                exp(-dt / 2.0 * inverse_gamma / taupi) +
                            piNS[i][j]);
#else
      if(taupi>0.5*dt)
       c->setpiH0(i, j, c->getpi(i, j) -
          (c->getpi(i, j) - piNS[i][j]) * dt / 2.0 * inverse_gamma / taupi);
      else
       c->setpiH0(i, j, piNS[i][j]);
#endif
      }
#ifdef FORMAL_SOLUTION
     c->setPiH0((c->getPi() - PiNS) * exp(-dt / 2.0 * inverse_gamma / tauPi) + PiNS);
#else
    if(tauPi>0.5*dt)
     c->setPiH0(c->getPi() - (c->getPi() - PiNS) * dt / 2.0 * inverse_gamma / tauPi);
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
       // now transversality and cross terms
       c->addpiH0(i, j, (- deltapipi * c->getpi(i, j) * du +
         lambdapiPi * c->getPi() * sigNS[i][j]) * inverse_gamma * 0.5 * dt);
       for (int k = 0; k < 4; k++) {
        // parts of terms with one internal summation index
        c->addpiH0(i, j, (phi7 * c->getpi(i,k) * c->getpi(j,k) * gmumu[k] - taupipi * 0.5 * (c->getpi(i,k) * sigNS[j][k] * gmumu[k] + c->getpi(j,k) * sigNS[i][k] * gmumu[k])) * inverse_gamma * 0.5 * dt);
        // parts of terms with two internal summation indexes
        for (int l = 0; l < 4; l++){
//         c->addpiH0(i, j, (-c->getpi(i, k) * u[j] - c->getpi(j, k) * u[i]) * u[l] * dmu[l][k] * gmumu[k] / gamma * 0.5 * dt // TADY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!???????????
           c->addpiH0(i, j, ((-c->getpi(i, k) * u_u_product_dmu_L[j][k][l] - c->getpi(j, k) * u_u_product_dmu_L[i][k][l]) * gmumu[k] * 0.5 * dt
          - 1. / 3. * Delta[index44(i,j)] * c->getpi(k, l) * ( phi7 * c->getpi(k, l) - taupipi * sigNS[k][l]) * gmumu[k] * gmumu[l] * 0.5 * dt)* inverse_gamma);//originaly u_u product [i/j][k] instead of [i/j][l]
         c->addPiH0(lamPipi * c->getpi(k, l) * sigNS[k][l] * inverse_gamma * 0.5 * dt);
        }
       }
      }
     c->addPiH0(-delPiPi * c->getPi() * du * inverse_gamma * 0.5 * dt);
     // 1) relaxation(piH)+source(piH) terms for full-step
     for (int i = 0; i < 4; i++)
      for (int j = 0; j <= i; j++) {
#ifdef FORMAL_SOLUTION
       c->setpi0(i, j,
                 (c->getpi(i, j) - piNS[i][j]) * exp(-dt * inverse_gamma / taupi) +
                     piNS[i][j]);
#else
      if(taupi>0.5*dt)
       c->setpi0(i, j, c->getpi(i, j) -
          (c->getpiH0(i, j) - piNS[i][j]) * dt * inverse_gamma / taupi);
      else
       c->setpi0(i, j, piNS[i][j]);
#endif
      }
#ifdef FORMAL_SOLUTION
     c->setPi0((c->getPi() - PiNS) * exp(-dt * inverse_gamma / tauPi) + PiNS);
#else
    if(tauPi>0.5*dt)
     c->setPi0(c->getPi() - (c->getPiH0() - PiNS) * dt * inverse_gamma / tauPi);
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
       // now transversality and cross terms
       c->addpi0(i, j, (- deltapipi * c->getpiH0(i, j) * du +
         lambdapiPi * c->getPiH0() * sigNS[i][j]) * inverse_gamma * dt);
       for (int k = 0; k < 4; k++) {
        // parts of terms with one internal summation index
        c->addpi0(i, j, (phi7 * c->getpi(i,k) * c->getpiH0(j,k) * gmumu[k] - taupipi * 0.5 * (c->getpiH0(i,k) * sigNS[j][k] * gmumu[k] + c->getpiH0(j,k) * sigNS[i][k] * gmumu[k])) * inverse_gamma * dt);
        for (int l = 0; l < 4; l++){
//         c->addpi0(i, j, ((-c->getpiH0(i, k) * u[j] - c->getpiH0(j, k) * u[i]) * u[l] * dmu[l][k] * gmumu[k]// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!????????????
         c->addpi0(i, j, ((-c->getpiH0(i, k) * u_u_product_dmu_L[j][k][l] - c->getpiH0(j, k) * u_u_product_dmu_L[i][k][l]) * gmumu[k]
          - 1. / 3. * Delta[index44(i,j)] * c->getpiH0(k, l) * ( phi7 * c->getpiH0(k, l) - taupipi * sigNS[k][l]) * gmumu[k] * gmumu[l])* inverse_gamma * dt);
         c->addPi0(lamPipi * c->getpiH0(k, l) * sigNS[k][l] * inverse_gamma * dt);
        }
       }
      }
     c->addPi0(-delPiPi * c->getPiH0() * du * inverse_gamma * dt);
    }  // end non-empty cell
   }   // end loop #1

 // 3) -- advection ---
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    Cell *c = f->getCell(ix, iy, iz);
    c->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vx, vy,
                         vz);  // getPrimVar() before
       
    double deltaVx = vx - vxH_const;
    double deltaVy = vy - vyH_const;
    double deltaVz = vz - vzH_const;
       
    if (e <= 0.) continue;
    double xm = -vx * dt / f->getDx();
    double ym = -vy * dt / f->getDy();
    double zm = -vz * dt / f->getDz() / tauMinusHalf;
    double xmH = -vx * dt / f->getDx() / 2.0;
    double ymH = -vy * dt / f->getDy() / 2.0;
    double zmH = -vz * dt / f->getDz() / tauMinusHalf / 2.0;
    double wx[2] = {(1. - fabs(xm)), fabs(xm)};
    double wy[2] = {(1. - fabs(ym)), fabs(ym)};
    double wz[2] = {(1. - fabs(zm)), fabs(zm)};
    double wxH[2] = {(1. - fabs(xmH)), fabs(xmH)};
    double wyH[2] = {(1. - fabs(ymH)), fabs(ymH)};
    double wzH[2] = {(1. - fabs(zmH)), fabs(zmH)};
       
    double xm0 = -vxH_const * dt / f->getDx();
    double ym0 = -vyH_const * dt / f->getDy();
    double zm0 = -vzH_const * dt / f->getDz() / tauMinusHalf;
    double xmH0 = -vxH_const * dt / f->getDx() / 2.0;
    double ymH0 = -vyH_const * dt / f->getDy() / 2.0;
    double zmH0 = -vzH_const * dt / f->getDz() / tauMinusHalf / 2.0;
         
    double deltaXm = -deltaVx * dt / f->getDx();
    double deltaYm = -deltaVy * dt / f->getDy();
    double deltaZm = -deltaVz * dt / f->getDz() / tauMinusHalf;
    double deltaXmH = -deltaVx * dt / f->getDx() / 2.0;
    double deltaYmH = -deltaVy * dt / f->getDy() / 2.0;
    double deltaZmH = -deltaVz * dt / f->getDz() / tauMinusHalf / 2.0;
         
         
         
    double www_product[2][2][2];
    www_product[0][0][0] = 1 - fabs(xm) - fabs(ym) - fabs(zm) + sgn(xm) * sgn(ym) * ( xm0 * ym0 + xm0 * deltaYm + ym0 * deltaXm ) + sgn(xm) * sgn(zm) * ( xm0 * zm0 + xm0 * deltaZm + zm0 * deltaXm ) + sgn(ym) * sgn(zm) * ( ym0 * zm0 + ym0 * deltaZm + zm0 * deltaYm ) - sgn(xm) * sgn(ym) * sgn(zm) * ( xm0 * ym0 * zm0 + xm0 * ym0 * deltaZm + xm0 * zm0 * deltaYm + ym0 * zm0 * deltaXm );
    www_product[0][0][1] = fabs(zm) - sgn(xm) * sgn(zm) * ( xm0 * zm0 + xm0 * deltaZm + zm0 * deltaXm) - sgn(ym) * sgn(zm) * ( ym0 * zm0 + ym0 * deltaZm + zm0 * deltaYm ) + sgn(xm) * sgn(ym) * sgn(zm) * ( xm0 * ym0 * zm0 + xm0 * ym0 * deltaZm + xm0 * zm0 * deltaYm + ym0 * zm0 * deltaXm );
    www_product[0][1][0] = fabs(ym) - sgn(xm) * sgn(ym) * ( xm0 * ym0 + xm0 * deltaYm + ym0 * deltaXm) - sgn(ym) * sgn(zm) * ( ym0 * zm0 + ym0 * deltaZm + zm0 * deltaYm ) + sgn(xm) * sgn(ym) * sgn(zm) * ( xm0 * ym0 * zm0 + xm0 * ym0 * deltaZm + xm0 * zm0 * deltaYm + ym0 * zm0 * deltaXm );
    www_product[1][0][0] = fabs(xm) - sgn(xm) * sgn(zm) * ( xm0 * zm0 + xm0 * deltaZm + zm0 * deltaXm) - sgn(xm) * sgn(ym) * ( xm0 * ym0 + xm0 * deltaYm + ym0 * deltaXm ) + sgn(xm) * sgn(ym) * sgn(zm) * ( xm0 * ym0 * zm0 + xm0 * ym0 * deltaZm + xm0 * zm0 * deltaYm + ym0 * zm0 * deltaXm );
    www_product[0][1][1] = sgn(ym) * sgn(zm) * ( ym0 * zm0 + ym0 * deltaZm + zm0 * deltaYm ) - sgn(xm) * sgn(ym) * sgn(zm) * ( xm0 * ym0 * zm0 + xm0 * ym0 * deltaZm + xm0 * zm0 * deltaYm + ym0 * zm0 * deltaXm );
    www_product[1][0][1] = sgn(xm) * sgn(zm) * ( xm0 * zm0 + xm0 * deltaZm + zm0 * deltaXm ) - sgn(xm) * sgn(ym) * sgn(zm) * ( xm0 * ym0 * zm0 + xm0 * ym0 * deltaZm + xm0 * zm0 * deltaYm + ym0 * zm0 * deltaXm );
    www_product[1][1][0] = sgn(xm) * sgn(ym) * ( xm0 * ym0 + xm0 * deltaYm + ym0 * deltaXm ) - sgn(xm) * sgn(ym) * sgn(zm) * ( xm0 * ym0 * zm0 + xm0 * ym0 * deltaZm + xm0 * zm0 * deltaYm + ym0 * zm0 * deltaXm );
    www_product[1][1][1] = sgn(xm) * sgn(ym) * sgn(zm) * ( xm0 * ym0 * zm0 + xm0 * ym0 * deltaZm + xm0 * zm0 * deltaYm + ym0 * zm0 * deltaXm );
       
    double www_product_H[2][2][2];
    www_product_H[0][0][0] = 1 - fabs(xmH) - fabs(ymH) - fabs(zmH) + sgn(xmH) * sgn(ymH) * ( xmH0 * ymH0 + xmH0 * deltaYmH + ymH0 * deltaXmH ) + sgn(xmH) * sgn(zmH) * ( xmH0 * zmH0 + xmH0 * deltaZmH + zmH0 * deltaXmH ) + sgn(ymH) * sgn(zmH) * ( ymH0 * zmH0 + ymH0 * deltaZmH + zmH0 * deltaYmH ) - sgn(xmH) * sgn(ymH) * sgn(zmH) * ( xmH0 * ymH0 * zmH0 + xmH0 * ymH0 * deltaZmH + xmH0 * zmH0 * deltaYmH + ymH0 * zmH0 * deltaXmH );
    www_product_H[0][0][1] = fabs(zmH) - sgn(xmH) * sgn(zmH) * ( xmH0 * zmH0 + xmH0 * deltaZmH + zmH0 * deltaXmH) - sgn(ymH) * sgn(zmH) * ( ymH0 * zmH0 + ymH0 * deltaZmH + zmH0 * deltaYmH ) + sgn(xmH) * sgn(ymH) * sgn(zmH) * ( xmH0 * ymH0 * zmH0 + xmH0 * ymH0 * deltaZmH + xmH0 * zmH0 * deltaYmH + ymH0 * zmH0 * deltaXmH );
    www_product_H[0][1][0] = fabs(ymH) - sgn(xmH) * sgn(ymH) * ( xmH0 * ymH0 + xmH0 * deltaYmH + ymH0 * deltaXmH) - sgn(ymH) * sgn(zmH) * ( ymH0 * zmH0 + ymH0 * deltaZmH + zmH0 * deltaYmH ) + sgn(xmH) * sgn(ymH) * sgn(zmH) * ( xmH0 * ymH0 * zmH0 + xmH0 * ymH0 * deltaZmH + xmH0 * zmH0 * deltaYmH + ymH0 * zmH0 * deltaXmH );
    www_product_H[1][0][0] = fabs(xmH) - sgn(xmH) * sgn(zmH) * ( xmH0 * zmH0 + xmH0 * deltaZmH + zmH0 * deltaXmH) - sgn(xmH) * sgn(ymH) * ( xmH0 * ymH0 + xmH0 * deltaYmH + ymH0 * deltaXmH ) + sgn(xmH) * sgn(ymH) * sgn(zmH) * ( xmH0 * ymH0 * zmH0 + xmH0 * ymH0 * deltaZmH + xmH0 * zmH0 * deltaYmH + ymH0 * zmH0 * deltaXmH );
    www_product_H[0][1][1] = sgn(ymH) * sgn(zmH) * ( ymH0 * zmH0 + ymH0 * deltaZmH + zmH0 * deltaYmH ) - sgn(xmH) * sgn(ymH) * sgn(zmH) * ( xmH0 * ymH0 * zmH0 + xmH0 * ymH0 * deltaZmH + xmH0 * zmH0 * deltaYmH + ymH0 * zmH0 * deltaXmH );
    www_product_H[1][0][1] = sgn(xmH) * sgn(zmH) * ( xmH0 * zmH0 + xmH0 * deltaZmH + zmH0 * deltaXmH ) - sgn(xmH) * sgn(ymH) * sgn(zmH) * ( xmH0 * ymH0 * zmH0 + xmH0 * ymH0 * deltaZmH + xmH0 * zmH0 * deltaYmH + ymH0 * zmH0 * deltaXmH );
    www_product_H[1][1][0] = sgn(xmH) * sgn(ymH) * ( xmH0 * ymH0 + xmH0 * deltaYmH + ymH0 * deltaXmH ) - sgn(xmH) * sgn(ymH) * sgn(zmH) * ( xmH0 * ymH0 * zmH0 + xmH0 * ymH0 * deltaZmH + xmH0 * zmH0 * deltaYmH + ymH0 * zmH0 * deltaXmH );
    www_product_H[1][1][1] = sgn(xmH) * sgn(ymH) * sgn(zmH) * ( xmH0 * ymH0 * zmH0 + xmH0 * ymH0 * deltaZmH + xmH0 * zmH0 * deltaYmH + ymH0 * zmH0 * deltaXmH );
       
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
//         pi[i][j] += wx[jx] * wy[jy] * wz[jz] * c1->getpi0(i, j); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!????????????????.....................
//         piH[i][j] += wxH[jx] * wyH[jy] * wzH[jz] * c1->getpiH0(i, j);
         pi[i][j] += www_product[jx][jy][jz] * c1->getpi0(i, j);
         piH[i][j] += www_product_H[jx][jy][jz] * c1->getpiH0(i, j);
        }
       Pi += www_product[jx][jy][jz] * c1->getPi0();
       PiH += www_product_H[jx][jy][jz] * c1->getPiH0();
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
//    double maxT0 = max((e + p) / (1. - vx * vx - vy * vy - vz * vz) - p,//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!??????????????????...................
//                       (e + p) * (vx * vx + vy * vy + vz * vz) /
//                               (1. - vx * vx - vy * vy - vz * vz) +
//                           p);
    double maxT0 = max((e + p) * (1. + vxH_const * vxH_const + vyH_const * vyH_const + vzH_const * vz_const) + 2. * ( eH_const + pH_const ) * ( vxH_const * deltaVx +                                   vyH_const * deltaVy + vzH_const * deltaVz ) - p,
                      (e + p) * ( vxH_const * vxH_const + vyH_const * vyH_const + vzH_const * vzH_const ) + 2. * ( eH_const + pH_const ) * ( vxH_const * deltaVx + vyH_const * deltaVy + vzH_const * deltaVz ) + p );
    // double maxpi = max(fabs(pi[1][1]),fabs(pi[2][2])) ;
    double maxpi = 0.;
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++)
      if (fabs(pi[i][j]) > maxpi) maxpi = fabs(pi[i][j]);
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
      c->setpi(i, j, pi[i][j]);
      c->setpiH(i, j, piH[i][j]);
     }
    c->setPi(Pi);
    c->setPiH(PiH);
    // source term  - (tau+dt)*delta_Q_(i+1)/delta_tau
    double gamma_const = 1. + 0.5 * (vxH_const * vxH_const + vyH_const * vyH_const + vzH_const * vzH_const);
    double gamma_delta = vxH_const * deltaVx + vyH_const * deltaVy + vzH_const * deltaVz;
    double gamma = gamma_const + gamma_delta;
    double u[4];
    u[0] = gamma;
    u[1] = u[0] * vxH_const + deltaVx * gamma_const;
    u[2] = u[0] * vyH_const + deltaVy * gamma_const;
    u[3] = u[0] * vzH_const + deltaVz * gamma_const;
//    double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
//    double u[4];
//    u[0] = gamma;
//    u[1] = u[0] * vx; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!??????????????????.......................
//    u[2] = u[0] * vy;
//    u[3] = u[0] * vz;
            
    double uu_product[4][4];
    //diagonal
    uu_product[0][0] = gamma_const * (u[0] + gamma_delta);
    uu_product[1][1] = gamma_const * ( (u[0] + gamma_delta) * vxH_const * vxH_const + 2 * gamma_const * vxH_const * deltaVx );
    uu_product[2][2] = gamma_const * ( (u[0] + gamma_delta) * vyH_const * vyH_const + 2 * gamma_const * vyH_const * deltaVy );
    uu_product[3][3] = gamma_const * ( (u[0] + gamma_delta) * vzH_const * vzH_const + 2 * gamma_const * vzH_const * deltaVz );
    //off-diagonal
    uu_product[0][1] = gamma_const * ( (u[0] + gamma_delta) * vxH_const + gamma_const * deltaVx );
    uu_product[0][2] = gamma_const * ( (u[0] + gamma_delta) * vyH_const + gamma_const * deltaVy );
    uu_product[0][3] = gamma_const * ( (u[0] + gamma_delta) * vzH_const + gamma_const * deltaVz );

    uu_product[1][0] = uu_product[0][1];
    uu_product[1][2] = gamma_const * ( (u[0] + gamma_delta) * vxH_const * vyH_const + gamma_const * vxH_const * deltaVy + gamma_const * vyH_const * deltaVx );
    uu_product[1][3] = gamma_const * ( (u[0] + gamma_delta) * vxH_const * vzH_const + gamma_const * vxH_const * deltaVz + gamma_const * vzH_const * deltaVx );

    uu_product[2][0] = uu_product[0][2];
    uu_product[2][1] = uu_product[1][2];
    uu_product[2][3] = gamma_const * ( (u[0] + gamma_delta) * vyH_const * vzH_const + gamma_const * vyH_const * deltaVz + gamma_const * vzH_const * deltaVy );
    
    uu_product[3][0] = uu_product[0][3];
    uu_product[3][1] = uu_product[1][3];
    uu_product[3][2] = uu_product[2][3];
            
    double flux[4];
    for (int i = 0; i < 4; i++)
     flux[i] = -tau * (c->getpi(0, i) + c->getPi() * uu_product[0][i]);
    flux[0] += tau * c->getPi();
    c->addFlux(flux[0], flux[1], flux[2], flux[3], 0., 0., 0.);
   }  // advection loop (all cells)
}

// this procedure explicitly uses T_==0, X_==1, Y_==2, Z_==3
void Hydro::visc_flux(Cell *left, Cell *right, int direction, double eH_const, double pH_const, double vxH_const, double vyH_const, double vzH_const) {
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
 // we need to know the velocities at both cell centers at (n+1/2) in order to
 // interpolate to
 // get the value at the interface
 left->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vxl, vyl, vzl);
 right->getPrimVarHCenter(eos, tauMinusHalf, e, p, nb, nq, ns, vxr, vyr, vzr);
 vxl = 0.5 * (vxl + vxr);
 vyl = 0.5 * (vyl + vyr);
 vzl = 0.5 * (vzl + vzr);
    
    double deltaVxl = vxl - vxH_const;
    double deltaVyl = vyl - vyH_const;
    double deltaVzl = vzl - vzH_const;
    
// double v = sqrt(vxl * vxl + vyl * vyl + vzl * vzl);
 double v = vxH_const + vyH_const + vzH_const + deltaVxl + deltaVyl + deltaVzl;
// if (v > 1.) {
//  vxl = 0.99 * vxl / v;
//  vyl = 0.99 * vyl / v;
//  vzl = 0.99 * vzl / v;
// }
    if (v > 1.) {
     vxl = 0.99;
     vyl = 0.99;
     vzl = 0.99;
    }
    
 double gamma_const = 1. + 0.5 * (vxH_const * vxH_const + vyH_const * vyH_const + vzH_const * vzH_const);
 double gamma_delta = vxH_const * deltaVxl + vyH_const * deltaVyl + vzH_const * deltaVzl;
 double gamma = gamma_const + gamma_delta;
// double gamma = 1. / sqrt(1. - v * v); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// double uuu[4] = {gamma, gamma * vxl, gamma * vyl, gamma * vzl};
 double uuu[4];
 uuu[0] = gamma;
 uuu[1] = gamma * vxH_const + gamma_const * deltaVxl;
 uuu[2] = gamma * vyH_const + gamma_const * deltaVyl;
 uuu[3] = gamma * vzH_const + gamma_const * deltaVzl;
 
 double uuu_product[4][4];
 //diagonal
 uuu_product[0][0] = gamma_const * (uuu[0] + gamma_delta);
 uuu_product[1][1] = gamma_const * ( (uuu[0] + gamma_delta) * vxH_const * vxH_const + 2 * gamma_const * vxH_const * deltaVxl );
 uuu_product[2][2] = gamma_const * ( (uuu[0] + gamma_delta) * vyH_const * vyH_const + 2 * gamma_const * vyH_const * deltaVyl );
 uuu_product[3][3] = gamma_const * ( (uuu[0] + gamma_delta) * vzH_const * vzH_const + 2 * gamma_const * vzH_const * deltaVzl );
 //off-diagonal
 uuu_product[0][1] = gamma_const * ( (uuu[0] + gamma_delta) * vxH_const + gamma_const * deltaVxl );
 uuu_product[0][2] = gamma_const * ( (uuu[0] + gamma_delta) * vyH_const + gamma_const * deltaVyl );
 uuu_product[0][3] = gamma_const * ( (uuu[0] + gamma_delta) * vzH_const + gamma_const * deltaVzl );

 uuu_product[1][0] = uuu_product[0][1];
 uuu_product[1][2] = gamma_const * ( (uuu[0] + gamma_delta) * vxH_const * vyH_const + gamma_const * vxH_const * deltaVyl + gamma_const * vyH_const * deltaVxl );
 uuu_product[1][3] = gamma_const * ( (uuu[0] + gamma_delta) * vxH_const * vzH_const + gamma_const * vxH_const * deltaVzl + gamma_const * vzH_const * deltaVxl );

 uuu_product[2][0] = uuu_product[0][2];
 uuu_product[2][1] = uuu_product[1][2];
 uuu_product[2][3] = gamma_const * ( (uuu[0] + gamma_delta) * vyH_const * vzH_const + gamma_const * vyH_const * deltaVzl + gamma_const * vzH_const * deltaVyl );

 uuu_product[3][0] = uuu_product[0][3];
 uuu_product[3][1] = uuu_product[1][3];
 uuu_product[3][2] = uuu_product[2][3];
    
 double gmumu[4] = {1., -1., -1., -1.};
 if (direction == X_)
  ind2 = 1;
 else if (direction == Y_)
  ind2 = 2;
 else if (direction == Z_)
  ind2 = 3;
 for (int ind1 = 0; ind1 < 4; ind1++) {
  flux[ind1] = 0.5 * (left->getpiH(ind1, ind2) + right->getpiH(ind1, ind2));
  if (ind1 == ind2)
   flux[ind1] += -0.5 * (left->getPiH() + right->getPiH()) *
                 gmumu[ind1];  // gmunu is diagonal
//  flux[ind1] +=
//      0.5 * (left->getPiH() + right->getPiH()) * uuu[ind1] * uuu[ind2];
    flux[ind1] +=
        0.5 * (left->getPiH() + right->getPiH()) * uuu_product[ind1][ind2];
 }
 for (int i = 0; i < 4; i++) flux[i] = flux[i] * tauMinusHalf * dt / dxa;
 left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
 right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
}

void Hydro::performStep(void) {
 // debugRiemann = false ; // turn off debug output

 f->updateM(tau, dt);

 tau_z = dt / 2. / log(1 + dt / 2. / tau);

 double e_const=0;
 double p_const=0;
 double vx_const=0;
 double vy_const=0;
 double vz_const=0;

 double eH_const=0;
 double pH_const=0;
 double vxH_const=0;
 double vyH_const=0;
 double vzH_const=0;
    
 double ePrev_const=0;
 double pPrev_const=0;
 double vxPrev_const=0;
 double vyPrev_const=0;
 double vzPrev_const=0;
    double vx_prev = 0;
    double p_prev = 0;
    
 //getting averaged value
 for (int iy = 0; iy < f->getNY(); iy++)
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++) {
       
   double e, p, nb, nq, ns, vx, vy, vz, vuv;
   Cell *c = f->getCell(ix, iy, iz);
       
   c -> getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
   e_const = e_const + e/((f->getNY())*(f->getNX())*(f->getNZ()));
   p_const = p_const + p/((f->getNY())*(f->getNX())*(f->getNZ()));
   vx_const = vx_const + vx/((f->getNY())*(f->getNX())*(f->getNZ()));
   vy_const = vy_const + vy/((f->getNY())*(f->getNX())*(f->getNZ()));
   vz_const = vz_const + vz/((f->getNY())*(f->getNX())*(f->getNZ()));
       
   c -> getPrimVarHCenter(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
   eH_const = eH_const + e/((f->getNY())*(f->getNX())*(f->getNZ()));
   pH_const = pH_const + p/((f->getNY())*(f->getNX())*(f->getNZ()));
   vxH_const = vxH_const + vx/((f->getNY())*(f->getNX())*(f->getNZ()));
   vyH_const = vyH_const + vy/((f->getNY())*(f->getNX())*(f->getNZ()));
   vzH_const = vzH_const + vz/((f->getNY())*(f->getNX())*(f->getNZ()));
       
   c -> getPrimVarPrev(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
   vx_prev = vx;
       p_prev = p;
   ePrev_const = ePrev_const + e/((f->getNY())*(f->getNX())*(f->getNZ()));
   pPrev_const = pPrev_const + p/((f->getNY())*(f->getNX())*(f->getNZ()));
   vxPrev_const = vxPrev_const + vx/((f->getNY())*(f->getNX())*(f->getNZ()));
   vyPrev_const = vyPrev_const + vy/((f->getNY())*(f->getNX())*(f->getNZ()));
   vzPrev_const = vzPrev_const + vz/((f->getNY())*(f->getNX())*(f->getNZ()));
   //cout << vx << endl;
   }
    

//    // vypis profilu energie
    for (int ix = 0; ix < f->getNX(); ix++) {
        //          for (int iy = 0; iy < f->getNY(); iy++){
        
        double e, p, nb, nq, ns, vx, vy, vz, cs, T, mub, muq, mus;
        Cell *c = f->getCell(ix, 0, 0);
        
        c -> getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
        cs = eos->cs();
        eos->eos(e, 0., 0., 0., T, mub, mus, muq, p);
        
        if(ix<=40){
            cout.precision(6);
        }
        else{
            cout.precision(7);
        }
        
        //              if(ix==20 || ix ==19 || ix ==21 || ix ==60 || ix ==61 || ix ==59){
//        if(ix==39 || ix== 40 || ix==41){
                  cout << vx << "   " << f->getX(ix) << "      " << e << "       " << p << "     " << T << endl;
//              }

          
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
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix + 1, iy, iz), X_, PREDICT, e_const, p_const, vx_const, vy_const, vz_const, ix);
   }
 //	cout << "predictor X done\n" ;
 // Y dir
 for (int iz = 0; iz < f->getNZ(); iz++)
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, PREDICT, e_const, p_const, vx_const, vy_const, vz_const, ix);
   }
 //	cout << "predictor Y done\n" ;
 // Z dir
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz + 1), Z_, PREDICT, e_const, p_const, vx_const, vy_const, vz_const, ix);
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
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix + 1, iy, iz), X_, CORRECT, e_const, p_const, vx_const, vy_const, vz_const, ix);
   }
 //	cout << "corrector X done\n" ;
 // Y dir
 for (int iz = 0; iz < f->getNZ(); iz++)
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, CORRECT, e_const, p_const, vx_const, vy_const, vz_const, ix);
   }
 //	cout << "corrector Y done\n" ;
 // Z dir
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz + 1), Z_, CORRECT, e_const, p_const, vx_const, vy_const, vz_const, ix);
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
  ISformal(e_const, p_const, vx_const, vy_const, vz_const, eH_const, pH_const, vxH_const, vyH_const, vzH_const, ePrev_const, pPrev_const, vxPrev_const, vyPrev_const, vzPrev_const);  // evolution of viscous quantities according to IS equations
  // X dir
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++)
    for (int ix = 0; ix < f->getNX(); ix++) {
     visc_flux(f->getCell(ix, iy, iz), f->getCell(ix + 1, iy, iz), X_, eH_const, pH_const, vxH_const, vyH_const, vzH_const);
    }
  //	cout << "visc_flux X done\n" ;
  // Y dir
  for (int iz = 0; iz < f->getNZ(); iz++)
   for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++) {
     visc_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy + 1, iz), Y_, eH_const, pH_const, vxH_const, vyH_const, vzH_const);
    }
  //	cout << "visc_flux Y done\n" ;
  // Z dir
  for (int ix = 0; ix < f->getNX(); ix++)
   for (int iy = 0; iy < f->getNY(); iy++)
    for (int iz = 0; iz < f->getNZ(); iz++) {
     visc_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz + 1), Z_, eH_const, pH_const, vxH_const, vyH_const, vzH_const);
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
