#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "rmn.h"
#include "cll.h"

using namespace std;

// slope limiter; chooses minimal abs of the neighbouring slopes

double minmod(double a, double b) {
 if (a * b <= 0.) return 0.;
 //	else return (a*a*b+a*b*b)/(a*a+b*b) ;
 if (fabs(a) > fabs(b))
  return b;
 else
  return a;
}

// index44: returns an index of pi^{mu nu} mu,nu component in a plain 1D array
int index44(const int &i, const int &j) {
 if (i > 3 || j > 3 || i < 0 || j < 0) {
  std::cout << "index44: i j " << i << " " << j << endl;
  exit(1);
 }
 if (j < i)
  return (i * (i + 1)) / 2 + j;
 else
  return (j * (j + 1)) / 2 + i;
}

Cell::Cell() {
 for (int i = 0; i < 7; i++) {
  Q[i] = 0.;
  Qh[i] = 0.;
  Qprev[i] = 0.;
  Q0[i] = 0.;
  Qh0[i] = 0.;
  Q0prev[i] = 0.;
  flux[i] = 0.;
 }
 viscCorrCut = 1.;
 for (int i = 0; i < 10; i++) {
  pi[i] = 0.0;
  piH[i] = 0.0;
  pi0[i] = 0.0;
  piH0[i] = 0.0;
  pi_bck[i] = 0.0;
  piH_bck[i] = 0.0;
  pi0_bck[i] = 0.0;
  piH0_bck[i] = 0.0;
  pi_bck_prev[i] = 0.0;
  piH_bck_prev[i] = 0.0;
  pi0_bck_prev[i] = 0.0;
  piH0_bck_prev[i] = 0.0;
 }
 Pi = 0.0;
 PiH = 0.0;
 Pi0 = 0.0;
 PiH0 = 0.0;
 Pi_bck = 0.0;
 PiH_bck = 0.0;
 Pi0_bck = 0.0;
 PiH0_bck = 0.0;
 Pi_bck_prev = 0.0;
 PiH_bck_prev = 0.0;
 Pi0_bck_prev = 0.0;
 PiH0_bck_prev = 0.0;
 setAllM(0.);
}

void Cell::importVars(Cell* c) {
 c->getQ(Q);
 c->getQh(Qh);
 c->getQprev(Qprev); // not necessary but for consistency
 Pi = c->getPi();
 PiH = c->getPiH();
 for(int i=0; i<4; i++)
 for(int j=0; j<4; j++) {
  pi[index44(i, j)] = c->getpi(i, j);
  piH[index44(i, j)] = c->getpiH(i, j); // not necessary either
 }
 for(int i=X_; i<=Z_; i++)
  m[i-1] = c->getM(i);
 viscCorrCut = c->getViscCorrCutFlag();
}

void Cell::updateByFlux() {
 if(Q0[0]+Q[0]+flux[0]<0.) // modified condition
  return;
 for (int i = 0; i < 7; i++) Q[i] += flux[i];
}

void Cell::updateByViscFlux() {
// if(fabs(flux[0]) <= fabs(0.5*Q[0])) { // modified condition - an absolute value was added
  for (int i = 0; i < 7; i++) Q[i] += flux[i];
// } else if (flux[0]!=0.){
//  double fac;
//  fac = fabs(0.5*Q[0]/flux[0]);
//  for (int i = 0; i < 7; i++) Q[i] += fac*flux[i];
// }
}

void Cell::updateQtoQhByFlux() {
 for (int i = 0; i < 7; i++) Qh[i] = Q[i] + flux[i];
}

void Cell::getPrimVar(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                      double &_nq, double &_ns, double &_vx, double &_vy,
                      double &_vz, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0) {
    double _Q[7];
//    double _Q0[7];
    for (int i = 0; i < 7; i++) {
        _Q[i] = Q[i] / tau;
//        _Q0[i] = Q0[i] / tau;
    }
// transformPVQ0(eos, _Q0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0);
 transformPV(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0 );
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVar:\n";
  Dump(tau);
 }
 #endif
 //------------------------------------------
}
void Cell::getPrimVarQ0(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                      double &_nq, double &_ns, double &_vx, double &_vy,
                      double &_vz) {
 double _Q0[7];
 for (int i = 0; i < 7; i++) _Q0[i] = Q0[i] / tau;
 transformPVQ0(eos, _Q0, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVar:\n";
  Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarLeft(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz, int dir, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0) {
    double Qr[7], Ql[7], dQ[7];
//    double Q0r[7], Ql0[7], dQ0[7];

    next[dir - 1]->getQ(Qr);
    prev[dir - 1]->getQ(Ql);
//    next[dir - 1]->getQ0(Q0r);
//    prev[dir - 1]->getQ0(Ql0);

    for (int i = 0; i < 7; i++){
        dQ[i] = minmod((Qr[i] - Q[i]) / 2., (Q[i] - Ql[i]) / 2.);
//        dQ0[i] = minmod((Q0r[i] - Q0[i]) / 2., (Q0[i] - Ql0[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
        Ql[i] = (Q[i] - dQ[i]) / tau;
//        Ql0[i] = (Q0[i] - dQ0[i]) / tau;
    }
// transformPVQ0(eos, Ql0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0 );
 transformPV(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0  );
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarLeft:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}
void Cell::getPrimVarLeftQ0(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz, int dir) {
 double Qr0[7], Ql0[7], dQ0[7];

 next[dir - 1]->getQ0(Qr0);
 prev[dir - 1]->getQ0(Ql0);

 for (int i = 0; i < 7; i++)
  dQ0[i] = minmod((Qr0[i] - Q0[i]) / 2., (Q0[i] - Ql0[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql0[i] = (Q0[i] - dQ0[i]) / tau;
 transformPVQ0(eos, Ql0, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarLeft:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarRight(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz, int dir, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0) {
    double Qr[7], Ql[7], dQ[7];
//    double Q0r[7], Q0l[7], dQ0[7];

    next[dir - 1]->getQ(Qr);
    prev[dir - 1]->getQ(Ql);
//    next[dir - 1]->getQ0(Q0r);
//    prev[dir - 1]->getQ0(Q0l);

    for (int i = 0; i < 7; i++){
        dQ[i] = minmod((Qr[i] - Q[i]) / 2., (Q[i] - Ql[i]) / 2.);
//        dQ0[i] = minmod((Q0r[i] - Q0[i]) / 2., (Q0[i] - Q0l[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
        Qr[i] = (Q[i] + dQ[i]) / tau;
//        Q0r[i] = (Q0[i] + dQ0[i]) / tau;
    }
//    transformPVQ0(eos, Q0r, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0 );
    transformPV(eos, Qr, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0  );
//-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarRight:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}
void Cell::getPrimVarRightQ0(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz, int dir) {
 double Qr0[7], Ql0[7], dQ0[7];

 next[dir - 1]->getQ0(Qr0);
 prev[dir - 1]->getQ0(Ql0);

 for (int i = 0; i < 7; i++)
  dQ0[i] = minmod((Qr0[i] - Q0[i]) / 2., (Q0[i] - Ql0[i]) / 2.);

 for (int i = 0; i < 7; i++) Qr0[i] = (Q0[i] + dQ0[i]) / tau;
 transformPVQ0(eos, Qr0, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarRight:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarHLeft(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz, int dir, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0) {
    double Qr[7], Ql[7], dQ[7];
//    double Q0r[7], Q0l[7], dQ0[7];

    next[dir - 1]->getQh(Qr);
    prev[dir - 1]->getQh(Ql);
//    next[dir - 1]->getQh0(Q0r);
//    prev[dir - 1]->getQh0(Q0l);
//    next[dir - 1]->getQ0(Q0r);
//    prev[dir - 1]->getQ0(Q0l);

    for (int i = 0; i < 7; i++){
        dQ[i] = minmod((Qr[i] - Qh[i]) / 2., (Qh[i] - Ql[i]) / 2.);
//        dQ0[i] = minmod((Q0r[i] - Qh0[i]) / 2., (Qh0[i] - Q0l[i]) / 2.);
//        dQ0[i] = minmod((Q0r[i] - Q0[i]) / 2., (Q0[i] - Q0l[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
        Ql[i] = (Qh[i] - dQ[i]) / tau;
//        Q0l[i] = (Qh0[i] - dQ0[i]) / tau;
//        Q0l[i] = (Q0[i] - dQ0[i]) / tau;
    }
//    transformPVQ0(eos, Q0l, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0 );
    transformPV(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0  ); //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarHLeft:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarHLeftQ0(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz, int dir) {
 double Q0r[7], Q0l[7], dQ0[7];

//    next[dir - 1]->getQh0(Q0r);
//    prev[dir - 1]->getQh0(Q0l);
    next[dir - 1]->getQ0(Q0r);
    prev[dir - 1]->getQ0(Q0l);

    for (int i = 0; i < 7; i++){
//        dQ0[i] = minmod((Q0r[i] - Qh0[i]) / 2., (Qh0[i] - Q0l[i]) / 2.);
        dQ0[i] = minmod((Q0r[i] - Q0[i]) / 2., (Q0[i] - Q0l[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
//        Q0l[i] = (Qh0[i] - dQ0[i]) / tau;
        Q0l[i] = (Q0[i] - dQ0[i]) / tau;
    }
    transformPVQ0(eos, Q0l, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz );
//-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarHLeft:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarHRight(EoS *eos, double tau, double &_e, double &_p,
                            double &_nb, double &_nq, double &_ns, double &_vx,
                            double &_vy, double &_vz, int dir, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0) {
    double Qr[7], Ql[7], dQ[7];
//    double Q0r[7], Q0l[7], dQ0[7];

    next[dir - 1]->getQh(Qr);
    prev[dir - 1]->getQh(Ql);
//    next[dir - 1]->getQh0(Q0r);
//    prev[dir - 1]->getQh0(Q0l);
//    next[dir - 1]->getQ0(Q0r);
//    prev[dir - 1]->getQ0(Q0l);

    for (int i = 0; i < 7; i++){
        dQ[i] = minmod((Qr[i] - Qh[i]) / 2., (Qh[i] - Ql[i]) / 2.);
//        dQ0[i] = minmod((Q0r[i] - Qh0[i]) / 2., (Qh0[i] - Q0l[i]) / 2.);
//        dQ0[i] = minmod((Q0r[i] - Q0[i]) / 2., (Q0[i] - Q0l[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
        Qr[i] = (Qh[i] + dQ[i]) / tau;
//        Q0r[i] = (Qh0[i] + dQ0[i]) / tau;
//        Q0r[i] = (Q0[i] + dQ0[i]) / tau;
    }
//    transformPVQ0(eos, Q0r, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0 );
    transformPV(eos, Qr, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0  );
//-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarHRight:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarHRightQ0(EoS *eos, double tau, double &_e, double &_p,
                            double &_nb, double &_nq, double &_ns, double &_vx,
                            double &_vy, double &_vz, int dir) {
 double Q0r[7], Q0l[7], dQ0[7];

//    next[dir - 1]->getQh0(Q0r);
//    prev[dir - 1]->getQh0(Q0l);
    next[dir - 1]->getQ0(Q0r);
    prev[dir - 1]->getQ0(Q0l);

    for (int i = 0; i < 7; i++){
//        dQ0[i] = minmod((Q0r[i] - Qh0[i]) / 2., (Qh0[i] - Q0l[i]) / 2.);
        dQ0[i] = minmod((Q0r[i] - Q0[i]) / 2., (Q0[i] - Q0l[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
//        Q0r[i] = (Qh0[i] + dQ0[i]) / tau;
        Q0r[i] = (Q0[i] + dQ0[i]) / tau;
    }
    transformPVQ0(eos, Q0r, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz );
//-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarHRight:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}


void Cell::getPrimVarHCenter(EoS *eos, double tau, double &_e, double &_p,
                             double &_nb, double &_nq, double &_ns, double &_vx,
                             double &_vy, double &_vz, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0) {
    double _Q[7];
//    double _Q0[7];
    for (int i = 0; i < 7; i++){
        _Q[i] = Qh[i] / tau;
//        _Q0[i] = Qh0[i] / tau;
//        _Q0[i] = Q0[i] / tau;
    }
//    transformPVQ0(eos, _Q0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0 );
    transformPV(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0  );
}

void Cell::getPrimVarHCenterQ0(EoS *eos, double tau, double &_e, double &_p,
                             double &_nb, double &_nq, double &_ns, double &_vx,
                             double &_vy, double &_vz) {
 double _Q0[7];
    for (int i = 0; i < 7; i++){
//        _Q0[i] = Qh0[i] / tau; -toto !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        _Q0[i] = Q0[i] / tau;
    }
    transformPVQ0(eos, _Q0, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz );
}

void Cell::getPrimVarPrev(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0) {
    double _Q[7];
//    double _Q0[7];
    for (int i = 0; i < 7; i++){
        _Q[i] = Qprev[i] / tau;
//        _Q0[i] = Q0prev[i] / tau;
//        _Q0[i] = Q0[i] / tau;
    }
//    transformPVQ0(eos, _Q0, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0 );
    transformPV(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz, e_0, p_0, nb_0, nq_0, ns_0, vx_0, vy_0, vz_0  );
}

void Cell::getPrimVarPrevQ0(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz) {
 double _Q0[7];
    for (int i = 0; i < 7; i++){
//        _Q0[i] = Q0prev[i] / tau; -toto !!!!!!!!!!!!!!!!!!!!!!
        _Q0[i] = Q0[i] / tau;
    }
    transformPVQ0(eos, _Q0, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz );
}
// setting of fluctuating variables
void Cell::setPrimVar(EoS *eos, double tau, double _e, double _nb, double _nq,
                      double _ns, double _vx, double _vy, double _vz, double e0, double vx0, double vy0, double vz0) {
 const double gamma0 = 1. / sqrt(1 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
 double p0 = eos->p(e0, _nb, _nq, _ns); // background pressure
 const double p = eos->p(_e+e0, _nb, _nq, _ns) - eos->p(e0, _nb, _nq, _ns); // fluctuation of pressure - will this work even for non-linear EOS?
 Q[T_] = tau * ( (_e + p) * gamma0 * gamma0 - p );
 Q[X_] = tau * ( (e0 + p0) * gamma0 * _vx + (_e + p) * gamma0 * gamma0 * vx0 );
 Q[Y_] = tau * ( (e0 + p0) * gamma0 * _vy + (_e + p) * gamma0 * gamma0 * vy0 );
 Q[Z_] = tau * ( (e0 + p0) * gamma0 * _vz + (_e + p) * gamma0 * gamma0 * vz0 );
 Q[NB_] = tau * _nb * gamma0; // not linearized
 Q[NQ_] = tau * _nq * gamma0; // not linearized
 Q[NS_] = tau * _ns * gamma0; // not linearized
 if (std::isinf(Q[NB_]) or std::isnan(Q[NB_])) {
  cout << "init error!\n";
  eos->p(_e, _nb, _nq, _ns);
  cout << "e = " << _e << " p = " << p << " vx = " << _vx << " vy = " << _vy
       << " vz = " << _vz << endl;
  //		exit(1) ;
  return;
 }
}

// function for setting of background variables
void Cell::setPrimVarQ0(EoS *eos, double tau, double _e, double _nb, double _nq,
                      double _ns, double _vx, double _vy, double _vz) {
 const double gamma2 = 1. / (1 - _vx * _vx - _vy * _vy - _vz * _vz);
 const double p = eos->p(_e, _nb, _nq, _ns);
 Q0[T_] = tau * (_e + p * (_vx * _vx + _vy * _vy + _vz * _vz)) * gamma2;
 Q0[X_] = tau * (_e + p) * _vx * gamma2;
 Q0[Y_] = tau * (_e + p) * _vy * gamma2;
 Q0[Z_] = tau * (_e + p) * _vz * gamma2;
 Q0[NB_] = tau * _nb * sqrt(gamma2);
 Q0[NQ_] = tau * _nq * sqrt(gamma2);
 Q0[NS_] = tau * _ns * sqrt(gamma2);
 if (std::isinf(Q[NB_]) or std::isnan(Q[NB_])) {
  cout << "init error!\n";
  eos->p(_e, _nb, _nq, _ns);
  cout << "e = " << _e << " p = " << p << " vx = " << _vx << " vy = " << _vy
       << " vz = " << _vz << " gamma2 = " << gamma2 << endl;
  //        exit(1) ;
  return;
 }
}

void Cell::Dump(double tau) {
 cout << "---------cell values dump-------\n";
 cout << setw(5) << ix << setw(5) << iy << setw(5) << iz << endl;
 cout << setw(14) << Q[0] / tau << setw(14) << Q[1] / tau << setw(14)
      << Q[2] / tau << setw(14) << Q[3] / tau << endl;
 cout << setw(14) << Q[4] / tau << setw(14) << Q[5] / tau << setw(14)
      << Q[6] / tau << endl;
 cout << setw(14) << Qh[0] / tau << setw(14) << Qh[1] / tau << setw(14)
      << Qh[2] / tau << setw(14) << Qh[3] / tau << endl;
 cout << setw(14) << Qh[4] / tau << setw(14) << Qh[5] / tau << setw(14)
      << Qh[6] / tau << endl;

 cout << "--------------------------------\n";
}

//void d_pi(int ix, int iy, int iz, double &dpi0[4][4][4]){
//#ifdef CARTESIAN
//    double tauPlusHalf = 1.0;
//#else
//    double tauPlusHalf = tau + 0.5 * dt;
//#endif
//    double dx = f->getDx(), dy = f->getDy(), dz = f->getDz();
//
//    for(int i = 0; i < 4; i++)
//        for(int j = 0; j < 4; j++){
//            dpi0[i][j][0] = (f->getCell(ix, iy, iz)->getpi_bck(i,j) - f->getCell(ix, iy, iz)->getpi_bck_prev(i,j)) / dt;
//            dpi0[i][j][1] = 0.5 * (f->getCell(ix + 1, iy, iz)->getpi_bck(i,j) - f->getCell(ix - 1, iy, iz)->getpi_bck(i,j)) / dx;
//            dpi0[i][j][2] = 0.5 * (f->getCell(ix, iy + 1, iz)->getpi_bck(i,j) - f->getCell(ix, iy - 1, iz)->getpi_bck(i,j)) / dy;
//            dpi0[i][j][3] = 0.5 * (f->getCell(ix, iy, iz + 1)->getpi_bck(i,j) - f->getCell(ix, iy, iz - 1)->getpi_bck(i,j)) / dz / tauPlusHalf;
//        }
////        } else {  // matter-vacuum
////            dpi0[i][j][1] = 0.;
////            dpi0[i][j][2] = 0.;
////            dpi0[i][j][3] = 0.;
////        }
////    if (e1+e_01 <= 0. || e0+e_0 <= 0.) {  // matter-vacuum
////        dpi0[i][j][0] = 0.;
////    }
//
//        // additional terms from Christoffel symbols :) ??????????????????????
//#ifndef CARTESIAN
////        dmu[3][0] += uuu[3] / (tau - 0.5 * dt);
////        dmu[3][3] += uuu[0] / (tau - 0.5 * dt);
////    dmu0[3][0] += uuu0[3] / (tau - 0.5 * dt);
////    dmu0[3][3] += uuu0[0] / (tau - 0.5 * dt);
//#endif
//}
