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
 //    else return (a*a*b+a*b*b)/(a*a+b*b) ;
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
  d_Q[i] = 0.;
  d_Qh[i] = 0.;
  d_Qprev[i] = 0.;
  Q_bck[i] = 0.;
  Qh_bck[i] = 0.;
  Qprev_bck[i] = 0.;
  d_flux[i] = 0.;
 }
 viscCorrCut = 1.;
 for (int i = 0; i < 10; i++) {
  d_pi[i] = 0.0;
  d_piH[i] = 0.0;
  d_pi0[i] = 0.0;
  d_piH0[i] = 0.0;
  pi_bck[i] = 0.0;
  piH_bck[i] = 0.0;
  pi0_bck[i] = 0.0;
  piH0_bck[i] = 0.0;
  pi_bck_prev[i] = 0.0;
  piH_bck_prev[i] = 0.0;
  pi0_bck_prev[i] = 0.0;
  piH0_bck_prev[i] = 0.0;
 }
 d_Pi = 0.0;
 d_PiH = 0.0;
 d_Pi0 = 0.0;
 d_PiH0 = 0.0;
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
 c->getQ(d_Q);
 c->getQh(d_Qh);
 c->getQprev(d_Qprev); // not necessary but for consistency
 d_Pi = c->getPi();
 d_PiH = c->getPiH();
 for(int i=0; i<4; i++)
 for(int j=0; j<4; j++) {
  d_pi[index44(i, j)] = c->getpi(i, j);
  d_piH[index44(i, j)] = c->getpiH(i, j); // not necessary either
 }
 for(int i=X_; i<=Z_; i++)
  m[i-1] = c->getM(i);
 viscCorrCut = c->getViscCorrCutFlag();
}

void Cell::updateByFlux() {
 if(Q_bck[0]+d_Q[0]+d_flux[0]<0.) // modified condition
  return;
 for (int i = 0; i < 7; i++) d_Q[i] += d_flux[i];
}

void Cell::updateByViscFlux() {
// if(fabs(flux[0]) <= fabs(0.5*Q[0])) { // modified condition - an absolute value was added
  for (int i = 0; i < 7; i++) d_Q[i] += d_flux[i];
// } else if (flux[0]!=0.){
//  double fac;
//  fac = fabs(0.5*Q[0]/flux[0]);
//  for (int i = 0; i < 7; i++) Q[i] += fac*flux[i];
// }
}

void Cell::updateQtoQhByFlux() {
 for (int i = 0; i < 7; i++) d_Qh[i] = d_Q[i] + d_flux[i];
}

void Cell::getPrimVar(EoS *eos, double tau, double &_d_e, double &_d_p, double &_d_nb,
                      double &_d_nq, double &_d_ns, double &_d_vx, double &_d_vy,
                      double &_d_vz, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck) {
    double _d_Q[7];
//    double _Q_bck[7];
    for (int i = 0; i < 7; i++) {
        _d_Q[i] = d_Q[i] / tau;
//        _Q_bck[i] = Q_bck[i] / tau;
    }
// transformPVQbck(eos, _Qbck, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck);
 transformPV(eos, _d_Q, _d_e, _d_p, _d_nb, _d_nq, _d_ns, _d_vx, _d_vy, _d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck );
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVar:\n";
  Dump(tau);
 }
 #endif
 //------------------------------------------
}
void Cell::getPrimVarQbck(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                      double &_nq, double &_ns, double &_vx, double &_vy,
                      double &_vz) {
 double _Q_bck[7];
 for (int i = 0; i < 7; i++) _Q_bck[i] = Q_bck[i] / tau;
 transformPVQbck(eos, _Q_bck, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVar:\n";
  Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarLeft(EoS *eos, double tau, double &_d_e, double &_d_p,
                          double &_d_nb, double &_d_nq, double &_d_ns, double &_d_vx,
                          double &_d_vy, double &_d_vz, int dir, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck) {
    double d_Qr[7], d_Ql[7], d_dQ[7];
//    double Qr_bck[7], Ql0[7], dQ_bck[7];

    next[dir - 1]->getQ(d_Qr);
    prev[dir - 1]->getQ(d_Ql);
//    next[dir - 1]->getQbck(Qr_bck);
//    prev[dir - 1]->getQbck(Ql0);

    for (int i = 0; i < 7; i++){
        d_dQ[i] = minmod((d_Qr[i] - d_Q[i]) / 2., (d_Q[i] - d_Ql[i]) / 2.);
//        dQ_bck[i] = minmod((Qr_bck[i] - Q_bck[i]) / 2., (Q_bck[i] - Ql0[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
        d_Ql[i] = (d_Q[i] - d_dQ[i]) / tau;
//        Ql0[i] = (Q_bck[i] - dQ_bck[i]) / tau;
    }
// transformPVQbck(eos, Ql0, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck );
 transformPV(eos, d_Ql, _d_e, _d_p, _d_nb, _d_nq, _d_ns, _d_vx, _d_vy, _d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck  );
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
void Cell::getPrimVarLeftQbck(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz, int dir) {
 double Qr_bck[7], Ql_bck[7], dQ_bck[7];

 next[dir - 1]->getQbck(Qr_bck);
 prev[dir - 1]->getQbck(Ql_bck);

 for (int i = 0; i < 7; i++)
  dQ_bck[i] = minmod((Qr_bck[i] - Q_bck[i]) / 2., (Q_bck[i] - Ql_bck[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql_bck[i] = (Q_bck[i] - dQ_bck[i]) / tau;
 transformPVQbck(eos, Ql_bck, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
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

void Cell::getPrimVarRight(EoS *eos, double tau, double &_d_e, double &_d_p,
                           double &_d_nb, double &_d_nq, double &_d_ns, double &_d_vx,
                           double &_d_vy, double &_d_vz, int dir, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck) {
    double d_Qr[7], d_Ql[7], d_dQ[7];
//    double Qr_bck[7], Ql_bck[7], dQ_bck[7];

    next[dir - 1]->getQ(d_Qr);
    prev[dir - 1]->getQ(d_Ql);
//    next[dir - 1]->getQbck(Qr_bck);
//    prev[dir - 1]->getQbck(Ql_bck);

    for (int i = 0; i < 7; i++){
        d_dQ[i] = minmod((d_Qr[i] - d_Q[i]) / 2., (d_Q[i] - d_Ql[i]) / 2.);
//        dQ_bck[i] = minmod((Qr_bck[i] - Q_bck[i]) / 2., (Q_bck[i] - Ql_bck[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
        d_Qr[i] = (d_Q[i] + d_dQ[i]) / tau;
//        Qr_bck[i] = (Q_bck[i] + dQ_bck[i]) / tau;
    }
//    transformPVQbck(eos, Qr_bck, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck );
    transformPV(eos, d_Qr, _d_e, _d_p, _d_nb, _d_nq, _d_ns, _d_vx, _d_vy, _d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck  );
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
void Cell::getPrimVarRightQbck(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz, int dir) {
 double Qr_bck[7], Ql_bck[7], dQ_bck[7];

 next[dir - 1]->getQbck(Qr_bck);
 prev[dir - 1]->getQbck(Ql_bck);

 for (int i = 0; i < 7; i++)
  dQ_bck[i] = minmod((Qr_bck[i] - Q_bck[i]) / 2., (Q_bck[i] - Ql_bck[i]) / 2.);

 for (int i = 0; i < 7; i++) Qr_bck[i] = (Q_bck[i] + dQ_bck[i]) / tau;
 transformPVQbck(eos, Qr_bck, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
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

void Cell::getPrimVarHLeft(EoS *eos, double tau, double &_d_e, double &_d_p,
                           double &_d_nb, double &_d_nq, double &_d_ns, double &_d_vx,
                           double &_d_vy, double &_d_vz, int dir, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck) {
    double d_Qr[7], d_Ql[7], d_dQ[7];
//    double Qr_bck[7], Ql_bck[7], dQ_bck[7];

    next[dir - 1]->getQh(d_Qr);
    prev[dir - 1]->getQh(d_Ql);
//    next[dir - 1]->getQh0(Qr_bck);
//    prev[dir - 1]->getQh0(Ql_bck);
//    next[dir - 1]->getQbck(Qr_bck);
//    prev[dir - 1]->getQbck(Ql_bck);

    for (int i = 0; i < 7; i++){
        d_dQ[i] = minmod((d_Qr[i] - d_Qh[i]) / 2., (d_Qh[i] - d_Ql[i]) / 2.);
//        dQ_bck[i] = minmod((Qr_bck[i] - Qh0[i]) / 2., (Qh0[i] - Ql_bck[i]) / 2.);
//        dQ_bck[i] = minmod((Qr_bck[i] - Q_bck[i]) / 2., (Q_bck[i] - Ql_bck[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
        d_Ql[i] = (d_Qh[i] - d_dQ[i]) / tau;
//        Ql_bck[i] = (Qh0[i] - dQ_bck[i]) / tau;
//        Ql_bck[i] = (Q_bck[i] - dQ_bck[i]) / tau;
    }
//    transformPVQbck(eos, Ql_bck, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck );
    transformPV(eos, d_Ql, _d_e, _d_p, _d_nb, _d_nq, _d_ns, _d_vx, _d_vy, _d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck  ); //-------------------- debug ---------------
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

void Cell::getPrimVarHLeftQbck(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz, int dir) {
 double Qr_bck[7], Ql_bck[7], dQ_bck[7];

//    next[dir - 1]->getQh0(Qr_bck);
//    prev[dir - 1]->getQh0(Ql_bck);
    next[dir - 1]->getQbck(Qr_bck);
    prev[dir - 1]->getQbck(Ql_bck);

    for (int i = 0; i < 7; i++){
//        dQ_bck[i] = minmod((Qr_bck[i] - Qh0[i]) / 2., (Qh0[i] - Ql_bck[i]) / 2.);
        dQ_bck[i] = minmod((Qr_bck[i] - Q_bck[i]) / 2., (Q_bck[i] - Ql_bck[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
//        Ql_bck[i] = (Qh0[i] - dQ_bck[i]) / tau;
        Ql_bck[i] = (Q_bck[i] - dQ_bck[i]) / tau;
    }
    transformPVQbck(eos, Ql_bck, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz );
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

void Cell::getPrimVarHRight(EoS *eos, double tau, double &_d_e, double &_d_p,
                            double &_d_nb, double &_d_nq, double &_d_ns, double &_d_vx,
                            double &_d_vy, double &_d_vz, int dir, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck) {
    double d_Qr[7], d_Ql[7], d_dQ[7];
//    double Qr_bck[7], Ql_bck[7], dQ_bck[7];

    next[dir - 1]->getQh(d_Qr);
    prev[dir - 1]->getQh(d_Ql);
//    next[dir - 1]->getQh0(Qr_bck);
//    prev[dir - 1]->getQh0(Ql_bck);
//    next[dir - 1]->getQbck(Qr_bck);
//    prev[dir - 1]->getQbck(Ql_bck);

    for (int i = 0; i < 7; i++){
        d_dQ[i] = minmod((d_Qr[i] - d_Qh[i]) / 2., (d_Qh[i] - d_Ql[i]) / 2.);
//        dQ_bck[i] = minmod((Qr_bck[i] - Qh0[i]) / 2., (Qh0[i] - Ql_bck[i]) / 2.);
//        dQ_bck[i] = minmod((Qr_bck[i] - Q_bck[i]) / 2., (Q_bck[i] - Ql_bck[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
        d_Qr[i] = (d_Qh[i] + d_dQ[i]) / tau;
//        Qr_bck[i] = (Qh0[i] + dQ_bck[i]) / tau;
//        Qr_bck[i] = (Q_bck[i] + dQ_bck[i]) / tau;
    }
//    transformPVQbck(eos, Qr_bck, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck );
    transformPV(eos, d_Qr, _d_e, _d_p, _d_nb, _d_nq, _d_ns, _d_vx, _d_vy, _d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck  );
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

void Cell::getPrimVarHRightQbck(EoS *eos, double tau, double &_e, double &_p,
                            double &_nb, double &_nq, double &_ns, double &_vx,
                            double &_vy, double &_vz, int dir) {
 double Qr_bck[7], Ql_bck[7], dQ_bck[7];

//    next[dir - 1]->getQh0(Qr_bck);
//    prev[dir - 1]->getQh0(Ql_bck);
    next[dir - 1]->getQbck(Qr_bck);
    prev[dir - 1]->getQbck(Ql_bck);

    for (int i = 0; i < 7; i++){
//        dQ_bck[i] = minmod((Qr_bck[i] - Qh0[i]) / 2., (Qh0[i] - Ql_bck[i]) / 2.);
        dQ_bck[i] = minmod((Qr_bck[i] - Q_bck[i]) / 2., (Q_bck[i] - Ql_bck[i]) / 2.);
    }

    for (int i = 0; i < 7; i++){
//        Qr_bck[i] = (Qh0[i] + dQ_bck[i]) / tau;
        Qr_bck[i] = (Q_bck[i] + dQ_bck[i]) / tau;
    }
    transformPVQbck(eos, Qr_bck, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz );
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


void Cell::getPrimVarHCenter(EoS *eos, double tau, double &_d_e, double &_d_p,
                             double &_d_nb, double &_d_nq, double &_d_ns, double &_d_vx,
                             double &_d_vy, double &_d_vz, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck) {
    double _d_Q[7];
//    double _Q_bck[7];
    for (int i = 0; i < 7; i++){
        _d_Q[i] = d_Qh[i] / tau;
//        _Q_bck[i] = Qh0[i] / tau;
//        _Q_bck[i] = Q_bck[i] / tau;
    }
//    transformPVQbck(eos, _Qbck, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck );
    transformPV(eos, _d_Q, _d_e, _d_p, _d_nb, _d_nq, _d_ns, _d_vx, _d_vy, _d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck  );
}

void Cell::getPrimVarHCenterQbck(EoS *eos, double tau, double &_e, double &_p,
                             double &_nb, double &_nq, double &_ns, double &_vx,
                             double &_vy, double &_vz) {
 double _Q_bck[7];
    for (int i = 0; i < 7; i++){
//        _Q_bck[i] = Qh0[i] / tau; -toto !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        _Q_bck[i] = Q_bck[i] / tau;
    }
    transformPVQbck(eos, _Q_bck, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz );
}

void Cell::getPrimVarPrev(EoS *eos, double tau, double &_d_e, double &_d_p,
                          double &_d_nb, double &_d_nq, double &_d_ns, double &_d_vx,
                          double &_d_vy, double &_d_vz, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck) {
    double _d_Q[7];
//    double _Q_bck[7];
    for (int i = 0; i < 7; i++){
        _d_Q[i] = d_Qprev[i] / tau;
//        _Q_bck[i] = Qbckprev[i] / tau;
//        _Q_bck[i] = Q_bck[i] / tau;
    }
//    transformPVQbck(eos, _Qbck, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck );
    transformPV(eos, _d_Q, _d_e, _d_p, _d_nb, _d_nq, _d_ns, _d_vx, _d_vy, _d_vz, e_bck, p_bck, nb_bck, nq_bck, ns_bck, vx_bck, vy_bck, vz_bck  );
}

void Cell::getPrimVarPrevQbck(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz) {
 double _Q_bck[7];
    for (int i = 0; i < 7; i++){
//        _Q_bck[i] = Qbckprev[i] / tau; -toto !!!!!!!!!!!!!!!!!!!!!!
        _Q_bck[i] = Q_bck[i] / tau;
    }
    transformPVQbck(eos, _Q_bck, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz );
}
// setting of fluctuating variables
void Cell::setPrimVar(EoS *eos, double tau, double _d_e, double _nb, double _nq,
                      double _ns, double _d_vx, double _d_vy, double _d_vz, double e_bck, double vx_bck, double vy_bck, double vz_bck) {
 const double gamma_bck = 1. / sqrt(1 - vx_bck * vx_bck - vy_bck * vy_bck - vz_bck * vz_bck);
 double p_bck = eos->p(e_bck, _nb, _nq, _ns); // background pressure
 const double d_p = eos->p(_d_e+e_bck, _nb, _nq, _ns) - eos->p(e_bck, _nb, _nq, _ns); // fluctuation of pressure - will this work even for non-linear EOS?
 d_Q[T_] = tau * ( (_d_e + d_p) * gamma_bck * gamma_bck - d_p );
 d_Q[X_] = tau * ( (e_bck + p_bck) * gamma_bck * _d_vx + (_d_e + d_p) * gamma_bck * gamma_bck * vx_bck );
 d_Q[Y_] = tau * ( (e_bck + p_bck) * gamma_bck * _d_vy + (_d_e + d_p) * gamma_bck * gamma_bck * vy_bck );
 d_Q[Z_] = tau * ( (e_bck + p_bck) * gamma_bck * _d_vz + (_d_e + d_p) * gamma_bck * gamma_bck * vz_bck );
 d_Q[NB_] = tau * _nb * gamma_bck; // not linearized
 d_Q[NQ_] = tau * _nq * gamma_bck; // not linearized
 d_Q[NS_] = tau * _ns * gamma_bck; // not linearized
 if (std::isinf(d_Q[NB_]) or std::isnan(d_Q[NB_])) {
  cout << "init error!\n";
  eos->p(_d_e, _nb, _nq, _ns);
  cout << "e = " << _d_e << " p = " << d_p << " vx = " << _d_vx << " vy = " << _d_vy
       << " vz = " << _d_vz << endl;
  //        exit(1) ;
  return;
 }
}

// function for setting of background variables
void Cell::setPrimVarQbck(EoS *eos, double tau, double _e, double _nb, double _nq,
                      double _ns, double _vx, double _vy, double _vz) {
 const double gamma2 = 1. / (1 - _vx * _vx - _vy * _vy - _vz * _vz);
 const double p = eos->p(_e, _nb, _nq, _ns);
 Q_bck[T_] = tau * (_e + p * (_vx * _vx + _vy * _vy + _vz * _vz)) * gamma2;
 Q_bck[X_] = tau * (_e + p) * _vx * gamma2;
 Q_bck[Y_] = tau * (_e + p) * _vy * gamma2;
 Q_bck[Z_] = tau * (_e + p) * _vz * gamma2;
 Q_bck[NB_] = tau * _nb * sqrt(gamma2);
 Q_bck[NQ_] = tau * _nq * sqrt(gamma2);
 Q_bck[NS_] = tau * _ns * sqrt(gamma2);
 if (std::isinf(Q_bck[NB_]) or std::isnan(Q_bck[NB_])) {
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
 cout << setw(14) << d_Q[0] / tau << setw(14) << d_Q[1] / tau << setw(14)
      << d_Q[2] / tau << setw(14) << d_Q[3] / tau << endl;
 cout << setw(14) << d_Q[4] / tau << setw(14) << d_Q[5] / tau << setw(14)
      << d_Q[6] / tau << endl;
 cout << setw(14) << d_Qh[0] / tau << setw(14) << d_Qh[1] / tau << setw(14)
      << d_Qh[2] / tau << setw(14) << d_Qh[3] / tau << endl;
 cout << setw(14) << d_Qh[4] / tau << setw(14) << d_Qh[5] / tau << setw(14)
      << d_Qh[6] / tau << endl;

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
////    if (e1+e_bck1 <= 0. || e0+e_bck <= 0.) {  // matter-vacuum
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
