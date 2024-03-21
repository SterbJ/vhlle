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
*                                                                             *
*******************************************************************************/

#pragma once
#include <iosfwd>
#include <algorithm>
#include <cmath>
#include "inc.h"
class EoS;

//#define NAN_DEBUG

// returns an index of pi^{mu nu} mu,nu component in a plain 1D array
int index44(const int &i, const int &j);

// this class stores the information about an individual hydro cell
class Cell {
private:
 // Q usually denotes the conserved quantities, T^{0i}
 // here, Q, Qh, Qprev etc. ~tau*T^{0i}, like Hirano'01
 double Q[7];      // final values at a given timestep
 double Q0[7];     // background values at a given timestep
 double Qh[7];     // half-step updated values
 double Qh0[7];     // half-step updated values
 double Qprev[7];  // values at the end of previous timestep
 double Q0prev[7];  // values at the end of previous timestep
 double pi[10], piH[10];  // pi^{mu nu}, WITHOUT tau factor, final (pi) and
                          // half-step updated (piH)
 double Pi, PiH;  // Pi, WITHOUT tau factor, final (Pi) and half-step updated (PiH)
 double pi0[10], piH0[10];  // // pi^{mu nu}, WITHOUT tau factor, auxiliary
 double pi_bck[10]={0,0,0,0,0,0,0,0,0,0}, piH_bck[10]={0,0,0,0,0,0,0,0,0,0};  // background
 double pi0_bck[10]={0,0,0,0,0,0,0,0,0,0}, piH0_bck[10]={0,0,0,0,0,0,0,0,0,0};  // background prev
 double pi_bck_prev[10]={0,0,0,0,0,0,0,0,0,0}, piH_bck_prev[10]={0,0,0,0,0,0,0,0,0,0};  // background prev
 double pi0_bck_prev[10]={0,0,0,0,0,0,0,0,0,0}, piH0_bck_prev[10]={0,0,0,0,0,0,0,0,0,0};  // background prev

 double Pi0, PiH0;          // viscous, WITHOUT tau factor, auxiliary
 double Pi_bck=0, PiH_bck=0;          // background
 double Pi0_bck=0, PiH0_bck=0;          // background
 double Pi_bck_prev=0, PiH_bck_prev=0;          // background prev
 double Pi0_bck_prev=0, PiH0_bck_prev=0;          // background prev
 double flux[7];            // cumulative fluxes
 Cell *next[3];             // pointer to the next cell in a given direction
 Cell *prev[3];             // pointer to the previous cell in a given direction
 double m[3];               // extend of matter propagation inside cell [0...1]
 double dm[3];              // auxiliary
 int ix, iy, iz;            // cell coordinate on the grid
 // viscCorrCut: flag if the viscous corrections are cut for this cell:
 // 1.0 = uncut, < 1 :  cut by this factor
 double viscCorrCut;

    
public:
 Cell();
 ~Cell(){};
 inline void setPos(int iix, int iiy, int iiz) {
  ix = iix;
  iy = iiy;
  iz = iiz;
 }
 inline int getX(void) { return ix; }
 inline int getY(void) { return iy; }
 inline int getZ(void) { return iz; }

 inline void setQ(double *_Q) {
  for (int i = 0; i < 7; i++) Q[i] = _Q[i];
//  if (Q[T_] < 0.) {
//   for (int i = 0; i < 7; i++) Q[i] = 0.;
//  }
 }
 inline void setQ0(double *_Q0) {
  for (int i = 0; i < 7; i++) Q0[i] = _Q0[i];
  if (Q0[T_] < 0.) {
  for (int i = 0; i < 7; i++) Q0[i] = 0.;
  }
 }
 inline void setQh(double *_Qh) {
  for (int i = 0; i < 7; i++) Qh[i] = _Qh[i];
//  if (Qh[T_] < 0.) {
//   for (int i = 0; i < 7; i++) Qh[i] = 0.;
//  }
 }
    inline void setQh0(double *_Qh0) {
     for (int i = 0; i < 7; i++) Qh0[i] = _Qh0[i];
     if (Qh0[T_] < 0.) {
      for (int i = 0; i < 7; i++) Qh0[i] = 0.;
     }
    }

 // getter and setter methods for the class members
    
 inline double getpi(const int &i, const int &j) { return pi[index44(i, j)]; }
 inline double getpiH(const int &i, const int &j) { return piH[index44(i, j)]; }
 inline double getpi0(const int &i, const int &j) { return pi0[index44(i, j)]; }
 inline double getpiH0(const int &i, const int &j) {
  return piH0[index44(i, j)];
 }
    inline double getpi_bck(const int &i, const int &j) { return pi_bck[index44(i, j)]; }
    inline double getpiH_bck(const int &i, const int &j) { return piH_bck[index44(i, j)]; }
    inline double getpi0_bck(const int &i, const int &j) { return pi0_bck[index44(i, j)]; }
    inline double getpiH0_bck(const int &i, const int &j) {
     return piH0_bck[index44(i, j)];
    }
    
    inline double getpi_bck_prev(const int &i, const int &j) { return pi_bck_prev[index44(i, j)]; }
    inline double getpiH_bck_prev(const int &i, const int &j) { return piH_bck_prev[index44(i, j)]; }
    inline double getpi0_bck_prev(const int &i, const int &j) { return pi0_bck_prev[index44(i, j)]; }
    inline double getpiH0_bck_prev(const int &i, const int &j) {
     return piH0_bck_prev[index44(i, j)];
    }
    
 inline double getPi(void) { return Pi; }
 inline double getPiH(void) { return PiH; }
 inline double getPi0(void) { return Pi0; }
 inline double getPiH0(void) { return PiH0; }
    inline double getPi_bck(void) { return Pi_bck; }
    inline double getPiH_bck(void) { return PiH_bck; }
    inline double getPi0_bck(void) { return Pi0_bck; }
    inline double getPiH0_bck(void) { return PiH0_bck; }
    
    inline double getPi_bck_prev(void) { return Pi_bck_prev; }
    inline double getPiH_bck_prev(void) { return PiH_bck_prev; }
    inline double getPi0_bck_prev(void) { return Pi0_bck_prev; }
    inline double getPiH0_bck_prev(void) { return PiH0_bck_prev; }

 inline void setpi(const int &i, const int &j, const double &val) {
  pi[index44(i, j)] = val;
 }
 inline void setpiH(const int &i, const int &j, const double &val) {
  piH[index44(i, j)] = val;
 }
 inline void setpi0(const int &i, const int &j, const double &val) {
  pi0[index44(i, j)] = val;
 }
 inline void setpiH0(const int &i, const int &j, const double &val) {
  piH0[index44(i, j)] = val;
 }
 inline void addpi0(const int &i, const int &j, const double &val) {
  pi0[index44(i, j)] += val;
 }
 inline void addpiH0(const int &i, const int &j, const double &val) {
  piH0[index44(i, j)] += val;
 }
 inline void setPi(const double &val) { Pi = val; }
 inline void setPiH(const double &val) { PiH = val; }
 inline void setPi0(const double &val) { Pi0 = val; }
 inline void setPiH0(const double &val) { PiH0 = val; }
 inline void addPi0(const double &val) { Pi0 += val; }
 inline void addPiH0(const double &val) { PiH0 += val; }

 inline void getQ(double *_Q) {
  for (int i = 0; i < 7; i++) _Q[i] = Q[i];
 }
 inline void getQ0(double *_Q0) {
  for (int i = 0; i < 7; i++) _Q0[i] = Q0[i];
 }
 inline void getQh(double *_Qh) {
  for (int i = 0; i < 7; i++) _Qh[i] = Qh[i];
 }
    inline void getQh0(double *_Qh0) {
     for (int i = 0; i < 7; i++) _Qh0[i] = Qh0[i];
    }
 inline void getQprev(double *_Qp) {
  for (int i = 0; i < 7; i++) _Qp[i] = Qprev[i];
 }
 inline void getQ0prev(double *_Q0p) {
     for (int i = 0; i < 7; i++) _Q0p[i] = Q0prev[i];
 }
 inline void saveQprev(void) {
  for (int i = 0; i < 7; i++) Qprev[i] = Q[i];
 }
 inline void saveQ0prev(void) {
     for (int i = 0; i < 7; i++) Q0prev[i] = Q0[i];
 }

 // imports Q, Qh, Qprev, Pi, pi, m, viscCorrCut from cell c
 void importVars(Cell* c);

 inline void setNext(int i, Cell *c) { next[i - 1] = c; }
 inline void setPrev(int i, Cell *c) { prev[i - 1] = c; }
 inline Cell *getNext(int i) { return next[i - 1]; }
 inline Cell *getPrev(int i) { return prev[i - 1]; }

 inline void setAllM(double value) { m[0] = m[1] = m[2] = value; }
 inline void addM(int dir, double inc) {
  m[dir - 1] += inc;
  if (m[dir - 1] > 0.9)
   for (int i = 0; i < 3; i++) m[i] = 1.;
 }
 inline double getM(int dir) { return m[dir - 1]; }
 inline double getMaxM(void) { return std::max(m[0], std::max(m[1], m[2])); }
 inline void setDM(int dir, double value) { dm[dir - 1] = value; }
 inline double getDM(int dir) { return dm[dir - 1]; }

 inline void setpi0(double values[4][4]) {
  for (int i = 0; i < 4; i++)
   for (int j = 0; j < 4; j++) pi0[index44(i, j)] = values[i][j];
 }
 inline void setpiH0(double values[4][4]) {
  for (int i = 0; i < 4; i++)
   for (int j = 0; j < 4; j++) piH0[index44(i, j)] = values[i][j];
 }

 // get the energy density, pressure, charge densities and flow velocity
 // components (e,p,n,v) from conserved quantities Q in the centre of the cell
 void getPrimVar(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                 double &_nq, double &_ns, double &_vx, double &_vy,
                 double &_vz, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0);
 void getPrimVarQ0(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                 double &_nq, double &_ns, double &_vx, double &_vy,
                 double &_vz);
 // (e,p,n,v) at cell's left boundary in a given direction dir
 void getPrimVarLeft(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                     double &_nq, double &_ns, double &_vx, double &_vy,
                     double &_vz, int dir, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0);
    void getPrimVarLeftQ0(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                        double &_nq, double &_ns, double &_vx, double &_vy,
                        double &_vz, int dir);
 // (e,p,n,v) at cell's right boundary in a given direction dir
 void getPrimVarRight(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                      double &_nq, double &_ns, double &_vx, double &_vy,
                      double &_vz, int dir, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0);
    void getPrimVarRightQ0(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                         double &_nq, double &_ns, double &_vx, double &_vy,
                         double &_vz, int dir);

 // (e,p,n,v) from half-step updated Qh at cell's left boundary in a given
 // direction
 void getPrimVarHLeft(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                      double &_nq, double &_ns, double &_vx, double &_vy,
                      double &_vz, int dir, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0);
 void getPrimVarHLeftQ0(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                         double &_nq, double &_ns, double &_vx, double &_vy,
                         double &_vz, int dir);
 // (e,p,n,v) from half-step updated Qh at cell's right boundary in a given
 // direction
 void getPrimVarHRight(EoS *eos, double tau, double &_e, double &_p,
                       double &_nb, double &_nq, double &_ns, double &_vx,
                       double &_vy, double &_vz, int dir, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0);
 void getPrimVarHRightQ0(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz, int dir);
 // (e,p,n,v) from half-step updated Qh at cell's centre
 void getPrimVarHCenter(EoS *eos, double tau, double &_e, double &_p,
                        double &_nb, double &_nq, double &_ns, double &_vx,
                        double &_vy, double &_vz, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0);
 void getPrimVarHCenterQ0(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz);
 // (e,p,n,v) at the previous timestep and cell's centre
 void getPrimVarPrev(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                     double &_nq, double &_ns, double &_vx, double &_vy,
                     double &_vz, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0);
 void getPrimVarPrevQ0(EoS *eos, double tau, double &_e, double &_p,
                              double &_nb, double &_nq, double &_ns, double &_vx,
                                double &_vy, double &_vz);
 // calculate and set Q from (e,n,v)
 void setPrimVar(EoS *eos, double tau, double _e, double _nb, double _nq,
                 double _ns, double _vx, double _vy, double _vz, double e0, double vx0, double vy0, double vz0);
 void setPrimVarQ0(EoS *eos, double tau, double _e, double _nb, double _nq,
                    double _ns, double _vx, double _vy, double _vz);

 // update the cumulative fluxes through the cell
 inline void addFlux(double Ft, double Fx, double Fy, double Fz, double Fnb,
                     double Fnq, double Fns) {
  #ifdef NAN_DEBUG
  if(std::isinf(Ft) or std::isnan(Ft)) {
   std::cout << "Cell::addFlux inf/nan\n";
  }
  #endif
  flux[T_] += Ft;
  flux[X_] += Fx;
  flux[Y_] += Fy;
  flux[Z_] += Fz;
  flux[NB_] += Fnb;
  flux[NQ_] += Fnq;
  flux[NS_] += Fns;
 }
 inline void clearFlux(void) {
  for (int i = 0; i < 7; i++) flux[i] = 0.;
 }
 void updateByFlux();       // Q = Q + flux
 void updateByViscFlux();   // this limits the update based on flux[0]/Q[0] ratio
 void updateQtoQhByFlux();  // Qh = Q + flux
 inline void setViscCorrCutFlag(double value) { viscCorrCut = value; }
 inline double getViscCorrCutFlag(void) { return viscCorrCut; }
 void Dump(double tau);  // dump the contents of the cell into dump.dat
// void d_pi(int ix, int iy, int iz, double dpi0[4][4][4]);
};
