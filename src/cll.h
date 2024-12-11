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
 double d_Q[7];      // final values at a given timestep
 double Q_bck[7];     // background values at a given timestep
 double d_Qh[7];     // half-step updated values
 double Qh_bck[7];     // half-step updated values
 double d_Qprev[7];  // values at the end of previous timestep
 double Qprev_bck[7];  // values at the end of previous timestep
 double d_pi[10], d_piH[10];  // pi^{mu nu}, WITHOUT tau factor, final (pi) and
                          // half-step updated (piH)
 double d_Pi, d_PiH;  // Pi, WITHOUT tau factor, final (Pi) and half-step updated (PiH)
 double d_pi0[10], d_piH0[10];  // // pi^{mu nu}, WITHOUT tau factor, auxiliary
 double pi_bck[10], piH_bck[10];  // background pi^munu
 double pi0_bck[10], piH0_bck[10];  // background pi^munu without tau factor
 double pi_bck_prev[10], piH_bck_prev[10];  // background pi^munu at (half)previous step
 double pi0_bck_prev[10], piH0_bck_prev[10];  // background pi^munu at previous step withou tau factor
 double d_Pi0, d_PiH0;          // viscous, WITHOUT tau factor, auxiliary
 double Pi_bck, PiH_bck;          // background
 double Pi0_bck, PiH0_bck;          // background
 double Pi_bck_prev, PiH_bck_prev;          // background prev
 double Pi0_bck_prev, PiH0_bck_prev;          // background prev
 double d_flux[7];            // cumulative fluxes
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

 inline void setQ(double *_Q) {//sets fluctuationg Q
  for (int i = 0; i < 7; i++) d_Q[i] = _Q[i];
//  if (Q[T_] < 0.) {
//   for (int i = 0; i < 7; i++) Q[i] = 0.;
//  }
 }
 inline void setQbck(double *_Q_bck) {//sets background Q
  for (int i = 0; i < 7; i++) Q_bck[i] = _Q_bck[i];
  if (Q_bck[T_] < 0.) {
  for (int i = 0; i < 7; i++) Q_bck[i] = 0.;
  }
 }
 inline void setQh(double *_Qh) {
  for (int i = 0; i < 7; i++) d_Qh[i] = _Qh[i];
//  if (Qh[T_] < 0.) {
//   for (int i = 0; i < 7; i++) Qh[i] = 0.;
//  }
 }
    inline void setQhbck(double *_Qh_bck) {
     for (int i = 0; i < 7; i++) Qh_bck[i] = _Qh_bck[i];
     if (Qh_bck[T_] < 0.) {
      for (int i = 0; i < 7; i++) Qh_bck[i] = 0.;
     }
    }

 // getter and setter methods for the class members
    
 inline double getpi(const int &i, const int &j) { return d_pi[index44(i, j)]; }
 inline double getpiH(const int &i, const int &j) { return d_piH[index44(i, j)]; }
 inline double getpi0(const int &i, const int &j) { return d_pi0[index44(i, j)]; }
 inline double getpiH0(const int &i, const int &j) {
  return d_piH0[index44(i, j)];
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
    
 inline double getPi(void) { return d_Pi; }
 inline double getPiH(void) { return d_PiH; }
 inline double getPi0(void) { return d_Pi0; }
 inline double getPiH0(void) { return d_PiH0; }
    inline double getPi_bck(void) { return Pi_bck; }
    inline double getPiH_bck(void) { return PiH_bck; }
    inline double getPi0_bck(void) { return Pi0_bck; }
    inline double getPiH0_bck(void) { return PiH0_bck; }
    
    inline double getPi_bck_prev(void) { return Pi_bck_prev; }
    inline double getPiH_bck_prev(void) { return PiH_bck_prev; }
    inline double getPi0_bck_prev(void) { return Pi0_bck_prev; }
    inline double getPiH0_bck_prev(void) { return PiH0_bck_prev; }

 inline void setpi(const int &i, const int &j, const double &val) {
  d_pi[index44(i, j)] = val;
 }
 inline void setpiH(const int &i, const int &j, const double &val) {
  d_piH[index44(i, j)] = val;
 }
 inline void setpi0(const int &i, const int &j, const double &val) {
  d_pi0[index44(i, j)] = val;
 }
 inline void setpiH0(const int &i, const int &j, const double &val) {
  d_piH0[index44(i, j)] = val;
 }
 inline void addpi0(const int &i, const int &j, const double &val) {
  d_pi0[index44(i, j)] += val;
 }
 inline void addpiH0(const int &i, const int &j, const double &val) {
  d_piH0[index44(i, j)] += val;
 }
 inline void setPi(const double &val) { d_Pi = val; }
 inline void setPiH(const double &val) { d_PiH = val; }
 inline void setPi0(const double &val) { d_Pi0 = val; }
 inline void setPiH0(const double &val) { d_PiH0 = val; }
 inline void addPi0(const double &val) { d_Pi0 += val; }
 inline void addPiH0(const double &val) { d_PiH0 += val; }

 inline void getQ(double *_Q) {
  for (int i = 0; i < 7; i++) _Q[i] = d_Q[i];
 }
 inline void getQbck(double *_Q_bck) {
  for (int i = 0; i < 7; i++) _Q_bck[i] = Q_bck[i];
 }
 inline void getQh(double *_Qh) {
  for (int i = 0; i < 7; i++) _Qh[i] = d_Qh[i];
 }
    inline void getQhbck(double *_Qh_bck) {
     for (int i = 0; i < 7; i++) _Qh_bck[i] = Qh_bck[i];
    }
 inline void getQprev(double *_Qp) {
  for (int i = 0; i < 7; i++) _Qp[i] = d_Qprev[i];
 }
 inline void getQprevbck(double *_Qp_bck) {
     for (int i = 0; i < 7; i++) _Qp_bck[i] = Qprev_bck[i];
 }
 inline void saveQprev(void) {
  for (int i = 0; i < 7; i++) d_Qprev[i] = d_Q[i];
 }
 inline void saveQprevbck(void) {
     for (int i = 0; i < 7; i++) Qprev_bck[i] = Q_bck[i];
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
   for (int j = 0; j < 4; j++) d_pi0[index44(i, j)] = values[i][j];
 }
 inline void setpiH0(double values[4][4]) {
  for (int i = 0; i < 4; i++)
   for (int j = 0; j < 4; j++) d_piH0[index44(i, j)] = values[i][j];
 }

 // get the energy density, pressure, charge densities and flow velocity
 // components (e,p,n,v) from conserved quantities Q in the centre of the cell
 void getPrimVar(EoS *eos, double tau, double &_d_e, double &_d_p, double &_d_nb,
                 double &_d_nq, double &_d_ns, double &_d_vx, double &_d_vy,
                 double &_d_vz, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck);
 void getPrimVarQbck(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                 double &_nq, double &_ns, double &_vx, double &_vy,
                 double &_vz);
 // (e,p,n,v) at cell's left boundary in a given direction dir
 void getPrimVarLeft(EoS *eos, double tau, double &_d_e, double &_d_p, double &_d_nb,
                     double &_d_nq, double &_d_ns, double &_d_vx, double &_d_vy,
                     double &_d_vz, int dir, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck);
 void getPrimVarLeftQbck(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                        double &_nq, double &_ns, double &_vx, double &_vy,
                        double &_vz, int dir);
 // (e,p,n,v) at cell's right boundary in a given direction dir
 void getPrimVarRight(EoS *eos, double tau, double &_d_e, double &_d_p, double &_d_nb,
                      double &_d_nq, double &_d_ns, double &_d_vx, double &_d_vy,
                      double &_d_vz, int dir, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck);
 void getPrimVarRightQbck(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                         double &_nq, double &_ns, double &_vx, double &_vy,
                         double &_vz, int dir);

 // (e,p,n,v) from half-step updated Qh at cell's left boundary in a given
 // direction
 void getPrimVarHLeft(EoS *eos, double tau, double &_d_e, double &_d_p, double &_d_nb,
                      double &_d_nq, double &_d_ns, double &_d_vx, double &_d_vy,
                      double &_d_vz, int dir, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck);
 void getPrimVarHLeftQbck(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                         double &_nq, double &_ns, double &_vx, double &_vy,
                         double &_vz, int dir);
 // (e,p,n,v) from half-step updated Qh at cell's right boundary in a given
 // direction
 void getPrimVarHRight(EoS *eos, double tau, double &_d_e, double &_d_p,
                       double &_d_nb, double &_d_nq, double &_d_ns, double &_d_vx,
                       double &_d_vy, double &_d_vz, int dir, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck);
 void getPrimVarHRightQbck(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz, int dir);
 // (e,p,n,v) from half-step updated Qh at cell's centre
 void getPrimVarHCenter(EoS *eos, double tau, double &_d_e, double &_d_p,
                        double &_d_nb, double &_d_nq, double &_d_ns, double &_d_vx,
                        double &_d_vy, double &_d_vz, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck);
 void getPrimVarHCenterQbck(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz);
 // (e,p,n,v) at the previous timestep and cell's centre
 void getPrimVarPrev(EoS *eos, double tau, double &_d_e, double &_d_p, double &_d_nb,
                     double &_d_nq, double &_d_ns, double &_d_vx, double &_d_vy,
                     double &_d_vz, double e_bck, double p_bck, double nb_bck, double nq_bck, double ns_bck, double vx_bck, double vy_bck, double vz_bck);
 void getPrimVarPrevQbck(EoS *eos, double tau, double &_e, double &_p,
                              double &_nb, double &_nq, double &_ns, double &_vx,
                                double &_vy, double &_vz);
 // calculate and set Q from (e,n,v)
 void setPrimVar(EoS *eos, double tau, double _d_e, double _d_nb, double _d_nq,
                 double _d_ns, double _d_vx, double _d_vy, double _d_vz, double e_bck, double vx_bck, double vy_bck, double vz_bck);
 void setPrimVarQbck(EoS *eos, double tau, double _e, double _nb, double _nq,
                    double _ns, double _vx, double _vy, double _vz);

 // update the cumulative fluxes through the cell
 inline void addFlux(double Ft, double Fx, double Fy, double Fz, double Fnb,
                     double Fnq, double Fns) {
  #ifdef NAN_DEBUG
  if(std::isinf(Ft) or std::isnan(Ft)) {
   std::cout << "Cell::addFlux inf/nan\n";
  }
  #endif
  d_flux[T_] += Ft;
  d_flux[X_] += Fx;
  d_flux[Y_] += Fy;
  d_flux[Z_] += Fz;
  d_flux[NB_] += Fnb;
  d_flux[NQ_] += Fnq;
  d_flux[NS_] += Fns;
 }
 inline void getFlux(double F[7]) {
     F[0] = d_flux[T_];
     F[1] = d_flux[X_];
     F[2] = d_flux[Y_];
     F[3] = d_flux[Z_];
     F[4] = d_flux[NB_];
     F[5] = d_flux[NQ_];
     F[6] = d_flux[NS_];
 }
 inline void setFlux(double Ft, double Fx, double Fy, double Fz, double Fnb,
                     double Fnq, double Fns) {
     d_flux[T_] = Ft;
     d_flux[X_] = Fx;
     d_flux[Y_] = Fy;
     d_flux[Z_] = Fz;
     d_flux[NB_] = Fnb;
     d_flux[NQ_] = Fnq;
     d_flux[NS_] = Fns;
 }
 inline void clearFlux(void) {
  for (int i = 0; i < 7; i++) d_flux[i] = 0.;
 }
 void updateByFlux();       // Q = Q + flux
 void updateByViscFlux();   // this limits the update based on flux[0]/Q[0] ratio
 void updateQtoQhByFlux();  // Qh = Q + flux
 inline void setViscCorrCutFlag(double value) { viscCorrCut = value; }
 inline double getViscCorrCutFlag(void) { return viscCorrCut; }
 void Dump(double tau);  // dump the contents of the cell into dump.dat
};
