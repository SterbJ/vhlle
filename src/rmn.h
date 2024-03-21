#include "eos.h"
#include "inc.h"

// transformation from conserved -> primitive quantities
// input: EoS instance and conserved quantities Q
// output: energy density e, pressure p, charge densities nb, nq, ns
// 3-velocity components vx, vy, vz
void transformPVQ0(EoS *eos, double Q[7], double &e, double &p, double &nb,
                 double &nq, double &ns, double &vx, double &vy, double &vz);

void transformPV(EoS *eos, double Q[7], double &e, double &p, double &nb,
                 double &nq, double &ns, double &vx, double &vy, double &vz, double e_0, double p_0, double nb_0, double nq_0, double ns_0, double vx_0, double vy_0, double vz_0);

// the same, except that known bulk pressure Pi is taken into account
//void transformPVBulk(EoS *eos, double Pi, double Q[7], double &e, double &p,
//                     double &nb, double &nq, double &ns, double &vx, double &vy,
//                     double &vz);

// backward transformation from primitive -> conserved quantities
// input: (e,p,nb,nq,ns,vx,vy,vz)
// output: Q[7]
void transformCV(double e, double p, double nb, double nq, double ns, double vx,
                 double vy, double vz, double Q[]);
