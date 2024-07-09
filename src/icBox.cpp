#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstring>

#include "eos.h"
#include "eoChiral.h"
#include "fld.h"
#include "icBox.h"
#include "rmn.h"
#include "s95p.h"
#include "cll.h"

using namespace std;

IcBox::IcBox(Fluid* f, const char* filename, double _tau0, const char* setup) {
// cout << "loading BOX IC\n";
 nx = f->getNX();
 ny = f->getNY();
 nz = f->getNZ();
 dx = f->getDx();
 dy = f->getDy();
 dz = f->getDz();
 xmin = f->getX(0);
 xmax = f->getX(nx - 1);
 ymin = f->getY(0);
 ymax = f->getY(ny - 1);
 zmin = f->getZ(0);
 zmax = f->getZ(nz - 1);

 tau0 = _tau0;
 
 if(strcmp(setup,"LHC276")==0) {
  sNN = 2760;
  eta0 = 2.3; // midrapidity plateau
  sigEta = 1.4; // diffuseness of rapidity profile
  etaM = 4;
  ybeam = 7.98; // beam rapidity
  alphaMix = 0.15; // WN/binary mixing
  Rg = 0.4; // Gaussian smearing in transverse dir
  //sNorm = 0.96; // normalization of initial entropy profile
  A = 0.0 ; /// 3.6e-5; // initial shear flow
  cout << "IcBox: setup for 2.76 TeV LHC\n";
 } else if(strcmp(setup,"RHIC200")==0) {
  sNN = 200;
  eta0 = 1.5; // midrapidity plateau
  sigEta = 1.4; // diffuseness of rapidity profile
  etaM = 3.36;
  ybeam = 5.36; // beam rapidity
  alphaMix = 0.145; // WN/binary mixing
  Rg = 0.4; // Gaussian smearing in transverse dir
  //sNorm = 0.56; // normalization of initial entropy profile
  A = 0.0 ; // 5e-4; // initial shear flow
  nsigma = 0.6;
  neta0 = 1.4;
//  cout << "IcBox: setup for 200 GeV RHIC\n";
 } else if(strcmp(setup,"LHC5020")==0) {
  sNN = 5020;
  eta0 = 3.7; //2.3 // midrapidity plateau
  sigEta = 1.4; // diffuseness of rapidity profile
  etaM = 4.5;
  ybeam = 8.585; // beam rapidity
  alphaMix = 0.15; // WN/binary mixing
  Rg = 0.4; // Gaussian smearing in transverse dir
  A = 0.0 ; // 5e-4; // initial shear flow
  cout << "IcBox: setup for 5.02 TeV LHC\n";
 } else if(strcmp(setup,"RHIC62")==0) {
  sNN = 62.4;
  etaM = 1.8;
  ybeam = 4.2;
  alphaMix = 0.132;
  Rg = 0.4;
  A = 0.0;
  eta0 = 1.8;
  sigEta = 0.7;
  nsigma = 1.0;
  neta0 = 2.2;
  cout << "IcBox: setup for 62.4 GeV RHIC\n";
 } else if(strcmp(setup,"AFTER72")==0) {
  sNN = 72;
  etaM = 1.8;
  ybeam = 4.2;
  alphaMix = 0.132;
  Rg = 0.4;
  A = 0.0;
  eta0 = 1.8;
  sigEta = 0.7;
  nsigma = 1.0;
  neta0 = 2.2;
  cout << "IcBox: setup for 72 GeV AFTER\n";
 } else if(strcmp(setup,"RHIC27")==0) {
  sNN = 27;
  etaM = 1.0;
  ybeam = 3.36; // beam rapidity
  alphaMix = 0.123; // 0.125 WN/binary mixing
  Rg = 0.4; // Gaussian smearing in transverse dir
  A = 0.0 ; // 5e-4; // initial shear flow
  cout << "IcBox: setup for 27 GeV RHIC\n";
 } else {
  cout << "IcBox: optional parameter LHC276 or RHIC200 is expected\n";
  exit(0);
 }

 nsmoothx = (int)(3.0 * Rg / dx);  // smoothly distribute to +- this many cells
 nsmoothy = nsmoothx;

 rho = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  rho[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   rho[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    rho[ix][iy][iz] = 0.0;
   }
  }
 }
 nrho = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  nrho[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   nrho[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    nrho[ix][iy][iz] = 0.0;
   }
  }
 }
}

IcBox::~IcBox() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] rho[ix][iy];
   delete[] nrho[ix][iy];
  }
  delete[] rho[ix];
  delete[] nrho[ix];
 }
 delete[] rho;
 delete[] nrho;
}


void IcBox::setIC(Fluid* f, EoS* eos) {

 double e, p, nb, vx, vy, vz;
    
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {

    Cell* c = f->getCell(ix, iy, iz);

       double x = f->getX(ix);
       double y = f->getY(iy);
       double z = f->getZ(iz);
       double x_maximal = f->getX(nx-1);
       double x_minimal = f->getX(0);
//       double iniTime = 0.;
       double cell_size_x = f->getDx();
       double L = x_maximal - x_minimal + cell_size_x;
       
       // x-direction
       vx = 0.; // initial velocity
       vy = vz = 0.;
       double deltaE = 0.;
//       double deltaE = 0.;
       
//        diagonal xy
//       vx = 0.;
//       vy = 0.;
//       double deltaE = 0.01*sin(2 * M_PI / L * (x+y) );
       
       // 2D wave
//       vx = 0.;
//       vy = vz = 0.;
//       double deltaE = 0.01*sin(2 * M_PI / L * x ) * sin(2 * M_PI / L * y );
       
//       e = e0 + deltaE; // initial energy density
       
//       if(iy == 0 && iz == 0){
//           cout.precision(8);
//           cout << e << endl;
//       }
//       // diagonal xy
//       vx = 0.;
//       e = e0 + 0.01*sin(2 * M_PI / L * (x+y) );
       

    nb = 0.;
    c->setPrimVar(eos, tau0, deltaE, nb, 0.4*nb, 0., vx, 0., 0., e_bck, vx, vy, vz);
    c->setPrimVarQbck(eos, tau0, e_bck, nb, 0.4*nb, 0., vx, vy, vz);
       c->setAllM(1.0);

   }

}


