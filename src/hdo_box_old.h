
class Cell;
class Fluid;
class EoS;
class TransportCoeff;

// this class implements the hydrodynamic evolution and
// contains the hydrodynamic algorithm
class Hydro {
private:
 Fluid *f;
 EoS *eos;
 TransportCoeff *trcoeff;
 double dt, tau;  // dt: timestep, tau: current value of the proper time
 #ifdef CARTESIAN
 double t;  // time in Cartesian frame
 #endif
 double tau_z;    // effective value of the proper time used in 1/tau factors in
                  // the fluxes. Used to increase the accuracy
public:
 Hydro(Fluid *_f, EoS *_eos, TransportCoeff *_trcoeff, double _t0, double _dt);
 ~Hydro();
 void setDtau(double deltaTau);  // change the timestep
 double getDtau() { return dt; }  // return current value of timestep
 void setFluid(Fluid *_f) { f = _f; }
 Fluid* getFluid() { return f; }

 // HLLE (ideal)flux between two neighbouring cells in a given direction
 // mode: PREDICT = used in predictor step; calculates fluxes for dt/2
 // CORRECT = used in corrector step, calculates fluxes based on predicted
 // half-step quantities
 void hlle_flux(Cell *left, Cell *right, int direction, int mode, int ix);
 // viscous flux \delta F
 void visc_flux(Cell *left, Cell *right, int direction, double ix, double iy, double iz);
 // viscous source step for a given cell (ix,iy,iz)
 void visc_source_step(int ix, int iy, int iz);
 void source(double tau, double x, double y, double z, double Q[7],
             double S[7]);
 // ideal source step for a given cell (ix,iy,iz)
 void source_step(int ix, int iy, int iz, int mode);
 // shear stress tensor and bulk pressure in Navier-Stokes (NS) limit
 // plus \partial_\mu u^\nu matrix (dmu) and
 // expansion scalar \partial_mu u^\mu (du)
 // for a given cell (ix,iy,iz)
 void NSquant(int ix, int iy, int iz, double pi[][4], double &Pi,
              double dmu[4][4], double dmu0[4][4], double &du, double &du0);
 // sets the values of shear stress/bulk pressure in NS limit in all hydro grid
// void setNSvalues();
 // advances numerical solution for shear/bulk in a whole grid over one
 // timestep
 void ISformal();
 // advances numerical solution for Q (including ideal and viscous fluxes and
 // source terms) over one timestep
 void performStep(double ctime, double Struc[1000][31][31][31], int event);
 // gets the current proper time
 inline double getTau() { return tau; }
 #ifdef CARTESIAN
 double time() { return t; }
 #endif
    
 int sgn(double v) {
  if (v < 0) return -1;
  if (v == 0) return 0;
  if (v > 0) return 1;
 }
    
};