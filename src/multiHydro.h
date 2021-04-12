
class Fluid;
class Hydro;
class EoS;
class TransportCoeff;
class CrossSections;

class MultiHydro {
 Fluid *f_p, *f_t, *f_f; // projectile and target fluids
 Hydro *h_p, *h_t, *h_f; // solutions of hydro eqs. for projectile and target
 EoS *eos;
 TransportCoeff *trcoeff;
 CrossSections *xsect;

public:
 MultiHydro(Fluid *f_p, Fluid *f_t, Fluid *f_f, Hydro *h_p, Hydro *h_t,
  Hydro* h_f, EoS *eos, TransportCoeff *trcoeff);
 void performStep();
 void frictionSubstep();
 void findFreezeout();
};