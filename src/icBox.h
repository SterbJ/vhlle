#include <unistd.h>

class Fluid;
class EoS;

class IcBox {
private:
 double eta0; // midrapidity plateau
 double sigEta; // diffuseness of rapidity profile
 double ybeam; // beam rapidity
 double alphaMix; // WN/binary mixing
 double Rg; // Gaussian smearing in transverse dir
 double sNorm; // normalization of initial entropy profile
 double nNorm; // normalization of baryon density
 double sNN; // energy
 double nsigma; // width of gaussian for baryon density
 double neta0; // mean of gaussian for baryon density
 double etaM;
 double A; // initial shear flow constant
 int nx, ny, nz, nevents;
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***rho;
 double ***nrho;
 static const int NP = 10000;  // dimension for particle arrays
 // auxiliary particle arrays
 double X[NP], Y[NP], W[NP];
 int C[NP];
 double e0 = 1.;
 double p0 = 1./3;
    

 double tau0;
 int nsmoothx;  // smoothly distribute to +- this many cells
 int nsmoothy;
 int nsmoothz;
 void makeSmoothTable(int npart);
    
    double calculateSin(double x) {
     // This function calculates the sine of a given angle using the provided expansion formula.

    // Initialize variables
    double sinValue = 0.0;
    double term = x;
    double precision = 0.000001;
    int n = 1;
    
    // Calculate the sine value using the expansion formula
    while (std::abs(term) >= precision) {
        sinValue += term;
        term = -term * x * x / ((2 * n) * (2 * n + 1));
        n++;
    }
    
    return 1. + 0.01*sinValue;
}

public:
 IcBox(Fluid *f, const char *filename, double tau0, const char* setup);
 ~IcBox();
 void setIC(Fluid *f, EoS *eos);
 double setNormalization(int npart);
 double setBaryonNorm(int npart);
};
