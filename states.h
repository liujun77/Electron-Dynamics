#ifndef STATE_H
#define STATE_H

class State
{
    public:
        double energy;
        int l;
        int n;
        int m;
        double spin;
        int Eind;
};

class Selection
{
    public:
        int Ei;
        int Ej;
        int mij;
        int k;
        int l;
};

double fermi_dirac(double energy, double temp);

int get_free_energies(std::string filename);

void ini_n_l_E(std::string filename);

double psiR(int i, double r);

void normalize_coeff(void);

void ini_gauss(int n_inte_gauss);

double scatter_ee_INT2(int i, int j, int k, int l);

double scatter_cos_INT2(int i, int j);

double scatter_rcos_INT2(int i, int j);

void generator_electron_photon_rate(double wp);

double exite_e(int i, int k, double temp, double wp);

double scatter_ep(int i, int k, double temp);

double total_energy(double p[], int n_p);

void ini_epsilon(void);

void ini_selection(void);

double evolution(double p[],
               int n_p,
               double dt,
               double t,
               double temp,
               double wp);

void debug_n(double p[]);

void ini_p(double p[], double temp);

void state_2_energy(double p[],
                    double x[],
                    double y[]);

void initial_excitation(double p[], double temp, double wp, double pht, double dt);

void hot_carrier(double x[],
                 double p[],
                 double temp,
                 double y[],
                 double yh[],
                 double ye[]);

double fermi_energy(double temp);

void delete_cache(void);

void save_electron_photon_matrix(void);

void ini_IJ_Photon(void);

void evolution_photon(double p[], double temp, double dt);

void save_w_photon(std::string filename);

void save_convolution(std::string infilename,
                      std::string outfilename,
                      double sigma,
                      int ntimes,
                      double xmin,
                      double xmax);

void save_convolution_photon(std::string infilename, double sigma, int ntimes);

void save_convolution_hotcarrier(std::string infilename, double sigma, int ntimes);

#endif
