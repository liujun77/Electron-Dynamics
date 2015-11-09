#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "const.h"
#include "save.h"
#include "states.h"
#include "readin.hpp"
#include "bessel_utilities.h"
#include "integration.hpp"
#include "legendre_rule_fast.h"
#include <gsl/gsl_sf_coupling.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/bessel.hpp>
#ifdef OPENMP
#include <omp.h>
#endif
#include "convolution.h"

//using declare
using namespace std;
using jun::sph_bessel_jp_j;
using jun::sph_hankel_ip_h;
using jun::sph_hankel_i;
using boost::math::sph_bessel;
using boost::math::spherical_harmonic;
using jun::I;
using jun::PI;
using jun::C_AU;
using jun::AU_EV;
using jun::Boltzman_k;
using jun::save;

//external variables
extern double V0;
extern double ex_val;
extern int n_states;
extern int n_ele;
extern double r0;
extern double radius;
extern double EF;
extern int ss;
extern double exite;
extern double ephonon;
extern double ee_a;
extern double ee_b;
extern double ee_c;
extern double ee_de;
extern State *state;
extern char dir[256];
extern int steps;

//*************************************************************************
//global var in this file
double *X;
double *W;
int n_inte;
double *COEF_J;
complex<double> *COEF_H;
double *Scatter_Electron_Photon;
double *Scatter_Exite;
double *Epsilon_r;
double *Epsilon_i;
double *IJ_Photon;
double *E_s;
int *d_s;
int *l_s;
int lmax;
double epsilon_gamma;

double fermi_dirac(double energy, double temp)
{
    if(energy<=EF)
        return 1;
    else {
        return 0;
    }
}

double s_func_energy(int l, double energy)
{
    double alfa, k;
    complex<double> val;
    alfa=sqrt(2*(energy-V0));
    k=sqrt(-2*energy);
    val=alfa*sph_bessel_jp_j(l,alfa*radius)-I*k*sph_hankel_ip_h(l,k*radius);
    return real(val);
}

double search_root(double x1, double x2, int l)
{
    //search the root between x1 &x2 of function func_energy
    double mid,fx1,fx2,fmid;
    do{
        mid=(x1+x2)/2;
        fx1=s_func_energy(l,x1);
        fx2=s_func_energy(l,x2);
        fmid=s_func_energy(l,mid);
        if(fx1*fx2>0||(fabs((fx1-fx2)/(x1-x2))>1e6&&fmid>10))
            return 1;   //avoid the divergence case means there is no root
        else{
            if(fx1*fmid>0)
                x1=mid;
            else x2=mid;
        }
    }while(fabs(fmid)>1e-10);
    return mid;
}

int get_free_energies(string filename)
{
    int l=-1,nnl=1;
    int i=0;
    double y1,y2,energy,root;
    double de=1e-5;
    ofstream filein(filename.c_str());
    while(nnl!=0){
        nnl=0;
        ++l;
        y1=s_func_energy(l,V0+de);
        for(energy=V0+de;energy<-de;energy+=de){
            y2=s_func_energy(l,energy+de);
            if(y1*y2<0){
                root=search_root(energy,energy+de,l);
                if(root<0){
                    filein<<l<<"\t"<<nnl<<"\t"<<root<<endl;
                    i++;
                    ++nnl;
                }
            }
            y1=y2;
        }
    }
    filein.close();
    lmax=l-1;
    cout<<"lmax = "<<lmax<<endl;
    return i;
}

int compare_State(const void* a, const void* b)
{
    double da,db;
    if(((const State*)a)->energy>((const State*)b)->energy)
        return 1;
    if(((const State*)a)->energy<((const State*)b)->energy)
        return -1;
    if(((const State*)a)->energy==((const State*)b)->energy){
        da=((const State*)a)->l;
        db=((const State*)b)->l;
        if(da>db)
            return 1;
        else
            return -1;
    }
    return 0;
}

void ini_n_l_E(string filename)
{
    int i;
    string line;
    state=new State [n_states];
    ifstream data_file(filename.c_str());
    for(i=0;i<n_states;++i){
        getline(data_file,line);
        istringstream line_str(line);
        line_str>>state[i].l>>state[i].n>>state[i].energy;
    }
    data_file.close();
    qsort(state,n_states,sizeof(State),compare_State);
    E_s=new double [n_states];
    l_s=new int [n_states];
    d_s=new int [n_states];
    for(i=0;i<n_states;i++){
        E_s[i]=state[i].energy;
        l_s[i]=state[i].l;
        d_s[i]=2*l_s[i]+1;
    }
#ifdef DEBUG
    ofstream debug_test("states");
    for(i=0;i<n_states;i++)
        debug_test<<state[i].l<<"\t"<<state[i].n<<"\t"<<state[i].energy<<endl;
    debug_test.close();
#endif
}

double f1norm(double r, int i)
{
    double alfa=sqrt(2*(state[i].energy-V0));
    int l=state[i].l;
    double re;
    if(r==0)
        return 0;
    else{
        re=r*r*sph_bessel(l,alfa*r)*sph_bessel(l,alfa*r);
        return re;
    }
}

double f2norm(double r, int i)
{
    double k=sqrt(-2*state[i].energy);
    int l=state[i].l;
    double auxr=real(sph_hankel_i(l,k*r));
    double auxi=imag(sph_hankel_i(l,k*r));
#ifdef DEBUG_STATES
    if(isnan(auxr)){
        cout<<"auxr is nan i="<<i
            <<" r="<<r<<endl;
        exit(0);
    }
    if(isnan(auxi)){
        cout<<"auxi is nan i="<<i
            <<" r="<<r<<endl;
        exit(0);
    }
#endif
    return r*r*(auxr*auxr+auxi*auxi);
}

void normalize_coeff(void)
{
    int i,n=n_states;
    int l;
    double inte1, inte2, inte;
    double k, alfa;
    COEF_J=new double [n_states];
    COEF_H=new complex<double> [n_states];
    for(i=0;i<n;++i){
        l=state[i].l;
        alfa=sqrt(2*(state[i].energy-V0));
        k=sqrt(-2*state[i].energy);
        inte1=jun::integ(0.,radius,f1norm,i,200);
        inte2=jun::integ(radius,1.5*radius,f2norm,i,100);
#ifdef DEBUG_STATES
        if(isnan(inte1)){
            cout<<"inte1 is nan i = "<<i<<endl;
            exit(0);
        }
        if(isnan(inte2)){
            cout<<"inte2 is nan i = "<<i<<endl;
            exit(0);
        }
#endif
        inte=inte1
            +norm(sph_bessel(l,alfa*radius)/sph_hankel_i(l,k*radius))*inte2;
        COEF_J[i]=sqrt(1/inte);
        COEF_H[i]=sph_bessel(l,alfa*radius)/sph_hankel_i(l,k*radius)*COEF_J[i];
  }
}

double psiR(int i, double r)
{
    double alfa,k;
    int l=state[i].l;
    alfa=sqrt(2*(state[i].energy-V0));
    k=sqrt(-2*state[i].energy);
    if(r<radius)
        return COEF_J[i]*sph_bessel(l,alfa*r);
    else
        return real(COEF_H[i]*sph_hankel_i(l,k*r));
}

void ini_gauss(int n_inte_gauss)
{
    n_inte=n_inte_gauss;
    X=new double[n_inte];
    W=new double[n_inte];
    legendre_compute_glr(n_inte,X,W);
    rescale(0,1.2*radius,n_inte,X,W);
}

double scatter_ee_INT2(int i, int j, int k, int l)
{
    int max_ij=max(l_s[i],l_s[j]);
    int min_ij=min(l_s[i],l_s[j]);
    int max_kl=max(l_s[k],l_s[l]);
    int min_kl=min(l_s[k],l_s[l]);
    if(abs(max_ij-max_kl)>1)
        return 0;
    if(abs(min_ij-min_kl)>1)
        return 0;
    if(max(i,j)>=max(k,l)&&min(i,j)>=min(k,l))
        return 0;
    if(max(i,j)<=max(k,l)&&min(i,j)<=min(k,l))
        return 0;
    double deltaE=max(E_s[i],E_s[j])-max(E_s[k],E_s[l]);
    double deltaEm=min(E_s[i],E_s[j])-min(E_s[k],E_s[l]);
    double x=0.5*(abs(deltaE)+abs(deltaEm));
    int mi=min(min_ij,min_kl);
    int ma=max(max_ij,max_kl);
    double factor=1./(2*ma+1)
                  -2.*mi*(2*mi+1)*(2*mi+2)/3./d_s[i]/d_s[j]/d_s[k]/d_s[l];
    double Tau_1=ee_de*1.7;
    //double f=Tau_1/n_ele/(pow(E_s[i]+E_s[j]-E_s[k]-E_s[l],2)+Tau_1*Tau_1);
    double f=1/Tau_1/n_ele*exp(-pow((E_s[i]+E_s[j]-E_s[k]-E_s[l])/Tau_1,2));
    if(abs(l_s[i]+l_s[j]-l_s[k]-l_s[l])==0)
        return ee_a*exp(-pow((x-ee_b)/ee_c,2))*factor*f;//XXX
    if(abs(l_s[i]+l_s[j]-l_s[k]-l_s[l])==2)
        return ee_a/20*factor*f;//XXX
    return 0;
}

double scatter_cos_m_INT2(int i, int j, int m)
{
    if(abs(state[i].l-state[j].l)!=1)
        return 0;
    if(abs(m)>min(state[i].l,state[j].l))
        return 0;
    double coeff_l_1_lp_000;
    coeff_l_1_lp_000=gsl_sf_coupling_3j(2*state[i].l,2,2*state[j].l,0,0,0);
    double coeffmmp=0;
    coeffmmp=gsl_sf_coupling_3j(2*state[i].l,2,2*state[j].l,2*m,0,-2*m);
    return (2*state[i].l+1)*(2*state[j].l+1)*coeff_l_1_lp_000*coeff_l_1_lp_000
           *coeffmmp*coeffmmp;

}

double scatter_cos_INT2(int i, int j)
{
    if(abs(state[i].l-state[j].l)!=1)
        return 0;
    int l=min(state[i].l,state[j].l);
    int lm=max(state[i].l,state[j].l);
    int m;
    double returnVal=0;
    for(m=-l;m<=l;m++)
        returnVal+=scatter_cos_m_INT2(i,j,m);
    return returnVal/(2*l+1)/(2*lm+1);
}

complex<double> epsilon(int i, int j)
{
    double wev=AU_EV*abs(state[i].energy-state[j].energy);
    return 4.6-9.06*9.06/wev/(wev+I*epsilon_gamma);//XXX
}

complex<double> epsilon(double w)
{
    double wev=w*AU_EV;
    return 4.6-9.06*9.06/wev/(wev+I*epsilon_gamma);//XXX
}

double scatter_rcos_INT2(int i, int j)
{
    // | /+oo            |2
    // | |   jirrjkrr^2dr|
    // | /0              |
    double re=0;
    int ir;
    double scatter_cos=scatter_cos_INT2(i,j);
    if(scatter_cos==0)
        return 0;
    for(ir=0;ir<n_inte;++ir)
            re+=psiR(i,X[ir])*psiR(j,X[ir])*X[ir]*X[ir]*X[ir]*W[ir];
    return re*re*scatter_cos;
}

int ep_selection(int i, int j)
{
    if(abs(state[i].l-state[j].l)!=1)
        return 0;
    return 1;
}

void generator_electron_photon_rate(double wp)
{
    Scatter_Electron_Photon=new double[n_states*n_states];
    Scatter_Exite=new double[n_states*n_states];
    int i,j;
    double w;
    double scatter;
#ifdef OPENMP
    #pragma omp parallel for private(j,w,scatter)
#endif
    for(i=0;i<n_states;++i)
        for(j=0;j<n_states;++j){
            w=abs(state[i].energy-state[j].energy);
            scatter=scatter_rcos_INT2(i,j);
            if(scatter!=0){
                Scatter_Electron_Photon[i*n_states+j]=scatter*pow(w/C_AU,3)*pow(abs(3./(2.+epsilon(i,j))-1.),2);//*pow(abs(epsilon(i,j)),1.5);
                Scatter_Exite[i*n_states+j]=scatter*2*jun::PI*pow(abs(3./(2.+epsilon(wp))),2);
            }
            else{
                Scatter_Electron_Photon[i*n_states+j]=0;
                Scatter_Exite[i*n_states+j]=0;
            }
        }
}

double scatter_ephoton(int i, int j, double temp)
{
    //double re=0;
    if(state[j].energy>EF||state[i].energy<EF)
        return 0;
    if(state[i].energy<=state[j].energy)
        return 0;
    else
        return Scatter_Electron_Photon[i*n_states+j];
}

double scatter_ep(int i, int k, double temp)
{
    double re=ephonon;//XXX
    if(re==0)
        return 0;
    if(abs(state[i].energy-state[k].energy)>0.30000001/AU_EV)
        return 0;
    if(state[i].energy<state[k].energy)
        return 0;
    else
        return re;
}

double exite_e(int i, int k, double temp, double wp)
{
    if(abs(state[i].l-state[k].l)!=1)
        return 0;
    double re=exite*Scatter_Exite[i*n_states+k];
    double Tau_1=ee_de;
    if(state[i].energy<state[k].energy)
        //return re*4/Tau_1*exp(-pow((state[k].energy-state[i].energy-wp)/Tau_1,2));
        return re*4*Tau_1/(pow(state[k].energy-state[i].energy-wp,2)+Tau_1*Tau_1);
    else
        return re*4*Tau_1/(pow(state[i].energy-state[k].energy-wp,2)+Tau_1*Tau_1);
        //return re*4/Tau_1*exp(-pow((state[i].energy-state[k].energy-wp)/Tau_1,2));
}

double total_energy(double p[], int n_p)
{
    int i;
    double energy=0;
    for(i=0;i<n_p;++i)
        energy=energy+state[i].energy*p[i]*(2*state[i].l+1);
    return energy;
}

double evolution(double p[], int n_p, double dt, double t, double temp, double wp)
{
    int i,j,k,l;
    double dp[n_p],dph[n_p];
    double sca_ee=0;
    double sca_eel;
    double sca;
    double return_val=0;
    double pij,pji;
#ifdef OPENMP
    #pragma omp parallel for private(j,k,l,sca,pij,pji,sca_ee,sca_eel)
#endif
    for(i=0;i<n_p;i++){
        dp[i]=0;
        dph[i]=0;
        for(j=0;j<n_p;j++){
            pij=p[i]*p[j];
            pji=(1-p[i])*(1-p[j]);
            for(k=0;k<n_p;k++)
                for(l=k;l<n_p;l++){
                    if(0//abs(E_s[i]+E_s[j]-E_s[k]-E_s[l])>ee_de
                           ||abs(l_s[i]+l_s[j]-l_s[k]-l_s[l])>2)
                        continue;
                    if(abs(l_s[i]+l_s[j]-l_s[k]-l_s[l])==1)
                        continue;
                    if((sca=scatter_ee_INT2(i,j,k,l))>0){
                        if(l!=k)
                            dp[i]+=2*2*sca*(-pij*(1-p[k])*(1-p[l])
                                    +p[k]*p[l]*pji)*d_s[k]*d_s[l]*d_s[j];
                        else
                            dp[i]+=2*sca*(-pij*(1-p[k])*(1-p[l])
                                    +p[k]*p[l]*pji)*d_s[k]*d_s[l]*d_s[j];
                    }
                }
        }
        for(j=0;j<n_p;j++){
            pij=p[i]*(1-p[j]);
            pji=p[j]*(1-p[i]);
            if(ep_selection(i,j)!=0){
                if(ss>0){
                    if(ss%steps==0)
                        sca=exite_e(i,j,temp,wp);
                    else
                        sca=0;
                    sca_ee=scatter_ephoton(i,j,temp)+sca;
                    sca_eel=scatter_ephoton(j,i,temp)+sca;
                }
                else{
                    sca_ee=scatter_ephoton(i,j,temp);
                    sca_eel=scatter_ephoton(j,i,temp);
                }
                dp[i]+=(-sca_ee*pij+sca_eel*pji)*d_s[j];
                //dph[i]+=(-scatter_ep(i,j,temp)*pij+scatter_ep(j,i,temp)*pji)*(2*l_s[j]+1);
            }
            dph[i]+=(-scatter_ep(i,j,temp)*pij+scatter_ep(j,i,temp)*pji)*d_s[j];
        }
        dp[i]+=dph[i];
        dp[i]*=dt;
        dph[i]*=dt;
    }
    if(ss>0)
        ss++;
    for(i=0;i<n_p;++i){
        p[i]+=dp[i];
        if(p[i]>1||p[i]<0||isnan(p[i])){
            cout<<"error! p is out of 0,1 or illegal!"<<endl;
            cout<<"p["<<i<<"]="<<p[i]<<"\tdp="<<dp[i]<<endl;
            exit(1);
        }
        return_val+=dp[i]*dp[i];
    }
    return return_val;
}

void debug_n(double p[])
{
    double n=0;
    int i;
    for(i=0;i<n_states;++i)
        n+=p[i]*d_s[i]*2;
    printf("test # of electrons %f\n",n);
}

void ini_p(double p[], double temp)
{
    int i;
    for(i=0;i<n_states;++i){
        p[i]=fermi_dirac(state[i].energy,temp);
        if(p[i]>1||p[i]<0||isnan(p[i])){
            cout<<"ini error! p is out of 0,1!"<<endl;
            cout<<"ini p["<<i<<"]="<<p[i]<<endl;
            exit(1);
        }
    }
}

void state_2_energy(double p[], double x[], double y[])
{
    int i;
    for(i=0;i<n_states;++i){
        x[i]=state[i].energy;
        y[i]=p[i]*d_s[i]*2;
    }
}

void initial_excitation(double p[], double temp, double wp, double pht, double dt)
{
    ini_p(p,temp);
    int i,k;
    double dp[n_states];
    for(i=0;i<n_states;++i){
        dp[i]=0;
        for(k=0;k<n_states;++k){
            if(state[i].energy<EF&&state[k].energy>EF)
                dp[i]+=(-exite_e(i,k,temp,wp)*p[i]*(1-p[k])*d_s[k]);
            if(state[i].energy>EF&&state[k].energy<EF)
                dp[i]+=exite_e(k,i,temp,wp)*p[k]*(1-p[i])*d_s[k];
        }
    }
    double n=0;
    for(i=0;i<n_states;++i){
        p[i]=p[i]+pht*dp[i]*dt;
        n=n+dp[i]*pht;
        if(p[i]>1||p[i]<0||isnan(p[i])){
            cout<<"ex error! p is out of 0,1!"<<endl;
            cout<<"ex p["<<i<<"]="<<p[i]<<endl;
            exit(1);
        }
    }
}

void hot_carrier(double x[], double p[], double temp, double y[], double yh[], double ye[])
{
    state_2_energy(p,x,y);
    double pfd[n_states];
    double yfd[n_states];
    ini_p(pfd,temp);
    state_2_energy(pfd,x,yfd);
#ifdef OPENMP
    #pragma omp parallel
#endif
    for(int i=0;i<n_states;++i){
        if(x[i]>EF){
            ye[i]=y[i]-yfd[i];
            yh[i]=0;
            //if(ye[i]<-0.01)
                //cout<<"ye<0"<<endl;
        }
        else{
            ye[i]=0;
            yh[i]=yfd[i]-y[i];
            //if(yh[i]<-0.01)
                //cout<<"yh<0"<<endl;
        }
    }
}

double fermi_energy(double temp)
{
    double ef= 0.5*pow(9.*PI/4/r0/r0/r0,2./3);
    epsilon_gamma=sqrt(2*ef)/radius*AU_EV+0.07;
    cout<<"fermi_energy = "<<ef<<endl;
    cout<<"epsilon_gamma = "<<epsilon_gamma<<endl;
    return ef+V0;
}

void delete_cache(void)
{
    delete [] state;
    delete [] COEF_H;
    delete [] COEF_J;
    delete [] Scatter_Electron_Photon;
    delete [] IJ_Photon;
    delete [] X;
    delete [] W;
    delete [] E_s;
}

void save_electron_photon_matrix(void)
{
    ofstream file_photon("matrix_photon");
    int i,j;
    for(i=0;i<n_states;++i)
        for(j=0;j<n_states;++j)
            if(Scatter_Electron_Photon[i*n_states+j]>0)
                file_photon<<i<<"\t"<<state[i].energy*27.211<<"\t"
                           <<j<<"\t"<<state[j].energy*27.211<<"\t"
                           <<(state[i].energy-state[j].energy)*27.211<<"\t"
                           <<Scatter_Electron_Photon[i*n_states+j]<<"\t"
                           <<exite_e(i,j,300,3.5/AU_EV)<<endl;;
    file_photon.close();
}

void save_electron_phonon_matrix(void)
{
    ofstream file_photon("matrix_phonon");
    int i,j;
    for(i=0;i<n_states;++i)
        for(j=0;j<n_states;++j)
            if(scatter_ep(i,j,300)>0)
                file_photon<<i<<"\t"<<j<<"\t"<<scatter_ep(i,j,300)<<endl;;
    file_photon.close();
}

class W_Photon
{
    public:
        double w;
        double n_photon;
};

void ini_IJ_Photon(void)
{
    IJ_Photon=new double[n_states*n_states];
    for(int i=0;i<n_states*n_states;++i)
        IJ_Photon[i]=0;
}

void evolution_photon(double p[], double temp, double dt)
{
    //double d_p[n_states*n_states];
    int i,j;
#ifdef OPENMP
    #pragma omp parallel for private(j)
#endif
    for(i=0;i<n_states;++i)
        for(j=0;j<n_states;++j){
            IJ_Photon[i*n_states+j]=scatter_ephoton(i,j,temp)*p[i]*(1-p[j])*d_s[i]*d_s[j]*2*dt;
            //d_p[i*n_states+j]=Scatter_Electron_Photon[i*n_states+j]*p[i]*(1-p[j])*dt;
            //IJ_Photon[i*n_states+j]=d_p[i*n_states+j];//XXX+=
        }
#ifdef DEBUGx
    ofstream debug_test("wphoton");
    for(i=0;i<n_states;++i)
        for(j=0;j<n_states;++j){
            if(scatter_ephoton(i,j,temp)>0)
            debug_test<<i<<"\t"
                      <<j<<"\t"
                      <<state[i].energy*27.211<<"\t"
                      <<state[j].energy*27.211<<"\t"
                      <<(state[i].energy-state[j].energy)*27.211<<"\t"
                      <<p[i]<<"\t"
                      <<p[j]<<"\t"
                      <<p[i]*(1-p[j])<<"\t"
                      <<scatter_ephoton(i,j,temp)<<"\t"
                      <<IJ_Photon[i*n_states+j]<<endl;
        }
    debug_test.close();
    save("testp",1,n_states,p);
#endif
}

int compare_W_Photon(const void* a, const void* b)
{
    if(((const W_Photon*)a)->w>((const W_Photon*)b)->w)
        return 1;
    if(((const W_Photon*)a)->w<((const W_Photon*)b)->w)
        return -1;
    return 0;
}

void save_w_photon(string filename)
{
    ofstream out(filename.c_str());
    W_Photon wphoton[n_states*n_states];
    int i,j;
    for(i=0;i<n_states;++i)
        for(j=0;j<n_states;++j){
            wphoton[i*n_states+j].w=state[i].energy-state[j].energy;
            wphoton[i*n_states+j].n_photon=IJ_Photon[i*n_states+j];
        }
    qsort(wphoton,n_states*n_states,sizeof(W_Photon),compare_W_Photon);
#ifdef DEBUGx
    ofstream debug_test("wphoton");
    for(i=0;i<n_states*n_states;++i){
        if(wphoton[i].w>0&&wphoton[i].n_photon>0)
            debug_test<<wphoton[i].w*27.211<<"\t"<<wphoton[i].n_photon<<endl;
    }
    debug_test.close();
#endif
    for(i=0;i<n_states*n_states;++i){
        if(wphoton[i].w<=0)
            continue;
        out<<wphoton[i].w<<"\t"<<wphoton[i].n_photon<<endl;
    }
    out.close();
}

void save_convolution(string infilename, string outfilename, double sigma, int num, double xmin, double xmax)
{
    int n;
    n=jun::count_lines(infilename);
    double x[n],y[n];
    ifstream in(infilename.c_str());
    for(int i=0;i<n;++i){
        string line;
        getline(in,line);
        istringstream line_str(line);
        line_str>>x[i]>>y[i];
    }
    in.close();
    double cx[num],cy[num];
    jun::convolution(x,y,cx,cy,xmin,xmax,n,num,sigma);
    ofstream out(outfilename.c_str());
    for(int i=0;i<num;++i)
        out<<cx[i]<<"\t"<<cy[i]<<endl;
    out.close();
}

void save_convolution_photon(string infilename, double sigma, int num)
{
    int n;
    n=jun::count_lines(infilename);
    double x[n],y[n];
    ifstream in(infilename.c_str());
    for(int i=0;i<n;++i){
        string line;
        getline(in,line);
        istringstream line_str(line);
        line_str>>x[i]>>y[i];
    }
    in.close();
    double cx[num],cy[num];
    jun::convolution(x,y,cx,cy,0,abs(V0),n,num,sigma);
    string outfilename=infilename+"conv";
    ofstream out(outfilename.c_str());
    for(int i=0;i<num;++i)
        out<<cx[i]<<"\t"<<cy[i]<<endl;
    out.close();
}

void save_convolution_hotcarrier(string infilename, double sigma, int num)
{
    int n;
    n=jun::count_lines(infilename);
    double x[n],ye[n],yh[n];
    ifstream in(infilename.c_str());
    for(int i=0;i<n;++i){
        string line;
        getline(in,line);
        istringstream line_str(line);
        line_str>>x[i]>>ye[i]>>yh[i];
    }
    in.close();
    double cx[num],cye[num],cyh[num];
    jun::convolution(x,ye,cx,cye,V0,0,n,num,sigma);
    jun::convolution(x,yh,cx,cyh,V0,0,n,num,sigma);
    string outfilename=infilename+"conv";
    ofstream out(outfilename.c_str());
    for(int i=0;i<num;++i){
        if(cx[i]<=EF)
            out<<cx[i]<<"\t"<<0<<"\t"<<cyh[i]<<endl;
        else
            out<<cx[i]<<"\t"<<cye[i]<<"\t"<<0<<endl;
    }
    out.close();
}
