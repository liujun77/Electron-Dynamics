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
#ifdef DEBUG_TIME
extern time_t start;
#endif
extern int select_status;
extern int reset;
extern double V0;
extern int n_states;
extern int n_energy;
extern int n_ele;
extern double r0;
extern double radius;
extern double EF;
extern int ss;
extern double exite;
extern double ephonon;
extern double ee_rate;
extern double ee_alpha;
extern double ee_de;
extern int epsilon_type;
extern State *state;
extern char dir[256];

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
int *l_s;
int lmax;
vector<Selection> selection;
long n_select;

double fermi_dirac(double energy, double temp)
{
    return 1/(exp((energy-EF)/Boltzman_k/temp)+1);
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
    int l=-1,nnl=1,m;
    int i=0;
    double y1,y2,energy,root;
    double de=1e-5;
    n_energy=0;
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
                    ++n_energy;
                    for(m=-l;m<=l;++m){
                        filein<<l<<"\t"<<nnl<<"\t"<<m<<"\t"<<-0.5<<"\t"<<root<<endl;
                        filein<<l<<"\t"<<nnl<<"\t"<<m<<"\t"<<0.5<<"\t"<<root<<endl;
                        i+=2;
                    }
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
        da=100*((const State*)a)->m+((const State*)a)->spin;
        db=100*((const State*)b)->m+((const State*)b)->spin;
        if(da>db)
            return 1;
        else
            return -1;
    }
    return 0;
}

void ini_n_l_m_spin_E(string filename)
{
    int i;
    string line;
    state=new State [n_states];
    ifstream data_file(filename.c_str());
    for(i=0;i<n_states;++i){
        getline(data_file,line);
        istringstream line_str(line);
        line_str>>state[i].l>>state[i].n>>state[i].m>>state[i].spin>>state[i].energy;
    }
    data_file.close();
    qsort(state,n_states,sizeof(State),compare_State);
    double energy_temp=0;
    int j=0;
    E_s=new double [n_energy];
    l_s=new int [n_energy];
    for(i=0;i<n_states;i++){
        if(state[i].energy!=energy_temp){
            energy_temp=state[i].energy;
            E_s[j]=energy_temp;
            l_s[j]=state[i].l;
            state[i].Eind=j;
            j++;
        }
        else
            state[i].Eind=j-1;
    }
    if(j!=n_energy){
        cout<<"error"<<endl;
        exit(1);
    }
#ifdef DEBUGx
    ofstream debug_test("states");
    for(i=0;i<n_states;i++)
        debug_test<<state[i].l<<"\t"<<state[i].energy<<"\t"<<state[i].Eind<<"\t"<<state[i].m<<"\t"<<state[i].spin<<endl;
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
        if((i>1)&&(state[i].energy==state[i-1].energy)){
            COEF_J[i]=COEF_J[i-1];
            COEF_H[i]=COEF_H[i-1];
            continue;
        }
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

double r_rp(double r, double rp, int l)
{
    if(r>rp)
        return pow(rp,l)/pow(r,l+1);
    else
        return pow(r,l)/pow(rp,l+1);
}

void ini_gauss(int n_inte_gauss)
{
    n_inte=n_inte_gauss;
    X=new double[n_inte];
    W=new double[n_inte];
    legendre_compute_glr(n_inte,X,W);
    rescale(0,1.2*radius,n_inte,X,W);
}

double scatter_ee_INT2(int ei, int ej, double spin, int mij, int k, int l)
{
    if(abs(E_s[ei]+E_s[ej]-state[k].energy-state[l].energy)>ee_de)
        return 0;
    if(mij!=state[k].m+state[l].m)
        return 0;
    if(spin!=state[k].spin+state[l].spin)
        return 0;
    double deltaE=max(E_s[ei],E_s[ej])-max(state[k].energy,state[l].energy);
    return ee_rate*exp(-abs(deltaE)/ee_alpha);
}

double scatter_ee_INT2(int i, int j, int k, int l)
{
    if(abs(state[i].energy+state[j].energy-state[k].energy-state[l].energy)>ee_de)
        return 0;
    if(state[i].m+state[j].m!=state[k].m+state[l].m)
        return 0;
    if(state[i].spin+state[j].spin!=state[k].spin+state[l].spin)
        return 0;
    double deltaE=max(state[i].energy,state[j].energy)-max(state[k].energy,state[l].energy);
    return ee_rate*exp(-abs(deltaE)/ee_alpha);
}

double scatter_cos_INT2(int i, int j)
{
    if(abs(state[i].l-state[j].l)!=1)
        return 0;
    if(state[i].m!=state[j].m)
        return 0;
    if(state[i].spin!=state[j].spin)
        return 0;
    double coeff_l_1_lp_000;
    coeff_l_1_lp_000=gsl_sf_coupling_3j(2*state[i].l,2,2*state[j].l,0,0,0);
    double coeffmmp=0;
    coeffmmp=gsl_sf_coupling_3j(2*state[i].l,2,2*state[j].l,2*state[i].m,0,-2*state[j].m);
    return (2*state[i].l+1)*(2*state[j].l+1)*coeff_l_1_lp_000*coeff_l_1_lp_000
           *coeffmmp*coeffmmp;

}

complex<double> epsilon(int i, int j)
{
    int ei=state[i].Eind;
    int ej=state[j].Eind;
    if(epsilon_type>1)
        return Epsilon_r[ei*n_energy+ej]+I*Epsilon_i[ei*n_energy+ej];
    else{
        double wev=AU_EV*abs(state[i].energy-state[j].energy);
        return 4.18-9.06*9.06/wev/(wev+I*0.66);
    }
}

complex<double> epsilon(double w)
{
    double wev=w*AU_EV;
    return 4.18-9.06*9.06/wev/(wev+I*0.66);
}

double scatter_rcos_INT2(int i, int j)
{
    // | /+oo            |2
    // | |   jirrjkrr^2dr|
    // | /0              |
    double re=0;
    int ir;
    double scatter_cos=scatter_cos_INT2(i,j);
    //if(state[i].energy=state[j].energy)
    //    return 0;
    if(scatter_cos==0)
        return 0;
    for(ir=0;ir<n_inte;++ir)
            re+=psiR(i,X[ir])*psiR(j,X[ir])*X[ir]*X[ir]*X[ir]*W[ir];
    //double w=abs(state[i].energy-state[j].energy);
    return re*re*scatter_cos;//*pow(w/C_AU,3)*pow(abs(3./(2.+epsilon(w))),2)*pow(abs(3.*epsilon(w)/(2.*epsilon(w)+1.)),2);
}

int ep_selection(int i, int j)
{
    if(abs(state[i].l-state[j].l)!=1)
        return 0;
    if(state[i].m!=state[j].m)
        return 0;
    if(state[i].spin!=state[j].spin)
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
    //cout<<"gener begin!"<<endl;
#ifdef OPENMP
    #pragma omp parallel for private(j,w,scatter)
#endif
    for(i=0;i<n_states;++i)
        for(j=0;j<n_states;++j){
            w=abs(state[i].energy-state[j].energy);
            //cout<<"gener#!i="<<i<<"\tj="<<j<<endl;
            scatter=scatter_rcos_INT2(i,j);
            if(scatter!=0){
                Scatter_Electron_Photon[i*n_states+j]=scatter*pow(w/C_AU,3)*pow(abs(3./(2.+epsilon(i,j))),2)*pow(abs(3.*epsilon(i,j)/(2.*epsilon(i,j)+1.)),2);
                Scatter_Exite[i*n_states+j]=scatter*2*jun::PI*pow(abs(3./(2.+epsilon(wp))),2);
            }
            else{
                Scatter_Electron_Photon[i*n_states+j]=0;
                Scatter_Exite[i*n_states+j]=0;
            }
        }
    //cout<<"gener done!"<<endl;
}

double scatter_ephoton(int i, int j, double temp)
{
    //double re=0;
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
    if(abs(state[i].energy-state[k].energy)>0.10000001/AU_EV)
        return 0;
    if(state[i].energy<state[k].energy)
        re=re*exp((state[i].energy-state[k].energy)/temp/Boltzman_k);
    else
        re=re;
    return re;
}

double exite_e(int i, int k, double temp, double wp)
{
    if(abs(state[i].l-state[k].l)!=1)
        return 0;
    if(state[i].m!=state[k].m)
        return 0;
    if(state[i].spin!=state[k].spin)
        return 0;
    double re=exite*Scatter_Exite[i*n_states+k];
    if(state[i].energy<state[k].energy)
        return re//*exp((state[i].energy-EF)*lamba)
               *exp(-50000*pow((state[k].energy-state[i].energy-wp),2));
    else
        return re//*exp((state[i].energy-EF)*lamba)
               *exp(-50000*pow((state[i].energy-state[k].energy-wp),2));
}

double total_energy(double p[], int n_p)
{
    int i;
    double energy=0;
    for(i=0;i<n_p;++i)
        energy=energy+state[i].energy*p[i];
    return energy;
}

void generate_eij(void)
{
    int ei,ej;
    char ff[256];
    sprintf(ff,"%s/energy_ij",dir);
    ofstream fileout(ff);
    for(ei=0;ei<n_energy;ei++)
        for(ej=ei+1;ej<n_energy;ej++)
            fileout<<ei<<"\t"<<ej<<"\t"<<abs(E_s[ei]-E_s[ej])<<endl;
    fileout.close();
}

void ini_epsilon(void)
{
    generate_eij();
    char ff[256];
    sprintf(ff,"nohup matlab -nodisplay -r \"epsilon(\'%s\');quit\" run.out",dir);
    system(ff);
    int ei,ej;
    double waste;
    Epsilon_r=new double[n_energy*n_energy];
    Epsilon_i=new double[n_energy*n_energy];
    for(ei=0;ei<n_energy;ei++)
        for(ej=0;ej<n_energy;ej++){
            Epsilon_r[ei*n_energy+ej]=0;
            Epsilon_i[ei*n_energy+ej]=0;
        }
    sprintf(ff,"%s/epsilon_ij",dir);
    ifstream filein(ff);
    if(epsilon_type==2)
        while(!filein.eof()){
            filein>>ei>>ej>>waste;
            filein>>Epsilon_r[ei*n_energy+ej]>>Epsilon_i[ei*n_energy+ej]>>waste>>waste;
            Epsilon_r[ej*n_energy+ei]=Epsilon_r[ei*n_energy+ej];
            Epsilon_i[ej*n_energy+ei]=Epsilon_i[ei*n_energy+ej];
        }
    else if(epsilon_type==3)
        while(!filein.eof()){
            filein>>ei>>ej>>waste;
            filein>>waste>>waste>>Epsilon_r[ei*n_energy+ej]>>Epsilon_i[ei*n_energy+ej];
            Epsilon_r[ej*n_energy+ei]=Epsilon_r[ei*n_energy+ej];
            Epsilon_i[ej*n_energy+ei]=Epsilon_i[ei*n_energy+ej];
        }
    filein.close();
    cout<<"Epsilon done!"<<endl;
}

void ini_selection(void)
{
    //ofstream file("eij.txt");
    //for(int i=0;i<n_states;i++)
    //    file<<state[i].l<<"\t"<<state[i].energy<<"\t"<<state[i].Eind<<"\t"<<state[i].m<<endl;
    //file.close();
    int ei,ej,mij,k,l;
    n_select=0;
    double sca_ee;
    char ff[n_energy][256];
    for(int i=0;i<n_energy;i++)
        sprintf(ff[i],"selection/selection_ee%d.dat",i);//XXX
    int n[n_energy];
    cout<<"select_status = "<<select_status<<endl;
    cout<<"reset = "<<reset <<endl;
    if(select_status!=0||reset==1){
        ofstream selection_file;
#ifdef OPENMP
        #pragma omp parallel for private(ej,mij,k,l,sca_ee,selection_file)
#endif
        for(ei=0;ei<n_energy;ei++){
            n[ei]=0;
            selection_file.open(ff[ei]);
            selection_file<<"                \n";
            //cout<<ei<<endl;
            for(ej=ei;ej<n_energy;ej++)
                for(mij=0;mij<=l_s[ei]+l_s[ej];mij++)
                    for(k=0;k<n_states;k=k+2)
                        for(l=k;l<n_states;l=l+2){
                            sca_ee=scatter_ee_INT2(ei,ej,-1,mij,k,l);
                            if(sca_ee>1e-8){
                                selection_file<<ei<<"\t"<<ej<<"\t"
                                    <<mij<<"\t"<<k<<"\t"<<l<<endl;
                                n[ei]++;
                            }
                        }
            selection_file.seekp(ios::beg);
            selection_file<<n[ei];;
            selection_file.close();
        }
    }
    for(int i=0;i<n_energy;i++){
        ifstream selection_file_in(ff[i]);
        selection_file_in>>n[i];
        for(int j=0;j<n[i];j++){
            Selection selec;
            selection_file_in>>selec.Ei
                             >>selec.Ej
                             >>selec.mij
                             >>selec.k
                             >>selec.l;
            selection.push_back(selec);
            n_select++;
        }
        selection_file_in.close();
    }
    cout<<"n_select = "<<n_select<<endl;
}

void debug_fij(int ei, int ej, int mij, double p[], double sf[])
{
    int k,l;
    double f,s;
    f=s=0;
    for(k=0;k<n_states;k+=2)
        for(l=0;l<n_states;l+=2){
            f+=scatter_ee_INT2(ei,ej,-1,mij,k,l)*(1-p[k])*(1-p[l]);
            s+=scatter_ee_INT2(ei,ej,-1,mij,k,l)*p[k]*p[l];
        }
    //cout<<"fij=\t"<<f<<"\tsij=\t"<<s<<endl;
    sf[0]=f;
    sf[1]=s;
}

void debug_fij(int i, int j, double p[], double sf[])
{
    int k,l;
    double f,s;
    f=s=0;
    for(k=0;k<n_states;k+=2)
        for(l=0;l<n_states;l+=2){
            f+=scatter_ee_INT2(i,j,k,l)*(1-p[k])*(1-p[l]);
            s+=scatter_ee_INT2(k,l,i,j)*p[k]*p[l];
        }
    //cout<<"fij=\t"<<f<<"\tsij=\t"<<s<<endl;
    sf[0]=f;
    sf[1]=s;
}

double evolution(double p[], int n_p, double dt, double t, double temp, double wp)
{
    int i,j,k,l;
    int n;
    int ei,ej;
    int mij;
    double dp[n_p],dph[n_p];
    double fij[n_energy][n_energy][2*lmax+1];
    double sij[n_energy][n_energy][2*lmax+1];
    double sca_ee=0;
    double sca_eel;
    double sca;
    double return_val=0;
    double pij,pji;
#ifdef OPENMP
    #pragma omp parallel for private(ej,mij)
#endif
    for(ei=0;ei<n_energy;ei++)
        for(ej=0;ej<n_energy;ej++)
            for(mij=0;mij<=l_s[ei]+l_s[ej];mij++){
                fij[ei][ej][mij]=0;
                sij[ei][ej][mij]=0;
            }
#ifdef OPENMP
    #pragma omp parallel for private(ei,ej,mij,k,l,sca_ee)
#endif
    for(n=0;n<n_select;n++){
        ei=selection[n].Ei;
        ej=selection[n].Ej;
        mij=selection[n].mij;
        k=selection[n].k;
        l=selection[n].l;
        sca_ee=scatter_ee_INT2(ei,ej,-1,mij,k,l);
        if(k!=l){
            fij[ei][ej][mij]+=(2*sca_ee*(1-p[k])*(1-p[l]));
            sij[ei][ej][mij]+=(2*sca_ee*p[k]*p[l]);
        }
        else{
            fij[ei][ej][mij]+=(sca_ee*(1-p[k])*(1-p[l]));
            sij[ei][ej][mij]+=(sca_ee*p[k]*p[l]);
        }
        /*
        if(ej==ei)
            continue;
        if(k!=l){
            fij[ej][ei][mij]+=(2*sca_ee*(1-p[k])*(1-p[l]));
            sij[ej][ei][mij]+=(2*sca_ee*p[k]*p[l]);
        }
        else{
            fij[ej][ei][mij]+=(sca_ee*(1-p[k])*(1-p[l]));
            sij[ej][ei][mij]+=(sca_ee*p[k]*p[l]);
        }
        */
    }

#ifdef OPENMP
    #pragma omp parallel for private(ei,mij)
#endif
    for(ej=0;ej<n_energy;ej++)
        for(ei=ej+1;ei<n_energy;ei++)
            for(mij=0;mij<2*lmax+1;mij++){
                fij[ei][ej][mij]=fij[ej][ei][mij];
                sij[ei][ej][mij]=sij[ej][ei][mij];
            }

#ifdef OPENMP
    #pragma omp parallel for private(j,k,l,mij,sca,sca_ee,sca_eel,pij,pji)
#endif
    for(i=0;i<n_p;i=i+2){
        dp[i]=0;
        dph[i]=0;
        for(j=0;j<n_p;j=j+2){
            if(1){
                mij=abs(state[i].m+state[j].m);
                dp[i]+=(2*(-fij[state[i].Eind][state[j].Eind][mij]*p[i]*p[j]
                       +sij[state[i].Eind][state[j].Eind][mij]*(1-p[i])*(1-p[j])));
            }
        }
        for(j=0;j<n_p;j=j+2){
            pij=p[i]*(1-p[j]);
            pji=p[j]*(1-p[i]);
            if(ep_selection(i,j)!=0){
                if(ss==1){
                    sca=exite_e(i,j,temp,wp);
                    sca_ee=scatter_ephoton(i,j,temp)+sca;
                    sca_eel=scatter_ephoton(j,i,temp)+sca;
                }
                else{
                    sca_ee=scatter_ephoton(i,j,temp);
                    sca_eel=scatter_ephoton(j,i,temp);
                }
                dp[i]+=(-sca_ee*pij+sca_eel*pji);
            }

            dph[i]+=(-scatter_ep(i,j,temp)*pij+scatter_ep(j,i,temp)*pji);

        }
        dp[i]+=dph[i];
        dp[i]*=dt;
        dph[i]*=dt;
        dp[i+1]=dp[i];
        dph[i+1]=dph[i];
    }
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
        n+=p[i];
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
    int i,j;
    x[0]=state[0].energy;
    y[0]=p[0];
    j=0;
    for(i=1;i<n_states;++i){
        if(state[i].energy==x[j])
            y[j]+=p[i];
        else{
            ++j;
            x[j]=state[i].energy;
            y[j]=p[i];
        }
    }
}

void initial_excitation(double p[], double temp, double wp, double pht)
{
    ini_p(p,temp);
    int i,k;
    double dp[n_states];
    for(i=0;i<n_states;++i){
        dp[i]=0;
        for(k=0;k<n_states;++k)
            dp[i]+=(-exite_e(i,k,temp,wp)*p[i]*(1-p[k])
                +exite_e(k,i,temp,wp)*p[k]*(1-p[i]));
    }
    double n=0;
    for(i=0;i<n_states;++i){
        p[i]=p[i]+pht*dp[i];
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
    double yfd[n_energy];
    ini_p(pfd,temp);
    state_2_energy(pfd,x,yfd);
#ifdef OPENMP
    #pragma omp parallel
#endif
    for(int i=0;i<n_energy;++i){
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

double delta_num_ini_state(double fermi, double temp)
{
    double p[n_states];
    int i;
    for(i=0;i<n_states;++i)
        p[i]=1/(exp((state[i].energy-fermi)/Boltzman_k/temp)+1);
    double n=0;
    for(i=0;i<n_states;++i)
        n=n+p[i];
    return n-n_ele;
}

double fermi_energy(double temp)
{
    double mid,fx1,fx2,fmid;
    double x1=-0.05;
    double x2=V0+0.05;
    do{
        mid=(x1+x2)/2;
        fx1=delta_num_ini_state(x1,temp);
        fx2=delta_num_ini_state(x2,temp);
        fmid=delta_num_ini_state(mid,temp);
        if(fx1*fx2>0)
            return 1;   //avoid the divergence case means there is no root
        else{
            if(fx1*fmid>0)
                x1=mid;
            else x2=mid;
        }
    }while(fabs(fmid)>1e-10);
    return mid;
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
                file_photon<<i<<"\t"<<j<<"\t"
                           <<state[i].energy-state[j].energy<<"\t"
                           <<Scatter_Electron_Photon[i*n_states+j]<<endl;;
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
            IJ_Photon[i*n_states+j]=Scatter_Electron_Photon[i*n_states+j]*p[i]*(1-p[j])*dt;
            //d_p[i*n_states+j]=Scatter_Electron_Photon[i*n_states+j]*p[i]*(1-p[j])*dt;
            //IJ_Photon[i*n_states+j]=d_p[i*n_states+j];//XXX+=
        }
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
    double wave=0;
    double n=0;
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
        if(wave!=wphoton[i].w){
            if(i>0&&n>0)//XXX
                out<<wphoton[i-1].w<<"\t"<<n<<endl;
            n=0;
            wave=wphoton[i].w;
            n=n+wphoton[i].n_photon;
        }
        else{
            n=n+wphoton[i].n_photon;
        }
    }
    if(i==n_states*n_states)//XXX
        out<<wphoton[i-1].w<<"\t"<<n<<endl;
    out.close();
}

void save_convolution(string infilename, string outfilename, double sigma, int ntimes, double xmin, double xmax)
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
    double cx[ntimes*n],cy[ntimes*n];
    jun::convolution(x,y,cx,cy,xmin,xmax,n,ntimes*n,sigma);
    ofstream out(outfilename.c_str());
    for(int i=0;i<ntimes*n;++i)
        out<<cx[i]<<"\t"<<cy[i]<<endl;
    out.close();
}

void save_convolution_photon(string infilename, double sigma, int ntimes)
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
    double cx[ntimes*n],cy[ntimes*n];
    jun::convolution(x,y,cx,cy,0,abs(V0),n,ntimes*n,sigma);
    string outfilename=infilename+"conv";
    ofstream out(outfilename.c_str());
    for(int i=0;i<ntimes*n;++i)
        out<<cx[i]<<"\t"<<cy[i]<<endl;
    out.close();
}

void save_convolution_hotcarrier(string infilename, double sigma, int ntimes)
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
    double cx[ntimes*n],cye[ntimes*n],cyh[ntimes*n];
    jun::convolution(x,ye,cx,cye,V0,0,n,ntimes*n,sigma);
    jun::convolution(x,yh,cx,cyh,V0,0,n,ntimes*n,sigma);
    string outfilename=infilename+"conv";
    ofstream out(outfilename.c_str());
    for(int i=0;i<ntimes*n;++i){
        if(cx[i]<=EF)
            out<<cx[i]<<"\t"<<0<<"\t"<<cyh[i]<<endl;
        else
            out<<cx[i]<<"\t"<<cye[i]<<"\t"<<0<<endl;
    }
    out.close();
}
