#include <iostream>
#include <complex>
#include <cmath>
#include <string>
#include <sstream>
#include <omp.h>
#include <ctime>
#include "states.h"
#include "const.h"
#include "readin.hpp"
#include "save.h"
#include <cstdlib>

using namespace std;
using jun::PI;
using jun::EC;
using jun::A0_NM;
using jun::AU_EV;
using jun::read;
using jun::save;

double ex_val;
double ee_b;
double ee_c;
double ee_a;
double ee_de;
double V0;
int n_states;
int n_ele;
double r0;
double density;        //in g*nm-3
double radius;
double EF;                      //fermi level
int ss;                      //fermi level
int steps;
double wp;            //excitation
double pht;
double exite;
double ephonon;
char dir[256];
State *state;

int main(int argc, char* argv[]){
    time_t start=time(NULL);
    if(argc!=2){
        cout<<"No input file!"<<endl;
        std::exit(0);
    }
//*************************************************************************
    double temp;
    double dt;
    double et;
    double t=0;
    string filename="data";
    int num=0;
    char ff[256];
//*************************************************************************
    read(V0,argv[1],"V0");
    read(n_ele,argv[1],"n_ele");
    read(radius,argv[1],"radius");
    read(r0,argv[1],"r0");
    read(density,argv[1],"density");
    read(wp,argv[1],"wp");
    read(pht,argv[1],"pht");
    read(temp,argv[1],"temp");
    read(dt,argv[1],"dt");
    read(et,argv[1],"et");
    read(exite,argv[1],"exite");
    read(ephonon,argv[1],"ephonon");
    read(ss,argv[1],"ss");
    read(ee_c,argv[1],"ee_c");
    read(steps,argv[1],"steps");
    read(ee_de,argv[1],"ee_de");
//*************************************************************************
    if(radius==0)
        radius=r0*pow(n_ele*1.,1./3);
    else{
        radius=radius/A0_NM;
        n_ele=int(pow(radius/r0,3));
    }
    if(ss==0)
        sprintf(dir,"de_%d_%.2f_%.1e",n_ele,wp,ee_de);
    else
        sprintf(dir,"ss_%d_%.2f_%d_%.3e_%.1e",n_ele,wp,steps,exite,ee_de);
    sprintf(ff,"mkdir %s",dir);
    system(ff);
    sprintf(ff,"rm %s/*",dir);
    system(ff);

    V0=V0/AU_EV;
    wp=wp/AU_EV;
    EF=fermi_energy(temp);
    ee_a=57.97e-4/n_ele;
    ee_b=(0.01495*n_ele+7.336)/(n_ele+106.6);

    ini_gauss(300);
    jun::CONST_dfactorial_generator(250);
    sprintf(ff,"%s/energy.txt",dir);
    n_states=get_free_energies(ff);
    ini_n_l_E(ff);
    normalize_coeff();
    generator_electron_photon_rate(wp);
    save_electron_photon_matrix();
//*************************************************************************
    double p[n_states];
    double x[n_states],y[n_states],ye[n_states],yh[n_states];
//*************************************************************************
#ifdef OPENMP
    int n_core;
    n_core=omp_get_num_procs();
    cout<<"==================================================="<<endl;
    omp_set_num_threads(n_core);
    cout<<n_core<<" processors will be used."<<endl;
#endif
    cout<<"radius = "<<radius<<endl;
    cout<<"ee_a = "<<ee_a<<endl;
    cout<<"ee_b = "<<ee_b<<endl;
    cout<<"ee_c = "<<ee_c<<endl;
    cout<<"ee_de = "<<ee_de<<endl;
//*************************************************************************
#ifdef DEBUG
    cout<<"==================================================="<<endl;
    ini_p(p,temp);
    debug_n(p);
    cout<<"==================================================="<<endl;
#endif
    cout<<"n_states = "<<n_states<<endl;

    initial_excitation(p,temp,wp,pht,dt);
//*************************************************************************
    ini_IJ_Photon();

    sprintf(ff,"%s/tandT",dir);
    ofstream tandT(ff);

    hot_carrier(x,p,temp,y,yh,ye);

    sprintf(ff,"%s/p0",dir);
    filename=ff;
    save(filename,1,n_states,p);

    sprintf(ff,"%s/data0",dir);
    save(ff,3,n_states,x,ye,yh);
    save_convolution_hotcarrier(ff,abs(V0)/100,1000);

    sprintf(ff,"%s/photon0",dir);
    save_w_photon(ff);
    save_convolution_photon(ff,abs(V0)/10,1000);

    tandT<<num<<"\t"<<t<<"\t"<<temp<<"\t"<<EF<<endl;
    cout<<"evolution start!"<<endl;

#ifdef DEBUG_TIME
        cout<<"time cost:"<<difftime(time(NULL),start)<<"s"<<endl;
#endif
    double sdp=0;
    double nnh=0,nne=0;
    for(t=dt;t<et;t+=dt){
        sdp=evolution(p,n_states,dt,t,temp,wp);//XXX
        ++num;
        int nnn=500;
        if(ss>0){
            nnn=50;
        }
        if(num%nnn==0){
            hot_carrier(x,p,temp,y,yh,ye);
            evolution_photon(p,temp,dt);

            sprintf(ff,"%s/p%d",dir,num/nnn);
            filename=ff;
            save(filename,1,n_states,p);

            //sprintf(ff,"%s/dos%d",dir,num/nnn);
            //filename=ff;
            //save(filename,2,n_states,x,y);

            cout<<"total energy = "<<total_energy(p,n_states)<<endl;;
            nnh=nne=0;
            for(int ii=0;ii<n_states;ii++){
                nnh+=yh[ii];
                nne+=ye[ii];
            }
            printf("nnh = %f, nne = %f   delta = %e sdp = %e ex_val = %e\n",nnh,nne,nnh-nne,sdp,ex_val);

            debug_n(p);

            sprintf(ff,"%s/data%d",dir,num/nnn);
            filename=ff;
            save(filename,3,n_states,x,ye,yh);
            save_convolution_hotcarrier(filename,abs(V0)/100,1000);

            sprintf(ff,"%s/photon%d",dir,num/nnn);
            filename=ff;
            save_w_photon(filename);
            save_convolution_photon(filename,abs(V0)/50,1000);

            tandT<<num<<"\t"<<t<<"\t"<<temp<<"\t"<<EF<<endl;
            cout<<"progress "<<t/et<<endl;
            cout<<"temperature = "<<temp<<"K"<<endl;
#ifdef DEBUG_TIME
            cout<<"time cost:"<<difftime(time(NULL),start)<<"s"<<endl;
#endif
        }
    }
    cout<<"Final temperature = "<<temp<<"K"<<endl;
    tandT.close();
#ifdef DEBUG_TIME
        cout<<"time cost:"<<difftime(time(NULL),start)<<"s"<<endl;
#endif
    delete_cache();
    return 0;
}
