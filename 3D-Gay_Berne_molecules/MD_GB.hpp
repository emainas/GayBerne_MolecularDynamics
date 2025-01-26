
#ifndef MD_GB_hpp
#define MD_GB_hpp

#include <stdio.h>
#include <string>
#include <memory>
#include <math.h>
#include <istream>
#include <ostream>
#include <iostream>
#include <Eigen/Dense>
#include "basics.h"


class GB {
  /*Assume k=1, m=1*/
protected:
  constexpr static int N = N_const;
  const int nu=1;
  const int mu=2;
  double dt = 0.001;
  double dphi;
  int steps;
  double totalTime;
  Vec3d pos[N];
  Vec3d vol[N];
  Vec3d force[N];
  Vec3d ort[N];
  Vec3d ortvol[N];
  Vec3d ortforce[N];
  Vec3d ortvolchg[N];
  Int3d periods[N];
  double rho;
  double T;
  int M; // number of cells on each edge
  double m = 1.;
  int I = 1.;
  constexpr static double sigma0 = 1.;
  constexpr static double epsilon0 = 1.;
  constexpr static double chi = 0.8;
  constexpr static double chiprime = 0.3819660112501052084965635913249570876359939575195312;
  double L;
  double space;
  double rcut;
  double pot;
  double Q[6];
  
  double KE[1000000];
  double PE[1000000];
  double S[1000000];
  double C[1000000];
  
  int initial_steps = 300000; //Initial steps, used to melt the lattice
  int compression_steps;// = 20000; //The number of steps during which to compress the system while keeping T constant
  int cool_steps = 300000; //The number of steps during which to cool the system down to T_target
  int init_rescale_steps = 300000; //The number of steps during which to rescale the temperature to T_target
  int init_equil_steps = 600000; //The number of 'initial' equilibration steps
  int rescale_steps = 300000; //The number of steps which to rescale the temperature
  int equil_steps = 3000000; //The number of steps during which to leave the system alone to 'equilibrate'
  int data_steps = 10000000; //The number of steps during which to collect data on the system
  
  double init_density = 0.05;
  double final_density;// = 0.3; //Isotropic Density
  double density_interval=0.001; //Interval in which to increase density
  double cool_interval;
  double T0 = 5.0; //The starting temperature
  double T_compress = 3.0; //The temperature during compression
  double EL;
  std::random_device rd;
  int seed = 4; // Use this instead of the above random device for reproducible results
  std::mt19937 generator;
  std::normal_distribution<double> normDistStd;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  double order;
  double chi1;
  Eigen::Vector3d eigenv;
  void _savedata(const Vec3d* var, const std::string& fn) const;
  void _loaddata(Vec3d* var, const std::string& fn);
  clock_t t1,t2,t3,t4, start, stop;
  
  
  
  
  
public:
  GB(double rho, double T):rho(rho),T(T),M(round(pow((N/4.),1./3.))),generator(seed), normDistStd(std::normal_distribution<>(0.0, 1))  {
    final_density = rho;
    L = pow(N/init_density,1.0/3.0);
    space = L/M;
    rcut = L/2;
    compression_steps = ceil((final_density - init_density)/density_interval*100);
    steps = initial_steps + compression_steps + cool_steps + init_rescale_steps + init_equil_steps + rescale_steps + equil_steps + data_steps;
    //    density_interval = (final_density - init_density)/compression_steps*100;
    cool_interval = (T_compress-T)/cool_steps*1000;
    totalTime = steps*dt;
  };
  GB(double rho, const Vec3d* position, const Vec3d* orientation):rho(rho),M(round(pow((N/4.),1./3.))), L(pow(N/rho, 1.0/3.0)), space(L/M), rcut(L/2){
    for (int i=0; i<N; i++) {
      pos[i] = position[i];
      ort[i] = orientation[i];
    }
  };
  void initialisation();
  void outputit() const;
  void savedata() const;
  void loaddata();
  double kineticEnergy() const;
  double rotKE() const;
  double transKE() const;
  double potentialEnergy() const;
  void calc_force();
  void calc_order();
  void susceptibility();
  void VV();
  void temp_control(double T_target);
  double run(int is_init_file);
  void g_corr();
};

#endif /* MD_GB_hpp */
