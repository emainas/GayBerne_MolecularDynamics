

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
  const int nu=1; //GB model
  const int mu=2; //GB model
  double dt = 0.001; //Timestep
  double dphi; 
  int steps; //Simulation steps
  double totalTime;
  
  Vec2d pos[N];
  Vec2d vel[N];
  Vec2d force[N];
  Vec2d ort[N];
  Vec2d ortvel[N];
  Vec2d ortforce[N];
  Vec2d ortvelchg[N];
  
  

  double rho;
  double T;
  int M; //Number of cells on each edge. Calculate this for 2D.
  double m = 1.; //Mass
  int I = 1.; //Inertia
  constexpr static double sigma0 = 1.;
  constexpr static double epsilon0 = 1.;
  constexpr static double chi = 0.8;
  constexpr static double chiprime = 0.3819660112501052084965635913249570876359939575195312; //(sqrt(5)-1)/(sqrt(5)+1)
  double L;
  double space;
  double rcut;
  double pot;
  double Q[4];
  
  double KE[1000]; 
  double PE[1000];
  double Total[1000];
  double S[1000];
  

  int initial_steps = 25000; //Initial steps, used to melt the lattice
  int compression_steps;// = 20000; //The number of steps during which to compress the system while keeping T constant
  int cool_steps = 20000; //The number of steps during which to cool the system down to T_target
  int init_rescale_steps = 20000; //The number of steps during which to rescale the temperature to T_target
  int init_equil_steps = 60000; //The number of 'initial' equilibration steps
  int rescale_steps = 20000; //The number of steps which to rescale the temperature
  int equil_steps = 5000; //The number of steps during which to leave the system alone to 'equilibrate'
  int data_steps = 10000; //The number of steps during which to collect data on the system
  
  double init_density = 0.05;
  double final_density;  // = 0.3; //Isotropic Density
  double density_interval=0.00015; //Interval in which to increase density
  double cool_interval;
  double T0 = 5.0; //The starting temperature
  double T_compress = 3.0; //The temperature during compression
  double EL;
  
  std::random_device rd; 
  std::mt19937 generator;
  std::normal_distribution<double> normDistStd;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es;
  
  double order;
  double order6;
  double order4;

  double chi1;
  

  Eigen::Vector2d eigenv;
  Vec2d eigenvv;
  double phi[N];
  
  void _savedata(const Vec2d* var, const std::string& fn) const;
  void _loaddata(Vec2d* var, const std::string& fn);
  
public:
  
  GB(double rho, double T):rho(rho),T(T),M(round(pow((N/2),1./2.))),generator(rd()), normDistStd(std::normal_distribution<>(0.0, 1))  {
    final_density = rho;
    L = pow(N/init_density,1.0/2.0);
    space = L/M;
    rcut = L/2;
    compression_steps = ceil((final_density - init_density)/density_interval*100);
    steps = initial_steps + compression_steps + cool_steps + init_rescale_steps + init_equil_steps + rescale_steps + equil_steps + data_steps;
    //    density_interval = (final_density - init_density)/compression_steps*100;
    cool_interval = (T_compress-T)/cool_steps*1000;
    totalTime = steps*dt;
  };
  GB(double rho, const Vec2d* position, const Vec2d* orientation):rho(rho),M(round(pow((N/2.),1./2.))), L(pow(N/rho, 1.0/2.0)), space(L/M), rcut(L/2){
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
  double potentialEnergy2(double V0) const; 
  void calc_force();
  void calc_force2(double V0);
  void calc_order();
  void calc_order6();
  void calc_order4();
  void VV();
  void VV2(double V0);
  void temp_control(double T_target);
  double run(int is_init_file);
  void g_corr();
  void g_orient();
  void eos();
  double f(double rscale, float thetaj, float thetak) const;
  void integration();

};

#endif /* MD_GB_hpp */
