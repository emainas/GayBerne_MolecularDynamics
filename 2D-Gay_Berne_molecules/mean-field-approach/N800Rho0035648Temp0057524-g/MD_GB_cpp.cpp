

#include "MD_GB_hpp.hpp"
#include <iostream>
#include <math.h>
#include <random>
#include <fstream>
#include <time.h>
#include <iomanip>

using namespace std;

void GB::outputit() const{
  cout<<"Density: "<<rho<<"\nTemperature:"<<T<<endl;
  cout<<"Initial box length: "<<L<<"\nFinal box length: " << pow(N/final_density, 1.0/2.0)<<endl;
}

void GB::savedata() const{
  _savedata(pos, "pos");
  _savedata(vel, "vel");
}

void GB::_savedata(const Vec2d* var, const string& fn) const{
  ofstream out("./N800Rho0035648Temp0057524-g/GB_"+fn+".txt");
  for (int i=0; i<N; i++) {
    out << var[i];
  }
  out.close();
}



void GB::loaddata(){
  _loaddata(pos, "pos");
  _loaddata(vel, "vel");
  L =  pow(N/final_density, 1.0/2.0);
  calc_force();
}

void GB::_loaddata(Vec2d* var, const string& fn) {
  ifstream in("./N800Rho0035648Temp0057524-g/GB_"+fn+".txt");
  for (int i=0; i<N; i++) {
    in >> var[i];
  }
  in.close();
}

double GB::kineticEnergy() const {
  return transKE();
}

double GB::transKE() const {
  double ke = 0;
  for (int i=0; i<N; i++) {
    ke += m*vel[i].norm2sq()/2;
  }
  return ke/N; 
}

double GB::potentialEnergy() const {
  double pe = 0;
  double a = 38.39547022;
  double b = 0.05162686;
  double rmin = 2.951;
  double r0 = 2.5;
  Vec2d r, rhat;
  double rscale;
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      r = pos[j]-pos[i];
      r.periodCheck(L); 
      rscale = sqrt(r.norm2sq());
      if (rscale > 2.951) continue;
      pe += a*exp(-((rscale-r0)/b)) - a*exp(-((rmin-r0)/b));
    }
  }
  return pe;
}

void GB::calc_force(){
  double a = 38.39547022;
  double b = 0.05162686;
  double r0 = 2.5;
  Vec2d r, rhat;
  Vec2d tmpf;
  double rscale;
  for (int i=0; i<N; i++) {
    force[i] = Vec2d(0.,0.);
  }
  //pot = 0;
  
  for (int i=0; i<N; i++) {
    
    for (int j=i+1; j<N; j++) {
      
      r = pos[j]-pos[i]; 
      r.periodCheck(L);
      rscale = sqrt(r.norm2sq()); 
      if (rscale > 2.951) continue;
      rhat = r/rscale; 
      
      tmpf = ((a/b)*exp(-(rscale - r0)/b))*rhat;
      force[i] -= tmpf; force[j] += tmpf;
      //pot += ;
    }
  }
  //pot /= N;
}

void GB::VV(){
  for (int i=0; i<N; i++) {
    vel[i] += force[i]*dt*0.5/m;
    pos[i] += vel[i]*dt;
    pos[i].periodCheck(L);
  }
  calc_force(); // Hence calculate a = F/m
  for (int i=0; i<N; ++i) {
    vel[i] += force[i]*dt*0.5/m;
  }
}

void GB::temp_control(double T_target){
  double T_cur_t = 0.;
  for (int i=0; i<N; i++) {
    T_cur_t += vel[i].norm2sq(); //KE_trans=3kT/2
  }
  T_cur_t *= 0.5*m/N*2./2.; //SOS Change this from 3-->2. There are two degrees of freedom in 2D.
  for (int i=0; i<N; i++) {
    vel[i] *= sqrt(T_target/T_cur_t);
  }
}


void GB::initialisation(){
  /* sampling position from fcc */
  int n = 0;
  for (int i=0; i<M; i++) {
    for (int j=0; j<M; j++) {
      pos[n] = Vec2d(i*space-L/2, j*space-L/2);
      pos[n+1] = Vec2d((i+0.5)*space-L/2, (j+0.5)*space-L/2);
      n += 2;
    }
  }
  
 
  //sampling momenta from exp(-0.5*m*v*v/(kT)) and normalize it 
  //std::normal_distribution<double> rand_init(0.0, sqrt(T0));
  
  Vec2d avevel; //Average trans velocity
  
  for (int i=0; i<N; i++) {
    vel[i] = Vec2d(normDistStd(generator),normDistStd(generator));
    avevel += vel[i]/N;
  }
  double tmp = 0;
  for (int i=0; i<N; i++) {
    vel[i] -= avevel;
    tmp += vel[i].norm2sq();
  }
  tmp /= N; 
  for (int i=0; i<N; i++) {
    vel[i] *= sqrt(2*T0/m)/sqrt(tmp);
  }
}

void GB::g_corr() {
  loaddata();
  double dr = 0.001;
  int maxbin = 9.0 / dr;
  double G[maxbin];

  for (int i = 0; i < maxbin; ++i) {
    G[i] = 0.0;
  }
  for (int i = 0; i < 50000000; ++i) {
    if (i % 10 != 0)
      continue;
    for (int j = 0; j < N; j++) {
      for (int k = j + 1; k < N; k++) {
        Vec2d dx = pos[k] - pos[j];
        dx.periodCheck(L);
        double rjk = sqrt(dx.norm2sq());
        int bin = int(rjk / dr);
        if (bin<maxbin)
        G[bin] += 1;
      }
    }
    VV();
  }

  ofstream g("./N800Rho0035648Temp0057524-g/g.txt");
  for (int bin = 0; bin < maxbin; ++bin) {
    double rlow = bin * dr;
    double rup = (bin + 1) * dr;
    double nideal = M_PI * rho * (rup * rup - rlow * rlow);
    G[bin] = G[bin] / 5000000 / N / nideal * 2;
    g << dr * bin + dr << " " << G[bin] << endl;
  }
  g.close();
}


double GB::run(int is_init_file){
  ofstream energy1("./N800Rho0035648Temp0057524-g/energy1.txt");

  int step_count = 0;
  int index;
  if (is_init_file==0) {
    cout<<"Initialising... "<<initial_steps<<" steps..."<<endl;
    for (int i=0; i<initial_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<" "<<endl;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
    savedata();                       //cnm
//    cout<<(float) t1/CLOCKS_PER_SEC<<endl;
//    exit(0);
    cout<<"Compressing... "<<compression_steps<<" steps..."<<endl;
    double density = init_density;
    for (int i=0; i<compression_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<" "<<endl;
      }
      temp_control(T_compress);
      if ((i+1)%100==0) {
        density += density_interval;
        L = pow(N/density, 1./2.);
        rcut = L/2;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Cooling... "<<cool_steps<<" steps..."<<endl;
    double T_cur = T_compress;
    for (int i=0; i<cool_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<" "<<endl;
      }
      if ((i+1)%1000==0) {
        T_cur -= cool_interval;
        temp_control(T_cur);
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Initial Rescaling... "<<init_rescale_steps<<" steps..."<<endl;
    for (int i=0; i<init_rescale_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<" "<<endl;
      }
      if ((i+1)%1000==0) {
        temp_control(T);
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Initial equilibrating... "<<init_equil_steps<<" steps..."<<endl;
    for (int i=0; i<init_equil_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<" "<<endl;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Rescaling... "<<rescale_steps<<" steps..."<<endl;
    for (int i=0; i<rescale_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<" "<<endl;
      }
      if ((i+1)%1000==0) {
        temp_control(T);
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Equilibrating... "<<equil_steps<<" steps..."<<endl;
    
    for (int i=0; i<equil_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<" "<<endl;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
    savedata();                       //cnm
  }else {
    
    //    space = L/M;
    loaddata();
  }
  
  energy1.close();

  cout<<"Collecting data... "<<data_steps<<" steps..."<<endl;
  
  return 0;
}


  
