

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
  _savedata(ort, "ort");
  _savedata(ortvel, "ortvel");
}

void GB::_savedata(const Vec2d* var, const string& fn) const{
  ofstream out("./N128Rho0285Temp1-eps/GB_"+fn+".txt");
  for (int i=0; i<N; i++) {
    out << var[i];
  }
  out.close();
}



void GB::loaddata(){
  _loaddata(pos, "pos");
  _loaddata(vel, "vel");
  _loaddata(ort, "ort");
  _loaddata(ortvel, "ortvel");
  L =  pow(N/final_density, 1.0/2.0);
  for (int i=0; i<N; i++) {
    ort[i].normalize();
    ortvelchg[i] =  (Vec2d::perp(ortforce[i], ort[i])*dt/I-2.*(ortvel[i]*ort[i])*ort[i])/2;
  }
  calc_force();
}

void GB::_loaddata(Vec2d* var, const string& fn) {
  ifstream in("./N128Rho0285Temp1-eps/GB_"+fn+".txt");
  for (int i=0; i<N; i++) {
    in >> var[i];
  }
  in.close();
}

double GB::kineticEnergy() const {
  return transKE()+rotKE();
}

double GB::transKE() const {
  double ke = 0;
  for (int i=0; i<N; i++) {
    ke += m*vel[i].norm2sq()/2;
  }
  return ke/N; //Calculating the intensive property
}

double GB::rotKE() const {
  double ke = 0;
  for (int i=0; i<N; i++) {
    ke += I*ortvel[i].norm2sq()/2;
  }
  return ke/N;
}

double GB::potentialEnergy() const {
  double pe = 0;
  double rinv;
  Vec2d r, rhat;
  double rscale, epsilon, epsilonprime, sigma, oij, roiaj, roimj;
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      r = pos[j]-pos[i];
      r.periodCheck(L); 
      rscale = sqrt(r.norm2sq());
      if (rscale>6) continue;
      rhat = r/rscale;
      oij = ort[i]*ort[j]; roiaj = rhat*(ort[i]+ort[j]); roimj = rhat*(ort[i]-ort[j]);
      epsilon = 1./sqrt(1-pow(chi*oij,2));
      epsilonprime = 1-chiprime/2.*(roiaj*roiaj/(1+chiprime*oij)+roimj*roimj/(1-chiprime*oij));
      sigma = sigma0/sqrt(1-chi/2.*(roiaj*roiaj/(1+chi*oij)+roimj*roimj/(1-chi*oij)));
      rinv = sigma0/(rscale-sigma+sigma0);// + 0.00000000001;
      pe += 4*epsilon0*epsilon*epsilonprime*epsilonprime*(pow(rinv, 12)-pow(rinv, 6));
    }
    
  }
  return pe;
}

void GB::calc_force(){
//  start = clock();
  double rinv;
  Vec2d r, rhat;
  double rscale, epsilon, epsilonprime, sigma, oij, roiaj, roimj;
  Vec2d pirr, pjrr, prr, tmpf;
  double pur, puir, pujr, r126, r7213, epsisum;
  double puij, pepij, peppij, psij;
  double rachipij, rmchipij, rachiij, rmchiij;
  for (int i=0; i<N; i++) {
    force[i] = Vec2d(0.,0.);
    ortforce[i] = Vec2d(0.,0.);
//    force[i].x=0.; force[i].y=0.; force[i].z=0.;
//    ortforce[i].x=0.; ortforce[i].y=0.; ortforce[i].z=0.;
  }
  pot = 0;
  
  for (int i=0; i<N; i++) {
    
    for (int j=i+1; j<N; j++) {
      // Force calculation
      r = pos[j]-pos[i]; //Distance between two nematogens
      r.periodCheck(L);
      rscale = sqrt(r.norm2sq()); //Modulus of r vector
      if (rscale>6) continue;
      rhat = r/rscale; //Unit vector
      oij = ort[i]*ort[j]; roiaj = rhat*ort[i]+rhat*ort[j]; roimj = rhat*ort[i]-rhat*ort[j];
      rachipij = chiprime*roiaj/(1.+chiprime*oij); rmchipij = chiprime*roimj/(1.-chiprime*oij);
      rachiij = chi*roiaj/(1.+chi*oij); rmchiij = chi*roimj/(1.-chi*oij);
      epsilon = 1./sqrt(1.-pow(chi*oij,2));
      epsilonprime = 1.-(roiaj*rachipij+roimj*rmchipij)/2.;
      
      sigma = sigma0/sqrt(1.-(roiaj*rachiij+roimj*rmchiij)/2.);
      rinv = sigma0/(rscale-sigma+sigma0);// + 0.00000000001;
      
      r126 = pow(rinv, 12) - pow(rinv, 6);
      r7213 = pow(rinv, 7) - 2.*pow(rinv, 13);
      epsisum = epsilon0*epsilon*epsilonprime*epsilonprime;
      pur = 24.*epsisum/sigma0*r7213; prr = rhat;
      pirr = 1./rscale*(ort[i]-rhat*(rhat*ort[i]));
      pjrr = 1./rscale*(ort[j]-rhat*(rhat*ort[j]));
      puir = 4.*epsisum/epsilonprime*2.*r126*(-rachipij-rmchipij)-24.*epsisum/sigma0*r7213*pow(sigma, 3)/sigma0/sigma0/2.*(rachiij+rmchiij);
      pujr = 4.*epsisum/epsilonprime*2.*r126*(-rachipij+rmchipij)-24.*epsisum/sigma0*r7213*pow(sigma, 3)/sigma0/sigma0/2.*(rachiij-rmchiij);
      tmpf = -pur*prr-puir*pirr-pujr*pjrr;
//      tmpf = 4.*epsilon0*epsilon*epsilonprime*(2.*r126/rscale*(rachipij*(ort[i]+ort[j]-rhat*roiaj)+rmchipij*(ort[i]-ort[j]-rhat*roimj))+6.*epsilonprime*r7213/sigma0*(0.5*sigma*sigma*sigma/(sigma0*sigma0)/rscale*(rachiij*(ort[i]+ort[j]-rhat*roiaj)+rmchiij*(ort[i]-ort[j]-rhat*roimj))))-pur*prr;
      force[i] -= tmpf; force[j] += tmpf;
      // Torque calculation
      pepij = pow(epsilon, 3)*chi*chi*oij;
      peppij = 1./2.*(rachipij*rachipij-rmchipij*rmchipij);
      psij = -pow(sigma, 3)/sigma0/sigma0/4.*(rachiij*rachiij-rmchiij*rmchiij);
      puij = 4.*epsilon0*r126*(epsilonprime*epsilonprime*pepij+2.*epsilon*epsilonprime*peppij)-24.*epsisum/sigma0*r7213*psij;
      ortforce[i] -= puir*rhat+puij*ort[j];
      ortforce[j] -= pujr*rhat+puij*ort[i];
      
      pot += 4.*epsisum*r126;
    }
    // Update the torques here perhaps.
  }
  pot /= N;
}


double GB::potentialEnergy2(double test) const {
  double pe2 = 0;
  double rinv;
  Vec2d r, rhat;

  Vec2d gamma;

  double gammax = 0.0;
  double gammay = 1.0;

  gamma = Vec2d(gammax, gammay);

  double rscale, epsilon, epsilonprime, sigma, oij, roiaj, roimj;

  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      r = pos[j]-pos[i];
      r.periodCheck(L); 
      rscale = sqrt(r.norm2sq());
      if (rscale>6) continue;
      rhat = r/rscale;
      oij = ort[i]*ort[j]; roiaj = rhat*(ort[i]+ort[j]); roimj = rhat*(ort[i]-ort[j]);
      epsilon = 1./sqrt(1-pow(chi*oij,2));
      epsilonprime = 1-chiprime/2.*(roiaj*roiaj/(1+chiprime*oij)+roimj*roimj/(1-chiprime*oij));
      sigma = sigma0/sqrt(1-chi/2.*(roiaj*roiaj/(1+chi*oij)+roimj*roimj/(1-chi*oij)));
      rinv = sigma0/(rscale-sigma+sigma0);// + 0.00000000001;
      pe2 += 4*epsilon0*epsilon*epsilonprime*epsilonprime*(pow(rinv, 12)-pow(rinv, 6));
    }
    pe2 -= test*(2*(ort[i]*gamma)*(ort[i]*gamma)-1);
  }
  return pe2;
}


void GB::calc_force2(double test){
//  start = clock();
  double rinv;
  Vec2d r, rhat;
  double rscale, epsilon, epsilonprime, sigma, oij, roiaj, roimj;
  Vec2d pirr, pjrr, prr, tmpf;
  double pur, puir, pujr, r126, r7213, epsisum;
  double puij, pepij, peppij, psij;
  double rachipij, rmchipij, rachiij, rmchiij;

  Vec2d gamma;

  double gammax = 0.0;
  double gammay = 1.0;

  gamma = Vec2d(gammax, gammay);

  for (int i=0; i<N; i++) {
    force[i] = Vec2d(0.,0.);
    ortforce[i] = Vec2d(0.,0.);
  }
  pot = 0;
  
  for (int i=0; i<N; i++) {
    
    for (int j=i+1; j<N; j++) {
      // Force calculation
      r = pos[j]-pos[i]; //Distance between two nematogens
      r.periodCheck(L);
      rscale = sqrt(r.norm2sq()); //Modulus of r vector
      if (rscale>6) continue;
      rhat = r/rscale; //Unit vector
      oij = ort[i]*ort[j]; roiaj = rhat*ort[i]+rhat*ort[j]; roimj = rhat*ort[i]-rhat*ort[j];
      rachipij = chiprime*roiaj/(1.+chiprime*oij); rmchipij = chiprime*roimj/(1.-chiprime*oij);
      rachiij = chi*roiaj/(1.+chi*oij); rmchiij = chi*roimj/(1.-chi*oij);
      epsilon = 1./sqrt(1.-pow(chi*oij,2));
      epsilonprime = 1.-(roiaj*rachipij+roimj*rmchipij)/2.;
      
      sigma = sigma0/sqrt(1.-(roiaj*rachiij+roimj*rmchiij)/2.);
      rinv = sigma0/(rscale-sigma+sigma0);// + 0.00000000001;
      
      r126 = pow(rinv, 12) - pow(rinv, 6);
      r7213 = pow(rinv, 7) - 2.*pow(rinv, 13);
      epsisum = epsilon0*epsilon*epsilonprime*epsilonprime;
      pur = 24.*epsisum/sigma0*r7213; prr = rhat;
      pirr = 1./rscale*(ort[i]-rhat*(rhat*ort[i]));
      pjrr = 1./rscale*(ort[j]-rhat*(rhat*ort[j]));
      puir = 4.*epsisum/epsilonprime*2.*r126*(-rachipij-rmchipij)-24.*epsisum/sigma0*r7213*pow(sigma, 3)/sigma0/sigma0/2.*(rachiij+rmchiij);
      pujr = 4.*epsisum/epsilonprime*2.*r126*(-rachipij+rmchipij)-24.*epsisum/sigma0*r7213*pow(sigma, 3)/sigma0/sigma0/2.*(rachiij-rmchiij);
      tmpf = -pur*prr-puir*pirr-pujr*pjrr;
//      tmpf = 4.*epsilon0*epsilon*epsilonprime*(2.*r126/rscale*(rachipij*(ort[i]+ort[j]-rhat*roiaj)+rmchipij*(ort[i]-ort[j]-rhat*roimj))+6.*epsilonprime*r7213/sigma0*(0.5*sigma*sigma*sigma/(sigma0*sigma0)/rscale*(rachiij*(ort[i]+ort[j]-rhat*roiaj)+rmchiij*(ort[i]-ort[j]-rhat*roimj))))-pur*prr;
      force[i] -= tmpf; force[j] += tmpf;
      // Torque calculation
      pepij = pow(epsilon, 3)*chi*chi*oij;
      peppij = 1./2.*(rachipij*rachipij-rmchipij*rmchipij);
      psij = -pow(sigma, 3)/sigma0/sigma0/4.*(rachiij*rachiij-rmchiij*rmchiij);
      puij = 4.*epsilon0*r126*(epsilonprime*epsilonprime*pepij+2.*epsilon*epsilonprime*peppij)-24.*epsisum/sigma0*r7213*psij;
      ortforce[i] -= puir*rhat+puij*ort[j];
      ortforce[j] -= pujr*rhat+puij*ort[i];
      
      pot += 4.*epsisum*r126;
    }
    ortforce[i] += 4*test*(ort[i]*gamma)*gamma;     
    pot -= test*(2*(ort[i]*gamma)*(ort[i]*gamma)-1);
  }
}

void GB::calc_order(){
  for (int i=0; i<4; i++) {
    Q[i] = 0;
  }
  for (int i=0; i<N; i++) {
    Q[0] +=  (2*ort[i].x*ort[i].x-1);
    Q[1] +=  (2*ort[i].x*ort[i].y);
    Q[2] +=  (2*ort[i].x*ort[i].y);
    Q[3] += -(2*ort[i].x*ort[i].x-1);
  }
  for (int i=0; i<4; i++) {
    Q[i] /= N;
  }
  Eigen::Matrix2d Qmat;
  Qmat << Q[0], Q[1],
  Q[2], Q[3];

  es.compute(Qmat);
  double maxeig = es.eigenvalues()[0];
  int maxind = 0;
  for (int i=0; i<2; i++) {
    if (es.eigenvalues()[i]>maxeig) {
      maxeig = es.eigenvalues()[i];
      maxind = i;
    }
  }

  eigenv = es.eigenvectors().col(0); //The eigenvector
  order = maxeig;
}


void GB::VV(){
  for (int i=0; i<N; i++) {
    vel[i] += force[i]*dt*0.5/m;
    ortvel[i] += ortvelchg[i]; 
    pos[i] += vel[i]*dt;
    ort[i] += ortvel[i]*dt;
    ort[i].normalize();
    pos[i].periodCheck(L);
  }
  calc_force(); // Hence calculate a = F/m
  for (int i=0; i<N; ++i) {
    vel[i] += force[i]*dt*0.5/m;
    ortvelchg[i] =  (Vec2d::perp(ortforce[i], ort[i])*dt/I - 2.*(ortvel[i]*ort[i])*ort[i])/2;
    ortvel[i] += ortvelchg[i];
  }
}

void GB::VV2(double test){
  for (int i=0; i<N; i++) {
    vel[i] += force[i]*dt*0.5/m;
    ortvel[i] += ortvelchg[i]; 
    pos[i] += vel[i]*dt;
    ort[i] += ortvel[i]*dt;
    ort[i].normalize();
    pos[i].periodCheck(L);
  }
  calc_force2(test);
  for (int i=0; i<N; ++i) {
    vel[i] += force[i]*dt*0.5/m;
    ortvelchg[i] =  (Vec2d::perp(ortforce[i], ort[i])*dt/I - 2.*(ortvel[i]*ort[i])*ort[i])/2;
    ortvel[i] += ortvelchg[i];
  }
}

void GB::temp_control(double T_target){
  double T_cur_t = 0., T_cur_r = 0.;
  for (int i=0; i<N; i++) {
    T_cur_t += vel[i].norm2sq(); //KE_trans=3kT/2
    T_cur_r += ortvel[i].norm2sq();
  }
  T_cur_t *= 0.5*m/N*2./2.; //SOS Change this from 3-->2. There are two degrees of freedom in 2D.
  T_cur_r *= 0.5*I/N*2.;
  for (int i=0; i<N; i++) {
    vel[i] *= sqrt(T_target/T_cur_t);
    ortvel[i] *= sqrt(T_target/T_cur_r);
  }
}


void GB::initialisation(){
  /* sampling position from fcc */
  int n = 0;
  double c = pow(1./2, 1./2); //Center of mass for orientations???
  for (int i=0; i<M; i++) {
    for (int j=0; j<M; j++) {
      pos[n] = Vec2d(i*space-L/2, j*space-L/2);
      pos[n+1] = Vec2d((i+0.5)*space-L/2, (j+0.5)*space-L/2);
      //pos[n+2] = Vec2d(i*space-L/2, (j+0.5)*space-L/2);
      //pos[n+3] = Vec2d((i+0.5)*space-L/2, j*space-L/2);
      ort[n] = Vec2d(c,c); ort[n+1] = Vec2d(c,-c);
      //ort[n+2] = Vec2d(-c,c); ort[n+3] = Vec2d(-c,-c);
      n += 2;
    }
  }
  
 
  //sampling momenta from exp(-0.5*m*v*v/(kT)) and normalize it 
  //std::normal_distribution<double> rand_init(0.0, sqrt(T0));
  
  Vec2d avevel; //Average trans velocity
  double ortscale;
  for (int i=0; i<N; i++) {
    vel[i] = Vec2d(normDistStd(generator),normDistStd(generator));
    avevel += vel[i]/N;
    ortvel[i] = Vec2d(normDistStd(generator),normDistStd(generator));
    ortscale = sqrt(ortvel[i].norm2sq());
    ortvel[i] = Vec2d::perp(ortvel[i], ort[i]);
    ortvel[i] *= ortscale;
  }
  double tmp = 0;
  double tmp2 = 0;
  for (int i=0; i<N; i++) {
    vel[i] -= avevel;
    tmp += vel[i].norm2sq();
    tmp2 += ortvel[i].norm2sq();
  }
  tmp /= N; tmp2 /= N;
  for (int i=0; i<N; i++) {
    vel[i] *= sqrt(2*T0/m)/sqrt(tmp);
    ortvel[i] *= sqrt(T0/I)/sqrt(tmp2);

  }
}

void GB::g_corr() {
  loaddata();
  double dr = 0.025;
  int maxbin = 9.0 / dr;
  double G[maxbin];

  for (int i = 0; i < maxbin; ++i) {
    G[i] = 0.0;
  }
  for (int i = 0; i < 5000000; ++i) {
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

  ofstream g("./N128Rho0285Temp1-eps/g.txt");
  for (int bin = 0; bin < maxbin; ++bin) {
    double rlow = bin * dr;
    double rup = (bin + 1) * dr;
    double nideal =
        M_PI * rho * (rup * rup - rlow * rlow);
    G[bin] = G[bin] / 500000 / N / nideal * 2;
    g << dr * bin + dr << " " << G[bin] << endl;
  }
  g.close();
}


void GB::g_orient() {
  loaddata();
  double dr = 0.025;
  int maxbin = 9.0 / dr;
  double G[maxbin];

  for (int i = 0; i < maxbin; ++i) {
    G[i] = 0.0;
  }
  for (int i = 0; i < 5000000; ++i) {
    if (i % 10 != 0)
      continue;
    for (int j = 0; j < N; j++) {
      for (int k = j + 1; k < N; k++) {
        Vec2d dx = pos[k] - pos[j];
        dx.periodCheck(L);
        double rjk = sqrt(dx.norm2sq());
        int bin = int(rjk / dr);
        if (bin<maxbin)
        G[bin] += 2*(ort[j].x*ort[k].x + ort[j].y*ort[k].y)*(ort[j].x*ort[k].x + ort[j].y*ort[k].y) - 1; 
      }
    }
    VV();
  }

  ofstream gor("./N128Rho0285Temp1-eps/gor.txt");
  for (int bin = 0; bin < maxbin; ++bin) {
    double rlow = bin * dr;
    double rup = (bin + 1) * dr;
    double nideal =
        M_PI * rho * (rup * rup - rlow * rlow);
    G[bin] = G[bin] / 500000 / N / nideal * 2;
    gor << dr * bin + dr << " " << G[bin] << endl;
  }
  gor.close();
}

void GB::eos(){
  loaddata();
  double Qeos[1500];
  double V0 = 0.0;
  double test = 0.0;
  double temp[1500];
  //double temp2[4000];
  //double potener[4000];
  
  Vec2d gamma;
  double gammax = 0.0;
  double gammay = 1.0;
  gamma = Vec2d(gammax, gammay);

  for (int m = 0; m < 1500; ++m) {
   Qeos[m] = 0.0;
   temp[m] = 0.0;
   //temp2[m] = 0.0;
   //potener[m] = 0.0;
   }

  ofstream eosfile("./N128Rho0285Temp1-eps/eosfile.txt");
  //ofstream potfile("./N128Rho0285Temp1-eps/pot.txt");

  for (int k = 0; k < 1500; ++k) {
    V0 = 0.0001;
    test = V0*k;
    
    for (int i = 0; i < 5000000; ++i) {
      
      VV2(test); 
      
      if (i % 10 != 0)
      continue;

      for (int j = 0; j < N; ++j) {
        temp[k] += 2*(ort[j]*gamma)*(ort[j]*gamma) - 1; 
       // temp2[k] += -test*(2*(ort[j]*gamma)*(ort[j]*gamma)-1);
      }
       
     // temp2[k] += potentialEnergy();

    }  
    
    Qeos[k] = temp[k] / N / 500000;
    // potener[k] = temp2[k] / 30000;

    // potfile << V0*k << " " << test << " " << potener[k] << " " << potentialEnergy() << endl;
    eosfile << V0*k << " " << test << " " << Qeos[k] << endl;
  }
  eosfile.close();
  // potfile.close();
}

double GB::f(double r, float thetaj, float thetak) const {

  double poter = 0.0;
  double rinv;
  double epsilon, epsilonprime, sigma, oij, roiaj, roimj;

  oij = cos(thetaj - thetak);
  double proj1 = cos(thetak);
  double proj2 = cos(thetaj);
  roiaj = proj1 + proj2;
  roimj = proj1 - proj2;
  epsilon = 1./sqrt(1-pow(chi*oij,2));
  epsilonprime = 1-chiprime/2.*(roiaj*roiaj/(1+chiprime*oij)+roimj*roimj/(1-chiprime*oij));
  sigma = sigma0/sqrt(1-chi/2.*(roiaj*roiaj/(1+chi*oij)+roimj*roimj/(1-chi*oij)));
  rinv = sigma0/(r-sigma+sigma0); // + 0.00000000001;
  poter = 4*epsilon0*epsilon*epsilonprime*epsilonprime*(pow(rinv, 12)-pow(rinv, 6)); 

  return poter;
}

double GB::h(double t, double r, float thetaj, float thetak) const {

  double uGB = 0.0;
  double rinv;
  double epsilon, epsilonprime, sigma, oij, roiaj, roimj;
  double chebj, chebk;

  chebj = cos(t*thetaj);
  chebk = cos(t*thetak);
  oij = cos(thetaj - thetak);
  double proj1 = cos(thetak);
  double proj2 = cos(thetaj);
  roiaj = proj1 + proj2;
  roimj = proj1 - proj2;
  epsilon = 1./sqrt(1-pow(chi*oij,2));
  epsilonprime = 1-chiprime/2.*(roiaj*roiaj/(1+chiprime*oij)+roimj*roimj/(1-chiprime*oij));
  sigma = sigma0/sqrt(1-chi/2.*(roiaj*roiaj/(1+chi*oij)+roimj*roimj/(1-chi*oij)));
  rinv = sigma0/(r-sigma+sigma0); // + 0.00000000001;
  uGB = 4*epsilon0*epsilon*epsilonprime*epsilonprime*(pow(rinv, 12)-pow(rinv, 6)); 

  return chebj*chebk*uGB;
}

void GB::integration() {
  
  const double PI = 3.141592653589793238463;
  double n = 20;
  double h = PI/n;
  double C2 = (h*h)/(PI*PI*4);
  double trap = 0.0;
  double c22 = 0.0;
  double c44 = 0.0;
  double c66 = 0.0;
  
  ofstream fcnn("./N128Rho0285Temp1-eps/fcnn.txt");

  for (int a = 1; a < 4001; ++a) { 

    double r = a*0.001;

    for (int k = 0; k < n; k++) {
      for (int j = 0; j < n; j++) {          
        //trap += f(r, j*h, k*h) + f(r, (j+1)*h, k*h) + f(r, j*h, (k+1)*h) + f(r, (j+1)*h, (k+1)*h);
        c22 += h(2, r, j*h, k*h) + h(2, r, (j+1)*h, k*h) + h(2, r, j*h, (k+1)*h) + h(2, r, (j+1)*h, (k+1)*h);
        c44 += h(4, r, j*h, k*h) + h(4, r, (j+1)*h, k*h) + h(4, r, j*h, (k+1)*h) + h(4, r, (j+1)*h, (k+1)*h);
        c66 += h(6, r, j*h, k*h) + h(6, r, (j+1)*h, k*h) + h(6, r, j*h, (k+1)*h) + h(6, r, (j+1)*h, (k+1)*h);
      } 
    } 
    fcnn<<r<<"  "<<C2*trap<<" "<<C2*c22<<" "<<C2*c44<<" "<<C2*c66 <<" "<<endl;
    //trap = 0.0;
    c22 = 0.0;
    c44 = 0.0;
    c66 = 0.0;
  } 
  fcnn.close();
}


double GB::run(int is_init_file){
  ofstream energy1("./N128Rho0285Temp1-eps/energy1.txt");

  int step_count = 0;
  int index;
  if (is_init_file==0) {
    cout<<"Initialising... "<<initial_steps<<" steps..."<<endl;
    for (int i=0; i<initial_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<" "<<endl;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
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
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<" "<<endl;
      }
      temp_control(T_compress);
      if ((i+1)%100==0) {
        density += density_interval;
        L = pow(N/density, 1./2.);
        rcut = L/2;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Cooling... "<<cool_steps<<" steps..."<<endl;
    double T_cur = T_compress;
    for (int i=0; i<cool_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<" "<<endl;
      }
      if ((i+1)%1000==0) {
        T_cur -= cool_interval;
        temp_control(T_cur);
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Initial Rescaling... "<<init_rescale_steps<<" steps..."<<endl;
    for (int i=0; i<init_rescale_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<" "<<endl;
      }
      if ((i+1)%1000==0) {
        temp_control(T);
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Initial equilibrating... "<<init_equil_steps<<" steps..."<<endl;
    for (int i=0; i<init_equil_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<" "<<endl;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Rescaling... "<<rescale_steps<<" steps..."<<endl;
    for (int i=0; i<rescale_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<" "<<endl;
      }
      if ((i+1)%1000==0) {
        temp_control(T);
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
//    savedata();                       //cnm
    cout<<"Equilibrating... "<<equil_steps<<" steps..."<<endl;
    
    for (int i=0; i<equil_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        energy1<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<" "<<endl;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
    savedata();                       //cnm
  }else {
    
    //    space = L/M;
    loaddata();
  }
  
  energy1.close();

  cout<<"Collecting data... "<<data_steps<<" steps..."<<endl;      
}


  
