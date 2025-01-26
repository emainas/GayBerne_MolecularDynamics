
#include "MD_GB.hpp"
#include <iostream>
#include <math.h>
#include <random>
#include <fstream>
#include <time.h>
#include <iomanip>

using namespace std;

void GB::outputit() const{
  cout<<"Density: "<<rho<<"\nTemperature:"<<T<<endl;
  cout<<"Initial box length: "<<L<<"\nFinal box length: " << pow(N/final_density, 1.0/3.0)<<endl;
}

void GB::savedata() const{
  _savedata(pos, "pos");
  _savedata(vol, "vol");
  _savedata(ort, "ort");
  _savedata(ortvol, "ortvol");
}

void GB::_savedata(const Vec3d* var, const string& fn) const{
  ofstream out("./N1372s4/GB_"+fn+".txt");
  for (int i=0; i<N; i++) {
    out << var[i];
  }
  out.close();
}

void GB::loaddata(){
  _loaddata(pos, "pos");
  _loaddata(vol, "vol");
  _loaddata(ort, "ort");
  _loaddata(ortvol, "ortvol");
  L =  pow(N/final_density, 1.0/3.0);
  for (int i=0; i<N; i++) {
    periods[i] = 0;
    ort[i].normalize();
    ortvolchg[i] =  (Vec3d::perp(ortforce[i], ort[i])*dt/I-2.*(ortvol[i]*ort[i])*ort[i])/2;
  }
  calc_force();
}

void GB::_loaddata(Vec3d* var, const string& fn) {
  ifstream in("./N1372s4/GB_"+fn+".txt");
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
    ke += m*vol[i].norm2sq()/2;
  }
  return ke/N;
}

double GB::rotKE() const {
  double ke = 0;
  for (int i=0; i<N; i++) {
    ke += I*ortvol[i].norm2sq()/2;
  }
  return ke/N;
}

double GB::potentialEnergy() const {
  double pe = 0;
  double rinv;
  Vec3d r, rhat;
  double rscale, epsilon, epsilonprime, sigma, oij, roiaj, roimj;
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      r = pos[j]-pos[i];
      r.periodCheck(L);
      rscale = sqrt(r.norm2sq());
      rhat = r/rscale;
      oij = ort[i]*ort[j]; roiaj = rhat*(ort[i]+ort[j]); roimj = rhat*(ort[i]-ort[j]);
      epsilon = 1./sqrt(1-pow(chi*oij,2));
      epsilonprime = 1-chiprime/2.*(roiaj*roiaj/(1+chiprime*oij)+roimj*roimj/(1-chiprime*oij));
      sigma = sigma0/sqrt(1-chi/2.*(roiaj*roiaj/(1+chi*oij)+roimj*roimj/(1-chi*oij)));
      rinv = sigma0/(rscale-sigma+sigma0);// + 0.00000000001;
      pe += 4*epsilon0*epsilon*epsilonprime*epsilonprime*(pow(rinv, 12)-pow(rinv, 6));
    }
  }
  return pe/N;
//  return pot;
}


void GB::calc_force(){
//  start = clock();
  double rinv;
  Vec3d r, rhat;
  double rscale, epsilon, epsilonprime, sigma, oij, roiaj, roimj;
  Vec3d pirr, pjrr, prr, tmpf;
  double pur, puir, pujr, r126, r7213, epsisum;
  double puij, pepij, peppij, psij;
  double rachipij, rmchipij, rachiij, rmchiij;
  for (int i=0; i<N; i++) {
    force[i] = Vec3d(0.,0.,0.);
    ortforce[i] = Vec3d(0.,0.,0.);
//    force[i].x=0.; force[i].y=0.; force[i].z=0.;
//    ortforce[i].x=0.; ortforce[i].y=0.; ortforce[i].z=0.;
  }
  pot = 0;
  
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      // Force calculation
      r = pos[j]-pos[i];
      r.periodCheck(L);
      rscale = sqrt(r.norm2sq());
      rhat = r/rscale;
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
  }
  pot /= N;
}

void GB::calc_order(){
  for (int i=0; i<6; i++) {
    Q[i] = 0;
  }
  for (int i=0; i<N; i++) {
    Q[0] += 0.5*(3*ort[i].x*ort[i].x-1);
    Q[1] += 0.5*(3*ort[i].x*ort[i].y);
    Q[2] += 0.5*(3*ort[i].x*ort[i].z);
    Q[3] += 0.5*(3*ort[i].y*ort[i].y-1);
    Q[4] += 0.5*(3*ort[i].y*ort[i].z);
    Q[5] += 0.5*(3*ort[i].z*ort[i].z-1);
  }
  for (int i=0; i<6; i++) {
    Q[i] /= N;
  }
  Eigen::Matrix3d Qmat;
  Qmat << Q[0], Q[1], Q[2],
  Q[1], Q[3], Q[4],
  Q[2], Q[4], Q[5];
  es.compute(Qmat);
  double maxeig = es.eigenvalues()[0];
  int maxind = 0;
  for (int i=1; i<3; i++) {
    if (es.eigenvalues()[i]>maxeig) {
      maxeig = es.eigenvalues()[i];
      maxind = i;
    }
  }
//  cout<<order<<endl;
  eigenv = es.eigenvectors().col(0);
  order = maxeig;
}

void GB::VV(){
  for (int i=0; i<N; i++) {
    vol[i] += force[i]*dt*0.5/m;
    ortvol[i] += ortvolchg[i];// wrong (Vec3d::perp(ortforce[i], ort[i])*dt/I-2.*(ortvol[i]*ort[i])*ort[i])/2;
    pos[i] += vol[i]*dt;
    ort[i] += ortvol[i]*dt;
    ort[i].normalize();
    periods[i] += pos[i].periodCheck(L);
  }
  calc_force();
  for (int i=0; i<N; ++i) {
    vol[i] += force[i]*dt*0.5/m;
    ortvolchg[i] =  (Vec3d::perp(ortforce[i], ort[i])*dt/I - 2.*(ortvol[i]*ort[i])*ort[i])/2;
    ortvol[i] += ortvolchg[i];
  }
}

void GB::temp_control(double T_target){
  double T_cur_t = 0., T_cur_r = 0.;
  for (int i=0; i<N; i++) {
    T_cur_t += vol[i].norm2sq(); //KE_trans=3kT/2
    T_cur_r += ortvol[i].norm2sq();
  }
  T_cur_t *= 0.5*m/N*2./3.;
  T_cur_r *= 0.5*I/N;
  for (int i=0; i<N; i++) {
    vol[i] *= sqrt(T_target/T_cur_t);
    ortvol[i] *= sqrt(T_target/T_cur_r);
  }
}


void GB::initialisation(){
  /* sampling position from fcc */
  int n = 0;
  double c = pow(1./3, 1./2);
  for (int i=0; i<M; i++) {
    for (int j=0; j<M; j++) {
      for (int k=0; k<M; k++) {
        pos[n] = Vec3d(i*space-L/2, j*space-L/2, k*space-L/2);
        pos[n+1] = Vec3d((i+0.5)*space-L/2, (j+0.5)*space-L/2, k*space-L/2);
        pos[n+2] = Vec3d(i*space-L/2, (j+0.5)*space-L/2, (k+0.5)*space-L/2);
        pos[n+3] = Vec3d((i+0.5)*space-L/2, j*space-L/2, (k+0.5)*space-L/2);
        ort[n] = Vec3d(c,c,c); ort[n+1] = Vec3d(c,-c,-c);
        ort[n+2] = Vec3d(-c,c,-c); ort[n+3] = Vec3d(-c,-c,c);
        n += 4;
      }
    }
  }
  /* sampling momentum from exp(-0.5*m*v*v/(kT)) and normalize it */
//  std::normal_distribution<double> rand_init(0.0, sqrt(T0));
  Vec3d avevol(0, 0, 0);
  double ortscale;
  for (int i=0; i<N; i++) {
    vol[i] = Vec3d(normDistStd(generator),normDistStd(generator),normDistStd(generator));
    avevol += vol[i]/N;
    ortvol[i] = Vec3d(normDistStd(generator),normDistStd(generator),normDistStd(generator));
    ortscale = sqrt(ortvol[i].norm2sq());
    ortvol[i] = Vec3d::perp(ortvol[i], ort[i]);
    ortvol[i] *= ortscale;
  }

  double tmp = 0;
  double tmp2 = 0;
  for (int i=0; i<N; i++) {
    vol[i] -= avevol;
    tmp += vol[i].norm2sq();
    tmp2 += ortvol[i].norm2sq();
  }
  tmp /= N; tmp2 /= N;
  for (int i=0; i<N; i++) {
    vol[i] *= sqrt(2*T0/m)/sqrt(tmp);
    ortvol[i] *= sqrt(T0/I)/sqrt(tmp2);

  }
}

void GB::susceptibility(){

  double susc = 0.0;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
        susc += 0.5*( 3*(ort[i]*ort[j])*(ort[i]*ort[j]) - 1);
    }
  }

  chi1 = susc/N;
}

void GB::g_corr(){
  loaddata();
  double dr=0.025;
  int maxbin=9.0/dr;
  double G[maxbin];
  
  for (int i=0; i<maxbin; ++i) {
    G[i] = 0.0;
  }
  for (int i=0; i<100000; ++i) {
    if (i%10 != 0) continue;
    for (int j=0; j<N; j++) {
      for (int k=j+1; k<N; k++) {
        Vec3d dx = pos[k] - pos[j];
        dx.periodCheck(L);
        double rjk = sqrt(dx.norm2sq());
        int bin=int(rjk/dr);
        G[bin] += 1;
      }
    }
    VV();
  }
  
  ofstream g("./N1372s4/g.txt");
  for (int bin=0; bin<maxbin; ++bin) {
    double rlow=bin*dr;
    double rup=(bin+1)*dr;
    double nideal=4*M_PI*rho*(rup*rup*rup-rlow*rlow*rlow)/3.0;
    G[bin]=G[bin]/10000/N/nideal*2;
    g << dr*bin+dr<<" "<<G[bin]<<endl;
  }
  g.close();
}


double GB::run(int is_init_file){
  ofstream energy("./N1372s4/energy.txt");

  int step_count = 0;
  int index;
  if (is_init_file==0) {
    cout<<"Initialising... "<<initial_steps<<" steps..."<<endl;
    for (int i=0; i<initial_steps; i++) {
      VV();
      step_count++;
      if ((step_count+1)%10==0) {
        index = round((step_count+1)/10)-1;
        energy<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
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
        energy<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
      }
      temp_control(T_compress);
      if ((i+1)%100==0) {
        density += density_interval;
        L = pow(N/density, 1./3.);
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
        energy<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
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
        energy<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
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
        energy<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
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
        energy<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
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
        energy<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
      }
    }
    cout<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
    cout<< pos[0]<<endl;              //cnm
    savedata();                       //cnm
  }else {
    
    //    space = L/M;
    loaddata();
  }
  cout<<"Collecting data... "<<data_steps<<" steps..."<<endl;
//  ofstream energy("./data/energy.txt");
  int index2;
  ofstream forder("./N1372s4/order.txt");
  ofstream fchi("./N1372s4/chi1.txt");
  for (int i=0; i<data_steps; i++) {
    VV();
    step_count++;
    if ((step_count+1)%10==0) {
      index = round((step_count+1)/10)-1;
      index2 = round((i+1)/10)-1;
      PE[index2] = potentialEnergy();
      KE[index2] = kineticEnergy();
      calc_order();
      S[index2] = order;
      susceptibility();
      C[index2] = chi1;
      energy<<index<<" "<< potentialEnergy()<<" "<<kineticEnergy()<<" "<<rotKE()<<" "<<transKE()<<endl;
      forder<<index2<<" "<< order <<endl;
      fchi<<index2<<" "<< chi1 <<endl;
    }
  }
  energy.close();
  double avePE = 0; double aveKE = 0;
  double aveTE = 0; double aveS = 0;
  double aveChi = 0;
  int K = 1000000;
  for (int i=0; i<K; i++) {
    avePE += PE[i];
    aveKE += KE[i];
    aveTE += KE[i] + PE[i];
    aveS += S[i];
    aveChi += C[i];
  }
  avePE /= K; aveKE /= K; aveTE /= K; aveS /= K; aveChi /= K;
  cout<<"Average PE, KE, TE, S are " << avePE <<" "<< aveKE<< " "<<aveTE<<" "<<aveS<<" "<<aveChi<<endl;
  return 0;
}




