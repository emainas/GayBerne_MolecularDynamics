
#ifndef basics_h
#define basics_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
constexpr static int N_const = 1372;
constexpr static int samples_const = 10;
class Vec3d;
class Int3d {
public:
  int x,y,z;
public:
  Int3d(int X=0, int Y=0, int Z=0):x(X),y(Y),z(Z) {}
  Int3d& operator=(const Int3d& src) {
    x=src.x; y=src.y; z=src.z;
    return *this;
  }
  ~Int3d(){};
  Int3d operator+(const Int3d& src) const{
    Int3d tmp = Int3d(x+src.x,y+src.y,z+src.z);
    return tmp;
  }
  Int3d operator-(const Int3d& src) const{
    Int3d tmp = Int3d(x-src.x,y-src.y,z-src.z);
    return tmp;
  }
  //
  Int3d& operator+=(const Int3d& src) {
    x += src.x; y += src.y; z +=src.z;
    return *this;
  }
  Int3d& operator-=(const Int3d& src) {
    x -= src.x; y -= src.y; z -=src.z;
    return *this;
  }
};


class Vec3d {
public:
  double x,y,z;
public:
  Vec3d(double X, double Y, double Z):x(X),y(Y),z(Z) {}
  Vec3d(){}
  ~Vec3d(){}
  double dist_from_2(const Vec3d& src) const {
    return (x-src.x)*(x-src.x)+(y-src.y)*(y-src.y)+(z-src.z)*(z-src.z);
  }
  double norm1() const{
    return abs(x)+abs(y)+abs(z);
  }
  double norm2sq() const{
    return x*x+y*y+z*z;
  }
  
  
  //
  Vec3d& operator=(const Vec3d& src) {
    x=src.x; y=src.y; z=src.z;
    return *this;
  }
  Vec3d operator+(const Vec3d& src) const {
    Vec3d tmp = Vec3d(x+src.x,y+src.y,z+src.z);
    return tmp;
  }
  Vec3d operator-(const Vec3d& src) const {
    Vec3d tmp = Vec3d(x-src.x,y-src.y,z-src.z);
    return tmp;
  }
  //
  Vec3d& operator+=(const Vec3d& src) {
    x += src.x; y += src.y; z +=src.z;
    return *this;
  }
  Vec3d& operator-=(const Vec3d& src) {
    x -= src.x; y -= src.y; z -=src.z;
    return *this;
  }
  //
  Vec3d operator+(double d) const {
    Vec3d tmp =  Vec3d(x+d, y+d, z+d);
    return tmp;
  }
  Vec3d operator-(double d) const {
    Vec3d tmp = Vec3d(x-d, y-d, z-d);
    return tmp;
  }
  Vec3d operator*(double d) const {
    Vec3d tmp = Vec3d(x*d,y*d,z*d);
    return tmp;
  }
  Vec3d operator/(double d) const {
    Vec3d tmp = Vec3d(x/d,y/d,z/d);
    return tmp;
  }
  Vec3d operator-() const {
    Vec3d tmp = Vec3d(-x,-y,-z);
    return tmp;
  }
  
  Vec3d& operator+=(double d) {
    x += d; y += d; z += d;
    return *this;
  }
  Vec3d& operator-=(double d) {
    x -= d; y -= d; z -= d;
    return *this;
  }
  Vec3d& operator*=(double d) {
    x *= d; y *= d; z *= d;
    return *this;
  }
  Vec3d& operator/=(double d) {
    x /= d; y /= d; z /= d;
    return *this;
  }
  
  friend Vec3d operator+(double d, const Vec3d& src){
    return src+d;
  }
  friend Vec3d operator-(double d, const Vec3d& src){
    return -src+d;
  }
  friend Vec3d operator*(double d, const Vec3d& src){
    return src*d;
  }
  friend double operator*(const Vec3d& a, const Vec3d& b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
  }
  friend std::ostream& operator<<(std::ostream& os, const Vec3d& src)
  {
    os << src.x << ' ' <<src.y <<' ' <<src.z <<'\n';
    return os;
  }
  friend std::istream& operator>>(std::istream& is, Vec3d& src)
  {
    is >> src.x >> src.y >> src.z;
    return is;
  }
  static Vec3d crosspdt(const Vec3d& a, const Vec3d& b){
    Vec3d tmp;
    tmp.x = a.y*b.z - a.z*b.y;
    tmp.y = a.z*b.x - a.x*b.z;
    tmp.z = a.x*b.y - a.y*b.x;
    return tmp;
  }
  static Vec3d perp(const Vec3d& a, const Vec3d& b){
    Vec3d tmp = a-b*(a*b);
    return tmp;
  }
  Int3d periodCheck(double L){
    Int3d tmp(0,0,0);
    tmp.x = rint(x/L); x -= L*rint(x/L);
    tmp.y = rint(y/L); y -= L*rint(y/L);
    tmp.z = rint(z/L); z -= L*rint(z/L);
    return tmp;
  }
  void normalize(){
    double n2 = sqrt(norm2sq());
    *this /= n2;
  }
  double dot(const Vec3d& src) const{
    return this->x*src.x+this->y*src.y+this->z*src.z;
  }
};




class config_Space {
public:
  constexpr static int N = N_const;
  Vec3d pos[N];
public:
  config_Space(){};
  config_Space(Vec3d position[N]){
    for (int i=0; i<N; i++) {
      pos[i] = position[i];
    }
  }
  config_Space(int flag){
    if (flag==0) {
      for (int i=0; i<N; i++) {
        pos[i] = Vec3d(0.0,0.0,0.0);
      }
    }else if (flag==1) {
      for (int i=0; i<N; i++) {
        pos[i] = Vec3d(1.0,1.0,1.0);
      }
    }
    
  }
  
  
  double dist_from_2(const config_Space& src) const{
    double dist2 = 0;
    for (int i=0; i<N; i++) {
      dist2 += pos[i].dist_from_2(src.pos[i]);
    }
    return dist2;
  }
  const double norm1() const{
    double n1 = 0;
    for (int i=0; i<N; i++) {
      n1 += pos[i].norm1();
    }
    return n1;
  }
  double norm2sq() const{
    double n2 = 0;
    for (int i=0; i<N; i++) {
      n2 += pos[i].norm2sq();
    }
    return n2;
  }
  //
  void assign_var(Vec3d* var){
    for (int i=0; i<N; i++) {
      var[i] = pos[i];
    }
  }
  void assign_spc(const Vec3d* var){
    for (int i=0; i<N; i++) {
      pos[i] = var[i];
    }
  }

  config_Space& operator=(const config_Space& src) {
    for (int i=0; i<N; i++) {
      pos[i] = src.pos[i];
    }
    return *this;
  }
  config_Space operator+(const config_Space& src) const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] += src.pos[i];
    }
    return tmp;
  }
  config_Space operator-(const config_Space& src) const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] -= src.pos[i];
    }
    return tmp;
  }
  config_Space operator+(const Vec3d& src) const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] += src;
    }
    return tmp;
  }
  config_Space operator-(const Vec3d& src) const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] -= src;
    }
    return tmp;
  }
  //
  config_Space& operator+=(const config_Space& src) {
    for (int i=0; i<N; i++) {
      pos[i] += src.pos[i];
    }
    return *this;
  }
  config_Space& operator-=(const config_Space& src) {
    for (int i=0; i<N; i++) {
      pos[i] -= src.pos[i];
    }
    return *this;
  }
  config_Space& operator+=(const Vec3d& src) {
    for (int i=0; i<N; i++) {
      pos[i] += src;
    }
    return *this;
  }
  config_Space& operator-=(const Vec3d& src) {
    for (int i=0; i<N; i++) {
      pos[i] -= src;
    }
    return *this;
  }
  //
  config_Space operator+(double d) const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] += d;
    }
    return tmp;
  }
  config_Space operator-(double d) const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] -= d;
    }
    return tmp;
  }
  config_Space operator*(double d) const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] *= d;
    }
    return tmp;
  }
  config_Space operator/(double d) const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] /= d;
    }
    return tmp;
  }
  config_Space operator-() const {
    config_Space tmp = *this;
    for (int i=0; i<N; i++) {
      tmp.pos[i] += -tmp.pos[i];
    }
    return tmp;
  }
  
  
  config_Space& operator+=(double d) {
    for (int i=0; i<N; i++) {
      pos[i] += d;
    }
    return *this;
  }
  config_Space& operator-=(double d) {
    for (int i=0; i<N; i++) {
      pos[i] -= d;
    }
    return *this;
  }
  config_Space& operator*=(double d) {
    for (int i=0; i<N; i++) {
      pos[i] *= d;
    }
    return *this;
  }
  config_Space& operator/=(double d) {
    for (int i=0; i<N; i++) {
      pos[i] /= d;
    }
    return *this;
  }
  friend double operator*(const config_Space& a, const config_Space& b){
    double tmp = 0;
    for (int i=0; i<N; i++) {
      tmp += a.pos[i]*b.pos[i];
    }
    return tmp;
  }
  
  friend config_Space operator+(double d, const config_Space& src){
    return src+d;
  }
  friend config_Space operator-(double d, const config_Space& src){
    return -src+d;
  }
  friend config_Space operator*(double d, const config_Space& src){
    return src*d;
  }
  friend std::ostream& operator<<(std::ostream& os, const config_Space& src)
  {
    for (int i=0; i<N; i++) {
      os << src.pos[i];
    }
    return os;
  }
  friend std::istream& operator>>(std::istream& is, config_Space& src)
  {
    for (int i=0; i<N; i++) {
      is >> src.pos[i];
    }
    return is;
  }
  //  void periodCheck(double L){
  //    for (int i=0; i<N; i++) {
  //      pos[i].periodCheck(L);
  //    }
  //  }
  void normalize(){
    double n2 = sqrt(norm2sq());
    *this /= n2;
  }
  
};

struct config_GB {
  config_Space pos;
  config_Space ort;
};

struct config_Pair {
  config_Space r1;
  config_Space r2;
};
#endif /* basics_h */
