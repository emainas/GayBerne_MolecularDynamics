


#ifndef basics_h
#define basics_h

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <random> 

constexpr static int N_const = 800;

class Vec2d {
public:
  double x,y;

public:
  Vec2d(double X=0.0, double Y=0.0):x(X),y(Y) {} 
  ~Vec2d(){}
  double dist_from_2(const Vec2d& src) const {
    return (x-src.x)*(x-src.x)+(y-src.y)*(y-src.y);
  }
  double norm1() const{
    return abs(x)+abs(y);
  }
  double norm2sq() const{
    return x*x+y*y;
  }
  
  
  //
  Vec2d& operator=(const Vec2d& src) {
    x=src.x; y=src.y;
    return *this;
  }
  Vec2d operator+(const Vec2d& src) const {
    Vec2d tmp = Vec2d(x+src.x,y+src.y);
    return tmp;
  }
  Vec2d operator-(const Vec2d& src) const {
    Vec2d tmp = Vec2d(x-src.x,y-src.y);
    return tmp;
  }
  //
  Vec2d& operator+=(const Vec2d& src) {
    x += src.x; y += src.y;
    return *this;
  }
  Vec2d& operator-=(const Vec2d& src) {
    x -= src.x; y -= src.y;
    return *this;
  }
  //
  Vec2d operator+(double d) const {
    Vec2d tmp =  Vec2d(x+d, y+d);
    return tmp;
  }
  Vec2d operator-(double d) const {
    Vec2d tmp = Vec2d(x-d, y-d);
    return tmp;
  }
  Vec2d operator*(double d) const {
    Vec2d tmp = Vec2d(x*d,y*d);
    return tmp;
  }
  Vec2d operator/(double d) const {
    Vec2d tmp = Vec2d(x/d,y/d);
    return tmp;
  }
  Vec2d operator-() const {
    Vec2d tmp = Vec2d(-x,-y);
    return tmp;
  }
  
  Vec2d& operator+=(double d) {
    x += d; y += d;
    return *this;
  }
  Vec2d& operator-=(double d) {
    x -= d; y -= d; 
    return *this;
  }
  Vec2d& operator*=(double d) {
    x *= d; y *= d; 
    return *this;
  }
  Vec2d& operator/=(double d) {
    x /= d; y /= d; 
    return *this;
  }
  
  friend Vec2d operator+(double d, const Vec2d& src){
    return src+d;
  }
  friend Vec2d operator-(double d, const Vec2d& src){
    return -src+d;
  }
  friend Vec2d operator*(double d, const Vec2d& src){
    return src*d;
  }
  friend double operator*(const Vec2d& a, const Vec2d& b){
    return a.x*b.x+a.y*b.y;
  }
  friend std::ostream& operator<<(std::ostream& os, const Vec2d& src)
  {
    os << src.x << ' ' <<src.y <<'\n';
    return os;
  }
  friend std::istream& operator>>(std::istream& is, Vec2d& src)
  {
    is >> src.x >> src.y;
    return is;
  }

  static Vec2d perp(const Vec2d& a, const Vec2d& b){
    Vec2d tmp = a-b*(a*b);
    return tmp;
  }
  void periodCheck(double L){
    x -= L*rint(x/L);
    y -= L*rint(y/L);
  }
  void normalize(){
    double n2 = sqrt(norm2sq());
    *this /= n2;
  }
  double dot(const Vec2d& src) const{
    return this->x*src.x+this->y*src.y;
  }
};


#endif /* basics_h */
