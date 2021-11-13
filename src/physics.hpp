#ifndef PHYSICS_H_
#define PHYSICS_H_

#include<vector>
#include<cmath>




// define the element of the grid
class U
{
public:
   double x1;
   double x2;
   double x3;
   double x4;
   
   U(){}
   
   U(double u1, double u2, double u3, double u4): x1(u1), x2(u2), x3(u3), x4(u4) {}
   
   //////////// Redefine useful operations on grid elements /////////////
   U& operator+=(const U& rhs)
   {
      this->x1 += rhs.x1;
      this->x2 += rhs.x2;
      this->x3 += rhs.x3;
      this->x4 += rhs.x4;
      return *this;
   }
   
   friend U operator+(U lhs, const U& rhs)
   {
     lhs += rhs;
     return lhs;
   }
  
   U& operator+=(const double rhs)
   {
      this->x1 += rhs;
      this->x2 += rhs;
      this->x3 += rhs;
      this->x4 += rhs;
      return *this;
   }
   
  friend U operator+(U lhs, const double rhs)
  {
    lhs += rhs;
    return lhs;
  }
   
   U& operator-=(const U& rhs)
   {
      this->x1 -= rhs.x1;
      this->x2 -= rhs.x2;
      this->x3 -= rhs.x3;
      this->x4 -= rhs.x4;
      return *this;
   }
   
  friend U operator-(U lhs, const U& rhs)
  {
    lhs -= rhs;
    return lhs;
  }
  
  U& operator*=(const double a)
   {
      this->x1 *= a;
      this->x2 *= a;
      this->x3 *= a;
      this->x4 *= a;
      return *this;
   }
  
  friend U operator*(U lhs, const double rhs)
  {
    lhs *= rhs;
    return lhs;
  }
  
  U& operator*=(const U a)
   {
      this->x1 *= a.x1;
      this->x2 *= a.x2;
      this->x3 *= a.x3;
      this->x4 *= a.x4;
      return *this;
   }
  
  friend U operator*(U lhs, const U rhs)
  {
    lhs *= rhs;
    return lhs;
  }
  
  U& operator/=(const double a)
   {
      this->x1 /= a;
      this->x2 /= a;
      this->x3 /= a;
      this->x4 /= a;
      return *this;
   }
   
   U& operator/=(const U v)
   {
      this->x1 /= v.x1;
      this->x2 /= v.x2;
      this->x3 /= v.x3;
      this->x4 /= v.x4;
      return *this;
   }
  
  friend U operator/(U lhs, const U rhs)
  {
    lhs /= rhs;
    return lhs;
  }
  
  bool operator!=(const U v)
   {
      return (this->x1 != v.x1 ||
      this->x2 != v.x2 ||
      this->x3 != v.x3 ||
      this->x4 != v.x4);
   }
   
   bool operator==(const U v)
   {
      return (this->x1 == v.x1 &
      this->x2 == v.x2 &
      this->x3 == v.x3 &
      this->x4 == v.x4);
   }
   
   U abs()
   {
      U Uabs;
      Uabs.x1 = fabs(this->x1);
      Uabs.x2 = fabs(this->x2);
      Uabs.x3 = fabs(this->x3);
      Uabs.x4 = fabs(this->x4);
      return Uabs;
   }
   
   double &operator[](int i)
   {
      if(i==0)
         return this->x1;
      if(i==1)
         return this->x2;
      if(i==2)
         return this->x3;
      else
         return this->x4;
   }
   
   double dot(const U &v)
   {
      return this->x1*v.x1 + this->x2*v.x2  + this->x3*v.x3 + this->x4*v.x4; 
   }
  
  void print(std::string s){std::cout<<s<<" is : ("<<x1<<";"<<x2<<";"<<x3<<";"<<x4<<")"<<std::endl;}
};




class PhysicModel
{
public:
   double _gamma;
   
   PhysicModel(){_gamma = 1.4;}
   
   inline double EOS(U &Uin)
   {
      double vx = Uin.x2/Uin.x1;
      double vy = Uin.x3/Uin.x1;
      return (_gamma-1.0)*(Uin.x4 - 0.5*Uin.x1*(vx*vx+vy*vy));
   }
   
   U Fx(U &Uin)
   {
      U Uout;
      double p = EOS(Uin);
      Uout.x1 = Uin.x2;
      Uout.x2 = Uin.x2*Uin.x2/Uin.x1 + p;
      Uout.x3 = Uin.x2*Uin.x3/Uin.x1;
      Uout.x4 = Uin.x2/Uin.x1*(Uin.x4 + p);
      return Uout;
   }
   
   U Gy(U &Uin)
   {
      U Uout;
      double p = EOS(Uin);
      Uout.x1 = Uin.x3;
      Uout.x2 = Uin.x2*Uin.x3/Uin.x1;
      Uout.x3 = Uin.x3*Uin.x3/Uin.x1 + p;
      Uout.x4 = Uin.x3/Uin.x1*(Uin.x4 + p);
      return Uout;
   }
   
   inline double gamma(){return _gamma;}
};

#endif

