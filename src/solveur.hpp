#ifndef SOLVER_H_
#define SOLVER_H_

#include <vector>
#include <iostream>
#include <cmath>

#include "physics.hpp"

////////////////////////////
class Solveur_HLLC
{
public:

  void solve_x(U &UL, U &UR, PhysicModel &model, U &Uhllc, U &Fhllc);
  
  void solve_y(U &UL, U &UR, PhysicModel &model, U &Uhllc, U &Fhllc);

};

///////////////////////////:
class Solveur_HLLE
{
public:

   void solve_x(U &UL, U &UR, PhysicModel &model, U &Fhlle);
  
   void solve_y(U &UL, U &UR, PhysicModel &model, U &Fhlle);
};

#endif
