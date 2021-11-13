#ifndef SCHEME_H_
#define SCHEME_H_

#include<string>
#include <tuple>

#include <omp.h>

#include "solveur.hpp"
#include "grid.hpp"
#include "physics.hpp"

typedef std::tuple<U,U,U> triple;

class MUSCL
{
   
public:

   void wavehllc(std::vector<triple> &wave, std::vector<U> Ui, PhysicModel &model);
   
   void reconstruct(std::vector<U> &Ui, std::vector<U> &UL, std::vector<U> &UR, PhysicModel &model);
   
   void evolution_x(std::vector<U> &UL, std::vector<U> &UR, double dt, double dx, PhysicModel &model);
   
   void evolution_y(std::vector<U> &UL, std::vector<U> &UR, double dt, double dy, PhysicModel &model);
};


class Simulation
{
protected:
  Grid<U> _grid;
  Solveur_HLLC _solveurHLLC;
  Solveur_HLLE _solveurHLLE;
  MUSCL _muscl;
  PhysicModel _model;
  
  std::string _splitting;
  std::string _scheme;
  std::string _solver;
  
  double _dt;
  
public:
  Simulation(double Lx, double Ly, double Nx, double Ny, std::string s): _grid(Lx,Ly,Nx,Ny) 
  {
     _scheme = "MUSCL";
     _splitting = "Strang";
     _solver = s;
  }
  
  void uniformboundaries(U &val){_grid.init_border_to_value(val);}
  
  void boundaries_extend();
  
  void setgridval(U & val, int i, int j) {_grid.at(i,j) = val;}
  
  void init_step_x(const U &v1, const U &v2);
  
  void init_step_y(const U &v1, const U &v2);
  
  void init_square(const U &v1, const U &v2);
  
  void init_circle(const U &v1, const U &v2, double r);
  
  void init_two_circles(const U &v1, const U &v2, const U &v3, double r);
  
  void init_two_flow(const U &v1, const U &v2);
  
  void init_diagonal(const U &v1, const U &v2);
  
  void init_jet_x(const U &v1, const U &v2);
  
  void force_jet_x(const U &v1);
  
  void print(int idx){_grid.print(idx);}
  
  void EDP1D_x(double dt);
  
  void EDP1D_y(double dt);
  
  void dimensionSplit(double dt);
   
  void EDP2D(double dt);
  
  Grid<U> getGrid(){return _grid;}

};



#endif

