#include "scheme.hpp"

///////////////////////////////////////////////////////

void MUSCL::wavehllc(std::vector<triple> &wave, std::vector<U> Ui, PhysicModel &model)
{
   double rhoL,factorL,uL,vL,pL,aL;
   double rhoR,factorR,uR,vR,pR,aR;
   rhoL = Ui[0].x1;
   factorL = 1./rhoL;
   uL = Ui[0].x2 * factorL;
   vL = Ui[0].x3 * factorL;
   pL = model.EOS(Ui[0]);
   aL = sqrt(model.gamma()*pL*factorL);
   
   for(int j=1; j<Ui.size(); j++)
   {
      rhoR = Ui[j].x1;
      factorR = 1./rhoR;
      uR = Ui[j].x2 * factorR;
      vR = Ui[j].x3 * factorR;
      pR = model.EOS(Ui[j]);
      aR = sqrt(model.gamma()*pR*factorR);
      
      
      double as,rs,sL,sR,ss,el,er;
      as = 0.5*(aL+aR);
      rs = 0.5*(rhoL+rhoR);
      sL = std::min(uL-aL,uR-aR);
      sR = std::max(uL+aL, uR+aR);
      ss = 0.5*((uL+uR)+(pL-pR)/(rs*as));
      el = 0.5*rhoL*(uL*uL+vL*vL)+pL/(model.gamma()-1.);
      er = 0.5*rhoR*(uR*uR+vR*vR)+pR/(model.gamma()-1.);
      
      U UL(rhoL,rhoL*ss,rhoL*vL,el+(ss-uL)*(ss*rhoL+pL/(sL-uL)));
      U UR(rhoR,rhoR*ss,rhoR*vR,er+(ss-uR)*(ss*rhoR+pR/(sR-uR)));
      U qL = UL*((sL-uL)/(sL-ss));
      U qR = UR*((sR-uR)/(sR-ss));
      
      std::get<0>(wave[j]) = qL - Ui[j-1];
      std::get<0>(wave[j]) = qR - qL;
      std::get<0>(wave[j]) = Ui[j] - qR;
      
      rhoL = rhoR;
      uL = uR;
      vL = vR;
      pL = pR;
      aL = aR;
   }

}


void MUSCL::reconstruct(std::vector<U> &Ui, std::vector<U> &UL, std::vector<U> &UR, PhysicModel &model)
{
   for(int j=1; j<Ui.size()-1; j++)
   {
      for(int k=0; k<4; k++){
      if(Ui[j][k] != Ui[j+1][k])
      {
         double r = (Ui[j][k]-Ui[j-1][k])/(Ui[j+1][k]-Ui[j][k]);
   	 // Ospre
         //double Phi = (r*r + r)/(r*r+r+1)*1.5; 
         // Albada 2 
         //double Phi = (r*2.)/(r*r+1.);
         // Minmod
         double Phi = std::max(0.,std::min(1.,r));
         double Dp = Ui[j+1][k] - Ui[j][k];
         double Dm = Ui[j][k] - Ui[j-1][k];
         UR[j][k] = Ui[j][k] + Dp*Phi*0.5;
         UL[j][k] = Ui[j][k] - Dp*Phi*0.5;
      }
      else
      {
         UL[j][k] = Ui[j][k];
         UR[j][k] = Ui[j][k];
      }
      }
   }

}

void MUSCL::evolution_x(std::vector<U> &UL,std::vector<U> &UR, double dt, double dx, PhysicModel &model)
{

   double alpha = 0.5*dt/dx;
   for(int j=1; j<UR.size()-1; j++)
   {
      U Ualpha = (model.Fx(UR[j]) - model.Fx(UL[j]))*alpha;
      UL[j] -= Ualpha;
      UR[j] -= Ualpha;
   }
}

void MUSCL::evolution_y(std::vector<U> &UL, std::vector<U> &UR, double dt, double dy, PhysicModel &model)
{
   double alpha = 0.5*dt/dy;
   for(int j=1; j<UR.size()-1; j++)
   {
      U Ualpha = (model.Gy(UR[j]) - model.Gy(UL[j]))*alpha;
      UL[j] -= Ualpha;
      UR[j] -= Ualpha;
   }
}


//////////////////////////////////////////////////////

void Simulation::init_step_x(const U &v1, const U &v2)
{
  for(int j=0; j<_grid.Ny(); j++)
   {
     for(int i=0; i<_grid.Nx()/2; i++)
     {
       _grid.at(i,j) = v1;
     }
     _grid.at(_grid.Nx()/2,j) = (v1+v2)*0.5;
     for(int i=_grid.Nx()/2+1; i<_grid.Nx()+1; i++)
     {
       _grid.at(i,j) = v2;
     }
   }
}

void Simulation::init_step_y(const U &v1, const U &v2)
{
  for(int j=0; j<_grid.Nx(); j++)
   {
     for(int i=0; i<_grid.Ny()/2; i++)
     {
       _grid.at(j,i) = v1;
     }
     _grid.at(j,_grid.Ny()/2) = (v1+v2)*0.5;
     for(int i=_grid.Ny()/2+1; i<_grid.Ny()+1; i++)
     {
       _grid.at(j,i) = v2;
     }
   }
}

void Simulation::init_square(const U &v1, const U &v2)
{
  for(int j=0; j<_grid.Ny(); j++)
   {
     for(int i=0; i<_grid.Nx(); i++)
     {
        if((j>_grid.Ny()/3 && j<_grid.Ny()*2./3) && (i>_grid.Nx()/3 && i<_grid.Nx()*2./3))
           _grid.at(i,j) = v1;
        else if((j>=_grid.Ny()/3 && j<=_grid.Ny()*2./3+1) && (i>=_grid.Nx()/3 && i<=_grid.Nx()*2./3+1))
           _grid.at(i,j) = (v2+v1)*0.5;
        else
           _grid.at(i,j) = v2;
     }
   }
}

void Simulation::init_circle(const U &v1, const U &v2, double r)
{
  double r2 = r*r;
  for(int j=0; j<_grid.Ny(); j++)
   {
     double d1 = pow(j-_grid.Ny()/2,2);
     for(int i=0; i<_grid.Nx(); i++)
     {
        double d2 = pow(i-_grid.Nx()/2,2);
        if(d1+d2 < r2)
           _grid.at(i,j) = v1;
        else
           _grid.at(i,j) = v2;
     }
   }
}

void Simulation::init_two_circles(const U &v1, const U &v2, const U &v3, double r)
{
  double r2 = r*r;
  for(int j=0; j<_grid.Ny(); j++)
   {
     double d3 = pow(j-_grid.Ny()/2,2);
     for(int i=0; i<_grid.Nx(); i++)
     {
        double d1 = pow(i-_grid.Nx()/4,2);
        double d2 = pow(i-_grid.Nx()*3/4,2);
        if(d1+d3 < r2)
           _grid.at(i,j) = v1;
        else if(d3+d2 < r2)
           _grid.at(i,j) = v2;
        else
           _grid.at(i,j) = v3;
     }
   }
}

void Simulation::init_diagonal(const U &v1, const U &v2)
{
  for(int j=0; j<_grid.Ny(); j++)
   {
     for(int i=0; i<_grid.Nx(); i++)
     {
        if(i+j<_grid.Nx())
           _grid.at(i,j) = v1;
        else
           _grid.at(i,j) = v2;
     }
   }
}

void Simulation::init_two_flow(const U &v1, const U &v2)
{
  for(int j=0; j<_grid.Ny(); j++)
   {
     if(j<_grid.Ny()/2)
     for(int i=0; i<_grid.Nx(); i++)   
           _grid.at(i,j) = v1;
     else
     for(int i=0; i<_grid.Nx(); i++) 
           _grid.at(i,j) = v2;
   }
}

void Simulation::init_jet_x(const U &v1, const U &v2)
{
  for(int j=0; j<_grid.Ny(); j++)
  {
     for(int i=0; i<_grid.Nx(); i++)
     if(i<_grid.Nx()/4 && (j>_grid.Ny()/2-5 && j<_grid.Ny()/2+5))
        _grid.at(i,j) = v1;
     else 
        _grid.at(i,j) = v2;
   }
}

void Simulation::force_jet_x(const U &v1)
{
  for(int j=_grid.Ny()/2-5; j<_grid.Ny()/2+5; j++)
  {
     for(int i=0; i<_grid.Nx()/4; i++)
        _grid.at(i,j) = v1;
   }
}

//////////////////////////////////////////////////////////////////


void Simulation::boundaries_extend()
{
   std::vector<U> vec;
    // Rows
   vec = _grid.row(0);
   _grid.setrow(-1,vec);
   _grid.setrow(-2,vec);
   vec = _grid.row(_grid.Ny()-1);
   _grid.setrow(_grid.Ny(),vec);
   // Columns
   vec = _grid.col(0);
   _grid.setcol(-1,vec);
   _grid.setcol(-2,vec);
   vec = _grid.col(_grid.Nx()-1);
   _grid.setcol(_grid.Nx(),vec);
  
}
void Simulation::EDP2D(double dt)
{
   dimensionSplit(dt);
}


void Simulation::dimensionSplit(double dt)
{
   if(_splitting.compare("Strang") == 0)
   {
      EDP1D_x(dt/2);
      boundaries_extend();
      EDP1D_y(dt);
      boundaries_extend();
      EDP1D_x(dt/2);
   }
   else
      std::cout<<"No splitting selected!"<<std::endl;
}
  

void Simulation::EDP1D_x(double dt)
{
   int gridsize = _grid.Nx()+3;
   
#pragma omp parallel for num_threads(10)
   for(int j=0; j<_grid.Ny(); j++)
   {
      std::vector<U> Unew(gridsize);
      std::vector<U> UL(gridsize);
      std::vector<U> UR(gridsize);
      std::vector<U> Uhllc(gridsize);
      std::vector<U> Fhllc(gridsize);
   
      std::vector<U> Ui = _grid.row(j);
   
      if(_scheme.compare("MUSCL") == 0)
      {
         _muscl.reconstruct(Ui, UL, UR, _model);
         _muscl.evolution_x(UL, UR, dt, _grid.dx(), _model);
      }
      else
         std::cout<<"No scheme selected!"<<std::endl;
  
      for(int i=2; i<Ui.size()-1; i++)
      {
         if(!_solver.compare("HLLC"))
           _solveurHLLC.solve_x(UL[i], UR[i-1], _model, Uhllc[i], Fhllc[i]);
         else
           _solveurHLLE.solve_x(UL[i], UR[i-1],  _model, Fhllc[i]);
      }
   
      // Update state
      for(int i=2; i<Ui.size()-2; i++)
      {
         Unew[i] = Ui[i] - (Fhllc[i+1] - Fhllc[i])*(dt/_grid.dx());
      }
      Unew[0] = Ui[0];
      Unew[1] = Ui[1];
      Unew[Ui.size()-2] = Ui[Ui.size()-2];
      Unew[Ui.size()-1] = Ui[Ui.size()-1];
   
      _grid.setrow(j,Unew);
   
   }
}


void Simulation::EDP1D_y(double dt)
{
   int gridsize = _grid.Ny()+3;
   
#pragma omp parallel for num_threads(10)
   for(int j=0; j<_grid.Nx(); j++)
   {
      std::vector<U> Unew(gridsize);
      std::vector<U> UL(gridsize);
      std::vector<U> UR(gridsize);
      std::vector<U> Uhllc(gridsize);
      std::vector<U> Fhllc(gridsize);
   
      std::vector<U> Ui = _grid.col(j);
   
      if(_scheme.compare("MUSCL") == 0)
      {
         _muscl.reconstruct(Ui, UL, UR, _model);
         _muscl.evolution_y(UL, UR, dt, _grid.dx(), _model);
      }
      else
         std::cout<<"No scheme selected!"<<std::endl;
  
      for(int i=2; i<Ui.size()-1; i++)
      {
         if(!_solver.compare("HLLC"))
           _solveurHLLC.solve_y(UL[i], UR[i-1], _model, Uhllc[i], Fhllc[i]);
         else
           _solveurHLLE.solve_y(UL[i], UR[i-1], _model, Fhllc[i]);
      }
   
      // Update state
      for(int i=2; i<Ui.size()-2; i++)
      {
         Unew[i] = Ui[i] - (Fhllc[i+1] - Fhllc[i])*(dt/_grid.dx());
      }
      Unew[0] = Ui[0];
      Unew[1] = Ui[1];
      Unew[Ui.size()-2] = Ui[Ui.size()-2];
      Unew[Ui.size()-1] = Ui[Ui.size()-1];
   
      _grid.setcol(j,Unew);
   
   }
}


