#include "solveur.hpp"
#include "physics.hpp"





void Solveur_HLLC::solve_x(U &UR, U &UL, PhysicModel &model, U &Uhllc, U &Fhllc)
{
   double gamma = model.gamma();
   U FL,FR;
   double pL = model.EOS(UL);
   double pR = model.EOS(UR);
   
   FL = model.Fx(UL);
   FR = model.Fx(UR);
   
   double aL = sqrt(gamma*pL/UL.x1);
   double aR = sqrt(gamma*pR/UR.x1);
   
   double as = 0.5*(aL+aR);
   double rhos = 0.5*(UL.x1+UR.x1);
   double uL = UL.x2/UL.x1;
   double uR = UR.x2/UR.x1;
   double sL = std::min(uL-aL,uR-aR);
   double sR = std::max(uL+aL,uR+aR);
   double ss = 0.5*(uL+uR)+0.5*(pL-pR)/(rhos*as);

   if (sL >= 0.)
   {
      Uhllc = UL;
      Fhllc = FL;         
      return;
   }
   else if ((sL < 0.) & (ss >= 0.))
   {
      double factor_l = (sL-uL)/(sL-ss);
      Uhllc.x1 = factor_l*UL.x1;
      Uhllc.x2 = factor_l*UL.x1*ss;
      Uhllc.x3 = factor_l*UL.x3;
      Uhllc.x4 = factor_l*(UL.x4+(ss-uL)*(ss*UL.x1+pL/(sL-uL)));        

      Fhllc = FL + (Uhllc - UL)*sL;
      return;
   }
   else if ((ss < 0.) & (sR >= 0.))
   {
      double factor_r = (sR-uR)/(sR-ss);
      Uhllc.x1 = factor_r*UR.x1;
      Uhllc.x2 = factor_r*UR.x1*ss;
      Uhllc.x3 = factor_r*UR.x3;
      Uhllc.x4 = factor_r*(UR.x4+(ss-uR)*(ss*UR.x1+pL/(sR-uR)));

      Fhllc = FR + (Uhllc - UR)*sR;
      return;
   }
   else if (sR < 0.)
   {
      Uhllc = UR;
      Fhllc = FR;
      return;
   }
   else
      std::cout<<"WRONG HLLC CASE"<< sL<<" "<<sR<<" "<<ss<<std::endl;
}


void Solveur_HLLC::solve_y(U &UR, U &UL, PhysicModel &model, U &Uhllc, U &Fhllc)
{
   double gamma = model.gamma();
   U FL,FR;
   double pL = model.EOS(UL);
   double pR = model.EOS(UR);
   
   FL = model.Gy(UL);
   FR = model.Gy(UR);
   
   double aL = sqrt(gamma*pL/UL.x1);
   double aR = sqrt(gamma*pR/UR.x1);
   
   double as = 0.5*(aL+aR);
   double rhos = 0.5*(UL.x1+UR.x1);
   double uL = UL.x3/UL.x1;
   double uR = UR.x3/UR.x1;
   double sL = std::min(uL-aL,uR-aR);
   double sR = std::max(uL+aL,uR+aR);
   double ss = 0.5*(uL+uR)+0.5*(pL-pR)/(rhos*as);

   if (sL >= 0.)
   {
      Uhllc = UL;
      Fhllc = FL;         
      return;
   }
   else if ((sL < 0.) & (ss >= 0.))
   {
      double factor_l = (sL-uL)/(sL-ss);
      Uhllc.x1 = factor_l*UL.x1;
      Uhllc.x3 = factor_l*UL.x1*ss;
      Uhllc.x2 = factor_l*UL.x2;
      Uhllc.x4 = factor_l*(UL.x4+(ss-uL)*(ss*UL.x1+pL/(sL-uL)));        

      Fhllc = FL + (Uhllc - UL)*sL;
      return;
   }
   else if ((ss < 0.) & (sR >= 0.))
   {
      double factor_r = (sR-uR)/(sR-ss);
      Uhllc.x1 = factor_r*UR.x1;
      Uhllc.x3 = factor_r*UR.x1*ss;
      Uhllc.x2 = factor_r*UR.x2;
      Uhllc.x4 = factor_r*(UR.x4+(ss-uR)*(ss*UR.x1+pL/(sR-uR)));

      Fhllc = FR + (Uhllc - UR)*sR;
      return;
   }
   else if (sR < 0.)
   {
      Uhllc = UR;
      Fhllc = FR;
      return;
   }
   else{
      std::cout<<"WRONG HLLC CASE"<< sL<<" "<<sR<<" "<<ss<<std::endl;
    }
}


////////////////////////////////////////////////////////////////////////////////////



void Solveur_HLLE::solve_x(U &UL, U &UR, PhysicModel &model, U &Fhlle)
{
   double rL,uL,vL,eL,hL;
   rL = UR.x1;
   double factorL = 1./rL;
   uL = UR.x2*factorL;
   vL = UR.x3*factorL;
   eL = UR.x4;
   double pL = model.EOS(UR);
   hL = (eL+pL)*factorL;
   double sqrL = sqrt(rL);
   U FL = U(rL*uL, rL*uL*uL+pL, rL*uL*vL, uL*(eL+pL));
   
   U diff = UL - UR;
   
   if(diff == U(0,0,0,0))
   {
      Fhlle = FL;
      return;
   }
   
   double gammar = model.gamma()-1.0;
   double rR,uR,vR,eR,hR;
   rR = UL.x1;
   double factorR = 1./rR;
   uR = UL.x2*factorR;
   vR = UL.x3*factorR;
   eR = UL.x4;
   double pR = model.EOS(UL);
   hR = (eR+pR)*factorR;
   double sqrR = sqrt(rR);
   U FR = U(rR*uR, rR*uR*uR+pR, rR*uR*vR, uR*(eR+pR));
   
   double factor_s = 1./(sqrR+sqrL);
   
   double um = (sqrL*uL+sqrR*uR)*factor_s;
   double vm = (sqrL*vL+sqrR*vR)*factor_s;
   double Vm2 = um*um + vm*vm;
   double hm = (sqrL*hL+sqrR*hR)*factor_s;
   double am = sqrt(gammar*(hm-0.5*Vm2));
   
   double sL = std::min(um-am,0.0);
   double sR = std::max(um+am,0.0);
   
   U L2 = U(am*am/gammar-0.5*Vm2, um, vm, -1.0)*(gammar/(am*am));
   U L3 = U(-vm, 0.0, 1.0, 0.0);
   U R2 = U(1.0, um, vm, 0.5*Vm2);
   U R3 = U(0.0, 0.0, 1.0, vm);
   
   double eta2 = L2.dot(diff);
   double eta3 = L3.dot(diff);
   
   Fhlle = (FL*sR-FR*sL + diff*sL*sR)*(1.0/(sR-sL)) 
         - (R2*eta2+R3*eta3)*(0.5*sR*sL/(sqrt(Vm2)+am));
         
}

void Solveur_HLLE::solve_y(U &UL, U &UR, PhysicModel &model, U &Fhlle)
{
   double rL,uL,vL,eL,hL;
   rL = UR.x1;
   double factorL = 1./rL;
   uL = UR.x2*factorL;
   vL = UR.x3*factorL;
   eL = UR.x4;
   double pL = model.EOS(UR);
   hL = (eL+pL)*factorL;
   double sqrL = sqrt(rL);
   U FL = U(rL*vL, rL*vL*uL, rL*vL*vL+pL, vL*(eL+pL));
   
   U diff = UL - UR;
   
   if(diff == U(0,0,0,0))
   {
      Fhlle = FL;
      return;
   }
   
   double gammar = model.gamma()-1.0;
   double rR,uR,vR,eR,hR;
   rR = UL.x1;
   double factorR = 1./rR;
   uR = UL.x2*factorR;
   vR = UL.x3*factorR;
   eR = UL.x4;
   double pR = model.EOS(UL);
   hR = (eR+pR)*factorR;
   double sqrR = sqrt(rR);
   U FR = U(rR*vR, rR*vR*uR, rR*vR*vR+pR, vR*(eR+pR));
   
   double factor_s = 1./(sqrR+sqrL);
   
   double um = (sqrL*uL+sqrR*uR)*factor_s;
   double vm = (sqrL*vL+sqrR*vR)*factor_s;
   double Vm2 = um*um + vm*vm;
   double hm = (sqrL*hL+sqrR*hR)*factor_s;
   double am = sqrt(gammar*(hm-0.5*Vm2));
   
   double sL = std::min(vm-am,0.0);
   double sR = std::max(vm+am,0.0);
   
   U L2 = U(am*am/gammar-0.5*Vm2, um, vm, -1.0)*(gammar/(am*am));
   U L3 = U(-um, 1., 0., 0.);
   U R2 = U(1.0, um, vm, 0.5*Vm2);
   U R3 = U(0., 1., 0., um);
 
   double eta2 = L2.dot(diff);
   double eta3 = L3.dot(diff);
   
   Fhlle = (FL*sR-FR*sL + diff*sL*sR)*(1.0/(sR-sL))
         - (R2*eta2+R3*eta3)*(0.5*sR*sL/(sqrt(Vm2)+am));
         
}


