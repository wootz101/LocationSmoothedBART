/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017-2018 Robert McCulloch, Rodney Sparapani
 *                          and Charles Spanbauer
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#ifndef GUARD_bart_h
#define GUARD_bart_h

#include <ctime>

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"

class bart {
public:
   //------------------------------
   //friends
   friend bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
		  std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen);
   //------------------------------
   //constructor/destructor
   bart();
   bart(size_t m);
   bart(const bart&);
   ~bart();
   //------------------------------
   //operators
   bart& operator=(const bart&);
   //------------------------------
   //get,set
   size_t getm() {return m;}
   void setm(size_t m);
   void setdata(size_t p, size_t n, double *x, double *y, size_t nc=100);
   void setdata(size_t p, size_t n, double *x, double *y, int* nc);


  // NEW SET DATA FOR FBART
   void setdata_fbart(size_t p, size_t n, size_t funcP, size_t intvls, double *x, double *y, size_t nc=100);
   void setdata_fbart(size_t p, size_t n, size_t funcP, size_t intvls, double *x, double *y, int* nc);

   void setpi(pinfo& pi) {this->pi = pi;}
   void setprior(double alpha, double beta, double tau)
      {pi.alpha=alpha; pi.mybeta = beta; pi.tau=tau;}

   void setdart(double _a, double _b, double _rho, bool _aug, bool _dart,
		double _theta=0., double _omega=1.) {
     this->a=_a; this->b=_b; this->rho=_rho; this->aug=_aug;
     this->dart=_dart; this->omega=_omega;
     if(_theta==0.){
       this->const_theta=false;
       this->theta=1.;
     }
     else{
       this->const_theta=true;
       this->theta=_theta;
     }
}

   void startdart() {this->dartOn=!(this->dartOn);}
   // add DART for functional predictors and function points

   void setdart_fp(double _a, double _b, double _rho, bool _aug, bool _dart,
                   double _theta=0., double _omega=1.) {
     this->a_fp=_a; this->b_fp=_b; this->rho_fp=_rho; this->aug_fp=_aug;
     this->dart_fp=_dart; this->omega_fp=_omega;
     if(_theta==0.){
       this->const_theta_fp=false;
       this->theta_fp=1.;
     }
     else{
       this->const_theta_fp=true;
       this->theta_fp=_theta;
     }
   }

   void setdart_int(double _a, double _b, double _rho, bool _aug, bool _dart,
                    double _theta=0., double _omega=1.) {
     this->a_int=_a; this->b_int=_b; this->rho_int=_rho; this->aug_int=_aug;
     this->dart_int=_dart; this->omega_int=_omega;
     if(_theta==0.){
       this->const_theta_int=false;
       this->theta_int=1.;
     }
     else{
       this->const_theta_int=true;
       this->theta_int=_theta;
     }
   }

   void startdart_fp() {this->dartOn_fp=!(this->dartOn_fp);}
   void startdart_int() {this->dartOn_int=!(this->dartOn_int);}

   ///////////////////////////////////////////////////////////////


   void settau(double tau) {pi.tau=tau;}
   tree& gettree(size_t i ) { return t[i];}
   xinfo& getxinfo() {return xi;}
   void setxinfo(xinfo& _xi);
   std::vector<size_t>& getnv() {return nv;}
   std::vector<double>& getpv() {return pv;}

   //New probabilites

   std::vector<size_t>& getnv_funcP() {return nv_funcP;}
   std::vector<size_t>& getnv_intvls() {return nv_intvls;}

   std::vector<double>& getpv_funcP() {return pv_funcP;}
   std::vector<double>& getpv_intvls() {return pv_intvls;}

   // Conditional Probs

   std::vector<size_t>& getnv_int_cond() {return  nv_intvls_cond; }
   std::vector<double>& getpv_int_cond() {return  pv_intvls_cond; }



   double gettheta() {return theta;}

   //new thetas
   double gettheta_fp() {return theta_fp;}
   double gettheta_inv() {return theta_int;}
   //------------------------------
   //public methods
   void birth(size_t i, size_t nid,size_t v, size_t c, double ml, double mr)
         {t[i].birth(nid,v,c,ml,mr);}
   void death(size_t i,size_t nid, double mu)
         {t[i].death(nid,mu);}
   void pr();
   void tonull() {for(size_t i=0;i!=t.size();i++) t[i].tonull();}
   void predict(size_t p, size_t n, double *x, double *fp);
   void draw(double sigma, rn& gen);

   //NEW Draw method

   void draw_global_conditional(size_t funcP, size_t intvls, double sigma,
                                rn& gen, rn& gen_fp, rn& gen_int);

   void draw_global(size_t funcP, size_t intvls, double sigma,
                          rn& gen, rn& gen_fp, rn& gen_int);

   void draw_global_neighborhood(size_t funcP, size_t intvls, double sigma,
                    rn& gen, rn& gen_fp, rn& gen_int, double sigma_int);



//   void draw_s(rn& gen);
   double f(size_t i) {return allfit[i];}
protected:
   size_t m;  //number of trees
   std::vector<tree> t; //the trees
   pinfo pi; //prior and mcmc info
   //data
   size_t p,n; //x has dim p, n obserations
   double *x,*y;  //x is column stack, pxn
   xinfo xi; //cutpoint info
   //working
   double *allfit; //if the data is set, should be f(x)
   double *r;
   double *ftemp;
   dinfo di;
   bool dart,dartOn,aug,const_theta;
   double a,b,rho,theta,omega;

   //dart for functions and intervals
   bool dart_fp,dartOn_fp,aug_fp,const_theta_fp;
   double a_fp,b_fp,rho_fp,theta_fp,omega_fp;

   bool dart_int,dartOn_int,aug_int,const_theta_int;
   double a_int,b_int,rho_int,theta_int,omega_int;



   std::vector<size_t> nv;
   std::vector<double> pv, lpv;

   //Add

   //creating new probabilities for interval and func_predictors

   std::vector<size_t> nv_funcP;
   std::vector<size_t> nv_intvls;

   std::vector<double> pv_funcP, lpv_funcP;
   std::vector<double> pv_intvls, lpv_intvls;

   //Creating conditional intervals
   std::vector<size_t>  nv_intvls_cond;

   std::vector<double>  pv_intvls_cond;
   std::vector<double>  lpv_intvls_cond;
};

#endif
