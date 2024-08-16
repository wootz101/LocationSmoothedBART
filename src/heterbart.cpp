/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
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

#include "heterbart.h"

//--------------------------------------------------
void heterbart::pr()
{
   cout << "+++++heterbart object:\n";
   bart::pr();
}
//--------------------------------------------------
void heterbart::draw(double *sigma, rn& gen)
{
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      heterbd(t[j],xi,di,pi,sigma,nv,pv,aug,gen);
      heterdrmu(t[j],xi,di,pi,sigma,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
   }
   if(dartOn) {
     draw_s(nv,lpv,theta,gen);
     draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
     for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
   }
}




void heterbart::draw_conditional(size_t funcP, size_t intvls, double *sigma,
                                 rn& gen, rn& gen_fp, rn& gen_int)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      r[k] = y[k]-allfit[k];
    }
    heterbd(t[j],xi,di,pi,sigma,nv,pv,aug,gen);
    heterdrmu(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
  }


  if(dartOn) {

    // ADD Fucntional Selection Portion 'draw_r'
    //DART for FunctionalPredictors FP
    //if(dartOn_fp){
    for(size_t j=0;j<funcP;j++){
      int temp_funcP = {0};
      for(size_t k=0; k<intvls; k++){
        temp_funcP += nv[(j)*intvls+k];
        //temp_funcP += nv[(j)*intvls+k];
      }
      nv_funcP[j] = temp_funcP;
    }

    draw_s(nv_funcP,lpv_funcP,theta_fp,gen_fp);
    draw_theta0(const_theta_fp,theta_fp,lpv_funcP,a_fp,b_fp,rho_fp,gen_fp);

    for(size_t j=0;j<funcP;j++) pv_funcP[j]=::exp(lpv_funcP[j]);

    //}

    //DART for intervals intvls


    for(size_t k=0; k<funcP; k++){
      for(size_t j=0;j<intvls;j++){
        nv_intvls[j] = nv[(j)+k*intvls];

      }
      //draw for intervals accross all K functions (CONDITIONAL PART) for function K

      draw_s(nv_intvls, lpv_intvls, theta_int, gen_int);
      draw_theta0(const_theta_int,theta_int, lpv_intvls, a_int,b_int,rho_int,gen_int);

      //updated probabilites for function K with J intervals
      //update the new values
      for(size_t j=0;j<intvls;j++){

        nv_intvls_cond[(j)+k*intvls] = nv_intvls[j];
        pv_intvls_cond[(j)+k*intvls] = ::exp(lpv_intvls[j]);
        //collect all pv's for the final probabilites
        pv[(k)*intvls+j] = ::exp( lpv_funcP[k] + lpv_intvls[j]);
      }
    }



  }



}




void heterbart::draw_conditional_neighborhood(size_t funcP, size_t intvls, double *sigma,
                                 rn& gen, rn& gen_fp, rn& gen_int, double sigma_int)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      r[k] = y[k]-allfit[k];
    }
    heterbd(t[j],xi,di,pi,sigma,nv,pv,aug,gen);
    heterdrmu(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
  }


  if(dartOn) {

    // ADD Fucntional Selection Portion 'draw_r'
    //DART for FunctionalPredictors FP
    //if(dartOn_fp){
    for(size_t j=0;j<funcP;j++){
      int temp_funcP = {0};
      for(size_t k=0; k<intvls; k++){
        temp_funcP += nv[(j)*intvls+k];
        //temp_funcP += nv[(j)*intvls+k];
      }
      nv_funcP[j] = temp_funcP;
    }

    draw_s(nv_funcP,lpv_funcP,theta_fp,gen_fp);
    draw_theta0(const_theta_fp,theta_fp,lpv_funcP,a_fp,b_fp,rho_fp,gen_fp);

    for(size_t j=0;j<funcP;j++) pv_funcP[j]=::exp(lpv_funcP[j]);

    //}

    //DART for intervals intvls


    for(size_t k=0; k<funcP; k++){

      for(size_t j=0;j<intvls;j++){
        nv_intvls[j] = nv[(j)+k*intvls];

      }
      //draw for intervals across all K functions (CONDITIONAL PART) for function K

      // WANT TO IMPLEMENT GAUSSIAN PROCESS HERE

      draw_s(nv_intvls, lpv_intvls, theta_int, gen_int);
      draw_theta0(const_theta_int,theta_int, lpv_intvls, a_int,b_int,rho_int,gen_int);


      //get neighborhood matrix
      std::vector<std::vector<double> > temp_mat(intvls, std::vector<double>(intvls));
      for(size_t j=0;j<intvls;j++){
        int temp_diff = {0};
        for(size_t l=0;l<intvls;l++){
          //temp_mat[j][k] = ::exp(lpv_intvls[j])*::exp( ( (-(k-j)^2 )/2 ) );
          temp_diff = l-j;
          temp_mat[j][l] =  ::exp(lpv_intvls[j])*::exp( -temp_diff*temp_diff/ (2*sigma_int*sigma_int) );

          //cout << "lpv_intvls[j]" << ::exp(lpv_intvls[j]) << '\n' << endl;
          //cout << "sigma_int " << sigma_int << '\n' << endl;
          //cout << "(2*sigma_int*sigma_int)" << (2*sigma_int*sigma_int) << '\n' << endl;

          //cout << "j k " << j << " "<< l << " "<< temp_diff << " "<< -temp_diff*temp_diff<< '\n' << endl;
          //cout << "Matrix Value " << temp_mat[j][k] << '\n' << endl;
        }
      }

    //Normalize the probabilities
      double temp_sum = {0};
      for(size_t j=0;j<intvls;j++){
        nv_intvls_cond[(j)+k*intvls] = nv_intvls[j];
        for(size_t l=0;l<intvls;l++){
          pv_intvls_cond[(j)+k*intvls] += temp_mat[l][j];
        }
        temp_sum += pv_intvls_cond[(j)+k*intvls];
      }

      //cout << "TEMP SUM VALUE " << temp_sum << '\n' << endl;

      /*
      cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: "
           << dart << ',' << theta << ',' << omega << ',' << a << ','
           << b << ',' << rho << ',' << aug << endl;
      printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n"
              );
       */


      for(size_t j=0;j<intvls;j++){
        pv_intvls_cond[(j)+k*intvls] /= temp_sum;
        //collect all pv's for the final probabilites
        pv[(k)*intvls+j] = ::exp( lpv_funcP[k])*pv_intvls_cond[(j)+k*intvls];
      }



      //update log probabilitiies
      /*
      for(size_t j=0;j<intvls;j++){
        for(size_t l=0;l<intvls;l++){
          lpv_intvls[j] = temp_mat[l][j] + lpv_intvls[j] ;
        }
      }
       */


      //update probabilities
     /* double temp_sum{0};
      for(size_t j=0;j<intvls;j++){

        pv_intvls[j]=::exp(lpv_intvls[j]);
        temp_sum += pv_intvls[j];
      }

      //normalize probabilities
      for(size_t j=0;j<intvls;j++){
        pv_intvls[j] /= temp_sum;
      }

      */

      //updated probabilites for function K with J intervals
      //update the new values

      /*
      double temp_sum{0};
      for(size_t j=0;j<intvls;j++){

        nv_intvls_cond[(j)+k*intvls] = nv_intvls[j];
        pv_intvls_cond[(j)+k*intvls] = ::exp(lpv_intvls[j]);
        for(size_t l=0;l<intvls;l++){
          pv_intvls_cond[(j)+k*intvls] += temp_mat[l][j];
        }
        temp_sum += pv_intvls_cond[(j)+k*intvls];
      }


      for(size_t j=0;j<intvls;j++){
        pv_intvls_cond[(j)+k*intvls] /= temp_sum;
        //collect all pv's for the final probabilites
        pv[(k)*intvls+j] = ::exp( lpv_funcP[k])*pv_intvls_cond[(j)+k*intvls];
      }
       */




    }
  }
}

