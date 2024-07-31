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






void heterbart::draw_GP(size_t funcP, size_t intvls, double *sigma,
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
      }
      nv_funcP[j] = temp_funcP;
    }

    draw_s(nv_funcP,lpv_funcP,theta_fp,gen_fp);
    draw_theta0(const_theta_fp,theta_fp,lpv_funcP,a_fp,b_fp,rho_fp,gen_fp);

    for(size_t j=0;j<funcP;j++) pv_funcP[j]=::exp(lpv_funcP[j]);

    //}

    //DART for intervals intvls


    for(size_t k=0; k<funcP; k++){

      VectorXd eta_l(intvls);
      VectorXd pseudo_probs(intvls);
      VectorXd locations = VectorXd::LinSpaced(intvls, 1.0 / intvls, 1.0);

      //cout << "Locations:" << endl;
      //cout << locations << endl;

      int temp_sum_1 = {0};
      for(size_t j=0;j<intvls;j++){
        nv_intvls[j] = nv[(j)+k*intvls];
        temp_sum_1 += nv_intvls[j];
      }

      //cout << "New all cuts :" << endl;
      //cout << temp_sum_1 << endl;

      double eps = 1e-6;  // A small constant to avoid numerical issues
      double total_prob = 0.0;
      // First, compute adjusted probabilities and total sum
      for(size_t j = 0; j < intvls; j++) {
        double prob = (nv_intvls[j] + 1.0) / (temp_sum_1 + intvls); // Adjust temp_sum_1 + intvls to ensure normalization
        prob = std::max(eps, std::min(1.0 - eps, prob));  // Ensure prob is within (0, 1)
        pseudo_probs[j] = prob;
        total_prob += prob; // Accumulate total probability
      }

      // Normalize the probabilities to ensure they sum to 1
      for(size_t j = 0; j < intvls; j++) {
        pseudo_probs[j] = pseudo_probs[j] / total_prob; // Normalize
        eta_l[j] = log(pseudo_probs[j] / (1 - pseudo_probs[j])); // Transform to logit scale
      }


     // for(size_t j=0;j<intvls;j++){
       // pseudo_probs[j] =  (nv_intvls[j] + 1) / (temp_sum_1 + 1); // one is added to avoid zero probability
      //  eta_l[j] = log(pseudo_probs[j] / (1 - pseudo_probs[j])); //map to real number line
      //}

      //cout << "New pseudo-prob :" << endl;
      //cout << pseudo_probs << endl;

      //cout << "New eta :" << endl;
      //cout << eta_l.head(1) << endl;

      // WANT TO IMPLEMENT GAUSSIAN PROCESS HERE

      // Gaussian Process Parameters
      double alpha = 4;
      double bandwidth = 0.05;
      double sigmasq_jitter = 0.01;
      double sigmasq_eta = 1;

      // Compute GP covariance matrix
      MatrixXd kernel_mat = covSE(locations, bandwidth, alpha, sigmasq_jitter);

      // Compute posterior covariance and mean
      MatrixXd gp_cov_post = (kernel_mat.inverse() + (1 / sigmasq_eta) * MatrixXd::Identity(intvls, intvls)).inverse();
      double prior_mean = log((1.0/intvls) / (1 - 1.0/intvls));
      VectorXd gp_mean_post = prior_mean * VectorXd::Ones(intvls) + gp_cov_post * ((1 / sigmasq_eta) * (eta_l.array() - prior_mean).matrix());

      //cout << "New gp_mean (first 2):" << endl;
      //cout << gp_mean_post.head(1) << endl;


      // Printing the determinant of kernel_mat
     // cout << "Determinant of kernel_mat: " << kernel_mat.determinant() << endl;

      // Printing the determinant of gp_cov_post
     // cout << "Determinant of gp_cov_post: " << gp_cov_post.determinant() << endl;

      // Sample from multivariate normal for the new logit probabilities
      VectorXd new_logit_probs = sampleMultivariateNormal(gp_mean_post, gp_cov_post);

      //cout << "New Logit (first 2):" << endl;
      //cout << new_logit_probs.head(2) << endl;


      // Convert back to probabilities ensuring they sum to 1
      VectorXd new_probs = new_logit_probs.array().exp() / (1 + new_logit_probs.array().exp());
      new_probs = new_probs / new_probs.sum();

      // Output some of the new probabilities
      //cout << "New Probabilities (first 2):" << endl;
      //cout << new_probs.head(2) << endl;


      // ensure we record all probabilities

      for(size_t j=0;j<intvls;j++){
        nv_intvls_cond[(j)+k*intvls] = nv_intvls[j];
        pv_intvls_cond[(j)+k*intvls] = new_probs[(j)];
        //collect all pv's for the final probabilites
        pv[(k)*intvls+j] = ::exp( lpv_funcP[k])*new_probs[(j)];
      }





    }
  }
}


