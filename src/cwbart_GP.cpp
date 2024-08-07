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

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cwbart_GP(
    SEXP _funcP,
    SEXP _intvls,




    SEXP _in,            //number of observations in training data
    SEXP _ip,		//dimension of x
    SEXP _inp,		//number of observations in test data
    SEXP _ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP _iy,		//y, train,  nx1
    SEXP _ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP _im,		//number of trees
    SEXP _inc,		//number of cut points
    SEXP _ind,		//number of kept draws (except for thinnning ..)
    SEXP _iburn,		//number of burn-in draws skipped
    SEXP _ipower,
    SEXP _ibase,
    SEXP _itau,
    SEXP _inu,
    SEXP _ilambda,
    SEXP _isigest,
    SEXP _iw,
    SEXP _idart,
    SEXP _itheta,
    SEXP _iomega,
    SEXP _igrp,
    SEXP _ia,
    SEXP _ib,
    SEXP _irho,
    SEXP _iaug,

    //NEW DART Variables
    SEXP _idart_fp,
    SEXP _itheta_fp,
    SEXP _iomega_fp,
    //SEXP _igrp,
    SEXP _ia_fp,
    SEXP _ib_fp,
    SEXP _irho_fp,
    SEXP _iaug_fp,

    SEXP _idart_int,
    SEXP _itheta_int,
    SEXP _iomega_int,
    //SEXP _igrp,
    SEXP _ia_int,
    SEXP _ib_int,
    SEXP _irho_int,
    SEXP _iaug_int,
    // For GP process
    //SEXP _gp_bool,




    SEXP _inkeeptrain,
    SEXP _inkeeptest,
    SEXP _inkeeptestme,
    SEXP _inkeeptreedraws,
    SEXP _inprintevery,
    //   SEXP _treesaslists,
    SEXP _Xinfo,
    SEXP _isigma_int


)
{

  //--------------------------------------------------
  //process args
  size_t funcP = Rcpp::as<int>(_funcP);
  size_t intvls = Rcpp::as<int>(_intvls);

  size_t n = Rcpp::as<int>(_in);
  size_t p = Rcpp::as<int>(_ip);
  size_t np = Rcpp::as<int>(_inp);
  Rcpp::NumericVector  xv(_ix);
  double *ix = &xv[0];
  Rcpp::NumericVector  yv(_iy);
  double *iy = &yv[0];
  Rcpp::NumericVector  xpv(_ixp);
  double *ixp = &xpv[0];
  size_t m = Rcpp::as<int>(_im);
  Rcpp::IntegerVector _nc(_inc);
  int *numcut = &_nc[0];
  //size_t nc = Rcpp::as<int>(_inc);
  size_t nd = Rcpp::as<int>(_ind);
  size_t burn = Rcpp::as<int>(_iburn);
  double mybeta = Rcpp::as<double>(_ipower);
  double alpha = Rcpp::as<double>(_ibase);
  double tau = Rcpp::as<double>(_itau);
  double nu = Rcpp::as<double>(_inu);
  double lambda = Rcpp::as<double>(_ilambda);
  double sigma=Rcpp::as<double>(_isigest);
  Rcpp::NumericVector  wv(_iw);
  double *iw = &wv[0];
  bool dart;
  if(Rcpp::as<int>(_idart)==1) dart=true;
  else dart=false;
  double a = Rcpp::as<double>(_ia);
  double b = Rcpp::as<double>(_ib);
  double rho = Rcpp::as<double>(_irho);
  bool aug;
  if(Rcpp::as<int>(_iaug)==1) aug=true;
  else aug=false;
  double theta = Rcpp::as<double>(_itheta);
  double omega = Rcpp::as<double>(_iomega);

  //New DART variables

  bool dart_int;
  if(Rcpp::as<int>(_idart_int)==1) dart_int=true;
  else dart_int=false;
  double a_int = Rcpp::as<double>(_ia_int);
  double b_int = Rcpp::as<double>(_ib_int);
  double rho_int = Rcpp::as<double>(_irho_int);
  bool aug_int;
  if(Rcpp::as<int>(_iaug_int)==1) aug_int=true;
  else aug_int=false;
  double theta_int = Rcpp::as<double>(_itheta_int);
  double omega_int = Rcpp::as<double>(_iomega_int);

  double sigma_int = Rcpp::as<double>(_isigma_int);


  bool dart_fp;
  if(Rcpp::as<int>(_idart_fp)==1) dart_fp=true;
  else dart_fp=false;
  double a_fp = Rcpp::as<double>(_ia_fp);
  double b_fp = Rcpp::as<double>(_ib_fp);
  double rho_fp = Rcpp::as<double>(_irho_fp);
  bool aug_fp;
  if(Rcpp::as<int>(_iaug_fp)==1) aug_fp=true;
  else aug_fp=false;
  double theta_fp = Rcpp::as<double>(_itheta_fp);
  double omega_fp = Rcpp::as<double>(_iomega_fp);

  //GP
 // bool gp_bool;
  //if(Rcpp::as<int>(_gp_bool)==1) gp_bool=true;
  //else gp_bool=false;



  Rcpp::IntegerVector _grp(_igrp);
  int *grp = &_grp[0];
  size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
  size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
  size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
  size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
  size_t printevery = Rcpp::as<int>(_inprintevery);
  //   int treesaslists = Rcpp::as<int>(_treesaslists);
  Rcpp::NumericMatrix Xinfo(_Xinfo);
  //   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);




  //return data structures (using Rcpp)
  Rcpp::NumericVector trmean(n); //train
  Rcpp::NumericVector temean(np);
  Rcpp::NumericVector sdraw(nd+burn);
  Rcpp::NumericMatrix trdraw(nkeeptrain,n);
  Rcpp::NumericMatrix tedraw(nkeeptest,np);
  //   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
  Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
  Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);


  Rcpp::NumericMatrix varprb_fp(nkeeptreedraws,funcP);
  Rcpp::IntegerMatrix varcnt_fp(nkeeptreedraws,funcP);


  Rcpp::NumericMatrix varprb_int(nkeeptreedraws, p);
  Rcpp::IntegerMatrix varcnt_int(nkeeptreedraws, p);


  //random number generation
  arn gen;

  arn gen_fp;
  arn gen_int;

  heterbart bm(m);

  if(Xinfo.size()>0) {
    xinfo _xi;
    _xi.resize(p);
    for(size_t i=0;i<p;i++) {
      _xi[i].resize(numcut[i]);
      //Rcpp::IntegerVector cutpts(Xinfo[i]);
      for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
    }
    bm.setxinfo(_xi);
  }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

  void cwbart(
      funcP,
      intvls,


      size_t n,            //number of observations in training data
      size_t p,		//dimension of x
      size_t np,		//number of observations in test data
      double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
      double* iy,		//y, train,  nx1
      double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
      size_t m,		//number of trees
      int* numcut,		//number of cut points
      size_t nd,		//number of kept draws (except for thinnning ..)
      size_t burn,		//number of burn-in draws skipped
      double mybeta,
      double alpha,
      double tau,
      double nu,
      double lambda,
      double sigma,
      double* iw,
      bool dart,
      double theta,
      double omega,
      int *grp,
      double a,
      double b,
      double rho,
      bool aug,

      //New DART values
      bool dart_fp,
      double theta_fp,
      double omega_fp,
      //int *grp,
      double a_fp,
      double b_fp,
      double rho_fp,
      bool aug_fp,

      bool dart_int,
      double theta_int,
      double omega_int,
      //int *grp,
      double a_int,
      double b_int,
      double rho_int,
      bool aug_int,

      //bool gp_bool;



  size_t nkeeptrain,
  size_t nkeeptest,
  size_t nkeeptestme,
  size_t nkeeptreedraws,
  size_t printevery,

  double sigma_int,
  //   int treesaslists,
  unsigned int n1, // additional parameters needed to call from C++
  unsigned int n2,
  double* trmean,
  double* temean,
  double* sdraw,
  double* _trdraw,
  double* _tedraw



  )
  {

    //return data structures (using C++)
    std::vector<double*> trdraw(nkeeptrain);
    std::vector<double*> tedraw(nkeeptest);

    for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
    for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

    std::vector< std::vector<size_t> > varcnt;
    std::vector< std::vector<double> > varprb;

    //conditional DART iterations

    //For conditional DART
    std::vector< std::vector<size_t>> varcnt_fp;
    std::vector< std::vector<double>> varprb_fp;

    std::vector< std::vector<size_t>> varcnt_int;
    std::vector< std::vector<double>> varprb_int;


    //random number generation
    arn gen(n1, n2);

    arn gen_fp(n1, n2);
    arn gen_int(n1, n2);

    heterbart bm(m);
#endif

    for(size_t i=0;i<n;i++) trmean[i]=0.0;
    for(size_t i=0;i<np;i++) temean[i]=0.0;

    printf("*****Into main of wbart\n");
    //-----------------------------------------------------------

    size_t skiptr,skipte,skipteme,skiptreedraws;
    if(nkeeptrain) {skiptr=nd/nkeeptrain;}
    else skiptr = nd+1;
    if(nkeeptest) {skipte=nd/nkeeptest;}
    else skipte=nd+1;
    if(nkeeptestme) {skipteme=nd/nkeeptestme;}
    else skipteme=nd+1;
    if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
    else skiptreedraws=nd+1;

    //--------------------------------------------------
    //print args
    printf("*****Data:\n");
    printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
    printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
    printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
    if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
    printf("*****Number of Trees: %zu\n",m);
    printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
    printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
    printf("*****Prior:beta,alpha,tau,nu,lambda: %lf,%lf,%lf,%lf,%lf\n",
           mybeta,alpha,tau,nu,lambda);
    printf("*****sigma: %lf\n",sigma);
    printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
    cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: "
         << dart << ',' << theta << ',' << omega << ',' << a << ','
         << b << ',' << rho << ',' << aug << endl;
    printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
           nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
    printf("*****printevery: %zu\n",printevery);
    printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);

    //--------------------------------------------------
    //heterbart bm(m);
    bm.setprior(alpha,mybeta,tau);
    //bm.setdata(p,n,ix,iy,numcut);
    bm.setdata_fbart(p,n, funcP, intvls, ix,iy,numcut);
    bm.setdart(a,b,rho,aug,dart,theta,omega);

    bm.setdart_fp(a_fp,b_fp,rho_fp,aug_fp,dart_fp,theta_fp,omega_fp);
    bm.setdart_int(a_int,b_int,rho_int,aug_int,dart_int,theta_int,omega_int);




    //--------------------------------------------------
    //sigma
    //gen.set_df(n+nu);
    double *svec = new double[n];
    for(size_t i=0;i<n;i++) svec[i]=iw[i]*sigma;

    //--------------------------------------------------

    std::stringstream treess;  //string stream to write trees to
    treess.precision(10);
    treess << nkeeptreedraws << " " << m << " " << p << endl;
    // dart iterations
    std::vector<double> ivarprb (p,0.);
    std::vector<size_t> ivarcnt (p,0);

    //Other DART
    std::vector<double> ivarprb_fp (funcP,0.);
    std::vector<size_t> ivarcnt_fp (funcP,0);

    std::vector<double> ivarprb_int (intvls,0.);
    std::vector<size_t> ivarcnt_int (intvls,0);

    //--------------------------------------------------
    //temporary storage
    //out of sample fit
    double* fhattest=0; //posterior mean for prediction
    if(np) { fhattest = new double[np]; }
    double restemp=0.0,rss=0.0;


    //--------------------------------------------------
    //mcmc
    printf("\nMCMC\n");
    //size_t index;
    size_t trcnt=0; //count kept train draws
    size_t tecnt=0; //count kept test draws
    size_t temecnt=0; //count test draws into posterior mean
    size_t treedrawscnt=0; //count kept bart draws
    bool keeptest,keeptestme,keeptreedraw;

    time_t tp;
    int time1 = time(&tp);
    xinfo& xi = bm.getxinfo();

    for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) {
        bm.startdart();
        bm.startdart_fp();
        bm.startdart_int();

      }
      //draw bart
      //bm.draw(svec,gen);

      //printf("sigma_int %zu ", sigma_int);
      //bm.draw_conditional_neighborhood(funcP, intvls, svec, gen, gen_fp, gen_int, sigma_int);


        bm.draw_GP(funcP, intvls, svec, gen, gen_fp, gen_int, sigma_int);



      //draw sigma
      rss=0.0;
      for(size_t k=0;k<n;k++) {restemp=(iy[k]-bm.f(k))/(iw[k]); rss += restemp*restemp;}
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
      for(size_t k=0;k<n;k++) svec[k]=iw[k]*sigma;
      sdraw[i]=sigma;
      if(i>=burn) {
        for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
        if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
          //index = trcnt*n;;
          //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
          for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
          trcnt+=1;
        }
        keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
        keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
        if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
        if(keeptest) {
          //index=tecnt*np;
          //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
          for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
          tecnt+=1;
        }
        if(keeptestme) {
          for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
          temecnt+=1;
        }
        keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
        if(keeptreedraw) {
          //	   #ifndef NoRcpp
          //	   Rcpp::List lists(m*treesaslists);
          //	   #endif

          for(size_t j=0;j<m;j++) {
            treess << bm.gettree(j);
            /*
#ifndef NoRcpp
             varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
             if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
#endif
             */
          }
#ifndef NoRcpp
          //	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
          ivarcnt=bm.getnv();
          ivarprb=bm.getpv();

          ivarcnt_fp =bm.getnv_funcP();
          ivarprb_fp =bm.getpv_funcP();

          ivarcnt_int = bm.getnv_int_cond();
          ivarprb_int = bm.getpv_int_cond();



          size_t k=(i-burn)/skiptreedraws;
          for(size_t j=0;j<p;j++){
            varcnt(k,j)=ivarcnt[j];
            //varcnt(i-burn,j)=ivarcnt[j];
            varprb(k,j)=ivarprb[j];
            //varprb(i-burn,j)=ivarprb[j];
            varcnt_int(k,j)=ivarcnt_int[j];
            //varcnt(i-burn,j)=ivarcnt[j];
            varprb_int(k,j)=ivarprb_int[j];
          }

          for(size_t j=0;j<funcP;j++){
            varcnt_fp(k,j)=ivarcnt_fp[j];
            //varcnt(i-burn,j)=ivarcnt[j];
            varprb_fp(k,j)=ivarprb_fp[j];
            //varprb(i-burn,j)=ivarprb[j];
          }


#else
          varcnt.push_back(bm.getnv());
          varprb.push_back(bm.getpv());

          varcnt_fp.push_back(bm.getnv_funcP());
          varprb_fp.push_back(bm.getpv_funcP());

          varcnt_int.push_back(bm.getnv_int_cond());
          varprb_int.push_back(bm.getpv_int_cond());
#endif

          treedrawscnt +=1;
        }
      }
    }
    int time2 = time(&tp);
    printf("time: %ds\n",time2-time1);
    for(size_t k=0;k<n;k++) trmean[k]/=nd;
    for(size_t k=0;k<np;k++) temean[k]/=temecnt;
    printf("check counts\n");
    printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
    //--------------------------------------------------
    //PutRNGstate();

    if(fhattest) delete[] fhattest;
    if(svec) delete [] svec;

    //--------------------------------------------------
    //return
#ifndef NoRcpp
    Rcpp::List ret;
    ret["sigma"]=sdraw;
    ret["yhat.train.mean"]=trmean;
    ret["yhat.train"]=trdraw;
    ret["yhat.test.mean"]=temean;
    ret["yhat.test"]=tedraw;
    //ret["varcount"]=varcount;
    ret["varcount"]=varcnt;
    ret["varprob"]=varprb;


    ret["varcount_fp"]=varcnt_fp;
    ret["varprob_fp"]=varprb_fp;

    ret["varcount_int"]=varcnt_int;
    ret["varprob_int"]=varprb_int;



    //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
    //}

    Rcpp::List xiret(xi.size());
    for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
    }

    Rcpp::List treesL;
    //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
    //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
    //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
    treesL["cutpoints"] = xiret;
    treesL["trees"]=Rcpp::CharacterVector(treess.str());
    //   if(treesaslists) treesL["lists"]=list_of_lists;
    ret["treedraws"] = treesL;

    return ret;
#else

#endif

  }
