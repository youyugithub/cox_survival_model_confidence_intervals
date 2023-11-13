# cox_survival_model_confidence_intervals
cox_survival_model_confidence_intervals

## fast_coxph.cpp
```
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

List fast_coxph(
    arma::vec time,
    arma::uvec event,
    arma::mat xx){
  
  List result_temp(10);
  
  int nind=xx.n_rows;
  int ndim=xx.n_cols;
  arma::vec SS(ndim);
  arma::mat II(ndim,ndim);
  double denominator;
  arma::vec numerator1(ndim);
  arma::mat numerator2(ndim,ndim);
  int iter,ii,jj,ll;
  
  arma::vec beta(ndim);
  
  arma::vec expr(nind);
  arma::field<arma::vec> exprx(nind);
  arma::field<arma::mat> exprx2(nind);
  
  arma::field<arma::vec> F_xx(nind);
  for(ii=0;ii<nind;ii++){
    F_xx(ii)=xx.row(ii).t();
  }

  for(iter=0;iter<300;iter++){
    SS.fill(0.0);
    II.fill(0.0);
    expr=arma::exp(xx*beta);
    for(ii=0;ii<nind;ii++){
      exprx(ii)=expr(ii)*F_xx(ii);
      exprx2(ii)=expr(ii)*F_xx(ii)*F_xx(ii).t();
    }
    for(ii=0;ii<nind;ii++){
      if(!event(ii))continue;
      SS=SS+F_xx(ii);
      denominator=0.0;
      numerator1.fill(0.0);
      numerator2.fill(0.0);
      for(ll=0;ll<nind;ll++){
        if(time(ll)<time(ii))continue;
        denominator=denominator+expr(ll);
        numerator1=numerator1+exprx(ll);
        numerator2=numerator2+exprx2(ll);
      }
      SS=SS-numerator1/denominator;
      II=II+numerator2/denominator-
        (numerator1/denominator)*(numerator1/denominator).t();
    }
    beta=beta+arma::solve(II,SS);
  }
  
  expr=arma::exp(xx*beta);
  arma::vec unique_time=time.elem(arma::find(event));
  unique_time=arma::sort(unique_time);
  int n_unique_time=unique_time.n_elem;
  arma::vec hazard0(n_unique_time);
  arma::vec hazard0_pulled(n_unique_time);
  arma::vec nevent(n_unique_time);
  arma::vec sum_expr(n_unique_time);
  for(jj=0;jj<n_unique_time;jj++){
    nevent(jj)=0.0;
    sum_expr(jj)=0.0;
    for(ii=0;ii<nind;ii++){
      if(time(ii)<unique_time(jj))continue;
      if(time(ii)==unique_time(jj))nevent(jj)=nevent(jj)+1.0;
      sum_expr(jj)=sum_expr(jj)+expr(ii);
    }
  }
  hazard0=1.0/sum_expr;
  hazard0_pulled=nevent/sum_expr;
  
  arma::mat IIbb(ndim,ndim);
  arma::mat IIbh(ndim,n_unique_time);
  arma::vec IIhh(n_unique_time);
  IIhh=nevent%sum_expr%sum_expr;
  
  IIbb.fill(0.0);
  for(jj=0;jj<n_unique_time;jj++){
    for(ii=0;ii<nind;ii++){
      if(time(ii)<unique_time(jj))continue;
      IIbb=IIbb+hazard0(jj)*exprx2(ii);
    }
  }
  
  IIbh.fill(0.0);
  for(jj=0;jj<n_unique_time;jj++){
    for(ii=0;ii<nind;ii++){
      if(time(ii)<unique_time(jj))continue;
      IIbh.col(jj)=IIbh.col(jj)-nevent(jj)*exprx(ii);
    }
  }
  
  List result(10);
  result(0)=beta;
  result(1)=II;
  result(2)=unique_time;
  result(3)=hazard0_pulled;
  result(4)=nevent;
  result(5)=IIhh;
  result(6)=IIbb;
  result(7)=IIbh;
  return(result);
}
```

## R
```
library(Rcpp)
sourceCpp("Rcpp/fast_coxph.cpp")
set.seed(0)
mytime<-rexp(100)
myevent<-sample(c(T,F),100,replace=T)
myxx<-matrix(rnorm(200),100,2)
x1<-myxx[,1]
x2<-myxx[,2]

a_fast_coxph<-fast_coxph(mytime,myevent,myxx)

calculate_var_b_H<-function(a_fast_coxph,a_time){
  unique_time<-c(a_fast_coxph[[3]])
  
  IIhh<-diag(c(a_fast_coxph[[6]]))
  IIbb<-a_fast_coxph[[7]]
  IIbh<-a_fast_coxph[[8]]
  dimb<-ncol(IIbb);idxb<-1:dimb
  dimh<-ncol(IIhh);idxh<-dimb+(1:dimh)
  II<-rbind(
    cbind(IIbb,IIbh),
    cbind(t(IIbh),IIhh))
  vv<-solve(II)
  vvbb<-vv[idxb,idxb,drop=F]
  vvbh<-vv[idxb,idxh,drop=F]
  vvhh<-vv[idxh,idxh,drop=F]
  idxt<-unique_time<=a_time
  
  vvbH<-rowSums(vvbh[,idxt,drop=F])
  vvHH<-sum(vvhh[idxt,idxt,drop=F])
  
  return(list(vvbb=vvbb,vvbH=vvbH,vvHH=vvHH))
}

a_fast_coxph<-fast_coxph(mytime,myevent,myxx)

beta<-c(a_fast_coxph[[1]])
haz<-c(a_fast_coxph[[4]])
cumhaz<-cumsum(haz)
unique_time<-c(a_fast_coxph[[3]])

# CI for exp(-H*exp(beta*x))

a_time<-0.5
a_x<-c(0.5,0.5)
a_1<-as.numeric(unique_time<=a_time)

var_b_H<-calculate_var_b_H(a_fast_coxph,a_time=a_time)
vvbb<-var_b_H$vvbb
vvbH<-var_b_H$vvbH
vvHH<-var_b_H$vvHH

est_logH_bx<-log(sum(a_1*haz))+sum(a_x*beta)
est_H<-sum(a_1*haz)
se_logH_bx<-
  sqrt(c(t(a_x)%*%vvbb%*%a_x)+
         c(2*t(a_x)%*%vvbH)/est_H+
         vvHH/est_H/est_H)

exp(-exp(est_logH_bx))
exp(-exp(est_logH_bx+c(1,-1)*qnorm(0.975)*se_logH_bx))
```
