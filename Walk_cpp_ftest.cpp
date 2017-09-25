#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <fstream>
#include <time.h>
using namespace Rcpp;

//Genera posicion al asar dimension d escala scale
NumericVector do_muller_c(int d,double scale =2.0){
  NumericVector v = rnorm(d);
  v = v/sqrt(sum(v*v))*scale;
  return(v);
}

//Genera posicion del optimo
NumericVector optimum_rn_c(NumericVector opt_i,double sg){
  int i = opt_i.size();
  double d=Rf_rnorm(0,sg);
  NumericVector v = do_muller_c(i,d);
  //NumericVector v = rnorm(i,sg)/sqrt(i);
  return(opt_i+v);
}

//Calcula el fitness
double fit_gauss(NumericVector vec) {
  return(exp((-sum(vec*vec))/2));
}

//Calcula el fitness relativo
double rel_fit(NumericVector vecwt,NumericVector vecmt){
  return(fit_gauss(vecmt)/fit_gauss(vecwt)-1);
}

//Calcula prob y decide si se fija N infinito
bool prob_fix_ndel(double s){
  
  double a = as<double>(runif(1,0.0,1.0));       
  double p = 1-exp(-2*s);
  
  if( a < p){
    return(1);
  }else{
    return(0);
  }
}

//Calcula prob fijacion y decide si se fija N finito
bool prob_fix_K(double s,int N){
  
  double a = as<double>(runif(1,0.0,1.0));       
  double p = (1-exp(-2*s))/(1-exp(-2*N*s));
  
  if( a < p){
    return(1);
  }else{
    return(0);
  }
}

//Si cae dentro de la esfera se fija, si no no. Orr
bool fix_sphere(NumericVector vec,double d=2){
  double n = sqrt(sum(vec*vec));
  
  if(n < d/2){
    
    return(1);
    
  }else{
    
    return (0);
  }
}

//Calcula el paso para el cambio ambiental. Devuelve el numero de paso
//Inputs paso actual y el Tau de la dist de poisson del ambiente
int env_chg_step (double tau = NA_REAL,double paso=0){
  
  if(!R_IsNA(tau)){
    return(as<int>(rexp(1,1/tau))+paso);
  }
  else{
    return(NA_REAL);
    }
}

//Funcion de walk 

// [[Rcpp::export]]
DataFrame walk_Nfix_f(List input){
  
  int times,timed,ti;
  times=time(NULL);
  
  int N = input["N"]; //Poblacion
  int n = input["n"]; //Num dimensiones
  double mms = input["mms"]; //Mutation medium size -Tamanio mutacional medio
  int Nm = input["Nm"]; //numero de mutaciones final
  double tau = input["tau"];  //Tau ambiental
  double sg = input["sg"]; //Sigma ambiental
  
  
  // std::vector<double> c;
  // std::vector<int> f ;
  // std::vector<double> cs ;
  // std::vector<double> r ;
  // std::vector<double> tm ;
  // std::vector<double> tf ;
  // std::vector<double> s ;
  // std::vector<double> w ;
  // std::vector<double> p ;
  // std::vector<double> p_Ni ;
  // std::vector<std::string> type ;
  
  std::cout<<"start, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<std::endl;
  
  NumericVector vec(n); 
  //vec = do_muller_c(n,as<double>(runif(1,0,mms)));
  vec[0]=-0.5;
  NumericVector opt(n,0.0);
  //opt= optimum_rn_c(opt,sg);
  
  NumericVector opt_in = opt;
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double chg_step = env_chg_step(tau);
  int i=0;
  int fi=0;
  double ps=chg_step;
  double fd=0;
  double fs=0;
 
  while(i<Nm){
    i=i+1;
    
    if(tau==tau && i==chg_step){
      
      opt = optimum_rn_c(opt_in,sg);
      chg_step = env_chg_step(tau,i);
      ps=chg_step-i;
      //std::cout<<ps<<" "<<i<<" "<<fi<<" "<<chg_step<<std::endl;
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,mms)));
    s_r = rel_fit(vec-opt,vec_az+vec-opt);
    n_p = sqrt(sum((vec-opt)*(vec-opt)));
    n_pn = sqrt(sum((vec+vec_az-opt)*(vec+vec_az-opt)));
    
    
    if(n_pn < n_p){
      // c.push_back(i);
      // f.push_back(fi);
      // cs.push_back(ps);
      // r.push_back(n_pn);
      // tm.push_back(sqrt(sum(vec_az*vec_az)));
      // s.push_back(s_r);
      // w.push_back(fit_gauss(vec-opt));
      // p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
      // p_Ni.push_back(1-exp(-2*s_r));
      // type.push_back("s");
      fs=fs+1-exp(-2*s_r);
      fd=fd+(1-exp(-2*s_r))/(1-exp(-2*N*s_r))-1+exp(-2*s_r);
    }
    else{
    //   c.push_back(i);
    //   f.push_back(fi);
    //   cs.push_back(ps);
    //   r.push_back(n_pn);
    //   tm.push_back(sqrt(sum(vec_az*vec_az)));
    //   s.push_back(s_r);
    //   w.push_back(fit_gauss(vec-opt));
    //   p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
    //   p_Ni.push_back(0);
    //   type.push_back("d");
      fd=fd+(1-exp(-2*s_r))/(1-exp(-2*N*s_r));
    }
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;
      vec = vec + vec_az;
      }
  }
  
  
  
  //double fs=0;
  //double fd=0;
  
  // for(int i=0;i<Nm;i++){
  //   
  //   fd=fd+p.at(i)-p_Ni.at(i);
  //   fs=fs+p_Ni.at(i);
  // 
  // }
  
  
  fs=fs*N/Nm;
  fd=fd*N/Nm;
  
  //std::cout<<fs<<" "<<fd<<std::endl;
  
  timed=time(NULL);
  ti=timed-times;
  
  
  
  DataFrame datos = DataFrame::create(Named("N")=N,Named("n")=n,Named("r")=mms,Named("tau")=tau,Named("sg")=sg,Named("fs")=fs,Named("fd")=(fd),Named("ft")=(fd+fs),Named("Nf")=fi);
  
  timed=time(NULL);
  ti=timed-times;
  
  std::cout<<"finish, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<" :"<<ti<<std::endl;
  
  return(datos);
  
}


/*** R

*/
