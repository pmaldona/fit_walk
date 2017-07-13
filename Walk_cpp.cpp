#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <fstream>
#include <time.h>
using namespace Rcpp;

/* Rutinas que generan caminatas evolutivas en un espacio geometrico de rasgos de n dimesiones
 * hacía el optimo mediante, llamado modelo geometrico de Fisher (Razeto 2011), estas tienen porbabilidad de 
 * fijación descritas por Kimura (kimura 1964).
 * 
 * A continuación se presenta la rutina descrita en el README.md, La idea de estas rutinas es de presentar
 * generar un input/output en forma de lista para se ejecutadas mediante algun lapply o mclapply por parLapply. 
 * Esto con el fin de generar una paralilización mediante alguna herramienta en R.
 */
 
// Función que genera un vector con dirección al azar y tamaño scala.
NumericVector do_muller_c(int d,double scale =2.0){
  NumericVector v = rnorm(d,0,scale);
  v = v/sqrt(sum(v*v))*scale;
  return(v);
}
//Función que genera la posición de del optimo 
NumericVector optimum_rn_c(NumericVector opt_i,double sg){
  int i = opt_i.size();
  NumericVector v = rnorm(i,sg)/sqrt(i);
  return(opt_i+v);
}

double fit_gauss(NumericVector vec) {
  return(exp((-sum(vec*vec))));
}

double rel_fit(NumericVector vecwt,NumericVector vecmt){
  return(fit_gauss(vecmt)/fit_gauss(vecwt)-1);
}

bool prob_fix_ndel(double s){
  
  double a = as<double>(runif(1,0.0,1.0));       
  double p = 1-exp(-2*s);
  
  if( a < p){
    return(1);
  }else{
    return(0);
  }
}

bool prob_fix_K(double s,int N){
  
  double a = as<double>(runif(1,0.0,1.0));       
  double p = (1-exp(-2*s))/(1-exp(-2*N*s));
  
  if( a < p){
    return(1);
  }else{
    return(0);
  }
}

bool fix_sphere(NumericVector vec,double d=2){
  double n = sqrt(sum(vec*vec));
  
  if(n < d/2){
    
    return(1);
    
  }else{
    
    return (0);
  }
}

// [[Rcpp::export]]
double env_chg_step (double tau = NA_REAL,double paso=0){
  
  if(!R_IsNA(tau)){
    return(as<double>(rexp(1,1/tau))+paso);
  }
  else{
    return(NA_REAL);
    }
}

// [[Rcpp::export]]
DataFrame walk_Nfix_f(List input){
  
  int times,timed,ti;
  times=time(NULL);
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double tau = input["tau"];
  double sg = input["sg"];
  
  
  std::vector<double> c (Nf);
  std::vector<int> f (Nf);
  std::vector<double> cs (Nf);
  std::vector<double> r (Nf);
  std::vector<double> tm (Nf);
  std::vector<double> s (Nf);
  std::vector<double> w (Nf);
  std::vector<double> p (Nf);
  std::vector<double> p_Ni (Nf);
  std::vector<std::string> type (Nf);
  
  std::cout<<"start, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<std::endl;
  
  NumericVector vec(n); 
  vec = do_muller_c(n,as<double>(runif(1,0,2*mms)));
  NumericVector opt(n);
  opt= optimum_rn_c(opt,sg);
  NumericVector opt_in = opt;
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double chg_step = env_chg_step(tau);
  double i=0;
  int fi=0;
  double ps=chg_step;
  
  while(fi<Nf){
    i=i+1;
    
    if(tau==tau && i>chg_step){
      
      opt = optimum_rn_c(opt_in,sg);
      chg_step = env_chg_step(tau,i);
      ps=chg_step-i;
      //std::cout<<ps<<" "<<i<<" "<<fi<<" "<<chg_step<<std::endl;
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec-opt,vec_az+vec-opt);
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;
      
      n_p = sqrt(sum((vec-opt)*(vec-opt)));
      vec = vec + vec_az;
      n_pn = sqrt(sum((vec-opt)*(vec-opt)));
      
      if(n_pn < n_p){
        c.at(fi-1)=i;
        f.at(fi-1)=fi;
        cs.at(fi-1)=ps;
        r.at(fi-1)=n_pn;
        tm.at(fi-1)=sqrt(sum(vec_az*vec_az));
        s.at(fi-1)=s_r;
        w.at(fi-1)=fit_gauss(vec-opt);
        p.at(fi-1)=(1-exp(-2*s_r))/(1-exp(-2*N*s_r));
        p_Ni.at(fi-1)=1-exp(-2*s_r);
        type.at(fi-1)="s";
      }
      else{
        c.at(fi-1)=i;
        f.at(fi-1)=fi;
        cs.at(fi-1)=ps;
        r.at(fi-1)=n_pn;
        tm.at(fi-1)=sqrt(sum(vec_az*vec_az));
        s.at(fi-1)=s_r;
        w.at(fi-1)=fit_gauss(vec-opt);
        p.at(fi-1)=(1-exp(-2*s_r))/(1-exp(-2*N*s_r));
        p_Ni.at(fi-1)=0;
        type.at(fi-1)="d";
      }
    }
  }
  
  
  double fs=0;
  double fd=0;
  
  for(int i=0;i<Nf;i++){
    
    fd=fd+p.at(i)-p_Ni.at(i);
    fs=fs+p_Ni.at(i);
  
  }
  
  
  fs=fs/c.at(Nf-1)*N*N;
  fd=fd/c.at(Nf-1)*N*N;
  
  //std::cout<<fs<<" "<<fd<<std::endl;
  
  timed=time(NULL);
  ti=timed-times;
  
  
  
  DataFrame datos = DataFrame::create(Named("N")=N,Named("n")=n,Named("r")=mms,Named("tau")=tau,Named("sg")=sg,Named("fs")=fs,Named("fd")=(fd),Named("ft")=(fd+fs));
  
  timed=time(NULL);
  ti=timed-times;
  
  std::cout<<"finish, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<" :"<<ti<<std::endl;
  
  return(datos);
  
}



// [[Rcpp::export]]
void walk_Nfix_v(List & input){
  
  int times,timed,ti;
  times=time(NULL);
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> cs;
  std::vector<double> r;
  std::vector<double> tm;
  std::vector<double> s;
  std::vector<double> w;
  std::vector<double> p;
  std::vector<double> p_Ni;
  std::vector<std::string> type;
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double tau = input["tau"];
  double sg = input["sg"];
  
  std::cout<<"start, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<std::endl;
  
  NumericVector vec(n); 
  vec = do_muller_c(n,as<double>(runif(1,0,2*mms)));
  NumericVector opt(n);
  opt= optimum_rn_c(opt,sg);
  NumericVector opt_in = opt;
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double chg_step = env_chg_step(tau);
  double i=0;
  int fi=0;
  double ps=chg_step;
  
  while(fi<Nf){
    i=i+1;
    
    if(tau==tau && i>chg_step){
      
      opt = optimum_rn_c(opt_in,sg);
      chg_step = env_chg_step(tau,i);
      ps=chg_step-i;
      std::cout<<ps<<" "<<i<<" "<<fi<<" "<<chg_step<<std::endl;
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec-opt,vec_az+vec-opt);
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;

      n_p = sqrt(sum((vec-opt)*(vec-opt)));
      vec = vec + vec_az;
      n_pn = sqrt(sum((vec-opt)*(vec-opt)));
      
      if(n_pn < n_p){
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec-opt));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(1-exp(-2*s_r));
        type.push_back("s");
      }
      
      else{
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec-opt));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(0);
        type.push_back("d");
      }
      
    }
  }
  
  timed=time(NULL);
  ti=timed-times;
  
  std::ostringstream oss;
  
  oss<<"walk_N"<<N<<"_n"<<n<<"_r"<<mms<<"_tau"<<tau<<"_sg"<<sg<<".csv";
  std::string nomarch = oss.str();
  
  std::ofstream arch(nomarch.c_str());
  
  arch<<"c,f,cs,r,tm,s,w,p,p_Ni,type"<<std::endl;
  for(int k=0; k<Nf ;k++){
    arch<<c.at(k)<<","<<f.at(k)<<","<<cs.at(k)<<","<<r.at(k)<<","<<tm.at(k)<<","<<s.at(k)<<","<<w.at(k)<<","<<p.at(k)<<","<<p_Ni.at(k)<<","<<type.at(k)<<std::endl;
  }
  
  arch.close();
  
  std::cout<<"finish, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<" :"<<ti<<std::endl;
  
}

// [[Rcpp::export]]
void walk_Nfix_v_T(List & input){
  
  int times,timed,ti;
  times=time(NULL);
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> cs;
  std::vector<double> r;
  std::vector<double> tm;
  std::vector<double> s;
  std::vector<double> w;
  std::vector<double> p;
  std::vector<double> p_Ni;
  std::vector<std::string> type;
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double sg = input["sg"];
  double nd = (double) n;
  double Nd = (double) N;
  double gl = pow(1-1/(2*Nd-1),nd/2);
  
  std::cout<<"start glt, N="<<N<<" n="<<n<<" r="<<mms<<" gl="<<gl<<std::endl;
  
  
  NumericVector vec(n); 
  vec = do_muller_c(n,1.0);
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double i=0;
  int fi=0;
  double ps=0;
  double pa;
  while(fi<Nf){
    i=i+1;
    
    if(gl<(fit_gauss(vec))){
      //std::cout<<ps<<" "<<i<<" "<<fi<<" "<<(fit_gauss(vec))<<" "<<gl<<std::endl;
      vec = do_muller_c(n,1.0);
      ps=pa-i;
      pa=i;
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec,vec_az+vec);
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;
      
      n_p = sqrt(sum((vec)*(vec)));
      vec = vec + vec_az;
      n_pn = sqrt(sum((vec)*(vec)));
      
      if(n_pn < n_p){
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(1-exp(-2*s_r));
        type.push_back("s");
      }
      
      else{
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(0);
        type.push_back("d");
      }
      
    }
  }
  
  timed=time(NULL);
  ti=timed-times;
  
  std::ostringstream oss;
  
  oss<<"walk_glt_N"<<N<<"_n"<<n<<"_r"<<mms<<"_sg"<<sg<<"_glt"<<gl<<".csv";
  std::string nomarch = oss.str();
  
  std::ofstream arch(nomarch.c_str());
  
  arch<<"c,f,cs,r,tm,s,w,p,p_Ni,type"<<std::endl;
  for(int k=0; k<Nf ;k++){
    arch<<c.at(k)<<","<<f.at(k)<<","<<cs.at(k)<<","<<r.at(k)<<","<<tm.at(k)<<","<<s.at(k)<<","<<w.at(k)<<","<<p.at(k)<<","<<p_Ni.at(k)<<","<<type.at(k)<<std::endl;
  }
  
  arch.close();
  
  std::cout<<"finish glt, N="<<N<<" n="<<n<<" r="<<mms<<" sg="<<sg<<" :"<<ti<<std::endl;
  
}

// [[Rcpp::export]]
void walk_Nfix_v_O(List & input){
  
  int times,timed,ti;
  times=time(NULL);
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> cs;
  std::vector<double> r;
  std::vector<double> tm;
  std::vector<double> s;
  std::vector<double> w;
  std::vector<double> p;
  std::vector<double> p_Ni;
  std::vector<std::string> type;
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double nd = (double) n;
  double Nd = (double) N;
  double gl = 1-nd/(nd+2*Nd);
  std::cout<<"start glo, N="<<N<<" n="<<n<<" r="<<mms<<" gl="<<gl<<std::endl;
  
  NumericVector vec(n); 
  vec = do_muller_c(n,1.0);
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double i=0;
  int fi=0;
  double ps=0;
  double pa;
  while(fi<Nf){
    i=i+1;
    
    if(gl<(fit_gauss(vec))){
      //std::cout<<ps<<" "<<i<<" "<<fi<<" "<<(fit_gauss(vec))<<" "<<gl<<std::endl;
      vec = do_muller_c(n,1.0);
      ps=pa-i;
      pa=i;
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec,vec_az+vec);
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;
      
      n_p = sqrt(sum((vec)*(vec)));
      vec = vec + vec_az;
      n_pn = sqrt(sum((vec)*(vec)));
      
      if(n_pn < n_p){
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(1-exp(-2*s_r));
        type.push_back("s");
      }
      
      else{
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(0);
        type.push_back("d");
      }
      
    }
  }
  
  timed=time(NULL);
  ti=timed-times;
  
  std::ostringstream oss;
  
  oss<<"walk_glo_N"<<N<<"_n"<<n<<"_r"<<mms<<"_glo"<<gl<<".csv";
  std::string nomarch = oss.str();
  
  std::ofstream arch(nomarch.c_str());
  
  arch<<"c,f,cs,r,tm,s,w,p,p_Ni,type"<<std::endl;
  for(int k=0; k<Nf ;k++){
    arch<<c.at(k)<<","<<f.at(k)<<","<<cs.at(k)<<","<<r.at(k)<<","<<tm.at(k)<<","<<s.at(k)<<","<<w.at(k)<<","<<p.at(k)<<","<<p_Ni.at(k)<<","<<type.at(k)<<std::endl;
  }
  
  arch.close();
  
  std::cout<<"finish glo, N="<<N<<" n="<<n<<" r="<<mms<<" gl="<<gl<<" :"<<ti<<std::endl;
  
}

// [[Rcpp::export]]
DataFrame walk_Nfix(List input){
  
  int times,timed,ti;
  times=time(NULL);
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> cs;
  std::vector<double> r;
  std::vector<double> tm;
  std::vector<double> s;
  std::vector<double> w;
  std::vector<double> p;
  std::vector<double> p_Ni;
  std::vector<std::string> type;
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double tau = input["tau"];
  double sg = input["sg"];
  
  std::cout<<"start, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<std::endl;
  
  NumericVector vec(n); 
  vec = do_muller_c(n,as<double>(runif(1,0,2*mms)));
  NumericVector opt(n);
  opt= optimum_rn_c(opt,sg);
  NumericVector opt_in = opt;
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double chg_step = env_chg_step(tau);
  double i=0;
  int fi=0;
  double ps=chg_step;
  
  while(fi<Nf){
    i=i+1;
    
    if(tau==tau && i>chg_step){
      
      opt = optimum_rn_c(opt_in,sg);
      chg_step = env_chg_step(tau,i);
      ps=chg_step-i;
      std::cout<<ps<<" "<<i<<" "<<fi<<" "<<chg_step<<std::endl;
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec-opt,vec_az+vec-opt);
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;
      
      n_p = sqrt(sum((vec-opt)*(vec-opt)));
      vec = vec + vec_az;
      n_pn = sqrt(sum((vec-opt)*(vec-opt)));
      
      if(n_pn < n_p){
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec-opt));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(1-exp(-2*s_r));
        type.push_back("s");
      }
      else{
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec-opt));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(0);
        type.push_back("d");
      }
      
    }
  }
  
  
  double fs=0;
  double fd=0;
  
  for(int i=1;i<Nf;i++){
    
    fs=fs+p.at(i);
    fd=fd+p_Ni.at(i);
    
  }
  
  
  fs=fs/c.at(Nf-1)*N*N;
  fd=fd/c.at(Nf-1)*N*N;
  
  std::cout<<fs<<" "<<fd<<std::endl;
  
  timed=time(NULL);
  ti=timed-times;
  
  
  
  DataFrame datos = DataFrame::create(Named("c")=c,Named("f")=f,Named("cs")=cs,Named("r")=r,Named("tm")=tm,Named("s")=s,Named("w")=w,Named("p")=p,Named("p_Ni")=p_Ni,Named("fs")=fs,Named("fd")=fd,Named("type")=type);
  
  timed=time(NULL);
  ti=timed-times;
  
  std::cout<<"finish, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<" :"<<ti<<std::endl;
  
  return(datos);
}
// [[Rcpp::export]]
DataFrame walk_Nfix_T(List & input){
  
  int times,timed,ti;
  times=time(NULL);
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> cs;
  std::vector<double> r;
  std::vector<double> tm;
  std::vector<double> s;
  std::vector<double> w;
  std::vector<double> p;
  std::vector<double> p_Ni;
  std::vector<std::string> type;
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double nd = (double) n;
  double Nd = (double) N;
  double gl = pow(1-1/(2*Nd-1),nd/2);
  
  std::cout<<"start glt, N="<<N<<" n="<<n<<" r="<<mms<<" gl="<<gl<<std::endl;
  
  NumericVector vec(n); 
  vec = do_muller_c(n,1.0);
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double i=0;
  int fi=0;
  double ps=0;
  double pa;
  while(fi<Nf){
    i=i+1;
    
    if(gl<(fit_gauss(vec))){
      //std::cout<<ps<<" "<<i<<" "<<fi<<" "<<(fit_gauss(vec))<<" "<<gl<<std::endl;
      vec = do_muller_c(n,1.0);
      ps=pa-i;
      pa=i;
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec,vec_az+vec);
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;
      
      n_p = sqrt(sum((vec)*(vec)));
      vec = vec + vec_az;
      n_pn = sqrt(sum((vec)*(vec)));
      
      if(n_pn < n_p){
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(1-exp(-2*s_r));
        type.push_back("s");
      }
      
      else{
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(0);
        type.push_back("d");
      }
      
    }
  }
  
  DataFrame datos = DataFrame::create(Named("c")=c,Named("f")=f,Named("cs")=cs,Named("r")=r,Named("tm")=tm,Named("s")=s,Named("w")=w,Named("p")=p,Named("p_Ni")=p_Ni,Named("type")=type);
  
  timed=time(NULL);
  ti=timed-times;
  
  std::cout<<"finish glt, N="<<N<<" n="<<n<<" r="<<mms<<" gl="<<gl<<" :"<<ti<<std::endl;
  
  return(datos);
}

// [[Rcpp::export]]
DataFrame walk_Nfix_O(List & input){
  
  int times,timed,ti;
  times=time(NULL);
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> cs;
  std::vector<double> r;
  std::vector<double> tm;
  std::vector<double> s;
  std::vector<double> w;
  std::vector<double> p;
  std::vector<double> p_Ni;
  std::vector<std::string> type;
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double nd = (double) n;
  double Nd = (double) N;
  double gl = 1-nd/(nd+2*Nd);
  std::cout<<"start glo, N="<<N<<" n="<<n<<" r="<<mms<<" gl="<<gl<<std::endl;
  
  NumericVector vec(n); 
  vec = do_muller_c(n,1.0);
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double i=0;
  int fi=0;
  double ps=0;
  double pa;
  while(fi<Nf){
    i=i+1;
    
    if(gl<(fit_gauss(vec))){
      //std::cout<<ps<<" "<<i<<" "<<fi<<" "<<(fit_gauss(vec))<<" "<<gl<<std::endl;
      vec = do_muller_c(n,1.0);
      ps=pa-i;
      pa=i;
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec,vec_az+vec);
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;
      
      n_p = sqrt(sum((vec)*(vec)));
      vec = vec + vec_az;
      n_pn = sqrt(sum((vec)*(vec)));
      
      if(n_pn < n_p){
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(1-exp(-2*s_r));
        type.push_back("s");
      }
      
      else{
        c.push_back(i);
        f.push_back(fi);
        cs.push_back(ps);
        r.push_back(n_pn);
        tm.push_back(sqrt(sum(vec_az*vec_az)));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(0);
        type.push_back("d");
      }
      
    }
  }
  
  DataFrame datos = DataFrame::create(Named("c")=c,Named("f")=f,Named("cs")=cs,Named("r")=r,Named("tm")=tm,Named("s")=s,Named("w")=w,Named("p")=p,Named("p_Ni")=p_Ni,Named("type")=type);
  
  timed=time(NULL);
  ti=timed-times;
  
  std::cout<<"finish glo, N="<<N<<" n="<<n<<" r="<<mms<<" gl="<<gl<<" :"<<ti<<std::endl;
  
  return(datos);
}




DataFrame walk_Nit(List input){
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> r;
  std::vector<double> s;
  std::vector<double> w;
  std::vector<double> p;
  std::vector<double> p_Ni;
  std::vector<std::string> type;
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double tau = input["tau"];
  double sg = input["sg"];
  
  NumericVector vec = do_muller_c(n,as<double>(runif(1,0,2*mms)));
  NumericVector opt(n);
  opt=opt*0;
  NumericVector opt_in = opt;
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double chg_step = env_chg_step(tau);
  double i=0;
  int fi=0;
  
  while(i<Nf){
    i=i+1;
    
    if(tau==tau && i>chg_step){
      
      opt = optimum_rn_c(opt_in,sg);
      chg_step = env_chg_step(tau,i);
      
    }
    
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec-opt,vec_az+vec-opt);
    
    if (prob_fix_K(s_r,N)){
      
      fi=fi+1;
      vec = vec + vec_az;
      n_p = sqrt(sum((vec-opt)*(vec-opt)));
      n_pn = sqrt(sum((vec_az+vec-opt)*(vec_az+vec-opt)));
      
      if(n_pn < n_p){
        c.push_back(i);
        f.push_back(fi);
        r.push_back(sqrt(sum((vec-opt)*(vec-opt))));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec-opt));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(1-exp(-2*s_r));
        type.push_back("s");
      }
      
      else{
        c.push_back(i);
        f.push_back(fi);
        r.push_back(sqrt(sum((vec-opt)*(vec-opt))));
        s.push_back(s_r);
        w.push_back(fit_gauss(vec-opt));
        p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
        p_Ni.push_back(1-exp(-2*s_r));
        type.push_back("d");
        // datos<-rbind(datos, data.frame(c=i, f=fi, r=sqrt(sum((vec-opt)^2)), s=s_r, w=fit_gauss(vec-opt),p=(1-exp(-2*s_r))/(1-exp(-2*input$N*s_r)),p_Ni=1-exp(-2*s_r),type="d"))  
      }
      
    }
  }
  
  std::cout<<"finish, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<std::endl;
  DataFrame datos = DataFrame::create(Named("c")=c,Named("f")=f,Named("r")=r,Named("s")=s,Named("w")=w,Named("p")=p,Named("p_Ni")=p_Ni,Named("type")=type);
  return(datos);
  
}


DataFrame walk_Nfix_test(List input){
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> r;
  std::vector<double> s;
  std::vector<double> w;
  std::vector<double> env;
  
  int N = input["N"];
  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double tau = input["tau"];
  double sg = input["sg"];
  
  NumericVector vec = do_muller_c(n,as<double>(runif(1,0,2*mms)));
  NumericVector opt(n);
  opt=opt*0;
  NumericVector opt_in = opt;
  NumericVector vec_az(n);
  
  double s_r,n_p,n_pn;
  double chg_step = env_chg_step(tau);
  double i=0;
  int fi=0;
  //std::cout<<chg_step<<std::endl;
  
  while(fi<Nf){
    i=i+1;
    
    if(tau==tau && i>chg_step){
      opt = optimum_rn_c(opt_in,sg);
      chg_step = env_chg_step(tau,i);
      //std::cout<<"entre"<<std::endl;
    }
    //std::cout<<i<<" "<<chg_step<<" "<<(tau==tau)<<" "<<(i==chg_step)<<" "<<(tau==tau && i==chg_step)<<std::endl;
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec-opt,vec_az+vec-opt);
    //std::cout<<i<<" "<<fi<<" "<<s_r<<" "<<prob_fix_K(s_r,N)<<std::endl;
    if (prob_fix_K(s_r,N)){
      fi=fi+1;
      vec = vec + vec_az;
      n_p = sqrt(sum((vec-opt)*(vec-opt)));
      n_pn = sqrt(sum((vec_az+vec-opt)*(vec_az+vec-opt)));
      /*if(n_pn < n_p){
      c.push_back(i);
      f.push_back(fi);
      r.push_back(sqrt(sum((vec-opt)*(vec-opt))));
      s.push_back(s_r);
      w.push_back(fit_gauss(vec-opt));
      p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
      p_Ni.push_back(1-exp(-2*s_r));
      type.push_back("s");
    }
      else{
      c.push_back(i);
      f.push_back(fi);
      r.push_back(sqrt(sum((vec-opt)*(vec-opt))));
      s.push_back(s_r);
      w.push_back(fit_gauss(vec-opt));
      p.push_back((1-exp(-2*s_r))/(1-exp(-2*N*s_r)));
      p_Ni.push_back(1-exp(-2*s_r));
      type.push_back("d");
      // datos<-rbind(datos, data.frame(c=i, f=fi, r=sqrt(sum((vec-opt)^2)), s=s_r, w=fit_gauss(vec-opt),p=(1-exp(-2*s_r))/(1-exp(-2*input$N*s_r)),p_Ni=1-exp(-2*s_r),type="d"))  
      }*/
      c.push_back(i);
      env.push_back(chg_step);
      f.push_back(fi);
      r.push_back(sqrt(sum((vec-opt)*(vec-opt))));
      s.push_back(s_r);
      w.push_back(fit_gauss(vec-opt));
      
  }
    /*c.push_back(i);
    env.push_back(chg_step);
    f.push_back(fi);
    r.push_back(sqrt(sum((vec-opt)*(vec-opt))));
    s.push_back(s_r);
    w.push_back(fit_gauss(vec-opt));
    */
    
  }
  
  std::cout<<"finish, N="<<N<<" n="<<n<<" r="<<mms<<" tau="<<tau<<" sg="<<sg<<std::endl;
  DataFrame datos = DataFrame::create(Named("c")=c,Named("env")=env,Named("f")=f,Named("r")=r,Named("s")=s,Named("w")=w);
  return(datos);
}


// [[Rcpp::export]]
DataFrame walk_sphere(List input){

  int n = input["n"];
  double mms = input["mms"];
  int Nf = input["Nf"];
  double d = input["d"];
  double rec = input["f"];
  
  NumericVector vec = do_muller_c(n,d/2);
  double xi = sqrt(n)/d;
  
  std::vector<double> c;
  std::vector<int> f;
  std::vector<double> r;
  std::vector<double> tm;
  std::vector<double> s;
  std::vector<double> w;

  double i=0;
  int fi=0;
  NumericVector vec_az(n);
  double s_r;
  double d_r=d;
  
  
  while(fi<Nf){
    i=i+1;
    vec_az = do_muller_c(n,as<double>(runif(1,0,2*mms)));
    s_r = rel_fit(vec,vec_az+vec);
    
    
    //out<<i<<" "<<sqrt(sum((vec+vec_az)*(vec+vec_az)))<<" "<<d_r<<" "<<fix_sphere(vec+vec_az,d_r)<<std::endl;
    //out<<i<<" "<<s_r<<" "<<prob_fix_ndel(s_r)<<std::endl;
    
    if (prob_fix_ndel(s_r) && fix_sphere(vec_az+vec,d_r)){
      
      fi=fi+1;
      
      vec = vec + vec_az;
      d_r=2*sqrt(sum((vec*vec)));
      c.push_back(i);
      f.push_back(fi);
      r.push_back(d_r/2);
      s.push_back(s_r);
      w.push_back(fit_gauss(vec));
      tm.push_back(sqrt(sum(vec_az*vec_az)));
    }  
    
    if (d_r<(1-rec)*d)
      {
      
      vec=do_muller_c(n,d/2);
      d_r=d;
    }
  }
  std::cout<<"finish, n="<<n<<" d="<<d<<" r="<<mms<<std::endl;
  DataFrame datos = DataFrame::create(Named("c")=c,Named("f")=f,Named("r")=r,Named("tm")=tm,Named("s")=s,Named("w")=w);
  return(datos);
}




/*** R

*/
