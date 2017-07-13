Rcpp::sourceCpp('Walk_cpp_ftest.cpp')
source('parallel_task_ftest.R')

N <- c(500,1000,2000,8000)
n <- c(2,5,10,20,25,30,40,50)
r <- c(1)
Nm <- 1e7
tau <-  seq(-6,-4,by = 0.25)
tau <- 1/10^(tau)
sg <- c(0.85)

force_r_8000 <- Seq_Nfix_f(N=N,n=n,r=r,tau=tau,sg=sg,Nf=Nm)


