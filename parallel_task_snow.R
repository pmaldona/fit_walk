require(compiler)
require(snow)
enableJIT(3) 
cl <- makeCluster(8,type="SOCK")



Seq_Nfix <- function(N=c(10,100,1000),n=c(2,10,20),r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=100,tau=1e6,sg=0.85){

  trails <-list()
  name <- character(0)
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
          for(m in tau){
            for(o in sg){

            trails[[length(trails)+1]]<-list(N=i, n=j, mms=k, Nf=l, tau=m, sg=o)

            }
          }
        }
      }
    }
  }

  trails <- parlapply(cl, trails,walk_Nfix)
  #trails <- lapply(trails,walk_Nfix)

  c=0
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
          for(m in tau){
            for(o in sg){

              c=c+1
              trails[[c]]<-list(N=i, n=j, mms=k, Nf=l, tau=m, sg=o,data=trails[[c]])

            }
          }
        }
      }
    }
  }


  return(trails)
}

Seq_Nfix_csv <- function(N=c(10,100,1000),n=c(2,10,20),r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=100,tau=1e6,sg=0.85){
  
  trails <-list()
  name <- character(0)
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
          for(m in tau){
            for(o in sg){
              
              trails[[length(trails)+1]]<-list(N=i, n=j, mms=k, Nf=l, tau=m, sg=o)
              
            }
          }
        }
      }
    }
  }
  
  parLapply(cl, trails,walk_Nfix_v)
}


Seq_Nfix_T <- function(N=c(10,100,1000),n=c(2,10,20),r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=100){
  
  trails <-list()
  name <- character(0)
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
              
          trails[[length(trails)+1]]<-list(N=i, n=j, mms=k, Nf=l)
              
          
        }
      }
    }
  }
  
  trails <- parLapply(cl, trails,walk_Nfix_T)
  #trails <- lapply(trails,walk_Nfix)
  
  c=0
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
            
          c=c+1
          trails[[c]]<-list(N=i, n=j, mms=k, Nf=l, data=trails[[c]])
              
          
        }
      }
    }
  }
  
  
  return(trails)
}

Seq_Nfix_T_csv <- function(N=c(10,100,1000),n=c(2,10,20),r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=100){
  
  trails <-list()
  name <- character(0)
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
            
          trails[[length(trails)+1]]<-list(N=i, n=j, mms=k, Nf=l)
            
          
        }
      }
    }
  }
  
  trails <- parLapply(cl, trails,walk_Nfix_v_T)
  #trails <- lapply(trails,walk_Nfix)
}


Seq_Nfix_O <- function(N=c(10,100,1000),n=c(2,10,20),r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=100){
  
  trails <-list()
  name <- character(0)
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
            
          trails[[length(trails)+1]]<-list(N=i, n=j, mms=k, Nf=l)
            
          
        }
      }
    }
  }
  
  trails <- parLapply(cl, trails,walk_Nfix_O)
  #trails <- lapply(trails,walk_Nfix)
  
  c=0
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
          
          c=c+1
          trails[[c]]<-list(N=i, n=j, mms=k, Nf=l, data=trails[[c]])
              
        }
      }
    }
  }
  
  return(trails)

}

Seq_Nfix_O_csv <- function(N=c(10,100,1000),n=c(2,10,20),r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=100){
  
  trails <-list()
  name <- character(0)
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
         
          trails[[length(trails)+1]]<-list(N=i, n=j, mms=k, Nf=l)
            
        }
      }
    }
  }
  
  parLapply(cl, trails,walk_Nfix_v_O)
  
}


 
Seq_sphere <- function(n=c(2,10,20),d=2.0,r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=100,f=0.9){

  trails <-list()
  for(i in n){
    for(j in d){
      for(k in r){
        for(l in Nf){
          for(m in f){

            trails[[length(trails)+1]]<-list(n=i, d=j, mms=k, Nf=l, f=m)

          }
        }
      }
    }
  }

  trails <- parLapply(cl, trails, walk_sphere)

  c=0
  for(i in n){
    for(j in d){
      for(k in r){
        for(l in Nf){
          for(m in f){

            c=c+1
            trails[[c]]<-list(n=i, d=j, mms=k, Nf=l, f=m,data=trails[[c]])

          }
        }
      }
    }
  }

  return(trails)
}

Seq_Nfix_f <- function(N=c(10,100,1000),n=c(2,10,20),r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=100,tau=1e6,sg=0.85){
  
  trails <-list()
  name <- character(0)
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
          for(m in tau){
            for(o in sg){
              
              trails[[length(trails)+1]]<-list(N=i, n=j, mms=k, Nf=l, tau=m, sg=o)
              
            }
          }
        }
      }
    }
  }
  
  trails <- parLapply(cl, trails,walk_Nfix_f)
  #trails <- lapply(trails,walk_Nfix)
  
  c=0
  for(i in N){
    for(j in n){
      for(k in r){
        for(l in Nf){
          for(m in tau){
            for(o in sg){
              
              c=c+1
              trails[[c]]<-list(N=i, n=j, mms=k, Nf=l, tau=m, sg=o,data=trails[[c]])
              
            }
          }
        }
      }
    }
  }
  
  
  return(trails)
}


# 
# Seq_all <- function(N=c(10,100,1000),n=c(2,10,20),r=c(0.001, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),Nf=25000,tau=NA,sg=0.85){
#   
#   trails <-list()
#   name <- character(0)
#   for(i in N){
#     for(j in n){
#       for(k in r){
#         for(l in Nf){
#           for(m in tau){
#             for(m in sg){
#               
#               trails[[length(trails)+1]]<-list(N=i, n=j, mms=k, Nf=l, tau=m, sg=n)
#               
#             }
#           }
#         }
#       }
#     }
#   }
#   
#   trails <- parLapply(cl, trails,walk_all)
#   
#   c=0
#   for(i in N){
#     for(j in n){
#       for(k in r){
#         for(l in Nf){
#           for(m in tau){
#             for(o in sg){
#               
#               c=c+1
#               trails[[c]]<-list(N=i, n=j, mms=k, Nf=l, tau=m, sg=o,data=trails[[c]])
#               
#             }
#           }
#         }
#       }
#     }
#   }
#   
#   return(trails)
# }

