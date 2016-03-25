simulator <- function(X,pw,net,u_ext,act,observ){
  
  # Reference: Imani, M., & Braga-Neto, U. "Optimal gene regulatory network inference using the Boolean Kalman  filter and 
  #  multiple model adaptive estimation". In 2015 49th Asilomar Conference on Signals, Systems and Computers (pp. 423-427).IEEE.
  
  # This algorithm obtains the optimal network inference of Partially-Observed Boolean Dynamical Systems
  # (e.g. genetic regulatory networs). The algorithm is based on the optimal MMSE state estimator for a Boolean 
  # dynamical system, known as the Boolean Kalman filter (BKF). In the presence of partial knowledge about the network, 
  # a bank of BKFs representing the candidate models is run in parallel in a framework known as Multiple Model Adaptive 
  # Estimation (MMAE). The method is able to infer any dicrete parameters of the process or observation model in an online manner.
  
  # For more information user is referred to the following articles or my personal website:
  # 1-  Imani, M., & Braga-Neto, U. "Optimal state estimation for boolean dynamical systems using a boolean Kalman smoother." 
  #  In 2015 IEEE Global Conference on Signal and Information Processing (GlobalSIP) (pp. 972-976). IEEE.
  # 2-  Imani, M., & Braga-Neto, U. "State-Feedback Control of Partially-Observed Boolean Dynamical Systems Using RNA-Seq Time Series Data "
  #  In 2016 IEEE American Control Conference (ACC2016). IEEE.
  # 3- my website: http://people.tamu.edu/~m.imani88/
  
  if(length(net)==1){  
  Xp <- ((net_model(X,net)+(runif(length(X)) < pw)*1)==1)*1
  }
  if(length(net)>1){  
    Xp <- ((c(((net%*%(c(X,1))))>0)*1+(runif(length(X)) < pw)*1)==1)*1
  }
  if(act<length(u_ext)){
    Xp[u_ext[act]] <- !Xp[u_ext[act]]*1
  }
  if(observ[1]=="direct"){
    return((list(X=Xp, Y=c())))
  }

  if(observ[1]==c("Bernoulli")){
    obs_model <- observ[1]
    obs_par <-  as.numeric(observ[2:(((length(observ)-1)/2)+1)])
    obs_var <-  as.numeric(observ[(((length(observ)-1)/2)+2):length(observ)])
    yy <- ((Xp[obs_var]+(runif(length(obs_par)) < obs_par)*1)==1)*1
    return((list(X=Xp, Y=yy)))
  }
  
  if(observ[1]==c("NB")){
    obs_model <- observ[1]
    obs_par <-  as.numeric(observ[2:(((length(observ)-3)/3)*2+3)])
    obs_var <-  as.numeric(observ[(((length(observ)-3)/3)*2+4):length(observ)])
    yy <- c()
    for(kkl in 1:length(obs_var)){
      ml <- obs_var[kkl]
      lam=obs_par[1]*exp(obs_par[2]+obs_par[2+kkl]*Xp[ml])
      phii <- obs_par[(3+length(obs_var))+kkl-1]
      yy  <- c(yy,rnbinom(1, size=phii, prob=phii/(phii+lam)))
    }
    return((list(X=Xp, Y=yy)))
  }
  
  if(observ[1]==c("Poisson")){
    obs_model <- observ[1]
    obs_par <-  as.numeric(observ[2:(((length(observ)-3)/2)*1+3)])
    obs_var <-  as.numeric(observ[(((length(observ)-3)/2)*1+4):length(observ)])
    yy <- c()
    for(kkl in 1:length(obs_var)){
      ml <- obs_var[kkl]
      lam=obs_par[1]*exp(obs_par[2]+obs_par[2+kkl]*Xp[ml])
       yy  <- c(yy,rpois(1, lam))
    }
    return((list(X=Xp, Y=yy)))
  }

  if(observ[1]==c("Gaussian")){
    obs_par <-  as.numeric(observ[2:(length(observ)-((length(observ)-1)/5))])
    obs_var <-  as.numeric(observ[(length(observ)-((length(observ)-1)/5)+1):length(observ)])
    
    obs_model <- observ[1]
    obs_par <-  as.numeric(observ[2:(((length(observ)-1)/5)*4+1)])
    obs_var <-  as.numeric(observ[(((length(observ)-1)/5)*4+2):length(observ)])
    mu0 <- obs_par[seq(1,length(obs_par)/2,2)]
    mu1 <- obs_par[seq(2,length(obs_par)/2,2)]
#      (length(obs_var)+1):(2*length(obs_var))]
    var0 <- obs_par[seq((length(obs_par)/2+1),length(obs_par),2)]
      #(2*length(obs_var)+1):(3*length(obs_var))]
    var1 <- obs_par[seq((length(obs_par)/2+2),length(obs_par),2)]
      #(3*length(obs_var)+1):(4*length(obs_var))]
    yy <- c()
    for(kkl in 1:length(obs_var)){
      ml <- obs_var[kkl]
      if(Xp[ml]==0){
        yy  <- c(yy,rnorm(1, mu0[kkl], var0[kkl]))
      }
      if(Xp[ml]==1){
        yy  <- c(yy,rnorm(1, mu1[kkl], var1[kkl]))
      }
    }
    return((list(X=Xp, Y=yy)))
  }

  
}
