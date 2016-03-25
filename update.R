update <- function(y,A,obs_par,obs_var,obs_model){
  
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
  
  if(obs_model=="Gaussian"){
    T <- c()
    mu0 <- obs_par[seq(1,length(obs_par)/2,2)]
    mu1 <- obs_par[seq(2,length(obs_par)/2,2)]
    var0 <- obs_par[seq((length(obs_par)/2+1),length(obs_par),2)]
    var1 <- obs_par[seq((length(obs_par)/2+2),length(obs_par),2)]
    for(kkl in 1:length(obs_var)){
      ml <- obs_var[kkl]
      Kc0 <- (1/sqrt(2*pi*var0[kkl]))*exp(-(y[kkl]-mu0[kkl])^2/(2*var0[kkl]))
      Kc1 <- (1/sqrt(2*pi*var1[kkl]))*exp(-(y[kkl]-mu1[kkl])^2/(2*var1[kkl]))
      T <- rbind(T,(A[,ml]==0)*Kc0+ (A[,ml]==1)*Kc1)
    }
    return(T)
  }
}
