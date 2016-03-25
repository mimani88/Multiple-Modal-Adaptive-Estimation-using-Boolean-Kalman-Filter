Prediction <- function(u,d,pw,A,net){
  
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
  # 2-  Imani, M., & Braga-Neto, U. "State-Feedback Control of Partially-Observed Boolean Dynamical Systems Using RNA-Seq Time Series Data,"
  #  In 2016 IEEE American Control Conference (ACC2016). IEEE.
  # 3- my website: http://people.tamu.edu/~m.imani88/
  
  Md=matrix(0,nrow=2^d,ncol=2^d)
  # Activation-Repression Model
  Xp <- ((t(net%*%t(cbind(A,1))))>0)*1
  for (j in 1:nrow(A)){
    ni <- abs(sweep(Xp, 2, A[j,], `-`))
    Md[j,]=apply(ni*c(pw)+(1-ni)*(1-pw),1,prod)
  }
  return(Md)
}
