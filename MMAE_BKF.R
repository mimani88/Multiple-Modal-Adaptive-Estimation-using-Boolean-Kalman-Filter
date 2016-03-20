MMAE_BKF <- function(net,observ,obs_var,correct_model,threshold,max_time) {
  
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
  # 2-  Imani, M., & Braga-Neto, U. "Optimal Intervention Strategy for Boolean Dynamical Systems Observed through RNA-Seq."
  #  In 2016 IEEE American Control Conference (ACC2016). IEEE.
  # 3- my website: http://people.tamu.edu/~m.imani88/
  
  # Number of state variables 
  d <- ncol(net_connection)
  
  # Create Boolean States
  A <- as.matrix(expand.grid(rep(list(0:1),d)))
  
  graphics.off()
  ########## Seperation of Observation Variables and Parameters 
  obs_model <- observ[1]
  obs_mean <-  as.numeric(observ[2:3])
  obs_variance <- as.numeric(observ[4:length(observ)])
  
  ########## Creating different possible models
  net_connection <- net[1:(nrow(net)-1),]
  net_bias <- net[nrow(net),]
  unk_bias <- which(c(net_bias)=="?")
  unk_connection <- which(c(t(net_connection))=="?")
  
  diff_model <- as.matrix(expand.grid( c(rep(list(-1:1),length(unk_connection)), rep(list(c(-1/2,1/2)),length(unk_bias)) , list((process_noise)), list((obs_variance))  ) ))
  
  # Computing the transition  for all models 
  M <- array(0,c(2^d,2^d,nrow(diff_model)))
  for(i in 1:nrow(diff_model)){
    net2_connection <- c(t(net_connection))
    net2_connection[unk_connection] <- (diff_model[i,1:length(unk_connection)])
    net2_bias <- c(net_bias)
    net2_bias[unk_bias] <- (diff_model[i,unk_bias])
    net2 <- cbind(t(matrix(as.numeric(net2_connection),nrow=ncol(net_connection))), as.numeric(net2_bias))
    pw2 <- rep(diff_model[i,ncol(diff_model)-1],ncol(net_connection))
    #obs_variance2 <- c(obs_variance2,diff_model[i,ncol(diff_model)],ncol(net_connection))
    M[,,i] <- Prediction(0,d,pw2,A,net2)
  }
  
  # Specifying the target model based on user 
  target <- rowmatch(diff_model,correct_model)
  target_connection <- c(t(net_connection))
  target_connection[unk_connection] <- (diff_model[target,1:length(unk_connection)])
  target_bias <- c(net_bias)
  target_bias[unk_bias] <- (diff_model[target,unk_bias])
  net_target <- cbind((t(matrix(as.numeric(target_connection),nrow=ncol(net_connection)))), as.numeric(target_bias))
  pw_target <- rep(correct_model[length(correct_model)-1],d)
  
  # Initial distribution of state
  Pi <- rep(1/2^d,2^d)
  #solve(rbind(rep(1,2^d),diag(2^d)[ c(1:(2^d-1)),]) %*% M[,,target] -  rbind(rep(0,2^d),diag(2^d)[ c(2:(2^d)),]))%*%matrix(c(1,rep(0,2^d-1)),nrow=2^d)
  PVD <- matrix( rep( t( Pi ) , nrow(diff_model) ) , ncol = nrow(diff_model) , byrow = FALSE ) # Posterior probability of different models
  nPVD <- PVD
  # initial state
  X <- A[which(rmultinom(1,1,Pi)>0),]
  #  sX <- matrix(X,nrow=1)
  
  stp <- 0
  #set the initial posterior probability for different models
  pr <- rep(1/nrow(diff_model),nrow(diff_model))
  # The initial probability of correct model
  pr_target <- pr[target]
  k=0
  
  observ_target <- c(obs_model, rep(obs_mean,length(obs_var)), rep(correct_model[length(correct_model)],2*length(obs_var)),obs_var)
  while (stp==0){
    k=k+1
    # Generate the next observation vector
    ndata=simulator(X,pw_target,net_target,0,0,observ_target)
    X <- ndata$X
    
    # Obtain the update matrix given the observed observation (ndata$Y) 
    
    if(length(obs_variance)>1){
      # Compute the posterior probability of Bank of BKF (all models)
      for(m in 1:nrow(diff_model)){
        mean <- rep(obs_mean,length(obs_var))
        var <- rep(diff_model[i,ncol(diff_model)],2*length(obs_var))
        To <- rep(1,2^d)
        To <- rbind(To,(update(ndata$Y,A,c(mean,var),obs_var,obs_model)))
        Tr <- apply(To, 2, prod)
        nPVD[,m] <- diag(Tr)%*%M[,,m]%*%PVD[,m]
      }
    }
    
    if(length(obs_variance)==1){
      # Compute the posterior probability of Bank of BKF (all models)
      mean <- rep(obs_mean,length(obs_var))
      var <- rep(obs_variance,2*length(obs_var))
      To <- rep(1,2^d)
      To <- rbind(To,(update(ndata$Y,A,c(mean,var),obs_var,obs_model)))
      Tr <- apply(To, 2, prod)
      for(m in 1:nrow(diff_model)){
        nPVD[,m] <- diag(Tr)%*%M[,,m]%*%PVD[,m]
      }
    }
    PVD <- sweep(nPVD, 2, colSums(nPVD), FUN="/")
    
    # Computes the unnormalized posterior probability for all models
    Beta <- colSums(nPVD)
    
    # Update the posterior probability of all models
    pr <- (Beta*pr)/sum(Beta*pr)
    
    # Stop the algorithm if the time steps is larger than the specified number defined by user
    if(k>max_time){
      stp <- 1
      cat("Posterior Probability of no model exceeded the defined threshold in the specified time interval.")
      cat("\n")
      if(which.max(pr)!=target){
        message("The correct model is not Inferred in the specified time interval. The inferred model is: ")
        message(diff_model[which.max(pr),])
        message("with Probability")
        message(max(pr))
      }
      if(which.max(pr)==target){
        cat("The correct model is inferred in the specified time interval")
        cat(" with Probability ")
        cat(max(pr))
      }
      plot(pr_target,type="l",lty=1,col=c("blue"),lwd=4.5,xlab="Time",ylab="Posterior Probability of The Correct Model",cex.axis = 1.3,pch=50, cex.lab = 1.3,ylim=c(0,1))
      abline(h=threshold)
      cat("\n")
      cat("Inferred model: ")
      return(diff_model[which.max(pr),])
    }
    
    # Stop the algorithm if posterior probability of any model exceeds the under defined threshold
    if(max(pr)>threshold){
      if(which.max(pr)!=target){
        message("The correct model is not Inferred. The inferred model is: ")
        cat(diff_model[which.max(pr),])
      }
      if(which.max(pr)==target){
        cat(" The correct model is inferred in ")
        cat(k)
        cat(" Time  Steps.")
      }
      stp <- 1
    }
    
    pr_target <- c(pr_target,pr[target])
  }
  
  plot(pr_target,type="l",lty=1,col=c("blue"),lwd=4.5,xlab="Time",ylab="Posterior Probability of The Correct Model",cex.axis = 1.3,pch=50, cex.lab = 1.3,ylim=c(0,1))
  abline(h=threshold)
  cat("\n")
  cat("Inferred model: ")
  return(diff_model[which.max(pr),])
}
