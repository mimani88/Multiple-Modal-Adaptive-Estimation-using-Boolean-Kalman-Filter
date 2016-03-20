
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

rm(list=ls(all=TRUE))
graphics.off()
source("MMAE_BKS.R")
source("update.R")
source("prediction.R")
source("simulator.R")
source("rowmatch.R")

# The coeficients of  activation/repression model for p53-Mdm2 network with the state Xk = (ATM, p53, Wip1, MDM2).
#a11 = 0, a12 = 0, a13 = -1, a14 = 0   # inputs to gene 1  
#a21 = +1, a22 = 0, a23 = +1, a24 = -1 # inputs to gene 2
#a31 = 0, a32 = +1, a33 = 0, a34 = 0   # inputs to gene 3
#a41 = -1, a42 = +1, a43 = +1, a44 = 0 # inputs to gene 4
net_connection <- rbind(c(0,   0,  -1,  0),
                        c(+1,  0, "?", -1),
                        c(0,  "?", 0, "?"),
                        c(-1, +1, "?",  0))
bias <- c(-1/2,-1/2,-1/2,-1/2) # bias units for genes  
process_noise <- c(0.1,0.15) # different discrete realizations of process noise


obs_mean0 <- c(0)           # mean of x=0
obs_mean1 <- c(1)           # mean of x=1 
obs_variance <- c(0.1)        # var of P(y|x)


correct_model <- c(+1, +1, 0, +1,      # correct elements of net_connection
                   c(),                # correct elements of bias units
                   process_noise[1],   # correct process noise
                   obs_variance[1])    # correct variance of observation noise

observ <- c("Gaussian",obs_mean0,obs_mean1,obs_variance)


obs_variables <- c(1:4)

net <- rbind(net_connection,bias)
threshold <- 0.90  # The threshold for stopping algorithm if one of posterior probabilities exceeds this threshold.  
max_time <- 200    # Maximum number of time steps for stopping the algorithm if the posterior probability of none of models exceeds the threshold.
MMAE_BKS(net,observ,obs_variables,correct_model,threshold,max_time) 
