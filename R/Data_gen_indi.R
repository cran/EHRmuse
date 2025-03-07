# ##
# # library(MASS)
# 
# K=3 ## Number of Cohorts
# set.seed(100)
# mean_w_p=0
# mean_z_1=0
# mean_z_2=0
# mean_z_3=0
# corr=0.5
# var_z_w_p=matrix(c(1,corr,corr,corr,
#                    corr,1,corr,corr,
#                    corr,corr,1,corr,
#                    corr,corr,corr,1),
#                  nrow=4,ncol=4)
# 
# theta=c(-2,0.35,0.45,0.25) ## Theta_Z vector
# N=5e4 ## Population size
# 
# ### selection models
# dw=1
# dwz1=c(1,0.8,0.6)
# dwz2=c(0.6,0.8,1)
# dwz3=rep(1,3)
# 
# gamma_ext=c(-0.6,1.2,0.4,-0.2,0.5)
# gamma_int_1=c(-1,1.5,0.2,0.8,-0.3)
# gamma_int_2=c(-1,1.25,0.4,0.6)
# gamma_int_3=c(-3,0.8,0.5)
# 
# ## Generation of population level data
# simu_popu<-function(N,mean_w_p,mean_z_1,mean_z_2,mean_z_3,
#                     var_z_w_p,theta,dw){
#   cov<- mvrnorm(n = N, mu = c(mean_w_p,mean_z_1,mean_z_2,mean_z_3), Sigma = var_z_w_p)
#   data <- data.frame(Z1 = cov[, 2], Z2 = cov[, 3], Z3=cov[,4])
#   W_p=cov[,1]
#   # Generate random uniforms
#   #set.seed(5678)
#   U1 <- runif(N)
#   #set.seed(4321)
#   # Generate Disease Status
#   DISEASE <- expit(theta[1] + theta[2] * data$Z1 + theta[3]*data$Z2 +theta[4]*data$Z3)
#   data$D   <- ifelse(DISEASE > U1, 1, 0)
#   # Relate W_p and D
#   data$W_1 <- W_p + dw* data$D + dwz1[1]*data$Z1 + 
#     dwz2[1]*data$Z2 + dwz3[1]*data$Z3 +
#     rnorm(n=N,0,1)
#   
#   data$W_2 <- W_p + dw* data$D + dwz1[2]*data$Z1 + 
#     dwz2[2]*data$Z2 + dwz3[2]*data$Z3 +
#     rnorm(n=N,0,1)
#   
#   data$W_3 <- W_p + dw* data$D + dwz1[3]*data$Z1 + 
#     dwz2[3]*data$Z2 + dwz3[3]*data$Z3 +
#     rnorm(n=N,0,1)
#   
#   data$id=c(1:N)
#   return(data)
# }
# ## Generation of external individual level data
# simu_ext<-function(data,gamma_ext){
#   U2e <- runif(N)
#   # Generate Sampling Status
#   SELECT <-0.75*expit(gamma_ext[1] +  
#                         gamma_ext[2]* data$D + 
#                         gamma_ext[3] * data$Z1 +
#                         gamma_ext[4]* data$Z2 +
#                         gamma_ext[5] * data$Z3)
#   S_e  <- ifelse(SELECT > U2e, T, F)
#   # Observed Data
#   data_e <- data[which(S_e==1),]
#   data_e$Select_Weights = 0.75*expit(gamma_ext[1] +  
#                                        gamma_ext[2]* data_e$D + 
#                                        gamma_ext[3] * data_e$Z1 +
#                                        gamma_ext[4]* data_e$Z2 +
#                                        gamma_ext[5] * data_e$Z3)
#   return(data_e)
# }
# ## Generation of internal data 1
# simu_int_1<-function(data,gamma_int_1){
#   U2i <- runif(N)
#   # Generate Sampling Status
#   SELECT <- expit(cbind(1,data$D,data$W_1,data$Z2,data$Z3)
#                   %*% gamma_int_1)
#   S_i  <- ifelse(SELECT > U2i, T, F)
#   # Observed Data
#   data_i <- data[which(S_i==1),]
#   return(data_i)
# }
# 
# ## Generation of internal data 2
# simu_int_2<-function(data,gamma_int_2){
#   U2i <- runif(N)
#   # Generate Sampling Status
#   SELECT <- expit(cbind(1,data$D,data$W_2,data$Z3)
#                   %*% gamma_int_2)
#   S_i  <- ifelse(SELECT > U2i, T, F)
#   # Observed Data
#   data_i <- data[which(S_i==1),]
#   return(data_i)
# }
# 
# 
# ## Generation of internal data 3
# simu_int_3<-function(data,gamma_int_3){
#   U2i <- runif(N)
#   # Generate Sampling Status
#   SELECT <- expit(cbind(1,data$W_3,data$Z2)
#                   %*% gamma_int_3)
#   S_i  <- ifelse(SELECT > U2i, T, F)
#   # Observed Data
#   data_i <- data[which(S_i==1),]
#   return(data_i)
# }
# 
# data=simu_popu(N,mean_w_p,mean_z_1,mean_z_2,mean_z_3,
#                var_z_w_p,theta,dw)
# 
# extdata=simu_ext(data,gamma_ext)
# 
# 
# intdata1=simu_int_1(data,gamma_int_1)
# intdata2=simu_int_2(data,gamma_int_2)
# intdata3=simu_int_3(data,gamma_int_3)
# 
# ## names of selection variables in each cohort
# select_var_list=list(c("D","W_1","Z2","Z3"),c("D","W_2","Z3"),c("W_3","Z2"))
# 
# ## names of auxiliary variables in each cohort
# aux_var_list=list(c("D","W_1","Z2","Z3"),c("D","W_2","Z3"),c("W_3","Z2"))
# 
# ## list of internal data
# intdata_list=list(intdata1,intdata2,intdata3)
# ## names of Z variables
# Z_names=c("Z1","Z2","Z3")

