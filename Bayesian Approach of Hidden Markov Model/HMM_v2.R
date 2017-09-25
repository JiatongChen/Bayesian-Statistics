################## setup ####################
set.seed(100)
library(DirichletReg)
library(invgamma)

######### generate simulated data ###########
## HMM generating function
HMM_generate <- function(trans,size,initial_state,sigma,mu){
  hidden_var = initial_state
  obs_var = rnorm(1, mean = mu[initial_state], sd = sigma)
  state = initial_state
  # Make n-1 steps through the markov chain
  for (i in 1:(size-1)){
    p = 0
    u = runif(1, 0, 1)
    #cat("> Dist:", paste(round(c(trans[state,]), 2)), "\n")
    #cat("> Prob:", u, "\n")
    # generate a new state from previous one
    for (j in 1:ncol(trans)){
      p = p + trans[state, j]
      if (p >= u){
        newState = j
        break
      }
    }
    #cat("*", state, "->", newState, "\n")
    hidden_var = c(hidden_var,newState)
    #make a emission from this new state
    obs = rnorm(1, mean = mu[newState], sd = sigma)
    obs_var = c(obs_var,obs)
    state = newState
  }
  return(list(hidden_var=hidden_var,obs_var=obs_var))
}

n = 1000 # size of the Markov chain
trans = matrix(c(1/3,1/3,1/3, # The probability transition matrix
                 0,2/3,1/3,
                 2/3,0,1/3), ncol=3, byrow=TRUE)
state = ceiling(3 * runif(1, 0, 1)) # initial prob for the 3 states are same here
cat("Starting state:", state, "\n")
# normal obs
sigma0 = 0.5
mu0 = c(-2,0,2)
list = HMM_generate(trans,n,state,sigma0,mu0)
hidden_var = list$hidden_var
obs_var = list$obs_var


############ sampling parameters ############
## set up the prior
xi = (min(obs_var)+max(obs_var))/2
R = max(obs_var)-min(obs_var)
alpha = 2
kapa = 1/(R^2)
g = 0.2
h = 10/(R^2)

## parameters
niters = 10000 # number of iterations
M = 3 # state space
mu= mat.or.vec(niters+1,M)
sigsq = mat.or.vec(niters+1,1)
beta = mat.or.vec(niters+1,1)
trans_est = mat.or.vec(niters+1,M^2) # each row represent the matrix ordered by row
rho = mat.or.vec(niters+1,M) # initial probability

# hidden chain
X = mat.or.vec(niters+1,n) 
cluster_size = mat.or.vec(M,1) #n_i
adj_size = mat.or.vec(M,M) #n_{ij}
sum_y = mat.or.vec(M,1) #S_i

# initialize parameters
mu[1,1] = -1
mu[1,2] = 0.5
mu[1,3] = 3
sigsq[1] = 0.4
trans_est[1,] = matrix(c(1/3+0.15,1/3-0.075,1/3-0.075,
                         0+0.075,2/3-0.15,1/3+0.075,
                         2/3-0.15,0+0.075,1/3+0.075), nrow=1)
rho[1,] =  matrix(c(1/3,1/3,1/3), nrow=1)
beta[1] = g/h

# initialize the chain
trans_initial = matrix(trans_est[1,],3,byrow = TRUE)
X[1,] = HMM_generate(trans_initial,n,state,sqrt(sigsq[1]),mu[1,])$hidden_var
#X[1,] = hidden_var #matrix(ceiling(3 * runif(n, 0, 1)),nrow=1)
for (i in 1:M){
  cluster_size[i] = sum(X[1,] == i)
}
for (i in 1:M){
  for(j in 1:M){
    for (k in 2:n){
      if(X[1,k-1] == i & X[1,k] == j){
        adj_size[i,j] = adj_size[i,j] +1
      }
    }
  }
}
for (i in 1:M){
  sum_y[i] = sum(obs_var[X[1,] == i])
}


########## Gibbs sampling iteration ###########
## backward recursion
backward_summation <- function(A,mu,sigsq,y){
  M = dim(A)[1]
  n = length(y)
  B = mat.or.vec(n,M) #backward variables
  B[n,] = 1
  for (m in (n-1):1){
    for (i in 1:M){
      #print(m+1)
      B[m,i]  = sum(A[i,]*dnorm(y[m+1],mean=mu,sd = sqrt(sigsq))*B[m+1,])
    }
    sumBrow = sum(B[m,])
    B[m,] = B[m,]/sumBrow # normalize the row to avoid decaying to 0
  }
  return(B)
}

## main iteration
for (t in 1:niters){
  
  if(t %% 100 == 0){
    cat("iter=",t,"\n")
  }
  
  # step 1: update mu, given hidden data, and sigma^2
  for (i in 1:M){
    mean_i = (sum_y[i] + kapa*xi*sigsq[t]) / (cluster_size[i]+kapa*sigsq[t])
    std_i = sqrt(sigsq[t] / (cluster_size[i]+kapa*sigsq[t]))
    mu[t+1,i] = rnorm(1,mean_i,std_i)
  }
  
  # step 2: update sigma^2, given hidden data, and beta
  shape1 = alpha + 0.5*n
  rate1 = beta + 0.5* sum((matrix(obs_var,nrow=1) - mu[t+1,X[t,]])^2)
  inv = rgamma(1,shape1,rate1)
  sigsq[t+1] = 1/inv
  
  # step 3: update beta, given sigma^2
  beta[t+1] = rgamma(1,g+alpha,h+inv)
  
  # step 4: update A, given sigma^2
  A = mat.or.vec(M,M)
  for (i in 1:M){
    A[i,] = rdirichlet(1, adj_size[i,]+1)
  }
  trans_est[t+1,] = matrix(t(A),nrow=1,byrow = FALSE)
  # matrix(trans_est[t+1,],3,byrow=TRUE)
  
  # step 5: update rho, given X1
  rho[t+1,] = rdirichlet(1, (X[t,1]==1:M)+1)
  
  # step 6: update X, given all parameters
  #X[t+1,1]
  B = backward_summation(A,mu[t+1,],sigsq[t+1],obs_var)
  w = mat.or.vec(n,M)
  for (i in 1:n){
    for (m in 1:M){
      if(i==1){
        w[i,m] =rho[m] * dnorm(obs_var[i],mean=mu[t+1,m],sd=sqrt(sigsq[t+1])) *B[i,m]
      }else{
        #print(i-1)
        w[i,m] =A[ X[t+1,i-1],m] * dnorm(obs_var[i],mean=mu[t+1,m],sd=sqrt(sigsq[t+1])) *B[i,m]
      }
    }
    
    w_sum = sum(w[i,])
    w_sample = runif(1, min = 0, max = w_sum)
    for (m in 1:M){
      if( w_sample<sum(w[i,1:m]) ){
        X[t+1,i] = m
        break
      }
    }
  }
  
  # step 7: update counting parameters
  for (i in 1:M){
    cluster_size[i] = sum(X[t+1,] == i)
  }
  for (i in 1:M){
    for(j in 1:M){
      for (k in 2:n){
        if(X[t+1,k-1] == i & X[t+1,k] == j){
          adj_size[i,j] = adj_size[i,j] +1
        }
      }
    }
  }
  for (i in 1:M){
    sum_y[i] = sum(obs_var[X[t+1,] == i])
  }
  
} 

#################### Plotting #####################
plot(1:(niters+1),mu[1:(niters+1),1],col = 'red',ylim = c(-3,3),xlab ="MCMC Iterations", ylab = "Normal Means", 
     main = "Posterior Samples of Emission Mean", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
points(1:(niters+1),mu[1:(niters+1),2],col = 'blue')
points(1:(niters+1),mu[1:(niters+1),3],col = 'black')


plot(1:(niters+1),sigsq[1:(niters+1)],col = 'blue',ylim = c(0,3.4),xlab ="MCMC Iterations", ylab = "Normal Variance",
     main ="Posterior Samples of Emission Variance", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )

######### calculate posterior statistics ##########
trans_mean = matrix(colMeans(trans_est[300:niters+1,]),M,byrow = TRUE)
sigsq_mean = mean(sigsq[300:niters+1])
mu_mean = colMeans(mu[300:niters+1,])

hidden_samples = X[300:niters+1,]
Xsumcount = mat.or.vec(M,n)
maj_vote_X = mat.or.vec(n,1)
for (j in 1:n){
  for (m in 1:M){
    Xsumcount[m,j] = sum(hidden_samples[,j] == m)
  }
  maj_vote_X[j] =which.max(Xsumcount[,j]) 
}
Accuracy_averaged = sum(maj_vote_X == hidden_var)/n
Accuracy_last_chain = sum(X[niters+1,] == hidden_var)/n