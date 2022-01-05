tm.adjacency <- function(T,RW.order) {
################################################################################
#    create the temporal adjacency matrix for the random walk prior
#    the codes were adapted from the GeoBUGS manual (V1.2 Example 9)
#
#    input:
#===================== 
#       T 			= number of time points
#       RW.order	= 1 or 2 (order of the random walk prior)
#    output: 
#===================== 
#       a list containing three components: adj, num and weights 
#
#       num:	 an array of size T with Entry i showing number of temporal neighbours
#            	 for time point i;
#       adj: 	 an array showing the neighbours of each time point;
#       weights: an array specifying the weights of the temporal 
#                neighbours.
################################################################################
  if (RW.order==1) {
  	####   for random walk order 1
    num <- numeric(T)
    weights <- numeric((T-2)*2+2)
    adj <- numeric((T-2)*2+2)
    for (t in 1:1) {
      weights[t] <- 1
      adj[t] <- t+1
      num[t] <- 1
    }
    for (t in 2:(T-1)) {
      weights[2+(t-2)*2] <- 1
      adj[2+(t-2)*2] <- t-1
      weights[3+(t-2)*2] <- 1
      adj[3+(t-2)*2] <- t+1
      num[t] <- 2
    }
    for (t in T:T) {
      weights[(T-2)*2+2] <- 1
      adj[(T-2)*2+2] <- t-1
      num[T] <- 1
    }
  }
  if (RW.order==2) {
  ####   for random walk order 2
    num <- numeric(T)
    weights <- numeric((T-4)*4+3*2+2*2)
    adj <- numeric((T-4)*4+3*2+2*2)
    for (t in 1:1) {
    	weights[t] <- 2
    	adj[t] <- t+1
    	weights[t+1] <- -1
    	adj[t+1] <- t+2
    	num[t] <- 2
    }
    for (t in 2:2) {
    	weights[t+1] <- 2
    	adj[t+1] <- t-1
    	weights[t+2] <- 4
    	adj[t+2] <- t+1
    	weights[t+3] <- -1
    	adj[t+3] <- t+2
    	num[t] <- 3
    }
    for (t in 3:(T-2)) {
    	weights[6+(t-3)*4] <- -1
    	adj[6+(t-3)*4] <- t-2
    	weights[7+(t-3)*4] <- 4
    	adj[7+(t-3)*4] <- t-1
    	weights[8+(t-3)*4] <- 4
    	adj[8+(t-3)*4] <- t+1
    	weights[9+(t-3)*4] <- -1
    	adj[9+(t-3)*4] <- t+2
    	num[t] <- 4
    }
    for (t in (T-1):(T-1)) {
    	weights[(T-4)*4+6] <- 2
    	adj[(T-4)*4+6] <- t+1
    	weights[(T-4)*4+7] <- 4
    	adj[(T-4)*4+7] <- t-1
    	weights[(T-4)*4+8] <- -1
    	adj[(T-4)*4+8] <- t-2
    	num[t] <- 3
    }
    for (t in T:T) {
    	weights[(T-4)*4+9] <- 2
    	adj[(T-4)*4+9] <- t-1
    	weights[(T-4)*4+10] <- -1
    	adj[(T-4)*4+10] <- t-2
    	num[t] <- 2
    }
  }  
  if (T==2) {
    adj <- c(2,1)
    weights <- c(1,1)
    num <- c(1,1)
  }
  adj.tm <- list()
  adj.tm$adj <- adj
  adj.tm$num <- num
  adj.tm$weights <- weights
  return(adj.tm)
}