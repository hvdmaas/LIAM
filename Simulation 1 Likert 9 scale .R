library(qgraph)
resample <- function(x, ...) x[sample.int(length(x), ...)]
hamiltonian <- function(x, n, t, w){
  -sum(t * x) - sum(w * x %*% t(x)/2)
}
glauber_step  <-  function(x, n, t, w, beta){
  i <- sample(1:n, size = 1) # take a random node
  # construct new state with flipped node:
  x_new <- x
  x_new[i] <- resample(spins[-which(abs(spins-x[i])<.000001)],1)
  # update probability
  p <- 1/(1 + exp(beta * (hamiltonian(x_new, n, t, w) - 
                            hamiltonian(x, n, t, w))))  
  if(runif(1) < p) x <- x_new # update state
  return(x)
}

layout(matrix(1:6,2,3,byrow=T))
epsilon <- .005; lambda <- .005 # low values = slow time scale
values = 9
spins=seq(-1,1,length=values)
spins=spins[-(values+1)/2] # delete zero

n <- 10; W <- matrix(rnorm(n^2, .0, .4), n, n)
#W <- matrix(0,n,n)
W <- (W + t(W)) / 2 # make symmetric
diag(W) <- 0
qgraph(W); title('before learning')
thresholds <- rep(0, n)
x=sample(spins,n,TRUE)

for(i in 1:500){
  x <- glauber_step(x, n, thresholds, W, beta = 2)
  # Hebbian learning:
  W <- W + epsilon * (1 - abs(W)) * outer(x, x, "*") - lambda * W 
  diag(W) <- 0
  if(i%%100==0) 
    {
    x_graph=x;W_graph= W 
    W_graph <- x * t(x * W); x_graph=abs(x) # relabeling
    color_func <- colorRamp(c("red", "yellow", "green"))
    normalized_values <- (x_graph + 1) / 2
    rgb_colors <- rgb(color_func(normalized_values) / 255)
    qgraph(W_graph,color=rgb_colors); title(paste("i = ",i))
  }
}
