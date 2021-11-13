##############################################
#########     The SpACE Function   ###########
##############################################

#Input: 
#x: The First Observed Vector; y: The Second Observed Vector; x and y should have the same dimention
#K: a integer, corresponding to the degree of approximate polynomials.

#Output: estimated T-score (numeric value)


#the hybrid thresholding estimator without sample splitting
TSK.F <- function(x,y,K=8){
  x = as.vector(x)
  y = as.vector(y)
  n = length(x)
  M = 8*sqrt(log(n))
  K = K
  
  g = rev(polyApprox(abs,-1,1,2*K)$p)
  
  S <- function(r){ 
    s=0
    for(k in 1:(K)){
      s = s + g[2*k+1]*M^(-2*k+1)*hermite(r,2*k)
    }
    return(s)
  }
  
  
  temp1 = rep(0,n)
  if(sum((abs(x)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(x)))>0){
    temp1[(abs(x)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(x))] = S(x[(abs(x)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(x))])
  }
  temp2 = rep(0,n)
  if(sum((abs(y)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(y)))>0){
    temp2[(abs(y)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(y))] = S(y[(abs(y)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(y))])
  }    
  total = sum((temp1 + abs(x) * (abs(x) > 2*sqrt(2*log(n)))) * (temp2 + abs(y) * (abs(y) > 2*sqrt(2*log(n)))))
  
  return(max(c(total,0)))
}

#the hybrid thresholding estimator
TSK <- function(x,y,K=8){
  x = as.vector(x)
  y = as.vector(y)
  n = length(x)
  M = 8*sqrt(log(n))
  K = K
  
  g = rev(polyApprox(abs,-1,1,2*K)$p)
  
  S <- function(r){ 
    s=0
    for(k in 0:(K)){
      s = s + g[2*k+1]*M^(-2*k+1)*hermite(r,2*k)
    }
    return(s)
  }
  
  z = rnorm(n,0,1)
  x1=(x-z)/sqrt(2)
  x2=(x+z)/sqrt(2)
  
  z = rnorm(n,0,1)
  y1=(y-z)/sqrt(2)
  y2=(y+z)/sqrt(2)
  
  
  temp1 = rep(0,n)
  if(sum((abs(x2)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(x2)))>0){
    temp1[(abs(x2)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(x2))] = S(x1[(abs(x2)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(x2))])
  }
  temp2 = rep(0,n)
  if(sum((abs(y2)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(y2)))>0){
    temp2[(abs(y2)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(y2))] = S(y1[(abs(y2)<= 2*sqrt(2*log(n)))&(sqrt(2*log(n))<abs(y2))])
  }    
  total = sum((temp1 + abs(x1) * (abs(x2) > 2*sqrt(2*log(n)))) * (temp2 + abs(y1) * (abs(y2) > 2*sqrt(2*log(n)))))
  
  return(max(c(2*total,0)))
}

#the simple thresholding estimator
Ttilde <- function(x,y){
  n = length(x)
  sum( (abs(x)*(2*sqrt(2*log(n))<abs(x))) * (abs(y)*(2*sqrt(2*log(n))<abs(y))) )
}

#the hybrid non-thresholding estimator
TK <- function(x,y,K=NULL){
  x = as.vector(x)
  y = as.vector(y)
  n = length(x)
  M = 8*sqrt(log(n))
  if(is.null(K)){
    K = round(log(n)/12) 
  }
  
  g = rev(polyApprox(abs,-1,1,2*K)$p)
  
  S <- function(r){ 
    s=0
    for(k in 1:K){
      s = s + g[2*k+1]*M^(-2*k+1)*hermite(r,2*k)
    }
    return(s)
  }
  
  z = rnorm(n,0,1)
  x1=(x-z)/sqrt(2)
  x2=(x+z)/sqrt(2)
  
  z = rnorm(n,0,1)
  y1=(y-z)/sqrt(2)
  y2=(y+z)/sqrt(2)
  
  
  temp1 = rep(0,n)
  if(sum(abs(x)<= 2*sqrt(2*log(n)))>0){
    temp1[abs(x)<= 2*sqrt(2*log(n))] = S(x[abs(x)<= 2*sqrt(2*log(n))])
  }
  temp2 = rep(0,n)
  if(sum((abs(y)<= 2*sqrt(2*log(n))))>0){
    temp2[(abs(y)<= 2*sqrt(2*log(n)))] = S(y[(abs(y)<= 2*sqrt(2*log(n)))])
  }    
  total = sum((temp1 + abs(x1) * (abs(x2) > 2*sqrt(2*log(n)))) * (temp2 + abs(y1) * (abs(y2) > 2*sqrt(2*log(n)))))
  
  
  return(max(c(2*total,0)))
}
