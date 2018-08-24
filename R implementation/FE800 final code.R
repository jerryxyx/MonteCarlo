# use sample data------------
S = 1
r = 0.06
dt = 1
K = 1.1
data = read.csv('/Users/wantengxi/Stevens/american option/example.csv',header = FALSE)
data
Ncol = Ncols = ncol(data)
Nrows = Nrow = nrow(data)
#class(data)
cashflow = matrix(0,nrow = Nrow,ncol = Ncol)
cashflow[,Ncol] = ifelse(K > data[,Ncol], K - data[,Ncol], 0)
#j =  Ncols

for(j in Ncols : 2){
  comb_table = cbind(data[,j - 1], cashflow[, j] * exp(-r * dt))
  # stock price at time = j - 1 (X) and value of time j discount back to j - 1 (Y)
  colnames(comb_table) = c('X','Y')
  itm_table = comb_table[comb_table[,'X'] < K,]
  # select only ITM option
  lr = lm(data = as.data.frame(itm_table),formula = Y~X + I(X ^ 2) + I(X ^ 3))
  # regression on constant, X and (X ^ 2) term
  EY  = lr$fitted.values
  # fitted.value = E[Y|X]
  Z = c()
  k = 0
  for(i in 1 : Nrow){
    # select in-the-money option paths
    if(comb_table[i,'X'] < K){
      k = k + 1
      Z[i] = EY[k]
      # assgin EY as Z (for only ITM)
    }
    else{
      Z[i] = comb_table[i,'Y']
      # for OTM, option will not be exercise, so assign continuation value to Z
    }
  }
  C = ifelse(K-comb_table[,'X'] > 0, K-comb_table[,'X'],0)
  # C is in the money exercise value
  comb_2 = cbind(comb_table,Z,C)
  # payoff_temp = cbind(ifelse(K > itm_table[,1],K - itm_table[,1],0),lr$fitted.values)
  cashflow[,j - 1] = ifelse(comb_2[,'C'] < comb_2[,'Z'],comb_2[,'Y'],comb_2[,'C'])
  # cashflow is the greater between C(exercise value) and Z (continuation value)
}
#plot(itm_table)
#points(itm_table[,1],EY,col = 'red')

# use self define parameters under Black-Scholes model------------
S = 36
sigma = 0.2
tau = 1
K = 40
r = 0.06
path = Nrow = 50000
n = 50
Ncol = n + 1
dt = tau / n

# simulation
simu_data = matrix(0,nrow = Nrow,ncol = Ncol)
colnames(simu_data) = paste("step", 0:n, sep="=")
head(simu_data)
# simulate Black Scholes ----------
simu_data[,1] = S
for(i in 2:Ncol){
  epsilon = rnorm(path,0,1)
  simu_data[,i] = simu_data[,i - 1] + r * simu_data[,i - 1] * dt + sigma * simu_data[,i - 1] * sqrt(dt) * epsilon
}
# some other parameter sets start--------
S = 100
K = 100
r = 0.04
tau = 3 / 12
kappa = 1.15
theta = 0.0348
sigma = 0.39
rho = -0.64
v_0 = 0.1866 ^ 2
# some other parameter sets--------
path = Nrow = 50000
n = 50
Ncol = n + 1
dt = tau / n

simu_data = simu_vol = matrix(0,nrow = Nrow,ncol = Ncol)
simu_data[,1] = S
simu_vol[,1] = v_0
v = v_bar = rep(v_0, Nrow)
full = function(v,v_bar,epsilon){
  # full truncation scheme
  v_bar = (v_bar - kappa * dt * (ifelse(v_bar > 0, v_bar, 0) - theta) +
             sigma * sqrt(ifelse(v_bar > 0, v_bar, 0)) * sqrt(dt) * epsilon)
  v = ifelse(v_bar > 0,v_bar,0)
  result = cbind(v, v_bar)
  return(result)
}
# stochasitic volatility regression----
sv_mc = function(data, v_data){
  t_1 = Sys.time()
  # we can add time consumption record
  Ncol = Ncols = ncol(data)
  Nrows = Nrow = nrow(data)
  option_value = matrix(0,nrow = Nrow,ncol = Ncol)
  option_value[,Ncol] = ifelse(K > data[,Ncol], K - data[,Ncol], 0)
  for(j in Ncols : 2){
    # change 3 to 2
    X = data[, j - 1]
    # stock price
    Y = option_value[, j] * exp(-r * dt)
    # discounted cashflow
    vol = v_data[, j - 1]
    # volatility
    comb_table = cbind(X, Y, vol)
    # stock price at time = j - 1 and value of time j discount back to j - 1 and volatility
    
    if(sum(comb_table[comb_table[,'X'] < K,]) < 2){
      EY = Y
    }else{
      itm_table = as.data.frame(comb_table[comb_table[,'X'] < K,])
      # In-the_Money option
      lr = lm(data = itm_table,formula = Y~ X + vol)
      # regression with quadratic term Y~ X + vol
      EY  = lr$fitted.values
      # regression fitted expectation of Y
    }
    Z = c()
    k = 0
    for(i in 1 : Nrow){
      if(comb_table[i,'X'] < K){
        k = k + 1
        Z[i] = EY[k]
      }
      else{
        Z[i] = comb_table[i,'Y']
      }
      
    }
    # Z is the option value combining 1. exercise, the expectation EY; 2 not exercise, the discounted value Y
    C = ifelse(K-comb_table[,'X'] > 0, K-comb_table[,'X'],0)
    # compare whether to exercise
    comb_2 = cbind(comb_table,Z,C)
    option_value[,j - 1] = ifelse(comb_2[,'C'] < comb_2[,'Z'],comb_2[,'Y'],comb_2[,'C'])
    # estimate option value one step backward
    
  }
  #price = mean(option_value[,2]) * exp(-r*dt)
  price = mean(option_value[,1])
  price = ifelse(price > K - S, price, K - S)
  # consider exercise at start
  se = sd(option_value[,2] * exp(-r*dt)) / sqrt(Nrows)
  t_2 = Sys.time()
  time = difftime(t_2,t_1,units = 'secs')
  # output regression_table = itm_table,EY = EY if needed
  result = list(price = price, se = se, time = time)
  return(result)
}
sv_mc_result = sv_mc(data = simu_data, v_data = simu_vol)
sv_mc_result

# check with EFD and crank nicolson method-------
crank.nicolson.method = function(S,K,tao,sigma,r,step,dx, first,
                                 div,type1=c('American','European'),type2=c('Call','Put')){
  payoff=function(S,K,expect,type1=c('American','Europrean'),type2=c('Call','Put')){
    
    temp1=ifelse(type1=='European',0,expect)
    temp2=ifelse(type2=='Call',1,-1)
    result=max(temp2*(S-K),temp1)  
    
    return(result)
  }
  v=r-div-0.5*sigma^2
  dt=tao/step
  pu=-0.25*dt*(sigma^2/dx^2+v/dx)
  pm=1+0.5*dt*sigma^2/dx^2+r*dt/2
  pd=-0.25*dt*(sigma^2/dx^2-v/dx)
  # First we calculate parameters we need
  # not the pm,pu,pd is different from implicit method
  firstRow = firstCol = 1
  nRows = lastRow = 2*step+1
  middleRow = step+1
  nCols = lastCol = step+1
  # Some variables we need to help us understand the position in tree.
  
  V.data = S.data = matrix(0, nrow=nRows, ncol=nCols, dimnames=list(
    paste("NumUps", step:-step, sep="="), paste("T", 0:step, sep="=")))
  S.data[step+1, 1] = S
  # Set the data table and initial stock value   
  
  for (j in 1:(nCols-1)) {
    for(i in (nCols-(j-1)):(nCols+(j-1))) {
      S.data [i-1, j+1] = S.data [i, j]*exp(dx)
      # up case
      S.data [i ,  j+1] = S.data [i, j] 
      # middle case
      S.data [i+1, j+1] = S.data [i, j]*exp(-dx)
      # down case
    }
  }
  # Calculating all stock prices.
  
  for (i in 1:nRows) {
    V.data[i, lastCol] = payoff(S=S.data[i,lastCol],K=K,type1 = 'European',type2 = type2)
  }
  # Calculating the option price at maturity.
  
  lambda.up = ifelse(type2=='Call',1 * (S.data[1, lastCol] - S.data[2,lastCol]),0)
  lambda.low = ifelse(type2=='Call',0,-1 * (S.data[lastRow-1, lastCol] - S.data[lastRow,lastCol]))
  # Boundary condition, same as in implicit method
  
  solve.crank.nicolson.tridiagnoal=function(V.data,pu,pm,pd,lambda.up,lambda.low,colI){
    lastRow = nrow(V.data)
    lastCol = ncol(V.data)
    p.prime = c()
    pm.prime = c()
    # we define p.prime and pm.prime for intermediate steps in the iterations
    pm.prime[lastRow-1] = pm + pd
    p.prime[lastRow-1]  = (-pu*V.data[lastRow-2,lastCol]
                           -(pm-2)*V.data[lastRow-1,lastCol]
                           -pd*V.data[lastRow,lastCol]+pd*lambda.low)
    
    # wo start from the last row (where the boundary took place)
    
    for (j in (lastRow-2):2) {
      pm.prime[j] = pm - pu*pd/pm.prime[j+1]
      p.prime[j] = (-pu*V.data[j-1,colI+1]
                    -(pm-2)*V.data[j,colI+1]
                    -pd*V.data[j+1,colI+1] 
                    - p.prime[j+1]*pd/pm.prime[j+1])
    }
    # solve all of the p.prime and pm.price
    
    V.data[1, colI] = (p.prime[2] + pm.prime[2]*lambda.up)/(pu + pm.prime[2])
    V.data[2, colI] = V.data[1,colI] - lambda.up
    # we get the first two option values
    
    # And then go back  the rest of them
    for(j in 3:(lastRow-1)) {
      V.data[j, colI] =  (p.prime[j] -pu*V.data[j-1, colI])/pm.prime[j]
    }
    V.data[lastRow, colI] = V.data[lastRow-1, colI] - lambda.low
    
    # Out put the V.data(option table)
    
    list(V.data=V.data) 
  }
  
  for(j in (nCols-1):1){
    V.data[, j] = solve.crank.nicolson.tridiagnoal(V.data,pu,pm,pd,lambda.up,lambda.low,colI=j)$V.data[,j]
    if(type1=='American'){
      for(i in 1:nRows){
        V.data[i, j] = payoff(S=S.data[i,lastCol],K=K,type1 = 'American',type2 = type2,
                              expect=V.data[i, j])
      }
      # consider American option can be exercised early
    }
  }
  list(Type = paste(type1,type2), probability=c(pu,pm,pd),
       Price = V.data[step+1,1],
       S.first.steps=S.data[(step+1-first):(step+1+first),1:(1+first)],
       V.first.steps=V.data[(step+1-first):(step+1+first),1:(1+first)]
  )
  # output result including Type, Option price, probability
  # and first steps of Stock and Opton.
}
cn_result = crank.nicolson.method(S = S, K = K, tao = tau, sigma = sigma,
                                  r = r, step = 50, dx = 0.05, first = 3, div = 0, type1 = 'American',type2 = 'Put')
cn_result$Price
sigma * sqrt(3 * dt)

explicit.method = function(S,K,tao,sigma,r,step,dx, first,
                           div,type1=c('American','European'),type2=c('Call','Put')){
  payoff=function(S,K,expect,type1=c('American','Europrean'),type2=c('Call','Put')){
    
    temp1=ifelse(type1=='European',0,expect)
    temp2=ifelse(type2=='Call',1,-1)
    result=max(temp2*(S-K),temp1)  
    
    return(result)
  }
  
  v=r-div-0.5*sigma^2
  dt=tao/step
  pu=0.5*dt*(sigma^2/dx^2+v/dx)
  pm=1-dt*sigma^2/dx^2-r*dt
  pd=0.5*dt*(sigma^2/dx^2-v/dx)
  # First we calculate parameters we need
  
  firstRow = firstCol = 1
  nRows = lastRow = 2*step+1
  middleRow = step+1
  nCols = lastCol = step+1
  # Some variables we need to help us understand the position in tree.
  
  V.data = S.data = matrix(0, nrow=nRows, ncol=nCols, dimnames=list(
    paste("NumUps", step:-step, sep="="), paste("T", 0:step, sep="=")))
  S.data[step+1, 1] = S
  # Set the data table and initial stock value    
  
  for (j in 1:(nCols-1)) {
    for(i in (nCols-(j-1)):(nCols+(j-1))) {
      S.data [i-1, j+1] = S.data [i, j]*exp(dx)
      # up case
      S.data [i ,  j+1] = S.data [i, j] 
      # middle case
      S.data [i+1, j+1] = S.data [i, j]*exp(-dx)
      # down case
    }
  }
  # Calculating all stock prices.
  
  for (i in 1:nRows) {
    V.data[i, lastCol] = payoff(S=S.data[i,lastCol],K=K,type1 = 'European',type2 = type2)
  }
  # Calculating the option price at maturity.
  
  for (j in (nCols-1):1) {
    for(i in (middleRow+(step-1)):(middleRow-(step-1))) {
      V.data[i, j] = (pu*V.data[i-1,j+1] + pm*V.data[i, j+1] + pd*V.data[i+1,j+1])
      
    }
    # Boundary Condition
    stockTerm = ifelse(type2=='Call', (S.data[1,lastCol]-S.data[2,lastCol]), 
                       (S.data[nRows-1,lastCol]-S.data[nRows,lastCol]))
    V.data[firstRow, j] = V.data[firstRow+1,j] + ifelse(type2=='Call', stockTerm, 0)
    V.data[lastRow , j] = V.data[lastRow-1, j] + ifelse(type2=='Call', 0, stockTerm)
    # That is for Call, when stock price is high, dV/dS = 1
    # when stock price is low, dV/dS = 0
    # For put, the when stock price is high, dV/dS = 0
    # when stock price is low, dV/dS = -1
    
    # Then we will add up American option case, deciding whether to exercise
    if(type1=='American') {
      for(i in lastRow:firstRow){
        V.data[i, j] = payoff(S=S.data[i,lastCol],K=K,type1 = 'American',type2 = type2,
                              expect=V.data[i, j])
      }
    }
  }
  ## Step backwards through the trinomial tree
  
  list(Type = paste(type1,type2), probability=c(pu,pm,pd),
       Price = V.data[step+1,1],
       S.first.steps=S.data[(step+1-first):(step+1+first),1:(1+first)],
       V.first.steps=V.data[(step+1-first):(step+1+first),1:(1+first)]
       ## output result including Type, Option price, probability
       ## and first steps of Stock and Opton.
  )
}
explicit.method(S = S, K = K, tao = tau, sigma = sigma,
                r = r, step = 100, dx = 0.05, first = 3, div = 0, type1 = 'American',type2 = 'Put')

# Bakshi parameters------
S = 100
K = 100
r = 0.04
tau = 3 / 12
kappa = 1.15
theta = 0.0348
sigma = 0.39
rho = -0.64
v_0 = 0.1866 ^ 2
path = Nrow = 1000
n = 100
Ncol = n + 1
dt = tau / n

# Simulation with QE scheme -------------

QE_scheme = function(v){
  #browser()
  len =  length(v)
  m = theta + (v - theta) * exp(-kappa * dt)
  s = sqrt(v * sigma ^ 2 * exp(-kappa * dt) / kappa * (1 - exp(-kappa * dt)) + theta * sigma ^ 2 / (2 * kappa) * (1 - exp(-kappa * dt)) ^ 2)
  
  psi = s ^ 2 / m ^ 2
  psi_c = 1.5
  u = runif(len)
  epsilon = qnorm(u)
  # for psi <= psi_c
  b = sqrt(2 / psi - 1 + sqrt(2 / psi * (2 / psi - 1)))
  a = m / (1 + b ^ 2)
  # for psi > psi_c
  p = (psi - 1) / (psi + 1)
  beta = (1 - p) / m
  psi_inv = ifelse((u <= p), 0, 1 / beta * log((1 - p) / (1 - u)))
  v_bar = ifelse(psi <= psi_c, a * (b + epsilon) ^ 2, psi_inv)
  return(v_bar)
}

QE_mc = function(data, v_data){
  t_1 = Sys.time()
  # we can add time consumption record
  Ncol = Ncols = ncol(data)
  Nrows = Nrow = nrow(data)
  option_value = matrix(0,nrow = Nrow,ncol = Ncol)
  option_value[,Ncol] = ifelse(K > data[,Ncol], K - data[,Ncol], 0)
  for(j in Ncols : 2){
    # change 3 to 2
    X = data[, j - 1]
    # stock price
    Y = option_value[, j] * exp(-r * dt)
    # discounted cashflow
    vol = v_data[, j - 1]
    # volatility
    comb_table = cbind(X, Y, vol)
    # stock price at time = j - 1 and value of time j discount back to j - 1 and volatility
    
    if(sum(comb_table[comb_table[,'X'] < K,]) < 2){
      EY = Y
    }else{
      itm_table = as.data.frame(comb_table[comb_table[,'X'] < K,])
      # In-the_Money option
      lr = lm(data = itm_table,formula = Y~ X + vol) 
      # regression with quadratic term Y~ X + vol
      EY  = lr$fitted.values
      # regression fitted expectation of Y
    }
    Z = c()
    k = 0
    for(i in 1 : Nrow){
      if(comb_table[i,'X'] < K){
        k = k + 1
        Z[i] = EY[k]
      }
      else{
        Z[i] = comb_table[i,'Y']
      }
      
    }
    # Z is the option value combining 1. exercise, the expectation EY; 2 not exercise, the discounted value Y
    C = ifelse(K-comb_table[,'X'] > 0, K-comb_table[,'X'],0)
    # compare whether to exercise
    comb_2 = cbind(comb_table,Z,C)
    option_value[,j - 1] = ifelse(comb_2[,'C'] < comb_2[,'Z'],comb_2[,'Y'],comb_2[,'C'])
    # estimate option value one step backward
    
  }
  #price = mean(option_value[,2]) * exp(-r*dt)
  price = mean(option_value[,1])
  price = ifelse(price > K - S, price, K - S)
  # consider exercise at start
  se = sd(option_value[,2] * exp(-r*dt)) / sqrt(Nrows)
  t_2 = Sys.time()
  time = difftime(t_2,t_1,units = 'secs')
  # output regression_table = itm_table,EY = EY if needed
  result = list(price = price, se = se,time = time)
  return(result)
}

QE_lsm = function(path){
  #browser()
  Nrow = path
  t_1 = Sys.time()
  simu_data = simu_vol = matrix(0,nrow = Nrow,ncol = Ncol)
  colnames(simu_data) = paste("step", 0:n, sep="=")
  simu_data[,1] = S
  simu_vol[,1] = v_0
  
  
  # Heston simulation with Euler Discretization
  for(i in 2 : Ncol){
    
    bm_s = rnorm(Nrow,0,1)
    bm_s2 = rnorm(Nrow,0,1)
    bm_v = rho * bm_s + sqrt(1 - rho ^ 2) * bm_s2
    drift = (r - 0.5 * simu_vol[,i - 1]) * dt
    diffusion = sqrt(simu_vol[,i - 1]) * sqrt(dt) * bm_s
    simu_data[,i] = exp(log(simu_data[,i - 1]) + drift + diffusion)
    
    simu_vol[,i] = QE_scheme(simu_vol[,i - 1])
  }
  
  
  lsmc_result = QE_mc(simu_data, simu_vol)
  price = lsmc_result[[1]]
  se = lsmc_result[[2]]
  t_2 = Sys.time()
  time = difftime(t_2,t_1,units = 'secs')
  result = list(price = price, se = se, time = time)
  return(result)
}

QE_lsm_result = QE_lsm(5000)
QE_lsm_result

# weighted Laguerre polynomials--------
laguerre_0 = function(x){
  y = exp(-x/2)
  return(y)
}
laguerre_1 = function(x){
  y = exp(-x/2) * (1 - x)
  return(y)
}
laguerre_2 = function(x){
  y = exp(-x/2) * (1 - 2 * x + x ^ 2 / 2)
  return(y)
}
laguerre_3 = function(x){
  y = exp(-x/2) * (2 * x + 3 / 2 * x ^ 2 - 1 / 6 * x ^ 3)
  # - 3*x + 1.5*x^2 - 0.1666667*x^3 
  return(y)
}
laguerre_4 = function(x){
  y = exp(-x/2) * (1 - 4 * x + 3 * x ^ 2 - 1 / 3 * x ^ 3 + 1 / 24 * x ^ 4)
  # - 3*x + 1.5*x^2 - 0.1666667*x^3 
  return(y)
}
laguerre.polynomials(3)

lsmc_laguerre = function(data){
  #browser()
  Ncol = Ncols = ncol(data)
  Nrows = Nrow = nrow(data)
  option_value = matrix(0,nrow = Nrow,ncol = Ncol)
  option_value[,Ncol] = ifelse(K > data[,Ncol], K - data[,Ncol], 0)
  for(j in Ncols : 2){
    # change 3 to 2
    X = data[,j - 1]
    Y = option_value[, j] * exp(-r * dt)
    
    
    # insert laguerre polynomial
    L_0 = laguerre_0(X)
    L_1 = laguerre_1(X)
    L_2 = laguerre_2(X)
    L_3 = laguerre_3(X)
    L_4 = laguerre_4(X)
    comb_table = cbind(X, L_0, L_1, L_2, L_3, L_4, Y)
    # stock price at time = j - 1 and value of time j discount back to j - 1
    
    if(nrow(comb_table[comb_table[,'X'] < K,]) < 2){
      EY = Y
    }else{
      itm_table = as.data.frame(comb_table[comb_table[,'X'] < K,])
      lr = lm(data = itm_table,formula = Y ~ L_0 + L_1 + L_2 + L_3  + L_4)
      # regression with quadratic term Y ~ L_0 + L_1 + L_2 + L_3  + L_4
      EY  = lr$fitted.values
      # regression fitted expectation of Y
    }
    Z = c()
    k = 0
    for(i in 1 : Nrow){
      if(comb_table[i,'X'] < K){
        k = k + 1
        Z[i] = EY[k]
      }
      else{
        Z[i] = comb_table[i,'Y']
      }
    }
    # Z is the option value combining 1. exercise, the expectation EY; 2 not exercise, the discounted value Y
    C = ifelse(K-comb_table[,'X'] > 0, K-comb_table[,'X'],0)
    # compare whether to exercise
    comb_2 = cbind(comb_table,Z,C)
    option_value[,j - 1] = ifelse(comb_2[,'C'] < comb_2[,'Z'],comb_2[,'Y'],comb_2[,'C'])
    # estimate option value one step backward
    
  }
  #price = mean(option_value[,2]) * exp(-r*dt)
  price = mean(option_value[,1])
  price = ifelse(price > K - S, price, K - S)
  # consider exercise at start
  se = sd(option_value[,2] * exp(-r*dt)) / sqrt(Nrows)
  # output regression_table = itm_table,EY = EY if needed
  result = c(price, se)
  return(result)
}


lsmc_laguerre_result = lsmc_laguerre(simu_data)
lsmc_laguerre_result

# path as input---------
BS_lsm_laguerre = function(path){
  Nrow = path
  t_1 = Sys.time()
  
  
  simu_data = matrix(0,nrow = Nrow,ncol = Ncol)
  colnames(simu_data) = paste("step", 0:n, sep="=")
  # simulate Black Scholes ----------
  simu_data[,1] = S
  for(i in 2:Ncol){
    epsilon = rnorm(path,0,1)
    simu_data[,i] = simu_data[,i - 1] + r * simu_data[,i - 1] * dt + sigma * simu_data[,i - 1] * sqrt(dt) * epsilon
  }
  
  lsmc_laguerre_result = lsmc_laguerre(simu_data)
  price = lsmc_laguerre_result[1]
  se = lsmc_laguerre_result[2]
  t_2 = Sys.time()
  time = difftime(t_2,t_1,units = 'secs')
  result = list(price = price, se = se, time = time)
  return(result)
}
BS_lsm_laguerre_result = BS_lsm_laguerre(2000)
BS_lsm_laguerre_result


