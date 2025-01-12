
## portfolio choice (Godley and Lavoie, ch 4) ####################################

# initial parameter values (to be calibrated)
ryear = 0.03
pc_model <- list()
pc_model$alpha1 = 0.5 #Propensity to consume out of income 
pc_model$alpha2 = 0.2933/365 #Propensity to consume out of wealth 
pc_model$lambda0 = 0.4/365 #Autonomous share of bills 
pc_model$lambda1 = 1 #Elasticity of bills demand to interest rate from mp
pc_model$lambda2 = 0.01 #Elasticity of bills demand to yd/v from mp
r = (1+ryear)^(1/365) - 1
pc_model$r = r #Interest rate as policy instrument 

# the names of the variables we solve for in the econ ODE model
pc_model$econvarnames = c('b_h','cons','v','yd')
pc_model$nEconODEs = length(pc_model$econvarnames)

# initial conditions
##!! need tax0 > g0

# using wb data
gdpdata <- wb_data("NY.GDP.MKTP.CN",country = country, start_date = 2018, end_date = 2024)
gdata <- wb_data("NE.CON.GOVT.CN",country = country, start_date = 2018, end_date = 2024)
tdata <- wb_data("GC.TAX.TOTL.CN",country = country, start_date = 2018, end_date = 2024)

tax = subset(tdata,date==2019)$GC.TAX.TOTL.CN/1e12
gment = subset(gdata,date==2019)$NE.CON.GOVT.CN/1e12
gdp = subset(gdpdata,date==2019)$NY.GDP.MKTP.CN/1e12

y0 = gdp/365 # 
g0 = gment/365 # 
tax0 = tax/365 # 
cons0 = y0 - g0
yd0 = cons0
bh0 = (yd0 - y0 + tax0) / r

theta = tax0 / (tax0 + yd0)
pc_model$theta = theta #Tax rate on income

param_to_initial <- function(par){
  
  alpha1 = par[1]
  alpha2 = par[2]
  lambda0 = par[3]
  lambda1 = par[4]
  lambda2 = par[5]
  
  denom = 1 + theta/(1-theta) - alpha1
  denom2 = 1 - (1 - alpha1 + r*365*((1 - alpha1)*(lambda0 + lambda1*r)/alpha2 - lambda2))/denom
  yd_est = y0/(1+denom*denom2)
  v_est = yd_est*(1 - alpha1)/alpha2
  bh_est = 365*(v_est*(lambda0 + lambda1*r) - lambda2*yd_est)
  cons_est = alpha1*yd_est + alpha2*v_est
  # b_s = v_est
  # bcb_est = b_s - bh_est
  # g_est = y0 - yd_est
  print(par)
  print(c(v_est,bh_est,bh0,cons_est ,cons0, yd_est,yd0))
  sum((c(bh0/365 ,cons0,yd0) - c(bh_est/365,cons_est, yd_est))^2) + abs(min(0,v_est-bh_est)/36.5)
}

par = with(pc_model, c(alpha1,alpha2,lambda0,lambda1,lambda2))
param_to_initial(par)
out <- optim(par     = with(pc_model, c(alpha1,alpha2,lambda0,lambda1,lambda2)),  # initial guess
             fn      = param_to_initial,
             method  = "L-BFGS-B",lower = 1e-6, upper=c(1,1/365,Inf,Inf,Inf)
            )
param_to_initial(out$par)
par = out$par

pc_model$alpha1 = par[1]
pc_model$alpha2 = par[2]
pc_model$lambda0 = par[3]
pc_model$lambda1 = par[4]
pc_model$lambda2 = par[5]


pc_model$econ_init = with(pc_model,{
  denom = 1 + theta/(1-theta) - alpha1
  denom2 = 1 - (1 - alpha1 + r*365*((1 - alpha1)*(lambda0 + lambda1*r)/alpha2 - lambda2))/denom
  yd0 = y0/(1+denom*denom2)
  v0 = yd0*(1 - alpha1)/alpha2
  bh0 = 365*(v0*(lambda0 + lambda1*r) - lambda2*yd0)
  cons0 = alpha1*yd0 + alpha2*v0
  b_s = v0
  bcb0 = b_s - bh0
  g0 = y0 - yd0
  c(bh0 ,cons0,v0,yd0)
})

pc_model$g = y0 - pc_model$econ_init[4] # y0 - yd0

pc_model$p_to_scale <- c('alpha1','alpha2')
pc_model$wealth = 'v'

##!! reuse consumption as we do not model workers in this economic model
pc_model$work_value = function(y,econ) {
  y[which(econ$econvarnames=='cons')]
}

pc_model$odes = function(y,econ){

  
  b_h = y[1]
  cons = y[2]
  v = y[3]
  yd = y[4]
  
  alpha1 = econ$alpha1
  alpha2 = econ$alpha2
  
  theta = econ$theta
  lambda0 = econ$lambda0
  lambda1 = econ$lambda1
  lambda2 = econ$lambda2
  r = econ$r
  g = econ$g
  # dot_g = 0
  
  # tax = yd*theta/(1-theta)
  dot_yd = (alpha2*v + g + r*b_h)/(1 + theta/(1-theta) - alpha1) - yd
  # gdp = cons + g
  dot_v = (yd + dot_yd)*(1 - alpha1) - alpha2*v
  dot_cons = alpha1*(yd + dot_yd) + alpha2*v - cons
  dot_b_h = (v + dot_v)*(lambda0 + lambda1*r) - lambda2*(yd + dot_yd) - b_h/365 ##!!
  # dot_b_s = g + dot_g - (yd + dot_yd)*theta/(1-theta) + r*b_h
  # dot_b_cb = dot_b_s - dot_b_h
  
  econ_derivs = c(dot_b_h, #Government bills held by households
                  dot_cons, #Consumption
                  dot_v, #Households wealth
                  dot_yd #Disposable income of households
  )
  
  return(econ_derivs)
  
}

