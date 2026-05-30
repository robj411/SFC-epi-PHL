
## Model 1, based on SIM (Godley and Lavoie, ch 3) ####################################

model1 <- list()
model1$model_name = 'Model 1'

# initial conditions

# get initial conditions using equations
y0 = ldata$gdp/365 # from annual to daily, in billions
tax0 = ldata$tax/365 # from annual to daily, in billions
theta = tax0 / y0 # from identity, in billions
yd0 = (1-theta)*y0 # from identity, in billions
cons0 = yd0 # from identity, in billions
g0 = y0 - cons0 # gment/365 #  # from identity, in billions
h0 = 19032.66/1e3 # in billions

# store parameters
model1$econ_init = h0
model1$G = g0 # we assume government spending is constant for now
model1$theta = theta #Tax rate on income
model1$alpha1 = 0.6839286/.9 # choices
model1$alpha2 = 0.005339403/365 # choices; goes from annual to daily
model1$alpha0 = with(model1, cons0 - h0*alpha2 - yd0*alpha1) 
print(model1$alpha0)
model1$gdp = y0*365 # for final comparisons
model1$y0 = model1$gdp/365
model1$lf = ldata$NNs[1]/1e6 # labour force in millions
model1$emp0 = model1$lf*ldata$employmentrate/100
model1$lambda = model1$y0/model1$emp0#model1$lf#
model1$cons0 = cons0

# get initial conditions
# model1$econ_init = with(model1,{
#   c( ( cons0*(1-alpha1*(1-theta)) - alpha0 - alpha1*(1-theta)*g0 )/alpha2)
# })


# variable names we need to know later
# which parameters we will scale
# model1$p_to_scale <- c('alpha1','alpha2')
model1$scalar <- 1
# labels for econ variables
model1$econvarlabels = c('Wealth')
# the names of the variables we solve for in the econ ODE model
model1$econvarnames = c('H_h')
# the number of odes
model1$nEconODEs = length(model1$econvarnames)

model1$get_cons_supply <- function(econ){
  # find Y given supply
  Y_s = econ$lf * econ$lambda
  G = econ$G
  cons_s = pmax(0, Y_s - G)
  cons_s
}

model1$get_cons_demand <- function(H_h, econ){
  scalar = econ$scalar
  alpha0 = econ$alpha0
  alpha1 = scalar*econ$alpha1
  alpha2 = scalar*econ$alpha2
  theta = econ$theta
  G = econ$G
  
  denom = 1 - alpha1*(1 - theta)
  cons_d = (alpha0 + alpha1*G*(1-theta) + alpha2*H_h) / denom
  cons_d
}


# function to get consumption from y and econ
model1$get_cons = function(H_h, econ){
  
  # find Y given supply
  cons_s = econ$get_cons_supply(econ)
  
  # compute cons_d (consumption assuming no reduction in labour supply)
  cons_d = econ$get_cons_demand(H_h, econ)
  
  # choose min
  cons = cons_d
  cons[cons_s < cons] <- cons_d[cons_s < cons]
  cons
}

model1$cons_link_fun = function(y, econ){
  H_h = y[1]
  # this is the first time cons is computed, so we save it for later
  econ$cons = econ$get_cons(H_h, econ)
  econ
}

# function to get amounts of economic activity, which will be compared to the counterfactual
model1$econ_to_epi = function(y,econ) {
  # cons has already been computed and saved
  cons_link = econ$cons
  work_link = (cons_link + econ$G)/econ$lambda
  list(cons_link=cons_link, work_link=work_link)
  
}

model1$odes = function(t, y, econ, confirmed, dot_confirmed){
  H_h = y[1] # household wealth
  
  G = econ$G # government spending
  theta = econ$theta # rate of tax
  
  cons = econ$cons # consumption
  Y = cons + G; # gdp
  YD = Y*(1 - theta); # disposable income
  S = YD - cons; # saving
  dot_H_h = S; # rate of wealth accumulation = change in money held by households
  
  econ_derivs = c(dot_H_h) 
  
  list(econ_derivs)
  
}

# function to get consumption from ode matrix output
model1$get_cons_from_timeseries = function(y, econ, epivar, data, integrate=1){
  H_h = y[,2]
  lf_index = 1 + econ$nEconODEs + (data$compindex$C_index-1)*data$nStrata + 1
  lf_confirmed = y[,lf_index]
  
  if (integrate==1){
    econ = epi_to_econ(epivar, lf_confirmed, econ, data)
  }
  cons = econ$get_cons(H_h, econ)
  cons_d = econ$get_cons_demand(H_h, econ)
  cons_s = econ$get_cons_supply(econ)
  list(cons, cons_s, cons_d)
}

# function to get gdp from ode matrix output
model1$get_cons_and_gdp_from_timeseries = function(y,econ,epivar,data,integrate=1){
  G = econ$G
  cons = econ$get_cons_from_timeseries(y,econ,epivar,data,integrate)
  GDP = cons[[1]] + G
  return(list(cons, GDP))
}

# demonstrate behaviour of econ model:
test_model = model1
test_model$alpha1 = model1$alpha1*(1+1/365)
# x <- deSolve::ode(times = 1:10000, y = test_model$econ_init, func = test_model$odes, parms=test_model, method='impAdams_d')
# plot(x[,2])


## Model 2 (Model 1 plus time-varying productivity) ####################################

model2 <- model1
model2$model_name = 'Model 2'

model2$lambda_0 = model1$lambda
model2$lambda_p0 = 0 # intrinsic growth of lambda
model2$lambda_p1 = 0.01 # rate that lambda tracks output

# get initial conditions
model2$econ_init = with(model2,{
  c( ( cons0*(1-alpha1*(1-theta)) - alpha0 - alpha1*(1-theta)*g0 )/alpha2, lambda_0)
})


# variable names we need to know later
# which parameters we will scale
# model2$p_to_scale <- c('alpha1','alpha2')
model2$scalar <- 1
# labels for econ variables
model2$econvarlabels = c('Wealth','Productivity')
# the names of the variables we solve for in the econ ODE model
model2$econvarnames = c('H_h','lambda')
# the number of odes
model2$nEconODEs = length(model2$econvarnames)

model2$get_cons_supply <- function(lambda, econ){
  # find Y given supply
  Y_s = econ$lf * lambda
  G = econ$G
  cons_s = pmax(0, Y_s - G)
  cons_s
}

# function to get consumption from y and econ
model2$get_cons = function(H_h, lambda, econ){
  
  # find Y given supply
  cons_s = econ$get_cons_supply(lambda, econ)
  # Y_s = econ$lf * lambda
  # G = econ$G
  # cons_s = pmax(0, Y_s - G)
  
  # compute cons_d (consumption assuming no reduction in labour supply)
  # cons_d = econ$get_cons_demand(H_h, econ)
  scalar = econ$scalar
  alpha0 = econ$alpha0
  alpha1 = scalar*econ$alpha1
  alpha2 = scalar*econ$alpha2
  theta = econ$theta
  G = econ$G
  
  denom = 1 - alpha1*(1 - theta)
  cons_d = (alpha0 + alpha1*G*(1-theta) + alpha2*H_h) / denom
  
  # choose min
  cons = cons_s
  if(sum(cons_s<cons_d)>1) print(which(cons_s<cons_d))
  cons[cons_d < cons] <- cons_d[cons_d < cons]
  cons
}

# function to get derivative of consumption function
model2$get_gdp_deriv = function(Y, H_h, dot_H_h, confirmed, dot_confirmed, econ, lambda, dot_lf){
  
  scalar = econ$scalar
  alpha0 = econ$alpha0
  alpha1 = econ$alpha1
  alpha2 = econ$alpha2
  alpha1t = scalar * alpha1
  alpha2t = scalar * alpha2
  theta = econ$theta
  G = econ$G
  
  # compute cons_d (consumption assuming no reduction in labour supply)
  denom = 1 - alpha1t*(1 - theta)
  cons_d = (alpha0 + alpha1t*G*(1-theta) + alpha2t*H_h) / denom

  # if consumption is supply determined
  if(econ$cons < cons_d) {
    lf = econ$lf
    # dot_gdp = (lambda * dot_lf + econ$lambda_p0 * lf) / (1 - econ$lambda_p1 * lf)
    dot_gdp = Y * lambda * (dot_lf + econ$lambda_p0 * lf) / 
      (Y - lambda * econ$lambda_p1 * lf)
    # print(c(0, Y, lambda, Y - lambda * econ$lambda_p1 * lf, Y, lambda * econ$lambda_p1 * lf))
    # derivative is not defined for econ$lambda_p1=1
  }else # if consumption is demand determined
  {
    
    dot_link_function = 0
    if(econ$integrate==1)
      dot_link_function = epi_to_econ_deriv(confirmed, dot_confirmed, ldata)
    
    f_term = alpha0 + alpha1t*G*(1-theta) + alpha2t*H_h
    f_prime = dot_link_function * (alpha1*G*(1-theta) + alpha2*H_h) + alpha2t * dot_H_h
    g_term = 1 - alpha1t*(1 - theta)
    g_prime = - dot_link_function * alpha1 * (1-theta)
    
    dot_cons = (g_term*f_prime - f_term*g_prime)/(g_term^2)
    dot_gdp = dot_cons
  }
  return(dot_gdp)
}

model2$cons_link_fun = function(y, econ){
  H_h = y[1]
  lambda = y[2]
  # this is the first time cons is computed, so we save it for later
  econ$cons = econ$get_cons(H_h, lambda, econ)
  econ
}

# function to get amounts of economic activity, which will be compared to the counterfactual
model2$econ_to_epi = function(y,econ) {
  # cons has already been computed and saved
  cons_link = econ$cons
  lambda = y[2]
  # work_link = (cons_link + econ$G)/econ$lambda
  work_link = (cons_link + econ$G)/lambda
  list(cons_link=cons_link, work_link=work_link)
  
}

model2$odes = function(t, y, econ, confirmed, dot_confirmed){
  H_h = y[1] # household wealth
  lambda = y[2] # productivity
  
  total_confirmed = sum(confirmed)
  sum_dot_confirmed = sum(dot_confirmed)
  dot_lf = - dot_confirmed[1] /1e6 # lf is counted in millions
  
  G = econ$G # government spending
  theta = econ$theta # rate of tax
  
  cons = econ$cons # consumption
  Y = cons + G; # gdp
  YD = Y*(1 - theta); # disposable income
  S = YD - cons; # saving
  dot_H_h = S; # rate of wealth accumulation = change in money held by households
  
  # dY/dt = dC/dt because Y=C+G and G is constant
  dot_gdp = econ$get_gdp_deriv(Y, H_h, dot_H_h, total_confirmed, sum_dot_confirmed, econ, lambda, dot_lf)
  
  dot_lambda = lambda * (econ$lambda_p0 + econ$lambda_p1 / Y * dot_gdp)
  # if(t<15)print(c(t,dot_lambda,lambda))
  econ_derivs = c(dot_H_h, dot_lambda) 
  list(econ_derivs)
  
}

# function to get consumption from ode matrix output
model2$get_cons_from_timeseries = function(y, econ, epivar, data, integrate=1){
  H_h = y[,2]
  lambda = y[,3]
  lf_index = 1 + econ$nEconODEs + (data$compindex$C_index-1)*data$nStrata + 1
  lf_confirmed = y[,lf_index]
  if (integrate==1){
    econ = epi_to_econ(epivar, lf_confirmed, econ, data)
  }
  cons = econ$get_cons(H_h, lambda, econ)
  cons_d = econ$get_cons_demand(H_h, econ)
  cons_s = econ$get_cons_supply(lambda, econ)
  list(cons, cons_s, cons_d)
}


