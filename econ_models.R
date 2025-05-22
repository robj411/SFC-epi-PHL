
## Model 1, based on SIM (Godley and Lavoie, ch 3) ####################################

model1 <- list()
model1$model_name = 'model1'

# initial conditions
# using wb data
gdpdata <- wb_data("NY.GDP.MKTP.CN",country = country, start_date = 2018, end_date = 2024)
# gdata <- wb_data("NE.CON.GOVT.CN",country = country, start_date = 2018, end_date = 2024)
tdata <- wb_data("GC.TAX.TOTL.CN",country = country, start_date = 2018, end_date = 2024)

ref_year = 2023
tax = c(subset(tdata,date==ref_year)$GC.TAX.TOTL.CN/1e12)
# gment = c(subset(gdata,date==ref_year)$NE.CON.GOVT.CN/1e12)
gdp = c(subset(gdpdata,date==ref_year)$NY.GDP.MKTP.CN/1e12)

# get initial conditions using equations
y0 = gdp/365 # from annual to daily
tax0 = tax/365 # from annual to daily
theta = tax0 / y0 # from identity
yd0 = (1-theta)*y0 # from identity
cons0 = yd0 # from identity
g0 = y0 - cons0 # gment/365 #  # from identity

# store parameters
model1$G = g0 # we assume government spending is constant for now
model1$theta = theta #Tax rate on income
model1$alpha1 = 0.9 # choices
model1$alpha2 = .2933/365 # choices; goes from annual to daily
model1$alpha0 = cons0*(1-model1$alpha1)/3 # choices
model1$gdp = y0*365 # for final comparisons

# get initial conditions
model1$econ_init = with(model1,{
  ( cons0*(1-alpha1*(1-theta)) - alpha0 - alpha1*(1-theta)*g0 )/alpha2
})


# variable names we need to know later
# which parameters we will scale
model1$p_to_scale <- c('alpha1','alpha2')
# which variable records wealth
model1$wealth = 'H_h'
# the names of the variables we solve for in the econ ODE model
model1$econvarnames = c('H_h')
# the number of odes
model1$nEconODEs = length(model1$econvarnames)


# function to get consumption from y and econ
model1$get_cons = function(y,econ){
  H_h = y
  
  alpha0 = econ$alpha0
  alpha1 = econ$alpha1
  alpha2 = econ$alpha2
  theta = econ$theta
  G = econ$G
  denom = 1 - alpha1*(1 - theta)
  cons = (alpha0 + alpha1*G*(1-theta) + alpha2*H_h) / denom
  cons
}

# function to get consumption from ode matrix output
model1$get_cons_from_out = function(y,econ,H,data,integrate=1){
  H_h = y[,2]
  if (integrate==1){
    econ = epi_to_econ(H,econ,ref_val=data$ref_val,baseline=data$baseline)
  }
  cons = econ$get_cons(H_h,econ)
  cons
}

# function to get gdp from ode matrix output
model1$get_gdp_from_out = function(y,econ,H,data,integrate=1){
  G = econ$G
  cons = econ$get_cons_from_out(y,econ,H,data,integrate)
  cons + G
}

##!! reuse consumption as we do not model workers in this economic model
model1$econ_to_epi = function(y,econ) {
  H_h = y[1]
  cons_link = econ$get_cons(H_h,econ)
  list(cons_link=cons_link, work_link=cons_link)
}

model1$odes = function(t,y,econ){
  H_h = y[1] # household wealth
  G = econ$G # government spending
  theta = econ$theta # rate of tax
  
  C = econ$get_cons(H_h,econ) # consumption
  Y = C + G; # gdp
  YD = Y*(1 - theta); # disposable income
  S = YD - C; # saving
  dot_H_h = S; # rate of wealth accumulation
  
  econ_derivs = dot_H_h # change in money held by households
  
  list(econ_derivs)
  
}

# demonstrate behaviour of econ model:
test_model = model1
test_model$alpha1 = model1$alpha1*(1+1/365)
x <- deSolve::ode(times = 1:10000, y = test_model$econ_init, func = test_model$odes, parms=test_model, method='impAdams_d')
plot(x[,2])

## model 2: online consumption ###############################

model2 = model1
model2$model_name = 'model2'
model2$q4 = 0.5;
model2$alpha1online = model2$q4 * model1$alpha1
model2$alpha2online = model2$q4 * model1$alpha2
model2$alpha1offline = (1 - model2$q4) * model1$alpha1
model2$alpha2offline = (1 - model2$q4) * model1$alpha2

model2$p_to_scale <- c('alpha1offline','alpha2offline')

## re-scale alpha1 and alpha2 within ODE model; hospitalisation affects the propensity to consume _in person_,
# but the economic impact is diminished because of online consumption. Say, half the transactions that would
# have been lost move online.

model2$get_cons = function(y,econ) {
  H_h = y
  alpha0 = econ$alpha0
  alpha1online = econ$alpha1online
  alpha2online = econ$alpha2online
  alpha1offline = econ$alpha1offline
  alpha2offline = econ$alpha2offline
  alpha1 = alpha1online + alpha1offline
  alpha2 = alpha2online + alpha2offline
  theta = econ$theta
  G = econ$G
  denom = 1 - alpha1*(1 - theta)
  cons = (alpha0 + alpha1*G*(1-theta) + alpha2*H_h) / denom
  cons
}


model2$econ_to_epi = function(y,econ) {
  H_h = y[1]
  alpha0 = econ$alpha0
  alpha1offline = econ$alpha1offline
  alpha2offline = econ$alpha2offline
  alpha1 = econ$alpha1online + alpha1offline
  alpha2 = econ$alpha2online + alpha2offline
  # in-person component of consumption
  # C = econ$get_cons(H_h,econ)
  # cons_link = (alpha0 + alpha1offline + alpha2offline) / (alpha0 + alpha1 + alpha2) * C
  G = econ$G
  denom = 1 - alpha1*(1 - theta)
  cons_link = (alpha0 + alpha1offline*G*(1-theta) + alpha2offline*H_h) / denom
  ##!! reuse consumption for work link as we do not model workers in this economic model
  list(cons_link=cons_link, work_link=cons_link)
}



## model 3: labour supply ###############################

model3 = model2
model3$model_name = 'model3'
model3$prop_to_work = 1
model3$y0 = model3$gdp/365
model3$lf = ldata$NNs[1]
model3$wage = model3$y0/model3$lf 

model3$p_to_scale <- c('alpha1offline','alpha2offline','prop_to_work')

model3$get_cons = function(y,econ) {
  H_h = y
  alpha0 = econ$alpha0
  alpha1online = econ$alpha1online
  alpha2online = econ$alpha2online
  alpha1offline = econ$alpha1offline
  alpha2offline = econ$alpha2offline
  alpha1 = alpha1online + alpha1offline
  alpha2 = alpha2online + alpha2offline
  theta = econ$theta
  G = econ$G
  prop_to_work = econ$prop_to_work
  
  # find out the labour supply
  wb_s = econ$lf * prop_to_work * econ$wage
  cons_s = pmax(0, wb_s - G)
  
  # compute cons_d assuming no reduction in labour supply
  denom = 1 - alpha1*(1 - theta)
  cons_d = (alpha0 + alpha1*G*(1-theta) + alpha2*H_h) / denom
  
  # choose min
  C = cons_s
  # if(is.na(cons)) browser()
  # print(format(c(cons,cons_d,cons_s),digits=20))
  C[cons_d <= C + 1e-10] <- cons_d[cons_d <= C + 1e-10]
  # if(cons_d <= cons_s + 1e-10){
  #   C = cons_d
  # }
  C
}

model3$econ_to_epi = function(y,econ) {
  H_h = y[1]
  alpha0 = econ$alpha0
  alpha1offline = econ$alpha1offline
  alpha2offline = econ$alpha2offline
  alpha1 = econ$alpha1online + alpha1offline
  alpha2 = econ$alpha2online + alpha2offline
  G = econ$G
  denom = 1 - alpha1*(1 - theta)
  ##!! unless in-person consumption is also limited by labour supply...
  cons_link = (alpha0 + alpha1offline*G*(1-theta) + alpha2offline*H_h) / denom
  # C = econ$get_cons(H_h,econ)
  # cons_link = (alpha1offline + alpha2offline) / (alpha1 + alpha2) * C
  work_link = econ$lf * econ$prop_to_work
  list(cons_link=cons_link, work_link=work_link)
  
}

## model 4 imports and exports ##############################################

# NE.IMP.GNFS.CN Imports of goods and services (current LCU)
# NE.EXP.GNFS.CN Exports of goods and services (current LCU)
idata <- wb_data("NE.IMP.GNFS.CN",country = country, start_date = 2018, end_date = 2024)
edata <- wb_data("NE.EXP.GNFS.CN",country = country, start_date = 2018, end_date = 2024)
imp = c(subset(idata,date==2019)$NE.IMP.GNFS.CN/1e12)
ex = c(subset(edata,date==2019)$NE.EXP.GNFS.CN/1e12)
bop = ex - imp


## portfolio choice (Godley and Lavoie, ch 4) ####################################

# initial parameter values (to be calibrated)
ryear = 0.03
pc_model <- model1
pc_model$model_name = 'PC'
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

tax = c(subset(tdata,date==2019)$GC.TAX.TOTL.CN/1e12)
gment = c(subset(gdata,date==2019)$NE.CON.GOVT.CN/1e12)
gdp = c(subset(gdpdata,date==2019)$NY.GDP.MKTP.CN/1e12)

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
  # print(par)
  # print(c(v_est,bh_est,bh0,cons_est ,cons0, yd_est,yd0))
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

pc_model$gdp = y0 * 365
pc_model$G = y0 - pc_model$econ_init[4] # y0 - yd0

pc_model$p_to_scale <- c('alpha1','alpha2')
pc_model$wealth = 'v'

##!! reuse consumption as we do not model workers in this economic model
pc_model$econ_to_epi = function(y,econ) {
  cons_link = y[which(econ$econvarnames=='cons')]
  list(cons_link=cons_link, work_link=cons_link)
}

pc_model$odes = function(t,y,econ){
  
  
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
  G = econ$G
  # dot_g = 0
  
  # tax = yd*theta/(1-theta)
  dot_yd = (alpha2*v + G + r*b_h)/(1 + theta/(1-theta) - alpha1) - yd
  # gdp = cons + G
  dot_v = (yd + dot_yd)*(1 - alpha1) - alpha2*v
  dot_cons = alpha1*(yd + dot_yd) + alpha2*v - cons
  dot_b_h = (v + dot_v)*(lambda0 + lambda1*r) - lambda2*(yd + dot_yd) - b_h/365 ##!!
  # dot_b_s = G + dot_g - (yd + dot_yd)*theta/(1-theta) + r*b_h
  # dot_b_cb = dot_b_s - dot_b_h
  
  econ_derivs = c(dot_b_h, #Government bills held by households
                  dot_cons, #Consumption
                  dot_v, #Households wealth
                  dot_yd #Disposable income of households
  )
  
  return(econ_derivs)
  
}



## model pc 2: online consumption ###############################

modelpc2 = pc_model
modelpc2$model_name = 'modelpc2'
modelpc2$q4 = 0.5;
modelpc2$alpha1online = modelpc2$q4 * pc_model$alpha1
modelpc2$alpha2online = modelpc2$q4 * pc_model$alpha2
modelpc2$alpha1offline = (1 - modelpc2$q4) * pc_model$alpha1
modelpc2$alpha2offline = (1 - modelpc2$q4) * pc_model$alpha2

modelpc2$econ_init = pc_model$econ_init[c(1,2,2,2,3,4)]
modelpc2$econ_init[2] = modelpc2$q4 * modelpc2$econ_init[2]
modelpc2$econ_init[3] = (1 - modelpc2$q4) * modelpc2$econ_init[3]

modelpc2$p_to_scale <- c('alpha1offline','alpha2offline')

modelpc2$econvarnames = c('b_h','cons_on','cons_off','cons','v','yd')
modelpc2$nEconODEs = length(modelpc2$econvarnames)

## re-scale alpha1 and alpha2 within ODE model; hospitalisation affects the propensity to consume _in person_,
# but the economic impact is diminished because of online consumption. Say, half the transactions that would
# have been lost move online.

modelpc2$odes = function(t,y,econ){
  
  
  b_h = y[1]
  cons_on = y[2]
  cons_off = y[3]
  cons = y[4]
  v = y[5]
  yd = y[6]
  
  q4 = econ$q4
  
  alpha1online = econ$alpha1online
  alpha2online = econ$alpha2online
  alpha1offline = econ$alpha1offline
  alpha2offline = econ$alpha2offline
  alpha1 = alpha1online + alpha1offline
  alpha2 = alpha2online + alpha2offline
  
  theta = econ$theta
  lambda0 = econ$lambda0
  lambda1 = econ$lambda1
  lambda2 = econ$lambda2
  r = econ$r
  G = econ$G
  # dot_g = 0
  
  # tax = yd*theta/(1-theta)
  dot_yd = (alpha2*v + G + r*b_h)/(1 + theta/(1-theta) - alpha1) - yd
  # gdp = cons + G
  dot_v = (yd + dot_yd)*(1 - alpha1) - alpha2*v
  dot_cons_on = alpha1online*(yd + dot_yd) + alpha2online*v - cons_on
  dot_cons_off = alpha1offline*(yd + dot_yd) + alpha2offline*v - cons_off
  dot_cons = dot_cons_on + dot_cons_off # alpha1*(yd + dot_yd) + alpha2*v - cons # 
  dot_b_h = (v + dot_v)*(lambda0 + lambda1*r) - lambda2*(yd + dot_yd) - b_h/365 ##!!
  # dot_b_s = G + dot_g - (yd + dot_yd)*theta/(1-theta) + r*b_h
  # dot_b_cb = dot_b_s - dot_b_h
  
  econ_derivs = c(dot_b_h, #Government bills held by households
                  dot_cons_on, #Consumption online
                  dot_cons_off, #Consumption offline
                  dot_cons, # consumption
                  dot_v, #Households wealth
                  dot_yd #Disposable income of households
  )
  
  return(econ_derivs)
  
}

modelpc2$econ_to_epi = function(y,econ) {
  ##!! reuse consumption as we do not model workers in this economic model
  cons_link = y[which(econ$econvarnames=='cons_off')]
  list(cons_link=cons_link, work_link=cons_link)
}


## model pc 3: labour supply ###############################

modelpc3 = modelpc2
modelpc3$model_name = 'modelpc3'
modelpc3$prop_to_work = 1
modelpc3$y0 = modelpc3$G + pc_model$econ_init[2]
modelpc3$lf = ldata$NNs[1]
modelpc3$wage = y0/modelpc3$lf 

modelpc3$p_to_scale <- c('alpha1offline','alpha2offline','prop_to_work')

modelpc3$econ_to_epi = function(y,econ) {
  lf = econ$lf
  cons_link = y[which(econ$econvarnames=='cons_off')]
  work_link = lf * econ$prop_to_work
  list(cons_link=cons_link, work_link=work_link)
  
}

modelpc3$get_cons_from_out = function(y,econ,H,data,integrate=1){
  y[,which(econ$econvarnames=='cons')+1]
}

## get cons_s and cons_d, assuming there was no hesitation to work.

modelpc3$odes = function(t,y,econ){
  
  
  b_h = y[1]
  cons_on = y[2]
  cons_off = y[3]
  cons = y[4]
  v = y[5]
  yd = y[6]
  
  alpha1online = econ$alpha1online
  alpha2online = econ$alpha2online
  alpha1offline = econ$alpha1offline
  alpha2offline = econ$alpha2offline
  alpha1 = alpha1online + alpha1offline
  alpha2 = alpha2online + alpha2offline
  prop_to_work = econ$prop_to_work
  
  theta = econ$theta
  lambda0 = econ$lambda0
  lambda1 = econ$lambda1
  lambda2 = econ$lambda2
  r = econ$r
  G = econ$G
  # dot_g = 0
  
  # find out the labour supply
  wb_s = econ$lf * prop_to_work * econ$wage
  cons_s = max(0, wb_s - G)
  dot_cons_s = cons_s - cons
  
  # compute cons_d assuming no reduction in labour supply
  
  dot_yd = (alpha2*v + G + r*b_h)/(1 + theta/(1-theta) - alpha1) - yd
  dot_cons_on = alpha1online*(yd + dot_yd) + alpha2online*v - cons_on
  dot_cons_off = alpha1offline*(yd + dot_yd) + alpha2offline*v - cons_off
  dot_cons_d = dot_cons_on + dot_cons_off # alpha1*(yd + dot_yd) + alpha2*v - cons # 
  cons_d = dot_cons_d + cons
  
  if(cons_off<0) browser()
  # choose min
  # if(is.na(cons)) browser()
  # print(format(c(cons,cons_d,cons_s),digits=20))
  if(cons_d <= cons_s + 1e-10){
    dot_cons = dot_cons_d
    dot_v = (yd + dot_yd)*(1 - alpha1) - alpha2*v
    dot_b_h = (v + dot_v)*(lambda0 + lambda1*r) - lambda2*(yd + dot_yd) - b_h/365 ##!!
    # dot_b_s = G + dot_g - (yd + dot_yd)*theta/(1-theta) + r*b_h
    # dot_b_cb = dot_b_s - dot_b_h
    # if(t<12) print(c(1,t))
  }else{
    # scale down all consumption
    offprop = alpha1offline/(alpha1online + alpha1offline)
    dot_cons = dot_cons_s
    dot_cons_off = dot_cons * offprop
    dot_cons_on = dot_cons - dot_cons_off
    # if(t<12) print(c(2,t))
    
    # print(c(2,offprop,cons_s,cons_d,prop_to_work,39377572-econ$lf))
    # tax = next_yd*theta/(1-theta)
    next_yd = (wb_s + r*b_h)/(1+theta/(1-theta))
    dot_yd = next_yd - yd
    dot_v = next_yd - cons_s
    
    # scale_alpha = (dot_cons_s + cons_s)/(dot_cons_d + cons_d)
    # alpha1 = scale_alpha*alpha1
    # alpha2 = scale_alpha*alpha2
    # alpha1offline = scale_alpha*alpha1offline
    # alpha2offline = scale_alpha*alpha2offline
    # # compute econ again
    # dot_yd = (alpha2*v + G + r*b_h)/(1 + theta/(1-theta) - alpha1) - yd
    # dot_v = (yd + dot_yd)*(1 - alpha1) - alpha2*v
    
    dot_b_h = (v + dot_v)*(lambda0 + lambda1*r) - lambda2*(yd + dot_yd) - b_h/365 ##!!
    
  }
  econ_derivs = c(dot_b_h, #Government bills held by households
                  dot_cons_on, #Consumption online
                  dot_cons_off, #Consumption offline
                  dot_cons, # consumption
                  dot_v, #Households wealth
                  dot_yd #Disposable income of households
  )
  
  return(econ_derivs)
  
}


