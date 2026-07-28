
## Model 1, based on SIM (Godley and Lavoie, ch 3) ####################################

model1 <- list()
model1$model_name = 'model1'

# get initial conditions using equations
y0    = ldata$gdp/365  # from annual to daily, in billions
tax0  = ldata$tax/365
theta = tax0 / y0
yd0   = (1-theta)*y0
cons0 = yd0
g0    = y0 - cons0
h0    = 19032.66/1e3   # in billions

# store parameters
model1$econ_init = h0
model1$G      = g0                     # government spending (constant)
model1$theta  = theta                  # tax rate on income
model1$alpha1 = 0.7695          # propensity to consume out of income
model1$alpha2 = 0.0105*4/365       # propensity to consume out of wealth (daily)
model1$alpha0 = with(model1, cons0 - h0*alpha2 - yd0*alpha1)
print(model1$alpha0)
model1$gdp    = y0*365
model1$y0     = model1$gdp/365
##!! change to people
model1$lf     = ldata$NNs[1]/1e6     # labour force in millions
model1$emp0   = model1$lf*ldata$employmentrate/100
##!! change to per million people
model1$lambda = model1$y0/model1$emp0
model1$cons0  = cons0

# variable names
model1$econvarlabels = c('Wealth')
model1$econvarnames  = c('H_h')
model1$nEconODEs     = length(model1$econvarnames)


## Model 2 (Model 1 plus time-varying productivity) ####################################

model2 <- model1
model2$model_name = 'model2'

model2$lambda_0  = model1$lambda
model2$lambda_p0 = 0     # intrinsic growth of lambda
model2$lambda_p1 = 0.5  # rate that lambda tracks output

model2$econ_init = c(h0, model2$lambda_0)

# variable names
model2$econvarlabels = c('Wealth','Productivity')
model2$econvarnames  = c('H_h','lambda')
model2$nEconODEs     = length(model2$econvarnames)
