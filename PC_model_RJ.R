# Model PC for R: main code
# from Wynne Godley and Marc Lavoie
# Monetary Economics
# Chapter 4 

# Version: 30 May 2019; revised: 6 November 2023
# setwd(getSrcDirectory(function(){})[1])
library(deSolve)


################################################################################
{
#STEP 1: PREPARE THE WORK-SPACE

#Clear environment
rm(list=ls(all=TRUE))

#Clear plots
# if(!is.null(dev.list())) dev.off()

#Clear console
# cat("\014")

#Number of periods
nPeriods = 90

#Number of scenarios
nScenarios=3 
}
################################################################################
{
#STEP 2: SET THE COEFFICIENTS (EXOGENOUS VARIABLES AND PARAMETERS)
alpha1=0.6 #Propensity to consume out of income
alpha2=0.4 #Propensity to consume out of wealth
lambda0 = 0.635 #Autonomous share of bills
lambda1 = 5 #Elasticity of bills demand to interest rate
lambda2 = 0.01 #Elasticity of bills demand to yd/v
theta=0.2 #Tax rate on income
r_bar=0.025 #Interest rate as policy instrument

# constants
r=r_bar #Interest rate on government bills
dot_r = 0
prev_r = r - dot_r
g=20; dot_g = 0 #Government expenditure

# end values
cons=64.87+21.62 #Consumption
# h_h=21.62 #Cash money held by households
# h_s=21.62 #Cash money supplied by central bank
# y=0 #Income = GDP
}
################################################################################
{
#STEP 3: CREATE VARIABLES AND ATTRIBUTE INITIAL VALUES
  b_cb=21.62 #Government bills held by Central Bank
  b_h=64.87 #Government bills held by households
  b_s=21.62+64.87 #Government bills supplied by government
  # tax=0 #Taxes
  v=64.87+21.62 #Households wealth
  yd=64.87+21.62 #Disposable income of households
}
################################################################################

#STEP 4: RUN THE MODEL

#STEP 4.A: DEFINE THE LOOPS

#Choose scenario

x = c(b_cb,b_h,b_s,cons,v,yd)


pc_model_ode = function(t,x,p){
  
  b_cb = x[1]
  b_h = x[2]
  b_s = x[3]
  cons = x[4]
  v = x[5]
  yd = x[6]
  
  
  # SIR
  
  Susc = x[7];
  I = x[8];
  beta = 1/2;
  recovery = 1/5;
  
  
  j = p$j
  integrate = p$integrate
  
  
  #Shock 1: higher interest rate
  if (t>=10 && j==2){r=0.035}
  
  #Shock 2: higher propensity to consume out of income
  if (t>=10 && j==3){alpha1=0.7}
  
  
  if (integrate==1){
    props = prop_to_consume(I,alpha1,alpha2,beta)
    alpha1 = props$alpha_1
    alpha2 = props$alpha_2
    beta = props$beta
  }
  
  Sdot = - beta*Susc*I/100000;
  Idot = - Sdot - recovery*I;
  
  #STEP 4.C: DEFINE THE SYSTEM OF EQUATIONS 
  
  #Tax payments - eq. 4.3
  # t[j,i+1] = theta*(y[j,i+1] + r*b_h)
  # t[j,i+1] = theta*(yd[j,i+1] + t[j,i+1])
  # t = theta*(yd + t)
  tax = yd*theta/(1-theta)
  
  #Disposable income - eq. 4.2
  # yd[j,i+1] = y[j,i+1] - t[j,i+1] + r*b_h
  dot_yd = (alpha2*v + g + r*b_h)/(1 + theta/(1-theta) - alpha1) - yd
  
  #Determination of output - eq. 4.1
  y = cons + g
  # y[j,i+1] = alpha1*yd[j,i+1] + alpha2*v + g
  # dot_y = alpha1*(yd + dot_yd) + alpha2*v + g + dot_g - y
  
  #Wealth accumulation - eq. 4.4
  # dot_v = yd[j,i+1] - cons[j,i+1]
  dot_v = (yd + dot_yd)*(1 - alpha1) - alpha2*v
  
  #Consumption function - eq. 4.5
  dot_cons = alpha1*(yd + dot_yd) + alpha2*v - cons
  
  #Demand for government bills - eq. 4.7
  dot_b_h = (v + dot_v)*(lambda0 + lambda1*r) - lambda2*(yd + dot_yd) - b_h
  
  #Supply of government bills - eq. 4.8
  # dot_b_s = g[j,i+1] + r*b_s - t[j,i+1] - r*b_cb
  dot_b_s = g + dot_g - (yd + dot_yd)*theta/(1-theta) + r*b_h
  
  #Government bills held by the central bank - eq. 4.10
  dot_b_cb = dot_b_s - dot_b_h
  
  dot = c(dot_b_cb, #Government bills held by Central Bank
          dot_b_h, #Government bills held by households
          dot_b_s, #Government bills supplied by government
          dot_cons, #Consumption
          # dot_g, #Government expenditure
          # dot_h_h, #Cash money held by households
          # dot_h_s, #Cash money supplied by central bank
          # dot_r, #Interest rate on government bills
          # dot_tax, #Taxes
          dot_v, #Households wealth
          # dot_y, #Income = GDP
          dot_yd, #Disposable income of households
          Sdot,
          Idot
  )
  
  
  x = x+dot
  list(dot)
  
}

prop_to_consume = function(I,alpha_1,alpha_2,beta){
  
  baseline = .5
  logistic = 1/(1+exp(-(10000-I)/1000));
  
  prop_t_consume = baseline + (1-baseline) * logistic;
  
  alpha_1 = prop_t_consume * alpha_1;
  alpha_2 = prop_t_consume * alpha_2;
  beta = prop_t_consume^2 * beta;
  print(c(alpha_1,alpha_2))
  
  list(alpha_1=alpha_1,alpha_2=alpha_2,beta=beta)
  
}

y0 = c(100000-1, 1) 

nScenarios = 3

b_cb=matrix(data=21.62,nrow=nScenarios,ncol=nPeriods) #Government bills held by Central Bank
b_h=matrix(data=64.87,nrow=nScenarios,ncol=nPeriods) #Government bills held by households
b_s=b_cb+b_h #Government bills supplied by government
cons=b_s #Consumption 
# g=matrix(data=20,nrow=nScenarios,ncol=nPeriods) #Government expenditure
h_h=b_cb #Cash money held by households
h_s=b_cb #Cash money supplied by central bank
# r=matrix(data=r_bar,nrow=nScenarios,ncol=nPeriods) #Interest rate on government bills
t=b_cb #Taxes
v=b_s #Households wealth
y=cons+g #Income = GDP
yd=b_s #Disposable income of households


for(j in 1:nScenarios){
  odeout = ode(y=c(x,y0), times=c(1:90), func=pc_model_ode,list(j=j,integrate=0));
  
  b_cb[j,] = odeout[,2]
  b_h[j,] = odeout[,3]
  b_s[j,] = odeout[,4]
  cons[j,] = odeout[,5]
  v[j,] = odeout[,6]
  yd[j,] = odeout[,7]
  
  #Supply of cash - eq. 4.9
  h_s[j,] = b_cb[j,]
  #Cash money - eq. 4.6
  h_h[j,] = v[j,] - b_h[j,]
  
  y[j,] = cons[j,] + g
}
#Interest rate as policy instrument - eq. 4.11
# r = r_bar

################################################################################

#STEP 5: CONSISTENCY CHECK (REDUNDANT EQUATION)

#Create consistency statement 
aerror=0
error=0
for (i in 2:(nPeriods-1)){error = error + (h_s[i]-h_h[i])^2}
aerror = error/nPeriods
if ( aerror<0.01 ){cat(" *********************************** \n Good news! The model is watertight! \n", "Average error =", aerror, "< 0.01 \n", "Cumulative error =", error, "\n ***********************************\n")} else{
  if ( aerror<1 && aerror<1 ){cat(" *********************************** \n Minor issues with model consistency \n", "Average error =", aerror, "> 0.01 \n", "Cumulative error =", error, "\n ***********************************\n")}
  else{cat(" ******************************************* \n Warning: the model is not fully consistent! \n", "Average error =", aerror, "> 1 \n", "Cumulative error =", error, "\n *******************************************\n")} }      

#Plot redundant equation
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(mar = c(5.1+1, 4.1+1, 4.1+1, 2.1+1))
plot(h_s[2:nPeriods]-h_h[2:nPeriods], type="l", col="green",lwd=3,lty=1,font.main=1,cex.main=1.5,
     main=expression("Consistency check (baseline scenario): " * italic(H[phantom("")]["s"]) - italic(H[phantom("")]["h"])),
     cex.axis=1.5,cex.lab=1.5,ylab = '£',
     xlab = 'Time',ylim = range(-1,1))

################################################################################

#STEP 6: CREATE AND DISPLAY GRAPHS

#Set layout
#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(mar = c(5.1+1, 4.1+1, 4.1+1, 2.1+1))

#Figure 1
plot(100*h_h[2,]/v[2,],type="l", col=1, lwd=3, lty=1, font.main=1,cex.main=1.5,
     main="Figure 1  Evolution of shares of bills and money following \n an increase in 100 points in the rate of interest on bills",
     ylab = '%',xlab = 'Time',ylim = range(19,26),cex.axis=1.5,cex.lab=1.5)
par(new=T)
plot(100*b_h[2,]/v[2,],type="l",lwd=3,lty=2,col=rgb(0,187,110,max=255),axes=F, ylab=NA,
     xlab=NA,ylim=range(74,81))
axis(side = 4,cex.axis=1,cex.lab=1.5)
legend("right",c("Share of money balances","Share of bills (right axis)"),  bty = "n",
       cex=1.5, lty=c(1,2), lwd=c(3,3), col = c(1,rgb(0,187,110,max=255)), box.lty=0)

#Figure 2
plot(yd[2,],type="l", col=1, lwd=3, lty=1, font.main=1,cex.main=1.5,
     main="Figure 2  Evolution of household disposable income and consumption \n following an increase in 100 points in the rate of interest on bills",
     ylab = '£',xlab = 'Time',ylim = range(60,110),cex.axis=1.5,cex.lab=1.5)
lines(cons[2,],type="l",lwd=3,lty=2,col=4)
legend("right",c("Disposable income","Consumption"),  bty = "n",
       cex=1.5, lty=c(1,2), lwd=c(3,3), col = c(1,4), box.lty=0)

#Figure 3 - Note: different from 4.5 as expectations are not included
plot(y[3,],type="l", col=2, lwd=3, lty=1, font.main=1,cex.main=1.5,
     main="Figure 3  Evolution of national income (GDP) following an increase in \n the propensity to consume out of expected disposable income",
     ylab = '£',xlab = 'Time',ylim = range(70,130),cex.axis=1.5,cex.lab=1.5)
abline(h=y[3,2],col="green",lty=2,lwd=2)
abline(h=y[3,60],col="blue",lty=2,lwd=2)
legend("right",c("National income","Old steady state","New steady state"),  bty = "n",
       cex=1.5, lty=c(1,2,2), lwd=c(3,2,2), col = c(2,"blue","green"), box.lty=0)
