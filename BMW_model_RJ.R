# Model BMW for R (from scratch): main code
# from Wynne Godley and Marc Lavoie
# Monetary Economics
# Chapter 7 

# Version: 30 May 2019; revised: 6 November 2022

library(deSolve)

################################################################################

#STEP 1: PREPARE THE WORK-SPACE

#Clear environment
rm(list=ls(all=TRUE))

#Clear plots
# if(!is.null(dev.list())) dev.off()

#Clear console
# cat("\014")

#Number of periods
nPeriods = 150

#Number of scenarios
nScenarios=6 

################################################################################

#STEP 2: SET THE COEFFICIENTS (EXOGENOUS VARIABLES AND PARAMETERS)
alpha2=0.1 #Propensity to consume out of wealth
delta=0.1 #Depreciation rate
gamma=0.15 #Reaction speed of adjustment of capital to its target value
alpha1r = 0.5 #Propensity to consume out of interests if c_d=c_d(r)
alpha1w = 0.5 #Propensity to consume out of wages if c_d=c_d(r)
alpha0=0 # matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Set autonomous component of consumption as a matrix (for shocks) 
alpha1=0.5 # matrix(data=0.75,nrow=nScenarios,ncol=nPeriods) #Set autonomous propensity to consume out of income as a matrix (for shocks) 
kappa=1 # matrix(data=1,nrow=nScenarios,ncol=nPeriods)  #Capital-output ratio
rl_bar=0.04 # matrix(data=0.04,nrow=nScenarios,ncol=nPeriods)  #Rate of interests on bank loans - exogenously set

################################################################################
c_d = 0
i_d = 0
k = 0
l_d = 0
mh = 0
yd = 0
x = c(c_d,i_d,k,l_d,mh,yd)

#STEP 3: CREATE VARIABLES AND ATTRIBUTE INITIAL VALUES
# af=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Create variables and attribute values to stocks
c_d=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Consumption goods demand by households
# c_s=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Consumption goods supply
# da=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Depreciation allowances
k=matrix(data=0,nrow=nScenarios,ncol=nPeriods)  #Stock of capital
# kt=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Target stock of capital
l_d=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Demans for bank loans
# l_s=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Supply of bank loans 
i_d=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Demand for Investment
# i_s=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Supply of Investment
mh=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Bank deposits held by households
ms=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Supply of bank deposits
# n_s=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Supply of labour
# n_d=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Demand for labour
# pr=matrix(data=1,nrow=nScenarios,ncol=nPeriods) #Labour productivity
# rl=matrix(data=rl_bar,nrow=nScenarios,ncol=nPeriods) #Rate of interests on banks loans
# w=matrix(data=0.86,nrow=nScenarios,ncol=nPeriods) #Wage rate
# wb_d=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Wage Bill
# wb_s=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Supply of Wages
# y=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Income
# y_star=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Steady-state income
yd=matrix(data=0,nrow=nScenarios,ncol=nPeriods) #Disposal Income of households
       
################################################################################

#STEP 4: RUN THE MODEL

#STEP 4.A: DEFINE THE LOOPS
j = 1



bmw_model_ode = function(t,x,p){
#Choose scenario
#
  j = p$j
  #Define time loop
  #for (i in 2:nPeriods-1){
    
    #Define iterations
    #for (iterations in 1:100){
      
      
      #STEP 4.B: ADD ALTERNATIVE SCENARIOS  
      
      #Initial shock to autonomous consumption
      if (t>=2){alpha0=25}    
      
      #Scenario 2: positive shock to target capital stock
      if (t>=52 & j==2){kappa=1.1}
      
      #Scenario 3: negative shock to target capital stock
      if (t>=52 & j==3){kappa=0.9}
      
      #Scenario 4: shock to autonomous consumption
      if (t>=52 & j==4){alpha0=28}
      
      #Scenario 5: shock to interest rate with c_d = c_d(r)
      if (t>=52 & j==5){rl_bar=0.05}
      
      #Scenario 6: increase in propensity to save
      if (t>=52 & j==6){alpha1=0.74}
      
      ri=rl_bar 
      
      #STEP 4.C: DEFINE THE SYSTEM OF EQUATIONS        
      
      
      
      c_d = x[1]
      i_d = x[2]
      k = x[3]
      l_d = x[4]
      mh = x[5]
      yd = x[6]
      
      integrate = p$integrate
      
      # SIR
      
      Susc = x[7];
      I = x[8];
      beta = 1/2;
      recovery = 1/5;
      
      
      #Behavioural equations and equilibrium condition
      
        #Supply of consumption goods - eq. 7.1
        c_s = c_d
        #Supply of investment goods - eq. 7.2
        i_s = i_d
        #Supply of labour - eq. 7.3
        # n_s = n_d
        #Supply of loans - eq. 7.4
        # dot_l_s = dot_l_d 
        
        
        #GDP - eq. 7.5
        y = c_s + i_s
        
        #The investment behaviour
        
        #Accumulation of capital - eq. 7.17
        dot_k = gamma*(kappa*y - k) 
        
        #Depreciation allowances - eq. 7.18
        # da[j,i+1] = delta*k
        # dot_da = delta*k - da
        
        #Capital stock target - eq. 7.19
        # kt[j,i+1] = kappa*y
        
        #Demand for investment goods - eq. 7.20
        # i_d[j,i+1] = gamma*(kappa*y - k) + delta*k
        dot_i_d = gamma*(kappa*y - k) + delta*k - i_d
        
        #Transactions of the banks  
        
        #Supply of deposits - eq. 7.11
        #   dot_ms = dot_l_s
        # 
        #   #Rate of interest on deposits - eq. 7.12
        rl = ri
        # 
      
      #Transactions of the firms
      
        #Wage bill - eq. 7.6
        # wb_d[j,i+1] = y[j,i+1] - rl*l_d - delta*k   
        
        #Depreciation allowances - eq. 7.7
        # af = da
        #Demand for bank loans - eq. 7.8
        dot_l_d = dot_k
      
      #Transactions of households
        
        #Disposable income - eq. 7.9
        # yd[j,i+1] = wb_d[j,i+1] + ri*mh
      
        #Bank deposits held by households - eq. 7.10
        # dot_mh = yd[j,i+1] - c_d[j,i+1]
        dot_mh = ri*mh + i_d + dot_i_d - rl*l_d - delta*k
        
        #Household behaviour
        
        #Demand for consumption goods - eq. 7.16 (and 7.16A)
        if (j==5){ 
          if (integrate==1){
            props = prop_to_consume(I,alpha1w,alpha2,beta)
            alpha1w = props$alpha_1
            alpha2 = props$alpha_2
            beta = props$beta
          }
          dot_yd = (dot_mh + alpha0 - alpha1w*ri*mh + alpha1r*ri*mh + alpha2*mh)/(1 - alpha1w) - yd
          dot_c_d = yd + dot_yd - dot_mh - c_d 
        }else{ 
          if (integrate==1){
            props = prop_to_consume(I,alpha1,alpha2,beta)
            alpha1 = props$alpha_1
            alpha2 = props$alpha_2
            beta = props$beta
          }
          dot_yd = (dot_mh + alpha0 + alpha2*mh)/(1 - alpha1) - yd
          # print(c(dot_yd,alpha1,alpha0))
          dot_c_d = yd + dot_yd - dot_mh - c_d 
        }
        
        Sdot = - beta*Susc*I/100000;
        Idot = - Sdot - recovery*I;
        
        # dot_wb_d = wb_d - ri*mh + yd + dot_yd     
        
        #The wage bill
        
        # dot_y = dot_i_d + dot_c_d
        # 
        # #Labour demand - eq. 7.14
        # # n_d = y/pr
        # dot_n_d = dot_y/pr
        # 
        # #Wage rate - eq. 7.15
        # # w = wb_d/n_d
        # dot_w = (dot_wb_d*n_d - wb_d*dot_n_d)/n_d^2
        
        #"Supply" of wages - eq. 7.13
        # wb_s = w*n_s
        # 
      # #The behaviour of banks
      # 
      #   #Interest rate on loans - eq. 7.21
      #   rl = rl_bar
      # 
      # #Additional equation for calculating y steady-state value
      #   y_star = alpha0/((1-alpha1)*(1-delta*kappa)-alpha2*kappa)       
      
        dot = c(dot_c_d,
                dot_i_d,
                dot_k,
                dot_l_d,
                dot_mh,
                dot_yd,
                Sdot,
                Idot#Disposable income of households
        )
        
        x = x+dot
        list(dot)
        
}

prop_to_consume = function(I,alpha_1,alpha_2,beta){
  
  baseline = .5;
  logistic = 1/(1+exp(-(10000-I)/1000));
  
  prop_t_consume = baseline + (1-baseline) * logistic;
  
  alpha_1 = prop_t_consume * alpha_1;
  alpha_2 = prop_t_consume * alpha_2;
  beta = prop_t_consume^2 * beta;
  
  list(alpha_1=alpha_1,alpha_2=alpha_2,beta=beta)
  
}

y0 = c(100000-1, 1) 

for (j in 1:nScenarios){
  odeout = ode(y=c(x,y0), times=c(1:nPeriods)-1, func=bmw_model_ode,list(j=j,integrate=1));
  
  c_d[j,] = odeout[,2]
  i_d[j,] = odeout[,3]
  k[j,] = odeout[,4]
  l_d[j,] = odeout[,5]
  mh[j,] = odeout[,6]
  yd[j,] = odeout[,7]
  ms[j,] = l_d[j,]
}
y = c_d + i_d
da = delta*k
y_star = 25/((1-alpha1)*(1-delta*kappa)-alpha2*kappa)  

################################################################################

#STEP 5: CONSISTENCY CHECK (REDUNDANT EQUATION)

#Create consistency statement 
aerror=0
error=0
for (i in 2:(nPeriods-1)){error = error + (ms[i+1]-mh[i+1])^2}
aerror = error/nPeriods
if ( aerror<0.01 ){cat(" *********************************** \n Good news! The model is watertight! \n", "Average error =", aerror, "< 0.01 \n", "Cumulative error =", error, "\n ***********************************\n")} else{
  if ( aerror<1 && aerror<1 ){cat(" *********************************** \n Minor issues with model consistency \n", "Average error =", aerror, "> 0.01 \n", "Cumulative error =", error, "\n ***********************************\n")}
  else{cat(" ******************************************* \n Warning: the model is not fully consistent! \n", "Average error =", aerror, "> 1 \n", "Cumulative error =", error, "\n *******************************************\n")} }      

#Plot redundant equation
layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow = TRUE))
par(mar = c(5.1+1, 4.1+1, 4.1+1, 2.1+1))
plot(ms[1,2:nPeriods]-mh[1,2:nPeriods], type="l", col="green",lwd=3,lty=1,font.main=1,cex.main=1.5,
     main=expression("Consistency check (baseline scenario): " * italic(M[phantom("")]["s"]) - italic(M[phantom("")]["h"])),
     cex.axis=1.5,cex.lab=1.5,ylab = '£',
     xlab = 'Time',ylim = range(-1,1))

################################################################################

#STEP 6: CREATE AND DISPLAY GRAPHS

#Change layout
#layout(matrix(c(1,2), 1, 2, byrow = TRUE))
par(mar = c(5.1+1, 4.1+1, 4.1+1, 2.1+1))

#Figure 1
plot(y[1,2:145],type="l",col=1,lwd=3,lty=1,font.main=1,cex.main=1.5,
     main="Figure 1  Evolution of national income \n following initial autonomous consumption",
     cex.axis=1.5,cex.lab=1.5,ylab = '£',xlab = 'Time',ylim=range(0,80))
abline(h=y_star,type="l",lwd=3,lty=2,col=4)
legend("bottomright",c("National income","Steady-state value"),  bty = "n", cex = 1.5,
       lty=c(1,2), lwd=c(3,3), col = c(1,4), box.lwd=0)

#Figure 2
plot(yd[4,48:140],type="l",col="aquamarine",lwd=3,lty=1,font.main=1,cex.main=1.5,
     main="Figure 2  Evolution of disposable income and consumption \n following shock to autonomous consumption",
     cex.axis=1.5,cex.lab=1.5,ylab = '£',xlab = 'Time',ylim=range(0,100))
lines(c_d[4,48:140],type="l",lwd=3,lty=3,col="aquamarine4")
legend("bottomright",c("Disposable income","Consumption"),  bty = "n", cex = 1.5,
       lty=c(1,3), lwd=c(3,3), col = c("aquamarine","aquamarine4"), box.lwd=0)

#Figure 3
plot(i_d[4,48:140],type="l",col="orchid",lwd=3,lty=1,font.main=1,cex.main=1.5,
     main="Figure 3  Evolution of investment and depreciation \n following shock to autonomous consumption",
     cex.axis=1.5,cex.lab=1.5,ylab = '£',xlab = 'Time',ylim=range(6,9))
lines(da[4,48:140],type="l",lwd=3,lty=2,col="gray50")
legend("right",c("Gross investment","Replacement investment \n (capital depreciation)"),  bty = "n", cex = 1.5,
       lty=c(1,2), lwd=c(3,3), col = c("orchid","gray50"), box.lwd=0)

#Figure 4
plot(yd[6,48:140],type="l",col="aquamarine",lwd=3,lty=1,font.main=1,cex.main=1.5,
     main="Figure 4  Evolution of disposable income and consumption \n following increase in propensity to save (paradox ot thrift!)",
     cex.axis=1.5,cex.lab=1.5,ylab = '£',xlab = 'Time',ylim=range(165,180))
lines(c_d[6,48:140],type="l",lwd=3,lty=3,col="aquamarine4")
legend("bottomright",c("Disposable income","Consumption"),  bty = "n", cex = 1.5,
       lty=c(1,3), lwd=c(3,3), col = c("aquamarine","aquamarine4"), box.lwd=0)

#Figure 5
plot(y[2,1:140],type="l",col=4,lwd=3,lty=1,font.main=1,cex.main=1.5,
     main="Figure 5  Evolution of national income \n following change in target capital stock",
     cex.axis=1.5,cex.lab=1.5,ylab = '£',xlab = 'Time',ylim=range(min(y[3,40:150]),max(y[2,40:150])))
mycol3 <- rgb(255,255,0, max = 255, alpha = 40, names = "myyellow")
lines(y[3,1:140],type="l",lwd=3,lty=1,col=3)
lines(y[1,1:140],type="l",lwd=3,lty=2,col=1)
#rect(xleft=12,xright=150,ybottom=30,ytop=130,col=mycol3,border=NA)
legend(40,105,c("Higher capital/output","Lower capital/output","Baseline value"),  bty = "n", cex = 1.5,
       lty=c(1,1,2), lwd=c(3,3,3), col = c(4,3,1), box.lwd=0)

#Figure 6
plot(y[5,48:140],type="l",col="coral3",lwd=3,lty=1,font.main=1,cex.main=1.5,
     main="Figure 6  Evolution of national income \n following an increase in the interest rate",
     cex.axis=1.5,cex.lab=1.5,ylab = '£',xlab = 'Time',ylim=range(60,80))
abline(h=y[5,48])
