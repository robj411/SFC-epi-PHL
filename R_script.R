library(pracma)
library(MASS)
library(abind)
library(Matrix)
library(haven)
library(ggplot2)
library(viridis)
library(latex2exp)

# set up basics
country <- 'Philippines'
iso3 = 'PHL'
set.seed(0)
setwd(getSrcDirectory(function(){})[1])

# load functions
source('R_functions.R');

# read or create epi and country structures
datafile <- 'data_file.Rds'
if(!file.exists(datafile)){
  library(Rilostat)
  library("conmat")
  library("wpp2024")
  library(wbstats)
  
  ## starting variables ###############################
  data <- data_start()
  
  ## disease variables ############################
  
  dis <- list()
  R0 <- 2.75 
  dis$TEtoI <- 4.6 
  dis$TItoR <- 4 
  dis$TItoC <- 2
  dis$TCtoR <- dis$TItoR - dis$TItoC
  dis$prob_detected <- 0.5
  dis$prob_isolated <- 1
  
  ## country data #################################
  
  ldata <- gather_country_data(data, country_name=country, iso3=iso3)
  # combine country and disease parameters
  CI <- get_candidate_infectees(dis, ldata)
  dis$beta <- R0 / CI
  
  ## complete p2, dis and data structs #################################
  
  ldata_dis <- list(data=ldata, dis=dis)
  # save:
  saveRDS(ldata_dis,datafile)
  
}else{
  ldata_dis <- readRDS(datafile)
}

## load stored objects
dis2 <- ldata_dis$dis
ldata <- ldata_dis$data

## load econ models, which are written into file econ_models.R
source('econ_models.R')
# choose an econ model
econ = model2

# plot response function to epidemic for reference
epivars = seq(0,3e5,by=1000)
eplot = epi_to_econ(epivars,model1)
# (plotresponse <- ggplot()  + 
#     geom_line(aes(x=epivars/1e3,y=eplot$scalar),linewidth=2,colour='midnightblue') +
#   theme_bw(base_size = 15) +
#   labs(x='Thousand confirmed cases',y='Relative propensity to consume') +
#   scale_y_continuous(limits=c(0,1)) +
#     geom_vline(xintercept=100,colour='grey',linewidth=1.5) +
#     geom_hline(yintercept=.5,colour='grey',linewidth=1.5) + 
#     geom_line(aes(x=epivars/1e3,y=eplot$scalar),linewidth=2,colour='midnightblue'))
# ggsave(plotresponse,filename='figures/response.png',width=5,height=5)

## simulate ###########################################################

q1_vals = c(0.86)
q2_vals = c(0.0001, 0.00001)
lambdap_vals = c(0.001, 0.01)
outtab = data.frame(expand.grid(lambdap_vals, q1_vals, q2_vals))
for(rv in 1:nrow(outtab))
{
  cat(paste0('Parameter combination ',rv,' out of ',nrow(outtab),'\n'))
  econ$lambda_p1 = outtab[rv,1]
  ldata$q1 = outtab[rv,2]
  ldata$q2 = outtab[rv,3]
  ## run model
  runlist <- simulate_epi_econ_model(data = ldata, dis = dis2, econ = econ)
  
  plotout = plot_trajectories(runlist, ldata)
  
  if(econ$model_name == 'Model 1'){
    label = TeX(paste0(econ$model_name,"; $q_1$ = ", ldata$q2))
  }else{
    label = TeX(paste0(econ$model_name,"; $\\gamma_1 = ", econ$lambda_p1,"$; $q_1$ = ", ldata$q2))
  }
  plotout = plotout + labs(title= label)
  ggsave(plotout,filename=paste0('figures/',econ$model_name,'_',econ$lambda_p1,'-',ldata$q2,'.png'),width=12,height=3.5)
  
  
  # gdpscen = trapz(Tout,econ$get_gdp_from_out(odevar,econ,C,ldata,1))
  # trapz(Tout0,cntr[,which(econ$econvarnames=='cons')+1] + econ$g)
}


