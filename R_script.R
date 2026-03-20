library(pracma)
library(MASS)
library(abind)
library(Matrix)
library(haven)
library(ggplot2)
library(viridis)
library(latex2exp)
library(data.table)

# set up basics

set.seed(0)
setwd(getSrcDirectory(function(){})[1])

# load functions
source('R_functions.R');

# read or create epi and country structures
datafile <- 'data_file.Rds'
if(!file.exists(datafile)){
  # library(Rilostat)
  library(conmat)
  library(wpp2024)
  library(wbstats)
  library(WDI)
  
  ## starting variables ###############################
  data <- data_start()
  
  ## country data #################################
  ldata <- gather_country_data(data)
  
  # save:
  saveRDS(ldata,datafile)
  
}else{
  ldata <- readRDS(datafile)
}

## load econ models, which are written into file econ_models.R
source('econ_models.R')
# choose an econ model
econ = model2

## simulate ###########################################################

nscen = 7
parameter_combinations = data.frame(q1 = rep(0.86, nscen), 
                                    q2 = rep(0.0001, nscen),
                                    lambdap = rep(0.01, nscen),
                                    cc = rep(1,nscen),
                                    isolate = rep(1,nscen))
parameter_combinations$q2[2] = 0.00001
parameter_combinations$q2[3] = 0.001
parameter_combinations$lambdap[4] = 0.001
parameter_combinations$lambdap[5] = 0.02
parameter_combinations$cc[6] = 0.5
parameter_combinations$isolate[7] = 0
for(rv in 1:nrow(parameter_combinations))
{
  cat(paste0('Parameter combination ',rv,' out of ',nrow(parameter_combinations),'\n'))
  econ$lambda_p1 = parameter_combinations$lambdap[rv]
  ldata$q1 = parameter_combinations$q1[rv]
  ldata$q2 = parameter_combinations$q2[rv]
  ldata$cc = parameter_combinations$cc[rv]
  ldata$epidemic$prob_isolated = parameter_combinations$isolate[rv]
  ldata = decompose_contacts(ldata, consumption_contact = parameter_combinations$cc[rv])
  
  ## run model
  runlist <- simulate_epi_econ_model(data = ldata, econ = econ)
  
  plotout = plot_trajectories(runlist, econ, ldata)
  
  if(econ$model_name == 'Model 1'){
    label = TeX(paste0(econ$model_name,"; $q_2$ = ", ldata$q2))
  }else{
    label = TeX(paste0(econ$model_name,"; $\\gamma_1 = ", econ$lambda_p1,
                       "$; $q_2$ = ", ldata$q2,
                       "$; $p_i$ = ", ldata$epidemic$prob_isolated,
                       "$; $\\eta$ = ", parameter_combinations$cc[rv]))
  }
  plotout = plotout + labs(title= label)
  ggsave(plotout,filename=paste0('figures/',econ$model_name,'_',
                                 econ$lambda_p1,'-',
                                 ldata$q2,'-',
                                 ldata$epidemic$prob_isolated,'-',
                                 parameter_combinations$cc[rv],'.png'),width=12,height=3.5)
  
}


