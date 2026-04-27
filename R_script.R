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
  library(furrr)
  
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

nscen = 9
parameter_combinations = data.frame(q1 = rep(0.862752, nscen), 
                                    q2 = rep(0.000237, nscen),
                                    lambdap = rep(0.75, nscen),
                                    cc = rep(0.4,nscen))
parameter_combinations$q1[2] = 0.73
parameter_combinations$q1[3] = 0.93
parameter_combinations$q2[4] = 0.0000001
parameter_combinations$q2[5] = 0.01
parameter_combinations$lambdap[6] = 0.5
parameter_combinations$lambdap[7] = 0.95
parameter_combinations$cc[8] = 0.2
parameter_combinations$cc[9] = 0.6
for(rv in 1:nrow(parameter_combinations))
{
  cat(paste0('Parameter combination ',rv,' out of ',nrow(parameter_combinations),'\n'))
  econ$lambda_p1 = parameter_combinations$lambdap[rv]
  ldata$q1 = parameter_combinations$q1[rv]
  ldata$q2 = parameter_combinations$q2[rv]
  ldata$cc = parameter_combinations$cc[rv]
  ldata = decompose_contacts(ldata, consumption_contact = parameter_combinations$cc[rv])
  print(ldata$epidemic$beta)
  
  ## run model
  runlist <- simulate_epi_econ_model(data = ldata, econ = econ)
  
  plotout = plot_trajectories(runlist, econ, ldata)
  
  if(econ$model_name == 'Model 1'){
    label = TeX(paste0(econ$model_name,"; $q_2$ = ", ldata$q2))
  }else{
    label = TeX(paste0(#econ$model_name,'; ',
                       "$\\gamma_1 = ", econ$lambda_p1,
                       "$; $q_1 = ", ldata$q1,
                       "$; $q_2 = ", ldata$q2,
                       "$; $\\eta = ", parameter_combinations$cc[rv],"$"))
  }
  plotout = plotout + labs(title= label)
  ggsave(plotout,filename=paste0('figures/',econ$model_name,'_',
                                 econ$lambda_p1,'-',
                                 ldata$q1,'-',
                                 ldata$q2,'-',
                                 parameter_combinations$cc[rv],'.png'),width=12,height=3.5)
  
}


