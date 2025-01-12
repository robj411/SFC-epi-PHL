library(readxl)
library(pracma)
library(MASS)
library(abind)
library(Matrix)
library(haven)
library(ggplot2)
library(wbstats)
library(viridis)

country <- 'Philippines'
set.seed(0)
setwd('~/projects/DAEDALUS/SFC-epi-PHL/')
source('R_functions.R');
 

datafile <- 'data_file.Rds'

if(!file.exists(datafile)){
  
  ## country variables ###############################
  
  CD_list <- load_country_data()
  CD <- CD_list[[1]]
  country_parameter_distributions <- CD_list[[2]]
  data <- data_start()
  
  ## disease variables ############################
  
  sevenpathogens <- read.csv('../Daedalus-P2-Dashboard/data/sevenpathogens.csv')
  
  dis <- list()
  dis$R0 <- sevenpathogens$R0[5]
  dis$asym_rel_transmission <- sevenpathogens$red[5]
  dis$TRtoS <- sevenpathogens$Ti[5] * 10e5
  dis$TEtoI <- sevenpathogens$Tlat[5]
  dis$THtoR <- sevenpathogens$Threc[5]
  dis$THtoD <- sevenpathogens$Thd[5]
  dis$TIatoR <- sevenpathogens$Tay[5]
  dis$TIstoR <- sevenpathogens$Tsr[5]
  dis$TIstoH <- sevenpathogens$Tsh[5]
  dis$prob_symp <- sevenpathogens$ps[5]
  
  dis$ihr <- as.numeric(sevenpathogens[5,sapply(paste0('ihr',1:17),
                                                function(x)which(names(sevenpathogens)==x))])
  dis$ifr <- as.numeric(sevenpathogens[5,sapply(paste0('ifr',1:17),
                                                function(x)which(names(sevenpathogens)==x))])
  
  R0_to_beta <- function(dis) c(dis$R0, dis$R0/dis$CI)
  
  ## country data #################################
  
  ldata1 <- p2RandCountry(data, CD, 'UMIC', country_parameter_distributions, country_name=country)
  # combine country and disease parameters
  dis1 <- population_disease_parameters(data=ldata1, dis, R0betafun=R0_to_beta)
  
  
  ## complete p2, dis and data structs #################################
  
  ldata_dis_p2 <- p2Params(ldata1, dis1)
  # save:
  saveRDS(ldata_dis_p2,datafile)
  
}else{
  ldata_dis_p2 <- readRDS('data_file.Rds')
}

## load stored objects
dis2 <- ldata_dis_p2$dis
ldata <- ldata_dis_p2$data
p2 <- ldata_dis_p2$p2

## load econ models, which are written into file econ_models.R
source('econ_models.R')
# choose a model
econ = pc_model

# plot response function to epidemic for reference
epivars = seq(0,1e6,by=1000)
eplot = prop_to_consume(epivars,econ,ref_val=400000)
plotresponse <- ggplot() + 
  geom_line(aes(x=epivars/1e6,y=eplot$alpha1/econ$alpha1),size=2,colour='midnightblue') +
  theme_bw(base_size = 15) +
  labs(x='Million hospital cases',y='Relative propensity to consume') +
  scale_y_continuous(limits=c(0,1))
ggsave(plotresponse,filename='figures/response.png',width=5,height=5)

#############################################################

ref_vals = seq(100000,800000,by=50000)
baselines = c(.2,.5,.8)
outtab = data.frame(expand.grid(ref_vals,baselines,0,0))
colnames(outtab) = c('Transition point','Baseline','GDP loss','Deaths averted')

for(rv in 1:nrow(outtab)){
  
  ldata$ref_val = outtab[rv,1]
  ldata$baseline = outtab[rv,2]
  ## run model
  runlist <- p2SimVax(ldata, dis2, p2, econ)
  odevar <- runlist$integrated
  
  ## OUTPUT VARIABLES
  
  # time
  Tout <- odevar[,1]
  # numerical variables
  out_without_time <- odevar[,-1]
  # indices for disease compartments
  compindex <- ldata$compindex
  # number of disease states
  nStates = max(unlist(compindex))
  
  # indices for disease states
  S_index <- compindex$S_index
  E_index <- compindex$E_index
  I_index <- compindex$I_index
  H_index <- compindex$H_index
  D_index <- compindex$D_index
  R_index <- compindex$R_index
  
  # number of population strata in epi model (age groups plus sectors = 4+1 = 5)
  nStrata <- ldata$nStrata
  # number of time points
  nTime = nrow(odevar)
  # array of epi variables
  epi_vars_mat = array(out_without_time[,-c(1:econ$nEconODEs)],dim=c(nTime,nStrata,nStates))
  # individual epi states
  Ia <- epi_vars_mat[,, I_index[1]]
  Is <- epi_vars_mat[,, I_index[2]]
  Iout <- Ia+Is
  Hout <-  epi_vars_mat[,, H_index[1]]
  Dout <- epi_vars_mat[,, D_index[1]]
  Sout <- epi_vars_mat[, , S_index[1]]
  
  # counterfactual outcome (no integration between models)
  cntr <- runlist$counterfactual
  
  # plot outcomes
  Tout0 <- cntr[,1]
  Hout0 <-  array(cntr[,-c(1:(econ$nEconODEs+1))],dim=c(length(Tout0),nStrata,nStates))[,, H_index[1]]
  Dout0 <-  array(cntr[,-c(1:(econ$nEconODEs+1))],dim=c(length(Tout0),nStrata,nStates))[,, D_index[1]]
  Sout0 <-  array(cntr[,-c(1:(econ$nEconODEs+1))],dim=c(length(Tout0),nStrata,nStates))[,, S_index[1]]
  Rout0 <-  array(cntr[,-c(1:(econ$nEconODEs+1))],dim=c(length(Tout0),nStrata,nStates))[,, R_index[1]]
  
  plotdata = data.frame(Day=c(Tout,Tout0),
                        Hospitalised=c(rowSums(Hout),rowSums(Hout0)),
                        Consumption=c(odevar[,which(econ$econvarnames=='cons')+1],cntr[,which(econ$econvarnames=='cons')+1]),
                        Wealth=c(odevar[,which(econ$econvarnames==econ$wealth)+1],cntr[,which(econ$econvarnames==econ$wealth)+1]),
                        GDP=c(odevar[,which(econ$econvarnames=='cons')+1],cntr[,which(econ$econvarnames=='cons')+1])+econ$g,
                        Integrated=c(rep('Integrated model',length(Tout)),rep('Counterfactual',length(Tout0))))
  
  plotout <- ggplot(reshape2::melt(plotdata,id.var=c('Day','Integrated'))) + 
    geom_line(aes(x=Day,y=value,colour=Integrated),linewidth=2) +
    facet_wrap(~variable,scales = 'free_y') +
    theme_bw(base_size=15) +
    theme(legend.position = 'top') + 
    labs(y='',colour='')
  # print(plotout)
  
  # ggsave(plotout,filename=paste0('figures/pc_plot',ldata$ref_val,'-',ldata$baseline,'.png'),width=7,height=6)
  
  
  gdpscen = trapz(Tout,odevar[,which(econ$econvarnames=='cons')+1] + econ$g)
  # trapz(Tout0,cntr[,which(econ$econvarnames=='cons')+1] + econ$g)
  deathsscen = max(Dout)
  deathscnt = max(Dout0)
  
  outtab[rv,3:4] = c((gdp-gdpscen)/gdp*100,deathscnt-deathsscen)
}

(plotout <- ggplot(outtab) + 
  geom_line(aes(x=`Deaths averted`/1e3,y=`GDP loss`,colour=`Transition point`/1e3,group=Baseline),size=2.5) +
  theme_bw(base_size = 15) +
  labs(x='Deaths averted, thousands',y='GDP loss, % of GDP',colour='Transition point') +
  scale_colour_viridis(discrete=F, name="Transition point (thousand hospital cases)",option='inferno',direction = -1) +
  theme(legend.position = 'top') +
    annotate('text',label='Baseline = 0.8',x=13,y=1.5,angle=15)+
    annotate('text',label='Baseline = 0.5',x=25,y=4.1,angle=20)+
    annotate('text',label='Baseline = 0.2',x=27,y=7,angle=25)
  
  )
  

# ggsave(plotout,filename=paste0('figures/losses.png'),width=7,height=6)
