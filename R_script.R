library(readxl)
library(pracma)
library(MASS)
library(abind)
library(Matrix)
library(haven)
library(ggplot2)
library(wbstats)
library(viridis)

# set up basics
country <- 'Philippines'
set.seed(0)
setwd(getSrcDirectory(function(){})[1])

# load functions
source('R_functions.R');

# read or create epi and country structures
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
  dis$R0 <- 1.5 # choosing something milder. sevenpathogens$R0[5]
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
  dis1 <- population_disease_parameters(data=ldata1, dis)
  
  
  ## complete p2, dis and data structs #################################
  
  ldata_dis_p2 <- p2Params(ldata1, dis1)
  # save:
  saveRDS(ldata_dis_p2,datafile)
  
}else{
  ldata_dis_p2 <- readRDS(datafile)
}

## load stored objects
dis2 <- ldata_dis_p2$dis
ldata <- ldata_dis_p2$data
p2 <- ldata_dis_p2$p2

## load econ models, which are written into file econ_models.R
source('econ_models.R')
# choose an econ model
econ = model1

# plot response function to epidemic for reference
epivars = seq(0,3e5,by=1000)
eplot = epi_to_econ(epivars,model1,ref_val=100000)
(plotresponse <- ggplot() + 
  geom_line(aes(x=epivars/1e3,y=eplot$alpha1/model1$alpha1),linewidth=2,colour='midnightblue') +
  theme_bw(base_size = 15) +
  labs(x='Thousand hospital cases',y='Relative propensity to consume') +
  scale_y_continuous(limits=c(0,1)) +
    geom_abline(slope=-1/85,intercept=1.925,colour='grey',linewidth=1.5,linetype=2) +
    geom_vline(xintercept=100,colour='grey',linewidth=1.5) +
    geom_hline(yintercept=.5,colour='grey',linewidth=1.5))
# ggsave(plotresponse,filename='figures/response.png',width=5,height=5)

## simulate ###########################################################

ref_vals = seq(25000,250000,by=25000)
baselines = c(.2,.5,.8)
outtab = data.frame(expand.grid(ref_vals,baselines,0,0))
colnames(outtab) = c('Transition point','Baseline','GDP loss','Deaths averted')
for(rv in 1:nrow(outtab))
{
  cat(paste0('Parameter combination ',rv,' out of ',nrow(outtab),'\n'))
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
  
  H = rowSums(Hout)
  Hc = rowSums(Hout0)
  
  as_pc <- function(scen,counter){
    100*c(scen/counter[1],counter/counter)
  }
  
  plotdata = data.frame(Day=c(Tout,Tout0),
                        Hospitalised=c(rowSums(Hout),rowSums(Hout0)),
                        Consumption=as_pc(econ$get_cons_from_out(odevar,econ,H,ldata,1),econ$get_cons_from_out(cntr,econ,Hc,ldata,0)),
                        Wealth=as_pc(odevar[,which(econ$econvarnames==econ$wealth)+1],cntr[,which(econ$econvarnames==econ$wealth)+1]),
                        GDP=as_pc(econ$get_gdp_from_out(odevar,econ,H,ldata,1),econ$get_gdp_from_out(cntr,econ,Hc,ldata,0)),
                        Integrated=c(rep('Integrated model',length(Tout)),rep('Counterfactual',length(Tout0))))
  colnames(plotdata)[3:5] <- c('Consumption, % of counterfactual','Wealth, % of counterfactual','GDP, % of counterfactual')
  plotout <- ggplot(reshape2::melt(plotdata,id.var=c('Day','Integrated'))) + 
    geom_line(aes(x=Day,y=value,colour=Integrated),linewidth=2) +
    facet_wrap(~variable,scales = 'free_y') +
    theme_bw(base_size=15) +
    theme(legend.position = 'top') + 
    labs(y='',colour='')
  print(plotout)
  
  ggsave(plotout,filename=paste0('figures/',econ$model_name,'_',ldata$ref_val,'-',ldata$baseline,'.png'),width=7,height=6)
  
  
  gdpscen = trapz(Tout,econ$get_gdp_from_out(odevar,econ,H,ldata,1))
  # trapz(Tout0,cntr[,which(econ$econvarnames=='cons')+1] + econ$g)
  deathsscen = max(Dout)
  deathscnt = max(Dout0)
  
  outtab[rv,3:4] = c((econ$gdp*2-gdpscen)/econ$gdp/2*100,deathscnt-deathsscen)
}



## plot results ###############################

mn = which(c('model1','model2','model3','PC','modelpc2','modelpc3')==econ$model_name)
xs = matrix(c(13.9,19.4,20.1, 15,19.5,20., 15,21,21.2,
              13.9,21.4,21.8, 9,21,21, 9,20,21),nrow=3)
ys = matrix(c(1.35,6,8.3, .95,3.4,4.8, 1.8,5.4,7,
              1.6,6.5,8.3, .45,3.7,4.6, 1,5.9,9.),nrow=3)
angles = matrix(c(18,57,59, 20,55,58, 22,50,52,
                  18,52,55, 11,51,53, 13,48,50),nrow=3)

(plotout <- ggplot(outtab) + 
    geom_line(aes(x=`Deaths averted`/1e3,y=`GDP loss`,colour=`Transition point`/1e3,group=Baseline),linewidth=2.5) +
    theme_bw(base_size = 15) +
    labs(x='Deaths averted, thousands',y='GDP loss, % of GDP',colour='Transition point') +
    scale_colour_viridis(discrete=F, name="Transition point (thousand hospital cases)",option='inferno',direction = -1) +
    theme(legend.position = 'top') )#+
    # annotate('text',label='Baseline = 0.8',x=xs[1,mn],y=ys[1,mn],angle=angles[1,mn])+
    # annotate('text',label='Baseline = 0.5',x=xs[2,mn],y=ys[2,mn],angle=angles[2,mn])+
    # annotate('text',label='Baseline = 0.2',x=xs[3,mn],y=ys[3,mn],angle=angles[3,mn]))
  
ggsave(plotout,filename=paste0('figures/losses_',econ$model_name,'.png'),width=7,height=6)


