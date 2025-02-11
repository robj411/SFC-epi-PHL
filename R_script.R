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
setwd(getSrcDirectory(function(){})[1])
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
econ = model2
econ = model3

# plot response function to epidemic for reference
epivars = seq(0,3e5,by=1000)
eplot = fear_of_infection(epivars,pc_model,ref_val=100000)
plotresponse <- ggplot() + 
  geom_line(aes(x=epivars/1e3,y=eplot$alpha1/pc_model$alpha1),linewidth=2,colour='midnightblue') +
  theme_bw(base_size = 15) +
  labs(x='Thousand hospital cases',y='Relative propensity to consume') +
  scale_y_continuous(limits=c(0,1))
ggsave(plotresponse,filename='figures/response.png',width=5,height=5)

## simulate ###########################################################

ref_vals = seq(50000,250000,by=25000)
baselines = c(.2,.5,.8)
outtab = data.frame(expand.grid(ref_vals,baselines,0,0))
colnames(outtab) = c('Transition point','Baseline','GDP loss','Deaths averted')
rv=18
for(rv in 1:nrow(outtab))
  {
  
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
  
  as_pc <- function(scen,counter){
    100*c(scen/counter[1],counter/counter)
  }
  
  plotdata = data.frame(Day=c(Tout,Tout0),
                        Hospitalised=c(rowSums(Hout),rowSums(Hout0)),
                        Consumption=as_pc(odevar[,which(econ$econvarnames=='cons')+1],cntr[,which(econ$econvarnames=='cons')+1]),
                        Wealth=as_pc(odevar[,which(econ$econvarnames==econ$wealth)+1],cntr[,which(econ$econvarnames==econ$wealth)+1]),
                        GDP=as_pc(odevar[,which(econ$econvarnames=='cons')+1]+econ$g,cntr[,which(econ$econvarnames=='cons')+1]+econ$g),
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
  
  
  gdpscen = trapz(Tout,odevar[,which(econ$econvarnames=='cons')+1] + econ$g)
  # trapz(Tout0,cntr[,which(econ$econvarnames=='cons')+1] + econ$g)
  deathsscen = max(Dout)
  deathscnt = max(Dout0)
  
  outtab[rv,3:4] = c((gdp*2-gdpscen)/gdp/2*100,deathscnt-deathsscen)
}

## plot results ###############################

mn = which(c('PC','model2','model3')==econ$model_name)
xs = matrix(c(14,21.5,21.8, 10,21,21, 9,20,21),nrow=3)
ys = matrix(c(1.5,6.5,8.3, .5,3.7,4.6, .9,5.9,9.),nrow=3)
angles = matrix(c(15,44,44, 10,45,45, 8,38,39),nrow=3)

(plotout <- ggplot(outtab) + 
    geom_line(aes(x=`Deaths averted`/1e3,y=`GDP loss`,colour=`Transition point`/1e3,group=Baseline),linewidth=2.5) +
    theme_bw(base_size = 15) +
    labs(x='Deaths averted, thousands',y='GDP loss, % of GDP',colour='Transition point') +
    scale_colour_viridis(discrete=F, name="Transition point (thousand hospital cases)",option='inferno',direction = -1) +
    theme(legend.position = 'top') +
    annotate('text',label='Baseline = 0.8',x=xs[1,mn],y=ys[1,mn],angle=angles[1,mn])+
    annotate('text',label='Baseline = 0.5',x=xs[2,mn],y=ys[2,mn],angle=angles[2,mn])+
    annotate('text',label='Baseline = 0.2',x=xs[3,mn],y=ys[3,mn],angle=angles[3,mn]))
  
ggsave(plotout,filename=paste0('figures/losses_',econ$model_name,'.png'),width=7,height=6)


## plot contacts ###########################

matrices = c('basic_contact_matrix','worker_worker_mat','community_worker_mat','CM_4','school_mat','hospitality_contacts')
matnames = c('Total','Workers','Community—worker','Base','School','Hospitality')
df = do.call(rbind,lapply(1:length(matrices),function(x)reshape2::melt(ldata$contacts[[matrices[x]]])))
colnames(df) = c('to','from','contacts')
df$matrix = rep(matnames,times=sapply(matrices,function(x)prod(dim(ldata$contacts[[x]]))))
fourbyfour <- c('School','Base','Hospitality')
df$to[df$matrix%in%fourbyfour] = df$to[df$matrix%in%fourbyfour] + 1 
df$from[df$matrix%in%fourbyfour] = df$from[df$matrix%in%fourbyfour] + 1 

groups = c('Workforce in place','Aged 0 to 4','Aged 5 to 14','Aged 15 to 69','Over 69')
cntcts = ggplot(df, aes(x=from, y=to, fill = contacts)) +
  geom_tile() +
  facet_wrap(~factor(matrix,levels=matnames), strip.position="bottom") + 
  scale_fill_viridis() +  
  theme_bw(base_size=18) +
  labs(title = "") +
  scale_x_continuous(name='',breaks=1:5,labels=groups,expand=c(0,0),position = "top")+
  scale_y_continuous(name='',breaks=1:5,labels=groups,expand=c(0,0),trans = "reverse")+
  theme(
    panel.margin=unit(.5, "lines"),
    panel.border = element_rect(color = "white", fill = NA), 
    strip.background = element_rect(fill="white",colour='white'),
    axis.text.x = element_text(angle=45,hjust=0,vjust=1),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank(),
    plot.title = element_text(hjust = 1),        
    legend.position="right",
    panel.grid.major = element_line(colour="white"),
    panel.grid.minor = element_line(colour="white")
  ) 
ggsave(cntcts,filename='figures/contacts.png',width=10,height=6.5)



# library(ggfortify)
library(ggdist)
library(gridExtra)

nSamples = 200000
qs <- c(.05,.5,.95)
df = data.frame(baseline = rbeta(nSamples,5,5),transition = rgamma(nSamples,2,.00001),
                work_frac = rbeta(nSamples,10,20), prop_to_consume = rbeta(nSamples,10,10))

p1 = ggplot(data=df,aes(x=baseline)) +
  stat_slab( aes(thickness = after_stat(ifelse(.width <= 0.9, pdf, NA))),fill = "gray85", .width = .9) +
  stat_slab(colour='midnightblue',linewidth=2,fill=NA) +
  theme_bw(base_size=15) +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks=quantile(df$baseline,qs),labels=round(quantile(df$baseline,qs),2)) +
  labs(   x = "Baseline behaviour",y = NULL)

p2 = ggplot(data=df,aes(x=transition)) +
  stat_slab( aes(thickness = after_stat(ifelse(.width <= 0.9, pdf, NA))),fill = "gray85", .width = .9) +
  stat_slab(colour='midnightblue',linewidth=2,fill=NA) +
  theme_bw(base_size=15) +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks=quantile(df$transition,qs),labels=round(quantile(df$transition,qs)/1000)) +
  labs(   x = "Transition point, thousand hospital cases",y = NULL)

p3 = ggplot(data=df,aes(x=work_frac)) +
  stat_slab( aes(thickness = after_stat(ifelse(.width <= 0.9, pdf, NA))),fill = "gray85", .width = .9) +
  stat_slab(colour='midnightblue',linewidth=2,fill=NA) +
  theme_bw(base_size=15) +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks=quantile(df$work_frac,qs),labels=round(quantile(df$work_frac,qs),2)) +
  labs(   x = "Fraction contacts from work",y = NULL) 

p4 = ggplot(data=df,aes(x=prop_to_consume)) +
  stat_slab( aes(thickness = after_stat(ifelse(.width <= 0.9, pdf, NA))),fill = "gray85", .width = .9) +
  stat_slab(colour='midnightblue',linewidth=2,fill=NA) +
  theme_bw(base_size=15) +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(breaks=quantile(df$prop_to_consume,qs),labels=round(quantile(df$prop_to_consume,qs),2)) +
  labs(   x = "Propensity to consume",y = NULL) 

params <- grid.arrange(p1,p2,p3,p4,nrow=2)
ggsave(params,file='figures/parameters.png')
