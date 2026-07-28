library(pracma)
library(MASS)
library(abind)
library(Matrix)
library(haven)
library(ggplot2)
library(viridis)
library(latex2exp)
library(data.table)
library(odin)

# set up basics

set.seed(0)
setwd(getSrcDirectory(function(){})[1])

# load functions
source('functions.R')
source('data_functions.R')

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
econ = model1
odinfile = paste0("odin_",econ$model_name,".R")
econ$gen = odin::odin(odinfile)

## simulate ###########################################################

nscen = 9
parameter_combinations = data.frame(q1 = rep(0.862752, nscen), 
                                    q2 = rep(0.000237, nscen),
                                    lambdap = rep(0.5, nscen),
                                    cc = rep(0.4,nscen))
parameter_combinations$q1[2] = 0.73
parameter_combinations$q1[3] = 0.93
parameter_combinations$q2[4] = 0.0000001
parameter_combinations$q2[5] = 0.001
parameter_combinations$lambdap[6] = 0.25
parameter_combinations$lambdap[7] = 0.75
parameter_combinations$cc[8] = 0.2
parameter_combinations$cc[9] = 0.6
sensvars = unique(colnames(parameter_combinations))
sensvarnames = c('q_1','q_2','\\gamma_1','\\eta')
plotdatalist = list()
for(rv in 1:nrow(parameter_combinations)){
  # cat(paste0('Parameter combination ',rv,' out of ',nrow(parameter_combinations),'\n'))
  econ$lambda_p1 = parameter_combinations$lambdap[rv]
  ldata$q1 = parameter_combinations$q1[rv]
  ldata$q2 = parameter_combinations$q2[rv]
  ldata$cc = parameter_combinations$cc[rv]
  ldata = decompose_contacts(ldata, consumption_contact = parameter_combinations$cc[rv])
  # print(ldata$epidemic$beta)
  if(rv==1) independent = simulate_epi_econ_model(data = ldata, econ = econ, integrate = 0)$trajectories
  
  ## run model
  runlist <- simulate_epi_econ_model(data = ldata, econ = econ, integrate = 1)
  runlist$independent = independent
  
  plotdata = plot_trajectories(runlist, econ, ldata, printflag=(rv==1))
  
  if(econ$model_name == 'model1'){
    label = TeX(paste0(econ$model_name,"; $q_2$ = ", ldata$q2))
  }else{
    label = TeX(paste0(#econ$model_name,'; ',
                       "$q_1 = ", signif(ldata$q1,3),
                       "$; $q_2 = ", signif(ldata$q2,3),
                       "$; $\\gamma_1 = ", signif(econ$lambda_p1,3),
                       "$; $\\eta = ", signif(parameter_combinations$cc[rv],3),"$"))
  }
  supplydemandcolour = 'midnightblue'
  indecolour = 'darkslategray3'
  if(rv==1){
    plotout <- ggplot(plotdata) +
      geom_point(data=data.frame(x=100,y=100),aes(x=x,y=y),colour='white') + 
      geom_line(aes(x=Day,y=value,colour=Integrated,linewidth=Integrated,linetype=Integrated)) +
      scale_linewidth_manual(values = c(`Independent models`=2,Demand=0.5,Supply=1,`Integrated model`=2),guide = 'none') +
      scale_linetype_manual(values = c(`Independent models`=1,Demand=1,Supply=3,`Integrated model`=1),guide = 'none') +
      scale_colour_manual(values = c(`Independent models`=indecolour,Demand=supplydemandcolour,Supply=supplydemandcolour,`Integrated model`='hotpink')) +
      facet_wrap(~variable,scales = 'free_y',nrow=1) +
      theme_bw(base_size=15) +
      theme(legend.position = 'top', 
            strip.text.y = element_blank(), 
            strip.background = element_blank(),
            # legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(-0,-10,-20,-10)) +
      guides(color = guide_legend(override.aes = list(linewidth = c(2,2,1,1),
                                                      linetype = c(1,1,3,1)) ) ) +
      labs(y='',colour='',x='Day') 
    
    # plotout = plotout + labs(title= label)
    ggsave(plotout,filename=paste0('figures/',econ$model_name,'_scenario_',rv,
                                   # econ$lambda_p1,'-',
                                   # ldata$q1,'-',
                                   # ldata$q2,'-',
                                   # parameter_combinations$cc[rv],
                                   '.png'),width=12,height=3)
  }else{
    sensind = floor(rv/2)
    plotdata$sensvar = parameter_combinations[[sensind]][rv]
    plotdatalist[[rv]] = plotdata
    if(rv%%2==1){
      l1 = paste0('$',sensvarnames[sensind],' = ',parameter_combinations[[sensind]][rv],'$')
      lm1 = paste0('$',sensvarnames[sensind],' = ',parameter_combinations[[sensind]][rv-1],'$')
      plotdata = rbind(plotdatalist[[rv]],plotdatalist[[rv-1]])
      appender <- function(string)  TeX(paste("$",sensvarnames[sensind],"\\,$ = ", string))
      # plotdata$variable <- fct_recode(plotdata$variable, "Confirmed cases" = "Confirmed cases (millions)")
      # plyr::revalue(plotdata$variable, c("Confirmed cases" = "Confirmed cases (millions)"))
      levels(plotdata$variable)[levels(plotdata$variable)=="Confirmed cases"] <- "Confirmed cases (millions)"
      
      plotdata$value[plotdata$variable=='Confirmed cases (millions)'] = plotdata$value[plotdata$variable=='Confirmed cases (millions)']/1e6
      plotout <- ggplot(plotdata) +
        # geom_point(data=data.frame(x=100,y=100),aes(x=x,y=y),colour='white') + 
        geom_line(aes(x=Day,y=value,colour=Integrated,linewidth=Integrated,linetype=Integrated)) +
        scale_linewidth_manual(values = c(`Independent models`=2,Demand=0.5,Supply=1,`Integrated model`=2),guide = 'none') +
        scale_linetype_manual(values = c(`Independent models`=1,Demand=1,Supply=3,`Integrated model`=1),guide = 'none') +
        scale_colour_manual(values = c(`Independent models`=indecolour,Demand=supplydemandcolour,Supply=supplydemandcolour,`Integrated model`='hotpink')) +
        facet_grid(variable~sensvar,labeller = labeller(
          sensvar = as_labeller(appender, default = label_parsed),  # transform + parse
          variable     =  label_wrap_gen(width = 10)                                         # no parsing
        ),scales = 'free_y',
                   switch='y') +
        theme_bw(base_size=15) +
        theme(strip.text.y.left = element_text(angle = 0),
              strip.switch.pad.grid = unit(-.5, "cm"),
              legend.position = 'top', 
              # strip.text.y = element_blank(),
              strip.background = element_blank(),
              strip.placement = 'outside',
              legend.margin=margin(0,100,0,0),
              legend.box.margin=margin(0,-10,-10,-10)) +
        guides(color = guide_legend(nrow=2,byrow=F,
                                    override.aes = list(linewidth = c(2,2,1,1),
                                                        linetype = c(1,1,3,1)) ) ) +
        labs(y='',colour='',x='Day') 
      
      print(plotout)
      # plotout = plotout + labs(title= label)
      ggsave(plotout,filename=paste0('figures/',econ$model_name,'_scenario_',rv,
                                     # econ$lambda_p1,'-',
                                     # ldata$q1,'-',
                                     # ldata$q2,'-',
                                     # parameter_combinations$cc[rv],
                                     '.png'),width=4.75,height=8)
      
    }
  }
}


