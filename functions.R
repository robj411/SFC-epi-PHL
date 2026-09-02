


simulate_epi_econ_model <- function(data, econ, integrate) {

  S0    <- data$NNs
  t0    <- data$tvec[1]
  tend  <- max(data$tvec)
  times <- seq(t0, tend, by = 1)

  imported <- 1 / sum(S0) * S0

  # parameters shared by all models
  shared <- list(
    alpha0        = econ$alpha0,
    alpha1        = econ$alpha1,
    alpha2        = econ$alpha2,
    theta         = econ$theta,
    G             = econ$G,
    cons0         = econ$cons0,
    emp0          = econ$emp0,
    lf            = econ$lf,
    q1            = data$q1,
    q2            = data$q2,
    integrate     = integrate,
    beta          = data$epidemic$beta,
    TEtoI         = data$epidemic$TEtoI,
    TItoR         = data$epidemic$TItoR,
    TItoC         = data$epidemic$TItoC,
    TCtoR         = data$epidemic$TCtoR,
    prob_detected = data$epidemic$prob_detected,
    NNs           = data$NNs,
    Mww           = data$contacts$Mww,
    Mcw           = data$contacts$Mcw,
    Mcc           = data$contacts$Mcc,
    Mcom          = data$contacts$Mcom,
    H_h_init      = econ$econ_init[[1]],
    S_init        = S0 - imported,
    E_init        = imported
  )

  # model-specific parameters and generator
  gen <- econ$gen
  if (econ$model_name == 'model1') {
    # gen         <- model1_gen
    model_extra <- list(lambda = econ$lambda)
  } else if (econ$model_name == 'model2') {
    # gen         <- model2_gen
    model_extra <- list(
      lambda_init = econ$econ_init[[2]],
      lambda_p0   = econ$lambda_p0,
      lambda_p1   = econ$lambda_p1
    )
  } else {
    stop(paste('Unknown model:', econ$model_name))
  }

  user_params <- c(shared, model_extra)

  mod <- gen$new(user = user_params)
  out <- mod$run(times)

  list(
    data           = data,
    trajectories     = out
  )
}


# get epi vars from ode output
get_epi_vars = function(mat, ldata, vars, nRemove = 0){
  
  # indices for disease compartments
  compindex <- ldata$compindex
  # number of disease states
  nStates = max(unlist(compindex))
  # number of population strata in epi model (age groups plus sectors = 4+1 = 5)
  nStrata <- ldata$nStrata
  # number of time points
  nTime = nrow(mat)
  
  # select only the epi state columns, ignoring any output() columns appended at the end
  epi_cols <- (nRemove + 1):(nRemove + nStrata * nStates)
  mat <- mat[, epi_cols]
  # array of epi variables
  epi_mat = array(mat,dim=c(nTime,nStrata,nStates))
  
  varlist = list()
  for(i in 1:length(vars)){
    firstletter = substr(vars[i],1,1)
    statename = paste0(firstletter, '_index')
    # determine which I if I
    statenumber = ifelse(vars[i]=='Id',2,1)
    statemat = epi_mat[,, compindex[[statename]][statenumber]]
    varlist[[vars[i]]] = rowSums(statemat)
  }
  varlist
}


get_results_df = function(runlist, epivars, econ, ldata, printflag=0){
  popsize = sum(ldata$Npop3)
  plotnames = c('Independent models','Integrated model')
  listnames = c('independent','trajectories')
  scen_df <- list()
  for(j in 1:2){
    mat = runlist[[listnames[j]]]
    integrate = j-1
    # time
    Tout <- mat[,1]
    # epi outcomes to plot
    epi_vars = get_epi_vars(mat=mat, ldata, vars=epivars, nRemove=econ$nEconODEs+1)
    
    # individual epi states
    Cout <-  epi_vars[[epivars]]
    
    cons_scen <- mat[, 'cons']
    gdp_scen  <- mat[, 'GDP']

    cons_list <- list(
      data.frame(Day = Tout, Integrated = 'Supply', variable = 'Consumption',
                 value = mat[, 'cons_s'] / cons_scen[1] * 100),
      data.frame(Day = Tout, Integrated = 'Demand', variable = 'Consumption',
                 value = mat[, 'cons_d'] / cons_scen[1] * 100)
    )
    cons_df <- do.call(rbind, cons_list)
                         
    econ_df = data.frame(Day = Tout,
               Consumption = cons_scen,
               GDP = gdp_scen,
               Integrated = plotnames[j])
    for(i in 1:length(econ$econvarlabels))
      econ_df[[econ$econvarlabels[i]]] = mat[,i+1]
    
    if('Productivity'%in%names(econ_df)){
      econ_df$Employment = with(econ_df, GDP/Productivity)/econ$lf * 100
    }else{
      econ_df$Employment = with(econ_df, GDP/econ$lambda)/econ$lf * 100
    }      
    
    scen_df[[j]] = cbind(econ_df, do.call(cbind, epi_vars))
    minemp = min(scen_df[[j]]$Employment)
    all_infected = Reduce('+',get_epi_vars(mat=mat, ldata, vars=c('C','Id','Iu'), nRemove=econ$nEconODEs+1))
    max_infected = max(all_infected)/1e6 # in millions
    recovered = get_epi_vars(mat=mat, ldata, vars='R', nRemove=econ$nEconODEs+1)
    cumulative_incidence_pc = max(recovered$R)/popsize*100
    if(printflag==1 | integrate==1)
      cat(paste0(plotnames[j],
        # ldata$q2,' & ',
        # econ$lambda_p1,' & ',
        # ldata$cc,' & ',
        ' & ',
        round(max_infected,1),' & ',
        scen_df[[j]]$Day[which.max(all_infected)],' & ',
        round(cumulative_incidence_pc),' & ',
        round(minemp,1),' & ',
        ifelse(j==1,'0 \\\\\n','')
      ))
  }
  # join, divide, and return
  alleconvars = c(econ$econvarlabels, 'Consumption', 'GDP')
  df = setDT(do.call(rbind,scen_df))
  # write as percent of counterfactual
  df[,(alleconvars) := lapply(.SD, function(x) x/x[1]*100),by=.(Day),.SDcols=alleconvars]
  
  ## get GDP loss
  cumulativegdp = df[,trapz(Day,GDP),by=Integrated]
  # GDP loss as percent of total (but could change to % of annual)
  gdplosspc = with(cumulativegdp, (V1[1]-V1[2])/V1[1])*100
  cat(paste0(round(gdplosspc,1),' \\\\\n '))
  
  colnames(df)[colnames(df)=='C'] <- 'Confirmed cases'
  
  melted_df = reshape2::melt(df,id.var=c('Day','Integrated'))
  melted_df = rbind(melted_df, cons_df)
  melted_df
}


plot_trajectories <- function(runlist, econ, ldata, printflag=0){
  
  ## OUTPUT VARIABLES
  
  plotdata = get_results_df(runlist, epivars='C', econ, ldata, printflag)
  plotdata$Integrated = factor(plotdata$Integrated, levels=c('Independent models','Integrated model','Supply','Demand'))
  
  plotdata
  
}


