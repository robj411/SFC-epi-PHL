
# data: list of general model parameters

data_start <- function(country_name = 'Philippines', iso3 = 'PHL') {
  
  data <- list()
  
  data$iso3 = iso3
  data$country_name = country_name
  data$nSectors <- 1
  data$tvec <- c(0, 300)
  
  data$ageindex <- list(1:3, 4:13, 14:21)
  data$adInd <- 2
  
  compindex <- list()
  
  compindex$S_index <- c(1) 
  compindex$E_index <- c(2)
  compindex$I_index <- c(3:4) # Iu, Id
  compindex$C_index <- c(5)
  compindex$R_index <- c(6)
  data$compindex <- compindex
  
  return(data)
}

# curate country data
#
# data: list of country-specific and general model parameters

gather_country_data <- function(data) {
  
  ref_year = 2023
  iso3 = data$iso3
  country_name = data$country_name
  
  # population
  # population by age
  data(popAge5dt)   # population by 5-year age group (long format)
  countrydata = subset(popAge5dt, name == country_name)
  maxyear = max(countrydata$year)
  pop_rec = subset(countrydata, year == maxyear, select = c(age, pop))
  Npop = pop_rec$pop*1000
  data$Npop <- Npop
  ageindex = data$ageindex
  Npop3 <- sapply(ageindex, function(x) sum(Npop[x]))
  data$Npop3 <- Npop3
  
  # employment rate
  dat <- get_ilostat(id = 'EAP_DWAP_SEX_AGE_RT_A', segment = 'indicator') 
  fiveyr <- subset(dat,ref_area==iso3 & sex=='SEX_T' & time==ref_year & grepl('5YR',classif1))
  fiveyr$pop <- c(sum(Npop[4:21]),Npop[4:13],sum(Npop[14:21]))
  
  # this was to add people aged over 65
  # shrinkage = with(subset(fiveyr,classif1%in%c('AGE_5YRBANDS_Y60-64','AGE_5YRBANDS_Y55-59','AGE_5YRBANDS_Y50-54')),obs_value[2:3]/obs_value[1:2])
  # shrinkagerate = shrinkage[2]/shrinkage[1]
  # lastval = fiveyr$obs_value[13-2]
  # workers = c()
  # for(i in 14:21){
  #   thispop = Npop[i]
  #   newval = lastval*shrinkage[2]*shrinkagerate^(i-14)
  #   workers[i-13] = newval*thispop/100
  #   lastval = newval
  # }
  # retiredage = workers/Npop[14:21]*100
  # percentages <- c(subset(fiveyr,!classif1%in%c('AGE_5YRBANDS_TOTAL','AGE_5YRBANDS_YGE65'))$obs_value, retiredage[1])
  
  adind = data$adInd
  percentages <- subset(fiveyr,!classif1%in%c('AGE_5YRBANDS_TOTAL','AGE_5YRBANDS_YGE65'))$obs_value
  populations <- Npop[ageindex[[adind]]]
  lfpr <- sum(percentages*populations/100)/sum(populations)
  data$employmentrate = lfpr
  
  workers_by_sector <- data$employmentrate*Npop3[adind] # sum(allsectors)
  # put into order: workers, then children, non-workers, and retired
  NNs <- as.numeric(c(workers_by_sector,
                      Npop3[1:(adind-1)],
                      Npop3[adind]-sum(workers_by_sector),
                      Npop3[(adind+1):length(Npop3)]))
  
  data$NNs <- NNs
  data$nStrata <- length(data$NNs)
  
  # generate basic contact components
  data$contacts <- get_basic_contacts(data)
  basic_contact_matrix <- get_scaled_contacts(data, data$NNs)
  data$contacts$basic_contact_matrix <- basic_contact_matrix
  
  # econ data
  # using wb data
  gdpdata <- wb_data("NY.GDP.MKTP.CN",country = country_name, start_date = 2018, end_date = 2024)
  tdata <- wb_data("GC.TAX.TOTL.CN",country = country_name, start_date = 2018, end_date = 2024)
  data$tax = c(subset(tdata,date==ref_year)$GC.TAX.TOTL.CN/1e12)
  data$gdp = c(subset(gdpdata,date==ref_year)$NY.GDP.MKTP.CN/1e12)
  
  return(data)
}

collapse_cm <- function(contact_matrix, age_groups16, age_counts16) {
  # collapse along 2nd dimension: sum all contacts made
  Ci = matrix(0, nrow = length(age_counts16), ncol = length(age_groups16))
  for(i in 1:ncol(Ci))
    Ci[,i] = rowSums(contact_matrix[, age_groups16[[i]]])
  
  # collapse along 1st dimension: weighted average over group making contacts
  cm = matrix(0, nrow = length(age_groups16), ncol = length(age_groups16))
  for(i in 1:nrow(cm)){
    indices = age_groups16[[i]]
    weights = age_counts16[indices]
    cm[i,] = colSums(weights * Ci[indices,]) / sum(weights) 
  }
  
  return(cm)
}

transpose_and_normalise_conmats = function(cms, age_counts16){
  # normalise matrix: each contact made should be reciprocated
  norm_cm = list()
  for(cmi in 1:length(cms)){
    cm = t(cms[[cmi]])
    contacts_made_subject = rowSums(cm)
    reciprocal_matrix = matrix(0, nrow = nrow(cm), ncol = ncol(cm))
    for(i in 1:nrow(cm)){
      for(j in 1:ncol(cm)){
        total_ij_contacts = cm[i,j] * age_counts16[i]
        contacts_per_j = total_ij_contacts / age_counts16[j]
        reciprocal_matrix[i,j] = contacts_per_j
      }
    }
    contacts_made_object = colSums(reciprocal_matrix)
    # print(names(cms)[cmi])
    # print(contacts_made_object-contacts_made_subject)
    # print(sum(abs(contacts_made_object-contacts_made_subject)))
    norm_cm[[names(cms)[cmi]]] = (cm + t(reciprocal_matrix)) / 2
  }
  norm_cm
}



four_to_three = function(cm){
  threebyfour = matrix(0, nrow = 3, ncol = 4)
  threebyfour[1,] = cm[2,]
  threebyfour[2,] = p*cm[1,] + q*cm[3,]
  threebyfour[3,] = cm[4,]
  newmat = matrix(0, nrow = 3, ncol = 3)
  newmat[,1] = threebyfour[,2]
  newmat[,2] = threebyfour[,1] + threebyfour[,3]
  newmat[,3] = threebyfour[,4]
  newmat
  # aleq = all.equal(newmat, reduced_cms$all)
}

# take population contact matrix and decompose into items that will be
# impacted differently by configurations
#
# data: list of general and country-specific model parameters
#
# contacts: list of contact parameters

get_basic_contacts <- function(data) {
  
  Npop <- data$Npop
  NN <- data$NNs
  ageindex = data$ageindex
  
  polymod_survey_data = get_polymod_population()
  polymod_survey_data$population <- Npop
  exp_contact = extrapolate_polymod(
    population = polymod_survey_data
  )
  
  # resize age groups (note that conmat matrices are transposed:
  # "The contact matrices created using this package are transposed when compared to the contact matrices discussed by Prem and Mossong. That is, the rows are “age group to”, and the columns are “age group from”."
  # we want the subject (contact from) on the rows and object (contact to) on the columns
  cm = t(exp_contact[['all']])
  maxindex = nrow(cm)
  age_counts16 = Npop[1:maxindex]
  age_counts16[maxindex] = sum(Npop[maxindex:length(Npop)])
  
  age_groups16 = ageindex
  age_groups16[[length(ageindex)]] = c(ageindex[[length(ageindex)]])[ageindex[[length(ageindex)]]<=16]
  
  normalised_cms = transpose_and_normalise_conmats(exp_contact, age_counts16)
  
  CM_16 <- unname(normalised_cms$all)
  s16 <- nrow(CM_16)
  
  reduced_cms = lapply(normalised_cms, function(x)  collapse_cm(contact_matrix = x, age_groups16, age_counts16))
  
  worker_index = data$adInd
  worker_age_contacts = sum(reduced_cms$all[worker_index,])
  worker_age_self_contacts = reduced_cms$all[worker_index,worker_index]
  frac_work_contacts = sum(reduced_cms$work[worker_index,]) / worker_age_contacts
  
  # workforce participation rate
  # LFPR 15–64 (total), Philippines
  lfpr_1564 = WDI(country=data$iso3, indicator="SL.TLF.ACTI.ZS", start=2010, end=2026)
  # latest observation
  latest = lfpr_1564[order(lfpr_1564$year, decreasing=TRUE), ][1, ]
  p = latest$SL.TLF.ACTI.ZS/100
  q = 1-p
  
  popsizes = sapply(age_groups16, function(x) sum(age_counts16[x]))
  popfracs = popsizes/sum(popsizes)
  
  # population size vector with working-age people split between two groups
  N_4 = c(p*popsizes[worker_index], popsizes)
  N_4[1+worker_index] = q*popsizes[worker_index]
  
  # x=\frac{pN_1+qN_3}{N_1}
  x = (p*N_4[1] + q*N_4[3])/N_4[1]
  
  # fractions of workplace contacts made with others by workers
  # need to use conmat, or some other source
  k_vec = reduced_cms$work[worker_index,] / sum(reduced_cms$work[worker_index,])
  # the fraction of adults making contact with workers as customers 
  kstar = popfracs[worker_index] * sum(k_vec[-worker_index]) / sum(popfracs[-worker_index])
  # the fraction of workplace contacts made with other workers
  fWW = k_vec[2] - kstar
  
  
  # the total number of workplace contacts between workers
  WW = fWW * frac_work_contacts * worker_age_contacts / p
  CW3 = kstar / (fWW * (x/2 +q)) * WW
  CW2_plus_CW4 = WW / fWW - WW - (x/2 - p + 1)*CW3 
  CW2 = CW2_plus_CW4 * k_vec[1] / (k_vec[1] + k_vec[3])
  CW4 = CW2_plus_CW4 * k_vec[3] / (k_vec[1] + k_vec[3])
  CC2 = reduced_cms$all[2,1] - p*CW2
  CC4 = reduced_cms$all[2,3] - p*CW4
  CC3 = reduced_cms$all[2,2] - p*WW - p*CW3 * (2*q + x)
  
  Mww <- Mcw <- Mcom <- Mcc <- matrix(0, nrow = 4, ncol = 4)
  Mww[1,1] = WW
  
  Mcw[1,] = c(x*CW3, CW2, q*CW3, CW4)
  Mcw[2:4,1] = c(N_4[1] / N_4[2] * CW2, p*CW3, N_4[1] / N_4[4] * CW4)
  
  Mcom[2:4, 2:4] = reduced_cms$all
  Mcom[1,] = c(p*CC3, CC2, q*CC3, CC4)
  Mcom[3,] = c(p*CC3, CC2, q*CC3, CC4)
  Mcom[2:4,1] = c(N_4[1] / N_4[2] * CC2, p*CC3, N_4[1] / N_4[4] * CC4)
  Mcom[2:4,3] = c(N_4[3] / N_4[2] * CC2, q*CC3, N_4[3] / N_4[4] * CC4)
  
  Mcc[2:4, 2:4] = reduced_cms$other
  Mcc[1,] = Mcc[3,] 
  Mcc[3,] = Mcc[3,] 
  Mcc[,1] = Mcc[,3] * p 
  Mcc[,3] = Mcc[,3] * q
  # four_to_three(Mcc)
  
  # subtract consumer--consumer contacts 
  Mcom = Mcom - Mcc
  
  Mtot = Mww + Mcw + Mcom + Mcc
  ## save
  
  contacts <- list()
  contacts$Mww <- Mww
  contacts$Mcw <- Mcw
  contacts$Mcc <- Mcc
  contacts$Mcom <- Mcom
  
  data$contacts = contacts
  # CM_4_orig - collapse_cm(contact_matrix, NN)
  
  return(contacts)
  
}



# construct contact matrix from component matrices
#
# data: list of general model parameters
# NN: population by stratum
# relative_consumption: consumption relative to counterfactual / normal times
# relative_work: number of workers relative to counterfactual / normal times

get_scaled_contacts <- function(data, NN, relative_consumption=1, relative_work=1) {
  
  ## variables to use
  contacts <- data$contacts
  Mcom <- contacts$Mcom
  Mww <- contacts$Mww
  Mcw <- contacts$Mcw
  Mcc <- contacts$Mcc
  
  contact_matrix = Mcom + 
    relative_work^2 * Mww +
    relative_work*relative_consumption * Mcw +
    relative_consumption^2 * Mcc
  
  return(contact_matrix)
  
}


## get candidate infectees (step before R0)
#
# dis: list of pathogen parameters
# data: list of general and country-specific model parameters
#
# CI: candidate infectees (R0 (reproduction number) / beta)

get_candidate_infectees <- function(dis, data) {
  
  
  S <- N <- data$NNs
  contact_matrix = data$contacts$basic_contact_matrix
  nStrata = length(N)
  
  # Rates
  prob_detected = dis$prob_detected
  prob_isolated = dis$prob_isolated
  sig1 <- prob_detected / dis$TEtoI
  sig2 <- (1 - prob_detected) / dis$TEtoI
  g1 <- 1 / dis$TItoR
  g2 <- 1 / dis$TItoC
  confirmed <- 1 / dis$TCtoR
  
  FOIin <- contact_matrix * t(pracma::repmat(S, nStrata, 1)) / pracma::repmat(N, nStrata, 1)
  
  # force of infection from groups E (none), Iu (undetected), Id (not detected yet), C (detected)
  Fmat <- matrix(0, 4 * nStrata, 4 * nStrata)
  Fmat[1:nStrata, (nStrata + 1):(4*nStrata)] <- cbind(FOIin, FOIin, prob_isolated*FOIin)
  
  # rate of exit from states
  ones <- matrix(1,nStrata,1)
  vvec <- c((sig1 + sig2)*ones, g1*ones, g2*ones, confirmed*ones)
  
  n <- length(vvec)
  V <- diag(vvec)
  nmat <- diag(rep(1,nStrata));
  V[(nStrata + 1):(2 * nStrata), 1:nStrata] <- -sig1*nmat
  V[(2 * nStrata + 1):(3 * nStrata), 1:nStrata] <- -sig2*nmat
  V[(3 * nStrata + 1):(4 * nStrata), (2 * nStrata + 1):(3 * nStrata)] <- -g2*nmat
  
  NGM <- Fmat %*% Matrix::solve(V)
  ev <- eigen(NGM,only.values = T,symmetric = F) #largest in magnitude (+/-) 
  d <- ev$values
  CI <- max(Re(d))
  
  return(CI)
}


simulate_epi_econ_model <- function(data, dis, econ) {
  
  ## PARAMETERS
  nStrata <- length(data$NNs)
  nSectors <- data$nSectors
  S0 <- data$NNs
  tend <- max(data$tvec)
  compindex <- data$compindex
  nStates <- max(unlist(compindex))
  
  ## copy over only required items
  
  rundata <- list(
    adInd = data$adInd,
    NNs = data$NNs,
    q1 = data$q1,
    q2 = data$q2,
    compindex = data$compindex,
    contacts = data$contacts
  )
  
  # initial conditions
  
  imported <- 100 / sum(S0) * S0
  t0 <- data$tvec[1]
  epi_init_mat <- matrix(0, nrow = nStrata, ncol = nStates)
  epi_init_mat[, compindex$S_index[1]] <- S0 - imported
  epi_init_mat[, compindex$E_index[1]] <- imported
  
  y0 <- c(econ$econ_init, as.vector(epi_init_mat))
  
  # run for counterfactual
  
  econ$integrate <- 0
  
  fun <- function(t, y, p) ODEs(data, t, dis, y, econ)
  tmpout <- deSolve::ode(times = seq(t0,tend,1), y = y0, func = fun,
                         parms=list(data=rundata, nStrata=nStrata, dis=dis, econ=econ),
                         method='impAdams_d')
  # plot(rowSums(tmpout[,12:23]))
  # store counterfactual values
  econ$counter_time <- tmpout[,1]
  linking_values = t(apply(tmpout[,-1],1,function(y){ econ = econ$cons_link_fun(y, econ);
    unlist(econ$econ_to_epi(y,econ))}))
  econ$counter_cons <- linking_values[,1]
  econ$counter_worker <- linking_values[,2]
  econ$integrate <- 1
  
  # solve ODEs
  out <- deSolve::ode(times = seq(t0, tend, by=1), y = y0, func = fun,
                      parms=list(data=rundata, nStrata=nStrata, dis=dis, econ=econ),
                      method='impAdams_d')
  
  ## OUTPUTS:  
  returnobject <- list(
    data=data,
    integrated=out,
    counterfactual=tmpout
  )
  
  return(returnobject)
}


epi_to_econ = function(epi_var,econ, q1 = 0.86, q2 = .0001){
  
  # make sure we start at one
  # ref_scalar = 1/(1+exp(-(ref_val)/gradient))
  # scalar = 1/(1+exp(-(ref_val-epi_var)/gradient)) / ref_scalar
  # prop_to = baseline + (1-baseline) * scalar
  
  prop_to = q1 + (1-q1) / (1+q2*epi_var)
  
  econ$scalar = prop_to
  return(econ)
  
}

epi_to_econ_deriv = function(epivar, dot_epivar, q1 = 0.86, q2 = .0001){
  deriv = (-(1-q1) * q2 * dot_epivar )/((1+q2*epivar )^2)
  deriv
}

ODEs <- function(data, t, dis, y, econ) {
  
  ## BLOCK 0: ACCESS VARIABLES ####################
  
  ## variables
  
  nEconODEs = econ$nEconODEs
  integrate = econ$integrate
  
  compindex <- data$compindex
  S_index <- compindex$S_index
  E_index <- compindex$E_index
  I_index <- compindex$I_index
  C_index <- compindex$C_index
  R_index <- compindex$R_index
  
  NN0 <- data$NNs
  nStrata <- length(NN0)
  nStates = max(unlist(compindex))
  
  epi_vars_mat <- matrix(y[-c(1:nEconODEs)], nrow = nStrata, ncol = nStates)
  
  # indices <- 1:nStrata
  S <- epi_vars_mat[,S_index[1]]
  E <- epi_vars_mat[,E_index[1]]
  Iu <- epi_vars_mat[,I_index[1]]
  Id <- epi_vars_mat[,I_index[2]]
  C <- epi_vars_mat[,C_index[1]]
  R <- epi_vars_mat[,R_index[1]]
  total_confirmed = sum(C)
  
  prob_isolated = dis$prob_isolated
  
  ## BLOCK 1: THE EPI->ECON LINK ####################
  ## response to pandemic / mandate
  
  if (integrate==1){
    econ = epi_to_econ(total_confirmed,econ,q1=data$q1,q2=data$q2)
    # labour force (lf) is the original lf minus those confirmed
    econ$lf = econ$lf - (prob_isolated * C[1] / 1e6) # in millions
  }

  ## BLOCK 2: THE ECON->EPI LINK ####################
  # first get consumption
  econ = econ$cons_link_fun(y, econ)
  
  # get relative values
  relative_consumption = 1
  relative_work = 1
  if (integrate==1){
    
    # "counter" = counterfactual; these are econ variables that came from a prior simulation with no
    # epidemic. In our simple models, they are just constant values. Here, we interpolate the values, 
    # in case we update (but it might save a lot of computation time to just use a constant)
    # time variables
    counter_time = econ$counter_time
    if(t<min(counter_time)) t = min(counter_time)
    if(t>max(counter_time)) t = max(counter_time)
    # interpolate counterfactuals
    counter_cons = pracma::interp1(x=counter_time,y=econ$counter_cons,xi = t)
    counter_worker = pracma::interp1(x=counter_time,y=econ$counter_worker,xi = t)
    # get work and consumption
    linked_vals = econ$econ_to_epi(y,econ)
    # take as proportion of counterfactual
    relative_consumption = linked_vals$cons_link/counter_cons
    relative_work = linked_vals$work_link/counter_worker #
  }
  
  ## BLOCK 3: EPI MODEL ######################
  
  ## FOI (force of infection)
  
  # we recompute the contact matrix at each time step based on the econ variables
  contact_matrix = get_scaled_contacts(data, NN0, relative_consumption = relative_consumption, 
                            relative_work = relative_work)
  I <- (1 - prob_isolated) * C + Id + Iu
  foi <- dis$beta * contact_matrix %*% (I/NN0)
  
  ## EQUATIONS
  prob_detected = dis$prob_detected
  TEtoI <- dis$TEtoI
  TItoR <- dis$TItoR
  TItoC <- dis$TItoC
  TCtoR <- dis$TCtoR
  
  new_infections = S * foi
  latent_det = E * prob_detected / TEtoI
  latent_undet = E * (1 - prob_detected) / TEtoI
  undet_recover = Iu / TItoR
  det_detected = Id / TItoC
  con_recover = C / TCtoR
  
  Sdot <- - new_infections #+ waning
  Edot <- new_infections - (latent_det + latent_undet)
  Iudot <- latent_undet - undet_recover
  Iddot <- latent_det - det_detected
  Cdot <- det_detected - con_recover
  Rdot <- undet_recover + con_recover
  
  epi_derivative_mat <- matrix(0, nrow = nrow(epi_vars_mat), ncol = ncol(epi_vars_mat))
  epi_derivative_mat[, S_index[1]] <- Sdot
  epi_derivative_mat[, E_index[1]] <- Edot
  epi_derivative_mat[, I_index[1]] <- Iudot
  epi_derivative_mat[, I_index[2]] <- Iddot
  epi_derivative_mat[, C_index[1]] <- Cdot
  epi_derivative_mat[, R_index[1]] <- Rdot
  
  epi_deriv <- as.vector(epi_derivative_mat)
  
  ## BLOCK 4: ECON MODEL ####################
  
  econ_derivs = econ$odes(t, y, econ, total_confirmed, sum(Cdot))
  
  ## return derivatives ############################
  
  # eps10 <- .Machine$double.eps * 1e12
  # epi_deriv[y < eps10] <- pmax(0, epi_deriv[y < eps10]) # exit wave was lost
  return(list(c(econ_derivs[[1]],epi_deriv)))
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
  
  # remove econ vars
  mat = if(nRemove>0) mat =  mat[ ,-c(1:nRemove)]
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


get_results_df = function(mat, epivars, econ, integrated=1){
  
  # time
  Tout <- mat[,1]
  # epi outcomes
  epi_vars = get_epi_vars(mat=mat, ldata, vars=epivars, nRemove=econ$nEconODEs+1)
  
  # individual epi states
  Cout <-  epi_vars[[epivars]]
  
  cons_scen = econ$get_cons_from_timeseries(mat,econ,Cout,ldata, integrated)
  gdp_scen = econ$get_gdp_from_timeseries(mat,econ,Cout,ldata, integrated)
  
  df = data.frame(Day = Tout,
             Consumption = cons_scen,
             GDP = gdp_scen,
             Integrated = c('Counterfactual','Integrated model')[integrated+1])
  for(i in 1:length(econ$econvarlabels))
    df[[econ$econvarlabels[i]]] = mat[,i+1]
  
  cbind(df, do.call(cbind, epi_vars))
  
}


plot_trajectories <- function(runlist,ldata){
  
  ## OUTPUT VARIABLES
  
  scen_df = get_results_df(runlist$integrated, epivars='C', econ, integrated=1)
  counter_df = get_results_df(runlist$counterfactual, epivars='C', econ, integrated=0)
  
  plotdata = rbind(scen_df, counter_df)
  
  plotout <- ggplot(reshape2::melt(plotdata,id.var=c('Day','Integrated'))) + 
    geom_line(aes(x=Day,y=value,colour=Integrated),linewidth=2) +
    facet_wrap(~variable,scales = 'free_y',nrow=1) +
    theme_bw(base_size=15) +
    theme(legend.position = 'top') + 
    labs(y='',colour='')
  
  print(plotout)
  plotout
  
}


