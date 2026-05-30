

# 117.800/5215.728*100 # march 7 -> april 5 2021


# data: list of general model parameters

data_start <- function(country_name = 'Philippines', iso3 = 'PHL') {
  
  data <- list()
  
  data$iso3 = iso3
  data$country_name = country_name
  data$nSectors <- 1
  data$tvec <- c(0, 365)
  
  data$ageindex <- list(1:3, 4:13, 14:21)
  data$adInd <- 2
  
  compindex <- list()
  
  compindex$S_index <- c(1) 
  compindex$E_index <- c(2)
  compindex$I_index <- c(3:4) # Iu, Id
  compindex$C_index <- c(5)
  compindex$R_index <- c(6)
  data$compindex <- compindex
  
  ## disease variables ############################
  
  epidemic <- list()
  epidemic$R0 <- 3
  epidemic$TEtoI <- 5
  epidemic$TItoR <- 13.5
  epidemic$TItoC <- 7.398
  epidemic$TCtoR <- epidemic$TItoR - epidemic$TItoC
  epidemic$prob_detected <- 0.604
  
  data$epidemic = epidemic
  
  return(data)
}

# curate country data
#
# data: list of country-specific and general model parameters

gather_country_data <- function(data) {
  
  ref_year = 2023
  iso3 = data$iso3
  country_name = data$country_name
  adind = data$adInd
  
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
  # dat <- get_ilostat(id = 'EAP_DWAP_SEX_AGE_RT_A', segment = 'indicator') 
  # fiveyr <- subset(dat,ref_area==iso3 & sex=='SEX_T' & time==ref_year & grepl('5YR',classif1))
  # fiveyr$pop <- c(sum(Npop[4:21]),Npop[4:13],sum(Npop[14:21]))
  
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
  
  # percentages <- subset(fiveyr,!classif1%in%c('AGE_5YRBANDS_TOTAL','AGE_5YRBANDS_YGE65'))$obs_value
  # populations <- Npop[ageindex[[adind]]]
  # lfpr <- sum(percentages*populations/100)/sum(populations)
  # data$employmentrate = lfpr
  
  # workforce participation rate
  # LFPR 15–64 (total), Philippines
  lfpr_1564 = WDI(country=data$iso3, indicator="SL.TLF.ACTI.ZS", start=2010, end=2026)
  latest = lfpr_1564[order(lfpr_1564$year, decreasing=TRUE), ][1, ]
  data$lfpr = latest$SL.TLF.ACTI.ZS/100
  
  unemployment = WDI(country=data$iso3, indicator="SL.UEM.TOTL.ZS", start=2010, end=2026)
  latest = unemployment[order(unemployment$year, decreasing=TRUE), ][1, ]
  data$employmentrate = (100 - latest$SL.UEM.TOTL.ZS)
  
  # generate basic contact components
  data <- get_basic_contacts(data)
  
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



four_to_three = function(cm,f1,f0){
  threebyfour = matrix(0, nrow = 3, ncol = 4)
  threebyfour[1,] = cm[2,]
  threebyfour[2,] = f1*cm[1,] + f0*cm[3,]
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
  
  data$full_cms = normalised_cms
  data$reduced_cms = lapply(normalised_cms, function(x)  collapse_cm(contact_matrix = x, age_groups16, age_counts16))
  
  # popsizes = sapply(age_groups16, function(x) sum(age_counts16[x]))
  
  data = decompose_contacts(data, consumption_contact=0.4, print_values=T)
  data$contacts$basic_contact_matrix <- get_scaled_contacts(data)
  
  # combine country and disease parameters
  CI <- get_candidate_infectees(data)
  data$epidemic$beta <- data$epidemic$R0 / CI
  
  data
}


##!! changing work-related contacts changes the counterfactual, because it changes the basic contact matrix
# changing consumption-related contacts should not change the counterfactual because the basic contact matrix is unchanged
decompose_contacts = function(data, consumption_contact=0.4, print_values=F){
  
  worker_index = data$adInd
  popsizes = data$Npop3
  workers_by_sector <- data$lfpr*popsizes[worker_index] # sum(allsectors)
  # put into order: workers, then children, non-workers, and retired
  NNs <- as.numeric(c(workers_by_sector,
                      popsizes[1:(worker_index-1)],
                      popsizes[worker_index]-sum(workers_by_sector),
                      popsizes[(worker_index+1):length(popsizes)]))
  
  data$NNs <- NNs
  data$nStrata <- length(data$NNs)
  
  reduced_cms = data$reduced_cms
  f1 = data$lfpr
  if(print_values) cat(paste0("Labour force participation rate: ",f1,'\n'))
  f0 = 1-f1
  
  worker_age_contacts = sum(reduced_cms$all[worker_index,])
  worker_age_self_contacts = reduced_cms$all[worker_index,worker_index]
  frac_work_contacts = sum(reduced_cms$work[worker_index,]) / worker_age_contacts
  if(print_values) cat(paste0("Fraction of workers' contacts made at work: ",frac_work_contacts,'\n'))
  
  popfracs = popsizes/sum(popsizes)
  
  # population size vector with working-age people split between two groups
  N_4 = c(f1*popsizes[worker_index], popsizes)
  N_4[1+worker_index] = f0*popsizes[worker_index]
  
  # mu=\frac{pN_1+qN_3}{N_1}
  mu = (f1*N_4[1] + f0*N_4[3])/N_4[1]
  
  # fractions of workplace contacts made with others by workers
  # need to use conmat, or some other source
  v_vec = reduced_cms$work[worker_index,] / sum(reduced_cms$work[worker_index,])
  # the fraction of adults making contact with workers as customers 
  vstar = popfracs[worker_index] * sum(v_vec[-worker_index]) / sum(popfracs[-worker_index])
  if(print_values) cat(paste0("Fraction of workplace contacts made with non-working adults: ",vstar,'\n'))
  # the fraction of workplace contacts made with other workers
  Omega = v_vec[2] - vstar
  if(print_values) cat(paste0("Fraction of workplace contacts made with other workers: ",Omega,'\n'))
  
  
  # the total number of workplace contacts between workers
  WW = Omega * frac_work_contacts * worker_age_contacts / f1
  CW3 = vstar / (Omega * (mu/2 +f0)) * WW
  CW2_plus_CW4 = WW / Omega - WW - (mu/2 - f1 + 1)*CW3 
  CW2 = CW2_plus_CW4 * v_vec[1] / (v_vec[1] + v_vec[3])
  CW4 = CW2_plus_CW4 * v_vec[3] / (v_vec[1] + v_vec[3])
  CC2 = reduced_cms$all[2,1] - f1*CW2
  CC4 = reduced_cms$all[2,3] - f1*CW4
  CC3 = reduced_cms$all[2,2] - f1*WW - f1*CW3 * (2*f0 + mu)
  
  Mww <- Mcw <- Mcom <- Mcc <- matrix(0, nrow = 4, ncol = 4)
  Mww[1,1] = WW
  
  Mcw[1,] = c(mu*CW3, CW2, f0*CW3, CW4)
  Mcw[2:4,1] = c(N_4[1] / N_4[2] * CW2, f1*CW3, N_4[1] / N_4[4] * CW4)
  
  Mcom[2:4, 2:4] = reduced_cms$all
  Mcom[1,] = c(f1*CC3, CC2, f0*CC3, CC4)
  Mcom[3,] = c(f1*CC3, CC2, f0*CC3, CC4)
  Mcom[2:4,1] = c(N_4[1] / N_4[2] * CC2, f1*CC3, N_4[1] / N_4[4] * CC4)
  Mcom[2:4,3] = c(N_4[3] / N_4[2] * CC2, f0*CC3, N_4[3] / N_4[4] * CC4)
  
  ##!! this is possibly not the right way to separate these components
  # what fraction of community contacts are associated with consumption?
  Mcc = consumption_contact * Mcom
  # Mcc[2:4, 2:4] = consumption_contact * reduced_cms$other
  # Mcc[1,] = Mcc[3,] 
  # Mcc[3,] = Mcc[3,] 
  # Mcc[,1] = Mcc[,3] * f1 
  # Mcc[,3] = Mcc[,3] * f0
  # four_to_three(Mcc,f1,f0)
  
  # subtract consumer--consumer contacts 
  Mcom = (1 - consumption_contact) * Mcom
  
  ## save
  contacts <- list()
  contacts$Mww <- Mww
  contacts$Mcw <- Mcw
  contacts$Mcc <- Mcc
  contacts$Mcom <- Mcom
  data$contacts = contacts
  # print(four_to_three(get_scaled_contacts(data),f1,f0))
  
  return(data)
}



# construct contact matrix from component matrices
#
# data: list of general model parameters
# relative_consumption: consumption relative to counterfactual / normal times
# relative_work: number of workers relative to counterfactual / normal times

get_scaled_contacts <- function(data, relative_consumption=1, relative_work=1) {
  
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
# data: list of general and country-specific model parameters
#
# CI: candidate infectees (R0 (reproduction number) / beta)

get_candidate_infectees <- function(data) {
  
  epidemic = data$epidemic
  
  S <- N <- data$NNs
  contact_matrix = data$contacts$basic_contact_matrix
  nStrata = length(N)
  
  # Rates
  prob_detected = epidemic$prob_detected
  prob_undetected = 1 - prob_detected
  sig1 <- prob_undetected / epidemic$TEtoI
  sig2 <- prob_detected / epidemic$TEtoI
  g1 <- 1 / epidemic$TItoR
  g2 <- 1 / epidemic$TItoC
  ##!! compute beta from R0 assuming no one is confirmed
  # g2 <- 1 / epidemic$TItoR
  # confirmed <- 1 / epidemic$TCtoR
  
  FOIin <- contact_matrix * t(pracma::repmat(S, nStrata, 1)) / pracma::repmat(N, nStrata, 1)
  
  # force of infection from groups E (none), Iu (undetected), Id (not detected yet), C (detected)
  Fmat <- matrix(0, 3 * nStrata, 3 * nStrata)
  Fmat[1:nStrata, (nStrata + 1):(3*nStrata)] <- cbind(FOIin, FOIin)
  
  # rate of exit from states
  ones <- matrix(1,nStrata,1)
  vvec <- c((sig1 + sig2)*ones, g1*ones, g2*ones)
  
  n <- length(vvec)
  V <- diag(vvec)
  nmat <- diag(rep(1,nStrata));
  V[(nStrata + 1):(2 * nStrata), 1:nStrata] <- -sig1*nmat
  V[(2 * nStrata + 1):(3 * nStrata), 1:nStrata] <- -sig2*nmat
  # V[(3 * nStrata + 1):(4 * nStrata), (2 * nStrata + 1):(3 * nStrata)] <- -g2*nmat
  
  NGM <- Fmat %*% Matrix::solve(V)
  ev <- eigen(NGM,only.values = T,symmetric = F) #largest in magnitude (+/-) 
  d <- ev$values
  CI <- max(Re(d))
  
  return(CI)
}


simulate_epi_econ_model <- function(data,  econ) {
  
  ## PARAMETERS
  nStrata <- length(data$NNs)
  # nSectors <- data$nSectors
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
    contacts = data$contacts,
    epidemic = data$epidemic
  )
  
  # initial conditions
  
  imported <- 200 / sum(S0) * S0
  t0 <- data$tvec[1]
  epi_init_mat <- matrix(0, nrow = nStrata, ncol = nStates)
  epi_init_mat[, compindex$S_index[1]] <- S0 - imported
  epi_init_mat[, compindex$E_index[1]] <- imported
  
  y0 <- c(econ$econ_init, as.vector(epi_init_mat))
  
  # run for counterfactual
  
  econ$integrate <- 0
  
  fun <- function(t, y, p) ODEs(data, t, y, econ)
  tmpout <- deSolve::ode(times = seq(t0,tend,1), y = y0, func = fun,
                         parms=list(data=rundata, nStrata=nStrata, econ=econ))
  # plot(rowSums(tmpout[,12:23]))
  # store counterfactual values
  # econ$counter_time <- tmpout[,1]
  # linking_values = t(apply(tmpout[,-1],1,function(y){ econ = econ$cons_link_fun(y, econ);
  #   unlist(econ$econ_to_epi(y,econ))}))
  # econ$counter_cons <- linking_values[,1]
  # econ$counter_worker <- linking_values[,2]
  econ$integrate <- 1
  
  # solve ODEs
  out <- deSolve::ode(times = seq(t0, tend, by=1), y = y0, func = fun,
                      parms=list(data=rundata, nStrata=nStrata, econ=econ))
  
  ## OUTPUTS:  
  returnobject <- list(
    data=data,
    integrated=out,
    counterfactual=tmpout
  )
  
  return(returnobject)
}


epi_to_econ = function(epi_var, lf_confirmed, econ, data){
  
  q1 = data$q1
  q2 = data$q2
  
  # make sure we start at one
  # ref_scalar = 1/(1+exp(-(ref_val)/gradient))
  # scalar = 1/(1+exp(-(ref_val-epi_var)/gradient)) / ref_scalar
  # prop_to = baseline + (1-baseline) * scalar
  
  # approximate notifications as total confirmed divided by duration confirmed
  epi_var_new = epi_var/data$epidemic$TCtoR
  
  prop_to = q1 + (1-q1) / (1+q2*epi_var_new)
  
  econ$scalar = prop_to
  
  # labour force (lf) is the original lf minus those confirmed
  econ$lf = econ$lf - (lf_confirmed / 1e6) # in millions
  
  return(econ)
  
}

epi_to_econ_deriv = function(epivar, dot_epivar, data){
  
  q1 = data$q1
  q2 = data$q2
  TCtoR = data$epidemic$TCtoR
  # approximate notifications as total confirmed divided by duration confirmed
  epivar_new = epivar/TCtoR
  dot_epivar_new = dot_epivar/TCtoR
  
  deriv = (-(1-q1) * q2 * dot_epivar_new )/((1+q2*epivar_new )^2)
  deriv
}

ODEs <- function(data, t, y, econ) {
  
  ## BLOCK 0: ACCESS VARIABLES ####################
  epidemic = data$epidemic
  
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
  
  
  prob_detected = epidemic$prob_detected
  TEtoI <- epidemic$TEtoI
  TItoR <- epidemic$TItoR
  TItoC <- epidemic$TItoC
  TCtoR <- epidemic$TCtoR
  
  latent_det = E * prob_detected / TEtoI
  latent_undet = E * (1 - prob_detected) / TEtoI
  undet_recover = Iu / TItoR
  det_detected = Id / TItoC
  con_recover = C / TCtoR
  
  total_notifications = sum(det_detected)
  
  ## BLOCK 1: THE EPI->ECON LINK ####################
  ## response to pandemic / mandate
  
  if (integrate==1){
    econ = epi_to_econ(epi_var = total_confirmed, lf_confirmed = C[1], econ, data)
  }

  ## BLOCK 2: THE ECON->EPI LINK ####################
  # first get consumption
  econ = econ$cons_link_fun(y, econ)
  
  # get relative values
  relative_consumption = 1
  relative_work = 1
  if (integrate==1){
    
    # get work and consumption
    linked_vals = econ$econ_to_epi(y,econ)
    # take as proportion of counterfactual
    relative_consumption = linked_vals$cons_link/econ$cons0
    relative_work = linked_vals$work_link/econ$emp0 #
  }
  
  ## BLOCK 3: EPI MODEL ######################
  
  ## FOI (force of infection)
  
  # we recompute the contact matrix at each time step based on the econ variables
  contact_matrix = get_scaled_contacts(data, relative_consumption = relative_consumption, 
                            relative_work = relative_work)
  I <- Id + Iu
  foi <- epidemic$beta * contact_matrix %*% (I/NN0)
  
  ## EQUATIONS
  
  new_infections = S * foi
  
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
  
  econ_derivs = econ$odes(t, y, econ, C, Cdot)
  
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


get_results_df = function(runlist, epivars, econ, ldata){
  popsize = sum(ldata$Npop3)
  plotnames = c('Counterfactual','Integrated model')
  listnames = c('counterfactual','integrated')
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
    
    consandgdp = econ$get_cons_and_gdp_from_timeseries(y=mat,econ,epivar=Cout,data=ldata, integrate)
    cons_scen = consandgdp[[1]][[1]]
    gdp_scen = consandgdp[[2]]
    
    cons_list <- c()
    for(i in 2:3) cons_list[[i-1]] = data.frame(Day = Tout,
                         Integrated = c('Supply','Demand')[i-1],
                         variable = 'Consumption',
                         value = consandgdp[[1]][[i]]/consandgdp[[1]][[1]][1]*100)
    cons_df = do.call(rbind,cons_list)
                         
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
    
    all_infected = Reduce('+',get_epi_vars(mat=mat, ldata, vars=c('C','Id','Iu'), nRemove=econ$nEconODEs+1))
    max_infected = max(all_infected)/1e6 # in millions
    recovered = get_epi_vars(mat=mat, ldata, vars='R', nRemove=econ$nEconODEs+1)
    cumulative_incidence_pc = max(recovered$R)/popsize*100
    # if(integrate==1)
      cat(paste0(plotnames[j],
        # ldata$q2,' & ',
        # econ$lambda_p1,' & ',
        # ldata$cc,' & ',
        ' & ',
        round(max_infected,1),' & ',
        scen_df[[j]]$Day[which.max(all_infected)],' & ',
        round(cumulative_incidence_pc),' & ',
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


plot_trajectories <- function(runlist, econ, ldata){
  
  ## OUTPUT VARIABLES
  
  plotdata = get_results_df(runlist, epivars='C', econ, ldata)
  plotdata$Integrated = factor(plotdata$Integrated, levels=c('Counterfactual','Integrated model','Supply','Demand'))
  
  plotdata
  
}


