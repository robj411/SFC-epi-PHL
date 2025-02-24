
# load in basic data and make general definitions
#
# data: struct of general model parameters

data_start <- function() {
  
  data <- list()
  # closurefile <- '../Daedalus-P2-Dashboard/data/closures.xlsx'
  # sheets <- excel_sheets(closurefile)
  # for (i in 1:length(sheets)) {
  #   thissheet <- sheets[i]
  #   closurei <- as.matrix(read_xlsx(closurefile, sheet = thissheet))
  #   data[[thissheet]] <- closurei
  # }
  
  data$adInd <- 3
  data$nSectors <- 1#nrow(data$x_elim)
  data$tvec <- c(0, 365*2)
  # data$EdInd <- 41 #education sector index
  # data$HospInd <- c(32,43,44) #hospitality sector indices
  
  contacts <- list()
  contacts$sectorcontacts <- as.numeric(read.csv('../Daedalus-P2-Dashboard/data/sectorcontacts.csv')$n_cnt)
  contacts$sectorcontactfracs <- read.csv('../Daedalus-P2-Dashboard/data/sec_contact_dist_UK.csv')
  data$contacts <- contacts
  
  data$ageindex <- list(1, 2:3, 4:14, 15:21)
  
  compindex <- list()
  
  compindex$S_index <- c(1) # S, Sn, S01, Sv, Sb, S02, S12
  compindex$E_index <- c(2)
  compindex$I_index <- c(3:4)
  compindex$H_index <- c(5)
  compindex$R_index <- c(6)
  compindex$D_index <- 7
  data$compindex <- compindex
  
  return(data)
}



load_country_data <- function() {
  CD <- read.csv('../Daedalus-P2-Dashboard/data/country_data.csv')
  CD$popsum <- rowSums(CD[, 4:24])
  
  CMcols <- grep('CM', names(CD))
  dim <- sqrt(length(CMcols))
  Npopcols <- grep('Npop', names(CD))
  average_contacts <- rep(NA, nrow(CD))
  for (i in 1:nrow(CD)) {
    cm <- as.matrix(CD[i, CMcols])
    cmm <- matrix(cm, nrow = dim, ncol = dim, byrow = TRUE)
    Npop <- as.numeric(CD[i, Npopcols])
    Npop[dim] <- sum(Npop[dim:length(Npop)])
    Npop <- Npop[1:dim]
    if (!is.na(cm[1])) {
      average_contacts[i] <- sum(rowSums(cmm) * Npop) / sum(Npop)
    }
  }
  CD$average_contacts <- average_contacts
  
  country_parameter_distributions <- read.csv('../Daedalus-P2-Dashboard/data/parameter_distributions.csv')
  country_parameter_distributions$Parameter.2[grep('gam',country_parameter_distributions$distribution)] <- 
    1/country_parameter_distributions$Parameter.2[grep('gam',country_parameter_distributions$distribution)]
  
  return(list(CD, country_parameter_distributions))
}



# simulate a random country by drawing from distributions and data
#
# data: struct of general model parameters
# CD: table of country data values
# income_level: string indicating income level (e.g. HIC)
# country_parameter_distributions: pre-specified, named distributions and
# parameters
# social_dist_coefs: table of parameters for social distancing function
#
# data: struct of general model parameters

p2RandCountry <- function(data, CD, income_level, country_parameter_distributions,country_name) {
  
  # start
  contacts <- data$contacts
  
  cdnames <- colnames(CD)
  
  # values from distributions
  pindices <- which(country_parameter_distributions$igroup == income_level | 
                      (country_parameter_distributions$igroup == 'all' & 
                         country_parameter_distributions$distribution != 'NA'))
  country_ind = which(CD$country==country_name)
  cpd <- country_parameter_distributions[pindices, ]
  cpd$distribution <- paste0('q',gsub('inv','',cpd$distribution))
  cpd$distribution <- gsub('gam','gamma',cpd$distribution)
  cpd$distribution <- gsub('logn','lnorm',cpd$distribution)
  for (i in 1:nrow(cpd)) {
    qname <- paste0(cpd$parameter_name[i], '_quantile')
    qexp <- paste0(' <- runif(1, 0, 1);')
    eval(parse(text = paste0(qname, qexp)))
    varname <- cpd$parameter_name[i]
    expression <- paste0(' <- ', cpd$distribution[i], '(', 
                         cpd$parameter_name[i], '_quantile, ', 
                         cpd$Parameter.1[i], ', ', 
                         cpd$Parameter.2[i], ');')
    eval(parse(text = paste0(varname, expression)))
  }
  
  # store values
  data$remote_quantile <- internet_coverage_quantile
  data$response_time_quantile <- runif(1, 0, 1)
  data$remote_teaching_effectiveness <- runif(1, 0, 1)
  
  # contacts.pt = pt
  contacts$school1_frac <- school1_frac
  contacts$school2_frac <- school2_frac
  contacts$hospitality_frac <- c(hospitality1_frac, hospitality2_frac, hospitality3_frac, hospitality4_frac)
  
  dindices <- which(grepl( 'hospitality_age',country_parameter_distributions$parameter_name))
  cpd <- country_parameter_distributions[dindices, ]
  contacts$hospitality_age <- matrix(0, nrow = 3, ncol = nrow(cpd))
  for (i in 1:nrow(cpd)) {
    p3 <- cpd$Parameter.1[i]
    p4 <- cpd$Parameter.2[i]
    contacts$hospitality_age[, i] <- c(1 - p3 - p4, p3, p4) # rdirichlet(1, c(1 - p3 - p4, p3, p4) * 10)
  }
  
  data$Hmax <- CD$Hmax[country_ind]
  labs <- read_dta(file.path('../Daedalus-P2-Dashboard/data/','lab_share_data.dta'))
  sublab <- subset(labs,countrycode=='MEX')
  labsh <- subset(sublab,year==max(year))$labsh
  data$labsh <- labsh
  
  # values by sampling
  
  if (income_level == 'LLMIC') {
    country_indices <- CD$igroup == 'LIC' | CD$igroup == 'LMIC'
  } else if (income_level == 'UMIC') {
    country_indices <- CD$igroup == 'UMIC'
  } else if (income_level == 'HIC') {
    country_indices <- CD$igroup == 'HIC'
  }
  
  # population
  # population by age
  cols <- grep('Npop', cdnames)
  Npop <- as.numeric(CD[country_ind, cols])
  data$Npop <- Npop
  Npop4 <- sapply(data$ageindex, function(x) sum(Npop[x]))
  data$Npop4 <- Npop4
  
  # population by stratum
  # sample workforce
  nonempind <- which(!is.na(CD$NNs1) & country_indices)
  colNNs <- grep('NNs', cdnames)
  allsectors <- as.matrix(CD[country_ind, colNNs])
  contacts$sectorcontacts = sum(allsectors*contacts$sectorcontacts)/sum(allsectors)
  contacts$sectorcontactfracs = c(colSums(allsectors*contacts$sectorcontactfracs)/sum(allsectors))[2:4]
  
  # correlation between workforce in place and contacts from work
  Z <- MASS::mvrnorm(1, c(0, 0), matrix(c(1, 0.8, 0.8, 1), nrow = 2))
  U <- pnorm(Z, 0, 1)
  
  wipindex <- which(country_parameter_distributions$parameter_name == 'workforce_in_place' & 
                      country_parameter_distributions$igroup == income_level)
  workers_by_sector <- sum(allsectors)
  # put into daedalus order: workers by sector, then infants, adolescents,
  # non-workers, and retired
  NNs <- as.numeric(c(workers_by_sector,Npop4[1],Npop4[2],Npop4[3]-sum(workers_by_sector),Npop4[4]))
  
  wfindex <- which(country_parameter_distributions$parameter_name == 'work_frac')
  cpd <- country_parameter_distributions[wfindex,]
  
  work_frac <- qbeta(U[1],cpd$Parameter.1,cpd$Parameter.2);
  # work contact fraction should not exceed worker fraction
  contacts$work_frac <- min(work_frac, sum(workers_by_sector) / Npop4[3])
  
  # contacts
  # workplace
  sectorcontacts <- contacts$sectorcontacts
  sB <- length(sectorcontacts)
  # s1 <- sB[1]
  # s2 <- sB[2]
  # contacts$sectorcontacts <- pmax((sectorcontacts + 1) * 2 ^ runif(sB, -1, 1) - 1, 0)
  uk_ptr <- 15.87574
  # contacts$sectorcontacts[data$EdInd] <- pupil_teacher_ratio / uk_ptr * contacts$sectorcontacts[data$EdInd]
  data$pupil_teacher_ratio <- pupil_teacher_ratio
  
  # matrix
  randvalue <- as.numeric(CD[country_ind, grep('CM', cdnames)])
  defivalue <- matrix(randvalue, nrow = 16, ncol = 16)
  contacts$CM <- defivalue
  
  # wfh = work from home
  mins <- apply(CD[country_ind, grep('wfhl', cdnames)], 2, min)
  maxs <- apply(CD[country_ind, grep('wfhu', cdnames)], 2, max)
  newprop <- sum(qunif(internet_coverage_quantile, mins, maxs)*allsectors)/sum(allsectors)
  data$wfh <- matrix(newprop, nrow = 2, ncol = length(mins))
  
  # date of importation
  data$t_import <- runif(1, 0, 20)
  
  # Hres = hospital occupancy at response time in origin country
  data$Hres <- qunif(data$response_time_quantile, 1, 20)
  
  # la = life expectancy
  cols <- grep('^la', cdnames)
  life_expectancy <- as.numeric(CD[country_ind, cols])
  data$life_exp <- life_expectancy
  
  data$NNs <- NNs
  
  # finish workers by sector
  nSectors = 1
  data$NNs[data$NNs == 0] <- 1
  data$nStrata <- length(data$NNs)
  data$employmentrate <- sum(data$NNs[1:nSectors]) / Npop4[3]
  
  # generate basic contact components
  data$contacts <- get_basic_contacts(data, contacts)
  basic_contact_matrix <- p2MakeDs(data, data$NNs)
  data$contacts$basic_contact_matrix <- basic_contact_matrix
  
  return(data)
}


# take population contact matrix and decompose into items that will be
# impacted differently by configurations
#
# data: struct of general model parameters
#
# contacts: struct of contact parameters

get_basic_contacts <- function(data, contacts) {
  
  NN <- data$NNs
  CM_16 <- contacts$CM
  s16 <- nrow(CM_16)
  
  Npop <- data$Npop
  Npop[s16] <- sum(Npop[s16:length(Npop)])
  Npop <- Npop[1:s16]
  ageindex <- data$ageindex
  ageindex[[4]] <- min(ageindex[[4]]):s16
  Npop4 <- data$Npop4
  pop_props <- Npop4 / sum(Npop)
  
  ## community contacts ##################################
  
  CM_164 <- sapply(ageindex, function(x) rowSums(CM_16[,x,drop=F]))
  CM_4 <- t(sapply(ageindex, function(x) (Npop[x] %*% CM_164[x,]) / sum(Npop[x]))  )
  CM_4_orig = CM_4
  CMav <- pop_props %*% colSums(CM_4_orig)
  contact_props <- CM_4_orig[3,] / sum(CM_4_orig[3,])
  workage_total <- sum(CM_4_orig[3,])
  
  ## indices
  
  adInd <- data$adInd #Adult index
  nSectors <- 1 #data$nSectors #Number of sectors
  nStrata <- length(NN)
  workage_indices <- c(1:nSectors,nSectors+adInd)
  
  NNrel <- NN[workage_indices] / sum(NN[workage_indices]) #adult population proportion vector
  
  ## work-related contacts ##########################################
  
  sectoragedist <- matrix(0, nSectors, nStrata)
  sectorcontactfracs <- contacts$sectorcontactfracs;
  # correct for age dist
  over65frac <- Npop4[4]/sum(Npop4);
  ukover65frac <- .25;
  sectorcontactfracs[['X65plus']] <- sectorcontactfracs[['X65plus']] * over65frac / ukover65frac;
  newtotal <- sectorcontactfracs[['workingage']]+sectorcontactfracs[['X65plus']]+sectorcontactfracs[['under18']];
  sectorcontactfracs[['under18']] <- sectorcontactfracs[['under18']] / newtotal;
  sectorcontactfracs[['workingage']] <- sectorcontactfracs[['workingage']] / newtotal;
  
  prop_working = NN[1]/(NN[1]+NN[4])
  target_work_contacts <- contacts$work_frac/2 * workage_total / prop_working
  worker_contacts_adults = sectorcontactfracs[['workingage']] * target_work_contacts
  total_nonworker_contacts = workage_total - prop_working * target_work_contacts
  total_worker_contacts = total_nonworker_contacts + target_work_contacts
  rel_worker = prop_working*total_worker_contacts / ((1-prop_working)*total_nonworker_contacts)
  w1 = worker_contacts_adults * rel_worker / (1 + rel_worker)
  w2 = worker_contacts_adults - w1
  s2 = prop_working*w2/(1-prop_working)
  
  sectoragedist[,nSectors+c(1,2)] <- t(repmat(sectorcontactfracs[['under18']],2,1)) * repmat(pop_props[1:2] / sum(pop_props[1:2]),nSectors,1) 
  sectoragedist[,nSectors+4] <- sectorcontactfracs[['X65plus']]
  sectoragedist[,workage_indices] <- t(repmat(sectorcontactfracs[['workingage']],length(workage_indices),1)) * c(w1,w2)/(w1+w2)
  
  
  community_to_worker_mat <- target_work_contacts * sectoragedist
  community_to_worker_mat <- rbind(community_to_worker_mat,matrix(0,nrow=4,ncol=nStrata))
  
  worker_worker_mat <- matrix(0,nrow=nStrata,ncol=nStrata)
  worker_worker_mat[1:nSectors,1:nSectors] <- community_to_worker_mat[1:nSectors,1:nSectors]
  
  contacts$worker_worker_mat <- worker_worker_mat
  print(worker_worker_mat)
  
  community_to_worker_mat[1:nSectors,1:nSectors] <- 0
  
  total_cn_workerage = NN[workage_indices] %*% community_to_worker_mat[workage_indices,]
  consumer_contacts <- total_cn_workerage / NN
  consumer_contacts[workage_indices] <- consumer_contacts[workage_indices[length(workage_indices)]]
  worker_to_community_mat <- cbind(t(consumer_contacts),matrix(0,ncol=4,nrow=nStrata))
  
  contacts$community_worker_mat <- worker_to_community_mat + community_to_worker_mat
  print(contacts$community_worker_mat)
  
  # worker_to_community_mat + community_to_worker_mat
  
  # get marginal contacts by age for workers
  av_workerage_contacts_full <- NNrel %*% community_to_worker_mat[workage_indices,]
  av_workerage_contacts_collapsed <- c(av_workerage_contacts_full[,nSectors+1], av_workerage_contacts_full[,nSectors+2], sum(av_workerage_contacts_full[,workage_indices]), av_workerage_contacts_full[,nSectors+4])
  
  c_to_w_back <- av_workerage_contacts_collapsed * Npop4[3] / Npop4
  
  ## get new contact rates
  contacts$school1 <- CM_4_orig[1,1] * contacts$school1_frac
  contacts$school2 <- CM_4_orig[2,2] * contacts$school2_frac
  
  ## school
  school_mat = matrix(0,nrow=4,ncol=4)
  school_mat[1,1] = contacts$school1
  school_mat[2,2] = contacts$school2
  contacts$school_mat = school_mat
  
  ## subtract contacts from C4
  # school
  CM_4 <- CM_4 - school_mat
  # customer to worker
  CM_4[3,] <- pmax(CM_4_orig[3,] - av_workerage_contacts_collapsed, 0)
  # worker to customer
  # CM_4[,3] <- pmax(CM_4[,3] - c_to_w_back, 0)
  CM_4[c(1,2,4),3] = CM_4[c(1,2,4),3] - c_to_w_back[c(1,2,4)]
  
  # hospitality
  
  ##!! too many work contacts to infants
  hospitality_age <- contacts$hospitality_age
  hospitality_age <- rbind(hospitality_age[1,], hospitality_age)
  total_contacts <- rowSums(CM_4)
  contacts$hospitality_contacts <- t(repmat(total_contacts * contacts$hospitality_frac,4,1)) %*% hospitality_age
  CM_4 <- pmax(CM_4 - contacts$hospitality_contacts, 0)
  
  ## save
  
  # contacts$sectorcontactfracs <- NULL
  contacts$CM_4 <- CM_4
  contacts$contact_props <- contact_props
  
  return(contacts)
}



# construct contact matrices from components for configurations
#
# data: struct of general model parameters
# NN: population by stratum
# x: economic configuration
# hw: proportion working from home by stratum
#
# D: contact matrix

p2MakeDs <- function(data, NN, relative_consumption=1, relative_work=1, home_working=0) {
  
  ## variables to use
  contacts <- data$contacts
  CM_4 <- contacts$CM_4
  contact_props <- contacts$contact_props
  
  workers_present <- relative_work
  # edInd <- data$EdInd
  adInd    <- data$adInd #Adult index
  # HospInd <- data$HospInd
  
  nSectors <- length(relative_work) #Number of sectors
  nStrata <- length(NN)
  workage_indices <- c(1:nSectors, nSectors + adInd)
  
  NNrel <- data$NNs[workage_indices] / sum(data$NNs[workage_indices]) #adult population proportion vector
  NNrepvecweighted <- rep(0, nStrata)
  NNrepvecweighted[workage_indices] <- NNrel * contact_props[3]
  NNrepvecweighted[nSectors + c(1, 2, 4)] <- contact_props[c(1, 2, 4)]
  NNrep <- matrix(rep(NNrepvecweighted, nStrata), nrow = nStrata, byrow = TRUE) #total population proportion matrix
  NNrea <- matrix(rep(data$NNs[1:nSectors] / sum(data$NNs[1:nSectors]), nSectors), nrow = nSectors, byrow = TRUE) #workforce population proportion matrix
  
  ## add school and hospitality contacts to CM_4
  
  # start with hospitality
  # hospitality_sectors <- NN[HospInd]
  # get weighted average
  # hospitality_sectors <- sum(hospitality_sectors * openness[HospInd]) / sum(hospitality_sectors) # constant from 0 to 1, weighted measure of how much sectors are open
  CM_4 <- CM_4 + relative_consumption^2 * contacts$hospitality_contacts
  workRow <- CM_4[adInd, ]
  
  # school
  CM_4 <- CM_4 + contacts$school_mat # workers_present[edInd]^2 * 
  
  ## Make community matrix
  community_mat <- matrix(0, nrow = nStrata, ncol = nStrata)
  community_mat[(nSectors + 1):nStrata, (nSectors + 1):nStrata] <- CM_4
  community_mat[1:nSectors, (nSectors + 1):nStrata] <- matrix(rep(workRow, nSectors), nrow = nSectors, byrow = TRUE)
  ad_row = community_mat[, nSectors + adInd]
  for(i in 1:length(workage_indices))
    community_mat[,workage_indices[i]] = ad_row * NNrel[i]
  # community_mat[, workage_indices] <- t(repmat(community_mat[, nSectors + adInd],nSectors+1,1)) * repmat(NNrel,nStrata,1)
  # browser()
  ## WORKER-WORKER AND COMMUNITY-WORKER MATRICES
  
  effective_openness <- pmax(0, relative_work - home_working)
  effective_openness <- c(effective_openness,rep(0,4))
  # effective_openness_mat <- t(repmat(effective_openness^2, nStrata, 1))
  
  # customer--worker
  communityworker_mat <- contacts$community_worker_mat * effective_openness * relative_consumption
  # browser()
  
  # worker--worker
  workerworker_mat <- contacts$worker_worker_mat * effective_openness^2
  
  ## add all together
  # contact_matrix <- community_mat + communitytoworker_mat + worker_back
  
  # return(contact_matrix)
  
  # worker--worker contacts
  # wwcontacts <- matrix(0, nrow = nStrata, ncol = nStrata)
  # wwcontacts[1:nSectors, 1:nSectors] <- (communitytoworker_mat[1:nSectors, 1:nSectors] + worker_back[1:nSectors, 1:nSectors]) 
  
  # community--worker contacts
  # cwcontacts <- matrix(0, nrow = nStrata, ncol = nStrata)
  # cwcontacts[1:nSectors, (nSectors + 1):nStrata] <- communitytoworker_mat[1:nSectors, (nSectors + 1):nStrata]
  # cwcontacts[(nSectors + 1):nStrata, 1:nSectors] <- worker_back[(nSectors + 1):nStrata, 1:nSectors]
  
  # add all together
  contact_matrix <- community_mat + workerworker_mat + communityworker_mat # wwcontacts + cwcontacts
  
  return(contact_matrix)
  
}

collapse_cm <- function(contact_matrix, NNs) {
  
  Ci <- cbind(contact_matrix[, 2], contact_matrix[, 3], rowSums(contact_matrix[, c(1, 4)]), contact_matrix[, 5])  # sum of the columns
  cm <- rbind(Ci[2, ],
              Ci[3, ],
              rowSums(NNs[c(1, 4)] * t(Ci[c(1, 4), ])) / sum(NNs[c(1, 4)]),
              Ci[5, ])
  
  return(cm)
}



# get parameters that characterise p2 (pandemic preparedness)
#
# data: struct of general model parameters
# dis: struct of pathogen parameters
# vaccine_day: 100 or 365, depending on SARS-X vaccine scenario
# bpsv: 0 or 1, depending on BPSV scenario (present or absent)
#
# data: struct of general model parameters
# dis: struct of pathogen parameters
# p2: struct of p2 intervention parameters

p2Params <- function(data, dis) {
  
  ## PREPAREDNESS PARAMETERS
  
  p2 <- list()
  
  # response time and testing start time set by global alert
  p2$Tres <- 30 # data$response_time
  # Hospital Capacity
  Hmax  <- data$Hmax * sum(data$Npop) / 10^5   #Hospital Capacity
  data$Hmax <- NULL
  p2$hosp_release_trigger   <- max(1, 0.25 * Hmax) #lower threshold can't be less than 1 occupant
  p2$Hmax  <- max(4 * p2$hosp_release_trigger, Hmax)
  # p2$SHmax <- 2 * p2$Hmax
  
  return(list(data = data, dis = dis, p2 = p2))
  
}



# function to combine population and disease parameters to get
# within-country pathogen parameters
#
# data: struct of general model parameters
# dis: struct of pathogen parameters
# R0betafun: the function that computes beta from R0
#
# dis: struct of pathogen parameters

population_disease_parameters <- function(data, dis, R0betafun) {
  
  ## COUNTRY PARAMETERS
  ## INITIAL DISEASE PARAMETERS
  # Population by Age
  nSectors <- data$nSectors
  adInd <- data$adInd
  subs0 <- 1:4
  lihr <- length(dis$ihr)
  Npop <- data$Npop
  Npop <- c(Npop[1:(lihr-1)], sum(Npop[lihr:length(Npop)])) # last age range for disease is 80+
  ageindex <- data$ageindex
  ageindex[[4]] <- min(ageindex[[4]]):lihr
  
  ranges <- sapply(ageindex, length)
  Npop4rep <- rep(data$Npop4, ranges)
  nnprop <- Npop / Npop4rep
  subs <- rep(subs0, ranges)
  
  # Probabilities
  probHgivenSym <- dis$ihr / dis$prob_symp
  probDgivenH <- dis$ifr / dis$ihr
  dis$hfr <- probDgivenH
  probHgivenSym <- as.numeric(tapply(probHgivenSym * nnprop, subs, sum))
  dis$prob_H <- c(rep(probHgivenSym[adInd], nSectors), probHgivenSym)
  nnh <- Npop * dis$ihr
  nnhtot <- sapply(ageindex, function(x) sum(nnh[x]))
  nnhtot <- rep(nnhtot, ranges)
  nnhprop <- nnh / nnhtot
  probDgivenH <- tapply(probDgivenH * nnhprop, subs, sum)
  dis$prob_D <- c(rep(probDgivenH[adInd], nSectors), probDgivenH)
  
  dis$rr_infection <- 1 # c(rep(data$bmi_rr[1, 1], nSectors), 1, 1, data$bmi_rr[1, 1], data$bmi_rr[2, 1])
  
  # Durations
  dis$TIs <- ((1 - dis$prob_H) * dis$TIstoR) + (dis$prob_H * dis$TIstoH)
  dis$TH <- ((1 - dis$prob_D) * dis$THtoR) + (dis$prob_D * dis$THtoD)
  
  NNs <- data$NNs
  zs <- rep(0, length(NNs))
  dis$CI <- get_candidate_infectees(nStrata=length(NNs), dis, S=NNs, N=NNs, contact_matrix=data$contacts$basic_contact_matrix)
  
  R0beta <- R0betafun(dis)
  dis$R0 <- R0beta[1]
  dis$beta <- R0beta[2]
  
  return(dis)
}



## get effective reproduction number
#
# nStrata: number of strata
# dis: struct of pathogen parameters
# S: susceptible unvaccinated
# Sv1: susceptible BPSV-vaccinated
# Sv2: susceptible SARS-X--vaccinated
# p3: fraction of asymptomatic infectious people's infectiousness averted
# p4: fraction of symptomatic infectious people's infectiousness averted
# N: population by stratum
# contact_matrix: contact matrix
#
# R: effective reproduction number

get_candidate_infectees <- function(nStrata, dis, S, N, contact_matrix) {
  
  
  # Rates
  sig1 <- (1 - dis$prob_symp) / dis$TEtoI
  sig2 <- dis$prob_symp / dis$TEtoI
  g1 <- 1 / dis$TIatoR
  g2 <- (1 - dis$prob_H) / dis$TIs
  g3 <- (1 - dis$pd) / dis$Th
  h <- dis$prob_H / dis$TIs
  
  red <- dis$asym_rel_transmission
  
  S_sum <- S 
  
  FOIin <- contact_matrix * t(repmat(S_sum * dis$rr_infection,nStrata,1)) / repmat(N, nStrata, 1)
  
  Fmat <- matrix(0, 3 * nStrata, 3 * nStrata)
  Fmat[1:nStrata, (nStrata + 1):(3*nStrata)] <- cbind(red * FOIin, FOIin)
  
  g2hweighted <- (g2 + h) 
  
  ones <- matrix(1,nStrata,1)
  vvec <- c((sig1 + sig2)*ones, (g1)*ones, g2hweighted*ones)
  
  n <- length(vvec)
  V <- diag(vvec)
  nmat <- diag(rep(1,nStrata));
  V[(nStrata + 1):(2 * nStrata), 1:nStrata] <- -sig1*nmat
  V[(2 * nStrata + 1):(3 * nStrata), 1:nStrata] <- -sig2*nmat
  
  NGM <- Fmat %*% Matrix::solve(V)
  ev <- eigen(NGM,only.values = T,symmetric = F) #largest in magnitude (+/-) 
  d <- ev$values
  CI <- max(Re(d))
  
  return(CI)
}


p2SimVax <- function(data, dis, p2, econ) {
  
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
    ref_val = data$ref_val,
    baseline = data$baseline,
    compindex = data$compindex,
    contacts = data$contacts
  )
  
  # initial conditions
  
  imported <- 5000 / sum(S0) * S0
  t0 <- data$tvec[1]
  epi_init_mat <- matrix(0, nrow = nStrata, ncol = nStates)
  epi_init_mat[, compindex$S_index[1]] <- S0 - imported
  epi_init_mat[, compindex$E_index[1]] <- imported
  
  y0 <- c(econ$econ_init, as.vector(epi_init_mat))
  
  # run for counterfactual
  
  econ$integrate <- 0
  
  fun <- function(t, y, p) ODEs(data, i, t, dis, y, p2, econ)
  tmpout <- deSolve::ode(times = seq(t0,tend,1), y = y0, func = fun,
                         parms=list(data=rundata, nStrata=nStrata, dis=dis, i=1, p2=p2, econ=econ),
                         method='impAdams_d')
  econ$counter_time <- tmpout[,1]
  linking_values = t(apply(tmpout[,-1],1,function(y)unlist(econ$epi_econ_link(y,econ))))
  econ$counter_cons <- linking_values[,1]
  econ$counter_worker <- linking_values[,2]
  econ$integrate <- 1
  
  
  # store outputs
  ode_out_list <- list()
  
  
  ## LOOP
  
  i <- 1
  isequence <- NULL
  
  while (t0 < tend) {
    
    isequence <- rbind(isequence, c(t0, i))
    
    # solve ODEs
    obligatory_times <- c(t0,seq(ceiling(t0),floor(tend),by=1),tend)
    if(p2$Tres > t0) obligatory_times <- sort(c(obligatory_times,p2$Tres))
    
    out <- deSolve::ode(times = obligatory_times, y = y0, func = fun,
                        parms=list(data=rundata, nStrata=nStrata, dis=dis, i=i, p2=p2, econ=econ),
                        rootfun = mitigate, method='impAdams_d')
    
    ode_out_list[[length(ode_out_list)+1]] <- out
    
    nTime = dim(out)[1]
    y0 = out[nTime, -1]
    t0 <- out[nTime, 1]
    i_event <- rev(which(attributes(out)$iroot==1))[1]
    if (!is.na(i_event)) {
      i <- data$inext[i_event]
    }
    
  }
  
  collapse_ode_out <- do.call(rbind, ode_out_list)
  max_times <- sapply(ode_out_list,function(x)max(x[,1]))
  data$tvec <- unique(c(data$tvec, max_times, tend))
  data$isequence = isequence
  # data$econ = econ
  
  ## OUTPUTS:  
  returnobject <- list(
    data=data,
    integrated=collapse_ode_out,
    counterfactual=tmpout
  )
  
  return(returnobject)
}


mitigate <- function(t, y, parms) {
  return(1)
  
}


fear_of_infection = function(epi_var,econ,
                           gradient = 10000 # small value => sharp corners
                           , ref_val = 200000 # value whereabouts change in behaviour occurs
                           , baseline = .5 # minimum value
){
  # baseline=.5
  scalar = 1/(1+exp(-(ref_val-epi_var)/gradient))
  prop_to = baseline + (1-baseline) * scalar
  
  for(val in econ$p_to_scale)
    econ[[val]] = prop_to * econ[[val]];
  
  return(econ)
  
}

ODEs <- function(data, i, t, dis, y, p2, econ) {
  
  ## BLOCK 0: ACCESS VARIABLES ####################
  
  ## variables
  
  nEconODEs = econ$nEconODEs
  
  compindex <- data$compindex
  S_index <- compindex$S_index
  E_index <- compindex$E_index
  I_index <- compindex$I_index
  H_index <- compindex$H_index
  R_index <- compindex$R_index
  D_index <- compindex$D_index
  
  NN0 <- data$NNs
  nStrata <- length(NN0)
  nStates = max(unlist(compindex))
  
  epi_vars_mat <- matrix(y[-c(1:nEconODEs)], nrow = nStrata, ncol = nStates)
  
  # indices <- 1:nStrata
  S <- epi_vars_mat[,S_index[1]]
  E <- epi_vars_mat[,E_index[1]]
  Ia <- epi_vars_mat[,I_index[1]]
  Is <- epi_vars_mat[,I_index[2]]
  H <- epi_vars_mat[,H_index[1]]
  R <- epi_vars_mat[,R_index[1]]
  D <- epi_vars_mat[,D_index[1]]
  
  
  
  ## BLOCK 1: THE LINK ####################
  ## response to pandemic / mandate
  
  integrate = econ$integrate
  if (integrate==1){
    econ = fear_of_infection(sum(H),econ,ref_val=data$ref_val,baseline=data$baseline)
    # lf is used in model 3 but not 2 or PC
    # lf is the original lf minus those dead and in hospital
    econ$lf = econ$lf - D[1] - H[1]
  }
  
  
  
  ## BLOCK 2: ECON MODEL ####################
  
  econ_vars = y[1:econ$nEconODEs]
  econ_derivs = econ$odes(t,y,econ)

  
    
  ## BLOCK 3: EPI MODEL ####################
  
  # get relative values
  relative_consumption = 1
  relative_work = 1
  if (integrate==1){
    
    # cons = y[which(econ$econvarnames=='cons')]
    
    counter_time = econ$counter_time
    if(t<min(counter_time)) t = min(counter_time)
    if(t>max(counter_time)) t = max(counter_time)
    counter_cons = interp1(x=counter_time,y=econ$counter_cons,xi = t)
    counter_worker = interp1(x=counter_time,y=econ$counter_worker,xi = t)
    
    linked_vals = econ$epi_econ_link(y,econ)
    relative_consumption = linked_vals$cons_link/counter_cons
    relative_work = linked_vals$work_link/counter_worker #
    # print(c(t,linked_vals$cons_link,counter_cons,relative_consumption,relative_work))
  }
  
  
  ## FOI
  
  contact_matrix = p2MakeDs(data, NN0, relative_consumption = relative_consumption, 
                            relative_work = relative_work, home_working = 0)
  I <- dis$asym_rel_transmission * Ia + Is
  foi <- dis$beta * contact_matrix %*% (I/NN0)
  
  ## EQUATIONS
  
  TEtoI <- dis$TEtoI
  TIs <- dis$TIs
  prob_H = dis$prob_H
  prob_D = dis$prob_D
  TH = dis$TH
  prob_symp = dis$prob_symp
  
  new_infections = S * foi
  waning = R / dis$TRtoS
  latent_sym = E * prob_symp / TEtoI
  latent_asym = E * (1 - prob_symp) / TEtoI
  asym_recover = Ia / dis$TIatoR
  sym_recover = Is * (1-prob_H) / TIs
  sym_hosp = Is * prob_H / TIs
  hosp_recover = H * (1 - prob_D) / TH
  hosp_death = H * prob_D / TH
  
  Sdot <- - new_infections + waning
  Edot <- new_infections - (latent_sym + latent_asym)
  Iadot <- latent_asym - asym_recover
  Isdot <- latent_sym - (sym_recover + sym_hosp)
  Hdot <- sym_hosp - (hosp_recover + hosp_death)
  Rdot <- asym_recover + sym_recover + hosp_recover - waning
  DEdot <- hosp_death
  
  ## return derivatives
  
  epi_derivative_mat <- matrix(0, nrow = nrow(epi_vars_mat), ncol = ncol(epi_vars_mat))
  epi_derivative_mat[, S_index[1]] <- Sdot
  epi_derivative_mat[, E_index[1]] <- Edot
  epi_derivative_mat[, I_index[1]] <- Iadot
  epi_derivative_mat[, I_index[2]] <- Isdot
  epi_derivative_mat[, H_index[1]] <- Hdot
  epi_derivative_mat[, R_index[1]] <- Rdot
  epi_derivative_mat[, D_index[1]] <- pmax(DEdot, 0)
  
  epi_deriv <- as.vector(epi_derivative_mat)
  eps10 <- .Machine$double.eps * 1e12
  # epi_deriv[y < eps10] <- pmax(0, epi_deriv[y < eps10]) # exit wave was lost
  return(list(c(econ_derivs,epi_deriv)))
}




