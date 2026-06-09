## Odin DSL: integrated epi-econ model, Model 2 (time-varying productivity)
##
## Extends Model 1: lambda (productivity) is a state variable that tracks output growth.
## State variables: H_h, lambda (econ), then S, E, Iu, Id, C, R (epi, each dim 4)
## Column order in output: t, H_h, lambda, S[1:4], E[1:4], Iu[1:4], Id[1:4], C[1:4], R[1:4]

## User parameters: econ #############################

alpha0    <- user()
alpha1    <- user()
alpha2    <- user()
theta     <- user()
G         <- user()
cons0     <- user()
emp0      <- user()
lf        <- user()
lambda_p0 <- user()
lambda_p1 <- user()

## User parameters: epi-econ coupling #############################

q1        <- user()
q2        <- user()
integrate <- user(0)

## User parameters: epidemic #############################

beta          <- user()
TEtoI         <- user()
TItoR         <- user()
TItoC         <- user()
TCtoR         <- user()
prob_detected <- user()

## User parameters: population and contact matrices #############################

NNs[] <- user()
dim(NNs) <- 4

Mww[,] <- user()
dim(Mww) <- c(4, 4)
Mcw[,] <- user()
dim(Mcw) <- c(4, 4)
Mcc[,] <- user()
dim(Mcc) <- c(4, 4)
Mcom[,] <- user()
dim(Mcom) <- c(4, 4)

## Initial conditions #############################

H_h_init    <- user()
lambda_init <- user()
S_init[]    <- user()
dim(S_init) <- 4
E_init[]    <- user()
dim(E_init) <- 4

initial(H_h)    <- H_h_init
initial(lambda) <- lambda_init
initial(S[])    <- S_init[i]
initial(E[])    <- E_init[i]
initial(Iu[])   <- 0
initial(Id[])   <- 0
initial(C[])    <- 0
initial(R[])    <- 0

dim(S)  <- 4
dim(E)  <- 4
dim(Iu) <- 4
dim(Id) <- 4
dim(C)  <- 4
dim(R)  <- 4

## Block 1: Epi -> Econ #############################

total_C       <- sum(C[])
notifications <- total_C / TCtoR
scalar        <- if (integrate == 1) q1 + (1 - q1) / (1 + q2 * notifications) else 1
available_lf        <- if (integrate == 1) lf - C[1] / 1e6 else lf

## Block 2: Econ -> Epi  #############################

# need consumption first
alpha1t <- scalar * alpha1
alpha2t <- scalar * alpha2
denom   <- 1 - alpha1t * (1 - theta)
cons_d  <- (alpha0 + alpha1t * G * (1 - theta) + alpha2t * H_h) / denom
cons_s  <- max(0, available_lf * lambda - G)
cons    <- min(cons_s, cons_d)

relative_consumption <- if (integrate == 1) cons / cons0 else 1
relative_work        <- if (integrate == 1) (cons + G) / lambda / emp0 else 1

## Block 3: Epi derivatives #############################
# Contact matrix and force of infection

contact[,] <- Mcom[i, j] +
  relative_work ^ 2 * Mww[i, j] +
  relative_work * relative_consumption * Mcw[i, j] +
  relative_consumption ^ 2 * Mcc[i, j]
dim(contact) <- c(4, 4)

I_total[]        <- Id[i] + Iu[i]
dim(I_total)     <- 4

foi_contrib[,]   <- contact[i, j] * I_total[j] / NNs[j]
dim(foi_contrib) <- c(4, 4)

foi[]    <- beta * sum(foi_contrib[i, ])
dim(foi) <- 4

## Transition rates #############################

latent_det[]    <- E[i] * prob_detected / TEtoI
latent_undet[]  <- E[i] * (1 - prob_detected) / TEtoI
undet_recover[] <- Iu[i] / TItoR
det_detected[]  <- Id[i] / TItoC
con_recover[]   <- C[i] / TCtoR

dim(latent_det)    <- 4
dim(latent_undet)  <- 4
dim(undet_recover) <- 4
dim(det_detected)  <- 4
dim(con_recover)   <- 4

dot_S[]  <- -S[i] * foi[i]
dot_E[]  <- S[i] * foi[i] - latent_det[i] - latent_undet[i]
dot_Iu[] <- latent_undet[i] - undet_recover[i]
dot_Id[] <- latent_det[i] - det_detected[i]
dot_C[]  <- det_detected[i] - con_recover[i]
dot_R[]  <- undet_recover[i] + con_recover[i]

dim(dot_S)    <- 4
dim(dot_E)  <- 4
dim(dot_Iu) <- 4
dim(dot_Id)  <- 4
dim(dot_C)   <- 4
dim(dot_R)   <- 4

## Block 4: Econ derivatives #############################

# Econ (Model 2: time-varying productivity)
Y       <- cons + G
YD      <- Y * (1 - theta)
dot_H_h <- YD - cons

## derivative of scalar with respect to time (needed for dot_lambda)
dot_total_C       <- sum(det_detected[]) - sum(con_recover[])
dot_notifications <- dot_total_C / TCtoR
dot_scalar        <- if (integrate == 1) -(1 - q1) * q2 * dot_notifications / (1 + q2 * notifications) ^ 2 else 0

# derivative of consumption when demand-determined (quotient rule)
f_term     <- alpha0 + alpha1t * G * (1 - theta) + alpha2t * H_h
f_prime    <- dot_scalar * (alpha1 * G * (1 - theta) + alpha2 * H_h) + alpha2t * dot_H_h
g_term     <- 1 - alpha1t * (1 - theta)
g_prime    <- -dot_scalar * alpha1 * (1 - theta)
dot_cons_d <- (g_term * f_prime - f_term * g_prime) / g_term ^ 2

# derivative of consumption when supply-determined: cons = available_lf * lambda - G
# substituting into dot_lambda and solving implicitly gives:
dot_lf     <- -(det_detected[1] - con_recover[1]) / 1e6
dot_cons_s <- Y * lambda * (dot_lf + lambda_p0 * available_lf) / (Y - lambda * lambda_p1 * available_lf)

# consumption derivative depends on whether consumption is demand or supply determined
dot_cons <- if (cons_s < cons_d) dot_cons_s else dot_cons_d

dot_lambda <- lambda * (lambda_p0 + lambda_p1 / Y * dot_cons)


## Odin outputs ############################################

deriv(H_h) <- dot_H_h
deriv(lambda) <- dot_lambda

deriv(S[])  <- dot_S[i]
deriv(E[])  <- dot_E[i]
deriv(Iu[]) <- dot_Iu[i]
deriv(Id[]) <- dot_Id[i]
deriv(C[])  <- dot_C[i]
deriv(R[])  <- dot_R[i]

## Outputs for plotting #############################

output(cons)   <- cons
output(cons_d) <- cons_d
output(cons_s) <- cons_s
output(GDP)    <- Y
