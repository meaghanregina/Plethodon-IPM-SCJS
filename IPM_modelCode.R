
setwd("C:/Users/gade.15/Box/Research/Mark_Recapture/IPM")
load("IPM_data_30Sept2020.RData")

ya <- countA_array
yj <- countJ_array

nsite <- dim(ya)[1]
nyear_obs <- dim(ya)[2]
nreps <- dim(ya)[3]

site_obs <- which(!(is.na(apply(ya, 1, mean, na.rm=T)))) # the sites with count
nsite_obs <- length(site_obs)

his <- CH_array
nind <- dim(his)[1]
nyear_cap <- dim(his)[2]
ncaps <- dim(his)[3]

site_cap <- as.numeric(sort(names(table(CH_array)))) # the sites with mrc
site_cap <- site_cap[which(site_cap != nsite+1)]
nsite_cap <- length(site_cap)
site_cap01 <- ifelse(1:nsite %in% site_cap, 1, 0)

elev <- (elev - mean(elev)) / sd(elev)
stream <- (stream - mean(stream)) / sd(stream)

## add offset
offset <- readxl::read_excel("data/combiningPlots/offsets.xlsx")

# Data0 bundle ------------------------------------------------------------
data0 = list(
  nsite = nsite, 
  nsite_obs = nsite_obs, 
  site_obs = site_obs, 
  nsite_cap = nsite_cap, 
  site_cap01 = site_cap01, 
  nyear_obs = nyear_obs, 
  nyear_cap = nyear_cap, 
  nreps = nreps, 
  nind = nind, 
  ncaps = ncaps, 
  ya = ya,
  yj = yj,
  his = his,
  elev = elev,
  stream = stream,
  dist = dist_mat,
  first = first, # year of first capture
  zf = zf # plot of first capture
)

# Jags Code ---------------------------------------------------------------
## obs is count, cap = MRC
writeLines(
  #"
  model{
# priors
lambda0_alpha ~ dnorm(0, .01)
lambda0_mean <- exp(lambda0_alpha) # mean initial population per site
lambda0_elev ~ dnorm(0, .01) #effect of elevation on initial population size
lambda0_stream ~ dnorm(0, .01)
lambda0_elev_stream ~ dnorm(0, .01)
lambda0_tau ~ dgamma(.01, .01)
lambda0_sd <- 1 / sqrt(lambda0_tau)
gamma_alpha ~ dnorm(0, .01)
gamma_mean <- exp(gamma_alpha) # mean per capita reproductionn
gamma_ddp ~ dnorm(0, .01) # density dependence of reproduction
gamma_elev ~ dnorm(0, .01) # effect of elevation on reproduction
gamma_stream ~ dnorm(0, .01)
gamma_elev_stream ~ dnorm(0, .01)
gamma_tau ~ dgamma(.01, .01)
gamma_sd <- 1 / sqrt(gamma_tau)
omega_alpha ~ dnorm(0, .01)
omega_mean <- ilogit(omega_alpha) #mean survival rate
omega_ddp ~ dnorm(0, .01) #density dependence of survival
omega_elev ~ dnorm(0, .01) #effect of elevation on survival 
omega_stream ~ dnorm(0, .01)
omega_elev_stream ~ dnorm(0, .01)
omega_tau ~ dgamma(.01, .01)
omega_sd <- 1 / sqrt(omega_tau)
kappa ~ dunif(0, 1) #emigration rate
theta_mean ~ dgamma(.01, .01) # mean decay parameter for dispersal 
theta_tau ~ dgamma(.01, .01)
theta_sd <- 1 / sqrt(theta_tau)
logit_pobs_mean ~ dnorm(0, .01)
pobs_mean <- ilogit(logit_pobs_mean) # mean detection probability
logit_pobs_tau ~ dgamma(.01, .01)
logit_pobs_sd <- sqrt(1 / logit_pobs_tau)
logit_pcap_mean ~ dnorm(0, .01)
pcap_mean <- ilogit(logit_pcap_mean) #mean capture probability 
logit_pcap_tau ~ dgamma(.01, .01)
logit_pcap_sd <- sqrt(1 / logit_pcap_tau)

for (i in 1:nsite) {
  for (j in 1:nsite) {
    w_ori_epsilon[i,j] ~ dnorm(0, theta_tau)
    w_ori[i,j] <- ifelse(i==j, 1e-6, exp(-1 * theta_mean * dist[i,j] + w_ori_epsilon[i,j]))
    w[i,j] <- w_ori[i,j] / sum(w_ori[i,1:nsite])
  } # j
} # i

for (t in 1:nyear_obs) {
  logit_pobs[t] ~ dnorm(logit_pobs_mean, logit_pobs_tau)
  pobs[t] <- ilogit(logit_pobs[t])
} # t

for (t in 1:nyear_cap) {
  logit_pcap[t] ~ dnorm(logit_pcap_mean, logit_pcap_tau)
  pcap[t] <- ilogit(logit_pcap[t])
  for (i in 1:nsite) {
    pcap_mat[i,t] <- ifelse(site_cap01[i]==0, 0, pcap[t])
  } # i
} # t

# likelihood
### population count
for (i in 1:nsite) {
  lambda0_epsilon[i] ~ dnorm(0, lambda0_tau)
  lambda0[i] <- exp(lambda0_alpha + # rate of initial pop size
                    lambda0_elev * elev[i] + 
                    lambda0_stream * stream[i] + 
                    lambda0_elev_stream * elev[i] * stream[i] + 
                    lambda0_epsilon[i]) 
  
  N[i,1] ~ dpois(lambda0[i]) # Initial local populations size
  
  for (t in 2:nyear_obs) {
    gamma_epsilon[i,t-1] ~ dnorm(0, gamma_tau)
    gamma[i,t-1] <- exp( ## per capita reproduction rate
      gamma_alpha + 
        gamma_ddp * (N[i,t-1] - lambda0[i]) / lambda0[i] + 
        gamma_elev * elev[i] +
        gamma_stream * stream[i] +
        gamma_elev_stream * elev[i] * stream[i] +
        gamma_epsilon[i,t-1])
    R[i,t-1] ~ dpois(N[i,t-1] * gamma[i,t-1]) # number of reproduced individuals
    
    omega_epsilon[i,t-1] ~ dnorm(0, omega_tau)
    omega[i,t-1] <- ilogit( #probability of survival
      omega_alpha + #survival
        omega_ddp * (N[i,t-1] - lambda0[i]) / lambda0[i] + #density dep. of survival
        omega_elev * elev[i] +
        omega_stream * stream[i] +
        omega_elev_stream * elev[i] * stream[i] +
        omega_epsilon[i,t-1])
    S[i,t-1] ~ dbin(omega[i,t-1], N[i,t-1]) #number of survived individuals
    
    E[i,t-1] ~ dbin(kappa, S[i,t-1]) # number f emigrants
    for (j in 1:nsite) {
      M[i,j,t-1] <- E[i,t-1] * w[i,j] #number of moved individuals
    } # j
    I[i,t-1] ~ dpois(sum(M[1:nsite,i,t-1])) #number of Immigrants
    
    N[i,t] <- R[i,t-1] + S[i,t-1] - E[i,t-1] + I[i,t-1] #Local pop size
  } # t
} # i

### capture-recapture
for (t in 2:nyear_cap) { #nyear
  for (i in 1:nsite) {
    for (j in 1:nsite) {
      phi[i,j,t-1] <- ifelse(i==j, omega[i,t-1]*(1-kappa), omega[i,t-1]*kappa*w[i,j])
    } # j
    phi[i,nsite+1,t-1] <- 1 - omega[i,t-1]
  } # i
  for (i in 1:nsite) {
    phi[nsite+1,i,t-1] <- 0
  } # i
  phi[nsite+1,nsite+1,t-1] <- 1
} # t

for (i in 1:nind) {
  z[i,first[i]] <- zf[i]
  for (t in (first[i]+1):nyear_cap) {
    z[i,t] ~ dcat(phi[z[i,t-1],,t-1])
  } # t
} # i

# observation
### population count
for (i in 1:nsite_obs) { #num of Count sites
  for (t in 1:nyear_obs) { #num years of Counts =3
    for (j in 1:nreps) { #replicates per year (4, 1, 1) ## orig was nreps
      ya[i,t,j] ~ dbin(pobs[t], N[site_obs[i],t])
    } # j
  } # t
} # i

for (i in 1:nsite_obs) {
  for (t in 1:(nyear_obs-1)) {
    for (j in 1:nreps) {
      yj[i,t,j] ~ dbin(pobs[t+1], R[site_obs[i],t])
    } # j
  } # t
} # i

### encounter history
for (t in 1:nyear_cap) {
  for (i in 1:nsite) {
    for (j in 1:nsite) {
      po[i,j,t] <- ifelse(i==j, pcap_mat[i,t], 0)
    } # j
    po[i,nsite+1,t] <- 1 - pcap_mat[i,t]
  } # i
  for (i in 1:nsite) {
    po[nsite+1,i,t] <- 0
  } # i
  po[nsite+1,nsite+1,t] <- 1
} # t

for (i in 1:nind) {
  for (t in (first[i]+1):nyear_cap) { 
    for (j in 1:ncaps) {
      his[i,t,j] ~ dcat(po[z[i,t],1:(nsite+1),t])
    } # j
  } # t
} # i

}
", con = "IPM_v2_int.txt")

# Model  Specifications ---------------------------------------------------

# #Initial values ---------------------------------------------------------
Ni <- matrix(max(countA_array,na.rm=T)*20+20, 50, 4)  
Si <- Ri <- matrix(max(countA_array,na.rm=T)*10+10, 50, 4-1) 
Ei <- Ii <- matrix(max(countA_array,na.rm=T)*5+5, 50, 4-1) 
Ni[,-1] <- NA

ch_mat <- matrix(, nind, nyear)
for (i in 1:nind) {
  for (t in 1:nyear) {
    tt <- his[i,t,]
    tt[which(tt==nsite+1)] <- NA
    if (length(which(!(is.na(tt))))==0) {
      ch_mat[i,t] <- nsite+1
    } else {
      ch_mat[i,t] <- mean(tt, na.rm=T)
    }
  }
}

zi <- matrix(,nind, nyear)
for (i in 1:nind) {
  tt <- ch_mat[i,]
  tt[which(tt == nsite+1)] <- zf[i] #plot of first capture
  zi[i, (first[i]+1):nyear] <- tt[(first[i]+1):nyear]
} # i

inits <- function() {list(
  lambda0_alpha=0, lambda0_elev=0, lambda0_stream=0, lambda0_elev_stream=0, lambda0_tau=1, 
  gamma_alpha=0, gamma_ddp=0, gamma_elev=0, gamma_stream=0, gamma_elev_stream=0, gamma_tau=1, 
  omega_alpha=0, omega_ddp=0, omega_elev=0, omega_stream=0, omega_elev_stream=0, omega_tau=1, 
  kappa=runif(1,0,1), theta_mean=1, theta_tau=1, 
  R=Ri, S=Si, E=Ei, I=Ii, N=Ni, 
  z=zi)}

# Params to Monitor
params0 <- c(
  'N', 'R', 'S', 'E', 'I', 'M', 
  'lambda0_mean', 'lambda0_elev', 'lambda0_stream', 'lambda0_elev_stream', 'lambda0_sd', 
  'gamma_mean', 'gamma_ddp', 'gamma_elev', 'gamma_stream', 'gamma_elev_stream', 'gamma_sd', 
  'omega_mean', 'omega_ddp', 'omega_elev', 'omega_stream', 'omega_elev_stream', 'omega_sd', 
  'kappa', 'theta_mean', 'theta_sd', 
  'pobs_mean', 'logit_pobs_sd', 
  'pcap_mean', 'logit_pcap_sd')

# Run Model ---------------------------------------------------------------
nc <- 2
nt <- 1
na <- 100
nb <- 100
ni <- 200

mod2 <- jagsUI::jags(
  data0,
  inits = inits, 
  parameters.to.save = params0,
  model.file = "IPM_v2_int.txt",
  parallel = T,
  n.cores = 15,
  n.adapt = 2000,
  n.chains = 15,
  n.thin = 15,
  n.iter = 35000,
  n.burnin = 2500
  ,codaOnly = c("S", "E", "M", "I", "N","R")
)

print(mod)
