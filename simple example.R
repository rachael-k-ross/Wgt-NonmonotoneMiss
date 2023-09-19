
###########################################
#
# R code for simple illustrative example
# from Ross et al. 
# Accounting for Nonmonotone Missing Data using Inverse Probability Weights
# https://onlinelibrary.wiley.com/doi/10.1002/sim.9860
#
##########################################

library("tidyverse")
library("geepack")
library("R2jags")

# For CBE estimator by GIBBS, must download JAGS here: https://sourceforge.net/projects/mcmc-jags/

######################
# Generate full data 
######################
set.seed(13)

# Full data
gendata <- function(n){  
  tibble( 
    id = c(1:n),
    z = rbinom(n, size=1, prob=0.5), #binary confounder
    x = rbinom(n, size=1, prob=0.2 + 0.2*z), #binary exposure
    py = 1/(1+exp(-1*(log(0.1/0.9) + log(2)*z))), #probability of y under no exposure
    y0 = rbinom(n, size=1, prob=py), #binary potential outcome under no exposure
    y1 = rbinom(n, size=1, prob=py + .1), #binary potential outcomes under exposure
    y = x*y1 + (1-x)*y0) #observed binary outcome
}

full <- gendata(20000)
apply(full, 2, mean)

##########################
# Generate missingness   
##########################

add_missing <- function(data,parm){
  
  pR2 <-  with(data,plogis(parm[1] + parm[2]*z  + parm[3]*x  + parm[4]*y)) #probability of being in pattern 2
  pR3 <-  with(data,plogis(parm[5] + parm[6]*z  + parm[7]*x  + parm[8]*y)) #probability of being in pattern 3
  pR4 <-  with(data,plogis(parm[9] + parm[10]*z + parm[11]*x + parm[12]*y)) #probability of being in pattern 4
  pR1 <- 1 - pR2 - pR3 - pR4 #probability of being in pattern 1
  probs <- cbind(pR1,pR2,pR3,pR4)
  colnames(probs) <- NULL
  
  # Generate R
  data %>%
    mutate(R = Hmisc::rMultinom(probs,1), #multinomial draw using the probabilities for each pattern
           z = ifelse(R %in% c(3,4),NA,z),
           y = ifelse(R %in% c(2,3),NA,y),
           R1 = ifelse(R==1,1,0),
           R2 = ifelse(R==2,1,0),
           R3 = ifelse(R==3,1,0),
           R4 = ifelse(R==4,1,0)) 
}

# Set parameters of missingness models
g20 <- -0.945 # selected to produce 25% pattern 2
g21 <- log(.5)
g22 <- log(1.7)
g23 <- 0
g30 <- -1.815 # selected to produce 15% pattern 3
g31 <- 0
g32 <- log(1.3)
g33 <- 0
g40 <- -2.41 # selected to produce 10% pattern 4
g41 <- 0
g42 <- log(2)
g43 <- log(0.8)
gams <- c(g20, g21, g22, g23,
          g30, g31, g32, g33,
          g40, g41, g42, g43)

withmiss <- add_missing(full, gams)
prop.table(table(withmiss$R))


##########################
# Analysis of full data   
##########################

# Fit propensity score model
full_ps <- glm(x ~ z, data=full, family='binomial')
  
# Obtain IPTW
full_wt <- full %>%
  mutate(pscore = predict(full_ps, ., type="response"),
           ipw = x/pscore + (1-x)/(1-pscore))
  
# Fit weighted linear-binomial outcome model
full_out <- geeglm(y ~ x, family=binomial("identity"), data=full_wt,
                 weights=full_wt$ipw, id=id, corstr="independence")
  
# Results
results_full <- tibble(rd = coef(full_out)[[2]],
                       se = coef(summary(full_out))[[2,2]],
                       lcl = tidy(full_out, conf.int=T, exp = F)[[2,"conf.low"]],
                       ucl = tidy(full_out, conf.int=T, exp = F)[[2,"conf.high"]])
results_full

##########################
# Complete case analysis   
##########################

# Fit propensity score model
cc_ps <- glm(x ~ z, data=withmiss[withmiss$R==1,], family='binomial')

# Create IPTW
cc_wt <- withmiss %>%
  filter(R==1) %>%
  mutate(pscore = predict(cc_ps, ., type="response"),
         ipw = x/pscore + (1-x)/(1-pscore))

# Fit weighted linear-binomial outcome model
cc_out <- geeglm(y ~ x, family=binomial("identity"), data=cc_wt,
                  weights=cc_wt$ipw, id=id, corstr="independence")

# Results
results_cc <- tibble(rd = coef(cc_out)[[2]],
                       se = coef(summary(cc_out))[[2,2]],
                       lcl = tidy(cc_out, conf.int=T, exp = F)[[2,"conf.low"]],
                       ucl = tidy(cc_out, conf.int=T, exp = F)[[2,"conf.high"]])
results_cc

##########################
# Implement ST IPW UMLE   
##########################

# Create function of negative joint log-likelihood of missingness models
umleLogL <- function(g, data) {
  p2 <- with(data,plogis(g[1] + g[2]*z + g[3]*x         )) #logistic models for each pattern R>1
  p3 <- with(data,plogis(g[4]          + g[5]*x         ))
  p4 <- with(data,plogis(g[6]          + g[7]*x + g[8]*y))
  sump <- rowSums(cbind(p2,p3,p4),na.rm=T)
 
  ilogL <- ifelse(data$R==1,log(1-sump),
                  ifelse(data$R==2,log(p2),
                         ifelse(data$R==3,log(p3),log(p4))))
  neglogL <- -1*sum(ilogL)
  neglogL
}

# Minimize negative log-likelihood
init <- c(-2,0,0,-2,0,-2,0,0) #starting values
umlefit <- nlm(umleLogL, init, data=withmiss) 
gumle <- umlefit$est #umle gamma estimates

# Obtain marginal pr R=1
mnum <- mean(withmiss$R1)

# Obtain probability of complete case
cc_umle <- withmiss %>%
  filter(R==1) %>%
  mutate(p2 = plogis(gumle[1] + gumle[2]*z + gumle[3]*x             ),
         p3 = plogis(gumle[4]              + gumle[5]*x             ),
         p4 = plogis(gumle[6]              + gumle[7]*x + gumle[8]*y),
         p1 = 1 - p2 - p3 - p4,
         ipmw = mnum/p1)

summary(cc_umle$ipmw)
sum(cc_umle$ipmw)

# Fit weighted propensity score model
upsmodel <- glm(x ~ z, data=cc_umle, family='binomial', weights = cc_umle$ipmw)

# Create iptw and combine weights
uwithweight <- cc_umle %>%
  mutate(pscore = predict(upsmodel, ., type="response"),
         iptw = x/pscore + (1-x)/(1-pscore),
         ipw = ipmw*iptw)


# Fit weighted linear-binomial outcome model
uoutmod <- geeglm(y ~ x, family=binomial("identity"), data=uwithweight, 
                 weights=uwithweight$ipw, id=id, corstr="independence")

# Results
results_umle <- tibble(rd = coef(uoutmod)[[2]],
                     se = coef(summary(uoutmod))[[2,2]],
                     lcl = tidy(uoutmod, conf.int=T, exp = F)[[2,"conf.low"]],
                     ucl = tidy(uoutmod, conf.int=T, exp = F)[[2,"conf.high"]])
results_umle


################################################################
# Implement ST IPW CBE by adaptive Gibbs sampling using R2jags   
################################################################

# standardize the data
scale_rmna <- function(data,x){ #function ignores missing data in standardizing
  mean <- mean(data[[x]],na.rm=TRUE)
  sd <- sd(data[[x]], na.rm=TRUE)
  (data[[x]]-mean)/sd
}

scaled <- withmiss %>%
  mutate(x = scale_rmna(.,"x"),
         y = scale_rmna(.,"y"),
         z = scale_rmna(.,"z"),
         R = as.numeric(R))

# Prepare data for JAGS
sorted <- as_tibble(scaled) %>% arrange(R) # Sort so that complete cases (R=1) are first
dat = as.list(sorted[,("R")])
dat$L = as.matrix(sorted[,c("z","x","y")]) 
dat$L[is.na(dat$L)] = -9999 # Replace NA to Inf
dat$N = length(dat$R)
dat$f = rep(1, dat$N) # Vector of 1s length N
dat$Nc = sum(dat$R==1) # Number of complete cases
dat$onesc = rep(1, dat$Nc) # Vector of 1s length of complete cases
dat$c = 10^-8 # As recommended by ST

# Function for random starting values
initialvals <- function(numcoeff){
  ints <- runif(3, min=-4,max=-1)
  lim <- .06
  c(ints[1],runif(numcoeff[1], min=-lim,max=lim),
    ints[2],runif(numcoeff[2], min=-lim,max=lim),
    ints[3],runif(numcoeff[3], min=-lim,max=lim))
}


# JAGS function
jmod <- function(){
  for(i in 1:N){
    f[i] ~ dbern(pmiss[i,R[i]]) # f = 1 for all obs 
    #
    logit(pmiss[i, 2]) <- g[1] + g[2]*L[i,1] + g[3]*L[i,2]         
    logit(pmiss[i, 3]) <- g[4]               + g[5]*L[i,2]           
    logit(pmiss[i, 4]) <- g[6]               + g[7]*L[i,2] + g[8]*L[i,3] 
    #
    pmiss[i,1] <- 1 - pmiss[i,2] - pmiss[i,3] - pmiss[i,4] 
  }
  # constraint
  for (j in 1:Nc){
    onesc[j] ~ dbern(C[j]) 
    C[j] <- step(pmiss[j,1]-c)
  }
  # priors
  for(k in 1:8){
    g[k] ~ dnorm(0,1/100) # diffuse prior as recommended by ST
  }
}
  
# Initial values for 3 chains
init <- list(list(g=initialvals(c(2,1,2))),
             list(g=initialvals(c(2,1,2))),
             list(g=initialvals(c(2,1,2))))
  
# Run jags
jagsfit <- R2jags::jags(data=dat, inits = init, n.chains = 3,
                          parameters.to.save='g',
                          n.iter=1000, n.burnin=500, n.thin=1,
                          model.file=jmod)

jagsfit 
  
# Store parameter estimates (medians)
gcbegibbs <- jagsfit$BUGSoutput$median$g
  
# Obtain marginal Pr(R=1)
mnum <- mean(scaled$R1)
  
# Obtain probability of complete case
cc_scaledgibbs <- scaled %>%
    filter(R==1) %>%
    mutate(p2 = plogis(gcbegibbs[1] + gcbegibbs[2]*z + gcbegibbs[3]*x                 ),
           p3 = plogis(gcbegibbs[4]                  + gcbegibbs[5]*x                 ),
           p4 = plogis(gcbegibbs[6]                  + gcbegibbs[7]*x + gcbegibbs[8]*y),
           p1 = 1 - p2 - p3 - p4,
           ipmw = mnum/p1)
  
# Merge with unscaled data
cc_gibbs <- cc_scaledgibbs %>%
  dplyr::select(id, ipmw) %>%
  left_join(withmiss, by = "id") 

summary(cc_gibbs$ipmw)
sum(cc_gibbs$ipmw)
  
# Fit weighted propensity score model
cgpsmodel <- glm(x ~ z, data=cc_gibbs, family='binomial', weights = cc_gibbs$ipmw)
  
# Create iptw and combine weights
cgwithweight <- cc_gibbs %>%
    mutate(pscore = predict(cgpsmodel, ., type="response"),
           iptw = x/pscore + (1-x)/(1-pscore),
           ipw = ipmw*iptw)
  
  
# Fit weighted linear-binomial outcome model
cgoutmod <- geeglm(y ~ x, family=binomial("identity"), data=cgwithweight, 
                    weights=cgwithweight$ipw, id=id, corstr="independence")
  
# Results
results_cbegibbs <- tibble(rd = coef(cgoutmod)[[2]],
                         se = coef(summary(cgoutmod))[[2,2]],
                         lcl = tidy(cgoutmod, conf.int=T, exp = F)[[2,"conf.low"]],
                         ucl = tidy(cgoutmod, conf.int=T, exp = F)[[2,"conf.high"]])
results_cbegibbs


########################################
# Implement ST IPW CBE by M-H Algorithm   
########################################

# standardize the data
scale_rmna <- function(data,x){ #function ignores missing data in standardizing
  mean <- mean(data[[x]],na.rm=TRUE)
  sd <- sd(data[[x]], na.rm=TRUE)
  (data[[x]]-mean)/sd
}

scaled <- withmiss %>%
  mutate(x = scale_rmna(.,"x"),
         y = scale_rmna(.,"y"),
         z = scale_rmna(.,"z"),
         R = as.numeric(R))

### Create functions to be used in MH
# Log of zero will resolve to -Inf
log_zero = function(x) ifelse(x==0, -Inf, log(x))

# Log likelihood fx
loglikf = function(data,g,c){
  p2 <- with(data,plogis(g[1] + g[2]*z + g[3]*x         ))
  p3 <- with(data,plogis(g[4]          + g[5]*x         ))
  p4 <- with(data,plogis(g[6]          + g[7]*x + g[8]*y))
  sump <- rowSums(cbind(p2,p3,p4),na.rm=TRUE)
  
  constrmet <- ifelse(sump < 1 - c,1,0) # check that constraint is met
  
  ilogL <- ifelse(data$R==1,log_zero((1-sump)*constrmet), 
                  ifelse(data$R==2,log(p2),
                         ifelse(data$R==3,log(p3),log(p4))))
  
  sum(ilogL)
} 

# Priors
priorf = function(parm){
  # diffuse priors var=100 for all parameters
  sum(dnorm(parm, mean = 0, sd=sqrt(100), log=TRUE))
}

### Metropolis-Hastings Algorithm
parnum <- 8
m <- 10000 # iterations
b <- 5000 # burn-in
c <- 1e-8
sigma = 0.015 # standard deviation of the proposal distribution  

# Make empty shells to store values
chain = matrix(nrow=m, ncol=parnum)
colnames(chain) = c('g20', 'g21', 'g22', 'g30', 'g32', 'g40', 'g42', 'g43')
accept = matrix(nrow=m, ncol=1)
  
# Starting values
chain[1,] = c(-2, 0, 0, -2, 0, -2, 0, 0)   
accept[1,1] = 1 # set first row to accept
  
# random walk
for(k in 2:m){
    oldp = chain[k-1,] # current parameter values
    old = loglikf(scaled, oldp, c) # log likelihood at current parameter values
    
    prop = rnorm(n = parnum, mean = 0, sd = sigma) # draws from proposal distribution 
    newp = oldp + prop # proposed parameter values
    new = loglikf(scaled, newp, c) # log likelihood at proposed values
    
    # Compare old and new log likelihoods, toss coin weighted by acceptance prob
    num = new + priorf(newp)
    den = old + priorf(oldp)
    acceptprob = exp(num-den)
    acc = (acceptprob > runif(1))
    
    if(acc){
      chain[k,] = newp # if accept, save proposed parameter values
      accept[k] = 1
    }else{
      chain[k,] = oldp # if reject, save previous parameter values
      accept[k] = 0
    }
    
    if(k%%500==0) cat(paste("iter",k))
}

mean(accept)

post <- chain[(b+1):m,]
# Store parameter estimates (medians)
gcbemh <- apply(post, 2, median)

# Obtain marginal Pr(R=1)
mnum <- mean(scaled$R1)

# Obtain probability of complete case
cc_scaledmh <- scaled %>%
  filter(R==1) %>%
  mutate(p2 = plogis(gcbemh[1] + gcbemh[2]*z + gcbemh[3]*x               ),
         p3 = plogis(gcbemh[4]                + gcbemh[5]*x              ),
         p4 = plogis(gcbemh[6]                + gcbemh[7]*x + gcbemh[8]*y),
         p1 = 1 - p2 - p3 - p4,
         ipmw = mnum/p1)

# Merge with unscaled data
cc_mh <- cc_scaledmh %>%
  dplyr::select(id, ipmw) %>%
  left_join(withmiss, by = "id") 

summary(cc_mh$ipmw)
sum(cc_mh$ipmw)

# Fit weighted propensity score model
cmhpsmodel <- glm(x ~ z, data=cc_mh, family='binomial', weights = cc_mh$ipmw)

# Create iptw and combine weights
cmhwithweight <- cc_mh %>%
  mutate(pscore = predict(cmhpsmodel, ., type="response"),
         iptw = x/pscore + (1-x)/(1-pscore),
         ipw = ipmw*iptw)


# Fit weighted linear-binomial outcome model
cmgoutmod <- geeglm(y ~ x, family=binomial("identity"), data=cmhwithweight, 
                   weights=cmhwithweight$ipw, id=id, corstr="independence")

# Results
results_cbemh <- tibble(rd = coef(cmgoutmod)[[2]],
                           se = coef(summary(cmgoutmod))[[2,2]],
                           lcl = tidy(cmgoutmod, conf.int=T, exp = F)[[2,"conf.low"]],
                           ucl = tidy(cmgoutmod, conf.int=T, exp = F)[[2,"conf.high"]])
results_cbemh


########################################
# Compare results   
########################################

mean(full$y1) - mean(full$y0)
results_full
results_cc
results_umle
results_cbegibbs
results_cbemh
