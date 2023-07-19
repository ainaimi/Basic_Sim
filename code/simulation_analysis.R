pacman::p_load(
  here, tidyverse, mvtnorm
               )

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

# create a function to simulate and analyze data:
cluster_sim <- function(sim_id=1){
    
  n = 1e3
  p = 5
  set.seed(sim_id)
  
  ## CONFOUNDERS
  sigma <- matrix(0,nrow=p,ncol=p); diag(sigma) <- 1
  c     <- rmvnorm(n, mean=rep(0,p), sigma=sigma)
  
  # DESIGN MATRIX FOR THE OUTCOME MODEL
  muMatT <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]", collapse="+"),")")))
  parms3 <- rep(log(1.5),5)
  beta   <- c(log(2), -.5, parms3)
  
  # DESIGN MATRIX FOR THE PROPENSITY SCORE MODEL
  piMatT <- model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]", collapse="+"),")")))
  parms4 <- rep(log(1.5),5)
  theta  <- c(-.5,parms4) # assumed coefficients of propensity score model
  
  # PROPENSITY SCORE MODEL AND EXPOSURE
  ## pi is an actual parameter in base R, so good practice not to overwrite
  pi_x <- expit(piMatT%*%theta)
  x    <- rbinom(n,1,pi_x)

  # OUTCOME MODEL AND OUTCOME
  mu_y <- expit(cbind(x,muMatT)%*%beta)
  y <- rbinom(n, 1, mu_y)
  
  gen_dat <- data.frame(y,x,c) # create the dataset that we will analyze
  names(gen_dat) <- c("y","x",paste0("c",1:5))
  
  ##### 
  ##### end data gen
  ##### 
    
  # analysis start
  # true model
  
  mod_true <- glm(y ~ ., data = gen_dat, family = binomial(link = "logit"))
  est_logOR1 <- coef(mod_true)[2]
  est_logOR.SE1 <- summary(mod_true)$coefficients[2,2]
  
  # what happens if we put in all 2 way interactions?
  
  mod_int <- glm(y ~ .^2, data = gen_dat, family = binomial(link = "logit"))
  est_logOR2 <- coef(mod_int)[2]
  est_logOR.SE2 <- summary(mod_int)$coefficients[2,2]
  
  # what happens if we put in all 3 way interactions?
  
  mod_int2 <- glm(y ~ .^3, data = gen_dat, family = binomial(link = "logit"))
  est_logOR3 <- coef(mod_int2)[2]
  est_logOR.SE3 <- summary(mod_int2)$coefficients[2,2]
  
  # store the point estimates and standard error

  sim_dat <- data.frame(true_effect = log(2),
                        est_logOR1, est_logOR.SE1,
                        est_logOR2, est_logOR.SE2,
                        est_logOR3, est_logOR.SE3)
  
  return(sim_dat)
    
}
    
# let's do a simulation with a single simulation id
sim_res <- cluster_sim(sim_id = 1)

# now let's repeat for 100 different ids
sim_res <- lapply(1:100, function(x) cluster_sim(sim_id = x))

# combine results into a dataframe
sim_res <- do.call(rbind, sim_res)

# here are the simulation results, stored to a dataset
write_csv(sim_res, file = here("data","sim_res.csv"))