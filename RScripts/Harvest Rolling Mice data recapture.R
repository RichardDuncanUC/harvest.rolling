
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(lubridate)
  library(jagsUI)

  remove(list = ls())
  unlogit <- function(x) exp(x) / (1 + exp(x))

# Fate 2 = means dead (not counted in pre-baiting as 'not part of "baited  population", none post-treat on control sites, some on post-cotol treated sites but ignored - died as part on experimental treatment.)
# Fate 3 = Escaped but alive so "1" really  
# Fate 4 = taken to lab (only on last day so treat as 1 for analysis) 
  
################################################################################
# read in trapping night data
  tn <- read_excel("data/Parkes Trapping Nights.xlsx", trim_ws = TRUE) %>%
        mutate(grid = as.numeric(factor(DataSiteName)),
                trip = as.numeric(factor(Trip)))
  glimpse(tn)

################################################################################
# read in capture data
  dat <- read_excel("data/CMR Parkes 19.03.21.xlsx", trim_ws = TRUE)
  glimpse(dat)

# remove dead individuals - only use Fate == 1, 3 or 4
  dat <- dat %>%
    filter(Fate != 2)
  
  table(dat$DataSiteName, dat$CaptureDate)
  
  dat$PitTag <- toupper(dat$PitTag)
  
  dat <- dat %>%
         mutate(grid = as.numeric(factor(DataSiteName)),
                date = CaptureDate,
                trip = as.numeric(factor(Trip))) %>%
    filter(DataType == 1)   #  I have changed some in excel worsheet if datatype=2 but is a full trapping night.
  
  table(dat$DataSiteName, dat$SessionName)
  table(dat$grid, dat$trip)
  
# set missing pit tags to unique individuals
  n_miss <- sum(is.na(dat$PitTag))
  table(is.na(dat$PitTag), dat$trip)
  
  dat$PitTag[is.na(dat$PitTag)] <- paste0("m", 1:n_miss)
  
# convert dates to days of trapping   ##  does not include days of 0 caught in the right Day #.   JLA TG Session 2 - Day 2 of 5 was 0 captures.  
  dat <- dat %>%
    group_by(SessionName, DataSiteName, grid, trip) %>%
    mutate(day = as.numeric(factor(rank(date))))
  
  table(dat$date, dat$day)
  table(dat$day, dat$SessionName, dat$DataSiteName)
  table(dat$day, dat$trip, dat$grid)

################################################################################    
# number of captures by day
  dd <- dat %>%
    group_by(DataSiteName, trip, day) %>%
    summarise(n = n())
  
  ggplot(dd, aes(y = n, x = day, colour = factor(DataSiteName))) +
    geom_point(size = 3) +
    geom_line() +
    facet_wrap(~ trip) +
    theme_bw()
  
# not an obvious tendency for capture rates to increase/decrease over time
# so could leave this out
  
  
################################################################################
# for each individual, get the days it was captured on each grid on each trip
# number individuals on each trip and grid
  ind <- dat %>%
         group_by(DataSiteName, trip, grid, day, PitTag) %>%
         summarise(n = n()) %>%
         mutate(n = ifelse(n > 1, 1, n)) %>%
         group_by(DataSiteName, trip, grid, PitTag) %>%
         summarise(n = n())
    
  ind
  
################################################################################
# the minimum number of animals caught on each grid on each trip
  min.mice <- ind %>%
              ungroup() %>%
              group_by(DataSiteName, grid, trip) %>%
              summarise(n = n()) 
  
# list of grids and trips from tn
  all <- tn %>%
    dplyr::select(DataSiteName, grid, trip, NightsTrapped)
  
  min.mice <- left_join(all, min.mice) %>%
    filter(NightsTrapped > 0)
  
  ggplot(min.mice, aes(y = n, x = trip)) +
    geom_point(size = 3) +
    geom_line() +
    facet_wrap(~DataSiteName) +
    theme_bw() 
    
################################################################################
# use the data frame ind to create an array with the number of captures for each individual on each grid

  n.grid <- max(tn$grid)
  n.trip <- max(tn$trip)
  n.day <- max(tn$NightsTrapped)
  
# set up array Y with rows i = indiviudals, columns = j = trip, k = grids
  Y <- array(0, dim = c(1000, n.trip, n.grid))
  z <- array(NA, dim = c(1000, n.trip, n.grid))

# number of nights trapping
  J <- array(0, dim = c(n.trip, n.grid))
  
    
  for(j in 1:n.trip) {
    for(k in 1:n.grid) {
    # mice caught on grid k during trip j 
      cap <- ind$n[ind$trip == j & ind$grid == k]
      if(length(cap) > 0) {
        n.cap <- length(cap)
        Y[1:n.cap, j, k] <- cap
        z[1:n.cap, j, k] <- 1
      }
      
      night <- tn$NightsTrapped[tn$trip == j & tn$grid == k]
      J[j, k] <- night
    }
  }
  
################################################################################
# row i is individual animal
# column j is trip
# sheet k is grid

# variables we need for the model

# Y = number of captures for each mouse (i) caught on each trip (j) on each grid (k) on each day (m)
# z = latent variable for whether pseudo mice were present or not

  N_grid <- dim(Y)[3]
  N_animals <- dim(Y)[1]
  N_trip <- dim(Y)[2]

################################################################################
# JAGS model
  
  # JAGS model
  
  CR_mice <- "model {

  for(k in 1:N_grid){
    for(j in 1:N_trip){
  		for(i in 1:N_animals)  {
  			Y[i, j, k] ~ dbin(mu[i, j, k], J[j, k])
		  	mu[i, j, k] <- z[i, j, k] * p[i, j, k]
	  		z[i, j, k] ~ dbern(psi[j, k])
	  		logit(p[i, j, k]) <- lp[i, j, k]
	  	# model individual-level heterogeneity in capture probability 	
	  		lp[i, j, k] ~ dnorm(mu.ind, tau.ind)
		  } #i

    # prior for presence probability differs betwen trips and grids
    # use scale prior as in Link (2013)
      psi[j, k] ~ dbeta(0.001, 1)

    # Derived parameters: the estimated number of individuals on each grid at each trip
      N[j, k] <- sum(z[, j, k])
      logN[j, k] <- log(N[j, k])
      
    }  #j
    
  } #k
  
  # priors
  
  mu.ind ~ dnorm(0, 0.001)
  tau.ind <- 1 / (sigma.ind * sigma.ind)
  sigma.ind ~ dunif(0, 10)

  }"  #model finish
  
  #model finish
  
  #RUN MODEL------------------------------------------------
  
  write(CR_mice, "CR_mice.txt")
  
  mod <- jags(model.file = "CR_mice.txt",
              data = list(Y = Y, J = J, z = z, N_trip = N_trip, N_grid = N_grid, N_animals = N_animals),
              parameters.to.save = c("mu.ind", "sigma.ind", "N", "logN"),
              n.chains = 3,
              n.iter = 15000,
              n.burnin = 5000,
              n.thin = 1,
              parallel = T)
  
  mod.sum <- mod$summary
  mod.sum[, c(1, 2, 3, 7, 8, 9)]
  
  
  N <- mod$sims.list$N
  
# where J = 0 set N = NA
  for(j in 1:n.trip) {
    for(k in 1:n.grid) {
      if(J[j, k] == 0) N[1:dim(N)[1], j, k] <- NA
    }
  }
  
  head(N)

# save output for later use  
  save.image("data/N estimates.RData")


# for each grid, get the population size per trip
  size <- function(grid) {
    out <- apply(N[, , grid], 2, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
    out <- data.frame(t(out))
    names(out) <- c("lcl", "med", "ucl")
    out$grid <- grid
    out$trip <- 1:5
    return(out)
  }
  
  Nout <- size(1)
  for(i in 2:6) {
    Nout <- rbind(Nout, size(i))
  }
  
  Nout <- full_join(all, Nout)

  min.mice <- min.mice %>%
    arrange(grid, trip)
  
  Nout <- full_join(Nout, min.mice)
  Nout
  
  ggplot(Nout, aes(y = med, x = trip)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2) +
    facet_wrap(~ DataSiteName) +
    ylab("Mark-recapture estimate of population size") +
    xlab("Trip") +
    theme_bw()
  
