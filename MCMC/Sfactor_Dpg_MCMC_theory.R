#
# Sfactor_Dpg_MCMC.R
#
# CONDITIONS:
# - relationship between x and y given by theory [table] 
# - theory relation from Marcucci et al., PRL 116, 102501 (2005)
#   [energy grid of Marcucci 2016 is too coarse] 
#
# PURPOSE: 
# - fit a theoretical model to data
#
# DATA:
# - includes 10 data sets
# - two kinds of data sets: those whose absolute normalization we trust
#   and those for which we trust their relative energy dependence only;
#   ABSOLUTE DATA:
#   sampled using extrinsic scatter and a highly informative lognormal 
#   prior density, where the lognormal parameters mu and sigma 
#   a chosen according to the systematic uncertainty reported for a 
#   given experiment
#   RELATIVE DATA:
#   sampled using extrinsic scatter and a very broad prior to describe
#   the relative normalization
#
# DETAILS:
# - x data have no error; only y data have errors
# - no outlier identification
# - likelihoods have a Gaussian shape
# - true relationship between variables is given by nuclear theory,
#   assuming two parameters:
#   -- a multiplicative factor, a.scale
#   -- an offset, b.offset
# - different data sets are statistically independent
# - includes systematic errors in all sets [y.norm...]
# - includes extrinsic scatter [unreported additional statistical 
#   scatter]
#
# CAREFUL:
# - change the name of the theory file both in the _MCMC and _PLOT
#   scripts
# 
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package 
library(rjags)
library(magicaxis)
## for block updating [we do not need to center predictor variables]
load.module("glm")  
# random number seed
set.seed(2021)
require(MASS)

######################################################################
## FUNCTIONS
######################################################################

######################################################################
# DATA INPUT
# [statistical uncertainties only]
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, S-factor in MeVb

###### DATA SETS WITH INFO ON BOTH STATISTICAL AMND SYSTEMATIC UNCERTAINTIES

## DATA SET 0: tur21
obsx0    <-  c(0.279, 0.400, 0.461, 0.535, 0.670, 0.694,
0.860, 1.033, 1.094, 0.265, 0.333, 0.400, 
0.468, 0.535)
obsy0    <- c(3.1e-6, 4.0e-6, 5.4e-6, 7.4e-6, 11.3e-6, 8.8e-6,
12.9e-6, 14.3e-6, 15.8e-6, 2.9e-6, 4.6e-6, 4.9e-6,
5.7e-6, 5.2e-6)
errobsy0 <- c(0.65e-6, 0.64e-6, 0.43e-6, 0.44e-6, 2.8e-6, 1.8e-6, 
3.7e-6, 2.6e-6, 3.3e-6, 0.38e-6, 0.37e-6, 0.39e-6,
0.46e-6, 0.42e-6)

## DATA SET 1: mos20 **
obsx1    <- c(0.0324,   0.0667,   0.0995,   0.1159,   0.1329,   0.1493, 
              0.1661,   0.1827,   0.1995,   0.2228,   0.2329,   0.2529,   
              0.2629)
obsy1    <- c(0.386e-6, 0.627e-6, 0.850e-6, 0.966e-6, 1.133e-6, 1.223e-6, 
              1.375e-6, 1.475e-6, 1.648e-6, 1.791e-6, 1.866e-6, 2.073e-6, 
              2.156e-6)
errobsy1 <- c(0.014e-6, 0.009e-6, 0.008e-6, 0.009e-6, 0.004e-6, 0.006e-6, 
              0.004e-6, 0.006e-6, 0.003e-6, 0.006e-6, 0.012e-6, 0.012e-6, 
              0.020e-6)
              
## DATA SET 2: tis19 ***
obsx2    <- c(0.097,          0.119,         0.170,         0.210)
obsy2    <- c(0.78111693e-6,  1.0104971e-6,  1.565835e-6,   1.8542257e-6)
errobsy2 <- c(0.13045707e-6,  0.16803717e-6, 0.25806931e-6, 0.26704461e-6)

## DATA SET 3: cas02
obsx3    <- c(2.600E-03, 3.200E-03, 3.200E-03, 3.800E-03, 4.200E-03,
              4.600E-03, 4.600E-03, 5.200E-03, 5.800E-03, 6.100E-03, 6.500E-03,
              6.700E-03, 7.200E-03, 7.800E-03, 8.400E-03, 8.500E-03, 9.100E-03,
              9.100E-03, 9.700E-03, 9.800E-03, 1.040E-02, 1.040E-02, 1.050E-02,
              1.060E-02, 1.120E-02, 1.120E-02, 1.190E-02, 1.190E-02, 1.230E-02,
              1.240E-02, 1.320E-02, 1.320E-02, 1.380E-02, 1.440E-02, 1.460E-02,
              1.500E-02, 1.580E-02, 1.590E-02, 1.630E-02, 1.720E-02, 1.720E-02,
              1.760E-02, 1.850E-02, 1.860E-02, 1.910E-02, 1.970E-02, 1.990E-02,
              1.980E-02, 2.040E-02, 2.120E-02)
obsy3    <- c(1.590E-07, 2.370E-07, 2.560E-07, 2.420E-07, 2.330E-07,
              2.290E-07, 2.510E-07, 2.490E-07, 2.430E-07, 2.520E-07, 2.620E-07,
              2.610E-07, 2.690E-07, 2.730E-07, 2.730E-07, 2.640E-07, 2.620E-07,
              2.710E-07, 2.700E-07, 2.790E-07, 3.010E-07, 2.850E-07, 2.900E-07,
              2.860E-07, 2.880E-07, 2.800E-07, 2.900E-07, 2.890E-07, 2.770E-07,
              2.890E-07, 3.000E-07, 2.680E-07, 3.110E-07, 2.900E-07, 3.460E-07,
              2.990E-07, 2.820E-07, 3.250E-07, 2.960E-07, 3.140E-07, 3.390E-07,
              3.020E-07, 3.240E-07, 3.550E-07, 3.250E-07, 3.280E-07, 3.280E-07,
              3.320E-07, 3.090E-07, 3.280E-07)
errobsy3 <- c(4.900E-08, 4.000E-08, 3.800E-08, 2.900E-08, 1.300E-08,
              1.500E-08, 1.500E-08, 1.800E-08, 1.500E-08, 1.300E-08, 1.100E-08,
              1.300E-08, 1.100E-08, 9.700E-09, 9.170E-08, 1.700E-08, 4.200E-08,
              1.900E-08, 1.800E-08, 1.600E-08, 6.200E-08, 1.900E-08, 5.800E-08,
              1.500E-08, 3.000E-08, 1.400E-08, 3.100E-08, 2.800E-08, 2.400E-08,
              2.500E-08, 2.500E-08, 2.700E-08, 2.200E-08, 1.900E-08, 2.400E-08,
              2.400E-08, 3.400E-08, 1.600E-08, 2.000E-08, 2.500E-08, 1.400E-08,
              1.800E-08, 2.200E-08, 1.400E-08, 1.900E-08, 1.800E-08, 2.200E-08,
              3.100E-08, 1.700E-08, 1.200E-08)

## DATA SET 4: sch97 ***
obsx4    <- c(1.000E-02, 1.670E-02, 2.330E-02, 3.000E-02, 3.670E-02, 4.330E-02,
              5.000E-02)
obsy4    <- c(2.425E-07, 2.740E-07, 3.452E-07, 3.974E-07, 4.452E-07, 4.738E-07,
              4.744E-07)
errobsy4 <- c(1.250E-08, 7.500E-09, 6.500E-09, 6.100E-09, 5.700E-09, 7.200E-09,
              6.400E-09)

## DATA SET 5: ma97 **
obsx5    <- c(7.490E-02, 1.070E-01, 1.330E-01, 1.730E-01)
obsy5    <- c(6.850E-07, 7.080E-07, 9.560E-07, 1.260E-06)
errobsy5 <- c(7.020E-08, 6.840E-08, 8.400E-08, 9.820E-08)


## DATA SET 6: war63 ***
obsx6    <- c(0.646,	0.64605,		  1.476,		1.586)
obsy6    <- c(6.0E-6,	5.6E-6,		15.0E-6,		15.6E-6)
errobsy6 <- c(5.4E-7,	3.4E-7,		 7.5E-7,		7.8E-7)
  
             
###### DATA SETS WITH INFO ON RELATIVE MEAN VALUES ONLY

## DATA SET 7: bai70 ***
obsx7    <- c(0.0351,		0.0563,		0.0850,		0.2542,
		    0.3058,		0.3359,		0.4047,		0.4262,
		    0.4578,		0.5237,		0.6571,		0.7245)
obsy7    <- c(0.5668E-6,		0.7140E-6,	0.8902E-6,	1.9624E-6,
		    2.3359E-6,		2.7242E-6,	3.2760E-6,	3.4517E-6,
		    3.7555E-6,		4.3403E-6,	5.9367E-6,	6.4757E-6)

## DATA SET 8: gri63 **
obsx8    <- c(1.500E-02,		1.600E-02,	1.800E-02,	2.000E-02,	
		    2.200E-02,		2.300E-02,	2.400E-02,	2.600E-02,
		    2.700E-02,		2.800E-02,	3.100E-02,	3.200E-02)
obsy8    <- c(4.300E-07,		4.200E-07,	3.900E-07,	3.800E-07,
		    3.900E-07,		4.100E-07,	4.100E-07,	4.700E-07,
		    4.400E-07,		4.400E-07,	4.700E-07,	4.100E-07)

## DATA SET 9: gri62 **
obsx9    <- c(0.183,	0.387,	0.503,	0.657,	1.167)
obsy9    <- c(1.18E-6,	3.13E-6,	4.28E-6,	6.25E-6,	12.16E-6)

## DATA SET 10: gri55 ***
obsx10   <- c(0.1100,		0.2926,		0.4719,		0.6479,	
		    0.8219,		0.8612,		0.9092,		0.9678,	
		    0.9978,		1.0318,		1.0938,		1.1525,		
		    1.1898,		1.2291)
obsy10   <- c(0.7633E-6,		2.5007E-6,	4.5011E-6,	6.1330E-6,
		    8.2896E-6,		8.2425E-6,	9.9648E-6,	9.7644E-6,
		    10.5937E-6,	10.4870E-6,	11.0952E-6,	12.1701E-6,
		    12.9806E-6,	13.0021E-6)

## DATA SET x: bys08
##obsxx    <- c(8.280E-03, 9.490E-03, 10.10E-03)
##obsyx    <- c(2.370E-07, 2.770E-07, 2.980E-07)
##errobsyx <- c(7.100E-08, 6.400E-08, 6.500E-08)

######################################################################
# INPUT OF THEORY MODEL [MARCUCCI et al.]
######################################################################
# E is in MeV, S is in eV b; convert to MeV and MeV b, respectively
theory <- read.table("Marcucci2005.dat", header=FALSE)

interp.x <- theory[,1] 
# convert to MeVb
interp.y <- theory[,2] * 1e-6

# we will use JAGS interp.lin function to use this theoretical S-factor:
#
# - the columns of this table define vectors x and y
# - a single point is given by x_i, y_i
# - interp.lin gives the y value for the x value provided as argument e,
#   interp.lin(e, x, y)

######################################################################                 
######################################################################                 
# rjags ----->

cat('model {

################################
# LIKELIHOODS
################################
# - careful: dnorm is differently defined in R and JAGS! 
# - precision=sigma^(-2)
# - in a for loop, make sure **all** variables on the LEFT of an 
#   expression has the index [i]
# - systematic error as normalization factor y.norm...

### ABSOLUTE DATA
#################

for (i in 1:length(obsx0)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy0[i] ~ dnorm(ya0[i], pow(errobsy0[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya0[i] ~ dnorm(ym0[i], pow(yscat0, -2))
  # ...subject to syst uncertainties: 
  ym0[i] <- y.norm0 * yt0[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt0[i] <- b.offset + a.scale * interp.lin(obsx0[i], interp.x, interp.y)
}    


for (i in 1:length(obsx1)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy1[i] ~ dnorm(ya1[i], pow(errobsy1[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya1[i] ~ dnorm(ym1[i], pow(yscat1, -2))
  # ...subject to syst uncertainties: 
  ym1[i] <- y.norm1 * yt1[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt1[i] <- b.offset + a.scale * interp.lin(obsx1[i], interp.x, interp.y)
}    

for (i in 1:length(obsx2)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy2[i] ~ dnorm(ya2[i], pow(errobsy2[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya2[i] ~ dnorm(ym2[i], pow(yscat2, -2))
  # ...subject to syst uncertainties: 
  ym2[i] <- y.norm2 * yt2[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt2[i] <- b.offset + a.scale * interp.lin(obsx2[i], interp.x, interp.y)
}    

for (i in 1:length(obsx3)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy3[i] ~ dnorm(ya3[i], pow(errobsy3[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya3[i] ~ dnorm(ym3[i], pow(yscat3, -2))
  # ...subject to syst uncertainties: 
  ym3[i] <- y.norm3 * yt3[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt3[i] <- b.offset + a.scale * interp.lin(obsx3[i], interp.x, interp.y)
}    

for (i in 1:length(obsx4)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy4[i] ~ dnorm(ya4[i], pow(errobsy4[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya4[i] ~ dnorm(ym4[i], pow(yscat4, -2))
  # ...subject to syst uncertainties: 
  ym4[i] <- y.norm4 * yt4[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt4[i] <- b.offset + a.scale * interp.lin(obsx4[i], interp.x, interp.y)
}    

for (i in 1:length(obsx5)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy5[i] ~ dnorm(ya5[i], pow(errobsy5[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya5[i] ~ dnorm(ym5[i], pow(yscat5, -2))
  # ...subject to syst uncertainties: 
  ym5[i] <- y.norm5 * yt5[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt5[i] <- b.offset + a.scale * interp.lin(obsx5[i], interp.x, interp.y)
}    

for (i in 1:length(obsx6)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy6[i] ~ dnorm(ya6[i], pow(errobsy6[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya6[i] ~ dnorm(ym6[i], pow(yscat6, -2))
  # ...subject to syst uncertainties: 
  ym6[i] <- y.norm6 * yt6[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt6[i] <- b.offset + a.scale * interp.lin(obsx6[i], interp.x, interp.y)
}    

### RELATIVE DATA
#################

for (i in 1:length(obsx7)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter:
  obsy7[i] ~ dnorm(ya7[i], pow(yscat7, -2))    
  # ...subject to syst shift: 
  ya7[i] <- yt7[i] * 10^y.norm7
  # true sigma [calculated from theory and then scaled]: 
  yt7[i] <- b.offset + a.scale * interp.lin(obsx7[i], interp.x, interp.y)
}    

for (i in 1:length(obsx8)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter:
  obsy8[i] ~ dnorm(ya8[i], pow(yscat8, -2))    
  # ...subject to syst shift: 
  ya8[i] <- yt8[i] * 10^y.norm8
  # true sigma [calculated from theory and then scaled]: 
  yt8[i] <- b.offset + a.scale * interp.lin(obsx8[i], interp.x, interp.y)
}    

for (i in 1:length(obsx9)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter:
  obsy9[i] ~ dnorm(ya9[i], pow(yscat9, -2))    
  # ...subject to syst shift: 
  ya9[i] <- yt9[i] * 10^y.norm9
  # true sigma [calculated from theory and then scaled]: 
  yt9[i] <- b.offset + a.scale * interp.lin(obsx9[i], interp.x, interp.y)
}    

for (i in 1:length(obsx10)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter:
  obsy10[i] ~ dnorm(ya10[i], pow(yscat10, -2))    
  # ...subject to syst shift: 
  ya10[i] <- yt10[i] * 10^y.norm10
  # true sigma [calculated from theory and then scaled]: 
  yt10[i] <- b.offset + a.scale * interp.lin(obsx10[i], interp.x, interp.y)
}    

################################
# PRIORS
################################

### MODEL PRIORS 
#  a.scale ~ dnorm(0.0, pow(10, -2))T(0,)        # "a.scale" cannot become negative!
  a.scale <- 10^a.fac
  a.fac ~ dunif(-1, 1)

  b.offset ~ dnorm(0.0, pow(1e-4, -2))          # additive offset

###################

### DATA SET PRIORS

# ABSOLUTE DATA:

# extrinsic scatter in S-factor
  yscat0 ~ dnorm(0.0, pow(1e-4, -2))T(0,) 
  yscat1 ~ dnorm(0.0, pow(1e-4, -2))T(0,)       # sigma assumed to be 1e-4 MeVb
  yscat2 ~ dnorm(0.0, pow(1e-4, -2))T(0,)
  yscat3 ~ dnorm(0.0, pow(1e-4, -2))T(0,)
  yscat4 ~ dnorm(0.0, pow(1e-4, -2))T(0,)
  yscat5 ~ dnorm(0.0, pow(1e-4, -2))T(0,)
  yscat6 ~ dnorm(0.0, pow(1e-4, -2))T(0,)

# systematic normalization factor for S-factor:
# log(): natural logarithm

  y.norm0 ~ dlnorm(logmu0, pow(logsigma0, -2))
  logmu0 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma0 <- log(1.14)  
  
  y.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
  logmu1 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma1 <- log(1.027) # factor uncertainty is 1.027, i.e. 2.7% for MOS20

  y.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
  logmu2 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma2 <- log(1.10)  # factor uncertainty is 1.10, i.e., 10% for TIS19

  y.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
  logmu3 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma3 <- log(1.045) # factor uncertainty is 1.045, i.e., 4.5% for CAS02

  y.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
  logmu4 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma4 <- log(1.09)  # factor uncertainty is 1.09, i.e., 9% for SCH97 

  y.norm5 ~ dlnorm(logmu5, pow(logsigma5, -2))
  logmu5 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma5 <- log(1.09)  # factor uncertainty is 1.09, i.e., 9% for MA97 

  y.norm6 ~ dlnorm(logmu6, pow(logsigma6, -2))
  logmu6 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma6 <- log(1.10)  # factor uncertainty is 1.10, i.e., 10% for WAR63 

# RELATIVE DATA:

# extrinsic scatter in S-factor
  yscat7 ~ dnorm(0.0, pow(1e-4, -2))T(0,)
  yscat8 ~ dnorm(0.0, pow(1e-4, -2))T(0,)
  yscat9 ~ dnorm(0.0, pow(1e-4, -2))T(0,)
  yscat10 ~ dnorm(0.0, pow(1e-4, -2))T(0,)

# chose very broad [weakly informative] for normalization...
# sampled values are logs of multiplication factors
  y.norm7 ~ dunif(-1, 1)
  y.norm8 ~ dunif(-1, 1)
  y.norm9 ~ dunif(-1, 1)
  y.norm10 ~ dunif(-1, 1)

}', file={f <- tempfile()})

# dunif(min=0, max=1): gives uniform density   
# pow(a,b) = a^b
######################################################################
# n.adapt:  number of iterations in the chain for adaptation  
#           [JAGS will use to choose the sampler and to assure optimum 
#           mixing of the MCMC chain; will be discarded] 
# update(): performs the burn-in on each chain by running the MCMC for 
#           n.burn iterations without saving any of the posterior samples
# coda.samples(): runs each MCMC chain for the number of iterations 
#           specified by n.iter, but it does not save every iteration; 
#           instead, it saves only ever nth iteration, where n is given 
#           by thin
# n.chains: number of mcmc chains

n.chains <- 3
n.adapt  <- 20000 
n.burn   <- 300000
n.iter   <- 1000000
thin     <- 10

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model
#
## jags wants all data in a list
ourmodel <- jags.model(f, data = list('obsx0'=obsx0,'obsy0'=obsy0,'errobsy0'=errobsy0,  
                 'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1,
                 'obsx2' = obsx2, 'obsy2' = obsy2, 'errobsy2' = errobsy2,
                 'obsx3' = obsx3, 'obsy3' = obsy3, 'errobsy3' = errobsy3,
                 'obsx4' = obsx4, 'obsy4' = obsy4, 'errobsy4' = errobsy4,
                 'obsx5' = obsx5, 'obsy5' = obsy5, 'errobsy5' = errobsy5,
                 'obsx6' = obsx6, 'obsy6' = obsy6, 'errobsy6' = errobsy6,
                 'obsx7' = obsx7, 'obsy7' = obsy7,
                 'obsx8' = obsx8, 'obsy8' = obsy8,
                 'obsx9' = obsx9, 'obsy9' = obsy9,
                 'obsx10' = obsx10, 'obsy10' = obsy10,                 
                 'interp.x' = interp.x, 'interp.y' = interp.y
                                     ),
#                inits = list(a.scale = 1.0), 
                 n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn) 
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                 variable.names=c(
		'a.scale', 
#		'a.fac',
		'b.offset','y.norm0',
		'y.norm1', 'y.norm2', 'y.norm3', 'y.norm4', 'y.norm5', 
		'y.norm6', 'y.norm7', 'y.norm8', 'y.norm9', 'y.norm10',
		'yscat0',
		'yscat1',  'yscat2',  'yscat3',  'yscat4',  'yscat5',  
		'yscat6',  'yscat7',  'yscat8',  'yscat9',  'yscat10'		
                                 ), 
                 n.iter=n.iter, thin=thin)

# <---- rjags
######################################################################
######################################################################
# OUTPUT RESULTS TO SCREEN
######################################################################
cat("", "\n")    # output empty line

# sample size adjusted for autocorrelation
effectiveChainLength = effectiveSize(mcmcChain) 
show(effectiveChainLength)

cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 

######################################################################
# TRACES AND DENSITIES
######################################################################
pdf("MCMC_Dpg_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# OUTPUT RESULTS TO FILES
######################################################################
# !!! make sure to check the order of the parameters in the MCMC output !!!

# matrix samplesmat contains the samples from all n.chains
# as.matrix() strips MCMC attributes from mcmc object 
samplesmat = as.matrix(mcmcChain)
write.matrix(samplesmat,"Dpg_SAMP")




