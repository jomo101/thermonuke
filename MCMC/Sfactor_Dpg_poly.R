###############    ##################################################     #############
#Author: Joseph Moscoso (DECEMBER 2020)
#######################################################################################
#
# Monte Carlo Markov Chain (MCMC) for S-factor for Deuterium-proton reaction
# D(p,g)3He
#
# Conditions:
# - relationship between x and y values given by polynomial fit with cubic fit parameters
# - In conjuction to Iliadis Marcucci theory fit, including new data values from Gran Sasso LUNA
# - Incorporate Data sets to compare S-factor fit with leading fit of 2020 (Mossa)
#
# Purpose:
# - learning experience for analyzing S-factor data using MCMC chain and comparing fit parameters
# - fit polynomial curve to DATA
#
# Details:
# - y data have provided errors for 6 data sets 
# - Data sets have published papers: ma97, sch97, cas02, bis08, tis19, mos20
# - Using RJAGS to perform the MCMC chain: https://cran.r-project.org/web/packages/rjags/index.html
# - different data sets are statistically independent (required)
# - A Cubic fit is assumed with fit parameters: alpha(x^0), beta(x^1), delta(x^2), gamma(x^3) 
# -* each data set fit also includes systematic errors (y.scale#)
# -* each data set fit includes extrinsic scatter 
# -* each data set fit includes stat uncertainties 
#
#######################################################################################
# REMOVA ALL VARIABLES BEFORE WE BEGIN
rm(list=ls())
# import packages we need
library(rjags)
load.module("glm")
#Tip. Use the set.seed function when running simulations to ensure all results, figures, etc are #reproducible.
set.seed(523)
require(MASS)
#######################################################################################
#

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
#######################################################################################
#######################################################################################    
#            RJAGS MODEL            #   


cat('model {

# ABSOLUTE DATA

	for (i in 1:length(x0)) {
		y0[i] ~ dnorm(g0[i],pow(s0[i],-2)) 
		g0[i] ~ dnorm(g0p[i], pow(yscat0, -2))
		g0p[i] = y.norm0 * z0[i]
		z0[i] = alpha + beta * x0[i] + delta * x0[i]^2 + gamma * x0[i]^3
}

	for (i in 1:length(x1)) {
		# S-FACTOR
		# take into account statistical uncertainties
		y1[i] ~ dnorm(g1[i],pow(s1[i],-2))
		# extrinsic scatter
		g1[i] ~ dnorm(g1p[i], pow(yscat1, -2))
		# systematic uncertainties
		g1p[i] = y.norm1 * z1[i]
		# Model true S-factor from the data
		z1[i] = alpha + beta * x1[i] + delta * x1[i]^2 +  gamma * x1[i]^3
}

	for (i in 1:length(x2)) {
		y2[i] ~ dnorm(g2[i],pow(s2[i],-2)) 
		g2[i] ~ dnorm(g2p[i], pow(yscat2, -2))
		g2p[i] = y.norm2 * z2[i]
		z2[i] = alpha + beta * x2[i] + delta * x2[i]^2 + gamma * x2[i]^3
}

	for (i in 1:length(x3)) {
		y3[i] ~ dnorm(g3[i],pow(s3[i],-2)) 
		g3[i] ~ dnorm(g3p[i], pow(yscat3, -2))
		g3p[i] = y.norm3 * z3[i]
		z3[i] = alpha + beta * x3[i] + delta * x3[i]^2 + gamma * x3[i]^3
}

	for (i in 1:length(x4)) {
		y4[i] ~ dnorm(g4[i],pow(s4[i],-2)) 
		g4[i] ~ dnorm(g4p[i], pow(yscat4, -2))
		g4p[i] = y.norm4 * z4[i]
		z4[i] = alpha + beta * x4[i] + delta * x4[i]^2 + gamma * x4[i]^3
}

	for (i in 1:length(x5)) {
		y5[i] ~ dnorm(g5[i],pow(s5[i],-2)) 
		g5[i] ~ dnorm(g5p[i], pow(yscat5, -2))
		g5p[i] = y.norm5 * z5[i]
		z5[i] = alpha + beta * x5[i] + delta * x5[i]^2 + gamma * x5[i]^3
}

	for (i in 1:length(x6)) {
		y6[i] ~ dnorm(g6[i],pow(s6[i],-2)) 
		g6[i] ~ dnorm(g6p[i], pow(yscat6, -2))
		g6p[i] = y.norm6 * z6[i]
		z6[i] = alpha + beta * x6[i] + delta * x6[i]^2 + gamma * x6[i]^3
}

## RELATIVE DATA 
	for (i in 1:length(x7)) {
		# extrinsic scatter
		y7[i] ~ dnorm(g7[i],pow(yscat7,-2)) 
		# systematic shift
		g7[i] = z7[i] * 10^y.norm7
		# true S-factor
		z7[i] = alpha + beta * x7[i]+ delta * x7[i]^2 + gamma * x7[i]^3
}
	for (i in 1:length(x8)) {
		y8[i] ~ dnorm(g8[i],pow(yscat8,-2)) 
		# systematic shift
		g8[i] = z8[i] * 10^y.norm8
		# true S-factor
		z8[i] = alpha + beta * x8[i]+ delta * x8[i]^2 + gamma * x8[i]^3
}
	for (i in 1:length(x9)) {
		y9[i] ~ dnorm(g9[i],pow(yscat9,-2)) 
		# systematic shift
		g9[i] = z9[i] * 10^y.norm9
		# true S-factor
		z9[i] = alpha + beta * x9[i]+ delta * x9[i]^2 + gamma * x9[i]^3
}
	for (i in 1:length(x10)) {
		y10[i] ~ dnorm(g10[i],pow(yscat10,-2)) 
		# systematic shift
		g10[i] = z10[i] * 10^y.norm10
		# true S-factor
		z10[i] = alpha + beta * x10[i]+ delta * x10[i]^2 + gamma * x10[i]^3
}


# PRIORS
# - precision=sigma^(-2)
# each parameter is picked from gaussian distributions

alpha ~ dnorm(0.0,pow(10,-2))
beta ~ dnorm(0.0,pow(10,-2))
delta ~ dnorm(0.0,pow(10,-2))
gamma ~ dnorm(0.0,pow(10,-2))


# ABSOLUTE DATA

# extrinsic scatter in S-factor
yscat0 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat1 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat2 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat3 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat4 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat5 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat6 ~ dnorm(0.0,pow(1e-4,-2))T(0,)

#y.norm is the systematic uncertainty for each data set

# these look completely out of order considering your data input, please
# check again:
#y.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
#logmu1 = log(1.0)
#logsigma1 = log(1.09) # uncertaint of 1.09 for Ma9, 9%

#y.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
#logmu2 = log(1.0)
#logsigma2 = log(1.08) #<1.08 for BYS08 <8%

#y.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
#logmu3 = log(1.0)
#logsigma3 = log(1.09) #9% for SCH97

#y.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
#logmu4 = log(1.0)
#logsigma4 = log(1.045) #4.5% for CAS02

#y.norm5 ~ dlnorm(logmu5, pow(logsigma5, -2))
#logmu5 = log(1.0)
#logsigma5 = log(1.027) # 2.7% for MOS20

#y.norm6 ~ dlnorm(logmu6, pow(logsigma6, -2))
#logmu6 = log(1.0)
#logsigma6 = log(1.1) #10% for TIS19

# here is my input; but please check: the systematic uncertainties for
# each data set are given in the appendix of our paper
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

# RELATIVE DATA

# extrinsic scatter in S-factor
yscat7 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat8 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat9 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
yscat10 ~ dnorm(0.0,pow(1e-4,-2))T(0,)
#
# chose very broad [weakly informative] for normalizaion 
y.norm7 ~ dunif(-1, 1)
y.norm8 ~ dunif(-1, 1)
y.norm9 ~ dunif(-1, 1)
y.norm10 ~ dunif(-1, 1)


}', file={f <- tempfile()})   

# dunif(min=0, max=1): gives uniform density   
# pow(a,b) = a^b
#######################################################################################
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

mimodela <- jags.model(f,
			data= list(
			'x0'=obsx0,'y0'=obsy0,'s0'=errobsy0,
			'x1'=obsx1,'y1'=obsy1,'s1'=errobsy1,
			'x2'=obsx2,'y2'=obsy2,'s2'=errobsy2, 
			'x3'=obsx3,'y3'=obsy3,'s3'=errobsy3,
			'x4'=obsx4,'y4'=obsy4,'s4'=errobsy4,
			'x5'=obsx5,'y5'=obsy5,'s5'=errobsy5,
			'x6'=obsx6,'y6'=obsy6,'s6'=errobsy6,
			'x7'=obsx7,'y7'=obsy7,
			'x8'=obsx8,'y8'=obsy8,
			'x9'=obsx9,'y9'=obsy9,
			'x10'=obsx10,'y10'=obsy10),
			n.chains = n.chains,
			n.adapt = n.adapt)
			
update(mimodela, n.burn)
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(mimodela,
			  variable.names=c('alpha','beta','delta','gamma',
			  'y.norm0', 'y.norm1', 'y.norm2', 'y.norm3', 'y.norm4', 'y.norm5', 
			  'y.norm6', 'y.norm7', 'y.norm8', 'y.norm9', 'y.norm10',
			  'yscat0', 'yscat1',  'yscat2',  'yscat3',  'yscat4',  'yscat5',  
			  'yscat6',  'yscat7',  'yscat8',  'yscat9',  'yscat10'),
			  n.iter=n.iter,
			  thin=thin)

# <---- rjags
#######################################################################################
#######################################################################################
# OUTPUT RESULTS TO SCREEN
#######################################################################################
cat("", "\n")    # output empty line

# sample size adjusted for autocorrelation
effectiveChainLength = effectiveSize(mcmcChain) 
show(effectiveChainLength)			  
 
cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 			  
			
#######################################################################################
# TRACES AND DENSITIES
#######################################################################################
pdf("MCMC_Dpg_Traces.pdf")
plot(mcmcChain)
dev.off()

#######################################################################################
# OUTPUT RESULTS TO FILES
#######################################################################################
samplesmat = as.matrix(mcmcChain)
gee = samplesmat[,1:4]
write.matrix(samplesmat,"Dpg_SAMPj")

