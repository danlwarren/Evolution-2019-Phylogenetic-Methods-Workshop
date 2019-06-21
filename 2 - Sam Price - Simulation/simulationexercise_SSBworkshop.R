# Simulation 

# There are three interelated questions. The first is what is the uncertainty associated with the estimated parameters for the models? The second is given the tree and the data do we have the power to distinguish between the different models we are fitting including the null? The third is do the models adequately represent the data?

# We simulate under our models with the estimated parameter values and re-estimate the model-fit and parameters. This can tell us if we have the power to distinguish between the models we are fitting including the null model, which is key as we want to know we can correctly reject the null model. If we do have the power we should be able to simulate data under the estimated parameters for each of the models and to be able to distinguish between them. If we don't then we make inaccurate biological statements based on estimates from analyses that have inadequate power. It can also tell us the uncertainty associated with the parameter estimates. For model adequacy we can identify important test statistics and fit them to our data, estimate model parameters and simulate under them 


## Two examples of estimating parameter uncertainty##

## 1) Using the integrated OUwie bootstrap function to calculate 95% Confidence Intervals around the parameter estimates from the best-fitting model 
# This will provide an estimate of the degree of uncertainty in our parameter estimates including the influence that the likelihood surface has on the parameter estimates and is a much more rigorous way to present the results compared to the single value estimated from each phylogeny. We are testing the hypothesis that diet influences body mass evolution across marsupials, the expectation is that herbivores (1) will have larger optima than omnivores (2) or carnivores (3).

library(OUwie)
dietmtr<-read.nexus("marsupialdiettree.nex")
marsupiald<-read.table("marsupialmass_diet.txt", header=T)

mods<-c("BM1", "OU1", "OUM", "OUMA", "OUMV") # pick all the models you want to fit - we are using all the OU models that allow different for diet as that matches our prediction and also the two simpler models

sizedietresults<-list()# setting up an empty list for the results
for(i in 1:length(mods)){#loop to run OUwie across all the models you want to fit
	sizedietresults[[i]]<-OUwie(dietmtr, marsupiald, model=mods[i], simmap.tree=FALSE, root.station=TRUE, diagn=T) #note you probably don't actually want to run this with root.station=TRUE for the BMS model
}
dietbestfit<-c(sizedietresults[[1]]$AICc, sizedietresults[[2]]$AICc, sizedietresults[[3]]$AICc, sizedietresults[[4]]$AICc, sizedietresults[[5]]$AICc)
names(dietbestfit)<-mods# 
dietbestfit-min(dietbestfit)

# We want to use the parameters from the best-fitting model
sizedietresults[[4]]

# This will take a long time to run but you can read in the results from the analyses of 1000 I ran earlier (takes about 8 hrs on my mac 3.1 GHz Intel Core i7) DON'T RUN DURING WORKSHOP

sizedietoumboot1000<-OUwie.boot(dietmtr, marsupiald, model=c("OUMA"), nboot=1000, alpha=sizedietresults[[4]]$solution[1,] , sigma.sq=sizedietresults[[4]]$solution[2,], theta=sizedietresults[[4]]$theta[,1], theta0= sizedietresults[[4]]$theta[2,1]) #When theta0 was dropped from the model (theta0=T) i.e. not estimated, then we need to set the root as the the value of the continuous trait for the selective regime mapped at the root which is omnivore (2)

write.table(sizedietoumboot1000, file="sizedietoumboot1000.txt")

sizedietoumboot1000 <-read.table("sizedietoumboot1000.txt")

head(sizedietoumboot1000)

# Plot the entire distribution of the alpha parameter (sigma was fixed between the different diets as the best fitting model was OUMA)
plot(density(sizedietoumboot1000[,1]), col="green", ylim=c(0,30), xlab="Marsupial Size Alpha Estimate") # distribution of the alpha estimates for herbivores
lines(density(sizedietoumboot1000[,2]), col="purple")# distribution of the alpha estimates for omnivores
lines(density(sizedietoumboot1000[,3]), col="blue")# distribution of the alpha estimates for carnivores
abline(v= sizedietresults[[4]]$solution[1,1], col="green")# This is the empirical estimate of alpha
abline(v= sizedietresults[[4]]$solution[1,2], col="purple")# This is the empirical estimate of alpha
abline(v= sizedietresults[[4]]$solution[1,3], col="blue")# This is the empirical estimate of alpha

# We can also calculate the 95% CI of the parameter estimates easily.
C195sizedietoumboot1000 <-apply(sizedietoumboot1000, 2, quantile, probs=c(0.025,0.975)) 

# What conclusions can you draw from this? Is it different for the theta parameter?

# Why might you not use bootstrapping?  If your dataset is very large then the time taken to run the simulations may be prohibitive, especially if you are running it on 1000 stochastically mapped trees.


## 2)Estimate the 95% CI interval around the estimate of Pagel's lambda for parrotfish. Pagel's lambda is frequently used to estimate the degree of phylogenetic signal within a continuous trait or of the residuals from a phylogenetic regression by estimating the optimal transformation of the internal branches. Thus lambda=0 you have a star phylogeny and no signal but when lambda=1 is the identical phylogeny and indicative of Brownian motion. Pagel's lambda is estimated in geiger's fitContinuous function model="lambda"


library(geiger)

ptr<-read.nexus("parrotfishtree.nex")
parrotfishmorph<-read.table("parrotfishesRawMorphology.txt", header=T, row.names=1)
rownames(parrotfishmorph)[1]<-"Bolbometopon_muricatum" #this is a misspelling in the dataset!

foo<-treedata(ptr, parrotfishmorph, sort=T) # match dataset and tree

prunedmorph<-foo[[2]]
prunedtr<-foo[[1]]

loglength<-log(prunedmorph[,1])
names(loglength)<-rownames(prunedmorph)

lambdalength<-fitContinuous(prunedtr, loglength, model="lambda") #Estimate the ML estimate of lambda, root state z0 and sigma^2

lambdatree<-rescale(prunedtr, model="lambda",  lambdalength$opt$lambda) # to simulate under the lambda model we transform the branch lengths of the tree by the estimated lambda parameter and then run a simple Brownian motion simulation on the transformed tree

par(mfrow=c(1,2))
plot(prunedtr)
plot(lambdatree)

lambdasim<-sim.char(lambdatree, lambdalength$opt$sigsq, nsim=500, model="BM", root= lambdalength$opt$z0)
# this takes a few minutes to run so I have provided the output
lambdasimres<-vector()
for(i in 1:dim(lambdasim)[3]){ # a loop that fits the lambda model to the 500 simulated characters using fitContinuous and puts the estimated lambda value in a vector called lambdasimres
	print(i)
	tmp<-fitContinuous(prunedtr , lambdasim[,,i], model="lambda")
	lambdasimres<-c(lambdasimres, tmp$opt$lambda)
}

write.csv(lambdasimres, file="ParrotfishSLlambdasim.csv")
lambdasimres <-read.csv("ParrotfishSLlambdasim.csv", row.names=1) # read in existing file to save time

plot(density(lambdasimres[,1]))# you only need the column if you are reading in the data from the csv, if not lambdasimres will suffice
abline(v= lambdalength$opt$lambda)
parrotfishlengthlambdaCI<-quantile(lambdasimres[,1], c(0.025, 0.975))# you only need the column if you are reading in the data from the csv, if not lambdasimres will suffice

## What conclusions do you draw about the uncertainty when estimating lambda for this dataset? What does this tell us about the phylogenetic signal within the data?



## One example of estimating power and the ability to distinguish between evolutionary models ##

# Constructing a parametric bootrap test to look at power and uncertainty using discrete trait model-fitting and Parrotfishes 

# reading in data and fitting it to tree
ptr<-read.nexus("parrotfishtree.nex")

pecol<-read.table("parrotfishecology.txt", header=T, row.names=1)
row.names(pecol)[1]<-"Bolbometopon_muricatum"
foo<-treedata(ptr, pecol, sort=T)

prunedecol<-foo[[2]]
prunedtr<-foo[[1]]

# generating dataset that distinguishes scapers from other types of feeder.
feeding<-rep(1, length(prunedecol[,3]))
feeding[prunedecol[,3]=="scrape"]<-2
names(feeding)<-row.names(prunedecol)

# testing to see whether we should reconstruct the history with equal or unequal rates between scraping and other modes of feeding - which one fits best?
equalfit<-fitDiscrete(prunedtr, feeding, model="ER")
diffit<-fitDiscrete(prunedtr, feeding, model="ARD")

# Set up Q matrices for the two models we fitted 

qER <- list(rbind(c(-0.007394146, 0.007394146), c(0.007394146, -0.007394146))) # this is the q matrix estimated from the data using fit discrete ER
qARD<- list(rbind(c(-0.008830519, 0.008830519), c(0.002747460, -0.002747460))) # this is the q matrix estimated from the data using fit discrete ARD

# Now simulate 1000 characters using those models 
simER <- sim.char(prunedtr, qER, model="discrete", n=100) # Note just run 100 due to time constraints but in reality you would run 1000+ 
simARD <- sim.char(prunedtr, qARD, model="discrete", n=100)
write.table(simER, file="SimulatedDiscreteTraitER.txt") # Always good practice to save the simulated data in case you need to go back and check something
write.table(simARD, file="SimulatedDiscreteTraitARD.txt")


# Now set up empty matrices for the results - think through all the information you want to save. If you are not sure I recommend you save it all as an rdata object and process the results after running the model-fitting. 

ERsimresultstable<-matrix(nrow=1000, ncol=6)
colnames(ERsimresultstable)<-c("ERsimdeltaAICc","ERestlnLER","ERestRateER", "ERestlnLARD","ERestRateARDq12", "ERestRateARDq21")
ARDsimresultstable<-matrix(nrow=1000, ncol=6)
colnames(ARDsimresultstable)<-c("ARDsimdeltaAICc","ARDestlnLER","ARDestRateER", "ARDestRatelnLARD", "ARDestRateARDq12", "ARDestRateARDq21")

# Now write loop to run fitDiscrete with model="ER" and "ARD" on the datasets simulated under the ER and ARD model THIS WILL TAKE TOO LONG TO RUN FOR THE EXERCISE READ IN THE OUTPUT BELOW
for(i in 1:dim(simER)[3]) { # again this should be done on all your simulated datasets 1000+ but that will take too long 
	print(i) #useful to see how far through you are
	tmp1<-vector();tmp3<-vector();tmp3<-vector();tmp4<-vector(); # set vectors to empty each time
	try(tmp1<-fitDiscrete(prunedtr, simER[,,i], model="ER")) # unfortunately due to the possibility that these models can generate datasets with all the species in the same category, which makes it impossible to fit these models leading it to error out of the loop we have to use a 'try catch' which means an exception can be caught and then the rest of the script will continue to run. You could also remove these instances prior to running the loop or write an if statement to only run the analyses when there are at least 2 states.
	try(tmp1<-fitDiscrete(prunedtr, simER[,,i], model="ER"))
	try(tmp2<-fitDiscrete(prunedtr, simER[,,i], model="ARD"))
	try(tmp3<-fitDiscrete(prunedtr, simARD[,,i], model="ER"))
	try(tmp4<-fitDiscrete(prunedtr, simARD[,,i], model="ARD"))
	try(ERsimresultstable[i,1]<-tmp1$opt$aicc-tmp2$opt$aicc )# Delta AICc. As we simulated under the equal rates model we would expect that the ER model would be the best-fitting model (lower AICc score) and thus this value should be negative and ideally >2 
	try(ERsimresultstable[i,2:6]<-c(tmp1$opt$lnL, tmp1$opt$q12, tmp2$opt$lnL, tmp2$opt$q12, tmp2$opt$q21)) # concatenating all the parameters we are interested in and putting them in the results table
	try(ARDsimresultstable[i,1]<-tmp4$opt$aicc-tmp3$opt$aicc )# Delta AICc. As we simulated under the different rates model we would expect that the ARD model would be the best-fitting model (lower AICc score) and thus this value should be negative and ideally >2 
	try(ARDsimresultstable[i,2:6]<-c(tmp3$opt$lnL, tmp3$opt$q12, tmp4$opt$lnL, tmp4$opt$q12, tmp4$opt$q21)) 
}

write.table(ERsimresultstable, file="DiscreteModelfitting_ERsimresults.txt")
write.table(ARDsimresultstable, file="DiscreteModelfitting_ARDsimresults.txt")

quit()

# Because we didn't have time to run this during the workshop I have provided the results, read them in now!

ERsimresultstable<-read.table("DiscreteModelfitting_ERsimresults.txt", na.strings="NA")

ARDsimresultstable<-read.table("DiscreteModelfitting_ARDsimresults.txt", na.strings="NA")

# So first of all look at the delta AIC, if we can distinguish the ER from the ARD models? (You could also consider doing this for the Likelihood Ratio Test as AIC is less reliable, see Boettiger et al. 2012 for an example). This tells you whether you have the power to estimate these models reliably.

par(mfrow=c(1,2))
hist(ERsimresultstable$ERsimdeltaAICc, main="Simulated ER data", xlab="Delta AICc ER - ARD")
abline(v=median(ERsimresultstable$ERsimdeltaAICc, na.rm=T), col="blue")
abline(v=-2, col="red", lty=2)
hist(ARDsimresultstable$ARDsimdeltaAICc, main="Simulated ARD data", xlab="Delta AICc ARD - ER")
abline(v=median(ARDsimresultstable$ARDsimdeltaAICc, na.rm=T), col="blue")
abline(v=-2, col="red", lty=2)


# Now we can also calculate the 95% Confidence Interval around the estimated parameter values and compare these to the values we estimated from the original data and that we simulated under. This tells you how much uncertainty there is in the parameter estimates, these can be reported alongside your empirical estimates. 

ER95CI<-apply(ERsimresultstable, 2, quantile, probs=c(0.025, 0.975), na.rm=T) # applying the quantile function to every column in the results table (margin =2) is a quick way to calculate the 95% CI
ARD95CI<-apply(ARDsimresultstable, 2, quantile, probs=c(0.025, 0.975), na.rm=T)

ER95CI[,3]
qER
ARD95CI[,5:6]
qARD

# There are a variety of factors which can determine whether you have the power to detect different models of evolution and the degree of uncertainty in the parameter estimates, this includes tree and dataset size, the magnitude of the difference and the distribution of the characters upon the tree.

## One example of estimating model adequacy using summary statistics ##

# Use arbutus to run a simple test of BM vs Early burst vs OU on the parrotfish standard length data (loglength created earlier in the script along with prunedtr). Arbutus uses 6 summary statistics based on independent contrasts, for convenience I have pasted the description of each stasticic from Pennell et al. 2015 below:

#MSIG The mean of the squared contrasts, is a metric of overall rate. Violations detected by MSIG indicate whether the overall rate of trait evolution is over- or underestimated.

#CVAR The coefficient of variation (standard deviation/mean) of the absolute value of the contrasts. If CVAR calculated from the observed contrasts is greater than that calculated from the simulated contrasts, it suggests that we are not properly accounting for rate heterogeneity across the phylogeny. If CVAR from the observed is smaller, it suggests that contrasts are even more than the model assumes. The coefficient of variation is used rather than the variance because the mean and variance of contrasts can be highly correlated.

#SVAR The slope of a linear model fitted to the absolute value of the contrasts against their expected variances. Each (standardized) contrast has an expected variance proportional to the sum of the branch lengths connecting the node at which it is computed to its daughter lineages. Under a model of Brownian motion, we expect no relationship between the contrasts and their variances. It is used to test whether contrasts are larger or smaller than we expect based on their branch lengths. If, for example, more evolution occurred per unit time on short branches than long branches, we would observe a negative slope. If SVAR calculated from the observed data deviates substantially from the expectations, a likely explanation is branch length error in the phylogenetic tree.

#SASR The slope of a linear model fitted to the absolute value of the contrasts against the ancestral state inferred at the corresponding node. The ancestral state is estimated using the least squares method suggested by Felsenstein (1985) for the calculation of contrasts - this is more properly thought of as a weighted average value for each node. This statistic evaluates whether there is variation in rates relative to the trait value. For example, do larger organisms evolve proportionally faster than smaller ones?

#SHGT The slope of a linear model fitted to the absolute value of the contrasts against node depth also known as the node-height test. This is used to capture variation relative to time and has been used to detect early bursts of trait evolution during adaptive radiations.

#DCDF The D statistic obtained from a Kolmolgorov-Smirnov test from comparing the distribution of contrasts to that of a normal distribution with mean 0 and standard deviation equal to the root of the mean of squared contrasts (the expected distribution of the contrasts under Brownian motion. This to captures deviations from normality. For example, if traits evolved via a “jump-diffusion”-type process (Landis et al. 2013) in which there were occasional bursts of rapid phenotypic evolution, the tip data would no longer be multivariate normal owing to a few contrasts throughout the tree being much larger than the rest (i.e., the distribution of contrasts would have heavy tails).

# Make sure you have the develtools installed then install arbutus from github

devtools::install_github("mwpennell/arbutus")
library(arbutus)

# Now fit a Brownian motion, OU and early burst model using fitContinuous in geiger (Slater & Pennell used the MCMC version) to the parrotfish length data (loglength created earlier in the script along with prunedtr)

bmlength<-fitContinuous(prunedtr, loglength, model="BM")
oulength<-fitContinuous(prunedtr, loglength, model="OU")
eblength<-fitContinuous(prunedtr, loglength, model="EB")


## Which model fits best?

# Check adequacy of the three models

bmlength_modeladeq<-arbutus(bmlength, nsim=1000) 
plot(bmlength_modeladeq)
eblength_modeladeq<-arbutus(eblength, nsim=1000) 
plot(eblength_modeladeq)
oulength_modeladeq<-arbutus(oulength, nsim=1000) 
plot(oulength_modeladeq)

# Is the best fitting model adequate? Are any of the models adequate?

## Can other models look like OU? yes, several things 1) a late burst (this is why fitContinuous prevents it from fitting a late burst unless you alter the bounds of the rate of decay) and 2) white noise which is Brownian motion without a phylogenetic covariance structure. We can try fitting these two models. If we try fitting it the decay parameter is at the + bounds indicating a very strong rate increase but what about model adequacy?

eblblength<-fitContinuous(prunedtr, loglength, model="EB", bounds=list(a=c(-5,20)))# this is at the bounds of 5
eblblength_modeladeq<-arbutus(eblblength, nsim=1000) 

wnlength<-fitContinuous(prunedtr, loglength, model="white")
wnlength_modeladeq<-arbutus(wnlength, nsim=1000) 

# Generate independent contrasts and visualize them on the phylogeny to see if there is anything that looks unusual

lengthpic<-pic(loglength, prunedtr, scaled=T)
hist(lengthpic)
plot(prunedtr)
nodelabels(round(lengthpic,2), cex=0.5)

# Try removing the very short node with the large contrast - does it make a difference? 

prunedtr1<-drop.tip(prunedtr, "Calotomus_spinidens")

bmlength1spdrop<-fitContinuous(prunedtr1, loglength, model="BM")
oulength1spdrop<-fitContinuous(prunedtr1, loglength, model="OU")
eblblength1spdrop<-fitContinuous(prunedtr1, loglength, model="EB", bounds=list(a=c(-5,20)))# this is at the bounds of 5
wnlength1spdrop<-fitContinuous(prunedtr1, loglength, model="white")

bmlength1spdrop_modeladeq<-arbutus(bmlength1spdrop, nsim=1000) 
plot(bmlength1spdrop_modeladeq)
eblength1spdrop_modeladeq<-arbutus(eblblength1spdrop, nsim=1000) 
plot(eblength1spdrop_modeladeq)
oulength1spdrop_modeladeq<-arbutus(oulength1spdrop, nsim=1000) 
plot(oulength1spdrop_modeladeq)
wnlength1spdrop_modeladeq<-arbutus(wnlength1spdrop, nsim=1000) 
plot(wnlength1spdrop_modeladeq)

## What do you conclude, are any of the models adequate? Do you think you have the power to distinguish between the models?  