library(unmarked)
library(AICcmodavg)
library(aod)
library(ggplot2)
library(Rcpp)
library(lme4)

#Importing Data into R ----
#First need to read in data from distance surveys (see files for how to set up), then separately read
#in transect covariates and transect lengths (if transect length varies)
#First file is just each observation in a two column format. 
#First column is perpendicular distance to transect
#Second column is the transect it was found on. For empty transects leave a blank in the first column
#See file for example on data structure

#Second file has one transect per row with its potential covariates in each column. Landcover data has already been standardized

#Left in the file paths I used as examples
#Using data from the inland habit of the Kennedy Space Center in this example.
inlandcovs <- read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/siteCov_Inland_Up.csv")
inlanddists<-read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/BurrowDist_Inland_Up.csv")


# Formating your data for Unmarked----
#Make a vector that contains the length of each transect, in the SAME order they are set up in the 'covs' file
#If transect length varies, can read lengths in from a file and use it that way (master length above is example of this)
#If all transects are same length can simply make a vector using the rep() command
#Below example would be for 50 transects, each 100 (distance unit) long
length<-rep(130,200)

# Checking covariates for correlations
#Have a priori standard for what is too high of a correlation
correlations<-cor(inlandcovs[,8:25])
View(correlations)

#Need to format data. 'breaks' is the bin range and spacing for distance models (so bins range from 0 to 30 meters in 0.8m increments)
#yDat and umf are from unmarked, preparing data for models
#Refer to unmarked help using ?unmarkedFrame() for double observer datasets
yDat<-formatDistData(inlanddists,distCol="dist", transectNameCol="transect",dist.breaks=breaks)
umf<-unmarkedFrameDS(y=as.matrix(iDat), siteCovs=inlandcovs,survey="line",
                      dist.breaks=breaks,tlength=length,unitsIn="m")

#Visualizing data, generates a histogram of frequency over distance. 
#If looks off, try adjusting bin size using the breaks variable
hist(umf, xlab="distance (m)", cex.lab=0.8, cex.axis=0.8,main="Burrow detection rates over distance")


#Selecting Detection Distributions ----
#Testing what of the four main potential distributions best fit the data using null models
haz_Null<-distsamp(~1 ~1,umf, keyfun="hazard",output="density",unitsOut="ha")
halfnorm_Null<-distsamp(~1 ~1,umf, keyfun="halfnorm",output="density",unitsOut="ha")
exp_Null<-distsamp(~1 ~1,umf, keyfun="exp",output="density",unitsOut="ha")
uni_Null<-distsamp(~1 ~1,umf, keyfun="uniform",output="density",unitsOut="ha")

#now comparing models to select best fit based on AIC values
modSel(fitList(haz_Null,halfnorm_Null,exp_Null,uni_Null))
#In this example, treating exponential as the best fit, will continue to use throughout models


#Burrow Detection Models ----
#Now need to start on detection covariate testing, formatting follows unmarked user manual.
#See ?distamp() for details on model construction. 
#'+' indicates additive effect, '*' for including interaction term in models
#no a priori reason for any interactions here, and they greatly increase number of parameters
#Not testing for effects of veg on density, since few landcover types found across majority of transects

#Below models are for the master dataset. Change 'umf' to appropriate model to test others.
det_m1<-distsamp(~1 ~1, umf, keyfun="exp",output="density",unitsOut="ha")
det_m2<-distsamp(~obs ~1, umf, keyfun="exp",output="density",unitsOut="ha")

#Now need to compare model set to pick top model. modSel function ranks them in ascending order
modSel(fitList(det_m1,det_m2))
#Change umf to others, and test per-site differences in detection covariates
#Be aware, not all models should be tested for detection


#Burrow Density Modelling ----
#Make sure when setting up models, each makes biological sense
fm1<-distsamp(~loc + obs ~1, umf, keyfun="exp",output="density",unitsOut="ha")
fm2<-distsamp(~obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, umf, keyfun="exp",output="density",unitsOut="ha")

#If models did not converge, need to adjust for that. Cannot compare non-converging models to others
#If models do not converge, or produce NAs for values, can try changing maxit, or giving starting values
#See below model for example with changed maxit
#May also get an an error 'Hessian is singular' this is likely due to one covariate containing a subset of information
#found in another covariate. Below is an example of how to change the maxit

fm3<-distsamp(~obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, umf, keyfun="exp",output="density",unitsOut="ha", control = list(maxit=100))

#For adding starting values, zero is a good value, 1 for potential positive effect, -1 for negative effect
#Need to know total number of parameters when inputing starts
#See below example

fm4<-distsamp(~obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
             , umf, keyfun="exp",output="density",unitsOut="ha", starts=rep(0,26))

modSel(fitList(fm1,fm2,fm3,fm4))

#Model 2 is best model
#now need to test fit
fitfm2 <- Nmix.gof.test(fm2,nsim=25,plot.hist=TRUE)
fitfm2
#c-hat is a measure of model dispersion. Want to be close to 1
#If much greater than 1, need to correct for error
#also simulates predicted values (change nsim to higher values for better precision, low value to run fast) to see if predicted matches observed. 
#Recommend using low value for nsim during initial testing for quick check of fit, 
#and increasing value for precision later
#Want a high p-value, indicating observed and predicted are not significantly different
#If have low p or high c-hat, consider alternative models. May be missing a source of variation
#Or if overall dataset has transects in very different habitat types, may be due to large number of zeros in dataset
#If that is the case, consider splitting up data into sub-models to better partition variance due to habitat differences

#Once we have our top models, need to look at parameter estimates for the covariates
#renaming them for ease of identification
gopherTortoiseBurrow_top<-fm2

#calculate asymtotic profile CI, more accurate than just using SE to estimate, will default to 95% CI
GTB_paraciint <- confint(gopherTortoiseBurrow_top, type='state', method = 'profile')

#This will plot out the detection function
#if you want to look at different detection functions, just change the model name
modelrate<-coef(gopherTortoiseBurrow_top, type='det')
#pulls out the estimates just for detection, below plot is based on the intercept alone. 
plot(function(x) gxexp(x, rate=modelrate[1]), 0, 30,
     xlab="Distance (m)", ylab="Detection prob.")

#if top model had a covariate for detection, need to combine modelrate as below
#will allow you to view alternative detection curves, should always be rate_index (1) + rate_index (n)
#unless contrast matrix was changed by user (outside the scope of this code)
plot(function(x) gxexp(x, rate=(modelrate[1]+modelrate[2])), 0, 30,
     xlab="Distance (m)", ylab="Detection prob.")
#Could use add line instead of plot to overlap detection curves, but not displayed here


#Burrow Occupancy Models ----------
#Here we will model the occupancy rates of burrow from data collected via scoping
#Burrow data has been coded as a binomial response where 1=Occupied, 0=Unoccupied.
#Burrows previously classifed as unknown (due to flooding or other issues) were excluded from analysis

#Reading in the data
occudata<-read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/burrow_occu_master_trimmed.csv")

#Using generalize linear models for occupancy data
#Testing two models, null model (occupancy is consistent across all sites)
#and model where occupancy varies across major habitat groups
occunull<-glm(occu~1,family="binomial", data=occudata)
occutype<- glm(occu~type, family ="binomial", data=occudata)
summary(occunull)
summary(occutype)

#By comparing delta AIC scores, occutype is superior model

#Now want to back transform parameter estimates and get 95% CIs
#Since model has multiple categorical parameters, need to build a dataframe 'predictdata' with them named
#If no covariates in the top model, can use 'backTransform' function instead
#Left in burrows from broader Gopher Tortoise surveys at the Kennedy Space Center here as examples of
#handling covariates in the data
predictdata<-with(occudata,data.frame(type=factor(c("corridor","inland","natural_dune","new_dune","old_dune"))))
#Generates predicted values
occupredict <- cbind(predictdata, predict(occutype, newdata = predictdata, type="link", se=TRUE))
#Plogis backtransforms from logistics to decimal to be able to interpret parameters
occupredict <- within(occupredict, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
occupredict

#At this point, we have top models for the following
#Burrow detection
#Burrow density
#Burrow occupancy
#All three models interact to drive overall Gopher Tortoise Density
#Now we need to unify the outputs of all three model sets


#Combining Models and Estimating Gopher Tortoise Density ----
#Function takes in the predictions for a given site, its occupancy distribution, and iteration number, 
#then generates a distribution for the site and for each iteration
#randomly samples for occupancy and density, then uses that to estimate tortoise density
#as a sum of n bernoulli trials
#after each iteration, adds sum to vector, and get distribution from the vector as estimate of 
#true tortoise density

#Now need to use those data frames and top models to predict burrow density in each site ----
Insite<-predict(gopherTortoiseBurrow_top,newdata=inlandcovs, type="state")

#And do the same for burrow occupancy in each major area
BOccuInland<-rlogis(10000,location = occupredict$fit[2], scale=((sqrt(3)/pi)*(occupredict$se.fit[2])))

TDensity<-function(Site,BOccu, n){
  result = rep(NA, n)
  #Here we used a normal distribution for tortoise burrows. Would recommend also trying a poisson distribution
  #Can evaluate your choice by doing a histogram of the number of burrows per transect
  BDensity<-rnorm(1000,Site$Predicted,Site$SE)
  for (i in 1:n) {
    trials <-as.integer(sample(BDensity,1))
    if (trials<0){
      trials<-0
      #This is needed since with standard error, density values may be estimated below zero, causing an error in later step
    }
    p<-plogis(sample(BOccu,1))
    boot.sample <- sum(sample(c(0,1), size=trials,replace=TRUE, prob=c(1-p,p)))
    result[i] <- boot.sample
  }
  return(result)
  hist(result)
}

#Now can generate distribution of possible gopher tortoise density per hectare
InlandGTD <- TDensity(Insite, BOccuInland, 10000)

#and calculate summary statistics
quantile(InlandGTD,c(0,0.025,0.5,0.975,1))
mean(InlandGTD)
sd(InlandGTD)
hist(InlandGTD)

#Can see from histogram, distribution is non-normal, thus why better to use quantiles to get intervals
#Below function takes in the area of tortoise habitat, the generation distribution from that area and number of iterations to run
#It returns the population distribution for that area
#Can change the number of iterations used to generate the population size, have 1000 set as the default
PopSize<- function(area, distri, iter=1000){
  popdist<-rep(NA,iter)
  for (i in 1:iter){
    poptrack<-rep(NA,area)
    for (j in 1:area)
    {
      poptrack[j]<-sample(distri,1)
    }
    popdist[i]<-sum(poptrack)
  }
  return(popdist)
}
#50 is just an example for area size here
InlandSize<-PopSize(50,InlandGTD)
quantile(InlandSize,c(0,0.025,0.5,0.975,1))
mean(InlandSize)
sd(InlandSize)
hist(InlandSize)

