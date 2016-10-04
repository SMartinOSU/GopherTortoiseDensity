library(unmarked)
library(AICcmodavg)
library(aod)
library(ggplot2)
library(Rcpp)
library(lme4)
# [Complete] Importing Data into R ----
#First need to read in data from distance surveys (see files for how to set up), then separately read
#in transect covariates

#Reading in data
#use getwd() to view folder R is importing from
#use setwd("folder address") to change target folder
#data needs to be in .csv (comma deliminated) or tab deliminated format. Can save from excel or other spreadsheet program
#Due to large separation of landcover types, splitting data into three main types (inland areas, corridor areas, beach areas)
masterdists <- read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/BurrowDist_Master_Up.csv")
mastercovs <- read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/siteCov_Master_Up.csv")
beachlength<-read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/Transects_beachLength.csv")
beachdists<-read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/BurrowDist_Beach_Up.csv")
inlanddists<-read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/BurrowDist_Inland_Up.csv")
corridordists<-read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/BurrowDist_Corridor_Up.csv")
inlandcovs<-mastercovs[1:200,]
inlandlength<-mastercovs[1:200,7]
corridorcovs<-mastercovs[201:300,]
corridorlength<-mastercovs[201:300,7]
beachcovs<-mastercovs[301:350,]
beachlength<-mastercovs[301:350,7]

#First file is just each observation in a two column format. 
#First column is perpendicular distance to transect
#Second column is the transect it was found on. For empty transects leave a blank in the first column
#See file for example on data structure

#Second file has one transect per row with its potential covariates in each column. Landcover data has already been standardized

# [Complete] Formating your data for Unmarked----
#Make a vector that contains the length of each transect, in the SAME order they are set up in the 'covs' file
#If transect length varies, can read lengths in from a file and use it that way (beach length is example of this)
#If all transects are same length can simply make a vector using the rep() command
length<-c(inlandlength, corridorlength, beachlength)

# Checking covariates for correlations
correlations<-cor(mastercovs[,8:25])
View(correlations)
#Open_H20 with mangrove and barren (use open h20), Ocean with beach (use beach) are correlated in the master dataset
#Have a priori standard for what is too high of a correlation, here used r >= 0.6

#Need to format data. 'breaks' is the bin range and spacing for distance models (so bins range from 0 to 30 meters in 0.8m increments)
#yDat and umf are from unmarked, preparing data for models
#umf (UnMarkedFunction) is master data, bumf (BeachUnMarkedFunction) is beach subset, 
#iumf (InlandUnMarkedFunction) is inland subset, cumf (CorridorUnMarkedFunction) is corridor subset
#Refer to unmarked help using ?unmarkedFrame() for double observer datasets
breaks<-seq(0,30,0.8)

yDat<-formatDistData(masterdists,distCol="dist", transectNameCol="transect",dist.breaks=breaks)
umf<-unmarkedFrameDS(y=as.matrix(yDat), siteCovs=mastercovs,survey="line",
                     dist.breaks=breaks,tlength=length,unitsIn="m")

bDat<-formatDistData(beachdists,distCol="dist", transectNameCol="transect",dist.breaks=breaks)
bumf<-unmarkedFrameDS(y=as.matrix(bDat), siteCovs=beachcovs,survey="line",
                     dist.breaks=breaks,tlength=beachlength,unitsIn="m")

iDat<-formatDistData(inlanddists,distCol="dist", transectNameCol="transect",dist.breaks=breaks)
iumf<-unmarkedFrameDS(y=as.matrix(iDat), siteCovs=inlandcovs,survey="line",
                     dist.breaks=breaks,tlength=inlandlength,unitsIn="m")

cDat<-formatDistData(corridordists,distCol="dist", transectNameCol="transect",dist.breaks=breaks)
cumf<-unmarkedFrameDS(y=as.matrix(cDat), siteCovs=corridorcovs,survey="line",
                     dist.breaks=breaks,tlength=corridorlength,unitsIn="m")

#Visualizing data, generates a histogram of frequency over distance. 
#If looks off, try adjusting bin size using the breaks variable
hist(umf, xlab="distance (m)", cex.lab=0.8, cex.axis=0.8,main="Master")
hist(iumf, xlab="distance (m)", cex.lab=0.8, cex.axis=0.8,main="Inland")
hist(cumf, xlab="distance (m)", cex.lab=0.8, cex.axis=0.8,main="Corridor")
hist(bumf, xlab="distance (m)", cex.lab=0.8, cex.axis=0.8,main="Coastal")

# [Complete] Selecting Detection Distributions ----
#Testing what of the four main potential distributions best fit the data using null models
haz_Null<-distsamp(~1 ~1,umf, keyfun="hazard",output="density",unitsOut="ha")
halfnorm_Null<-distsamp(~1 ~1,umf, keyfun="halfnorm",output="density",unitsOut="ha")
exp_Null<-distsamp(~1 ~1,umf, keyfun="exp",output="density",unitsOut="ha")
uni_Null<-distsamp(~1 ~1,umf, keyfun="uniform",output="density",unitsOut="ha")

#now comparing models to select best fit based on AIC values
modSel(fitList(haz_Null,halfnorm_Null,exp_Null,uni_Null))

#With our data, exponential best fit, will continue to use throughout models

# [Complete] Burrow Detection Models ----
#Now need to start on detection covariate testing, formatting follows unmarked user manual.
#See ?distamp() for details on model construction. 
#'+' indicates additive effect, '*' for including interaction term in models
#no a priori reason for any interactions here, and they greatly increase number of parameters
#Not testing for effects of veg on density, since few landcover types found across majority of transects

#Below models are for the master dataset. Change 'umf' to appropriate model to test others.
#Be aware of parameters when testing bumf/cumf/iumf
det_m1<-distsamp(~1 ~1, umf, keyfun="exp",output="density",unitsOut="ha")
det_m2<-distsamp(~obs ~1, umf, keyfun="exp",output="density",unitsOut="ha")
det_m3<-distsamp(~loc ~1, umf, keyfun="exp",output="density",unitsOut="ha")
det_m4<-distsamp(~season ~1, umf, keyfun="exp",output="density",unitsOut="ha")
det_m5<-distsamp(~type ~1, umf, keyfun="exp",output="density",unitsOut="ha")
det_m6<-distsamp(~obs+loc ~1, umf, keyfun="exp",output="density",unitsOut="ha")
det_m7<-distsamp(~obs+season ~1, umf, keyfun="exp",output="density",unitsOut="ha")

#Now need to compare model set to pick top model. modSel function ranks them in ascending order
modSel(fitList(det_m1,det_m2,det_m3,det_m4,det_m5,det_m6, det_m7))
#Change umf to others, and test per-site differences in detection covariates
#Be aware, not all models should be tested for each region; some covariates only make sense for master list
#Beach data has only a single observer for example, so using 'obs' covariate will cause an error
#Model 4 fits best for inland (season)
#Model 2 fit best for corridor (obs)
#model 1 (null) is best for beach
#Model 6 fits best for master data

# [Complete] Burrow Density Modelling ----
#Make sure when setting up models, each makes biological sense
#For example, should not use both 'loc' and 'type' in a model,  
#since 'loc' is a more specific version of  'type' in this data. 
#Doing so generates 'NAs' for parameter estimates if we do this
#these try to model burrow density using only the master dataset

fm1<-distsamp(~loc + obs ~1, umf, keyfun="exp",output="density",unitsOut="ha")
fm2<-distsamp(~loc + obs ~PER_HAMMOCK + PER_STRAND+PER_OCEAN+PER_INFRA+PER_PEPPER+PER_MARSH+
               PER_OAK_SCRUB+PER_OPEN_H20+PER_PALM_SCRUB+
              PER_PINE_FLAT+PER_CONSTRUCT+PER_RUD_HERB+PER_RUD_WOOD
             +PER_SCRUB_FLAT+PER_WETLAND_SCRUB, umf, keyfun="exp",output="density",unitsOut="ha")
fm3<-distsamp(~loc + obs ~PER_WETLAND_SCRUB, umf, keyfun="exp",output="density",unitsOut="ha")
fm4<-distsamp(~loc + obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB
             +PER_OAK_SCRUB+loc+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
            +PER_SCRUB_FLAT, umf, keyfun="exp",output="density",unitsOut="ha")
fm5<-distsamp(~loc + obs ~PER_BEACH*PER_CONSTRUCT+PER_OPEN_H20+PER_STRAND
              +PER_PEPPER+PER_RUD_HERB+PER_PALM_SCRUB, umf, keyfun="exp",output="density",unitsOut="ha")
fm6<-distsamp(~loc + obs ~PER_BEACH+PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB+
                PER_OAK_SCRUB+loc+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT
              +PER_HAMMOCK+PER_SCRUB_FLAT+PER_CONSTRUCT+PER_OPEN_H20
             +PER_STRAND, umf, keyfun="exp",output="density",unitsOut="ha")
fm7<-distsamp(~loc + obs ~PER_BEACH+PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB+
                PER_OAK_SCRUB+type+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT
              +PER_HAMMOCK+PER_SCRUB_FLAT+PER_CONSTRUCT+PER_OPEN_H20
              +PER_STRAND, umf, keyfun="exp",output="density",unitsOut="ha")
fm8<-distsamp(~loc + obs ~PER_BEACH+PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB+
                PER_OAK_SCRUB+loc+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT
              +PER_HAMMOCK+PER_SCRUB_FLAT+PER_CONSTRUCT+PER_OPEN_H20
              +PER_STRAND, umf, keyfun="exp",output="density",unitsOut="ha")

#If models did not converge, need to adjust for that. Cannot compare non-converging models to others
#If models do not converge, or produce NAs for values, can try changing maxit, or giving starting values
#See below model for example with changed maxit
#May also get an an error 'Hessian is singular' this is likely due to one covariate containing a subset of information found in another covariate

fm9<-distsamp(~loc + obs  ~PER_BEACH+PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB+
                PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT
             +PER_HAMMOCK+PER_SCRUB_FLAT+PER_CONSTRUCT+PER_OPEN_H20
              +PER_STRAND, umf, keyfun="exp",output="density",unitsOut="ha", control = list(maxit=100))

#For adding starting values, zero is a good value, 1 for potential positive effect, -1 for negative effect
#Need to know total number of parameters when inputing starts
#See below example

fm10<-distsamp(~loc + obs  ~PER_BEACH+PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB+
                 PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT
              +PER_HAMMOCK+PER_SCRUB_FLAT+PER_CONSTRUCT+PER_OPEN_H20
               +PER_STRAND, umf, keyfun="exp",output="density",unitsOut="ha", starts=rep(0,26))

modSel(fitList(fm1,fm2,fm3,fm4,fm5,fm6,fm7,fm9,fm10))

#Model 6 is best model
#now need to test fit
fitfm6 <- Nmix.gof.test(fm6,nsim=25,plot.hist=TRUE)
fitfm6
#c-hat is a measure of model dispersion. Want to be close to 1
#If much greater than 1, need to correct for error
#also simulates predicted values (change nsim to higher values for better precision, low value to run fast) to see if predicted matches observed. 
#Recommend using low value for nsim during initial testing for quick check of fit, 
#and increasing value for precision later
#Want a high p-value, indicating observed and predicted are not significantly different
#If have low p or high c-hat, consider alternative models. May be missing a source of variation
#master top model shows bad fit (p<0.05), due to large number of zeros in land cover data 
#from variety of habitat areas. Will test sub groups of inland/corridor and beach to improve fit

#The top beach model is from prior research
beach_burrow_density<-distsamp(~1 ~PER_BEACH*PER_CONSTRUCT+PER_OPEN_H20+PER_STRAND
                               +PER_PEPPER+PER_RUD_HERB+PER_PALM_SCRUB, bumf, 
                               keyfun="exp", output="density", unitsOut="ha")
#Inland Model Construction
im1<-distsamp(~season ~PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB
             +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
             +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
im2<-distsamp(~season ~PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
im3<-distsamp(~season ~PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
im4<-distsamp(~season ~PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
im5<-distsamp(~season ~PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
im6<-distsamp(~season ~PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha", starts=rep(0,50))
im7<-distsamp(~season ~PER_RUD_WOOD+season+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
im8<-distsamp(~season ~PER_RUD_WOOD+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
im9<-distsamp(~season+obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
im10<-distsamp(~obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_PINE_FLAT+PER_HAMMOCK
              +PER_SCRUB_FLAT, iumf, keyfun="exp",output="density",unitsOut="ha")
#im6 gives Hessian error, may have too many covariates. Cannot compare non-converging models to converging models so we dropped it
modSel(fitList(im1,im2,im3,im4,im5,im7,im8,im9,im10))
#im9 & im10 both potentially top models. Only differ by detection covariates. Will assess fit for each model before continuing. 
fitim10 <- Nmix.gof.test(im10,nsim=50,plot.hist=TRUE)
fitim9 <- Nmix.gof.test(im10,nsim=50,plot.hist=TRUE)
#red line indicates where observed values are, want p-value > 0.05 since that means observed data is not statistically 
#different than predicted data
#im10 has better fit, and lower c-hat than im9, so going to use it for later analysis

#Corridor model construction
cm1<-distsamp(~obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+season+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_HAMMOCK
              , cumf, keyfun="exp",output="density",unitsOut="ha")
cm2<-distsamp(~obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_PALM_SCRUB+PER_MARSH+PER_HAMMOCK
              , cumf, keyfun="exp",output="density",unitsOut="ha")
cm3<-distsamp(~obs ~PER_WETLAND_SCRUB+PER_RUD_WOOD+PER_RUD_HERB
              +PER_OAK_SCRUB+PER_MARSH+PER_HAMMOCK
              , cumf, keyfun="exp",output="density",unitsOut="ha")
cm4<-distsamp(~obs ~PER_MARSH
              , cumf, keyfun="exp",output="density",unitsOut="ha")
cm5<-distsamp(~obs ~PER_RUD_HERB
              , cumf, keyfun="exp",output="density",unitsOut="ha")
cm6<-distsamp(~obs ~1
              , cumf, keyfun="exp",output="density",unitsOut="ha")
cm7<-distsamp(~obs ~season
              , cumf, keyfun="exp",output="density",unitsOut="ha")
cm8<-distsamp(~obs ~PER_RUD_WOOD+PER_RUD_HERB
              +PER_HAMMOCK, cumf, keyfun="exp",output="density",unitsOut="ha")
cm9<-distsamp(~obs ~ +PER_MARSH, cumf, keyfun="exp",output="density",unitsOut="ha")
cm10<-distsamp(~obs ~PER_WETLAND_SCRUB
               , cumf, keyfun="exp",output="density",unitsOut="ha")
cm11<-distsamp(~obs ~PER_RUD_WOOD
              , cumf, keyfun="exp",output="density",unitsOut="ha")

modSel(fitList(cm1,cm2,cm3,cm4,cm5,cm6,cm7,cm8,cm9,cm10,cm11))
fitcm10 <- Nmix.gof.test(cm10,nsim=500,plot.hist=TRUE)

#Can see in each case top model has predicted values overlapping with observed values after separating out 
#main habitat areas. Importance of noticing completely separated covariates

#Once we have our top models, need to look at parameter estimates for the covariates
#renaming them for ease of identification
beach_top<-beach_burrow_density
corridor_top <-cm10
inland_top <- im10

#calculate asymtotic profile CI, more accurate than just using SE to estimate, will default to 95% CI
b_paraciint <- confint(beach_top, type='state', method = 'profile')
c_paraciint <- confint(corridor_top, type='state', method = 'profile')
i_paraciint <- confint(inland_top, type='state', method = 'profile')

#This will plot out the detection function
#if you want to look at different detection functions, just change the model name
modelrate<-coef(inland_top, type='det')
#pulls out the estimates just for detection, below plot is based on the intercept alone. 
plot(function(x) gxexp(x, rate=modelrate[1]), 0, 30,
     xlab="Distance (m)", ylab="Detection prob.")

#if top model had a covariate for detection, need to combine modelrate as below
#will allow you to view alternative detection curves, should always be rate_index (1) + rate_index (n)
#unless contrast matrix was changed by user (outside the scope of this code)
plot(function(x) gxexp(x, rate=(modelrate[1]+modelrate[2])), 0, 30,
     xlab="Distance (m)", ylab="Detection prob.")
#Could use add line instead of plot to overlap detection curves, but not displayed here

# [Complete] Burrow Occupancy Models ----------
#Here we will model the occupancy rates of burrow from data collected via scoping
#Burrow data has been coded as a binomial response where 1=Occupied, 0=Unoccupied.
#Burrows previously classifed as unknown (due to flooding or other issues) were excluded from analysis

#Reading in the data
occudata<-read.csv("C:/Users/smart/Dropbox/KSC/GopherTortoiseSurveyManuscript/burrow_occu_master_trimmed.csv")

#Using generalize linear models for occupancy data
#Testing two models, null model (occupancy is consistent across all sites)
#and model where occupancy varies across major habitat groups
#See (thesis paper) for more details on beach dune data
occunull<-glm(occu~1,family="binomial", data=occudata)
occutype<- glm(occu~type, family ="binomial", data=occudata)
summary(occunull)
summary(occutype)

#By comparing delta AIC scores, occutype is superior model

#Now want to back transform parameter estimates and get 95% CIs
#Since model has multiple categorical parameters, need to build a dataframe 'predictdata' with them named
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

# Review:
#At this point, we have top models for the following
#Burrow detection
#Burrow density
#Burrow occupancy
#All three models interact to drive overall Gopher Tortoise Density
#Now we need to unify the outputs of all three model sets

# [Complete] Combining Models and Estimating Gopher Tortoise Density ----
#Function takes in the predictions for a given site, its occupancy distribution, and iteration number, 
#then generates a distribution for the site and for each iteration
#randomly samples for occupancy and density, then uses that to estimate tortoise density
#as a sum of n bernoulli trials
#after each iteration, adds sum to vector, and get distribution from the vector as estimate of 
#true tortoise density

#First need to contruct new data frame representing average of habitat coverage for each major site

#Need to deal with different sized transects along the beach 
#Will split up into n 100m segements (where n = original distance of transect/100) to have representation equal to total distance of original transects

northnat<-beachcovs[1:2,]
northnat<-northnat[rep(seq_len(nrow(northnat)), each=6),]
northcon<-beachcovs[11:12,]
northcon<-northcon[rep(seq_len(nrow(northcon)),each=2),]
conold<-beachcovs[21:22,]
southcon<-beachcovs[31:32,]
southcon<-southcon[rep(seq_len(nrow(southcon)),each=9),]
southnat<-beachcovs[41:42,]
southnat<-southnat[rep(seq_len(nrow(southnat)), each=6),]
beachpredictor<-rbind(northnat,northcon,conold,southcon,southnat)

#Now need to use those data frames and top models to predict burrow density in each site ----
Beach_Site <- predict(beach_top,newdata=beachpredictor, type = "state")
Insite<-predict(inland_top,newdata=inlandcovs, type="state")
CSite<-predict(corridor_top,newdata=corridorcovs,type="state")

#And do the same for burrow occupancy in each major area
BOccuCorr<-rlogis(10000,location = occupredict$fit[1], scale=((sqrt(3)/pi)*(occupredict$se.fit[1])))
BOccuInland<-rlogis(10000,location = occupredict$fit[2], scale=((sqrt(3)/pi)*(occupredict$se.fit[2])))
BOccuBeach<-rlogis(10000,location = ((occupredict$fit[3]*2+occupredict$fit[4]*2+occupredict$fit[5])/5), scale=((sqrt(3)/pi)*((occupredict$fit[3]*2+occupredict$fit[4]*2+occupredict$fit[5])/5)))

TDensity<-function(Site,BOccu, n){
  result = rep(NA, n)
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
BeachGTD <- TDensity(Beach_Site,BOccuBeach,10000)
CorridorGTD <- TDensity(CSite,BOccuCorr,10000)
InlandGTD <- TDensity(Insite, BOccuInland, 10000)
#and calculate summary statistics
##aiming for mean around 3.28 according to thesis data
quantile(BeachGTD,c(0,0.025,0.5,0.975,1))
mean(BeachGTD)
sd(BeachGTD)
hist(BeachGTD)

quantile(InlandGTD,c(0,0.025,0.5,0.975,1))
mean(InlandGTD)
sd(InlandGTD)
hist(InlandGTD)

quantile(CorridorGTD,c(0,0.025,0.5,0.975,1))
mean(CorridorGTD)
sd(CorridorGTD)
hist(CorridorGTD)
#Can see from histogram, distribution is non-normal, thus why better to use quantiles to get intervals
#Note: These are density of tortoises per hectare, need to multiple by area of site for final population size estimates!!
#To get range for population size of gopher tortoises at size, would then just multiple density per hectare by size of site