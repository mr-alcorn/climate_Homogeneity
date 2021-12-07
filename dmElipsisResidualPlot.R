library(SWTools)
library(tidyverse)
library(lubridate)
library(zoo)
library(hydroTSM)
library(dygraphs)
# library(plotrix)
source("dyplot_calibrations.R")
climatePath <- "C:/Source/R/shDoubleMass/data"
sDate <- as.Date("1900-01-01")
eDate <- as.Date("2016-12-31")

dmSite <- c("23820") # The site to check. A dowloaded Silo txt file must exist in Climate Path
refSites <- c("23003","23756","23737","23752") # Some Comparison Sites

# Create Table of z Values for ellipsis
Zp <- data.frame(p=c(50,60,70,75,80,85,90,95),zp=c(0.00,0.25,0.52,0.67,0.84,1.04,1.28,1.64))

dates <- seq(sDate, eDate, "years")

# Function to return an annual mean of the compariosn sites
averageRefSites <- function(refSites,climatePath){
  y <- SILOLoad(refSites,climatePath,startdate =  sDate,enddate = eDate) 
  y <- lapply(y,function(x) daily2annual(x$tsd$Rain,FUN="sum"))
  y <- data.frame(matrix(unlist(y),nrow=length(y[[1]]),byrow=FALSE))
  y <- rowMeans(y)
  #y <- zoo(y,order.by = dates)
  return(y)
}
#Load the data: Note for app: Have already loaded the data here in initial Cumulative Deviation Plot
X <- SILOLoad(dmSite,climatePath,startdate =  sDate,enddate = eDate)
Xann <- lapply(X,function(x) daily2annual(x$tsd$Rain,FUN="sum")) 
Xann <- matrix(unlist(Xann),nrow=length(Xann[[1]]),byrow=FALSE)
XannZ <- zoo(Xann,order.by = dates)
Yann <- averageRefSites(refSites,climatePath)
YannZ <- zoo(Yann,order.by = dates)
n <- length(Xann)

dat <- xts(merge(XannZ,YannZ))
#Calculate the recession Y~X
# Where Y is the average of surrounding sites
#       X is the current site being assessed
linReg <- lm(Xann~Yann)
summary(linReg)
#par(col=c("blue","red","green"))
yhigh <- round(max(max(Yann),max(Xann))*1.2,-2)
plot(Yann,Xann,type = "p",main = paste0("Scatter Plot (Ave of Sites vs. ",dmSite,")"))
lines(Yann,linReg$fitted.values)


# Plot the two annual series:
dygraph(dat,main = paste("Annual Rainfall Series") ) %>%
  dyOptions(labelsUTC = TRUE, 
            fillGraph=TRUE, 
            fillAlpha=0.1, 
            drawGrid = TRUE) %>%
  dyAxis("y", label = "Cumulative residual (mm)") %>%
  dyRangeSelector() %>%
  dyCrosshair(direction = "vertical") %>%
  dyRoller(rollPeriod = 1)


eps1 <- linReg$residuals # the residuals of Y ~ X
cumRes <- cumsum(eps1) # Cumulative Sum Residuals
df <- data.frame(Year=year(dates))
df$cumRes <- cumsum(eps1)


#Set up the elipsis
p <- 80 # p is the level of confidence required
pZ <- Zp[match(p,Zp[,1]),2] # Z value at p
alpha <- n/2
beta <- (n/(n-1)^0.5)*pZ*sd(eps1)
df$y1 <- sqrt(abs((beta^2)*(1-((index(df)-alpha)/alpha)^2)))
df$y2 <- -df$y1

dpdata <- xts(df[2:4],order.by = dates)

# plot the Cumulative Residuals including the ellipsis
dygraph(dpdata,main = paste0("Cumulative Residuals: M0",dmSite," vs Surrounding stations") ) %>%
  dyOptions(labelsUTC = TRUE,
            fillGraph=TRUE,
            fillAlpha=0.1,
            drawGrid = TRUE) %>%
  dyAxis("y", label = "Cumulative Residual") %>%
  dyRangeSelector() %>%
  dyCrosshair(direction = "vertical") %>%
  dyRoller(rollPeriod = 1)


#-- 8. CHeck that the cumulative residuals lie within the elipsis.
# If they do: No further work required
 
 badResidualCount =0
   for(i in 1:n){
     
     ires <- df$cumRes[i]
     
     if(ires>df[i,3] | ires<df[i,4]){
       badResidualCount = badResidualCount +1
     }
   }
 badResPercent <- round((badResidualCount/n)*100,1)
 rsdMsg <- paste("This test shows that approx",badResPercent,
                 "% of data falls outside the computed range with a p-value of",
                 p)
 print(rsdMsg)

 # Plot the double mass to check
 plotDm <- function(dmAnn, refAnn) {
   dmX <- cumsum(refAnn)
   dmY <- cumsum(dmAnn)
   plot(dmX,dmY)
   abline(0,max(dmY)/max(dmX))
 }
 plotDm(Xann,Yann)
# If they don't: Continue on..........

#-- 9. Select break point(s) where it;
#     ceases to ceases to increase or decrease i.e.  changes direction >^<
# -- This Break point is K=i
j <- 1 # Start point for Correction (year)
k <- 68 # End Point For Correction (year)


#-- 10. Assume the line <(k=i ) is Non-Homogeneous and;
#       the line > ( k = i) is Homogeneous.
#     This yields two new regressions:
#
# nonHomog regression: yi= aNh + bNh.Xi

nonhLinReg <- lm( Xann[j:k] ~ Yann[j:k])
homogLinReg <- lm( Xann[k+1:n] ~ Yann[k+1:n])

nhYann <- data.frame(Yann[1:k])
bnh <- nonhLinReg$coefficients[2]
anh <- nonhLinReg$coefficients[1]
bh <- homogLinReg$coefficients[2]
ah <- homogLinReg$coefficients[1]

ynonHom <- anh[1] + bnh[1]*nhYann
yHom <- ah[1] + bh[1]*nhYann

df2 <- data.frame(cbind(Year=year(dates[1:k]),ynonHom,yHom))

# Plot the predicted annual totals
plot(df2$Yann.1.k.)
points(df2$Yann.1.k..1,col="green")


df2$diff <- df2$Yann.1.k..1-df2$Yann.1.k.
df2$cf <- df2$Yann.1.k..1/df2$Yann.1.k.


summary(homogLinReg)
corrFactor <- homogLinReg$coefficients[2]/nonhLinReg$coefficients[2]

nonFitLine <- nonhLinReg$model
hFitLine <- homogLinReg$model
# This doesnt work and not sure what its for
# plot(nonFitLine$`Yann[1:k]`,nonFitLine$`Xann[1:k]`)
#   abline(nonhLinReg)
#   points(hFitLine$`Yann[k + 1:n]`,hFitLine$`Xann[k + 1:n]`,col="red")
#   abline(homogLinReg,col="red")

# Homog regression: yi= ah + bh.Xi

#-- 11. Compute the difference between the two regressions
# i.e  delta Yi = {ah+bh.Xi} - {aNh+bNh.Xi} for the non-homogenous set. i<=k

plot(nonFitLine, col = "red")
 points(hFitLine,col="blue")
  
 

#--- 12. Correct the non-homgenous subset portion of data set.
#       Yci = Yi + deltaYi for i <= k  (Yci = the corrected Y values)
#     
corXann <- Xann
for (i in 1:k){
  corXann[i] <- Xann[i]+df2$diff[i]
}
#CorrectedSet[k+1:n] <- Xann[k+1:n]
plot(Yann, type = "l") # All Original Annual Totals
  lines(Xann, col="red")
  lines(corXann,col="blue")

#--- Lucky 13! At this point can run the homogeneity check again to confirm that 
#               the test is satisfied post correction.

 #  This section doesnt work.....
# dmDaily <- X$`23878`$tsd$Rain
# corDmDaily <- dmDaily
# corYears <- df2$Year
# #length(dmDaily
# for (i in 1:length(dmDaily)) {
#   iYear <- year(index(dmDaily[i]))
#   if(iYear %in% corYears){
#     cf <- df2$cf[index(match(iYear,df2$Year))]
#     corDmDaily[[i]] <- dmDaily[[i]]*cf
#   }
# }


# Function to Correct the daily data using Annual Factor
correctDailyData <- function(X, df2) {
  dmDaily <- X[[1]]$tsd$Rain
  corDmDaily <- dmDaily
  corYears <- df2$Year
  #length(dmDaily
  for (i in 1:length(dmDaily)) {
    iYear <- year(index(dmDaily[i]))
    if(iYear %in% corYears){
      cf <- df2$cf[index(match(iYear,df2$Year))]
      corDmDaily[[i]] <- dmDaily[[i]]*cf
    }
  }
  return(corDmDaily)
}
# New Daily Data
d <- correctDailyData(X,df2)
# New Annual Data
corDmAnn <- as.vector(daily2annual(d,FUN = "sum"))
# Plot the Double Mass of the Corrected Data as a  Check
plotDm(corDmAnn,Yann)
# Write the Double Mass corrected data to file
write.zoo(d,paste0("Corrected-",dmSite,"-Rain.csv"),sep=",")

