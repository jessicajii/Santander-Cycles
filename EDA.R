## Empty list object
weekly_data = list()
## List all the .csv file names (and sort them in case they are not ordered by number)
file_names = sort(list.files('~/Documents/github/Santander-Cycles-Network/santander_bikes_data/', pattern='.csv'))
## Total number of files
n_weeks = length(file_names)
## Import the last 6 weeks of data
for(week in (n_weeks-5):n_weeks){
  weekly_data[[week]] = read.table(paste('~/Documents/github/Santander-Cycles-Network/santander_bikes_data/',file_names[week],sep=''), 
                                   sep=',', header=FALSE, col.names=c('start_id','end_id','start_time','duration'))
}
## Create a unique dataframe
library(dplyr)
df = dplyr::bind_rows(weekly_data)
str(df)

# locations
stations = read.table('~/Documents/github/Santander-Cycles-Network/santander_locations.csv', sep=',', header=TRUE)

# Distance matrix in km
library(geodist)
distances = geodist(stations[,c('longitude','latitude')], measure = 'vincenty')/1000

which(distances == max(distances), arr.ind = TRUE) 

## Most popular stations
require(data.table)
start_st <- as.data.frame(table(df$start_id))
end_st <- as.data.frame(table(df$end_id))

# Back to the starting station
x = 130
a <- which(df$start_id == x)
a1 <- as.data.frame(table(df[a,]$end_id))
a2 <- which(df$start_id == x & df$end_id == x)
isTRUE(all.equal(length(a2), max(a1$Freq)))
# FALSE:[IDs]=[3,5,19,20,22,23,24,27,30,40,41,46,(50),53,54,55,69,72,73,74,82,84,
#              102,119,120,127,128,]
# no journey start or end at ID=21,44,59,65,...

#
library(MASS)
library(survival)
library(fitdistrplus)
library(logspline)
library(plyr)
library(flexsurv) 
library(reldist) 

### Hyde Park Corner, Hyde Park
Hyde <- which(df$start_id == 191 & df$end_id == 191)
hist(df[Hyde,]$duration)

# a plot visualizing the data together with some possible 
# theoretical distributions in a skewness-kurtosis space
descdist(as.numeric(df[Hyde,]$duration), discrete = FALSE, boot=1000) # Cullen and Frey-plot
# detailed comparisons
weibullfit  <-  fitdistrplus::fitdist(as.numeric(df[Hyde,]$duration), "weibull")
lnormfit  <-  fitdistrplus::fitdist(as.numeric(df[Hyde,]$duration), "lnorm")   
gengammafit  <-  fitdistrplus::fitdist(as.numeric(df[Hyde,]$duration), "gengamma",
                                       start=function(d) list(mu=mean(d),sigma=sd(d),Q=0))
# compare the fits with qqplots
qqcomp(list( weibullfit, lnormfit, gengammafit),
       legendtext=c("weibull","lnorm", "gengamma"), fitpch = 16)
# the comparison is a relative density plot, let us use the best fitting generalized gamma 
# distribution as reference distribution. Then if the data exactly follows the reference 
# distribution, the plot will be a uniform density
N  <-  length(as.numeric(df[Hyde,]$duration))
q_rel  <-  qgengamma(ppoints(N), mu=gengammafit$estimate["mu"],
                     sigma=gengammafit$estimate["sigma"],
                     Q=gengammafit$estimate["Q"])
reldist(as.numeric(df[Hyde,]$duration), q_rel, main="Relative distribution comp to generalized gamma")

# Start from Hyde Park Corner, Hyde Park
Hyde_Park_Corner <- which(df$start_id == 191)
summary(df[Hyde_Park_Corner,])
hist(log((df[Hyde_Park_Corner,]$duration)/60), breaks = 15) 

# a plot visualizing the data together with some possible 
# theoretical distributions in a skewness-kurtosis space
descdist(log((df[Hyde_Park_Corner,]$duration)/60), discrete = FALSE, boot=1000) # Cullen and Frey-plot
# detailed comparisons
weibullfit1  <-  fitdistrplus::fitdist(log((df[Hyde_Park_Corner,]$duration)), "weibull")
lnormfit1  <-  fitdistrplus::fitdist(log((df[Hyde_Park_Corner,]$duration)), "lnorm")   
gengammafit1  <-  fitdistrplus::fitdist(log((df[Hyde_Park_Corner,]$duration)), "gengamma",
                                       start=function(d) list(mu=mean(d),sigma=sd(d),Q=0))
# compare the fits with qqplots
qqcomp(list( gammafit1, weibullfit1, lnormfit1, gengammafit1),
       legendtext=c("gamma", "lnorm", "weibull", "gengamma"), fitpch = 16)
# the comparison is a relative density plot, let us use the best fitting generalized gamma 
# distribution as reference distribution. Then if the data exactly follows the reference 
# distribution, the plot will be a uniform density
N1  <-  length(log((df[Hyde_Park_Corner,]$duration)))
q_rel1  <-  qgengamma(ppoints(N), mu=gengammafit1$estimate["mu"],
                     sigma=gengammafit1$estimate["sigma"],
                     Q=gengammafit1$estimate["Q"])
reldist(log((df[Hyde_Park_Corner,]$duration)), q_rel1, main="Relative distribution comp to generalized gamma")

end_1 <- df[Hyde_Park_Corner,]$end_id
end_1 <- as.data.frame(end_1)
dis_1 <- function(i) {geodist(stations[stations$Station.Id == 191, c('longitude','latitude')], 
                              stations[stations$Station.Id == end_1[i,] , c('longitude','latitude')], 
                              measure = 'vincenty')/1000}
Hyde_df <- data.frame(df[Hyde_Park_Corner,])
Hyde_df$distance <- lapply(seq(1:nrow(Hyde_df)), dis_1)
plot(Hyde_df$end_id, Hyde_df$distance)
plot(log(as.numeric(Hyde_df$distance) +1), log(Hyde_df$duration), xlab = "distance", ylab = "duration") # Distance vs. duration

hist(as.numeric(Hyde_df$distance), prob = TRUE)
hist(log(as.numeric(Hyde_df$distance)), prob = TRUE)
lines(density(as.numeric(Hyde_df$distance), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(as.numeric(Hyde_df$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density

as.data.frame(table(df[Hyde_Park_Corner,]$end_id)) 
summary(as.data.frame(table(df[Hyde_Park_Corner,]$end_id)))
summary(as.data.frame(table(df[Hyde_Park_Corner,]$duration)))
hist((df[Hyde_Park_Corner,]$duration)/60, col = 'skyblue3', breaks = 20)
descdist(as.numeric(df[Hyde_Park_Corner,]$duration), discrete = FALSE)

# a plot visualizing the data together with some possible 
# theoretical distributions in a skewness-kurtosis space
descdist((df[Hyde_Park_Corner,]$duration)/60, discrete = FALSE, boot=1000) # Cullen and Frey-plot
# detailed comparisons
weibullfit2  <-  fitdistrplus::fitdist((df[Hyde_Park_Corner,]$duration), "weibull")
lnormfit2  <-  fitdistrplus::fitdist((df[Hyde_Park_Corner,]$duration), "lnorm")   
gengammafit2  <-  fitdistrplus::fitdist((df[Hyde_Park_Corner,]$duration), "gengamma",
                                        start=function(d) list(mu=mean(d),sigma=sd(d),Q=0))
# compare the fits with qqplots
qqcomp(list( weibullfit2, lnormfit2, gengammafit2),
       legendtext=c( "weibull", "lnorm","gengamma"), fitpch = 16)
# the comparison is a relative density plot, let us use the best fitting generalized gamma 
# distribution as reference distribution. Then if the data exactly follows the reference 
# distribution, the plot will be a uniform density
N2 <-  length((df[Hyde_Park_Corner,]$duration))
q_rel2  <-  qweibull(ppoints(N), shape=weibullfit2$estimate["shape"],
                      scale=weibullfit2$estimate["scale"])
reldist((df[Hyde_Park_Corner,]$duration), q_rel2, main="Relative distribution comp to generalized gamma")

# End at Hyde Park Corner, Hyde Park
Hyde_Park_Corner1 <- which(df$end_id == 191)
summary(df[Hyde_Park_Corner1,])
hist((df[Hyde_Park_Corner1,]$duration), col = 'skyblue3', breaks = 20) 

start_1 <- df[Hyde_Park_Corner1,]$start_id
start_1 <- as.data.frame(start_1)
dis_11 <- function(i) {geodist(stations[stations$Station.Id == start_1[i,] , c('longitude','latitude')], 
                               stations[stations$Station.Id == 191, c('longitude','latitude')], 
                              measure = 'vincenty')/1000}
Hyde_df1 <- data.frame(df[Hyde_Park_Corner1,])
Hyde_df1$distance <- lapply(seq(1:nrow(Hyde_df1)), dis_11)
plot(Hyde_df1$start_id, Hyde_df1$distance)
plot(Hyde_df1$distance, Hyde_df1$duration) 
plot(log(as.numeric(Hyde_df1$distance)+1), log(Hyde_df1$duration)) # Distance vs. duration

hist(as.numeric(Hyde_df1$distance), prob = TRUE)
lines(density(as.numeric(Hyde_df1$distance), na.rm= TRUE), col="blue", lwd=2) # add a kernel density estimates
lines(density(as.numeric(Hyde_df1$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density

as.data.frame(table(df[Hyde_Park_Corner1,]$start_id)) # 191 2292
summary(as.data.frame(table(df[Hyde_Park_Corner1,]$start_id)))
summary(as.data.frame(table(df[Hyde_Park_Corner1,]$duration)))
hist((df[Hyde_Park_Corner1,]$duration)/3600, col = 'skyblue3', breaks = 20)
descdist(as.numeric(df[Hyde_Park_Corner1,]$duration), discrete = FALSE)

# a plot visualizing the data together with some possible 
# theoretical distributions in a skewness-kurtosis space
descdist(as.numeric(df[Hyde_Park_Corner1,]$duration), discrete = FALSE, boot=1000) # Cullen and Frey-plot
# detailed comparisons
weibullfit3  <-  fitdistrplus::fitdist(as.numeric(df[Hyde_Park_Corner1,]$duration), "weibull")
lnormfit3  <-  fitdistrplus::fitdist(as.numeric(df[Hyde_Park_Corner1,]$duration), "lnorm")   
gengammafit3  <-  fitdistrplus::fitdist(as.numeric(df[Hyde_Park_Corner1,]$duration), "gengamma",
                                        start=function(d) list(mu=mean(d),sigma=sd(d),Q=0))
# compare the fits with qqplots
qqcomp(list( weibullfit3, lnormfit3, gengammafit3),
       legendtext=c("lnorm", "weibull", "gengamma"), fitpch = 16)
# the comparison is a relative density plot, let us use the best fitting generalized gamma 
# distribution as reference distribution. Then if the data exactly follows the reference 
# distribution, the plot will be a uniform density
N3 <-  length((as.numeric(df[Hyde_Park_Corner1,]$duration)))
q_rel3  <-  qgengamma(ppoints(N), mu=gengammafit3$estimate["mu"],
                      sigma = gengammafit3$estimate["sigma"],
                     Q=gengammafit3$estimate["Q"])
reldist(as.numeric(df[Hyde_Park_Corner1,]$duration), q_rel3, main="Relative distribution comp to generalized gamma")

# Start from Queen's Gate North
Queens_Gate_North <- which(df$start_id == 266)
summary(df[Queens_Gate_North,])

end <- df[Queens_Gate_North,]$end_id
end <- as.data.frame(end)
dis <- function(i) {geodist(stations[stations$Station.Id == 266, c('longitude','latitude')], 
              stations[stations$Station.Id == end[i,] , c('longitude','latitude')], 
              measure = 'vincenty')/1000}
Queens_df <- data.frame(df[Queens_Gate_North,])
Queens_df$distance <- lapply(seq(1:nrow(Queens_df)), dis)
plot(Queens_df$end_id, Queens_df$distance)
plot(Queens_df$distance, Queens_df$duration) 
plot(log(as.numeric(Queens_df$distance)+1), log(Queens_df$duration), xlab = "distance")

hist(as.numeric(Queens_df$distance), prob = TRUE)
lines(density(as.numeric(Queens_df$distance), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(as.numeric(Queens_df$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density
  
as.data.frame(table(df[Queens_Gate_North,]$end_id))#266 53
summary(as.data.frame(table(df[Queens_Gate_North,]$end_id)))
summary(as.data.frame(table(df[Queens_Gate_North,]$duration)))
hist((df[Queens_Gate_North,]$duration)/3600, col = 'skyblue3', breaks = 20)
descdist(as.numeric(df[Queens_Gate_North,]$duration), discrete = FALSE)

g <- which(df$start_id == 266 & df$end_id == 266)
hist(df[g,]$duration)
descdist(as.numeric(df[g,]$duration), discrete = FALSE)


# Start from Imperial College, Knightsbridge
IC <- which(df$start_id == 392)
summary(df[IC,])
hist(df[IC,]$duration)

end_2 <- df[IC,]$end_id
end_2 <- as.data.frame(end_2)
dis_2 <- function(i) {geodist(stations[stations$Station.Id == 392, c('longitude','latitude')], 
                              stations[stations$Station.Id == end_2[i,] , c('longitude','latitude')], 
                              measure = 'vincenty')/1000}
IC_df <- data.frame(df[IC,])
IC_df$distance <- lapply(seq(1:1443), dis_2)
plot(IC_df$end_id, IC_df$distance)
plot(IC_df$distance, IC_df$duration) 
plot(log(as.numeric(IC_df$distance)+1), log(IC_df$duration))

hist(as.numeric(IC_df$distance), prob = TRUE)
lines(density(as.numeric(IC_df$distance), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(as.numeric(IC_df$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density

as.data.frame(table(df[IC,]$end_id)) # 392 97
summary(as.data.frame(table(df[IC,]$end_id)))
summary(as.data.frame(table(df[IC,]$duration)))
hist((df[IC,]$duration)/3600, col = 'skyblue3', breaks = 20)
descdist(as.numeric(df[IC,]$duration), discrete = FALSE)

j <- which(df$start_id == 392 & df$end_id == 392)
j
df[j,]
hist(df[j,]$duration)
descdist(as.numeric(df[j,]$duration), discrete = FALSE)

## Plots
par(mfrow = c(1, 3))
plot(log(as.numeric(Hyde_df$distance) +1), log(Hyde_df$duration), xlab = "log(distance+1) (km)", ylab = "log(duration) (seconds)", main = "Hyde Park Corner") # Distance vs. duration
plot(log(as.numeric(Queens_df$distance)+1), log(Queens_df$duration), xlab = "log(distance+1) (km)", ylab = "log(duration) (seconds)", main = "Queen's Gate North")
plot(log(as.numeric(IC_df$distance)+1), log(IC_df$duration), xlab = "log(distance+1) (km)", ylab = "log(duration) (seconds)", main = "Imperial College, Knightsbridge")
mtext("Distance vs. Duration", side = 3, line = -3, outer = TRUE)

hist(as.numeric(Hyde_df$distance), prob = TRUE, xlab = "distance (km)", main = "Hyde Park Corner")
lines(density(as.numeric(Hyde_df$distance), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(as.numeric(Hyde_df$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density

hist(as.numeric(Queens_df$distance), prob = TRUE,xlab = "distance (km)", 
     main = "Start at Queen's Gate North", breaks = 50)
lines(density(as.numeric(Queens_df$distance), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(as.numeric(Queens_df$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density

hist(as.numeric(IC_df$distance), prob = TRUE, xlab = "distance (km)", 
     main = "Start at Imperial College, Knightsbridge",breaks = 50)
lines(density(as.numeric(IC_df$distance), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(as.numeric(IC_df$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density
mtext("Distributions of Distance", side = 3, line = -3, outer = TRUE)
dev.off()

#### Hyde Park Corner
par(mfrow = c(1, 2))
plot(log(as.numeric(Hyde_df$distance) +1), log(Hyde_df$duration), xlab = "log(distance+1) (km)", ylab = "log(duration) (seconds)", main = "Start at Hyde Park Corner") # Distance vs. duration
plot(log(as.numeric(Hyde_df1$distance) +1), log(Hyde_df1$duration), xlab = "log(distance+1) (km)", ylab = "log(duration) (seconds)", main = "End at Hyde Park Corner") # Distance vs. duration
mtext("Distance vs. Duration", side = 3, line = -3, outer = TRUE)

hist(as.numeric(Hyde_df$distance), prob = TRUE, xlab = "distance (km)", 
     main = "Start at Hyde Park Corner", breaks = 50)
lines(density(as.numeric(Hyde_df$distance), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(as.numeric(Hyde_df$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density
hist(as.numeric(Hyde_df1$distance), prob = TRUE, xlab = "distance (km)", 
     main = "End at Hyde Park Corner", breaks = 50)
lines(density(as.numeric(Hyde_df1$distance), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(as.numeric(Hyde_df1$distance), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density

hist(log(as.numeric(Hyde_df$distance) +1), prob = TRUE, xlab = "log(distance+1) (km)", 
     main = "Start at Hyde Park Corner",breaks = 50)
lines(density(log(as.numeric(Hyde_df$distance) +1), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(log(as.numeric(Hyde_df$distance) +1), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density
hist(log(as.numeric(Hyde_df1$distance) +1), prob = TRUE, xlab = "log(distance+1) (km)", 
     main = "End at Hyde Park Corner",breaks = 50)
lines(density(log(as.numeric(Hyde_df1$distance) +1), na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults
lines(density(log(as.numeric(Hyde_df1$distance) +1), adjust=2, na.rm= TRUE), lty="dotted", col="darkgreen", lwd=2)  # add another "smoother" density
dev.off()
