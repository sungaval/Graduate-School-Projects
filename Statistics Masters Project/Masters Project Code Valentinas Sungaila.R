# Valentinas Sungaila ##
## Masters Project ##

#---------------------------------------------------------------------------------------------------------------------------
########################################
### trying to work with metabolites data###
#########################################


# loading my metabolites csv file

metabolites_dat <- read.csv("C:/Users/sunga/Desktop/stats thesis project/metabolites_multi_testing.csv")


# lets look at some descriptive stats
summary(metabolites_dat)
head(metabolites_dat, n = 2)
str(metabolites_dat)

# my most useful function
library(psych)
describe(metabolites_dat)

# I seem to have 40 observations with 4897 variables

# let me standerdize my OG data set to be able to apply the GAo's method

for(k in 1: ncol(metabolites_dat)){metabolites_dat[,k] <- as.numeric(as.character(metabolites_dat[,k]))}

metabolites_dat_std<-(metabolites_dat-rep(apply(metabolites_dat,2,mean,na.rm=T),
                                          each=nrow(metabolites_dat)))/rep(apply(metabolites_dat,2,sd,na.rm=T),each=nrow(metabolites_dat))

## check mean 0, var 1
summary(apply(metabolites_dat_std,2, mean))
summary(apply(metabolites_dat_std,2,var))
# looks good

# now I need to remove all the non-numeric non metabolites columns from my dataset

drops = c('QC.num','sample','id.basic','batch','order','subject','pre.post','intervention','id')
just_metabolites_std = metabolites_dat_std[ , !(names(metabolites_dat_std) %in% drops)]
View(just_metabolites_std)

summary(apply(just_metabolites_std,2, mean))
summary(apply(just_metabolites_std,2,var))

# looking great, but I still have 71 NA's which I need to remove

# one way of removing columns with NA values in R
just_metabolites_std = just_metabolites_std[,colSums(is.na(just_metabolites_std)) == 0]

# lets check now 
summary(apply(just_metabolites_std,2, mean))
summary(apply(just_metabolites_std,2,var))
# great, no more NA values

# now to be able to get the eigenvalues, lets make a covariance matrix
metabolites_cov = cov(just_metabolites_std)
head(metabolites_cov)


# ---------------------------------------------------------------------------------------------------------

### Lets create some density plots of the first nine metabolites to see the distribution

library(ggplot2)



# graph of first metabolite in the data
plot1 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$PC.34.2..)) + 
  stat_density() + 
  xlab("Metabolite PC.34.2..") +
  ylab("Density") +
  ggtitle("Density for Metabolite PC.34.2..")


# graph of second metabolite in the data
plot2 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$PC.22.6.4Z.7Z.10Z.13Z.16Z.19Z..16.0.)) + 
  stat_density() + 
  xlab("Metabolite PC.22.6.4Z.7Z.10Z.13Z.16Z.19Z..16.0.") +
  ylab("Density") +
  ggtitle("Density for Metabolite PC.22.6.4Z.7Z.10Z.13Z.16Z.19Z..16.0.")


# graph of third metabolite in the data
plot3 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$LysoPC.16.0..)) + 
  stat_density() + 
  xlab("Metabolite LysoPC.16.0..") +
  ylab("Density") +
  ggtitle("Density for Metabolite LysoPC.16.0..")


# graph of fourth metabolite in the data
plot4 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$PC.16.0.18.1.9Z...S.)) + 
  stat_density() + 
  xlab("Metabolite PC.16.0.18.1.9Z...S.") +
  ylab("Density") +
  ggtitle("Density for Metabolite PC.16.0.18.1.9Z...S.")


# graph of fifth metabolite in the data
plot5 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$Dioleoylphosphatidylcholine)) + 
  stat_density() + 
  xlab("Metabolite Dioleoylphosphatidylcholine") +
  ylab("Density") +
  ggtitle("Density for Metabolite Dioleoylphosphatidylcholine")


# graph of sixth metabolite in the data
plot6 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$PC.16.0.20.4.8Z.11Z.14Z.17Z..)) + 
  stat_density() + 
  xlab("Metabolite PC.16.0.20.4.8Z.11Z.14Z.17Z..") +
  ylab("Density") +
  ggtitle("Density for Metabolite PC.16.0.20.4.8Z.11Z.14Z.17Z..")


# graph of 7th metabolite in the data
plot7 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$SM.d18.0.16.1.9Z..)) + 
  stat_density() + 
  xlab("Metabolite SM.d18.0.16.1.9Z..") +
  ylab("Density") +
  ggtitle("Density for Metabolite SM.d18.0.16.1.9Z..")



# graph of 8th metabolite in the data
plot8 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$LysoPC.18.0..)) + 
  stat_density() + 
  xlab("Metabolite LysoPC.18.0..") +
  ylab("Density") +
  ggtitle("Density for Metabolite LysoPC.18.0..")


# graph of 9th metabolite in the data
plot9 = ggplot(data = just_metabolites_std, aes(x = just_metabolites_std$PC.14.0.22.1.13Z..)) + 
  stat_density() + 
  xlab("Metabolite PC.14.0.22.1.13Z..") +
  ylab("Density") +
  ggtitle("Density for Metabolite PC.14.0.22.1.13Z..")



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



### code for the function can be found here http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/


# saving it as a png picture so I can use it in my paper
png(paste(format(Sys.time(), "AAA"), 'png', sep = '.'),
    width = 12,
    height = 10,
    units = 'in',
    res = 300)

# the place the picture is being saved at
setwd("C:/Users/sunga/Desktop/stats thesis project/density graphs")

multiplot(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, cols=3)


dev.off()



#----------------------------------------------------------------------------------------------------------
# # lets create a heat map to see how corelated are my metabolites with one another. 
# 
metabolites_cov_matrix = as.matrix(metabolites_cov)


library(gplots)

# saving it as a png picture so I can use it in my paper
png(paste(format(Sys.time(), "%Y_%m-%d 11"), 'png', sep = '.'),
    width = 12,
    height = 10,
    units = 'in',
    res = 300)

# the place the picture is being saved at
setwd("C:/Users/sunga/Desktop/stats thesis project")

heatmap.2(metabolites_cov,
          Colv = NA,
          Rowv = NA,
          #lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 ), # Place the color key to the top right of the image plot
          #keysize = 1.5,
          margins = c(4,4), # change how close the heatmap is to border
          cexRow = 0.80, #change the font size of row variables
          cexCol = 0.80, #change the font size of column variable
          main="Metabolites Covariance Heatmap",
          scale = "none",
          col = bluered(100),
          trace = "none",
          density.info = "none",
          key.ylab = 'NA')

dev.off()

#---------------------------------------------------------------------------------------------------------------

# now lets calculate the eigenvalues using eigen function with the covariance matrix from above
eigenvalues_cov = eigen(metabolites_cov, symmetric=T, only.values=T)
summary(eigenvalues_cov$values)

##  Gao's Method  ##
sum_eigen = sum(eigenvalues_cov$values)

e2 = cbind(eigenvalues_cov$values, prop=eigenvalues_cov$values/sum_eigen)

prop_total = e2[1,2]

if(prop_total >= .995){
  prop_total = prop_total
  
}else{
  
  for (i in 2:length(e2[,2])) {prop_total = c(prop_total, prop_total[i-1]+e2[i,2])
  
  if(prop_total[i] >= .995) {break}}
  
}


M_gao_cov = length(prop_total)

M_gao_cov




# ______________________________________________________________________________________________________________________






#############################################
### trying gaos method on simulated data ###
############################################


# start by creating a function where I can input different number of vairbles, observations, and correlations and run the function
# for gaos method

set.seed(1)
library(MASS)

gao.sim.func<-function(N, p, mu, r){
  
  x <- matrix(nrow=p, ncol=p, r)
  
  diag(x)<-1
  
  dat.r1 <- scale(mvrnorm(N, mu, x))
  
  cor_dat_1 = cor(dat.r1)
  
  eigenvalues_cov = eigen(cor_dat_1, symmetric=T, only.values=T)
  
  ##  Gao's Method  ##
  sum_eigen = sum(eigenvalues_cov$values)
  
  e2 = cbind(eigenvalues_cov$values, prop=eigenvalues_cov$values/sum_eigen)
  
  prop_total = e2[1,2]
  
  if(prop_total >= .995){
    prop_total = prop_total
    
  }else{
    
    for (i in 2:length(e2[,2])) {prop_total = c(prop_total, prop_total[i-1]+e2[i,2])
    
    if(prop_total[i] >= .995) {break}}
    
  }
  
  
  M_gao_cov = length(prop_total)
  
  return(M_gao_cov)
  
}



# let me create a loop for running my function 10 times for each set of N and p because since you are sampling from a distribution
# and gao method method has a value of .995 value we might get different results for our Gao method depending on the simulation so I will
# run each number of variables and the number of predictors 10 times with different correlation values and average them to see what I get



# star by simulating gao with 0 correlation
# trying to see how many 
# number of variables = 10 changing number of observations and cor = 0

# N = 10 P = 10 r = 0
sim_cor0_n10_p10 = numeric(10)
for (i in 1:10){
  sim_cor0_n10_p10[i] = gao.sim.func(N=10, p=10, mu=rep(0,10), r=0)
}
sim_cor0_n10_p10

# N = 100 P = 10 r = 0
sim_cor0_n100_p10 = numeric(10)
for (i in 1:10){
  sim_cor0_n100_p10[i] = gao.sim.func(N=100, p=10, mu=rep(0,10), r=0)
}
sim_cor0_n100_p10

# N = 1000 P = 10 r = 0
sim_cor0_n1000_p10 = numeric(10)
for (i in 1:10){
  sim_cor0_n1000_p10[i] = gao.sim.func(N=1000, p=10, mu=rep(0,10), r=0)
}
sim_cor0_n1000_p10

# N = 10000 P = 10 r = 0
sim_cor0_n10000_p10 = numeric(10)
for (i in 1:10){
  sim_cor0_n10000_p10[i] = gao.sim.func(N=10000, p=10, mu=rep(0,10), r=0)
}
sim_cor0_n10000_p10



# number of variables = 100 changing number of observations and cor = 0

# N = 10 P = 100 r = 0
sim_cor0_n10_p100 = numeric(10)
for (i in 1:10){
  sim_cor0_n10_p100[i] = gao.sim.func(N=10, p=100, mu=rep(0,100), r=0)
}
sim_cor0_n10_p100

# N = 100 P = 100 r = 0
sim_cor0_n100_p100 = numeric(10)
for (i in 1:10){
  sim_cor0_n100_p100[i] = gao.sim.func(N=100, p=100, mu=rep(0,100), r=0)
}
sim_cor0_n100_p100

# N = 1000 P = 100 r = 0
sim_cor0_n1000_p100 = numeric(10)
for (i in 1:10){
  sim_cor0_n1000_p100[i] = gao.sim.func(N=1000, p=100, mu=rep(0,100), r=0)
}
sim_cor0_n1000_p100

# N = 10000 P = 100 r = 0
sim_cor0_n10000_p100 = numeric(10)
for (i in 1:10){
  sim_cor0_n10000_p100[i] = gao.sim.func(N=10000, p=100, mu=rep(0,100), r=0)
}
sim_cor0_n10000_p100



# number of variables = 1000 changing number of observations and cor = 0

# N = 10 P = 1000 r = 0
sim_cor0_n10_p1000 = numeric(10)
for (i in 1:10){
  sim_cor0_n10_p1000[i] = gao.sim.func(N=10, p=1000, mu=rep(0,1000), r=0)
}
sim_cor0_n10_p1000

# N = 100 P = 1000 r = 0
sim_cor0_n100_p1000 = numeric(10)
for (i in 1:10){
  sim_cor0_n100_p1000[i] = gao.sim.func(N=100, p=1000, mu=rep(0,1000), r=0)
}
sim_cor0_n100_p1000

# N = 1000 P = 1000 r = 0
sim_cor0_n1000_p1000 = numeric(10)
for (i in 1:10){
  sim_cor0_n1000_p1000[i] = gao.sim.func(N=1000, p=1000, mu=rep(0,1000), r=0)
}
sim_cor0_n1000_p1000

# N = 10000 P = 1000 r = 0
sim_cor0_n10000_p1000 = numeric(10)
for (i in 1:10){
  sim_cor0_n10000_p1000[i] = gao.sim.func(N=10000, p=1000, mu=rep(0,1000), r=0)
}
sim_cor0_n10000_p1000



# number of variables = 4000 changing number of observations and cor = 0

# N = 10 P = 4000 r = 0
sim_cor0_n10_p4000 = numeric(10)
for (i in 1:10){
  sim_cor0_n10_p4000[i] = gao.sim.func(N=10, p=4000, mu=rep(0,4000), r=0)
}
sim_cor0_n10_p4000

# N = 100 P = 4000 r = 0
sim_cor0_n100_p4000 = numeric(10)
for (i in 1:10){
  sim_cor0_n100_p4000[i] = gao.sim.func(N=100, p=4000, mu=rep(0,4000), r=0)
}
sim_cor0_n100_p4000

# N = 1000 P = 4000 r = 0
sim_cor0_n1000_p4000 = numeric(10)
for (i in 1:10){
  sim_cor0_n1000_p4000[i] = gao.sim.func(N=1000, p=4000, mu=rep(0,4000), r=0)
}
sim_cor0_n1000_p4000

# N = 10000 P = 4000 r = 0
sim_cor0_n10000_p4000 = numeric(10)
for (i in 1:10){
  sim_cor0_n10000_p4000[i] = gao.sim.func(N=10000, p=4000, mu=rep(0,4000), r=0)
}
sim_cor0_n10000_p4000





gao_effective_number_cor0 = cbind("effective number cor0 n10 p10"      = (mean(sim_cor0_n10_p10)),
                                  "effective number cor0 n100 p10"     = (mean(sim_cor0_n100_p10)),
                                  "effective number cor0 n1000 p10"    = (mean(sim_cor0_n1000_p10)),
                                  "effective number cor0 n10000 p10"   = (mean(sim_cor0_n10000_p10)),
                                  "effective number cor0 n10 p100"     = (mean(sim_cor0_n10_p100)),
                                  "effective number cor0 n100 p100"    = (mean(sim_cor0_n100_p100)),
                                  "effective number cor0 n1000 p100"   = (mean(sim_cor0_n1000_p100)),
                                  "effective number cor0 n10000 p100"  = (mean(sim_cor0_n10000_p100)),
                                  "effective number cor0 n10 p1000"    = (mean(sim_cor0_n10_p1000)),
                                  "effective number cor0 n100 p1000"   = (mean(sim_cor0_n100_p1000)),
                                  "effective number cor0 n1000 p1000"  = (mean(sim_cor0_n1000_p1000)),
                                  "effective number cor0 n10000 p1000" = (mean(sim_cor0_n10000_p1000)),
                                  "effective number cor0 n10 p4000"    = (mean(sim_cor0_n10_p4000)),
                                  "effective number cor0 n100 p4000"   = (mean(sim_cor0_n100_p4000)),
                                  "effective number cor0 n1000 p4000"  = (mean(sim_cor0_n1000_p4000)),
                                  "effective number cor0 n10000 p4000" = (mean(sim_cor0_n10000_p4000)))



gao_effective_number_cor0



# ---------------------------------------------------------------------------------------------------------------------------------



# start with 100% correlated matrix and simulate my data with number of variables = 100,
# changing number of predictors and see what values I get from my Gao's method

# N = 100 P = 10 r = 1
sim_cor1_n100_p10 = numeric(10)
for (i in 1:10){
  sim_cor1_n100_p10[i] = gao.sim.func(N=100, p=10, mu=rep(0,10), r=1)
}
sim_cor1_n100_p10

# N = 100 P = 100 r = 1
sim_cor1_n100_p100 = numeric(10)
for (i in 1:10){
  sim_cor1_n100_p100[i] = gao.sim.func(N=100, p=100, mu=rep(0,100), r=1)
}
sim_cor1_n100_p100

# N = 100 P = 4000 r = 1
sim_cor1_n100_p4000 = numeric(10)
for (i in 1:10){
  sim_cor1_n100_p4000[i] = gao.sim.func(N=100, p=4000, mu=rep(0,4000), r=1)
}
sim_cor1_n100_p4000


gao_effective_number_cor1 = cbind("effective number cor1 n100 p10" = (mean(sim_cor1_n100_p10)),
                                  "effective number cor1 n100 p100" = (mean(sim_cor1_n100_p100)),
                                  "effective number cor1 n100 p4000" = (mean(sim_cor1_n100_p4000)))
gao_effective_number_cor1
# looks great. It seems like the function did what I wanted it to do and got 1 effective number of e values for
# gao's method when correlation is 100%


###---------------------------------------------------------------------------------------------------------------------------------

### now let me do this with 80% correlated matix

# number of observations = 100 changing number of vairables and cor = .8
# N = 100 P = 10 r = .8
sim_cor.8_n100_p10 = numeric(10)
for (i in 1:10){
  sim_cor.8_n100_p10[i] = gao.sim.func(N=100, p=10, mu=rep(0,10), r=.8)
}
sim_cor.8_n100_p10

# N = 100 P = 20 r = .8
sim_cor.8_n100_p20 = numeric(10)
for (i in 1:10){
  sim_cor.8_n100_p20[i] = gao.sim.func(N=100, p=20, mu=rep(0,20), r=.8)
}
sim_cor.8_n100_p20

# N = 100 P = 100 r = .8
sim_cor.8_n100_p100 = numeric(10)
for (i in 1:10){
  sim_cor.8_n100_p100[i] = gao.sim.func(N=100, p=100, mu=rep(0,100), r=.8)
}
sim_cor.8_n100_p100

# N = 100 P = 200 r = .8
sim_cor.8_n100_p200 = numeric(10)
for (i in 1:10){
  sim_cor.8_n100_p200[i] = gao.sim.func(N=100, p=200, mu=rep(0,200), r=.8)
}
sim_cor.8_n100_p200

# N = 100 P = 1000 r = .8
sim_cor.8_n100_p1000 = numeric(10)
for (i in 1:10){
  sim_cor.8_n100_p1000[i] = gao.sim.func(N=100, p=1000, mu=rep(0,1000), r=.8)
}
sim_cor.8_n100_p1000

# N = 100 P = 4000 r = .8
sim_cor.8_n100_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.8_n100_p4000[i] = gao.sim.func(N=100, p=4000, mu=rep(0,4000), r=.8)
}
sim_cor.8_n100_p4000

gao_effective_number_cor.8_different_p = cbind("effective number cor.8 n100 p10" = (mean(sim_cor.8_n100_p10)),
                                               "effective number cor.8 n100 p20" = (mean(sim_cor.8_n100_p20)),
                                               "effective number cor.8 n100 p100" = (mean(sim_cor.8_n100_p100)),
                                               "effective number cor.8 n100 p200" = (mean(sim_cor.8_n100_p200)),
                                               "effective number cor.8 n100 p1000" = (mean(sim_cor.8_n100_p1000)),
                                               "effective number cor.8 n100 p4000" = (mean(sim_cor.8_n100_p4000)))
gao_effective_number_cor.8_different_p


# number of variables = 100 changing number of observations and cor = .8

# N = 10 P = 100 r = .8
sim_cor.8_n10_p100 = numeric(10)
for (i in 1:10){
  sim_cor.8_n10_p100[i] = gao.sim.func(N=10, p=100, mu=rep(0,100), r=.8)
}
sim_cor.8_n10_p100

# N = 20 P = 100 r = .8
sim_cor.8_n20_p100 = numeric(10)
for (i in 1:10){
  sim_cor.8_n20_p100[i] = gao.sim.func(N=20, p=100, mu=rep(0,100), r=.8)
}
sim_cor.8_n20_p100

# N = 100 P = 100 r = .8
sim_cor.8_n100_p100 = numeric(10)
for (i in 1:10){
  sim_cor.8_n100_p100[i] = gao.sim.func(N=100, p=100, mu=rep(0,100), r=.8)
}
sim_cor.8_n100_p100

# N = 200 P = 100 r = .8
sim_cor.8_n200_p100 = numeric(10)
for (i in 1:10){
  sim_cor.8_n200_p100[i] = gao.sim.func(N=200, p=100, mu=rep(0,100), r=.8)
}
sim_cor.8_n200_p100

# N = 200 P = 100 r = .8
sim_cor.8_n1000_p100 = numeric(10)
for (i in 1:10){
  sim_cor.8_n1000_p100[i] = gao.sim.func(N=1000, p=100, mu=rep(0,100), r=.8)
}
sim_cor.8_n1000_p100

gao_effective_number_cor.8_different_n = cbind("effective number cor.8 n10 p100" = (mean(sim_cor.8_n10_p100)),
                                               "effective number cor.8 n20 p100" = (mean(sim_cor.8_n20_p100)),
                                               "effective number cor.8 n100 p100" = (mean(sim_cor.8_n100_p100)),
                                               "effective number cor.8 n200 p100" = (mean(sim_cor.8_n200_p100)),
                                               "effective number cor.8 n1000 p100" = (mean(sim_cor.8_n1000_p100)))
gao_effective_number_cor.8_different_n


# number of variables = 4000 changing number of observations and cor = .8

# N = 10 P = 4000 r = .8
sim_cor.8_n10_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.8_n10_p4000[i] = gao.sim.func(N=10, p=4000, mu=rep(0,4000), r=.8)
}
sim_cor.8_n10_p4000

# N = 20 P = 4000 r = .8
sim_cor.8_n20_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.8_n20_p4000[i] = gao.sim.func(N=20, p=4000, mu=rep(0,4000), r=.8)
}
sim_cor.8_n20_p4000

# N = 100 P = 4000 r = .8
sim_cor.8_n100_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.8_n100_p4000[i] = gao.sim.func(N=100, p=4000, mu=rep(0,4000), r=.8)
}
sim_cor.8_n100_p4000

# N = 200 P = 4000 r = .8
sim_cor.8_n200_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.8_n200_p4000[i] = gao.sim.func(N=200, p=4000, mu=rep(0,4000), r=.8)
}
sim_cor.8_n200_p4000

# N = 200 P = 4000 r = .8
sim_cor.8_n1000_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.8_n1000_p4000[i] = gao.sim.func(N=1000, p=4000, mu=rep(0,4000), r=.8)
}
sim_cor.8_n1000_p4000

gao_effective_number_cor.8_different_n_p4000 = cbind("effective number cor.8 n10 p4000" = (mean(sim_cor.8_n10_p4000)),
                                                     "effective number cor.8 n20 p4000" = (mean(sim_cor.8_n20_p4000)),
                                                     "effective number cor.8 n100 p4000" = (mean(sim_cor.8_n100_p4000)),
                                                     "effective number cor.8 n200 p4000" = (mean(sim_cor.8_n200_p4000)),
                                                     "effective number cor.8 n1000 p4000" = (mean(sim_cor.8_n1000_p4000)))
gao_effective_number_cor.8_different_n_p4000

###---------------------------------------------------------------------------------------------------------------------------------
### now let me do this with 50% correlated matix

# number of observations = 100 changing number of vairables and cor = .5
# N = 100 P = 10 r = .5
sim_cor.5_n100_p10 = numeric(10)
for (i in 1:10){
  sim_cor.5_n100_p10[i] = gao.sim.func(N=100, p=10, mu=rep(0,10), r=.5)
}
sim_cor.5_n100_p10

# N = 100 P = 20 r = .5
sim_cor.5_n100_p20 = numeric(10)
for (i in 1:10){
  sim_cor.5_n100_p20[i] = gao.sim.func(N=100, p=20, mu=rep(0,20), r=.5)
}
sim_cor.5_n100_p20

# N = 100 P = 100 r = .5
sim_cor.5_n100_p100 = numeric(10)
for (i in 1:10){
  sim_cor.5_n100_p100[i] = gao.sim.func(N=100, p=100, mu=rep(0,100), r=.5)
}
sim_cor.5_n100_p100

# N = 100 P = 200 r = .5
sim_cor.5_n100_p200 = numeric(10)
for (i in 1:10){
  sim_cor.5_n100_p200[i] = gao.sim.func(N=100, p=200, mu=rep(0,200), r=.5)
}
sim_cor.5_n100_p200

# N = 100 P = 1000 r = .5
sim_cor.5_n100_p1000 = numeric(10)
for (i in 1:10){
  sim_cor.5_n100_p1000[i] = gao.sim.func(N=100, p=1000, mu=rep(0,1000), r=.5)
}
sim_cor.5_n100_p1000

# N = 100 P = 4000 r = .5
sim_cor.5_n100_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.5_n100_p4000[i] = gao.sim.func(N=100, p=4000, mu=rep(0,4000), r=.5)
}
sim_cor.5_n100_p4000

gao_effective_number_cor.5_different_p = cbind("effective number cor.5 n100 p10" = (mean(sim_cor.5_n100_p10)),
                                               "effective number cor.5 n100 p20" = (mean(sim_cor.5_n100_p20)),
                                               "effective number cor.5 n100 p100" = (mean(sim_cor.5_n100_p100)),
                                               "effective number cor.5 n100 p200" = (mean(sim_cor.5_n100_p200)),
                                               "effective number cor.5 n100 p1000" = (mean(sim_cor.5_n100_p1000)),
                                               "effective number cor.5 n100 p4000" = (mean(sim_cor.5_n100_p4000)))
gao_effective_number_cor.5_different_p

# number of variables = 100 changing number of observations and cor = .5

# N = 10 P = 100 r = .5
sim_cor.5_n10_p100 = numeric(10)
for (i in 1:10){
  sim_cor.5_n10_p100[i] = gao.sim.func(N=10, p=100, mu=rep(0,100), r=.5)
}
sim_cor.5_n10_p100

# N = 20 P = 100 r = .5
sim_cor.5_n20_p100 = numeric(10)
for (i in 1:10){
  sim_cor.5_n20_p100[i] = gao.sim.func(N=20, p=100, mu=rep(0,100), r=.5)
}
sim_cor.5_n20_p100

# N = 100 P = 100 r = .5
sim_cor.5_n100_p100 = numeric(10)
for (i in 1:10){
  sim_cor.5_n100_p100[i] = gao.sim.func(N=100, p=100, mu=rep(0,100), r=.5)
}
sim_cor.5_n100_p100

# N = 200 P = 100 r = .5
sim_cor.5_n200_p100 = numeric(10)
for (i in 1:10){
  sim_cor.5_n200_p100[i] = gao.sim.func(N=200, p=100, mu=rep(0,100), r=.5)
}
sim_cor.5_n200_p100

# N = 200 P = 100 r = .5
sim_cor.5_n1000_p100 = numeric(10)
for (i in 1:10){
  sim_cor.5_n1000_p100[i] = gao.sim.func(N=1000, p=100, mu=rep(0,100), r=.5)
}
sim_cor.5_n1000_p100

gao_effective_number_cor.5_different_n = cbind("effective number cor.5 n10 p100" = (mean(sim_cor.5_n10_p100)),
                                               "effective number cor.5 n20 p100" = (mean(sim_cor.5_n20_p100)),
                                               "effective number cor.5 n100 p100" = (mean(sim_cor.5_n100_p100)),
                                               "effective number cor.5 n200 p100" = (mean(sim_cor.5_n200_p100)),
                                               "effective number cor.5 n1000 p100" = (mean(sim_cor.5_n1000_p100)))
gao_effective_number_cor.5_different_n


# number of variables = 4000 changing number of observations and cor = .5

# N = 10 P = 4000 r = .5
sim_cor.5_n10_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.5_n10_p4000[i] = gao.sim.func(N=10, p=4000, mu=rep(0,4000), r=.5)
}
sim_cor.5_n10_p4000

# N = 20 P = 4000 r = .5
sim_cor.5_n20_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.5_n20_p4000[i] = gao.sim.func(N=20, p=4000, mu=rep(0,4000), r=.5)
}
sim_cor.5_n20_p4000

# N = 100 P = 4000 r = .5
sim_cor.5_n100_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.5_n100_p4000[i] = gao.sim.func(N=100, p=4000, mu=rep(0,4000), r=.5)
}
sim_cor.5_n100_p4000

# N = 200 P = 4000 r = .5
sim_cor.5_n200_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.5_n200_p4000[i] = gao.sim.func(N=200, p=4000, mu=rep(0,4000), r=.5)
}
sim_cor.5_n200_p4000

# N = 200 P = 4000 r = .5
sim_cor.5_n1000_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.5_n1000_p4000[i] = gao.sim.func(N=1000, p=4000, mu=rep(0,4000), r=.5)
}
sim_cor.5_n1000_p4000

gao_effective_number_cor.5_different_n_p4000 = cbind("effective number cor.5 n10 p4000" = (mean(sim_cor.5_n10_p4000)),
                                                     "effective number cor.5 n20 p4000" = (mean(sim_cor.5_n20_p4000)),
                                                     "effective number cor.5 n100 p4000" = (mean(sim_cor.5_n100_p4000)),
                                                     "effective number cor.5 n200 p4000" = (mean(sim_cor.5_n200_p4000)),
                                                     "effective number cor.5 n1000 p4000" = (mean(sim_cor.5_n1000_p4000)))
gao_effective_number_cor.5_different_n_p4000

###---------------------------------------------------------------------------------------------------------------------------------
### now let me do this with 20% correlated matix

# number of observations = 100 changing number of vairables and cor = .2
# N = 100 P = 10 r = .2
sim_cor.2_n100_p10 = numeric(10)
for (i in 1:10){
  sim_cor.2_n100_p10[i] = gao.sim.func(N=100, p=10, mu=rep(0,10), r=.2)
}
sim_cor.2_n100_p10

# N = 100 P = 20 r = .2
sim_cor.2_n100_p20 = numeric(10)
for (i in 1:10){
  sim_cor.2_n100_p20[i] = gao.sim.func(N=100, p=20, mu=rep(0,20), r=.2)
}
sim_cor.2_n100_p20

# N = 100 P = 100 r = .2
sim_cor.2_n100_p100 = numeric(10)
for (i in 1:10){
  sim_cor.2_n100_p100[i] = gao.sim.func(N=100, p=100, mu=rep(0,100), r=.2)
}
sim_cor.2_n100_p100

# N = 100 P = 200 r = .2
sim_cor.2_n100_p200 = numeric(10)
for (i in 1:10){
  sim_cor.2_n100_p200[i] = gao.sim.func(N=100, p=200, mu=rep(0,200), r=.2)
}
sim_cor.2_n100_p200

# N = 100 P = 1000 r = .2
sim_cor.2_n100_p1000 = numeric(10)
for (i in 1:10){
  sim_cor.2_n100_p1000[i] = gao.sim.func(N=100, p=1000, mu=rep(0,1000), r=.2)
}
sim_cor.2_n100_p1000

# N = 100 P = 4000 r = .2
sim_cor.2_n100_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.2_n100_p4000[i] = gao.sim.func(N=100, p=4000, mu=rep(0,4000), r=.2)
}
sim_cor.2_n100_p4000

gao_effective_number_cor.2_different_p = cbind("effective number cor.2 n100 p10" = (mean(sim_cor.2_n100_p10)),
                                               "effective number cor.2 n100 p20" = (mean(sim_cor.2_n100_p20)),
                                               "effective number cor.2 n100 p100" = (mean(sim_cor.2_n100_p100)),
                                               "effective number cor.2 n100 p200" = (mean(sim_cor.2_n100_p200)),
                                               "effective number cor.2 n100 p1000" = (mean(sim_cor.2_n100_p1000)),
                                               "effective number cor.2 n100 p4000" = (mean(sim_cor.2_n100_p4000)))
gao_effective_number_cor.2_different_p

# number of variables = 100 changing number of observations and cor = .2

# N = 10 P = 100 r = .2
sim_cor.2_n10_p100 = numeric(10)
for (i in 1:10){
  sim_cor.2_n10_p100[i] = gao.sim.func(N=10, p=100, mu=rep(0,100), r=.2)
}
sim_cor.2_n10_p100

# N = 20 P = 100 r = .2
sim_cor.2_n20_p100 = numeric(10)
for (i in 1:10){
  sim_cor.2_n20_p100[i] = gao.sim.func(N=20, p=100, mu=rep(0,100), r=.2)
}
sim_cor.2_n20_p100

# N = 100 P = 100 r = .2
sim_cor.2_n100_p100 = numeric(10)
for (i in 1:10){
  sim_cor.2_n100_p100[i] = gao.sim.func(N=100, p=100, mu=rep(0,100), r=.2)
}
sim_cor.2_n100_p100

# N = 200 P = 100 r = .2
sim_cor.2_n200_p100 = numeric(10)
for (i in 1:10){
  sim_cor.2_n200_p100[i] = gao.sim.func(N=200, p=100, mu=rep(0,100), r=.2)
}
sim_cor.2_n200_p100

# N = 200 P = 100 r = .2
sim_cor.2_n1000_p100 = numeric(10)
for (i in 1:10){
  sim_cor.2_n1000_p100[i] = gao.sim.func(N=1000, p=100, mu=rep(0,100), r=.2)
}
sim_cor.2_n1000_p100

gao_effective_number_cor.2_different_n = cbind("effective number cor.2 n10 p100" = (mean(sim_cor.2_n10_p100)),
                                               "effective number cor.2 n20 p100" = (mean(sim_cor.2_n20_p100)),
                                               "effective number cor.2 n100 p100" = (mean(sim_cor.2_n100_p100)),
                                               "effective number cor.2 n200 p100" = (mean(sim_cor.2_n200_p100)),
                                               "effective number cor.2 n1000 p100" = (mean(sim_cor.2_n1000_p100)))
gao_effective_number_cor.2_different_n


# number of variables = 4000 changing number of observations and cor = .2

# N = 10 P = 4000 r = .2
sim_cor.2_n10_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.2_n10_p4000[i] = gao.sim.func(N=10, p=4000, mu=rep(0,4000), r=.2)
}
sim_cor.2_n10_p4000

# N = 20 P = 4000 r = .2
sim_cor.2_n20_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.2_n20_p4000[i] = gao.sim.func(N=20, p=4000, mu=rep(0,4000), r=.2)
}
sim_cor.2_n20_p4000

# N = 100 P = 4000 r = .2
sim_cor.2_n100_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.2_n100_p4000[i] = gao.sim.func(N=100, p=4000, mu=rep(0,4000), r=.2)
}
sim_cor.2_n100_p4000

# N = 200 P = 4000 r = .2
sim_cor.2_n200_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.2_n200_p4000[i] = gao.sim.func(N=200, p=4000, mu=rep(0,4000), r=.2)
}
sim_cor.2_n200_p4000

# N = 200 P = 4000 r = .2
sim_cor.2_n1000_p4000 = numeric(10)
for (i in 1:10){
  sim_cor.2_n1000_p4000[i] = gao.sim.func(N=1000, p=4000, mu=rep(0,4000), r=.2)
}
sim_cor.2_n1000_p4000

gao_effective_number_cor.2_different_n_p4000 = cbind("effective number cor.2 n10 p4000" = (mean(sim_cor.2_n10_p4000)),
                                                     "effective number cor.2 n20 p4000" = (mean(sim_cor.2_n20_p4000)),
                                                     "effective number cor.2 n100 p4000" = (mean(sim_cor.2_n100_p4000)),
                                                     "effective number cor.2 n200 p4000" = (mean(sim_cor.2_n200_p4000)),
                                                     "effective number cor.2 n1000 p4000" = (mean(sim_cor.2_n1000_p4000)))
gao_effective_number_cor.2_different_n_p4000

#-----------------------------------------------------------------------------------------------------------------------------------------------------


gao_effective_simulated_data = cbind(gao_effective_number_cor0,
                                     gao_effective_number_cor1, 
                                     gao_effective_number_cor.8_different_p, 
                                     gao_effective_number_cor.8_different_n,
                                     gao_effective_number_cor.8_different_n_p4000,
                                     gao_effective_number_cor.5_different_p,
                                     gao_effective_number_cor.5_different_n,
                                     gao_effective_number_cor.5_different_n_p4000,
                                     gao_effective_number_cor.2_different_p,
                                     gao_effective_number_cor.2_different_n,
                                     gao_effective_number_cor.2_different_n_p4000)

gao_effective_simulated_data

setwd("C:/Users/sunga/Desktop/stats thesis project")
write.table(gao_effective_simulated_data, file = "gao simulated data.txt")

#------------------------------------------------------------------------------------------------------------------------------------------------------




########################################################################################
### time to create blocks of data and see how well gaos method performs with them ###
#######################################################################################

set.seed(1)
library(magic)
library(MASS)


# ---------------------------------------------------------------------------------------------------------------

# this function is for 2 blocks 

gao.2_bloc.sim.func<-function(N, p, b, mu, r){
  
  # creating a matrix of my correlations
  x <- matrix(nrow=p, ncol=p, r)
  # creating a block matrix
  x2<-adiag(x,x)
  
  
  dat.r1 <- scale(mvrnorm(N, mu, x2))
  # now I am creating a correlation matrix for the simulated block data
  cor_dat_1 = cor(dat.r1)
  
  
  # now lets calculate the eigenvalues using eigen function with the covariance matrix from above
  eigenvalues_cov = eigen(cor_dat_1, symmetric=T, only.values=T)
  summary(eigenvalues_cov$values)
  
  ##  Gao's Method  ##
  sum_eigen = sum(eigenvalues_cov$values)
  
  e2 = cbind(eigenvalues_cov$values, prop=eigenvalues_cov$values/sum_eigen)
  
  prop_total = e2[1,2]
  
  if(prop_total >= .995){
    prop_total = prop_total
    
  }else{
    
    for (i in 2:length(e2[,2])) {prop_total = c(prop_total, prop_total[i-1]+e2[i,2])
    
    if(prop_total[i] >= .995) {break}}
    
  }
  
  
  M_gao_cov = length(prop_total)
  
  return(M_gao_cov)
  
}


#-------

# mu=rep(0,p*b) I need to input the p*b because the function does not seem to take in 

# 2 blocks p = 10 r = 1 n = 100
block_gao_p10_b2_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p10_b2_r1_n100[i] = gao.2_bloc.sim.func(N=100, p=10, b=2, mu=rep(0,20), r=1)
}
(block_gao_p10_b2_r1_n100)


# 2 blocks p = 100 r = 1 n = 100
block_gao_p100_b2_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p100_b2_r1_n100[i] = gao.2_bloc.sim.func(N=100, p=100, b=2, mu=rep(0,200), r=1)
}
(block_gao_p100_b2_r1_n100)


# 2 blocks p = 1000 r = 1 n = 100
block_gao_p1000_b2_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p1000_b2_r1_n100[i] = gao.2_bloc.sim.func(N=100, p=1000, b=2, mu=rep(0,2000), r=1)
}
(block_gao_p1000_b2_r1_n100)

# 2 blocks p = 1000 r = 1 n = 1000
block_gao_p1000_b2_r1_n1000 = numeric(10)
for (i in 1:10){
  block_gao_p1000_b2_r1_n1000[i] = gao.2_bloc.sim.func(N=1000, p=1000, b=2, mu=rep(0,2000), r=1)
}
(block_gao_p1000_b2_r1_n1000)



gao_effective_number_2.blocks = cbind("effective number 2 blocks p = 10 r = 1 n = 100"    = ((block_gao_p10_b2_r1_n100)),
                                      "effective number 2 blocks p = 100 r = 1 n = 100"   = ((block_gao_p100_b2_r1_n100)),
                                      "effective number 2 blocks p = 1000 r = 1 n = 100"  = ((block_gao_p1000_b2_r1_n100)),
                                      "effective number 2 blocks p = 1000 r = 1 n = 1000" = ((block_gao_p1000_b2_r1_n1000)))

gao_effective_number_2.blocks

# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------

# this function is for 3 blocks 

gao.3_bloc.sim.func<-function(N, p, b, mu, r){
  
  # creating a matrix of my correlations
  x <- matrix(nrow=p, ncol=p, r)
  # creating a block matrix
  x2<-adiag(x,x,x)
  
  
  dat.r1 <- scale(mvrnorm(N, mu, x2))
  # now I am creating a correlation matrix for the simulated block data
  cor_dat_1 = cor(dat.r1)
  
  
  # now lets calculate the eigenvalues using eigen function with the covariance matrix from above
  eigenvalues_cov = eigen(cor_dat_1, symmetric=T, only.values=T)
  summary(eigenvalues_cov$values)
  
  ##  Gao's Method  ##
  sum_eigen = sum(eigenvalues_cov$values)
  
  e2 = cbind(eigenvalues_cov$values, prop=eigenvalues_cov$values/sum_eigen)
  
  prop_total = e2[1,2]
  
  if(prop_total >= .995){
    prop_total = prop_total
    
  }else{
    
    for (i in 2:length(e2[,2])) {prop_total = c(prop_total, prop_total[i-1]+e2[i,2])
    
    if(prop_total[i] >= .995) {break}}
    
  }
  
  
  M_gao_cov = length(prop_total)
  
  return(M_gao_cov)
  
}


#-------

# mu=rep(0,p*b) I need to input the p*b because the function does not seem to take in 

# 3 blocks p = 10 r = 1 n = 100
block_gao_p10_b3_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p10_b3_r1_n100[i] = gao.3_bloc.sim.func(N=100, p=10, b=3, mu=rep(0,30), r=1)
}
(block_gao_p10_b3_r1_n100)


# 3 blocks p = 100 r = 1 n = 100
block_gao_p100_b3_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p100_b3_r1_n100[i] = gao.3_bloc.sim.func(N=100, p=100, b=3, mu=rep(0,300), r=1)
}
(block_gao_p100_b3_r1_n100)


# 3 blocks p = 1000 r = 1 n = 100
block_gao_p1000_b3_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p1000_b3_r1_n100[i] = gao.3_bloc.sim.func(N=100, p=1000, b=3, mu=rep(0,3000), r=1)
}
(block_gao_p1000_b3_r1_n100)

# 3 blocks p = 1000 r = 1 n = 1000
block_gao_p1000_b3_r1_n1000 = numeric(10)
for (i in 1:10){
  block_gao_p1000_b3_r1_n1000[i] = gao.3_bloc.sim.func(N=1000, p=1000, b=3, mu=rep(0,3000), r=1)
}
(block_gao_p1000_b3_r1_n1000)



gao_effective_number_3.blocks = cbind("effective number 3 blocks p = 10 r = 1 n = 100"   = ((block_gao_p10_b3_r1_n100)),
                                      "effective number 3 blocks p = 100 r = 1 n = 100"  = ((block_gao_p100_b3_r1_n100)),
                                      "effective number 3 blocks p = 1000 r = 1 n = 100" = ((block_gao_p1000_b3_r1_n100)),
                                      "effective number 3 blocks p = 1000 r = 1 n = 1000" = ((block_gao_p1000_b3_r1_n1000)))

gao_effective_number_3.blocks

# ---------------------------------------------------------------------------------------------------------------





# ---------------------------------------------------------------------------------------------------------------

# this function is for 5 blocks 

gao.5_bloc.sim.func<-function(N, p, b, mu, r){
  
  # creating a matrix of my correlations
  x <- matrix(nrow=p, ncol=p, r)
  # creating a block matrix
  x2<-adiag(x,x,x,x,x)
  
  
  dat.r1 <- scale(mvrnorm(N, mu, x2))
  # now I am creating a correlation matrix for the simulated block data
  cor_dat_1 = cor(dat.r1)
  
  
  # now lets calculate the eigenvalues using eigen function with the covariance matrix from above
  eigenvalues_cov = eigen(cor_dat_1, symmetric=T, only.values=T)
  summary(eigenvalues_cov$values)
  
  ##  Gao's Method  ##
  sum_eigen = sum(eigenvalues_cov$values)
  
  e2 = cbind(eigenvalues_cov$values, prop=eigenvalues_cov$values/sum_eigen)
  
  prop_total = e2[1,2]
  
  if(prop_total >= .995){
    prop_total = prop_total
    
  }else{
    
    for (i in 2:length(e2[,2])) {prop_total = c(prop_total, prop_total[i-1]+e2[i,2])
    
    if(prop_total[i] >= .995) {break}}
    
  }
  
  
  M_gao_cov = length(prop_total)
  
  return(M_gao_cov)
  
}


#-------

# mu=rep(0,p*b) I need to input the p*b because the function does not seem to take in 

# 5 blocks p = 10 r = 1 n = 100
block_gao_p10_b5_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p10_b5_r1_n100[i] = gao.5_bloc.sim.func(N=100, p=10, b=5, mu=rep(0,50), r=1)
}
(block_gao_p10_b5_r1_n100)


# 5 blocks p = 100 r = 1 n = 100
block_gao_p100_b5_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p100_b5_r1_n100[i] = gao.5_bloc.sim.func(N=100, p=100, b=5, mu=rep(0,500), r=1)
}
(block_gao_p100_b5_r1_n100)


# 5 blocks p = 1000 r = 1 n = 100
block_gao_p1000_b5_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p1000_b5_r1_n100[i] = gao.5_bloc.sim.func(N=100, p=1000, b=5, mu=rep(0,5000), r=1)
}
(block_gao_p1000_b5_r1_n100)


# 5 blocks p = 1000 r = 1 n = 1000
block_gao_p1000_b5_r1_n1000 = numeric(10)
for (i in 1:10){
  block_gao_p1000_b5_r1_n1000[i] = gao.5_bloc.sim.func(N=1000, p=1000, b=5, mu=rep(0,5000), r=1)
}
(block_gao_p1000_b5_r1_n1000)


gao_effective_number_5.blocks = cbind("effective number 5 blocks p = 10 r = 1 n = 100"   = ((block_gao_p10_b5_r1_n100)),
                                      "effective number 5 blocks p = 100 r = 1 n = 100"  = ((block_gao_p100_b5_r1_n100)),
                                      "effective number 5 blocks p = 1000 r = 1 n = 100" = ((block_gao_p1000_b5_r1_n100)),
                                      "effective number 5 blocks p = 1000 r = 1 n = 1000" = ((block_gao_p1000_b5_r1_n1000)))

gao_effective_number_5.blocks

# ---------------------------------------------------------------------------------------------------------------





# ---------------------------------------------------------------------------------------------------------------

# this function is for 10 blocks 

gao.10_bloc.sim.func<-function(N, p, b, mu, r){
  
  # creating a matrix of my correlations
  x <- matrix(nrow=p, ncol=p, r)
  # creating a block matrix
  x2<-adiag(x,x,x,x,x,x,x,x,x,x)
  
  
  dat.r1 <- scale(mvrnorm(N, mu, x2))
  # now I am creating a correlation matrix for the simulated block data
  cor_dat_1 = cor(dat.r1)
  
  
  # now lets calculate the eigenvalues using eigen function with the covariance matrix from above
  eigenvalues_cov = eigen(cor_dat_1, symmetric=T, only.values=T)
  summary(eigenvalues_cov$values)
  
  ##  Gao's Method  ##
  sum_eigen = sum(eigenvalues_cov$values)
  
  e2 = cbind(eigenvalues_cov$values, prop=eigenvalues_cov$values/sum_eigen)
  
  prop_total = e2[1,2]
  
  if(prop_total >= .995){
    prop_total = prop_total
    
  }else{
    
    for (i in 2:length(e2[,2])) {prop_total = c(prop_total, prop_total[i-1]+e2[i,2])
    
    if(prop_total[i] >= .995) {break}}
    
  }
  
  
  M_gao_cov = length(prop_total)
  
  return(M_gao_cov)
  
}


#-------

# mu=rep(0,p*b) I need to input the p*b because the function does not seem to take in 

# 10 blocks p = 10 r = 1 n = 100
block_gao_p10_b10_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p10_b10_r1_n100[i] = gao.10_bloc.sim.func(N=100, p=10, b=10, mu=rep(0,100), r=1)
}
(block_gao_p10_b10_r1_n100)


# 10 blocks p = 100 r = 1 n = 100
block_gao_p100_b10_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p100_b10_r1_n100[i] = gao.10_bloc.sim.func(N=100, p=100, b=10, mu=rep(0,1000), r=1)
}
(block_gao_p100_b10_r1_n100)


# 10 blocks p = 1000 r = 1 n = 100
block_gao_p1000_b10_r1_n100 = numeric(10)
for (i in 1:10){
  block_gao_p1000_b10_r1_n100[i] = gao.10_bloc.sim.func(N=100, p=1000, b=10, mu=rep(0,10000), r=1)
}
(block_gao_p1000_b10_r1_n100)


# 10 blocks p = 1000 r = 1 n = 1000
block_gao_p1000_b10_r1_n1000 = numeric(10)
for (i in 1:10){
  block_gao_p1000_b10_r1_n1000[i] = gao.10_bloc.sim.func(N=1000, p=1000, b=10, mu=rep(0,10000), r=1)
}
(block_gao_p1000_b10_r1_n1000)


gao_effective_number_10.blocks = cbind("effective number 10 blocks p = 10 r = 1 n = 100"   = ((block_gao_p10_b10_r1_n100)),
                                       "effective number 10 blocks p = 100 r = 1 n = 100"  = ((block_gao_p100_b10_r1_n100)),
                                       "effective number 10 blocks p = 1000 r = 1 n = 100" = ((block_gao_p1000_b10_r1_n100)),
                                       "effective number 10 blocks p = 1000 r = 1 n = 1000" = ((block_gao_p1000_b10_r1_n1000)))

gao_effective_number_10.blocks

# ---------------------------------------------------------------------------------------------------------------

gao_effective_simulated_data = cbind(gao_effective_number_2.blocks,
                                     gao_effective_number_3.blocks,
                                     gao_effective_number_5.blocks,
                                     gao_effective_number_10.blocks)

gao_effective_simulated_data

















##############################################################
### calculating the standard deviations for each simulation###
##############################################################
# standard deviation for  simulation with cor 0

###--------------------------------------------------------------------------------------------
# p = 10 changing n
sd_sim_cor0_n10_p10 = sd(sim_cor0_n10_p10)
sd_sim_cor0_n10_p10

sd_sim_cor0_n100_p10 = sd(sim_cor0_n100_p10)
sd_sim_cor0_n100_p10

sd_sim_cor0_n1000_p10 = sd(sim_cor0_n1000_p10)
sd_sim_cor0_n1000_p10

sd_sim_cor0_n10000_p10 = sd(sim_cor0_n10000_p10)
sd_sim_cor0_n10000_p10


# p = 100 changin n
sd_sim_cor0_n10_p100 = sd(sim_cor0_n10_p100)
sd_sim_cor0_n10_p100

sd_sim_cor0_n100_p100 = sd(sim_cor0_n100_p100)
sd_sim_cor0_n100_p100

sd_sim_cor0_n1000_p100 = sd(sim_cor0_n1000_p100)
sd_sim_cor0_n1000_p100

sd_sim_cor0_n10000_p100 = sd(sim_cor0_n10000_p100)
sd_sim_cor0_n10000_p100

# p =1000 changin n
sd_sim_cor0_n10_p1000 = sd(sim_cor0_n10_p1000)
sd_sim_cor0_n10_p1000

sd_sim_cor0_n100_p1000 = sd(sim_cor0_n100_p1000)
sd_sim_cor0_n100_p1000

sd_sim_cor0_n1000_p1000 = sd(sim_cor0_n1000_p1000)
sd_sim_cor0_n1000_p1000

sd_sim_cor0_n10000_p1000 = sd(sim_cor0_n10000_p1000)
sd_sim_cor0_n10000_p1000


# p = 4000 changin n
sd_sim_cor0_n10_p4000 = sd(sim_cor0_n10_p4000)
sd_sim_cor0_n10_p4000

sd_sim_cor0_n100_p4000 = sd(sim_cor0_n100_p4000)
sd_sim_cor0_n100_p4000

sd_sim_cor0_n1000_p4000 = sd(sim_cor0_n1000_p4000)
sd_sim_cor0_n1000_p4000

sd_sim_cor0_n10000_p4000 = sd(sim_cor0_n10000_p4000)
sd_sim_cor0_n10000_p4000


sd_gao_effective_number_cor0 = cbind("sd for effective number cor0 n10 p10"      = sd_sim_cor0_n10_p10,
                                     "sd for effective number cor0 n100 p10"     = sd_sim_cor0_n100_p10,
                                     "sd for effective number cor0 n1000 p10"    = sd_sim_cor0_n1000_p10,
                                     "sd for effective number cor0 n10000 p10"   = sd_sim_cor0_n10000_p10,
                                     "sd for effective number cor0 n10 p100"     = sd_sim_cor0_n10_p100,
                                     "sd for effective number cor0 n100 p100"    = sd_sim_cor0_n100_p100,
                                     "sd for effective number cor0 n1000 p100"   = sd_sim_cor0_n1000_p100,
                                     "sd for effective number cor0 n10000 p100"  = sd_sim_cor0_n10000_p100,
                                     "sd for effective number cor0 n10 p1000"    = sd_sim_cor0_n10_p1000,
                                     "sd for effective number cor0 n100 p1000"   = sd_sim_cor0_n100_p1000,
                                     "sd for effective number cor0 n1000 p1000"  = sd_sim_cor0_n1000_p1000,
                                     "sd for effective number cor0 n10000 p1000" = sd_sim_cor0_n10000_p1000,
                                     "sd for effective number cor0 n10 p4000"    = sd_sim_cor0_n10_p4000,
                                     "sd for effective number cor0 n100 p4000"   = sd_sim_cor0_n100_p4000,
                                     "sd for effective number cor0 n1000 p4000"  = sd_sim_cor0_n1000_p4000,
                                     "sd for effective number cor0 n10000 p4000" = sd_sim_cor0_n10000_p4000)



sd_gao_effective_number_cor0

###--------------------------------------------------------------------------------------------

# standard deviation for  simulation with cor = 1

sd_sim_cor1_n100_p10 = sd(sim_cor1_n100_p10)
sd_sim_cor1_n100_p10

sd_sim_cor1_n100_p100 = sd(sim_cor1_n100_p100)
sd_sim_cor1_n100_p100

sd_sim_cor1_n100_p4000 = sd(sim_cor1_n100_p4000)
sd_sim_cor1_n100_p4000


sd_gao_effective_number_cor1 = cbind("sd for effective number cor1 n100 p10" = sd_sim_cor1_n100_p10,
                                     "sd for effective number cor1 n100 p100" = sd_sim_cor1_n100_p100,
                                     "sd for effective number cor1 n100 p4000" = sd_sim_cor1_n100_p4000)


sd_gao_effective_number_cor1

###--------------------------------------------------------------------------------------------

# standard deviation for  simulation with cor = 0.8 changing p 
sd_sim_cor.8_n100_p10 = sd(sim_cor.8_n100_p10)
sd_sim_cor.8_n100_p10

sd_sim_cor.8_n100_p20 = sd(sim_cor.8_n100_p20)
sd_sim_cor.8_n100_p20

sd_sim_cor.8_n100_p100 = sd(sim_cor.8_n100_p100)
sd_sim_cor.8_n100_p100

sd_sim_cor.8_n100_p200 = sd(sim_cor.8_n100_p200)
sd_sim_cor.8_n100_p200

sd_sim_cor.8_n100_p1000 = sd(sim_cor.8_n100_p1000)
sd_sim_cor.8_n100_p1000

sd_sim_cor.8_n100_p4000 = sd(sim_cor.8_n100_p4000)
sd_sim_cor.8_n100_p4000

sd_gao_effective_number_cor.8_different_p = cbind("sd for effective number cor.8 n100 p10" = sd_sim_cor.8_n100_p10,
                                                  "sd for effective number cor.8 n100 p20" = sd_sim_cor.8_n100_p20,
                                                  "sd for effective number cor.8 n100 p100" = sd_sim_cor.8_n100_p100,
                                                  "sd for effective number cor.8 n100 p200" = sd_sim_cor.8_n100_p200,
                                                  "sd for effective number cor.8 n100 p1000" = sd_sim_cor.8_n100_p1000,
                                                  "sd for effective number cor.8 n100 p4000" = sd_sim_cor.8_n100_p4000)
sd_gao_effective_number_cor.8_different_p

# standard deviation for  simulation with cor = 0.8 changing n
sd_sim_cor.8_n10_p100 = sd(sim_cor.8_n10_p100)
sd_sim_cor.8_n10_p100

sd_sim_cor.8_n20_p100 = sd(sim_cor.8_n20_p100)
sd_sim_cor.8_n20_p100

sd_sim_cor.8_n100_p100 = sd(sim_cor.8_n100_p100)
sd_sim_cor.8_n100_p100

sd_sim_cor.8_n200_p100 = sd(sim_cor.8_n200_p100)
sd_sim_cor.8_n200_p100

sd_sim_cor.8_n1000_p100 = sd(sim_cor.8_n1000_p100)
sd_sim_cor.8_n1000_p100



sd_gao_effective_number_cor.8_different_n = cbind("sd for effective number cor.8 n10 p100" = sd_sim_cor.8_n10_p100,
                                                  "sd for effective number cor.8 n20 p100" = sd_sim_cor.8_n20_p100,
                                                  "sd for effective number cor.8 n100 p100" = sd_sim_cor.8_n100_p100,
                                                  "sd for effective number cor.8 n200 p100" = sd_sim_cor.8_n200_p100,
                                                  "sd for effective number cor.8 n1000 p100" = sd_sim_cor.8_n1000_p100)
sd_gao_effective_number_cor.8_different_n


# standard deviation for  simulation with cor = 0.8 changing n p = 4000
sd_sim_cor.8_n10_p4000 = sd(sim_cor.8_n10_p4000)
sd_sim_cor.8_n10_p4000

sd_sim_cor.8_n20_p4000 = sd(sim_cor.8_n20_p4000)
sd_sim_cor.8_n20_p4000

sd_sim_cor.8_n100_p4000 = sd(sim_cor.8_n100_p4000)
sd_sim_cor.8_n100_p4000

sd_sim_cor.8_n200_p4000 = sd(sim_cor.8_n200_p4000)
sd_sim_cor.8_n200_p4000

sd_sim_cor.8_n1000_p4000 = sd(sim_cor.8_n1000_p4000)
sd_sim_cor.8_n1000_p4000


sd_gao_effective_number_cor.8_different_n_p4000 = cbind("sd for effective number cor.8 n10 p4000" = sd_sim_cor.8_n10_p4000,
                                                        "sd for effective number cor.8 n20 p4000" = sd_sim_cor.8_n20_p4000,
                                                        "sd for effective number cor.8 n100 p4000" = sd_sim_cor.8_n100_p4000,
                                                        "sd for effective number cor.8 n200 p4000" = sd_sim_cor.8_n200_p4000,
                                                        "sd for effective number cor.8 n1000 p4000" = sd_sim_cor.8_n1000_p4000)
sd_gao_effective_number_cor.8_different_n_p4000

###--------------------------------------------------------------------------------------------

# standard deviation for  simulation with cor = 0.5 changing p 
sd_sim_cor.5_n100_p10 = sd(sim_cor.5_n100_p10)
sd_sim_cor.5_n100_p10

sd_sim_cor.5_n100_p20 = sd(sim_cor.5_n100_p20)
sd_sim_cor.5_n100_p20

sd_sim_cor.5_n100_p100 = sd(sim_cor.5_n100_p100)
sd_sim_cor.5_n100_p100

sd_sim_cor.5_n100_p200 = sd(sim_cor.5_n100_p200)
sd_sim_cor.5_n100_p200

sd_sim_cor.5_n100_p1000 = sd(sim_cor.5_n100_p1000)
sd_sim_cor.5_n100_p1000

sd_sim_cor.5_n100_p4000 = sd(sim_cor.5_n100_p4000)
sd_sim_cor.5_n100_p4000

sd_gao_effective_number_cor.5_different_p = cbind("sd for effective number cor.5 n100 p10" = sd_sim_cor.5_n100_p10,
                                                  "sd for effective number cor.5 n100 p20" = sd_sim_cor.5_n100_p20,
                                                  "sd for effective number cor.5 n100 p100" = sd_sim_cor.5_n100_p100,
                                                  "sd for effective number cor.5 n100 p200" = sd_sim_cor.5_n100_p200,
                                                  "sd for effective number cor.5 n100 p1000" = sd_sim_cor.5_n100_p1000,
                                                  "sd for effective number cor.5 n100 p4000" = sd_sim_cor.5_n100_p4000)
sd_gao_effective_number_cor.5_different_p

# standard deviation for  simulation with cor = 0.5 changing n
sd_sim_cor.5_n10_p100 = sd(sim_cor.5_n10_p100)
sd_sim_cor.5_n10_p100

sd_sim_cor.5_n20_p100 = sd(sim_cor.5_n20_p100)
sd_sim_cor.5_n20_p100

sd_sim_cor.5_n100_p100 = sd(sim_cor.5_n100_p100)
sd_sim_cor.5_n100_p100

sd_sim_cor.5_n200_p100 = sd(sim_cor.5_n200_p100)
sd_sim_cor.5_n200_p100

sd_sim_cor.5_n1000_p100 = sd(sim_cor.5_n1000_p100)
sd_sim_cor.5_n1000_p100



sd_gao_effective_number_cor.5_different_n = cbind("sd for effective number cor.5 n10 p100" = sd_sim_cor.5_n10_p100,
                                                  "sd for effective number cor.5 n20 p100" = sd_sim_cor.5_n20_p100,
                                                  "sd for effective number cor.5 n100 p100" = sd_sim_cor.5_n100_p100,
                                                  "sd for effective number cor.5 n200 p100" = sd_sim_cor.5_n200_p100,
                                                  "sd for effective number cor.5 n1000 p100" = sd_sim_cor.5_n1000_p100)
sd_gao_effective_number_cor.5_different_n


# standard deviation for  simulation with cor = 0.5 changing n p = 4000
sd_sim_cor.5_n10_p4000 = sd(sim_cor.5_n10_p4000)
sd_sim_cor.5_n10_p4000

sd_sim_cor.5_n20_p4000 = sd(sim_cor.5_n20_p4000)
sd_sim_cor.5_n20_p4000

sd_sim_cor.5_n100_p4000 = sd(sim_cor.5_n100_p4000)
sd_sim_cor.5_n100_p4000

sd_sim_cor.5_n200_p4000 = sd(sim_cor.5_n200_p4000)
sd_sim_cor.5_n200_p4000

sd_sim_cor.5_n1000_p4000 = sd(sim_cor.5_n1000_p4000)
sd_sim_cor.5_n1000_p4000


sd_gao_effective_number_cor.5_different_n_p4000 = cbind("sd for effective number cor.5 n10 p4000" = sd_sim_cor.5_n10_p4000,
                                                        "sd for effective number cor.5 n20 p4000" = sd_sim_cor.5_n20_p4000,
                                                        "sd for effective number cor.5 n100 p4000" = sd_sim_cor.5_n100_p4000,
                                                        "sd for effective number cor.5 n200 p4000" = sd_sim_cor.5_n200_p4000,
                                                        "sd for effective number cor.5 n1000 p4000" = sd_sim_cor.5_n1000_p4000)
sd_gao_effective_number_cor.5_different_n_p4000


# -----------------------------------------------------------------------------------------------------------

# standard deviation for  simulation with cor = 0.2 changing p 
sd_sim_cor.2_n100_p10 = sd(sim_cor.2_n100_p10)
sd_sim_cor.2_n100_p10

sd_sim_cor.2_n100_p20 = sd(sim_cor.2_n100_p20)
sd_sim_cor.2_n100_p20

sd_sim_cor.2_n100_p100 = sd(sim_cor.2_n100_p100)
sd_sim_cor.2_n100_p100

sd_sim_cor.2_n100_p200 = sd(sim_cor.2_n100_p200)
sd_sim_cor.2_n100_p200

sd_sim_cor.2_n100_p1000 = sd(sim_cor.2_n100_p1000)
sd_sim_cor.2_n100_p1000

sd_sim_cor.2_n100_p4000 = sd(sim_cor.2_n100_p4000)
sd_sim_cor.2_n100_p4000

sd_gao_effective_number_cor.2_different_p = cbind("sd for effective number cor.2 n100 p10" = sd_sim_cor.2_n100_p10,
                                                  "sd for effective number cor.2 n100 p20" = sd_sim_cor.2_n100_p20,
                                                  "sd for effective number cor.2 n100 p100" = sd_sim_cor.2_n100_p100,
                                                  "sd for effective number cor.2 n100 p200" = sd_sim_cor.2_n100_p200,
                                                  "sd for effective number cor.2 n100 p1000" = sd_sim_cor.2_n100_p1000,
                                                  "sd for effective number cor.2 n100 p4000" = sd_sim_cor.2_n100_p4000)
sd_gao_effective_number_cor.2_different_p

# standard deviation for  simulation with cor = 0.2 changing n
sd_sim_cor.2_n10_p100 = sd(sim_cor.2_n10_p100)
sd_sim_cor.2_n10_p100

sd_sim_cor.2_n20_p100 = sd(sim_cor.2_n20_p100)
sd_sim_cor.2_n20_p100

sd_sim_cor.2_n100_p100 = sd(sim_cor.2_n100_p100)
sd_sim_cor.2_n100_p100

sd_sim_cor.2_n200_p100 = sd(sim_cor.2_n200_p100)
sd_sim_cor.2_n200_p100

sd_sim_cor.2_n1000_p100 = sd(sim_cor.2_n1000_p100)
sd_sim_cor.2_n1000_p100



sd_gao_effective_number_cor.2_different_n = cbind("sd for effective number cor.2 n10 p100" = sd_sim_cor.2_n10_p100,
                                                  "sd for effective number cor.2 n20 p100" = sd_sim_cor.2_n20_p100,
                                                  "sd for effective number cor.2 n100 p100" = sd_sim_cor.2_n100_p100,
                                                  "sd for effective number cor.2 n200 p100" = sd_sim_cor.2_n200_p100,
                                                  "sd for effective number cor.2 n1000 p100" = sd_sim_cor.2_n1000_p100)
sd_gao_effective_number_cor.2_different_n


# standard deviation for  simulation with cor = 0.2 changing n p = 4000
sd_sim_cor.2_n10_p4000 = sd(sim_cor.2_n10_p4000)
sd_sim_cor.2_n10_p4000

sd_sim_cor.2_n20_p4000 = sd(sim_cor.2_n20_p4000)
sd_sim_cor.2_n20_p4000

sd_sim_cor.2_n100_p4000 = sd(sim_cor.2_n100_p4000)
sd_sim_cor.2_n100_p4000

sd_sim_cor.2_n200_p4000 = sd(sim_cor.2_n200_p4000)
sd_sim_cor.2_n200_p4000

sd_sim_cor.2_n1000_p4000 = sd(sim_cor.2_n1000_p4000)
sd_sim_cor.2_n1000_p4000


sd_gao_effective_number_cor.2_different_n_p4000 = cbind("sd for effective number cor.2 n10 p4000" = sd_sim_cor.2_n10_p4000,
                                                        "sd for effective number cor.2 n20 p4000" = sd_sim_cor.2_n20_p4000,
                                                        "sd for effective number cor.2 n100 p4000" = sd_sim_cor.2_n100_p4000,
                                                        "sd for effective number cor.2 n200 p4000" = sd_sim_cor.2_n200_p4000,
                                                        "sd for effective number cor.2 n1000 p4000" = sd_sim_cor.2_n1000_p4000)
sd_gao_effective_number_cor.2_different_n_p4000

### make a table of all the standard deviations for my simulation

sd_gao_effective_simulated_data = cbind(sd_gao_effective_number_cor0,
                                        sd_gao_effective_number_cor1, 
                                        sd_gao_effective_number_cor.8_different_p, 
                                        sd_gao_effective_number_cor.8_different_n,
                                        sd_gao_effective_number_cor.8_different_n_p4000,
                                        sd_gao_effective_number_cor.5_different_p,
                                        sd_gao_effective_number_cor.5_different_n,
                                        sd_gao_effective_number_cor.5_different_n_p4000,
                                        sd_gao_effective_number_cor.2_different_p,
                                        sd_gao_effective_number_cor.2_different_n,
                                        sd_gao_effective_number_cor.2_different_n_p4000)

sd_gao_effective_simulated_data
###-------------------------------------------------------------------------------------------------------------------------------

# time to find the standard deviation for my block simulations

### -----------------------------------------------------------------------------------------------------------------------------



# standard deviation for 2 block simulation 
sd_gao_effective_number_2.blocks = cbind("effective number 2 blocks p = 10 r = 1 n = 100"    = sd((block_gao_p10_b2_r1_n100)),
                                         "effective number 2 blocks p = 100 r = 1 n = 100"   = sd((block_gao_p100_b2_r1_n100)),
                                         "effective number 2 blocks p = 1000 r = 1 n = 100"  = sd((block_gao_p1000_b2_r1_n100)),
                                         "effective number 2 blocks p = 1000 r = 1 n = 1000" = sd((block_gao_p1000_b2_r1_n1000)))

sd_gao_effective_number_2.blocks



# standard deviation for 3 block simulation 
sd_gao_effective_number_3.blocks = cbind("effective number 3 blocks p = 10 r = 1 n = 100"    = sd((block_gao_p10_b3_r1_n100)),
                                         "effective number 3 blocks p = 100 r = 1 n = 100"   = sd((block_gao_p100_b3_r1_n100)),
                                         "effective number 3 blocks p = 1000 r = 1 n = 100"  = sd((block_gao_p1000_b3_r1_n100)),
                                         "effective number 3 blocks p = 1000 r = 1 n = 1000" = sd((block_gao_p1000_b3_r1_n1000)))

sd_gao_effective_number_3.blocks


# standard deviation for 5 block simulation 
sd_gao_effective_number_5.blocks = cbind("effective number 5 blocks p = 10 r = 1 n = 100"    = sd((block_gao_p10_b5_r1_n100)),
                                         "effective number 5 blocks p = 100 r = 1 n = 100"   = sd((block_gao_p100_b3_r1_n100)),
                                         "effective number 5 blocks p = 1000 r = 1 n = 100"  = sd((block_gao_p1000_b5_r1_n100)),
                                         "effective number 5 blocks p = 1000 r = 1 n = 1000" = sd((block_gao_p1000_b5_r1_n1000)))

sd_gao_effective_number_5.blocks


# standard deviation for 5 block simulation 
sd_gao_effective_number_5.blocks = cbind("effective number 5 blocks p = 10 r = 1 n = 100"    = sd((block_gao_p10_b5_r1_n100)),
                                         "effective number 5 blocks p = 100 r = 1 n = 100"   = sd((block_gao_p100_b3_r1_n100)),
                                         "effective number 5 blocks p = 1000 r = 1 n = 100"  = sd((block_gao_p1000_b5_r1_n100)),
                                         "effective number 5 blocks p = 1000 r = 1 n = 1000" = sd((block_gao_p1000_b5_r1_n1000)))

sd_gao_effective_number_5.blocks


# standard deviation for 10 block simulation 
sd_gao_effective_number_10.blocks = cbind("effective number 10 blocks p = 10 r = 1 n = 100"    = sd((block_gao_p10_b10_r1_n100)),
                                          "effective number 10 blocks p = 100 r = 1 n = 100"   = sd((block_gao_p100_b3_r1_n100)),
                                          "effective number 10 blocks p = 1000 r = 1 n = 100"  = sd((block_gao_p1000_b10_r1_n100)),
                                          "effective number 10 blocks p = 1000 r = 1 n = 1000" = sd((block_gao_p1000_b10_r1_n1000)))

sd_gao_effective_number_10.blocks



#### time to make a nice table for my standard deviations of my block simulations
sd_gao_effective_simulated_data = cbind(sd_gao_effective_number_2.blocks,
                                        sd_gao_effective_number_3.blocks,
                                        sd_gao_effective_number_5.blocks,
                                        sd_gao_effective_number_10.blocks)

sd_gao_effective_simulated_data


# -------------------------------------------------------------------------------------------------------------------------------------------------

# time to find the range of all my simulated values

# cor 0
range_gao_effective_number_cor0 = cbind("effective number cor0 n10 p10"      = (range(sim_cor0_n10_p10)),
                                  "effective number cor0 n100 p10"     = (range(sim_cor0_n100_p10)),
                                  "effective number cor0 n1000 p10"    = (range(sim_cor0_n1000_p10)),
                                  "effective number cor0 n10000 p10"   = (range(sim_cor0_n10000_p10)),
                                  "effective number cor0 n10 p100"     = (range(sim_cor0_n10_p100)),
                                  "effective number cor0 n100 p100"    = (range(sim_cor0_n100_p100)),
                                  "effective number cor0 n1000 p100"   = (range(sim_cor0_n1000_p100)),
                                  "effective number cor0 n10000 p100"  = (range(sim_cor0_n10000_p100)),
                                  "effective number cor0 n10 p1000"    = (range(sim_cor0_n10_p1000)),
                                  "effective number cor0 n100 p1000"   = (range(sim_cor0_n100_p1000)),
                                  "effective number cor0 n1000 p1000"  = (range(sim_cor0_n1000_p1000)),
                                  "effective number cor0 n10000 p1000" = (range(sim_cor0_n10000_p1000)),
                                  "effective number cor0 n10 p4000"    = (range(sim_cor0_n10_p4000)),
                                  "effective number cor0 n100 p4000"   = (range(sim_cor0_n100_p4000)),
                                  "effective number cor0 n1000 p4000"  = (range(sim_cor0_n1000_p4000)),
                                  "effective number cor0 n10000 p4000" = (range(sim_cor0_n10000_p4000)))



range_gao_effective_number_cor0


# start with 100% correlated matrix and simulate my data with number of variables = 100,


range_gao_effective_number_cor1 = cbind("effective number cor1 n100 p10" = (range(sim_cor1_n100_p10)),
                                  "effective number cor1 n100 p100" = (range(sim_cor1_n100_p100)),
                                  "effective number cor1 n100 p4000" = (range(sim_cor1_n100_p4000)))
range_gao_effective_number_cor1


# number of observations = 100 changing number of vairables and cor = .8


range_gao_effective_number_cor.8_different_p = cbind("effective number cor.8 n100 p10" = (range(sim_cor.8_n100_p10)),
                                               "effective number cor.8 n100 p20" = (range(sim_cor.8_n100_p20)),
                                               "effective number cor.8 n100 p100" = (range(sim_cor.8_n100_p100)),
                                               "effective number cor.8 n100 p200" = (range(sim_cor.8_n100_p200)),
                                               "effective number cor.8 n100 p1000" = (range(sim_cor.8_n100_p1000)),
                                               "effective number cor.8 n100 p4000" = (range(sim_cor.8_n100_p4000)))
range_gao_effective_number_cor.8_different_p


# number of variables = 100 changing number of observations and cor = .8



range_gao_effective_number_cor.8_different_n = cbind("effective number cor.8 n10 p100" = (range(sim_cor.8_n10_p100)),
                                               "effective number cor.8 n20 p100" = (range(sim_cor.8_n20_p100)),
                                               "effective number cor.8 n100 p100" = (range(sim_cor.8_n100_p100)),
                                               "effective number cor.8 n200 p100" = (range(sim_cor.8_n200_p100)),
                                               "effective number cor.8 n1000 p100" = (range(sim_cor.8_n1000_p100)))
range_gao_effective_number_cor.8_different_n


# number of variables = 4000 changing number of observations and cor = .8



range_gao_effective_number_cor.8_different_n_p4000 = cbind("effective number cor.8 n10 p4000" = (range(sim_cor.8_n10_p4000)),
                                                     "effective number cor.8 n20 p4000" = (range(sim_cor.8_n20_p4000)),
                                                     "effective number cor.8 n100 p4000" = (range(sim_cor.8_n100_p4000)),
                                                     "effective number cor.8 n200 p4000" = (range(sim_cor.8_n200_p4000)),
                                                     "effective number cor.8 n1000 p4000" = (range(sim_cor.8_n1000_p4000)))
range_gao_effective_number_cor.8_different_n_p4000


# number of observations = 100 changing number of vairables and cor = .5


range_gao_effective_number_cor.5_different_p = cbind("effective number cor.5 n100 p10" = (range(sim_cor.5_n100_p10)),
                                               "effective number cor.5 n100 p20" = (range(sim_cor.5_n100_p20)),
                                               "effective number cor.5 n100 p100" = (range(sim_cor.5_n100_p100)),
                                               "effective number cor.5 n100 p200" = (range(sim_cor.5_n100_p200)),
                                               "effective number cor.5 n100 p1000" = (range(sim_cor.5_n100_p1000)),
                                               "effective number cor.5 n100 p4000" = (range(sim_cor.5_n100_p4000)))
range_gao_effective_number_cor.5_different_p

# number of variables = 100 changing number of observations and cor = .5



range_gao_effective_number_cor.5_different_n = cbind("effective number cor.5 n10 p100" = (range(sim_cor.5_n10_p100)),
                                               "effective number cor.5 n20 p100" = (range(sim_cor.5_n20_p100)),
                                               "effective number cor.5 n100 p100" = (range(sim_cor.5_n100_p100)),
                                               "effective number cor.5 n200 p100" = (range(sim_cor.5_n200_p100)),
                                               "effective number cor.5 n1000 p100" = (range(sim_cor.5_n1000_p100)))
range_gao_effective_number_cor.5_different_n


# number of variables = 4000 changing number of observations and cor = .5


range_gao_effective_number_cor.5_different_n_p4000 = cbind("effective number cor.5 n10 p4000" = (range(sim_cor.5_n10_p4000)),
                                                     "effective number cor.5 n20 p4000" = (range(sim_cor.5_n20_p4000)),
                                                     "effective number cor.5 n100 p4000" = (range(sim_cor.5_n100_p4000)),
                                                     "effective number cor.5 n200 p4000" = (range(sim_cor.5_n200_p4000)),
                                                     "effective number cor.5 n1000 p4000" = (range(sim_cor.5_n1000_p4000)))
range_gao_effective_number_cor.5_different_n_p4000


# number of observations = 100 changing number of vairables and cor = .2


range_gao_effective_number_cor.2_different_p = cbind("effective number cor.2 n100 p10" = (range(sim_cor.2_n100_p10)),
                                               "effective number cor.2 n100 p20" = (range(sim_cor.2_n100_p20)),
                                               "effective number cor.2 n100 p100" = (range(sim_cor.2_n100_p100)),
                                               "effective number cor.2 n100 p200" = (range(sim_cor.2_n100_p200)),
                                               "effective number cor.2 n100 p1000" = (range(sim_cor.2_n100_p1000)),
                                               "effective number cor.2 n100 p4000" = (range(sim_cor.2_n100_p4000)))
range_gao_effective_number_cor.2_different_p

# number of variables = 100 changing number of observations and cor = .2


range_gao_effective_number_cor.2_different_n = cbind("effective number cor.2 n10 p100" = (range(sim_cor.2_n10_p100)),
                                               "effective number cor.2 n20 p100" = (range(sim_cor.2_n20_p100)),
                                               "effective number cor.2 n100 p100" = (range(sim_cor.2_n100_p100)),
                                               "effective number cor.2 n200 p100" = (range(sim_cor.2_n200_p100)),
                                               "effective number cor.2 n1000 p100" = (range(sim_cor.2_n1000_p100)))
range_gao_effective_number_cor.2_different_n


# number of variables = 4000 changing number of observations and cor = .2



range_gao_effective_number_cor.2_different_n_p4000 = cbind("effective number cor.2 n10 p4000" = (range(sim_cor.2_n10_p4000)),
                                                     "effective number cor.2 n20 p4000" = (range(sim_cor.2_n20_p4000)),
                                                     "effective number cor.2 n100 p4000" = (range(sim_cor.2_n100_p4000)),
                                                     "effective number cor.2 n200 p4000" = (range(sim_cor.2_n200_p4000)),
                                                     "effective number cor.2 n1000 p4000" = (range(sim_cor.2_n1000_p4000)))
range_gao_effective_number_cor.2_different_n_p4000



### BLOCKS RANGE ###



range_gao_effective_number_2.blocks = cbind("effective number 2 blocks p = 10 r = 1 n = 100"    = (range(block_gao_p10_b2_r1_n100)),
                                      "effective number 2 blocks p = 100 r = 1 n = 100"   = (range(block_gao_p100_b2_r1_n100)),
                                      "effective number 2 blocks p = 1000 r = 1 n = 100"  = (range(block_gao_p1000_b2_r1_n100)),
                                      "effective number 2 blocks p = 1000 r = 1 n = 1000" = (range(block_gao_p1000_b2_r1_n1000)))

range_gao_effective_number_2.blocks




range_gao_effective_number_3.blocks = cbind("effective number 3 blocks p = 10 r = 1 n = 100"   = (range(block_gao_p10_b3_r1_n100)),
                                      "effective number 3 blocks p = 100 r = 1 n = 100"  = (range(block_gao_p100_b3_r1_n100)),
                                      "effective number 3 blocks p = 1000 r = 1 n = 100" = (range(block_gao_p1000_b3_r1_n100)),
                                      "effective number 3 blocks p = 1000 r = 1 n = 1000" = (range(block_gao_p1000_b3_r1_n1000)))

range_gao_effective_number_3.blocks




range_gao_effective_number_5.blocks = cbind("effective number 5 blocks p = 10 r = 1 n = 100"   = (range(block_gao_p10_b5_r1_n100)),
                                      "effective number 5 blocks p = 100 r = 1 n = 100"  = (range(block_gao_p100_b5_r1_n100)),
                                      "effective number 5 blocks p = 1000 r = 1 n = 100" = (range(block_gao_p1000_b5_r1_n100)),
                                      "effective number 5 blocks p = 1000 r = 1 n = 1000" = (range(block_gao_p1000_b5_r1_n1000)))

range_gao_effective_number_5.blocks



range_gao_effective_number_10.blocks = cbind("effective number 10 blocks p = 10 r = 1 n = 100"   = (range(block_gao_p10_b10_r1_n100)),
                                       "effective number 10 blocks p = 100 r = 1 n = 100"  = (range(block_gao_p100_b10_r1_n100)),
                                       "effective number 10 blocks p = 1000 r = 1 n = 100" = (range(block_gao_p1000_b10_r1_n100)),
                                       "effective number 10 blocks p = 1000 r = 1 n = 1000" = (range(block_gao_p1000_b10_r1_n1000)))

range_gao_effective_number_10.blocks

