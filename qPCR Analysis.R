###################################################
# STATISTICAL ANALYSIS                            #
#                                                 #
#                                                 #
# 96-flat bottomed well plate                     #
# Values: qPCR 2^ddCp values                      #
#                                                 #
# Format:                                         #
#     treatments  reporter1 reporter2 ...         #
#     t1                                          #
#     t1                                          #
#     t1                                          #
#     t2                                          #
#     t2                                          #
#     ...                                         #
#                                                 #
# Contrasts:                                      #
#     Set them where indicated below              #
#                                                 #
#                                                 #
###################################################

rm (list = ls())
setwd("Set_Your_WorkingDirectory")


library("reshape")
library("ggplot2")
library("pastecs")
library("car")
library("plyr")

## Prepare data ########################
#get data into R
fname <- file.choose()
data <- read.csv(fname,header=T)

#convert to long-format
dataL <- melt(id="treatment",data=data)
  #keep data in order
  dataL$treatment <- factor(dataL$treatment, levels=unique(dataL$treatment))
  dataL$variable <- factor(dataL$variable, labels=c("29a","29b1","NQO1"))

#split data into single experiments
d1 <- data[,c(1,2)]
d2 <- data[,c(1,3)]
d3 <- data[,c(1,4)]

p63 <- factor(c(rep(0,6),rep(1,6)),labels=c("-","+"))
Nrf2 <- factor(c(0,0,0,1,1,1,0,0,0,1,1,1),labels=c("-","+"))
d1 <- cbind(d1,p63,Nrf2)
d2 <- cbind(d2,p63,Nrf2)
d3 <- cbind(d3,p63,Nrf2)


## Plot the data ########################
#calculate standard positions for errorbars
len <- rep(stack(cast(dataL, variable~treatment, mean))$values,each=3)
sdd <- rep(stack(cast(dataL, variable~treatment, sd))$values,each=3)

low <- len-sdd
high <- len+sdd

#plot the graph
box.plot <- ggplot(dataL,aes(variable,value,fill=variable))
box.plot + stat_summary(fun.y="mean", geom="bar", position="dodge") +
  labs(x="",y="RNA levels, fold change over RPL27",title="pri-miRNA expression") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin=low, ymax=high), width=.2) +
  facet_wrap(~treatment) +
  theme(legend.text=element_text(size=10), legend.title=element_blank())

## test robustness  ########################
#look at data: do they make sense, are the variances equal?
#29a
by(d1$X29a,d1$treatment,stat.desc)
leveneTest(d1$X29a,d1$treatment,center=mean)
#29b1
by(d2$X29b1,d2$treatment,stat.desc)
leveneTest(d2$X29b1,d2$treatment,center=mean)
#29b1a
by(d3$NQO1,d3$treatment,stat.desc)
leveneTest(d3$NQO1,d3$treatment,center=mean)




## perform two-way ANOVA  ########################
#29a
  #calculate ANOVA according to your fit
    contrasts(d1$Nrf2) <- c(-1,1)
    contrasts(d1$p63) <- c(-1,1)
    
    fit1 <- aov(X29a ~ p63*Nrf2, data=d1)
    Anova(fit1, type="III")
    summary.lm(fit1)

#29b1
  #calculate ANOVA according to your fit
    contrasts(d2$Nrf2) <- c(-1,1)
    contrasts(d2$p63) <- c(-1,1)
    
    fit2 <- aov(X29b1 ~ p63*Nrf2, data=d2)
    Anova(fit2, type="III")
    summary.lm(fit2)

#29b1a
  #calculate ANOVA according to your fit
    contrasts(d3$Nrf2) <- c(-1,1)
    contrasts(d3$p63) <- c(-1,1)
    
    fit3<- aov(X29b1a ~ p63*Nrf2, data=d3)
    Anova(fit3, type="III")
    summary.lm(fit3)

