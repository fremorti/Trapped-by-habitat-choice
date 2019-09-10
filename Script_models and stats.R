library(car)
library(tidyr)
library(plyr)
library(ggplot2)
library(rethinking)
library(gridExtra)

#################################
# Habitat choice each selection #
#################################

#import dataset and resolve some issues with dataframe
choice <- read.csv("choice.csv", header = 1)
choice$treatment <- factor(choice$treatment, c('R', 'T', 'C'))
colnames(choice)[1] <- 'code'
choice$code <-as.factor(choice$code)
choice$n <-as.factor(choice$n)
#convert selection regime to 2 binary variables
choice$treatT <- choice$treatment == 'T'
choice$treatC <- choice$treatment == 'C'
str(choice)


#model the expected number of mites on tomato each selection round (selection)
#all other model descriptions follow following template.
m_ch3 <- map2stan(
  alist(
    tomato ~ dbetabinom(total, p, theta),                                                     #model tomato as distributed according to a betabinomial distribution
    logit(p) <- a[code] + aC*treatC + aT*treatT +(b[code] + bC*treatC + bT*treatT)*selection, #we use a logit function as link function on the proportion (p), which depends on the selection round (selection), the selection regime (treaT, treatC) and its interaction
                                                                                              #there is a varying intercept and varying slope for each replica
    c(a, b)[code] ~ dmvnorm2(c(a_0, b_0), sigma_code, rho),                                   #multivariate normal distribution of varying intercepts and slopes
    a_0 ~ dnorm(0,2),                                                              #prior distribution of varying intercepts
    c(aC, aT, b_0, bC, bT) ~ dnorm(0, 0.5),                                                   #prior distributions of cucumber choice parameter, tomato choice parameter, varying slopes, cucumber slope parameter and tomato choice slope parameter
    sigma_code ~ dcauchy(0, 2),                                                               #prior distribution for standard deviation of varying slopes and varying intercepts
    theta ~ dexp(1),                                                                          #prior distribution of betabinomial shape (theta)
    rho ~ dlkjcorr(2)                                                                         #prior of covariance of multivariate normal distribution
  ) ,
  data =  choice[!is.na(choice$tomato),],                                                      #omit NA's due to lost data
  constraints = list(theta="lower=0"),                                                        #contrain theta to positive values 
  chains = 2, cores = 2, iter = 10000, warmup = 5000)                                         #define the number of markov chains, cpu cores to use, the number of iterations to be run and the number of warmup chain iterations before sampling

#check markov chain mixing
plot(m_ch3)
par(mfrow = c(1,1))
#check model estimates and Gelman-Rubin convergence criterion (Rhat, should be 1.00)
precis(m_ch3, depth = 2)

#sample posterior slopes (similar to link)
post <- extract.samples(m_ch3)
pbR <- post$b_0
pbC <- post$b_0 + post$bC
pbT <- post$b_0 + post$bT

#plot estimated slopes
p <- data.frame('R'= pbR, 'C' = pbC, 'T' = pbT)
p_ <- gather(p, key = 'treatment', value = 'post')
(pred_slopes <- ggplot(data = p_, aes(x = treatment, y = post, color = treatment))+
    geom_violin(alpha = 0.3, draw_quantiles = c(0.09, 0.5, 0.91), aes(fill = treatment))+
    ggtitle('posterior predicted log-odds slope')+
    ylab('predicted log-odds slope')+
    scale_color_brewer( palette = "Dark2")+
    scale_fill_brewer( palette = "Dark2")+
    theme(legend.position='none'))

#plot estimated differences in slopes
pdiff <- data.frame('C_R'= pbC-pbR, 'C_T' = pbC-pbT, 'T_R' = pbT-pbR)
pdiff_ <- gather(pdiff, key = 'pairwise_difference', value = 'post')
(pred_diff_slopes <- ggplot(data = pdiff_, aes(x = pairwise_difference, y = post))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91))+
    ggtitle('posterior predicted differences in log-odds slope')+
    ylab('predicted differences')+
    xlab('pairs of treatments'))

#plot data on tomato preference
choice$treatment <- factor(choice$treatment, levels = c('C', 'R', 'T'))
(choice_sep2 <- ggplot(data = choice, aes(x=selection, y=proptomato, color = treatment, fill = treatment))+
    geom_point()+
    geom_smooth(method = 'lm')+ #smooth using a linear trend (great to sea increase/decrease + good approximation for log odds around 0.5)
    facet_grid(treatment~.)+
    ylab('proportion on tomato')+
    xlab('selection round')+
    scale_fill_brewer( palette = "Dark2")+
    scale_color_brewer( palette = "Dark2")+
    theme(legend.position="top"))

#combine plots
grid.arrange(grobs = list(choice_sep2, pred_slopes, pred_diff_slopes), widths = c(1,2), layout_matrix = rbind(c(1, 2),c(1, 3)))

######################
# Habitat imprinting #
######################
# visualisation #

#induced habitat choice (habitat imprinting) data after experiment
#import dataset and resolve some issues with dataframe
hi <- read.csv("HI_post.csv", header = 1)
hi$rep <- as.factor(hi$rep)
colnames(hi)[1] <- 'treatment'
str(hi)

#convert maternal habitat (mat, which is also the developmental habitat) and selection regime (treatment) to binary variables
hi$matC <- hi$mat == 'C'
hi$treatC <- hi$treatment == 'C'
hi$treatT <- hi$treatment == 'T'

#model the expected number of mites on tomato after a day (tomato2) depending on developmental habitat
m_hi1 <- map2stan(
  alist(
    tomato2 ~ dbetabinom(total2, p, theta), 
    logit(p) <- a[treatment] + b1*matC,     
    a[treatment] ~ dnorm(a, sigma_t),       
    a ~ dnorm(0, 2),                      
    b1 ~ dnorm(0, 2),                       
    theta ~ dexp(1),                        
    sigma_t ~ dcauchy(0, 1)                 
    
  ) ,
  data =  hi, chains = 2, cores = 2, iter = 10000, warmup = 1000) #define data set that is used, the number of markov chains, cpu cores to use, the number of iterations to be run and the number of warmup chain iterations before sampling

#check estimated model parameters
precis(m_hi1, depth = 2)
plot(precis(m_hi1, depth = 2))

#sample posterior 
l_after <- link(m_hi1, data_ = list(matC = c(0,1)))
#construct a dataframe with estimated expected tomato2 values (a) and estimated expected differences between both developmental habitats (b)
a_after <- gather(data.frame('T' = l_after[,1], 'C' = l_after[,2]), key = 'treatment', value = 't') 
b_after <- data.frame('delta_after' = l_after[, 1]-l_after[,2])
#proportion choosing tomato according to developmental environment
(bp_hi_post <- ggplot(data = hi, aes(x = mat, y = rel_tomato2))+
    geom_boxplot()+
    geom_violin(data = a_after, draw_quantiles = c(0.09, 0.5, 0.91), aes(x = treatment, y = t, width = 0.5))+
    scale_y_continuous(limits = c(0, 1))+
    ggtitle('tomato preference end')+
    xlab('developmental host')+
    ylab('proportion on tomato')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))
  
#induced habitat choice (habitat imprinting) data before experiment
#import dataset and resolve some issues with dataframe
hi_pre <- read.csv("HI_pre.csv", header = 1)
hi_pre$rep <- as.factor(hi_pre$rep)
colnames(hi_pre)[1] <- 'mat'
str(hi_pre)


hi_pre$total <- hi_pre$tomato+hi_pre$cucumber #total of mites tested
hi_pre$matC <- hi_pre$mat == 'C'

#model the expected number of mites on tomato after a day (tomato) depending on developmental habitat
m_hipre <- map2stan(
  alist(
    tomato ~ dbetabinom(total, p, theta),
    logit(p) <- a + b1*matC,
    a <- dnorm(0,2),
    b1 <- dnorm(0, 2),
    theta ~ dexp(1)
  ) ,
  data =  hi_pre, chains = 2, cores = 2, iter = 10000, warmup = 5000)

#check estimated variables
precis(m_hipre, depth = 2)

#sample posterior
l <- link(m_hipre, data = list(matC = c(0,1)))

#estimated proportion choosing tomato according to developmental environment
a_pre <- gather(data.frame('T' = l[,1], 'C' = l[,2]), key = 'treatment', value = 't')
(bp_hi_pre <- ggplot(data = hi_pre, aes(x = mat, y = rel_tomato))+
    geom_boxplot()+
    geom_violin(data = a_pre, draw_quantiles = c(0.09, 0.5, 0.91), aes(x = treatment, y=t, width = 0.5))+
    scale_y_continuous(limits = c(0, 1))+
    ggtitle('tomato preference start')+
    xlab('developmental host')+
    ylab('proportion on tomato')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#estimated difference in proportion choosing tomato between developmental environment at start and end of experiment
b <- gather(data.frame('start' = l[, 1]-l[, 2], 'end' = l_after[, 1]-l_after[, 2]), key = 'time', value = 't')
b$time <- factor(b$time, c('start', 'end'))
(delta_hi <- ggplot(data = b, aes(x = time, y = t))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91))+
    scale_y_continuous(limits = c(-0.35, 1))+
    ggtitle('posterior predicted differences')+
    xlab('time')+
    ylab('difference in proportion on tomato')+
    geom_vline(xintercept=1.5, linetype= 'dotted'))

#combine graphs
grid.arrange(bp_hi_pre, delta_hi, bp_hi_post, nrow = 1)



##########
#LH data #
##########

#life-history data at end
#import dataset and resolve some issues with dataframe
LH <- read.table("LH.txt", header = 1)
LH$treatment = factor(LH$treatment,c("T","R","C", "S"))
LH$treplicate = as.factor(LH$treplicate)
LH$patchC <- LH$patch == 'C'
LH$patchT <- LH$patch == 'T'
LH$patchB <- LH$patch == 'B'
LH$trep <- paste0(LH$treatment, LH$treplicate)
str(LH)



#life-history data at start
#import dataset and resolve some issues with dataframe
LH_pre_ <- read.csv('LH_pre.csv', header = 1, sep = ';')
colnames(LH_pre_)[1] = 'cg'
LH_pre_$nr <- as.factor(LH_pre_$nr)
str(LH_pre_)
LH_pre <- spread(LH_pre_, key = 'life.stage', value = count) #convert dataframe so that every life stage has its column
LH_pre$total <- LH_pre$af + LH_pre$am + LH_pre$d             #column of all adults and deutonymphs
LH_pre$nomales <- LH_pre$af + LH_pre$d                       #column of all deutonymphs and female adults
LH_pre$d6 <- LH_pre$e + LH_pre$l + LH_pre$p                  #column of the total of eggs, larvae and protonymphs: the total count for counts at day6
#tested host as binary variables
LH_pre$patchC <- LH_pre$plant == 'C'
LH_pre$patchT <- LH_pre$plant == 'T'

#data on reproductive success (day 12)

#model females and deutonymphs (reproductive success, X12nomales) dependent on tested host at end of experiment
m_LH <- map2stan(
  alist(
    X12nomales ~ dgampois(lambda, theta),
    log(lambda) <- a_trep[trep] + bC * patchC + bT * patchT,
    a_trep[trep] ~ dnorm(a, sigma_trep),
    a ~ dnorm(0, 5),
    c(bC, bT) ~ dnorm(0, 2),
    theta ~ dexp(1),
    sigma_trep ~ dcauchy(0, 2)
  ) ,
  data =  LH, chains = 2, cores = 2, iter = 10000, warmup = 1000)

#check model parameters
precis(m_LH, depth = 2)
#sample posterior
post <- extract.samples(m_LH)
simB <- exp(post$a)
simC <- exp(post$a+post$bC)
simT <- exp(post$a+post$bT)

#plot estimated reproductive success on different hosts at end
a <- gather(data.frame('B' = simB, 'T' = simT, 'C' = simC), key = 'treatment', value = 't')
(fit_after <- ggplot(data = LH, aes(x = patch, y = X12nomales))
  +geom_boxplot()+
    expand_limits(y = 0)+
    geom_violin(data = a, draw_quantiles = c(0.09, 0.5, 0.91), aes(x = treatment, y = t, width = 0.5))+
    ggtitle('reproductive success end')+
    xlab('host')+
    ylab('deutonymphs produced in 12 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot estimated differences of reproductive success between different hosts at end
b <- gather(data.frame('B_C'=simB-simC, 'B_T'=simB-simT, 'T_C'=simT-simC), key = 'pairs', value = 'pairwise_differences')
(fit_diff_after <- ggplot(data = b, aes(x = pairs, y = pairwise_differences))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91))+
  expand_limits(y = 0)+
  ggtitle('posterior predicted differences end')+
  xlab('pairs of treatments'))


#model females and deutonymphs (reproductive success, X12nomales) dependent on tested host at start of experiment
m_LHpre <- map2stan(
  alist(
    nomales ~ dgampois(lambda, theta),
    log(lambda) <- a + bC * patchC + bT * patchT,
    a ~ dnorm(0,5),
    theta ~ dexp(1),
    c(bC, bT) ~ dnorm(0, 1)
  ) ,
  data =  LH_pre[LH_pre$day == 12 & LH_pre$cg == 'tomato',], chains = 2, cores = 2, iter = 10000, warmup = 1000)

#check model parameters
precis(m_LHpre)
#sample posterior
post <- extract.samples(m_LHpre)
simB <- exp(post$a)
simC <- exp(post$a+post$bC)
simT <- exp(post$a+post$bT)

#plot estimated reproductive success on different hosts at start
a_ <- gather(data.frame('B' = simB, 'T' = simT, 'C' = simC), key = 'treatment', value = 't')
(fit_before <- ggplot(data = LH_pre[LH_pre$day == 12 & LH_pre$cg == 'tomato',], aes(x = plant, y = nomales))+
    geom_boxplot()+
    geom_violin(data = a_, draw_quantiles = c(0.09, .5, 0.91), aes(x = treatment, y = t, width = 0.5, trim = 1))+
    ggtitle('reproductive success start')+
    xlab('host')+
    ylab('deutonymphs produced in 12 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot estimated differences of reproductive success between different hosts at start
b_ <- gather(data.frame('B_C'=simB-simC, 'B_T'=simB-simT, 'T_C'=simT-simC), key = 'pairs', value = 'pairwise_differences')
(fit_diff_before <- ggplot(data = b_, aes(x = pairs, y = pairwise_differences))+
  geom_violin(draw_quantiles = c(0.09, .5, 0.91))+
  expand_limits(y = 0)+
  ggtitle('posterior predicted differences start')+
  xlab('pairs of treatments'))

#combine plot
grid.arrange(fit_before, fit_after, fit_diff_before, fit_diff_after, nrow = 2)

#Appendix 1: data on fertility (day6)

#model fertility at end
m_LH6 <- map2stan(
  alist(
    X6total ~ dgampois(lambda, theta),
    log(lambda) <- a_trep[trep] + bC * patchC + bT * patchT,
    a_trep[trep] ~ dnorm(a, sigma_trep),
    a ~ dnorm(0, 5),
    c(bC, bT) ~ dnorm(0, 2),
    theta ~ dexp(1),
    sigma_trep ~ dcauchy(0, 2)
  ) ,
  data =  LH, chains = 2, cores = 2, iter = 10000, warmup = 1000)


post <- extract.samples(m_LH6)
#sample posterior
simB <- exp(post$a)
simC <- exp(post$a+post$bC)
simT <- exp(post$a+post$bT)

#plot fertility at end
a <- gather(data.frame('B' = simB, 'T' = simT, 'C' = simC), key = 'treatment', value = 't')
(fit6_after <- ggplot(data = LH, aes(x = patch, y = X6total))
  +geom_boxplot()+
    expand_limits(y = 0)+
    geom_violin(data = a, aes(x = treatment, y = t, width = 0.5))+
    ggtitle('fertility end')+
    xlab('host')+
    ylab('eggs produced in 6 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot differences in fertility at end
b <- gather(data.frame('B_C'=simB-simC, 'B_T'=simB-simT, 'T_C'=simT-simC), key = 'pairs', value = 'pairwise_differences')
(fit6_diff_after <- ggplot(data = b, aes(x = pairs, y = pairwise_differences))+
    geom_violin()+
    expand_limits(y = 0)+
    ggtitle('posterior predicted differences end')+
    xlab('pairs of treatments'))


#model fertility at start
m_LHpre6 <- map2stan(
  alist(
    d6 ~ dgampois(lambda, theta),
    log(lambda) <- a + bC * patchC + bT * patchT,
    a ~ dnorm(0,5),
    theta ~ dexp(1),
    c(bC, bT) ~ dnorm(0, 1)
  ) ,
  data =  LH_pre[LH_pre$day == 6 & LH_pre$cg == 'tomato',], chains = 2, cores = 2, iter = 10000, warmup = 1000)


post <- extract.samples(m_LHpre6)
#sample distribution
simB <- exp(post$a)
simC <- exp(post$a+post$bC)
simT <- exp(post$a+post$bT)

#plot fertility at start
a_ <- gather(data.frame('B' = simB, 'T' = simT, 'C' = simC), key = 'treatment', value = 't')
(fit6_before <- ggplot(data = LH_pre[LH_pre$day == 6 & LH_pre$cg == 'tomato',], aes(x = plant, y = d6))+
    geom_boxplot()+
    geom_violin(data = a_, aes(x = treatment, y = t, width = 0.5, trim = 1))+
    ggtitle('fertility start')+
    xlab('host')+
    ylab('eggs produced in 6 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot differences in fertility at start
b_ <- gather(data.frame('B_C'=simB-simC, 'B_T'=simB-simT, 'T_C'=simT-simC), key = 'pairs', value = 'pairwise_differences')
(fit6_diff_before <- ggplot(data = b_, aes(x = pairs, y = pairwise_differences))+
    geom_violin()+
    expand_limits(y = 0)+
    ggtitle('posterior predicted differences start')+
    xlab('pairs of treatments'))

grid.arrange(fit_before, fit6_after, fit_diff_before, fit_diff_after, nrow = 2)


#Appendix 2: effect of selection regime on LH-data
#plot reproductive success by host tested and selection regime
(S21 <- ggplot(data = LH, aes(x = treatment, y = X12nomales))
  +geom_boxplot()+
    expand_limits(y = 0)+
    ggtitle('reproductive success end')+
    xlab('host')+
    ylab('deutonymphs produced in 12 days')+
    facet_grid(.~patch)+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot fertility by host tested and selection regime
(S22 <- ggplot(data = LH, aes(x = treatment, y = X6total))
  +geom_boxplot()+
    expand_limits(y = 0)+
    ggtitle('fertility end')+
    xlab('host')+
    ylab('eggs produced in 6 days')+
    facet_grid(.~patch)+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))
