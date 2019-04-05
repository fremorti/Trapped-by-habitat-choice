library(car)
library(tidyr)
library(plyr)
library(ggplot2)
library(rethinking)
library(gridExtra)

######################
# Habitat imprinting #
######################
# visualisation #


hi <- read.csv("HI_post.csv", header = 1)
hi$rep <- as.factor(hi$rep)
colnames(hi)[1] <- 'treatment'
str(hi)

hi$matC <- hi$mat == 'C'
hi$treatC <- hi$treatment == 'C'
hi$treatT <- hi$treatment == 'T'

m_hi1 <- map2stan(
  alist(
    tomato2 ~ dbetabinom(total2, p, theta),
    logit(p) <- a[treatment] + b1*matC,
    a[treatment] ~ dnorm(a, sigma_t),
    a ~ dnorm(0, 0.5),
    b1 ~ dnorm(0, 1),
    theta ~ dexp(1),
    sigma_t ~ dcauchy(0, 1)
  ) ,
  data =  hi, chains = 2, cores = 2, iter = 20000, warmup = 1000)

precis(m_hi1, depth = 2)
plot(precis(m_hi1, depth = 2))
l_after <- link(m_hi1, data_ = list(matC = c(0,1)))

a <- gather(data.frame('T' = l_after[,1], 'C' = l_after[,2]), key = 'treatment', value = 't')
b_after <- data.frame('delta_after' = l_after[, 1]-l_after[,2])
#proportion choosing tomato according to developmental environment
(bp_hi_post <- ggplot(data = hi, aes(x = mat, y = rel_tomato2))+
    geom_boxplot()+
    geom_violin(data = a, aes(x = treatment, y = t, width = 0.5))+
    scale_y_continuous(limits = c(0, 1))+
    ggtitle('tomato preference end')+
    xlab('developmental host')+
    ylab('proportion on tomato')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))
  


hi_pre <- read.csv("HI_pre.csv", header = 1)
hi_pre$rep <- as.factor(hi_pre$rep)
colnames(hi_pre)[1] <- 'mat'
str(hi_pre)

hi_pre$total <- hi_pre$tomato+hi_pre$cucumber
hi_pre$matC <- hi_pre$mat == 'C'
m_hipre <- map2stan(
  alist(
    tomato ~ dbetabinom(total, p, theta),
    logit(p) <- a + b1*matC,
    a <- dnorm(0,0.5),
    b1 <- dnorm(0, 1),
    theta ~ dexp(1)
  ) ,
  data =  hi_pre, chains = 2, cores = 2, iter = 4000, warmup = 1000)

precis(m_hipre)
l <- link(m_hipre, data = list(matC = c(0,1)))
a_pre <- gather(data.frame('T' = l[,1], 'C' = l[,2]), key = 'treatment', value = 't')
b_pre <- data.frame('delta' = l[, 1]-l[, 2])

#proportion choosing tomato according to developmental environment
(bp_hi_pre <- ggplot(data = hi_pre, aes(x = mat, y = rel_tomato))+
    geom_boxplot()+
    geom_violin(data = a_pre, aes(x = treatment, y=t, width = 0.5))+
    scale_y_continuous(limits = c(0, 1))+
    ggtitle('tomato preference start')+
    xlab('developmental host')+
    ylab('proportion on tomato')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

b <- gather(data.frame('start' = l[, 1]-l[, 2], 'end' = l_after[, 1]-l_after[, 2]), key = 'time', value = 't')
b$time <- factor(b$time, c('start', 'end'))
(delta_hi <- ggplot(data = b, aes(x = time, y = t))+
    geom_violin()+
    scale_y_continuous(limits = c(-0.35, 1))+
    ggtitle('posterior predicted differences')+
    xlab('time')+
    ylab('difference in proportion on tomato')+
    geom_vline(xintercept=1.5, linetype= 'dotted'))

grid.arrange(bp_hi_pre, delta_hi, bp_hi_post, nrow = 1)


#################################
# Habitat choice each selection #
#################################

choice <- read.csv("choice.csv", header = 1)
choice$treatment <- factor(choice$treatment, c('R', 'T', 'C'))
colnames(choice)[1] <- 'code'
choice$code <-as.factor(choice$code)
choice$n <-as.factor(choice$n)
choice$treatT <- choice$treatment == 'T'
choice$treatC <- choice$treatment == 'C'
str(choice)


m_ch3 <- map2stan(
  alist(
    tomato ~ dbetabinom(total, p, theta),
    logit(p) <- a[code] + aC*treatC + aT*treatT +(b[code] + bC*treatC + bT*treatT)*selection,
    c(a, b)[code] ~ dmvnorm2(c(a_0, b_0), sigma_code, rho),
    a_0 ~ dnorm(0.5,0.5),
    c(aC, aT, b_0, bC, bT) ~ dnorm(0, 0.5),
    sigma_code ~ dcauchy(0, 2),
    theta ~ dexp(1),
    rho ~ dlkjcorr(2)
  ) ,
  data =  choice[!is.na(choice$tomato),],
  constraints = list(theta="lower=0"), chains = 2, cores = 2, iter = 10000, warmup = 1000)

precis(m_ch3, depth = 2)
post <- extract.samples(m_ch3)
pbR <- post$b_0
pbC <- post$b_0 + post$bC
pbT <- post$b_0 + post$bT

p <- data.frame('R'= pbR, 'C' = pbC, 'T' = pbT)
p_ <- gather(p, key = 'treatment', value = 'post')
(pred_slopes <- ggplot(data = p_, aes(x = treatment, y = post, color = treatment))+
  geom_violin()+
  ggtitle('posterior predicted log-odds slope')+
  ylab('predicted log-odds slope')+
    theme(legend.position='none'))


pdiff <- data.frame('C_R'= pbC-pbR, 'C_T' = pbC-pbT, 'T_R' = pbT-pbR)
pdiff_ <- gather(pdiff, key = 'pairwise_difference', value = 'post')
(pred_diff_slopes <- ggplot(data = pdiff_, aes(x = pairwise_difference, y = post))+
  geom_violin()+
  ggtitle('posterior predicted differences in log-odds slope')+
  ylab('predicted differences')+
  xlab('pairs of treatments'))

choice$treatment <- factor(choice$treatment, levels = c('C', 'R', 'T'))
(choice_sep2 <- ggplot(data = choice, aes(x=selection, y=proptomato, color = treatment))+
    geom_point()+
    geom_smooth(method = 'lm')+
    facet_grid(treatment~.)+
    ylab('proportion on tomato')+
    xlab('selection round')+
    theme(legend.position="top"))

grid.arrange(grobs = list(choice_sep2, pred_slopes, pred_diff_slopes), widths = c(1,2), layout_matrix = rbind(c(1, 2),c(1, 3)))


##########
#LH data #
##########
#day12
LH <- read.table("LH.txt", header = 1)
LH$treatment = factor(LH$treatment,c("T","R","C", "S"))
str(LH)

LH_pre_ <- read.csv('LH_pre.csv', header = 1, sep = ';')
colnames(LH_pre_)[1] = 'cg'
LH_pre_$nr <- as.factor(LH_pre_$nr)
str(LH_pre_)
LH_pre <- spread(LH_pre_, key = 'life.stage', value = count)
LH_pre$total <- LH_pre$af + LH_pre$am + LH_pre$d 
LH_pre$nomales <- LH_pre$af + LH_pre$d
LH_pre$d6 <- LH_pre$e + LH_pre$l + LH_pre$p
LH_pre$eggs <- LH_pre$e+LH_pre$p+LH_pre$l+LH_pre$d
LH_pre$patchC <- LH_pre$plant == 'C'
LH_pre$patchT <- LH_pre$plant == 'T'


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


precis(m_LH, depth = 2)
post <- extract.samples(m_LH)
simB <- exp(post$a)
simC <- exp(post$a+post$bC)
simT <- exp(post$a+post$bT)

a <- gather(data.frame('B' = simB, 'T' = simT, 'C' = simC), key = 'treatment', value = 't')
(fit_after <- ggplot(data = LH, aes(x = patch, y = X12nomales))
  +geom_boxplot()+
    expand_limits(y = 0)+
    geom_violin(data = a, aes(x = treatment, y = t, width = 0.5))+
    ggtitle('reproductive success end')+
    xlab('host')+
    ylab('deutonymphs produced in 12 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

b <- gather(data.frame('B_C'=simB-simC, 'B_T'=simB-simT, 'T_C'=simT-simC), key = 'pairs', value = 'pairwise_differences')
(fit_diff_after <- ggplot(data = b, aes(x = pairs, y = pairwise_differences))+
  geom_violin()+
  expand_limits(y = 0)+
  ggtitle('posterior predicted differences end')+
  xlab('pairs of treatments'))



m_LHpre <- map2stan(
  alist(
    nomales ~ dgampois(lambda, theta),
    log(lambda) <- a + bC * patchC + bT * patchT,
    a ~ dnorm(0,5),
    theta ~ dexp(1),
    c(bC, bT) ~ dnorm(0, 1)
  ) ,
  data =  LH_pre[LH_pre$day == 12 & LH_pre$cg == 'tomato',], chains = 2, cores = 2, iter = 10000, warmup = 1000)


precis(m_LHpre)
post <- extract.samples(m_LHpre)
simB <- exp(post$a)
simC <- exp(post$a+post$bC)
simT <- exp(post$a+post$bT)

a_ <- gather(data.frame('B' = simB, 'T' = simT, 'C' = simC), key = 'treatment', value = 't')
(fit_before <- ggplot(data = LH_pre[LH_pre$day == 13 & LH_pre$cg == 'tomato',], aes(x = plant, y = nomales))+
    geom_boxplot()+
    geom_violin(data = a_, aes(x = treatment, y = t, width = 0.5, trim = 1))+
    ggtitle('reproductive success start')+
    xlab('host')+
    ylab('deutonymphs produced in 12 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

b_ <- gather(data.frame('B_C'=simB-simC, 'B_T'=simB-simT, 'T_C'=simT-simC), key = 'pairs', value = 'pairwise_differences')
(fit_diff_before <- ggplot(data = b_, aes(x = pairs, y = pairwise_differences))+
  geom_violin()+
  expand_limits(y = 0)+
  ggtitle('posterior predicted differences start')+
  xlab('pairs of treatments'))

grid.arrange(fit_before, fit_after, fit_diff_before, fit_diff_after, nrow = 2)

#day6
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
simB <- exp(post$a)
simC <- exp(post$a+post$bC)
simT <- exp(post$a+post$bT)

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

b <- gather(data.frame('B_C'=simB-simC, 'B_T'=simB-simT, 'T_C'=simT-simC), key = 'pairs', value = 'pairwise_differences')
(fit6_diff_after <- ggplot(data = b, aes(x = pairs, y = pairwise_differences))+
    geom_violin()+
    expand_limits(y = 0)+
    ggtitle('posterior predicted differences end')+
    xlab('pairs of treatments'))



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
simB <- exp(post$a)
simC <- exp(post$a+post$bC)
simT <- exp(post$a+post$bT)

a_ <- gather(data.frame('B' = simB, 'T' = simT, 'C' = simC), key = 'treatment', value = 't')
(fit6_before <- ggplot(data = LH_pre[LH_pre$day == 6 & LH_pre$cg == 'tomato',], aes(x = plant, y = d6))+
    geom_boxplot()+
    geom_violin(data = a_, aes(x = treatment, y = t, width = 0.5, trim = 1))+
    ggtitle('fertility start')+
    xlab('host')+
    ylab('eggs produced in 6 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

b_ <- gather(data.frame('B_C'=simB-simC, 'B_T'=simB-simT, 'T_C'=simT-simC), key = 'pairs', value = 'pairwise_differences')
(fit6_diff_before <- ggplot(data = b_, aes(x = pairs, y = pairwise_differences))+
    geom_violin()+
    expand_limits(y = 0)+
    ggtitle('posterior predicted differences start')+
    xlab('pairs of treatments'))

grid.arrange(fit_before, fit6_after, fit_diff_before, fit_diff_after, nrow = 2)


#supplement 2
(S21 <- ggplot(data = LH, aes(x = treatment, y = X12nomales))
  +geom_boxplot()+
    expand_limits(y = 0)+
    ggtitle('reproductive success end')+
    xlab('host')+
    ylab('deutonymphs produced in 12 days')+
    facet_grid(.~patch)+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

(S22 <- ggplot(data = LH, aes(x = treatment, y = X6total))
  +geom_boxplot()+
    expand_limits(y = 0)+
    ggtitle('fertility end')+
    xlab('host')+
    ylab('eggs produced in 6 days')+
    facet_grid(.~patch)+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))
