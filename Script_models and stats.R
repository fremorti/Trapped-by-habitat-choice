library(car)
library(tidyverse)
library(ggplot2)
library(rethinking)
library(brms)
library(gridExtra)

#################################
# Habitat choice each selection #
#################################

#import dataset and resolve some issues with dataframe
choice <- read.csv("choice.csv", header = 1)
colnames(choice)[1] <- 'code'
choice$code <-as.factor(choice$code)
choice$n <-as.factor(choice$n)
#convert selection regime to 2 binary variables
choice$treatT <- choice$treatment == 'T'
choice$treatC <- choice$treatment == 'C'
str(choice)

#define a beta binomial distribution with logit link

betabinom <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)
stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

#model the expected number of mites on tomato each selection round (selection)
#all other model descriptions follow following template.
stan(control = list(adapt_delta = 0.99))
fithc1 <- brm(data = choice[!is.na(choice$tomato),], family = betabinom ,
              tomato | vint(total) ~ 1 + treatment*selection + (1 + selection|code),
              prior = c(prior(normal(0,2), class = Intercept),
                        prior(normal(0,2), class = b),
                        prior(cauchy(0,1), class = sd),
                        prior(lkj(2), class = cor),
                        prior(exponential(1), class = phi)),
              iter = 5000, warmup = 2000, chains = 2, cores = 2, stanvars = stanvars)
summary(fithc1)
fithc1$fit

#define post-processing functions for the beta binomial distribution
expose_functions(fithc1, vectorize = TRUE)

log_lik_beta_binomial2 <- function(i, draws) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  trials <- draws$data$vint1[i]
  y <- draws$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

predict_beta_binomial2 <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  trials <- draws$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}

fitted_beta_binomial2 <- function(draws) {
  mu <- draws$dpars$mu
  trials <- draws$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

#sample posterior
post <- fithc1 %>% posterior_samples()
pbC <- post$b_selection
pbR <- post$b_selection + post$'b_treatmentR:selection'
pbT <- post$b_selection + post$'b_treatmentT:selection'


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
    geom_smooth(method = 'glm', se = 0)+
    facet_grid(treatment~.)+
    ylab('proportion on tomato')+
    xlab('selection round')+
    scale_fill_brewer( palette = "Dark2")+
    scale_color_brewer( palette = "Dark2")+
    theme(legend.position="top"))

#combine plots
grid.arrange(grobs = list(choice_sep2, pred_slopes, pred_diff_slopes), widths = c(1,2), layout_matrix = rbind(c(1, 2),c(1, 3)))


#Alternative model: temproal effect removed
fithc2 <- brm(data = choice[!is.na(choice$tomato),], family = betabinom ,
              tomato | vint(total) ~ 1 + treatment + (1|code),  
              prior = c(prior(normal(0,2), class = Intercept),
                        prior(normal(0,2), class = b),
                        prior(cauchy(0,1), class = sd),
                        prior(exponential(1), class = phi)),
              iter = 5000, warmup = 2000, chains = 2, cores = 2, stanvars = stanvars)
summary(fithc2)
fithc2$fit

#Weigh the alternative to the first model
waic(fithc1)
waic(fithc2)
model_weights(fithc1, fithc2) %>% round(digits = 2)

#Sample posterior
post2 <- fithc2 %>% posterior_samples()
pC2 <- logistic(post2$b_Intercept)
pR2 <- logistic(post2$b_Intercept+post2$b_treatmentR)
pT2 <- logistic(post2$b_Intercept+post2$b_treatmentT)

#Plot the estimated tomato preference for every treatment independent of time
g <- data.frame('C' = pC2, 'R'= pR2, 'T' = pT2)
g_ <- gather(g, key = 'treatment', value = 'post')
(pred_pref <- ggplot(data = g_, aes(x = treatment, y = post))+
    geom_boxplot(data = choice[!is.na(choice$tomato),], aes(x=treatment, y=proptomato ))+
    geom_violin(alpha = 0.3, draw_quantiles = c(0.09, 0.5, 0.91), aes(fill = treatment, color = treatment))+
    ggtitle('observed and estimated tomato preference')+
    ylab('tomato preference')+
    scale_color_brewer( palette = "Dark2")+
    scale_fill_brewer( palette = "Dark2")+
    theme(legend.position='none'))



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

#model the expected number of mites on tomato after a day (tomato2) depending on developmental habitat
fithi1 <- brm(data = hi, family = betabinom ,
              tomato2 | vint(total2) ~ 1 + mat + treatment,  
              prior = c(prior(normal(0,2), class = Intercept),
                        prior(normal(0,2), class = b),
                        prior(exponential(1), class = phi)),
              iter = 5000, warmup = 2000, chains = 2, cores = 2, stanvars = stanvars)
summary(fithi1)
fithi1$fit

#sample posterior
posthi <- fithi1 %>% posterior_samples()
hC <- logistic(posthi$b_Intercept)
hT <- logistic(posthi$b_Intercept + posthi$b_matT)

#construct a dataframe with estimated expected tomato2 values (a) and estimated expected differences between both developmental habitats (b)
a_after <- gather(data.frame('T' = hT, 'C' = hC), key = 'treatment', value = 't')
b_after <- data.frame('end' = hT-hC)
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

fithi2 <- brm(data = hi_pre, family = betabinom ,
              tomato | vint(total) ~ 1 + mat,  
              prior = c(prior(normal(0,2), class = Intercept),
                        prior(normal(0,2), class = b),
                        prior(exponential(1), class = phi)),
              iter = 5000, warmup = 2000, chains = 2, cores = 2, stanvars = stanvars)
summary(fithi2)
fithi2$fit
#sample posterior
posthi2 <- fithi2 %>% posterior_samples()
hC2 <- logistic(posthi2$b_Intercept)
hT2 <- logistic(posthi2$b_Intercept + posthi2$b_matT)

#construct a dataframe with estimated expected tomato2 values (a) and estimated expected differences between both developmental habitats (b)
a_pre <- gather(data.frame('T' = hT2, 'C' = hC2), key = 'treatment', value = 't')
b_pre <- data.frame('start' = hT2-hC2)

#estimated proportion choosing tomato according to developmental environment
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
b <- gather(data.frame(cbind(b_pre, b_after)), key = 'time', value = 't')
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

#model females and deutonymphs (reproductive success, nomales12) dependent on tested host at end of experiment
fitLH <- brm(data = LH, family = 'negbinomial',
             nomales12 ~ 1 + patch + (1|trep),
             prior = c(prior(normal(0,4), class = Intercept),
                       prior(normal(0,4), class = b),
                       prior(exponential(1), class = shape),
                       prior(cauchy(0,2), class = sd)),
             iter = 5000, warmup = 2000, chains = 2, cores = 2)

summary(fitLH)
fitLH$fit
#sample posterior
postLH <- fitLH %>% posterior_samples()
sucC <- exp(postLH$b_Intercept + postLH$b_patchC)
sucT <- exp(postLH$b_Intercept + postLH$b_patchT)
sucB <- exp(postLH$b_Intercept)


#plot estimated reproductive success on different hosts at end
a <- gather(data.frame('B' = sucB, 'T' = sucT, 'C' = sucC), key = 'treatment', value = 't')
(fit_after <- ggplot(data = LH, aes(x = patch, y = nomales12))
  +geom_boxplot()+
    expand_limits(y = 0)+
    geom_violin(data = a, draw_quantiles = c(0.09, 0.5, 0.91), aes(x = treatment, y = t, width = 0.5))+
    ggtitle('reproductive success end')+
    xlab('host')+
    ylab('deutonymphs produced in 12 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot estimated differences of reproductive success between different hosts at end
b <- gather(data.frame('B_C'=sucB-sucC, 'B_T'=sucB-sucT, 'T_C'=sucT-sucC), key = 'pairs', value = 'pairwise_differences')
(fit_diff_after <- ggplot(data = b, aes(x = pairs, y = pairwise_differences))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91))+
    expand_limits(y = 0)+
    ggtitle('posterior predicted differences end')+
    xlab('pairs of treatments'))


#model females and deutonymphs (reproductive success, nomales12) dependent on tested host at start of experiment
fitLHpre <- brm(data = LH_pre[LH_pre$day == 12 & LH_pre$cg == 'tomato',], family = 'negbinomial',
             nomales ~ 1 + patch,
             prior = c(prior(normal(0,4), class = Intercept),
                       prior(normal(0,4), class = b),
                       prior(exponential(1), class = shape)),
             iter = 5000, warmup = 2000, chains = 2, cores = 2)

summary(fitLHpre)
fitLHpre$fit
#sample posterior
postLHpre <- fitLHpre %>% posterior_samples()
sucCp <- exp(postLHpre$b_Intercept + postLHpre$b_patchC)
sucTp <- exp(postLHpre$b_Intercept + postLHpre$b_patchT)
sucBp <- exp(postLHpre$b_Intercept)


#plot estimated reproductive success on different hosts at start
a_ <- gather(data.frame('B' = sucBp, 'T' = sucTp, 'C' = sucCp), key = 'treatment', value = 't')
(fit_before <- ggplot(data = LH_pre[LH_pre$day == 12 & LH_pre$cg == 'tomato',], aes(x = patch, y = nomales))+
    geom_boxplot()+
    geom_violin(data = a_, draw_quantiles = c(0.09, .5, 0.91), aes(x = treatment, y = t, width = 0.5, trim = 1))+
    ggtitle('reproductive success start')+
    xlab('host')+
    ylab('deutonymphs produced in 12 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot estimated differences of reproductive success between different hosts at start
b_ <- gather(data.frame('B_C'=sucBp-sucCp, 'B_T'=sucBp-sucTp, 'T_C'=sucTp-sucCp), key = 'pairs', value = 'pairwise_differences')
(fit_diff_before <- ggplot(data = b_, aes(x = pairs, y = pairwise_differences))+
    geom_violin(draw_quantiles = c(0.09, .5, 0.91))+
    expand_limits(y = 0)+
    ggtitle('posterior predicted differences start')+
    xlab('pairs of treatments'))

#combine plots
grid.arrange(fit_before, fit_after, fit_diff_before, fit_diff_after, nrow = 2)

#Appendix 1: data on fertility (day6)

#model fertility at end
fitLH3 <- brm(data = LH, family = 'negbinomial',
             total6 ~ 1 + patch + (1|trep),
             prior = c(prior(normal(0,4), class = Intercept),
                       prior(normal(0,4), class = b),
                       prior(exponential(1), class = shape),
                       prior(cauchy(0,2), class = sd)),
             iter = 5000, warmup = 2000, chains = 2, cores = 2)

summary(fitLH3)
fitLH3$fit
#sample posterior
postLH3 <- fitLH3 %>% posterior_samples()
ferC <- exp(postLH3$b_Intercept + postLH3$b_patchC)
ferT <- exp(postLH3$b_Intercept + postLH3$b_patchT)
ferB <- exp(postLH3$b_Intercept)


#plot estimated reproductive fertility on different hosts at end
a <- gather(data.frame('B' = ferB, 'T' = ferT, 'C' = ferC), key = 'treatment', value = 't')
b <- gather(data.frame('B_C'=ferB-ferC, 'B_T'=ferB-ferT, 'T_C'=ferT-ferC), key = 'pairs', value = 'pairwise_differences')

(fit6_after <- ggplot(data = LH, aes(x = patch, y = total6))
  +geom_boxplot()+
    expand_limits(y = 0)+
    geom_violin(data = a, aes(x = treatment, y = t, width = 0.5))+
    ggtitle('fertility end')+
    xlab('host')+
    ylab('eggs produced in 6 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot differences in fertility at end
(fit6_diff_after <- ggplot(data = b, aes(x = pairs, y = pairwise_differences))+
    geom_violin()+
    expand_limits(y = 0)+
    ggtitle('posterior predicted differences end')+
    xlab('pairs of hosts'))


#model fertility at start
fitLH4 <- brm(data = LH_pre[LH_pre$day == 6 & LH_pre$cg == 'tomato',], family = 'negbinomial',
                d6 ~ 1 + patch,
                prior = c(prior(normal(0,4), class = Intercept),
                          prior(normal(0,4), class = b),
                          prior(exponential(1), class = shape)),
                iter = 5000, warmup = 2000, chains = 2, cores = 2)

summary(fitLH4)
fitLH4$fit
#sample posterior
postLH4 <- fitLH4 %>% posterior_samples()
ferCp <- exp(postLH4$b_Intercept + postLH4$b_patchC)
ferTp <- exp(postLH4$b_Intercept + postLH4$b_patchT)
ferBp <- exp(postLH4$b_Intercept)


#plot estimated reproductive fertility on different hosts at start
a_ <- gather(data.frame('B' = ferBp, 'T' = ferTp, 'C' = ferCp), key = 'treatment', value = 't')
b_ <- gather(data.frame('B_C'=ferBp-ferCp, 'B_T'=ferBp-ferTp, 'T_C'=ferTp-ferCp), key = 'pairs', value = 'pairwise_differences')

(fit6_before <- ggplot(data = LH_pre[LH_pre$day == 6 & LH_pre$cg == 'tomato',], aes(x = patch, y = d6))+
    geom_boxplot()+
    geom_violin(data = a_, aes(x = treatment, y = t, width = 0.5, trim = 1))+
    ggtitle('fertility start')+
    xlab('host')+
    ylab('eggs produced in 6 days')+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=3,show_guide = FALSE))

#plot differences in fertility at start
(fit6_diff_before <- ggplot(data = b_, aes(x = pairs, y = pairwise_differences))+
    geom_violin()+
    expand_limits(y = 0)+
    ggtitle('posterior predicted differences start')+
    xlab('pairs of hosts'))

#combine plots
grid.arrange(fit6_before, fit6_after, fit6_diff_before, fit6_diff_after, nrow = 2)


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
