# LOAD PACKAGES ####
library(lme4) # upload package to use lmer function and glmer function
library(car) #upload package to use Anova function for Analysis of Deviance
library(MuMIn) # to calculate marginal R^2
library(ggplot2) # to produce figures
library(gridExtra) # to combine figures
library(multcomp) # to make post-hoc mean comparisons
library(bbmle) # to use AICc model selection


# OYSTER SURVIVAL AFTER ONE MONTH (JUVENILES) #### 
data = read.csv("juv.abundance_8May2017.csv")
data$cage.number=as.factor(data$cage.number) # convert cage.number to factor to use as random effect

# FAR site analysis of juvenle survival with glmer function
SIN = subset(data, site == "sin") # SIN = FAR site
model.surv.sin = glmer(cbind(live,dead)~ treatment*tile.type + (1|cage.number), data = SIN, family = "binomial", REML = TRUE)
Anova(model.surv.sin, REML = FALSE)

#run simpler model based on Analysis of deviance to facilitate calculating marginal R^2
model.surv.sin.c = glmer(cbind(live,dead)~ tile.type + (1|cage.number), data = SIN, family = "binomial")
r.squaredGLMM(model.surv.sin.c) # to get marginal R2 of the significant factor of tile.type

library(multcomp) # load the  multcomp package to conduct post-hoc comparison of means
summary(glht(model.surv.sin.c, mcp(tile.type="Tukey"))) # compare means within factor "tile.type"

# repeat same process for middle site = ML
ML = subset(data, site == "ml")
model.surv.ml = glmer(cbind(live,dead)~ treatment*tile.type + (1|cage.number), data = ML, family = "binomial")
Anova(model.surv.ml)
model.surv.ml.c = glmer(cbind(live,dead)~ tile.type + (1|cage.number), data = ML, family = "binomial")
r.squaredGLMM(model.surv.ml.c)
summary(glht(model.surv.ml.c, mcp(tile.type="Tukey")))

# repeat same process for CLOSE site = PF
PF = subset(data, site == "pf")
model.surv.pf = glmer(cbind(live,dead)~ treatment*tile.type + (1|cage.number), data = PF, family = "binomial")
Anova(model.surv.pf) # Has singular fit but, It only affects the random effect estimate; the fixed effects are not different if you just run a GLM. Probably not enough differences among some of the cages.
r.squaredGLMM(model.surv.pf.c)
model.surv.pf.c = glmer(cbind(live,dead)~ tile.type + (1|cage.number), data = PF, family = "binomial")
summary(glht(model.surv.pf.c, mcp(tile.type="Tukey")))

# Model selection approach on results only from no culling and no predator cue treatments to evaluate best explanation of spatial variation in juvenile oyster survival
juv.surv.mod.data = read.csv("FigS5_data.csv")
juv.surv.mod.data$cage.number = as.factor(juv.surv.mod.data$cage.number)
model1.jsurv = lmer(surv~1 + (1|cage.number), data = juv.surv.mod.data)
model2.jsurv = lmer(surv~mean.flow + (1|cage.number), data = juv.surv.mod.data)
model3.jsurv = lmer(surv~mean.sal + (1|cage.number), data = juv.surv.mod.data)
model4.jsurv = lmer(surv~mean.temp + (1|cage.number), data = juv.surv.mod.data)
model5.jsurv = lmer(surv~mean.prop +(1|cage.number), data = juv.surv.mod.data)
model6.jsurv = lmer(surv~mean.chl + (1|cage.number), data = juv.surv.mod.data)
AICctab(logLik(model1.jsurv),logLik(model2.jsurv),logLik(model3.jsurv),logLik(model4.jsurv), logLik(model5.jsurv),logLik(model6.jsurv),weights=TRUE,delta=TRUE)
r.squaredGLMM(model5.jsurv) # to get R^2

# OYSTER GROWTH AFTER ONE MONTH (JUVENILES) ####
juv.data=read.csv("juvenile_growth.csv")
juv.data$cage.number=factor(juv.data$cage.number) #convert this to a factor to use as random effect
model.juv.growth.all = lmer(growth~ site*treatment*tile.type + (1|cage.number), data = juv.data)
Anova(model.juv.growth.all)
#re-run with only significant factors to facilitate posthoc mean comparisons
model.juv.growth.all.e = lmer(growth~ site + treatment + (1|cage.number), data = juv.data)
summary(glht(model.juv.growth.all.e, mcp(treatment="Tukey")))
# re-run with single-factor models that allow calculation of marginal R^2 for each significant main effect
model.juv.growth.all.h = lmer(growth~ site  + (1|cage.number), data = juv.data)
r.squaredGLMM(model.juv.growth.all.h)
model.juv.growth.all.i = lmer(growth~ treatment  + (1|cage.number), data = juv.data)
r.squaredGLMM(model.juv.growth.all.i)

# Model selection approach on results only from no culling and no predator cue treatments to evaluate best explanation of spatial variation in juvenile oyster survival
juv_growth = read.csv("FigS6_data.csv")
juv_growth$cage.number = as.factor(juv_growth$cage.number)
model1 = lmer(growth~1 + (1|cage.number), data = juv_growth)
model2 = lmer(growth~mean.flow + (1|cage.number), data = juv_growth)
model3 = lmer(growth~mean.sal + (1|cage.number), data = juv_growth)
model4 = lmer(growth~mean.temp + (1|cage.number), data = juv_growth)
model5 = lmer(growth~mean.prop +(1|cage.number), data = juv_growth)
model6 = lmer(growth~mean.chl + (1|cage.number), data = juv_growth)
AICctab(logLik(model1),logLik(model2),logLik(model3),logLik(model4), logLik(model5),logLik(model6),weights=TRUE,delta=TRUE)
r.squaredGLMM(model2) # to get pseudo R^2

# OYSTER CONDITION INDEX AFTER ONE MONTH (JUVENILES) ####
library(dplyr) # to allow data processing step
juv_ci_data = read.csv("juv_ci_data_28July2017.csv") 
juv_ci_data$cage=as.factor(juv_ci_data$cage) #convert cage.number to a factor to use as random effect
juv_ci_data2 = subset(juv_ci_data, tissue.mass >0) #remova all data that lacks tissue mass
ci_avg <- juv_ci_data2 %>%
  group_by(site,treatment,tile.treatment,cage,tile) %>%
  summarize(mean_ci = mean(ci, na.rm=TRUE))
model.juv.ci.all = lmer(mean_ci ~ site*treatment*tile.treatment + (1|cage), data = ci_avg)
Anova(model.juv.ci.all)

# FIGURE 2 ILLUSTRATION OF ONE MONTH RESULTS (JUVENILES) ####
#rename and re-order sites
data$site=factor(data$site, levels=c("sin", "ml", "pf"))
levels(data$site)[levels(data$site)=="sin"] <- "Far"
levels(data$site)[levels(data$site)=="ml"] <- "Mid"
levels(data$site)[levels(data$site)=="pf"] <- "Close"
#rename and re-order cue treatments
data$treatment=factor(data$treatment, levels=c("np", "mc",  "cc"))
levels(data$treatment)[levels(data$treatment)=="np"] <- "1a"
levels(data$treatment)[levels(data$treatment)=="mc"] <- "2b"
levels(data$treatment)[levels(data$treatment)=="cc"] <- "3c"
#rename and re-order CE stage treatments
data$tile.type=factor(data$tile.type, levels=c("no.cull", "juv.cull"))
levels(data$tile.type)[levels(data$tile.type)=="no.cull"] <- "1None"
levels(data$tile.type)[levels(data$tile.type)=="juv.cull"] <- "2Juvenile"

#pd <- position_dodge(0.50)
# produce Figure 2a-c. This was further modified in Adobe 
juv_abund_fig <-ggplot(data=data,
                       mapping = aes(x=treatment, 
                                     y = surv)) +
  facet_grid(rows = vars(site),cols=vars(tile.type)) +
  geom_boxplot(aes(fill=treatment),outlier.shape = NA,alpha=0.5)+
  scale_fill_manual(values=c("black", "dark grey", "light grey"), 
                    name=NULL, labels=c("1np", "2mc", "3cc"))+
  geom_jitter(color='black',width=0.2,size=0.4,height=0.05) +
  ylab("Survival after one month (proportional)")+
  xlab("Oyster life-stage exposed to consumptive effect")+
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()) +
  theme(legend.position = "none")

#produce Figure 2d-f
juv.data$site=factor(juv.data$site, levels=c("SIN", "ML", "PF"))
levels(juv.data$site)[levels(juv.data$site)=="1SIN"] <- "Far"
levels(juv.data$site)[levels(juv.data$site)=="2ML"] <- "Mid"
levels(juv.data$site)[levels(juv.data$site)=="3PF"] <- "Close"

juv.data$treatment=factor(juv.data$treatment, levels=c("NP", "MC",  "CC"))
levels(juv.data$treatment)[levels(juv.data$treatment)=="NP"] <- "1NP"
levels(juv.data$treatment)[levels(juv.data$treatment)=="MC"] <- "2MC"
levels(juv.data$treatment)[levels(juv.data$treatment)=="CC"] <- "3CC"

juv.data$tile.type=factor(juv.data$tile.type, levels=c("No.cull", "Juv.cull"))
levels(juv.data$tile.type)[levels(juv.data$tile.type)=="No.cull"] <- "1No.cull"
levels(juv.data$tile.type)[levels(juv.data$tile.type)=="Juv.cull"] <- "2Juv.cull"

juv_growth_fig <-ggplot(data=juv.data,
                        mapping = aes(x=treatment, 
                                      y = growth)) +
  facet_grid(rows = vars(site),cols=vars(tile.type)) +
  geom_boxplot(aes(fill=treatment),outlier.shape = NA,alpha=0.5)+
  scale_fill_manual(values=c("black", "dark grey", "light grey"), 
                    name=NULL, labels=c("1NP", "2MC", "3CC"))+
  geom_jitter(color='black',width=0.2,size=0.4,height=0.05) +
  ylab("Growth increment (mm) after one month")+
  xlab("Oyster life-stage exposed to consumptive effect")+
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()) +
  theme(legend.position = "none")

# combine survival and growth plots into Figure 2
quartz(width=8.75,height=5.25)
grid.arrange(juv_abund_fig, juv_growth_fig, ncol = 2)
quartz.save('Fig2.pdf',type="pdf", dpi = 600)

# OYSTER SURVIVAL AFTER FOUR MONTHS (ADULTS) ####

adult_abundance = read.csv("adult_abundance_29Jan2020.csv")

# subset data for just site FAR from freshwater input = SIN
sin_adult_all = subset(adult_abundance, site == "SIN")
gC = glmerControl(check.conv.grad = .makeCC("warning",tol=1e-2)) # relax convergence tolerance somewhat
model.adult.surv.sin.all = glmer(cbind(live,dead)~ treatment*tile.type + (1|cage), data = sin_adult_all, family = "binomial",control=gC)
Anova(model.adult.surv.sin.all) # FAILED to converge
model.adult.surv.sin.all.c = glmer(cbind(live,dead)~ tile.type + (1|cage), data = sin_adult_all, family = "binomial") # run simpler model to support post-hoc comparison of means
summary(glht(model.adult.surv.sin.all.c, mcp(tile.type="Tukey"))) # comprae means for CE factor
r.squaredGLMM(model.adult.surv.sin.all.c) # marginal R^2

# Repeat analyses for MID site = ML
ml_adult = subset(adult_abundance, site == "ML")
model.adult.surv.ml = glmer(cbind(live,dead)~ treatment*tile.type + (1|cage), data = ml_adult, family = "binomial",control=gC)
Anova(model.adult.surv.ml) # significant interaction but model FAILED to converge
# To facilitate targeted mean comparisions, we created new models for each of the four CE or culling groups and then compare predator cue means within each group
ml_adult_nocull = subset(ml_adult, tile.type == "No.cull") # CE or culling group #1
no.cull.model = glmer(cbind(live,dead)~ treatment + (1|cage), data = ml_adult_nocull, family = "binomial")
summary(glht(no.cull.model, mcp(treatment="Tukey")))
ml_adult_juvcull = subset(ml_adult, tile.type == "Juv.cull") # CE or culling group #2
juv.cull.model = glmer(cbind(live,dead)~ treatment + (1|cage), data = ml_adult_juvcull, family = "binomial")
summary(glht(juv.cull.model, mcp(treatment="Tukey")))
ml_adult_adultcull = subset(ml_adult, tile.type == "Adult.cull") # CE or culling group #3
adult.cull.model = glmer(cbind(live,dead)~ treatment + (1|cage), data = ml_adult_adultcull, family = "binomial")
summary(glht(adult.cull.model, mcp(treatment="Tukey")))
ml_adult_juvadultcull = subset(ml_adult, tile.type == "Juv.Adult.cull") # CE or culling group #4
juv.adult.cull.model = glmer(cbind(live,dead)~ treatment + (1|cage), data = ml_adult_juvadultcull, family = "binomial")
summary(glht(juv.adult.cull.model, mcp(treatment="Tukey")))
# Becuse all mean comparisions within each CE group were not significant, we focused on main effect of CE
# So run a model with just CE or culling treatment factor () averaged across predator cue treatments)
model.adult.surv.ml.c = glmer(cbind(live,dead)~ tile.type + (1|cage), data = ml_adult, family = "binomial")
summary(glht(model.adult.surv.ml.c, mcp(tile.type="Tukey")))
r.squaredGLMM(model.adult.surv.ml.c) # calculate marginal R^2

# Repeat analyses for CLOSE site = PF
pf_adult = subset(adult_abundance, site == "PF")
model.adult.surv.pf = glmer(cbind(live,dead)~ treatment*tile.type + (1|cage), data = pf_adult, family = "binomial")
model.adult.surv.pf2 = glm(cbind(live,dead)~ treatment*tile.type, data = pf_adult, family = "binomial")
Anova(model.adult.surv.pf) # Has singular fit but, It only affects the random effect estimate; the fixed effects are not different if you just run a GLM. Probably not enough differences among some of the cages.
model.adult.surv.pf.c = glmer(cbind(live,dead)~ tile.type + (1|cage), data = pf_adult, family = "binomial")
Anova(model.adult.surv.pf.c)
summary(glht(model.adult.surv.pf.c, mcp(tile.type="Tukey")))

# Model selection on results from no predator cue and no CE treatments...this is of survival/duration
adult.surv.mod.data = read.csv("FigS7_data.csv")
adult.surv.mod.data$cage = as.factor(adult.surv.mod.data$cage)
model1.asurv = lmer(daily.surv~1 + (1|cage), data = adult.surv.mod.data)
model2.asurv = lmer(daily.surv~flow + (1|cage), data = adult.surv.mod.data)
model3.asurv = lmer(daily.surv~salinity + (1|cage), data = adult.surv.mod.data)
model4.asurv = lmer(daily.surv~temp + (1|cage), data = adult.surv.mod.data)
model5.asurv = lmer(daily.surv~exposure +(1|cage), data = adult.surv.mod.data)
model6.asurv = lmer(daily.surv~chl + (1|cage), data = adult.surv.mod.data)
model7.asurv = lmer(daily.surv~waterht + (1|cage), data = adult.surv.mod.data)
AICctab(logLik(model1.asurv),logLik(model2.asurv),logLik(model3.asurv),logLik(model4.asurv), logLik(model5.asurv),logLik(model6.asurv,logLik(model7.asurv)),weights=TRUE,delta=TRUE)

# OYSTER GROWTH AFTER FOUR MONTHS (ADULTS) ####
adult_growth = read.csv("adult_growth_29Jan2020.csv")
model.adult.growth.all = lmer(trans.daily.growth~ site*treatment*tile.type + (1|cage), data = adult_growth)
Anova(model.adult.growth.all)
model.adult.growth.all.d = lmer(trans.daily.growth~ site + tile.type + treatment + (1|cage), data = adult_growth)
Anova(model.adult.growth.all.d) # run simpler model to facilitate post-hoc mean comparisions within each factor
summary(glht(model.adult.growth.all.d, mcp(site="Tukey")))
summary(glht(model.adult.growth.all.d, mcp(treatment="Tukey")))
summary(glht(model.adult.growth.all.d, mcp(tile.type="Tukey")))
# run single factor models to calculat marginal R^2 of each significant main effect
model.adult.growth.all.h = lmer(trans.daily.growth~ site + (1|cage), data = adult_growth)
r.squaredGLMM(model.adult.growth.final.h) # just site
model.adult.growth.all.i = lmer(trans.daily.growth~ treatment + (1|cage), data = adult_growth)
r.squaredGLMM(model.adult.growth.final.i) # just cue treatment
model.adult.growth.all.j = lmer(trans.daily.growth~tile.type + (1|cage), data = adult_growth)
r.squaredGLMM(model.adult.growth.final.j) # just culling treatment

# OYSTER CONDITION INDEX AFTER FOUR MONTHS (ADULTS) ####
adult_ci_data = read.csv("adult_ci_30Jan2020.csv")
adult_ci_data$cage = as.factor(adult_ci_data$cage)
model.adult.ci.all = lmer(daily.ci ~ site*treatment*tile.type + (1|cage), data = adult_ci_data)
Anova(model.adult.ci.all)

# FIGURE 3 ILLUSTRATION OF FOUR MONTH RESULTS (ADULTS) ####
adult_abundance$site=factor(adult_abundance$site, levels=c("SIN", "ML", "PF"))
levels(adult_abundance$site)[levels(adult_abundance$site)=="SIN"] <- "1Far"
levels(adult_abundance$site)[levels(adult_abundance$site)=="ML"] <- "2Mid"
levels(adult_abundance$site)[levels(adult_abundance$site)=="PF"] <- "3Close"

adult_abundance$treatment=factor(adult_abundance$treatment, levels=c("NP", "MC", "CC", "MC+CC"))
levels(adult_abundance$treatment)[levels(adult_abundance$treatment)=="NP"] <- "1Control"
levels(adult_abundance$treatment)[levels(adult_abundance$treatment)=="MC"] <- "2Crab"
levels(adult_abundance$treatment)[levels(adult_abundance$treatment)=="CC"] <- "3Conch"
levels(adult_abundance$treatment)[levels(adult_abundance$treatment)=="MC+CC"] <- "4MP"

adult_abundance$tile.type=factor(adult_abundance$tile.type, levels=c("No.cull", "Juv.cull", "Adult.cull", "Juv.Adult.cull"))
levels(adult_abundance$tile.type)[levels(adult_abundance$tile.type)=="No.cull"] <- "1None"
levels(adult_abundance$tile.type)[levels(adult_abundance$tile.type)=="Juv.cull"] <- "2Juvenile"
levels(adult_abundance$tile.type)[levels(adult_abundance$tile.type)=="Adult.cull"] <- "3Adult"
levels(adult_abundance$tile.type)[levels(adult_abundance$tile.type)=="Juv.Adult.cull"] <- "4both"

# Plot of Adult survival
adult_abund_fig <-ggplot(data=adult_abundance,
                         mapping = aes(x=treatment, 
                                       y = survival)) +
  facet_grid(rows = vars(site),cols=vars(tile.type)) +
  geom_boxplot(aes(fill=treatment),outlier.shape = NA,alpha=0.5)+
  scale_fill_manual(values=c("black", "dark grey", "light grey", "white"), 
                    name=NULL, labels=c("1Control", "2Crab", "3Conch", "4MP"))+
  geom_jitter(color='black',width=0.2,size=0.4,height=0.0005) +
  ylab("Survival after four months (Proportion)") +
  xlab("Oyster life-stage exposed to consumptive effect")+
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()) +
  theme(legend.position = "none")


adult_growth$site=factor(adult_growth$site, levels=c("SIN", "ML", "PF"))
levels(adult_growth$site)[levels(adult_growth$site)=="SIN"] <- "1Far"
levels(adult_growth$site)[levels(adult_growth$site)=="ML"] <- "2Mid"
levels(adult_growth$site)[levels(adult_growth$site)=="PF"] <- "3Close"

adult_growth$treatment=factor(adult_growth$treatment, levels=c("NP", "MC", "CC", "MC+CC"))
levels(adult_growth$treatment)[levels(adult_growth$treatment)=="NP"] <- "1NP"
levels(adult_growth$treatment)[levels(adult_growth$treatment)=="MC"] <- "2MC"
levels(adult_growth$treatment)[levels(adult_growth$treatment)=="CC"] <- "3CC"
levels(adult_growth$treatment)[levels(adult_growth$treatment)=="MC+CC"] <- "4MP"

adult_growth$tile.type=factor(adult_growth$tile.type, levels=c("No.cull", "Juv.cull", "Adult.cull", "Juv.Adult.cull"))
levels(adult_growth$tile.type)[levels(adult_growth$tile.type)=="No.cull"] <- "1None"
levels(adult_growth$tile.type)[levels(adult_growth$tile.type)=="Juv.cull"] <- "2Juv"
levels(adult_growth$tile.type)[levels(adult_growth$tile.type)=="Adult.cull"] <- "3Adult"
levels(adult_growth$tile.type)[levels(adult_growth$tile.type)=="Juv.Adult.cull"] <- "4Both"

adult_growth_fig <-ggplot(data=adult_growth,
                          mapping = aes(x=treatment, 
                                        y = daily.growth)) +
  facet_grid(rows = vars(site),cols=vars(tile.type)) +
  geom_boxplot(aes(fill=treatment),outlier.shape = NA,alpha=0.5)+
  scale_fill_manual(values=c("black", "dark grey", "light grey", "white"), 
                    name=NULL, labels=c("1NP", "2MC", "3CC", "4Both"))+
  geom_jitter(color='black',width=0.2,size=0.4,height=0.05) +
  ylab("Growth increment (mm/d) over four months") +
  xlab("Oyster life-stage exposed to consumptive effect")+
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()) +
  theme(legend.position = "none")

quartz(width=8.75,height=5.25)
grid.arrange(adult_abund_fig, adult_growth_fig, ncol = 2)
quartz.save('Fig3.pdf',type="pdf", dpi = 600)

# OYSTER RECRUITMENT AFTER FOUR MONTHS ####
recruit = read.csv("recruitment_29Jan2020.csv")
recruit = subset(recruit_data, harvest != "Juvenile")
recruit$cage.number=factor(recruit$cage.number) 
model.recruit.all = lmer(std.recruits~ site*treatment*tile.treatment + (1|cage.number), data = recruit)
Anova(model.recruit.all)
# Given output, will now run simpler models to support targeted post-hoc comparison of means
  # For each site, look at tile treatment comparisons
  # For all three sites together, look at treatment comparisons

sin = subset(recruit, site == "SIN") # The FAR Site
model.recruit.sin = lmer(std.recruits~ tile.treatment + (1|cage.number), data = sin)
summary(glht(model.recruit.sin, mcp(tile.treatment="Tukey")))
ml = subset(recruit, site == "ML")# The MID Site
model.recruit.ml = lmer(std.recruits~ tile.treatment + (1|cage.number), data = ml)
summary(glht(model.recruit.ml, mcp(tile.treatment="Tukey")))
pf = subset(recruit, site == "PF") # The CLOSE site
model.recruit.pf = lmer(std.recruits~ tile.treatment + (1|cage.number), data = pf)
summary(glht(model.recruit.pf, mcp(tile.treatment="Tukey")))
model.recruit.all.i = lmer(std.recruits~ treatment  + (1|cage.number), data = recruit)
summary(glht(model.recruit.all.i, mcp(treatment="Tukey")))
r.squaredGLMM(model.recruit.all.c) # site*culling + cue
model.recruit.site.cull = lmer(std.recruits~ site * tile.treatment  + (1|cage.number), data = recruit)
r.squaredGLMM(model.recruit.site.cull) # just site
r.squaredGLMM(model.recruit.all.i) # just cue treatment

#FIGURE 4 ILLUSTRATION OF RECRUITMENT RESULTS AFTER FOUR MONTHS####
recruit$site=factor(recruit$site, levels=c("SIN", "ML", "PF"))
levels(recruit$site)[levels(recruit$site)=="SIN"] <- "1Far"
levels(recruit$site)[levels(recruit$site)=="ML"] <- "2Mid"
levels(recruit$site)[levels(recruit$site)=="PF"] <- "3Close"

recruit$treatment=factor(recruit$treatment, levels=c("NP", "MC", "CC", "MC+CC"))
levels(recruit$treatment)[levels(recruit$treatment)=="NP"] <- "1NP"
levels(recruit$treatment)[levels(recruit$treatment)=="MC"] <- "2MC"
levels(recruit$treatment)[levels(recruit$treatment)=="CC"] <- "3CC"
levels(recruit$treatment)[levels(recruit$treatment)=="MC+CC"] <- "4MP"

recruit$tile.treatment=factor(recruit$tile.treatment, levels=c("No.cull", "Juv.cull", "Adult.cull", "Juv.Adult.cull"))
levels(recruit$tile.treatment)[levels(recruit$tile.treatment)=="No.cull"] <- "1None"
levels(recruit$tile.treatment)[levels(recruit$tile.treatment)=="Juv.cull"] <- "2Juv"
levels(recruit$tile.treatment)[levels(recruit$tile.treatment)=="Adult.cull"] <- "3Adult"
levels(recruit$tile.treatment)[levels(recruit$tile.treatment)=="Juv.Adult.cull"] <- "4Both"

quartz(width=8.75,height=5.25)
ggplot(data=recruit,
                     mapping = aes(x=treatment, 
                                   y = std.recruits)) +
  facet_grid(rows = vars(site),cols=vars(tile.treatment)) +
  geom_boxplot(aes(fill=treatment),outlier.shape = NA,alpha=0.5)+
  scale_fill_manual(values=c("black", "dark grey", "light grey", "white"), 
                    name=NULL, labels=c("1NP", "2MC", "3CC", "4MP"))+
  geom_jitter(color='black',width=0.2,size=0.4,height=0.05) +
  ylab("Recruitment after four months (tile/d)")+
  xlab("Oyster life-stage exposed to consumptive effect")+
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()) +
  theme(legend.position = "none")
quartz.save('Fig4.pdf',type="pdf", dpi = 600)
