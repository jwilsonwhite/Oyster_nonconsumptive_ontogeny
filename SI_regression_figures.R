# Figures for SI of Kimbro et al. NCE ontogeny paper
# Figures showing results of AIC selection for physical variables predicting oyster growth & survival

library(ggplot2)
library(lme4)
library(bbmle)
#library(MCMCglmm)
#library(RColorBrewer)
library(colorspace)
library(MuMIn)


#----------------------------------------------------------------------
# Fig S5: juvenile oyster survival
D = read.csv('FigS5_data.csv',header=TRUE)
D$surv.asin <- asin(sqrt(D$surv))
D$cage.number = as.factor(D$cage.number)

m0 = lmer(surv.asin~1+(1|cage.number),data=D)
m1 = lmer(surv.asin~mean.flow+(1|cage.number),data=D)
m2 = lmer(surv.asin~mean.sal+(1|cage.number),data=D)
m3 = lmer(surv.asin~mean.temp+(1|cage.number),data=D)
m4 = lmer(surv.asin~mean.prop+(1|cage.number),data=D)
m5 = lmer(surv.asin~mean.chl+(1|cage.number),data=D)
AICctab(m0,m1,m2,m3,m4,m5,weights=TRUE)
# m4 is best (mean.prop)

# Get predicted line
Dtmp = data.frame(mean.prop=seq(from=min(D$mean.prop),to=max(D$mean.prop),length.out=100))
c = fixef(m4)
Dtmp$surv.asin = c[1]+Dtmp$mean.prop*c[2]
Dtmp$surv = (sin(Dtmp$surv.asin))^2


# Create colorscale (using colorspace)
reds2 = rainbow_hcl(n=9)

quartz(width=4,height=4)
ggplot(data=D,aes(x=mean.prop,y=surv))+
  geom_jitter(aes(color=cage.number),width=0.002,height=0.01,shape=21,stroke=1)+
  scale_color_manual(values=reds2)+
  geom_line(data=Dtmp,aes(x=Dtmp$mean.prop,y=Dtmp$surv))+
  xlab('Proportion of day exposed at low tide')+
  ylab('Oyster survival after one month')+
  theme_bw()+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title = element_text(size=14))
quartz.save('FigS5.pdf',type="pdf")

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Fig S6: juvenile  oyster growth
D = read.csv('FigS6_data.csv',header=TRUE)
D$cage.number = as.factor(D$cage.number)

m0 = lmer(growth~1+(1|cage.number),data=D)
m1 = lmer(growth~mean.flow+(1|cage.number),data=D)
m2 = lmer(growth~mean.sal+(1|cage.number),data=D)
m3 = lmer(growth~mean.temp+(1|cage.number),data=D)
m4 = lmer(growth~mean.prop+(1|cage.number),data=D)
m5 = lmer(growth~mean.chl+(1|cage.number),data=D)
AICctab(m0,m1,m2,m3,m4,m5,weights=TRUE)

# Get predicted line
Dtmp = data.frame(mean.flow=seq(from=min(D$mean.flow),to=max(D$mean.flow),length.out=100))
c = fixef(m1)
Dtmp$growth <- c[1]+Dtmp$mean.flow*c[2]

# Model 4 is the best based on DIC, but no clear winner. None have signif regression coefficients.

# Create colorscale (using colorspace)
reds2 = rainbow_hcl(n=9)

quartz(width=4,height=4)
ggplot(data=D,aes(x=mean.flow,y=growth))+
  geom_jitter(aes(color=cage.number),width=0.2,height=0,shape=21,stroke=1)+
  scale_color_manual(values=reds2)+
  #geom_ribbon(data=Dtmp,aes(x=Dtmp$mean.flow,ymin=Pl[,2],ymax=Pl[,3]),alpha=0.3)+
  geom_line(data=Dtmp,aes(x=mean.flow,y=growth))+
  xlab('Bulk water flow (µg/d dissolution)')+
  ylab('Oyster growth after one month (mm)')+
  theme_bw()+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title = element_text(size=14))
quartz.save('FigS6.pdf',type="pdf")

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Fig S7: adult oyster survival
D = read.csv('FigS7_data.csv',header=TRUE)
D$surv.asin <- asin(sqrt(D$survival))
D$cage = as.factor(D$cage)

m0 = lmer(surv.asin~1+(1|cage),data=D)
m1 = lmer(surv.asin~flow+(1|cage),data=D)
m2 = lmer(surv.asin~salinity+(1|cage),data=D)
m3 = lmer(surv.asin~temp+(1|cage),data=D)
m4 = lmer(surv.asin~exposure+(1|cage),data=D)
AICctab(m0,m1,m2,m3,m4,weights=TRUE)
# m4 is best (exposure)

# Get predicted line
Dtmp = data.frame(exposure=seq(from=min(D$exposure),to=max(D$exposure),length.out=100))
c = fixef(m4)
Dtmp$surv.asin = c[1]+Dtmp$exposure*c[2]
Dtmp$surv = (sin(Dtmp$surv.asin))^2
c = fixef(m0)
Dtmp$null = Dtmp$surv
Dtmp$null = (sin(c[1]))^2

# Create colorscale (using colorspace)
reds2 = rainbow_hcl(n=9)

quartz(width=4,height=4)
ggplot(data=D,aes(x=exposure,y=survival))+
  geom_jitter(aes(color=cage),width=0.002,height=0.0001,shape=21,stroke=1)+
  scale_color_manual(values=reds2)+
  #geom_ribbon(data=Dtmp,aes(x=Dtmp$mean.chl,ymin=Pl[,2],ymax=Pl[,3]),alpha=0.3)+
  geom_line(data=Dtmp,aes(x=exposure,y=surv),linetype=2)+
  geom_line(data=Dtmp,aes(x=exposure,y=null))+
  xlab('Proportion of day exposed at low tide')+
  ylab('Oyster survival after four months')+
  theme_bw()+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title = element_text(size=14))
quartz.save('FigS7.pdf',type="pdf")

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Fig S8: adult oyster growth
D = read.csv('FigS8_data.csv',header=TRUE)
D$cage = as.factor(D$cage)

m0 = lmer(daily.growth~1+(1|cage),data=D)
m1 = lmer(daily.growth~flow+(1|cage),data=D)
m2 = lmer(daily.growth~salinity+(1|cage),data=D)
m3 = lmer(daily.growth~temp+(1|cage),data=D)
m4 = lmer(daily.growth~exposure+(1|cage),data=D)
m5 = lmer(daily.growth~chl+(1|cage),data=D)
AICctab(m0,m1,m2,m3,m4,m5,weights=TRUE)

# Get predicted line
Dtmp = data.frame(exposure=seq(from=min(D$exposure),to=max(D$exposure),length.out=100))
c = fixef(m4)
Dtmp$growth <- c[1]+Dtmp$exposure*c[2]
Dtmp$flow=seq(from=min(D$flow),to=max(D$flow),length.out=100)
c = fixef(m1)
Dtmp$growth.flow <- c[1]+Dtmp$flow*c[2]
c = fixef(m0)
Dtmp$growth.null = Dtmp$growth
Dtmp$growth.null = c[1]


# Create colorscale (using colorspace)
reds2 = rainbow_hcl(n=9)

quartz(width=4,height=4)
ggplot(data=D,aes(x=exposure,y=daily.growth))+
  geom_jitter(aes(color=cage),width=0.002,height=0.01,shape=21,stroke=1)+
  scale_color_manual(values=reds2)+
  #geom_ribbon(data=Dtmp,aes(x=Dtmp$mean.flow,ymin=Pl[,2],ymax=Pl[,3]),alpha=0.3)+
  geom_line(data=Dtmp,aes(x=exposure,y=growth),linetype=2)+
  geom_line(data=Dtmp,aes(x=exposure,y=growth.null))+
  xlab('Proportion of day exposed at low tide')+
  ylab('Oyster growth over five months (mm/d)')+
  theme_bw()+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title = element_text(size=14))
quartz.save('FigS8.pdf',type="pdf")

# flow instead
ggplot(data=D,aes(x=flow,y=daily.growth))+
  geom_jitter(aes(color=cage),width=0.002,height=0.01,shape=21,stroke=1)+
  scale_color_manual(values=reds2)+
  #geom_ribbon(data=Dtmp,aes(x=Dtmp$mean.flow,ymin=Pl[,2],ymax=Pl[,3]),alpha=0.3)+
  geom_line(data=Dtmp,aes(x=flow,y=growth.flow),linetype=2)
  xlab('Bulk water flow (µg/d dissolution)')+
  ylab('Oyster growth over five months (mm/d)')+
  theme_bw()+
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title = element_text(size=14))
quartz.save('FigS8.pdf',type="pdf")

#----------------------------------------------------------------------
