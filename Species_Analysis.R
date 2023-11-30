
library(terra)
library(sf)
library(lubridate)
library(sp)
library(magrittr)
library(ggplot2)
library(dplyr)  
library(reshape2)
library(patchwork)
library(sjPlot)
library(vegan)
library(performance)

### -------------load counting data

df<-read.csv("C:/OnboardSeabirdCounts/DrakePassage/nearestECA57-59c.csv")

summary(as.factor(df$Species))

df$Species[df$Species=="AntTern"]<-"ANT"
df$Species[df$Species=="PachBelch"]<-"SBPR"
df$Species[df$Species=="PachDesol"]<-"ANPR"
df$Species[df$Species=="FurSeal"]<-"AFS"
df$Species[df$Species=="Dolphin"]<-"PED"
df$Species[df$Species=="Whale"]<-"NIDWH"
df$Species[df$Species=="Tern"]<-"SAMT"
df$Species[df$Species=="Prion"]<-"NIDPR"
df$Species[df$Species=="PilotWhale"]<-"PWH"

spf<-plyr::ddply(df, c("Species"), summarise,
                 totcount=length(Year))




ggplot(spf,aes(reorder(Species,+totcount),totcount))+
  geom_bar(stat="identity")+theme_bw()+
  xlab("Taxa")+ylab("Total count")+coord_flip()





head(df)
head(trsm)

envi<-data.frame(trsm[1:5],trsm[10:11],trsm[15:18],trsm[21:22])
head(envi)

envim<-plyr::ddply(envi, c("Transect","Year"), summarise,
                   lon=mean(lon),lat=mean(lat),
                   chl=mean(na.omit(chl)),
                   sst=mean(na.omit(sst)),
                   ws=mean(na.omit(windspeed)),
                   mld=mean(na.omit(mld)),
                   swv=mean(na.omit(swv)),
                   zoo=mean(na.omit(zooc)))

envim$TRY<-paste(envim$Transect,envim$Year)

tail(envim)

enviS <- subset(envim,TRY!="1 2020" &TRY!="1 2021"&TRY!="2 2021"&TRY!="12 2021") # transects in the magellanic and bransfield strait were suppressed


dfnn<-merge(df[1:13],enviS[1:10], by=c("Transect","Year"))  # for nearest neighbor analysis later


df2<-data.frame(df[1:7],pres=c(1))


dfm<-varcast<-dcast(data=df2,formula=Transect+Year+Month+Day+Hour+Min~Species,value.var = "pres",sum)

dfmm<-merge(dfm,enviS, by=c("Transect","Year"))


head(dfmm)




dfom<-na.omit(dfmm)

head(dfom)

ggplot(dfom,aes(Prion,BBA))+geom_point()

### ----- environmental and species relations (interspecific associations)----
### canonical correspondence analysis

sp<-dfom[7:38]
hab<-data.frame(scale(dfom[41:46]))

head(sp)
head(hab)


rd1<-rda(sp,hab,scale=T)

rd1
rd1$CCA$tot.chi

ordiplot(rd1, choices = c(1, 2), display = c("sp", "wa", "cn"),
     scaling = "symmetric",  correlation = F, hill =F)

text(rd1, display = "species", choices = c(1, 2),
     scaling = "symmetric", head.arrow = 0.05,arrow.mul=1,
     axis.bp = FALSE, correlation = FALSE, hill = FALSE)

anova(rd1)
rd1$tot.chi

rd1$CCA$eig
rd1$CA$eig


sps<-data.frame(scores(rd1, display = 'species',scaling = "symmetric",  correlation = F, hill =F))
scr<-data.frame(scores(rd1, display = 'sites',scaling = "symmetric",  correlation = F, hill =F))

vars<-data.frame(scores(rd1, display = 'bp',scaling = "symmetric",  correlation = F, hill =F))

summary(sps)
sps$sp<-row.names(sps)    


summary(sps)

sps$RDA1<-sps$RDA1/2.6372
sps$RDA2<-sps$RDA2/1.246

vars$vars<-row.names(vars) 


summary(vars)

vars$RDA1<-vars$RDA1/0.3253
vars$RDA2<-vars$RDA2/0.2695

summary(scr)

scr$RDA1<-scr$RDA1/2.14
scr$RDA2<-scr$RDA2/2.662

scores<-data.frame(dfom[1:2],scr)



head(scores)
ggplot()+
  geom_point(data=scores,aes(x=RDA1,y=RDA2,shape=as.factor(Year)),alpha=0.5)+
  geom_hline(yintercept = 0,linewidth=1,linetype="dashed",colour="grey50")+
  geom_vline(xintercept = 0,linewidth=1,linetype="dashed",colour="grey50")+
  geom_text(data=sps,aes(x=RDA1,y=RDA2,label=sp),colour="red2",size=3)+
  geom_text(data=vars,aes(x=RDA1,y=RDA2,label=vars),colour="blue2",size=5)+
  theme_bw()+xlab("RDA 1 (58.6%)")+ylab("RDA 2 (31.9%)")


### --------conspecific associations---------
### --- test the variables more relevant to axis 1 and 2
### ------ CHL, SST and WS 
head(dfnn)

(ggplot(na.omit(dfnn),aes(chl,Consp))+stat_smooth(method="glm", method.args = list(family = "binomial"))+
  theme_bw()+xlab("mg/m3")+ylab("Conspecific probability")+ggtitle(label="a. Chlrophyll-a concentration"))+

  (ggplot(na.omit(dfnn),aes(zoo,Consp))+stat_smooth(method="glm", method.args = list(family = "binomial"))+
     theme_bw()+xlab("g/m2")+ylab("Conspecific probability")+ggtitle(label="b. Mass content of zooplankton"))+

(ggplot(na.omit(dfnn),aes(sst,Consp))+stat_smooth(method="glm", method.args = list(family = "binomial"))+
  theme_bw()+xlab("°C")+ylab("Conspecific probability")+ggtitle(label="c. Sea surface temperature"))+


(ggplot(na.omit(dfnn),aes(mld,Consp))+stat_smooth(method="glm", method.args = list(family = "binomial"))+
  theme_bw()+xlab("meters")+ylab("Conspecific probability")+ggtitle(label="d. Thickness of the mixed layer"))+

(ggplot(na.omit(dfnn),aes(ws,Consp))+stat_smooth(method="glm", method.args = list(family = "binomial"))+
  theme_bw()+xlab("m/s")+ylab("Conspecific probability")+ggtitle(label="e. Surface Wind Speed"))+

(ggplot(na.omit(dfnn),aes(swv,Consp))+stat_smooth(method="glm", method.args = list(family = "binomial"))+
  theme_bw()+xlab("m/s")+ylab("Conspecific probability")+ggtitle(label="f. Sea water velocity"))



(ggplot(subset(dfnn,Species=="SBPR"|Species=="BPT"|Species=="ANPR"),aes(swv,Consp))+stat_smooth(method="glm", method.args = list(family = "binomial"))+
    theme_bw()+xlab("m/s")+ylab("Conspecific probability")+ggtitle(label="f. Sea water velocity"))+
  facet_wrap(Species~.)

(ggplot(subset(dfnn,Species=="BBA"|Species=="WA"|Species=="LMA"|Species=="SRA"|Species=="GHA"),aes(chl,Consp))+stat_smooth(method="glm", method.args = list(family = "binomial"))+
    theme_bw()+xlab("m/s")+ylab("Conspecific probability")+ggtitle(label="f. Sea water velocity"))+
  facet_wrap(Species~.)


library(lmerTest)

dfnn$zchl<-scale(dfnn$chl)
spf
dfN<-merge(dfnn,spf)

dfN5<-subset(dfN,totcount>5)

lmer1<-glmer(Consp~chl+(chl|Species),data=dfN5,family="binomial")
lmer01<-glm(Consp~chl,data=dfN5,family="binomial")
summary(lmer1)

anova(lmer1,lmer01)

lmer2<-glmer(Consp~zoo+(zoo|Species),data=dfN5,family="binomial")
lmer02<-glm(Consp~zoo,data=dfN5,family="binomial")
summary(lmer2)

anova(lmer2,lmer02)

lmer3<-glmer(Consp~sst+(sst|Species),data=dfN5,family="binomial")
lmer03<-glm(Consp~sst,data=dfN5,family="binomial")
summary(lmer3)

anova(lmer3,lmer03)


lmer4<-glmer(Consp~mld+(mld|Species),data=dfN5,family="binomial")
lmer04<-glm(Consp~mld,data=dfN5,family="binomial")
summary(lmer4)

anova(lmer4,lmer04)

lmer5<-glmer(Consp~ws+(ws|Species),data=dfN5,family="binomial")
lmer05<-glm(Consp~ws,data=dfN5,family="binomial")
summary(lmer5)

anova(lmer5,lmer05)


lmer6<-glmer(Consp~swv+(swv|Species),data=dfN5,family="binomial")
lmer06<-glm(Consp~swv,data=dfN5,family="binomial")
summary(lmer6)

anova(lmer6,lmer06)



plot_model(lmer1,type="emm",terms="chl[all]",pred.type="re")+theme_bw()+
  xlab("mg/m3")+ylab("Conspecific probability")+ggtitle(label="a. Chlrophyll-a concentration")+
#plot_model(lmer2,type="emm",terms="zoo[all]",pred.type="re")+theme_bw()+
#  xlab("g/m2")+ylab("Conspecific probability")+ggtitle(label="b. Mass content of zooplankton")+
plot_model(lmer3,type="emm",terms="sst[all]",pred.type="re")+theme_bw()+
  xlab("°C")+ylab("Conspecific probability")+ggtitle(label="b. Sea surface temperature")+
plot_model(lmer4,type="emm",terms="mld[all]",pred.type="re")+theme_bw()+
  xlab("meters")+ylab("Conspecific probability")+ggtitle(label="c. Thickness of the mixed layer")+
plot_model(lmer5,type="emm",terms="ws[all]",pred.type="re")+theme_bw()+
  xlab("m/s")+ylab("Conspecific probability")+ggtitle(label="d. Surface Wind Speed")
#plot_model(lmer6,type="emm",terms="swv[all]",pred.type="re")+theme_bw()+
 # xlab("m/s")+ylab("Conspecific probability")+ggtitle(label="f. Sea water velocity")
  
plot_model(lmer1,type="re",sort.est="sort.all",grid=F)[[1]]+
  ggtitle(label="a.Chlrophyll-a concentration")+theme_bw()+ylab("odds ratio")+xlab("taxa")+
#plot_model(lmer2,type="re",sort.est="sort.all",grid=F)[[2]]+
 # ggtitle(label="b. Mass content of zooplankton")+theme_bw()+ylab("odds ratio")+xlab("taxa")+
plot_model(lmer3,type="re",sort.est="sort.all",grid=F)[[2]]+
  ggtitle(label="b. Sea surface temperature")+theme_bw()+ylab("odds ratio")+xlab("taxa")+
plot_model(lmer4,type="re",sort.est="sort.all",grid=F)[[1]]+ylim(0.95,1.05)+
  ggtitle(label="c. Thickness of the mixed layer")+theme_bw()+ylab("odds ratio")+xlab("taxa")+
plot_model(lmer5,type="re",sort.est="sort.all",grid=F)[[2]]+
  ggtitle(label="d. Surface Wind Speed")+theme_bw()+ylab("odds ratio")+xlab("taxa")
#plot_model(lmer6,type="re",sort.est="sort.all",grid=F,ci.lvl=0.25)[[2]]+
 # ggtitle(label="f. Sea water velocity")+theme_bw()+ylab("odds ratio")+xlab("taxa")



(plot_model(lmer1,type="diag")[[1]]+
plot_model(lmer2,type="diag")[[1]]+
plot_model(lmer3,type="diag")[[1]]+

plot_model(lmer4,type="diag")[[1]]+
plot_model(lmer5,type="diag")[[1]]+
plot_model(lmer6,type="diag")[[1]])







