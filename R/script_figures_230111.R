setwd("C:/Users/cmc/Desktop/R_lec/") # error. Use Project
getwd()
install.packages(readxl)
library(readxl)
library(ggplot2)

a=read_excel("scatter_data.xlsx") # error

a=read_excel(here::here("data/scatter_data.xlsx")) # error free

head(a)
b=na.omit(a)

Stage=factor(b$stage,levels=c("1","2"),labels=c("Stage I/II","Stage III/IV"))
c=ggplot(b,aes(x=year_s,y=percent_s,color=Stage,shape=Stage))+geom_point(na.rm=T,size=3)
c
d=c+scale_shape_manual(values=c(1,4))
d

e=d+stat_smooth(method=lm, se=F)+scale_x_continuous(breaks=seq(2008,2017,3))+
  xlab("Year")+ylab("Percent (%)")+
  annotate("text",label="R^2== 0.02",parse=T,x=2013,y=75,color="Dodger Blue")+
  annotate("text",label="R^2== 0.02",parse=T,x=2013,y=27,color="tomato")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

e


a=read_excel("scatter_data.xlsx") # error 

a=read_excel(here::here("data/scatter_data.xlsx")) # error free

RISK=factor(a$risk, levels=c("1","2","3"),labels=c("Low risk","Intermediate risk","High risk"))
b=ggplot(a,aes(x=year_r,y=percent_r,color=RISK,shape=RISK))+geom_point(size=3)
c=b+scale_shape_manual(values=c(1,2,4))
d=c+stat_smooth(method=lm, se=F)
d
d+scale_x_continuous(breaks=seq(2008,2017,3))
e=d+
  ylab("Percent (%)")+xlab("Year")+
  annotate("text",label="R^2== - 0.04",parse=T,x=2013,y=57,color="tomato")+
  annotate("text",label="R^2== - 0.03",parse=T,x=2013,y=8,color="Medium Sea Green")+
  annotate("text",label="R^2== - 0.02",parse=T,x=2013,y=20,color="Dodger Blue")

e+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                   
                  
 
#2.Multiple ROC curves
install.packages("pROC")
library(pROC)


a=read_excel("roc_data.xlsx") # error
a=read_excel(here::here("data/roc_data.xlsx")) # error free
head(a)

r1=roc(result_group1~test_group1,data=a)
r2=roc(result_group2~test_group2,data=a)
r3=roc(result_group3~test_group3,data=a)

plot(r1,legacy.axes=T)
plot(r2,add=T,col="red")
plot(r3,add=T,col="blue")

##color and legend
plot(r1,col="gray 40",legacy.axes=T)
plot(r2,add=T,col="Tomato")
plot(r3,add=T,col="Dodger Blue")

legend("bottomright",legend=c("Group 1: AUC=0.908","Group 2: AUC=0.775",'Group 3: AUC=0.605')
       ,fill=c("Tomato","gray 40","Dodger Blue"),border="white",box.lty=0,cex=1.5)



#3. Survival graph_number at risk
install.packages("survminer")
library(survival)
library(survminer)
library(ggplot2)


a=read_excel("survival_data.xlsx") # error
a=read_excel(here::here("data/survival_data.xlsx"))

fit=survfit(Surv(f_u_months,survival==1)~expression,data=a)
ggsurvplot(fit,data=a)
ggsurvplot(fit,risk.table=T,fun='pct',linetype=c(2,1),palette=c("black","red"),
           pval="p = 0.132",pval.coord=c(1,10),
           xlab="Months",
           legend.title="Expression",legend.labs=c("Low","High"),
           legend=c(0.9,0.9),font.legend="bold")

# package_version(R.version) 
