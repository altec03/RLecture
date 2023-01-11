:/Users/cmc/Desktop/R_lec/")
getwd() # error. Use Project
install.packages(readxl)
library(readxl)

#1.?Ñº??kages("ggplot2")
library(ggplot2)
#ggplot(dataa)
b=na.omit(a)

Stage=factor(b$st # errora
a=read_excel(here::here("data/scatter_data.xlsx")) # error free

ge,levels=c("1","2"),labels=c("Stage I/II","Stage III/IV"))
c=ggplot(b,aes(x=year_s,y=percent_s,color=Stage,shape=Stage))+geom_point(na.rm=T,size=3)
c
d=c+scale_shape_manual(values=c(1,4))
d
##?ß¼??? ?×¸??? 
e=d+stat_smooth(method=x_continuous(breaks=seq(2008,2017,3))+
  xlab("Year")+ylab("Percent (%)")+
  annotate("text",label="R^2== 0.02",parse=T,x=2013,y=75,color="Dodger Blue")+
  annotate("text",label="R^2== 0.02",parse=T,x=2013,y=27,color="tomato")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

e

##??Á¦ 2) three groups
a=read_excel("scaor(a$risk, levels=c("1","2","3"),l # error 
aa=read_excel(here::here("data/scatter_data.xlsx")) # error free

bels=c("Low risk","Intermediate risk","High risk"))
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

##??Á¦
a=read_excel("roc_data.xlsx")
head(a)esult_group1~("st_group1,data= # erroraa=read_excel(here::here("data/roc_data.xlsx")) # error free
)
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

##??Á¦
a=read_excel("survival_data.xlsx")
fit=su(f_u_months,survival==1)~expression # error
a=read_excel(here::here("data/survival_data.xlsx")),
data=a)
ggsurvplot(fit,data=a)
ggsurvplot(fit,risk.table=T,fun='pct',linetype=c(2,1),palette=c("black","red"),
           pval="p = 0.132",pval.coord=c(1,10),
           xlab="Months",
           legend.title="Expression",legend.labs=c("Low","High"),
           legend=c(0.9,0.9),font.legend="bold")

package_version(R.version) 
