rm(list=ls())
library(tidyverse)
library(readxl)
library(spBayesSurv)
library(survival)
library(coda)
library(car)
library(mnorm)
library(ldsr)
library(LaplacesDemon)
library(SurvMetrics)
library(latex2exp)
library(survminer)
library(gridExtra)
library(RColorBrewer)
library(wesanderson)

load("~/2024/Eliana/MyeloEli/Databases/DatosTotal.RData")
Arkansas <- read_excel("Databases/Arkansas.xlsx")


#### Graficas

tabA = summary(modelo)$table %>% data.frame()
rownames(tabA) = c("(Intercept)","PROTTT3","AGE","SEXmale","RACEwhite","ISOIgA",
                "ISOIgD","ISOIgG","ISONS","B2M","CRP","CREAT","LDH","ALB",
                "HGB","ASPC","BMPC","MRI","Cyto_Abn","Log(scale)")
tabA$Var = rownames(tabA)
tabA$Impor = -log10(tabA$p)
tabA$Model = "CoxComp"

tabB = summary(modeloCoxNA)$table %>% data.frame()
tabB$Var = rownames(tabA)
tabB$Impor = -log10(tabB$p)
tabB$Model = "CoxImp"


tabC = t(t(rowSums(modeloBayes$gamma)/2000)) %>% data.frame()
colnames(tabC) = "Impor"

tabC$Var =rownames(tabA)[2:19]
tabC$Model = "BayesComp"

tabD = t(t(rowSums(modeloBayesNAMediana$gamma)/2000)) %>% data.frame()
colnames(tabD) = "Impor"
tabD$Var =rownames(tabA)[2:19]
tabD$Model = "BayesImp"

Impor = rbind(tabA[2:19,] %>% dplyr::select(Impor, Var, Model),
              tabB[2:19,] %>% dplyr::select(Impor, Var, Model),
              tabC, tabD)

p1 = Impor %>% dplyr::filter(Model %in% c("BayesComp", "BayesImp")) %>%
  ggplot(aes(x=fct_reorder(Var, Impor), y=Impor, fill=Model))+
  geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  scale_fill_manual(name="Bayes Model",
                     labels=c("Complete", "Imputed"),
                     values=c("gray1", "gray50"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Variable") + ylab("Importance of Variable")



p2= Impor %>% dplyr::filter(!Model %in% c("BayesComp", "BayesImp")) %>%
  ggplot(aes(x=fct_reorder(Var, Impor), y=Impor, fill=Model))+
  geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  scale_fill_manual(name="Cox Model",
                    labels=c("Complete", "Imputed"),
                    values=c("gray1", "gray50"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Variable") + ylab(TeX(r"($-log_{10}(P-Value)$)"))

grid.arrange(p1, p2)                    

ggsave("Graphics/Importance.pdf")

#####  Four models of survival


p1 = datosTotal %>%
  group_by(IDX) %>% 
  mutate(n=n()) %>% filter(n==1600)%>% 
  ggplot(aes(x=Mes, y = Sobrevida, colour=TM)) +
  geom_smooth(se=FALSE) + theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="Complete",values=c("red", "orange", "blue", "lightblue"))+
  theme_bw() +
  coord_cartesian(xlim=c(0,150), ylim=c(0.3,1))

p2 = datosTotal %>% 
  group_by(IDX) %>% 
  mutate(n=n()) %>% filter(n==800)%>% 
  ggplot(aes(x=Mes, y = Sobrevida, colour=TM)) +
  geom_smooth( se=FALSE) + theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="Imputed",
                     values=c("orange", "lightblue","blue","red"))+
  theme_bw()+
  coord_cartesian(xlim=c(0,150),ylim=c(0.3,1))

grid.arrange(p1, p2)
p1




p3 = datosTotal %>% 
  group_by(IDX) %>% 
  mutate(n=n(), Lost= ifelse(n()==800, "Missing data", "No missing data"))%>% 
  filter(n==800| (Lost=="No missing data" & TM %in% c("B-Comp", "C-Comp"))) %>% 
  group_by(Mes, TM, Lost) %>% 
  summarise(SD = mean(Sobrevida),
         SDS = sd(Sobrevida),
         n=n()) %>% 
  mutate(linf =SD-1.96*SDS/(sqrt(n)),
         lsup = SD+1.96*SDS/(sqrt(n)), Sobrevida=SD) %>% 
  ggplot(aes(x=Mes, y = Sobrevida, group = TM, fill=TM)) +
  facet_grid(Lost~.)+
  geom_ribbon(aes(ymin=linf, ymax=lsup), alpha=0.2) +
  geom_line(aes(group = TM, color=TM), lwd=1.2)+
  xlab("Months") + ylab("Predicted survival")+
  theme_bw()                             
p3

ggsave("Graphics/Curvas_4_Models.pdf")

#### 
datosTotal = datosTotal %>%  left_join(Arkansas %>% dplyr::select(ISSR2, PATID), by=c("IDX"="PATID"))

p4 = datosTotal %>% 
  filter(TM=="B-Comp") %>% 
  group_by(IDX) %>% 
  mutate(n=n(), Lost= ifelse(n()==800, "Missing data", "No missing data"))%>% 
  ggplot(aes(x=Mes, y = Sobrevida, color=ISSR2)) +
  facet_wrap(.~Lost)+
  geom_smooth(se=TRUE, level=0.99) + theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="Model",
                     values=c("orange", "lightblue","blue","red")) +
  theme_bw()                             
p4



p5 = datosTotal %>% 
  filter(TM=="B-Imp") %>% 
  group_by(IDX) %>% 
  mutate(n=n(), Lost= ifelse(n()==800, "Missing data", "No missing data"))%>% 
  ggplot(aes(x=Mes, y = Sobrevida, color=ISSR2)) +
  facet_wrap(.~Lost)+
  geom_smooth(se=TRUE, level=0.99) + theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="Model",
                     values=c("orange", "lightblue","blue","red")) +
  theme_bw()                             
p5


##### Graphs of Kaplan Meier

fit = survfit(Surv(OS_Time, OS_Censor)~ISSR2, data=Arkansas)
pKM = ggsurvplot(fit)

pKM$plot  =pKM$plot+theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="ISSR2",
                     values=c("gray10", "gray40","gray60"), guide="none")+
  theme_bw()
pKM = pKM$plot

pKM + geom_smooth(data=datosTotal %>% 
                    group_by(IDX) %>% 
                    mutate(n=n(), Lost= ifelse(n()==800, "Missing data", "No missing data"),
                           ISSR2 = paste("ISSR2=",ISSR2, sep="" )),
                  aes(x=Mes, y = Sobrevida, fill=TM, colour = ISSR2),
                  se=TRUE, level=0.99) + theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="ISSR2",
                     values=c("gray10", "gray40","gray60","gray10", "gray40","gray60"))+
  theme_bw()+
  coord_cartesian(ylim=c(0.35,1), xlim=c(0,100), expand = FALSE)+
  scale_fill_manual(name="Model",values = c("skyblue", "skyblue3", "chocolate", "orange"))

