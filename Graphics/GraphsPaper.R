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
library(gridExtra)

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

#####  Four models of survival

p1 = datosTotal %>%
  group_by(IDX) %>% 
  mutate(n=n()) %>% filter(n==1600)%>% 
  ggplot(aes(x=Mes, y = Sobrevida, colour=TM)) +
  geom_smooth(se=FALSE) + theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="Complete",values=c("red", "orange", "blue", "lightblue"))+
  theme_bw()

p2 = datosTotal %>% 
  group_by(IDX) %>% 
  mutate(n=n()) %>% filter(n==800)%>% 
  ggplot(aes(x=Mes, y = Sobrevida, colour=TM)) +
  geom_smooth(se=FALSE) + theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="Imputed",
                     values=c("orange", "lightblue","blue","red"))+
  theme_bw()

grid.arrange(p1, p2)

p3 = datosTotal %>% 
  group_by(IDX) %>% 
  mutate(n=n(), Lost= ifelse(n()==800, "Missing data", "No missing data"))%>% 
  ggplot(aes(x=Mes, y = Sobrevida, colour=TM)) +
  facet_wrap(.~Lost)+
  geom_smooth(se=TRUE, level=0.99) + theme_bw()+
  xlab("Months") + ylab("Predicted survival")+
  scale_color_manual(name="Model",
                     values=c("orange", "lightblue","blue","red")) +
  theme_bw()                             
p3



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
