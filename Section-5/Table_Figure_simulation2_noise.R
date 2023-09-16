#################
#### Table for result according to noise
### Figure main text 
## 2 tables:  MSE appendix and sensibility and specificity
################
library(latex2exp)
library(ggplot2)
library(ggpubr)
load("result_SIMU2new_noise1.5.Rdata")


MSE <- cbind(result$MSE.spls,result$MSE.Bss)
colnames(MSE) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")

sens <- cbind(result$sens.spls,result$sens.bss)
colnames(sens) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")

spe <- cbind(result$spe.spls,result$spe.bss)
colnames(spe) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")



MSEP <- matrix(c(MSE),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(MSE),each=100)
data.MSE <- data.frame(MSEP,model)
MSEPplot <- ggplot(data.MSE, aes(x=model, y=MSEP)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=1.5$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  

Sensitivity <- matrix(c(sens),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(sens),each=100)
data.sens <- data.frame(Sensitivity,model)
sens.plot <- ggplot(data.sens, aes(x=model, y=Sensitivity)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=1.5$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  

Specificity <- matrix(c(spe),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(spe),each=100)
data.spe <- data.frame(Sensitivity,model)
speplot <- ggplot(data.spe, aes(x=model, y=Specificity)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=1.5$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  

ggarrange(MSEPplot, sens.plot,speplot,  ncol=3,nrow=1,common.legend = TRUE)

ggsave("figsimu2newnoise1.5.pdf", width = 10, height = 5)

load("result_SIMU2new_noise3.Rdata")


MSE <- cbind(result$MSE.spls,result$MSE.Bss)
colnames(MSE) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")

sens <- cbind(result$sens.spls,result$sens.bss)
colnames(sens) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")

spe <- cbind(result$spe.spls,result$spe.bss)
colnames(spe) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")

MSEP <- matrix(c(MSE),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(MSE),each=100)
data.MSE <- data.frame(MSEP,model)
MSEPplot <- ggplot(data.MSE, aes(x=model, y=MSEP)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=3$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  


Sensitivity <- matrix(c(sens),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(sens),each=100)
data.sens <- data.frame(Sensitivity,model)
sens.plot <- ggplot(data.sens, aes(x=model, y=Sensitivity)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=3$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  

Specificity <- matrix(c(spe),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(spe),each=100)
data.spe <- data.frame(Sensitivity,model)
speplot <- ggplot(data.spe, aes(x=model, y=Specificity)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=3$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  

ggarrange(MSEPplot, sens.plot,speplot,  ncol=3,nrow=1,common.legend = TRUE)
ggsave("figsimu2newnoise3.pdf", width = 10, height = 5)


load("result_SIMU2new_noise6.Rdata")


MSE <- cbind(result$MSE.spls,result$MSE.Bss)
colnames(MSE) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")

sens <- cbind(result$sens.spls,result$sens.bss)
colnames(sens) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")

spe <- cbind(result$spe.spls,result$spe.bss)
colnames(spe) <- c("spls(H=1)","spls(H=2)","BSS(H=1)","BSS(H=2)")

MSEP <- matrix(c(MSE),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(MSE),each=100)
data.MSE <- data.frame(MSEP,model)
MSEPplot <- ggplot(data.MSE, aes(x=model, y=MSEP)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=6$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  


Sensitivity <- matrix(c(sens),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(sens),each=100)
data.sens <- data.frame(Sensitivity,model)
sens.plot <- ggplot(data.sens, aes(x=model, y=Sensitivity)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=6$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  

Specificity <- matrix(c(spe),ncol=1,nrow=400,byrow=F)
model <- rep(colnames(spe),each=100)
data.spe <- data.frame(Sensitivity,model)
speplot <- ggplot(data.spe, aes(x=model, y=Specificity)) + 
  geom_boxplot(fill="gray")+
  labs(title=(TeX("noise: $\\sigma=6$")),x="")+
  theme_classic()+ theme(text = element_text(size = 17),axis.text.x = element_text(angle = -90))  

ggarrange(MSEPplot, sens.plot,speplot,  ncol=3,nrow=1,common.legend = TRUE)
ggsave("figsimu2newnoise6.pdf", width = 10, height = 5)

result$SNR




