#################
#### Table for result according to noise
## 3 tables: first one frequency in main text, then MSE appendix and sensibility and specificity
################

load("result_SIMU1_NOISE1.5.Rdata")

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.nois1.5 <- cbind(res.bss,res.spls)
result$SNR

load("result_SIMU1_NOISE3.Rdata")

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.nois3 <- cbind(res.bss,res.spls)

result$SNR
load("result_SIMU1_NOISE6.Rdata")

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.nois6 <- cbind(res.bss,res.spls)
res.nois6
result$SNR
load("result_SIMU1_NOISE8.Rdata")

result$SNR
res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.nois8 <- cbind(res.bss,res.spls)



tab.Freq <- cbind(res.nois1.5[,c(1)],res.nois3[,c(1)],res.nois6[,c(1)],res.nois8[,c(1)],
      res.nois1.5[,c(5)],res.nois3[,c(5)],res.nois6[,c(5)],res.nois8[,c(5)])

tab.MSE <- cbind(res.nois1.5[,c(2)],res.nois3[,c(2)],res.nois6[,c(2)],res.nois8[,c(2)],
                  res.nois1.5[,c(6)],res.nois3[,c(6)],res.nois6[,c(6)],res.nois8[,c(6)])


tab.spe <- cbind(res.nois1.5[,c(3,4)],res.nois3[,c(3,4)],res.nois6[,c(3,4)],res.nois8[,c(3,4)],
                 res.nois1.5[,c(7,8)],res.nois3[,c(7,8)],res.nois6[,c(7,8)],res.nois8[,c(7,8)])


library(xtable)

tab.MSE <- xtable(tab.spe,digits = 2,caption = "Number of times sparse PLS and BSS-PLS retrieve the true best subset for different model size  over 100 runs according to different noisel level")
print(tab.MSE, caption.placement = "top")
