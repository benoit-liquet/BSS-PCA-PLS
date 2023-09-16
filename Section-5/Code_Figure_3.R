pdf("dimension-p-noise5.pdf",width=10,height=10)
load("result_SIMU1_high_noise550.Rdata")
result.int <- result
result.int$tabcrit.cpls <- - result.int$tabcrit.cpls
result.int$tabcrit.spls <- - result.int$tabcrit.spls
result.int$tabcrit.true <- - result.int$tabcrit.true
xmin <- min(c(result.int$tabcrit.true,result.int$tabcrit.cpls))
ymax <- max(c(result.int$tabcrit.true,result.int$tabcrit.cpls))


res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.50 <- cbind(res.bss,res.spls)

xlim <- c(xmin,ymax)
par(mfrow=c(2,2))
plot(result.int$tabcrit.true,result.int$tabcrit.cpls,xlab=expression(paste(delta[1]," for generated subset")),ylab=expression(paste(delta[1], " for BSS-PLS subset")),main=expression(p==50),xlim=xlim,ylim=xlim)
abline(0,1)
sum(result.int$tabcrit.cpls>result.int$tabcrit.spls)
sum(result.int$tabcrit.spls==result.int$tabcrit.cpls)
load("result_SIMU1_high_noise3100.Rdata")
load("result_SIMU1_high_noise5100.Rdata")
result.int <- result
result.int$tabcrit.cpls <- - result.int$tabcrit.cpls
result.int$tabcrit.spls <- - result.int$tabcrit.spls
result.int$tabcrit.true <- - result.int$tabcrit.true
xmin <- min(c(result.int$tabcrit.true,result.int$tabcrit.cpls))
ymax <- max(c(result.int$tabcrit.true,result.int$tabcrit.cpls))
xlim <- c(xmin,ymax)
plot(result.int$tabcrit.true,result.int$tabcrit.cpls,xlab=expression(paste(delta[1]," for generated subset")),ylab=expression(paste(delta[1], " for BSS-PLS subset")),main=expression(p==100),xlim=xlim,ylim=xlim)
abline(0,1)
#plot(result.int$tabcrit.spls,result.int$tabcrit.cpls,ylab="BSS-PLS",xlab="sparse PLS",main="p=100, model size=90")
#abline(0,1)
sum(result.int$tabcrit.cpls>result.int$tabcrit.spls)
sum(result.int$tabcrit.spls==result.int$tabcrit.cpls)

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.100 <- cbind(res.bss,res.spls)




load("result_SIMU1_high_noise5200.Rdata")
result.int <- result
result.int$tabcrit.cpls <- - result.int$tabcrit.cpls
result.int$tabcrit.spls <- - result.int$tabcrit.spls
result.int$tabcrit.true <- - result.int$tabcrit.true
xmin <- min(c(result.int$tabcrit.true,result.int$tabcrit.cpls))
ymax <- max(c(result.int$tabcrit.true,result.int$tabcrit.cpls))
xlim <- c(xmin,ymax)
plot(result.int$tabcrit.true,result.int$tabcrit.cpls,xlab=expression(paste(delta[1]," for generated subset")),ylab=expression(paste(delta[1], " for BSS-PLS subset")),main=expression(p==200),xlim=xlim,ylim=xlim)
abline(0,1)
#plot(result.int$tabcrit.cpls~result.int$tabcrit.spls,ylab="BSS-PLS",xlab="sparse PLS",main="p=200,model size=190")
#abline(0,1)
sum(result.int$tabcrit.cpls>result.int$tabcrit.spls)
sum(result.int$tabcrit.spls==result.int$tabcrit.cpls)

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.200 <- cbind(res.bss,res.spls)



load("result_SIMU1_high_noise5500.Rdata")
result.int <- result
result.int$tabcrit.cpls <- - result.int$tabcrit.cpls
result.int$tabcrit.spls <- - result.int$tabcrit.spls
result.int$tabcrit.true <- - result.int$tabcrit.true
xmin <- min(c(result.int$tabcrit.true,result.int$tabcrit.cpls))
ymax <- max(c(result.int$tabcrit.true,result.int$tabcrit.cpls))

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.500 <- cbind(res.bss,res.spls)



xlim <- c(xmin,ymax)
plot(result.int$tabcrit.true,result.int$tabcrit.cpls,xlab=expression(paste(delta[1]," for generated subset")),ylab=expression(paste(delta[1], " for BSS-PLS subset")),main=expression(p==500),xlim=xlim,ylim=xlim)
abline(0,1)
#plot(result.int$tabcrit.cpls~result.int$tabcrit.spls,ylab="BSS-PLS",xlab="sparse PLS",main="p=200,model size=190")
#abline(0,1)
sum(result.int$tabcrit.cpls>result.int$tabcrit.spls)
sum(result.int$tabcrit.spls==result.int$tabcrit.cpls)
dev.off()


res.tab <- rbind(res.50,res.100,res.200,res.500)
row.names(res.tab) <- c(50,100,200,500)

library(xtable)

tab <- xtable(res.tab,digits = c(0,3,3,3,3,3,3),caption = "Number of times sparse PLS and BSS-PLS retrieve the true best subset for different model size  over 100 runs according to different noisel level")
print(tab, caption.placement = "top",include.rownames = TRUE)



