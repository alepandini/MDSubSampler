traj5_full <- read.csv("full.csv", header=FALSE)
traj5_0.05 <- read.csv("0.05.csv", header=FALSE)
traj5_0.1 <- read.csv("0.1.csv", header=FALSE)
traj5_0.5 <- read.csv("0.5.csv", header=FALSE)
traj5_1 <- read.csv("1.csv", header=FALSE)
traj5_5 <- read.csv("5.csv", header=FALSE)
traj5_10 <- read.csv("10.csv", header=FALSE)


plot(density(traj5_full$V2), xlab="Distance between G55 and P127", main="", ylim = c(0,1))
lines(density(traj5_full$V2), type = 'h', col = rgb(0.4,0.4,0.4,0.3), lty=7)
lines(density(traj5_0.05$V2), col = 'black', lty=2)
lines(density(traj5_0.05$V2), type = 'h', col = rgb(1,0,0,0.2))
lines(density(traj5_0.1$V2), col = 'black', lty=3)
lines(density(traj5_0.1$V2), type = 'h', col = rgb(0,1,0,0.2))
lines(density(traj5_0.5$V2), col = 'black', lty=6)
lines(density(traj5_0.5$V2), type = 'h', col = rgb(0,0,1,0.2))

par(mfrow=c(1,2))
par(las = 2) 
par(mar = c(4,4,3,1))
plot(density(traj5_full$V2), xaxt = 'n', xlab="Distance between G55 and P127", main="", ylim = c(0,1))
axis(1, seq(0,5,0.25))
lines(density(traj5_full$V2), type = 'h', col = rgb(0.4,0.4,0.4,0.2), lty=7)
lines(density(traj5_0.05$V2), col = 'blue', lty=2)
lines(density(traj5_0.1$V2), col = 'red', lty=3)
lines(density(traj5_0.5$V2), col = 'green', lty=4)
legend(x = "topright",        
       legend = c("Full trajectory", "0.05% sample", "0.1% sample", "0.5% sample"), 
       lty = c(7,2,3,4),          
       col = c("black","blue","red","green"),          
       lwd = 2,
       bty = 'n')                

plot(density(traj5_full$V2),xaxt = 'n', xlab="Distance between G55 and P127", main="", ylim = c(0,1))
axis(1, seq(0,5,0.25))
lines(density(traj5_full$V2), type = 'h', col = rgb(0.4,0.4,0.4,0.2), lty=7)
lines(density(traj5_1$V2), col = 'blue', lty=2)
lines(density(traj5_5$V2), col = 'red', lty=3)
lines(density(traj5_10$V2), col = 'green', lty=4)
legend(x = "topright",        
       legend = c("Full trajectory", "1% sample", "5% sample", "10% sample"), 
       lty = c(7,2,3,4),          
       col = c("black","blue","red","green"),          
       lwd = 2,
       bty = 'n')                                


