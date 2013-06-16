out<-read.table("miapacaS.txt")
matplot(out[,1],out[,2:3],type="l",xlab="time (hours)",ylab="count",main="cells in S-phase after 7 nM paclitaxel (0-24h)",log="y",lwd=2, xaxs="i")
legend(40,1e4, c("tumour","normal"), col = 1:2, lty = 1:2,bty="n")
