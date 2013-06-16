out<-read.table("table4.txt")
matplot(out[,1],out[,2:7],type="l",xlab="time (hours)",ylab="count",main="7 nM paclitaxel (0-24h)",log="y",lwd=2, xaxs="i")
legend(20,1e4, c("tumour","tumour M","tumour G1","normal","normal M","Normal G1" ), col = 1:7, lty = 1:7,bty="n")
 