out<-read.table("myelo.dat")
matplot(out[,1],out[,2:6],type="l",xlab="time (hours)",ylab="count",main="Simulation 3: Continuous exposure to 50 nM interleukin 8",log="y",lwd=2, xaxs="i")
legend(40,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")
