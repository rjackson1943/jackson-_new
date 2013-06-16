# dose-response curve plotting program
#
string1<-scan("headers.txt",what="",nmax=3,sep=",")
labelX<-string1[1]
labelY<-string1[2]
X<-matrix(scan("myelo2.dat",0),ncol=3, byrow=TRUE)

matplot(X[ , 1], X[ , 2:3], type = "l", xlab = "Seliciclib (micromolar)", ylab = "Percent Control",
		main = "Effect of Seliciclib on ROS and Circulating Neutrophils", lwd = 2)
 legend(60,100, c("ROS", "neutrophils"), col = 1:2, lty = 1:2, bty="n")
