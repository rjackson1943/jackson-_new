# isobologram plotting program
model<-scan("grecfit.ini",what="",nmax=1)
string1<-scan(model,what="",nmax=3,sep=",")
drug1<-string1[1]
drug2<-string1[2]
isobol<-read.table("isobol.dat", header=TRUE)
plot(isobol$Drug_A, isobol$combination, xlab=drug1, ylab=drug2, type="n")
abline(1,-1, col=8)
lines(isobol$Drug_A, isobol$Upper_95, col=2)
lines(isobol$Drug_A, isobol$Lower_95, col=2)
lines(isobol$Drug_A, isobol$combination)
text(.8,1,string1[3])

