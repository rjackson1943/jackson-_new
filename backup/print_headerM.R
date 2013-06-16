# print headers for a myelo output table
#
string1<-scan("headerM.txt",what="",nmax=7,sep="")
stringX<-c(string1[1],'     ',string1[2],'     ',string1[3],'     ',string1[4],'     ',string1[5],'     ',string1[6])
write.table(stringX,"R_output.txt",eol="",quote=F,row.names=F,col.names=F)
