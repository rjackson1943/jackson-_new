# print headers for an output table
#
string1<-scan("headers.txt",what="",nmax=3,sep=",")
stringX<-c(string1[1],'     ',string1[2],'     ',string1[3])
write.table(stringX,"R_output.txt",eol="",quote=F,row.names=F,col.names=F)
