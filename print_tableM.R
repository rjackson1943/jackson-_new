# print headers for a myelo output table
#
str1<-scan("headerM.txt",what="",nmax=7,sep="")
stringX<-c(str1[1],' ',str1[2],' ',str1[3],' ',str1[4],' ',str1[5],' ', str1[6], ' ', str1[7])
write.table(stringX,"R_output.txt",eol="",quote=F,row.names=F,col.names=F)
write.table(" ","R_output.txt",append=T,quote=F,row.names=F,col.names=F)
Mdata<-matrix(scan("Mtable1.txt",0),ncol=7,byrow=T)
write.table(Mdata,"R_output.txt",append=T,sep=" ",eol="\n",row.names=F,col.names=F)
