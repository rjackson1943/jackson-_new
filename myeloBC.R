# myeloBC.R  8 May 13
# modelling progression of CML
#
PS<-c(2e7,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
P<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
N<-c(4.55e6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DT<-c(1.1,0.9,2.7,3.1,3.0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
k<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
frac<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
m<-c(1e-5, 1e-6, 1e-7, 3e-7, 1e-9)           # mutation rates
apo=-.126                                    # rate of neutrophil apoptosis
df=.0419                                     # differentiation rate
Kn1<-5e6                                     # Ki for neutrophil feedback of proliferation
nc<-7                                        # number of cell compartments
volP<-2900                                   # plasma volume
volM<-1400                                   # marrow volume
time<-0
Ccapacity<-2e7                               # carrying capacity
for(i in seq(1,33,1))P[i]<-PS[i]
#
for(i in seq(1,nc,1))k[i]<-log(2)/DT[i]
k[1]<-log(2)/DT[1]/(1+(N[1]+N[2])/Kn1)              # feedback
k[2]<-log(2)/DT[2]/(1+(N[1]+N[2])/Kn1)
# P[n] are progenitor cells
# PS[n] are starting values of P
# N[n] are circulating myeloid cells
dx<-1.0
tlimit<-365
#
# proliferation
#
for(i in seq(0,tlimit-1,1)) {
time=time+dx
for(i2 in seq(1,nc,1))P[i2]<-P[i2]*exp(k[i2]*dx)
#
# compartment transitions
for(i2 in seq(3,nc,1))P[i2]<-P[i2]+m[i2-2]*P[2]*dx
P[2]=P[2]-(m[1]+m[2]+m[3]+m[4]+m[5])*P[2]*dx
#
# apoptosis
N[1]<-N[1]*exp(apo*dx)
#
# differentiation
P[1]<-P[1]*(1-df)
P[2]<-P[2]*(1-df)
N[1]<-N[1]+P[1]*df*volM/volP
N[2]<-N[2]+P[2]*df*volM/volP
#
# objective function
for(j in seq(1,nc,1))frac[j]=P[j]/sum(P)
if(sum(P)>Ccapacity){for(j1 in seq(1,nc,1))P[j1]<-frac[j1]*Ccapacity}
#
# print(P[1:nc])
# print(N[1:2])
}
time
print(P[1:nc])
print(N[1:2])
#Lst<-list(days=time, wbc=P[1])
#print(Lst)
#sink("C:/R/myeloBC.dat")
#print(P[1:nc])
#sink()
#sum(P)



