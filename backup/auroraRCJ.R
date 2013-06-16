#(drug=exp(seq(-2.5,1,.1)))
(drug=10^seq(-2.5,1,.1))
ATP=1000
KaA=22
KaB=94
KiA=.008
KiB=.009
A=(1+ATP/KaA)/(1+drug/KiA+ATP/KaA)
B=(1+ATP/KaB)/(1+drug/KiB+ATP/KaB)
logD=log10(drug)
plot(logD,A,type="l",col="red")
lines(logD,B,col="blue")

## same lines
#IC50A=KiA*(1+ATP/KaA)
#IC50B=KiB*(1+ATP/KaB)
#A=1-drug/(drug+IC50A)
#B=1-drug/(drug+IC50B)
#plot(drug,A,log="x",type="l",col="red")
#lines(drug,B,col="blue")


KiA=.008
KiA=.00043  # MLN8237
(IC50A=KiA*(1+ATP/KaA)) # 6.7 nM observed, 20 nM expected here

#IC50B=KiB*(1+ATP/KaB)
IC50B=1.53
(KiB=IC50B/(1+ATP/KaB)) # 131 nM; 131/.43 = 300 fold tighter for A

A=1-drug/(drug+IC50A)
B=1-drug/(drug+IC50B)

plot(drug,A,log="x",type="l",col="red")
lines(drug,B,col="blue")

library(deSolve)

F=0.343 #bioavailabily
V=3100  #Volume of distribution in mL/Kg
ka= .1517 # Rate constant for oral uptake in 1/min
ke=0.0339 #Rate constant for elimination
t=10^seq(0,2.5,.1)
D=10e6 # ng/kg
y=F*D*ka*(exp(-ke*t)-exp(-ka*t))/(V*(ka-ke)) # ng/mL as in paper
plot(log10(t),y,type="l",ylab="ng/ml")
# redo using concentrations in uM
y=y/368.5 # ng/nmoles => nmoles/ml = uM
plot(log10(t),y,type="l",ylab="uM")

library(odesolve)
# inverse minutes
k01=0.04 #U to ME state transition rate [21]
k02=0.032 #U to S state transition rate Set here
k03=0.178 # U to C state transition rate [21]
k04=0.13 #S to C state transition rate Set here
k05=0.3   #ME to C state transition rate Set here

# in cubic um per min
k0k=1200 #SAC protein-kinetochore interaction rate Set here
k=0.2 #SAC protein-protein interaction rate Set here
alph=0.2 #Degradn. rate of inhibited SAC complexes
vc=6000 # fL= um^3 = Cell volume, i.e. 6 pL
N=800000 # CDC20 proteins per cell



# define the ODE right hand side
faurora <- function(t, X, p)
{ # first define the fluxes
  U=X[1]; S=X[2];  ME=X[3]
  vU2ME =k01*A * U
  vU2S =k02*A * U
  vU2C =k03*A * U
  vS2C =k04*A*B * S
  vME2C =k05*A*B * ME
  Ne=X[4];  Nes=X[5];  Nec=X[6];  Neg=X[7]
#  vE2es=(k0k*B/vc)*(round(X[1])+round(X[2]))*Ne
#  vE2ec=(k*A*B/vc)*(round(X[1])+round(X[2]))*Ne*Nes
  vE2es=(k0k*B/vc)*(X[1]+X[2])*Ne
  vE2ec=(k*A*B/vc)*(X[1]+X[2])*Ne*Nes
  vDes=alph*Nes
  vDec=alph*Nec
  vDeg=alph*Neg
  # now define dX/dt (Xi prime) as fluxes into a node minus fluxes out
  X1p = -vU2ME-vU2S-vU2C
  X2p = vU2S-vS2C
  X3p = vU2ME-vME2C
  X4p = -vE2es-vE2ec+vDes+vDec+vDeg
  X5p = vE2es-vDes
  X6p = vE2ec-vDec
  X7p = -vDeg
  XP = c(X1p, X2p, X3p, X4p, X5p, X6p, X7p);
  V=c(C=46-sum(X[1:3]),vU2ME,vU2S,vU2C,vS2C,vME2C,vE2es,vE2ec,vDes,vDec,vDeg)
  names(V)<-c("C","vU2ME","vU2S","vU2C","vS2C","vME2C","vE2es","vE2ec","vDes","vDec","vDeg")
  list(XP,V)}
y0=c(U=46,S=0,ME=0,Ne=0,Nes=0,Nec=0,Neg=800000);
A=0.5; B=0.5
A=0.75; B=0.75
A=1; B=1

out1=lsoda(y=y0,times=0:80,faurora, parms=c(test1=1), rtol=1e-4, atol= rep(1e-4,7))
#ny0=out1[nrow(out1),2:8]
(outs=data.frame(out1))
par(mfrow=c(2,1))
with(outs[1:20,],{plot(time,U,type="l",xlab="minutes",ylab="# of Kinetochores");lines(time,S);lines(time,ME);lines(time,C)})
attach(outs)
plot(time,Neg,type="l",xlab="minutes",ylab="# of CDC20s")
lines(time,Nes)
lines(time,Nec)
lines(time,Ne)
#par(mfrow=c(1,1))
detach(outs)
