# myeloL2.R Model of CML 23 Nov 2012
#
library(deSolve)
# state variables are capitalized
# P   = progenitors (old CFU-GM)
# Q   = quiescent progenitors
# N   = Neutrophils
#TN   = Tissue Neutrophils
#AN   = Activated Neutrophils
# D   = dead cell
#GM   = GM-CSF   
#NAP2 = neutrophil activating peptide-2 
#RAC2 = G protein associated with increased ROS production
#Mcl1 = mantle cell lymphoma 1 (antiapoptosis factor)
#ROS  = reactive oxygen species

# Drugs used to manipulate the system
#bu = Busulfan = DNA crosslinker  (1)  drop state variable P by 10^(2*bu/20)
#mo = MOR103 = Ab blocker of GM   (2)
#das= dasatinib                   (3)
#sb = SB272844 = IL8 inhibitor    (4)
#cn = CNDAC = sapacitabine becomes this = SSB producer = DSBs in S via HR  (5)
#seli = seliciclib = cdk9 inhibitor, blocks transcription of Mcl-1 (6)
#imab = imatinib = anti-TNF monoclonal antibody (7)
#ifna = interferon alpha (8)
#antiox = ascorbic acid, antoxidant (9)


fmyelo <- function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				vX.P =Vx.p*(1+il3/Kil3)
				vP.P =P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*(1-Qf)
				vP.N =P*Vp.n*g/(g+Kg)
				vN.TN =N*Vn.tn*il8/(il8+Kil8)/(1+sb/Ksb)
				vTN.X =TN*Vtn.x/(1+Mcl1/Kmcl1)
				vTN.AN =TN*Vtn.an*NAP2/(NAP2+Knap2)
				vAN.TN =AN*Van.tn
				vP.Q =P*Vp.p*GM/(GM+Kgm)/(1+cn/Kcn)*Qf
				vQ.P =Q*Vq.p/(1+P/Kp)
				vX.GM =Vx.gm/(1+N/Kn)/(1+mo/Kmo)
				vGM.X=GM*Vgm.x
				vX.NAP2=Vx.nap2*(1+BcrAbl/(Kabl*(1+das/Kdas)))
				vNAP2.X=NAP2*Vnap2.x
				vX.RAC2=Vx.rac2*(1+AN/Kan)
				vRAC2.X=RAC2*Vrac2.x
				vN.X=N*Vn.x/(1+Mcl1/Kmcl1)
				vX.MCL1=Vx.mcl1*(1+BcrAbl/Kabl2)/(1+seli/Kseli)/(1+ifna/Kifna)
				vMCL1.X=Mcl1*Vmcl1.x
				vX.ROS=RAC2*Vx.ros/(RAC2+Krac2)
				vGSH.GSSG=ROS*Vgsh.gssg/(ROS+Kros)*(gsh/Kgsh/(1+gsh/Kgsh))
				vN.M=P*Vn.m
				vGSSG.GSH=GSSG*Vgssg.gsh/(GSSG+Kgssg)
				vROS.X=ROS*antiox*Kros2
				vBACT.X=bact*Vbact.x*ROS^nHr/(ROS^nHr+kkmax)
				vX.ABL=0
				vABL.X=0
				
				dP =vX.P+vP.P-vP.N-vP.Q+vQ.P-vN.M;
				dQ =vP.Q-vQ.P;
				dN =vP.N*volM/volP-vN.TN-vN.X
				dTN=vN.TN*volP/volT-vTN.X-vTN.AN+vAN.TN
				dAN=vTN.AN-vAN.TN
				dGM =vX.GM-vGM.X
				dNAP2   =vX.NAP2-vNAP2.X
				dRAC2   =vX.RAC2-vRAC2.X
				dMcl1   =vX.MCL1-vMCL1.X
				dROS    =vX.ROS-vGSH.GSSG-vROS.X
				dBACT   =0
				dabl    =vX.ABL-vABL.X
				dgsh	=vGSSG.GSH-vGSH.GSSG
				dGSSG   =vGSH.GSSG-vGSSG.GSH
				return(list(c(dP,dQ,dN,dTN,dAN,dGM,dNAP2,dRAC2,dMcl1,dROS,dBACT,dabl,dgsh,dGSSG),
   						    c(v0=vX.P,vP.P=vP.P,vP.N=vP.N,vN.TN=vN.TN,vTN.X=vTN.X,vTN.AN=vTN.AN,vAN.TN=vAN.TN,vP.Q=vP.Q,vQ.P=vQ.P,
vX.GM=vX.GM,vGM.X=vGM.X,vX.NAP2=vX.NAP2,vNAP2.X=vNAP2.X,vX.RAC2=vX.RAC2,vN.X=vN.X,
vX.MCL1=vX.MCL1,vMCL1.X=vMCL1.X,vN.M=vN.M,vBACT.X=vBACT.X,v14=vRAC2.X,v18=vX.ROS,
v19=vGSH.GSSG,v21=vGSSG.GSH,v22=vROS.X,v23=vX.ABL,v24 =vABL.X))
					   )
			})
}


pars=c(Vx.p=4.25e3, Kil3=1, Vp.p=.0203, Kgm=2, Kcn=100, Vp.n=.026, Kg=1, Kmo=1,
	  Vn.tn=.781, Kil8=396, Ksb=250, Vtn.x=.038, Kmcl1=2.438, Vtn.an=.001, Knap2=.1, 
	  Van.tn=.09, Vq.p=.53778, Kp=7e4, Vx.gm=2.007, Kn=5e5, Vgm.x=.05,
      Vx.nap2=.002, Kabl=1, Kdas=5, Vnap2.x=.2, Vx.rac2=1.5, Kan=100, Vrac2.x=.1, 
	  Vn.x=2e-3, Vx.mcl1=1, Kabl2=1, Vmcl1.x=.5, Vx.ros=980, Krac2=250,
      Vgsh.gssg=382.8, Vgssg.gsh=155.8,Kgssg=30,Kros=30, Kgsh=9.5e5, Vn.m=1.309e-3,Kros2=.104,
          Kseli=.3, Kifna=1, Vbact.x=5.756e-2, nHr=2.5,kkmax=320,Vx.abl=.22,Ktnfa=1,Kimab=1, Vabl.x=.0145,
	  g=11.0, il3=51.0, il8=4,                             
	  Qf=.0373, volP=2900, volM=1400, volT=65700,        # these are boundary  conditions
	  bu=0,mo=0,das=0,sb=0,cn=0,seli=0,imab=0,ifna=0,antiox=7)      # Drug concentrations

X0=c(P=1.222e7, Q=1.094e6, N=1.727e7, TN=9.296e5, AN=3875,
	 GM=1.129,NAP2=.01004, RAC2=26.59, Mcl1=2.021, ROS=1.967, bact=0, BcrAbl=5, gsh=1e6, GSSG=0)                                                                                                  # Starting values for variables
                                                                           # When bact>0 initialise il8 to bact*5e-4 +4
									   # When bact>0 initialise tnfa to bact*5e-6
									   # When bact>0 initialise abl to tnfa
times <- seq(1,80, by = 1)
out   <- ode(X0, times, fmyelo, pars)
head(out)
tail(out)  
X0=out[80,2:19]

# DNA cross-linking
#  bu=20 # set to MTD
#  (eventdat <- data.frame(var = "P", time = 0, value = 1/10^(2*bu/20), method = "mult"))
# Bob, this is a one time 100 fold drop in progenitor numbers due to busulfan
#  out   <- ode(X0, -20:600, fmyelo, pars,events = list(data = eventdat))
# tail(out) 

 matplot(out[ , 1], out[ , 2:6], type = "l", xlab = "time", ylab = "count",
		main = "Simulation 1: Uninhibited control", log="y", lwd = 2)
 legend(20,1e4, c("Progenitor Cells", "G0 Progenitor Cells","Blood Neutrophils","Tissue Neutrophils",
				"activated Neutrophils"), col = 1:5, lty = 1:5,bty="n")


