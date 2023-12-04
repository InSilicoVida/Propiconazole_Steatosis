########################################################
#Install all these packages
#############################################################
library(deSolve)
library(FME)
library(coda)
library(rootSolve)
# Propiconazole PBPK model  file for rats 
#Date 1 Dec 2023
#IF the dose is in ug/Kg BW/day, then concentration will be ug/L. If dose is in mg/Kg/day, concentration will be in mg/L and maount in mg.
#######################################################################


States = unlist(c(data.frame(
  Agut= 0,		 #amount of Propiconazole gut
  Aliver= 0,		 #amount of Propiconazole liver
  Abrain= 0,        #amount of Propiconazole in brain
  Akidney= 0,      #amount of Propiconazole in kidney
  Afat= 0,		 #amount ofPropiconazole fat
  Agonads= 0,      #amount of Propiconazole gonads
  Arestbody= 0,    #amount of Propiconazole rest of the body
  Aplasma= 0,		 #amount of Propiconazole blood plasma
  Aurine= 0
)))

States
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Physiological parameters 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#subject weight

BW = 0.25     #kg, you need to vary it based on experimental data.
MW=249.67     #Variable parameters, for propiconazole

#constant Fraction of blood flows to organs (blood flow rate)
QCC = 18.7                    #Total Cardiac blood output (L/h/kg)  
HCT = 0.45  				         #hematocrit percentage
FQliver = 0.174 		         #Fraction cardiac output going to liver 
FQbrain = 0.02  			       #Fraction cardiac output going to lung
FQkidney = 0.141  		       #Fraction cardiac output going to kidney                                                       
FQfat =  0.07  		         #Fraction cardiac output going to fat  
FQgonads = 0.0005            #fraction cardiac output going to gonads

#constant organ volume as a fraction of total body weight
Fliver = 0.034 		           #Fraction liver volume
Fbrain    = 0.006             #fraction lung volume  		         		                    
Fkidney = 0.0073 			       #Fraction kidney volume  
Fgonads = 0.0063             #fractional volume of gonads
Ffat = 0.07                 #fractional volume of fat  
Fplasma = 0.074  	         #fractional volume of plasma

#change the paramerers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# partition coefficient parameter For propiconazole
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# @Melina, you need to change these parameters.
k_liver_plasma = 1.78		 	   #Liver/blood partition coefficient  
k_brain_plasma = 1          #brain/blood partition coefficient
k_kidney_plasma = 2.63		 #kidney/blood partition coefficient
k_gonads_plasma = 1.7		 #gonads/blood partition coefficient  
k_fat_plasma =  2.78	     #Fat/blood partition coefficient  
k_restbody_plasma = 5.57     #Rest of the body/blood partition coefficient	

# absorption and elimination parameters need to scale per body weight   
kgut = 0.6     
fu = 1
Cl=1
frac = 1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#initialize parameters 
############################################################################
QCblood = QCC*((BW^0.74)) 		 			                    #Initial cardiac output for blood L/h  
QCplasma =QCblood *(1-HCT)   					              		#Adjust initial cardiac output for plasma flow  
Qliver= FQliver* QCplasma   	 				              		#Plasma flow to liver 
Qbrain= FQbrain*QCplasma 									  		#Plasma flow to brain
Qkidney= FQkidney*QCplasma 								  		#Plasma flow to kidney
Qgonads = FQgonads*QCplasma         				          		#plasma flow to gonads  
Qfat = FQfat*QCplasma                     				      		#plasma flow to fat  
Qrestbody = QCplasma  -(Qliver + Qbrain + Qkidney +  Qgonads + Qfat )  #plasma flow to rest of the body   (brain compartment removed)
vliver = Fliver * BW  									 	  		#Liver Volume  
vbrain = Fbrain*BW  											  	#volume of brain
vkidney = Fkidney*BW												#volume of kidney
vgonads = Fgonads * BW 	 							 	  		#gonads Volume 
vfat = Ffat*BW                   							  		#Volume of  fat    
vplasma = Fplasma*BW       							      		#volume of plasma  
vrestbody = 0.84*BW - vliver -vbrain- vkidney- vgonads- vfat - vplasma   #


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#    Dosing Parameters  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Differenzierbare Versionen einer Rechtecksfunktion - Feeding function

anaus <- function(t,t0,t1)
{
  y <- (tanh(100*(t-t0)) - tanh(100*(t-t1)))/2
  return(y)}

#Oral Dosing
multiple <- FALSE                                           # TRUE for multiple dosing, FALSE for single dose
#case     <- factor(c("avg","high"))[2]                       #             |Select the case
#D.o      <- c(avg=1204, high= 250000)[case]      # (micromole)   |oral   dose   #sometimes this messes up... check if neccessary
# @Melina, change the dose here.
dose.O   <- 1                                        # (ug/kgBW/day) |oral   dose, change the dose here.
EoA.O    <- 1
ifelse (multiple==TRUE, n.O <- 3, n.O <- 1)                 #             |number of dosing per day, currently its 1 dose per day.
uptake.O <- EoA.O*dose.O/n.O                             # (nmol)      |amount of uptake
period.O <- 3/60                                            # (h)         |uptake period
koa      <- uptake.O/period.O                               # (nmol/h)    |uptake rate
ifelse (multiple==TRUE,t0 <- c(0, 6, 12) + rep((0:10)*24, each=n.O),t0<-0)
t1 <- t0 + period.O
OD <- data.frame(t0, t1)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compile parameters to be computed  in initialized
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para <- unlist(c(data.frame( 
  
  QCblood,
  QCplasma,
  Qliver,
  Qbrain,
  Qkidney,
  Qgonads,
  Qfat,
  Qrestbody,
  vliver,
  vbrain,
  vkidney,
  vgonads,
  vfat,
  vplasma,
  vrestbody,
  k_liver_plasma,        
  k_brain_plasma,          
  k_kidney_plasma,
  k_gonads_plasma,
  k_fat_plasma,  
  k_restbody_plasma,    
  kgut,
  fu,
  Cl,
  frac
)))
para



########################################################################
#PBPK function model equations 
#########################################################################

PBTKmod <- function(t,y,para) 
  
{
  with (as.list(c(y, para)),
        {
          
          
          ANAUS.O   <- sum(apply(OD, 1, function(x) anaus(t, x[1], x[2])))
          
          Input     <- koa*ANAUS.O                          # Dosing (oral)
          
          cliver = Aliver/vliver
          cbrain = Abrain/vbrain
          ckidney = Akidney/vkidney
          cfat = Afat/vfat
          cgonads = Agonads/vgonads
          crestbody = Arestbody/vrestbody
          cplasma = Aplasma/vplasma

          kurine=Cl*BW
          
          dAgut = -kgut*Agut      + Input*frac             #amount of chemical in gut tissue          
          dAliver = kgut*Agut + Qliver*(cplasma*fu - cliver*(fu/k_liver_plasma)) 
          dAbrain = Qbrain*(cplasma*fu - cbrain*(fu/k_brain_plasma)) 
          dAkidney = Qkidney *(cplasma*fu - ckidney*(fu/k_kidney_plasma))- kurine*ckidney
          dAfat = Qfat *(cplasma*fu - cfat*(fu/k_fat_plasma))                     						 	
          dAgonads = Qgonads *(cplasma*fu - cgonads*(fu/k_gonads_plasma))      					 	# Amount of chemical in gonads
          dArestbody = Qrestbody *(cplasma*fu - crestbody*(fu/k_restbody_plasma))
          dAplasma = (Qliver*cliver*(fu/ k_liver_plasma)) + (Qbrain*cbrain*(fu/k_brain_plasma)) + (Qkidney*ckidney*(fu/k_kidney_plasma))+ 
            (Qfat * cfat*(fu/k_fat_plasma))+ (Qgonads *cgonads*(fu/k_gonads_plasma))+ (Qrestbody *crestbody*(fu/k_restbody_plasma)) - 
            (QCplasma* cplasma*fu) #
          dAurine = kurine*ckidney
    
          dydt = c(dAgut,dAliver,dAbrain, dAkidney,dAfat,dAgonads,dArestbody,dAplasma,dAurine)
          
          conc <- c(cplasma=cplasma, cliver = cliver, cbrain = cbrain, ckidney = ckidney,cfat = cfat,cgonads = cgonads, 
                    crestbody=crestbody)
          
          res  <- list(dydt, conc)
          return(res)
          
        })}

times= seq(0, 1*24*60, 2)/60                    #hour (time)

v<-ode(y=States, func=PBTKmod, times=times, parms=para, method="lsoda")
head(v)
plot(v)

v1=data.frame(v)

par(mfrow=c(2,2))
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))
plot(v1$time,v1$cplasma,type="l", xlab="Time (hours)", ylab="Plasma concentration (µg/L)",pch=19,lwd=2 )
plot(v1$time,v1$cliver,type="l", xlab="Time (hours)", ylab="Liver concentration (µg/L)",pch=19 ,lwd=2)
plot(v1$time,v1$cprostate,type="l", xlab="Time (hours)", ylab="Brain concentration (µg/L)",pch=19,lwd=2 )


