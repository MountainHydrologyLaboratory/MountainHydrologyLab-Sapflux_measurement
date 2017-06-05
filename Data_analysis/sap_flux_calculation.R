# SAP Flux calculation from Heat pulse method, based on McJannet et al. 2004
# R Code By Rodrigo Fernandez and Luis Andres Guillen
# 06/01/2017


# Prepare workspace}
#Clear workspace and read variables}

setwd("~/Desktop/Data_analysis/") #set working directory
rm(list=ls())     #clear workspace
inmat<- as.matrix(read.csv("mc_visit6.csv", header=FALSE)[47:57278,6:23]) #reads datafile
n<-nrow(inmat) #read the sample size
trees<-ncol(inmat)/3 #number of trees
## Declare Datalogger Variables

ms<-98 #number of measurements per pulse
sec<-2 # seconds per interval
isec<-6 # initial seconds without reading

#### Declare other variables
pul <-n/ms
ndays<-pul/48
stvar<- matrix(data=-9999,nrow=pul,ncol=trees) #should be NA to save space
temp<- matrix(data=-9999,nrow=pul,ncol=trees)
w<-1 #wound width (mm) of the borehole
d<-2.5 #determine the spacing between the heater and the middle of the two sensors
vwd<-1 #volume fraction of wood
vwt<-1 #volume fraction of water
error<- -9999


for (t in 1:trees){ #tree loop
  
  df <- as.numeric(inmat[,t*3])# read current tree sensor differences
  pos<-1 #position counter
  
  for (p in 1:pul) { #pulse loop
    
    pre<- df[pos]   # temperature difference at the beginning of the measurement
    con<-0 # condition
    temp[p,t]<-as.numeric(inmat[pos,t*3-2])
    
    for (s in 10:ms){ #time step loop
      
      post<-df[pos+s-1] #temperature post pulse
    
      if (con==0) { #condition
        
        if (post<=pre){ # condition to test of post temperature is back to pre temperature
          
          stvar[p,t] <- s*sec + isec #store the number of seconds for post=pre
         
           con<-1 #change condition
          
        } #end of condition
        
      } #end of condition
      
      
    } #end of step loop
    
    pos <- pos +ms #add number of time steps to position counter
    
  }# end of pulse loop
  
}# end of tree loop


vh<-matrix(data = -9999,nrow=pul,ncol=6)

#Calculate heat pulse velocity}

for (t in 1:6){ #tree loop
  
  for(p in 1:pul){#pulse cycle 
    
    if (stvar[p,t]==-9999) { #condition
   
     vh[p,t]<-error
   
    } else {
    
  
    vh[p,t]<-d/stvar[p,t]*360 #Equation 1 from McJannet etal 2004: initialize vector to store heat pulse velocity (vh)
   
     }#end of condition
    
  }# end loop pulse
  
  
}# end loop tree

vc<-matrix(data=-9999,nrow=pul,ncol=6) #initialize vector to store corrected vh for wounding effects (vc)


#Correct heat pulse for wounding effect

for (t in 1:6){ #tree loop
  
  for(p in 1:pul){#pulse cycle 
    
    if (vh[p,t]>3){ #conditional for vh in relation to 3cms/h (>)
      
      b1<--0.1175*w^2+1.46*w-1.6432 #Equation 3 McJannet 2004
      b2<-0.721*w^2-0.644*w+2.2024 #Equation 4 McJannet 2004
      b3<-0.0239*w^2-0.0319*w+0.0259 #Equation 5 McJannet 2004
      
    } else if(vh[p,t]>0) { #conditional for vh in relation to 3cms/h (<)
      
      b1<-0.0259*w^2-0.0397*w-0.4409 #Equation 6 McJannet 2004
      b2<-0.081*w^2-0.1298*w+1.3316 #Equation 7 McJannet 2004
      b3<-0.051*w^2-0.0003*w+0.0166 #Equation 8 McJannet 2004
      
    } else {
      
      b1<-error
      b2<-0
      b3<-0
          
    } #end of conditional
    
    vc[p,t]<-b1+b2*vh[p,t]+b3*(vh[p,t])^2  #Equation 2 McJannet 2004
    
  }# end loop pulse
  
  
}# end loop tree



#Converting heat flux to sap velocity


#### NEEDS TO BE REVISED FOR A TEMPERATURE CYCLE ACCORDING TO TEMPERATURE INPUT
k<-0.4+0.00214*temp-0.000006*temp^2 #Equation 10 McJannet 2004
vl<-matrix(data=-9999,nrow = pul, ncol = 6)
for (t in 1:6){
  for (p in 1:pul){
    if (vh[p,t]>0){
      vl[p,t]<-vc[p,t]*(k[p,t]*vwd+vwt) #Equation 9 McJannet 2004
    }
  }
}

write.csv(vl,"sapveloity.csv")


