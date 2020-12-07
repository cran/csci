# MPF edit: I could not find where the ssanv and survival packages were used, so I removed 
#      them from the 
#      required package list in the DESCRIPTION file
#Required packages: survival, exactci, and ssanv
#library(survival)
#library(ssanv)
#library(exactci)
#
#CSCI: the Confidence Intervals for Current Status data
#           1) type_CSCI: a) "VALID": the valid confidence interval. If the Lower limit is greater than the Upper, then the Upper and Lower are NPMLE.;
#                         b) "ABA": the approximate binomial appoach confidence interval
#                                  It produces 1) ABA CI and 2) ABA CI (mid-p version). Both of them contain NPMLE.
#                         c) "LIKELIHOOD": the likelihood ratio based confidence interval
#           2) C: the assessment times
#           3) D: indicators (0 or 1)
#           4) times: unobsered event times
#           5) Confidence.level: the confidence level

# MPF edit: I changed the order of the arguments and some names to be more in line with R convensions
#      I added the control argument to allow changing some parameters in the different algorithms 
CSCI<-function(C,D,times=NULL,type=c("VALID","ABA","LIKELIHOOD"),conf.level=0.95, control=controlCSCI()){
		if (any(C<=0) | any(C==Inf)) stop("must have 0<C[i]<Inf for all i")
	   # changed Confidence.level to conf.level to match common R argument
     Confidence.level<- conf.level
     type_CSCI<- match.arg(type)
     if (is.null(times)) times<- sort(unique(C))
	   #INDEX
	   #1) Valid CI
	   #     a) CSCI_Valid function
	   #2) ABA CI
	   #     a) kernel_function_density
	   #     b) kernel_function_derivative
	   #     c) kernel_function_distribution
	   #     d) ghat_function
	   #     e) NPT_function
	   #     f) CSCI_ABA function
	   #3) Likeihood based CI
	   #     a) CSCI_Likelihood function
	   #Selection

#1-a) CSCI_Valid function ##################################################### 
CSCI_Valid<-function(C,D,times,Confidence.level){
	alp<- 1-Confidence.level
  o<-order(C)
  C<- C[o]
  D<- D[o]   

  uTimes<- sort(unique(c(0,Inf,C)))
  k<- length(uTimes)
  n<- length(C)
  uN<-uD<- rep(0,k)
  
  # may be a more efficient way to do this with table?
  for (i in 1:k){
    I<- C==uTimes[i]
    uN[i]<- length(C[I])
    uD[i]<- sum(D[I])
  }
    
  Nat<-Ntb<-Yat<-Ytb<- rep(0,(k-1)) # Note: We need (k-1) intervals.
  # MPF edit: allow the power parameter to be changed in the control function
  #power=2/3
  power<- control$power
  m<- ceiling(n^power)
 
 
  is.even<-function(x){ round(x) %% 2 ==0 }  
  rcumsum<-function(x){ rev(cumsum(rev(x)))   }
    
  
  for (i in 1:(k-1)){
    
    # left end => lower limit 
    AT<- uTimes<=uTimes[i] #less than or equal to t
    if (sum(uN[AT])<=m){
      # at the lower end of F, less than or equal to m assessments 
      # at or before times[i]
      Nat[i]<- sum(uN[AT])
      Yat[i]<- sum(uD[AT])
      } else {  
             rc<- rcumsum(uN[AT])
               # sometimes we can get exactly m values or greater than m values less than or equal to t
                 if (any(rc==m)){
                         I<- rc<=m
                         Nat[i]<- sum(uN[AT][I])
                         Yat[i]<- sum(uD[AT][I])
                          } else {
                          # Let  rc[h] be the largest rc value < m
                          # then we pick values from (h-1):length(rc)
                          # c(1:length(rc))[rc<m] can be interger(0). We set rc_2.
                          rc_2<-c(1:length(rc))[rc<m]
                               if(length(rc_2)==0){
        	                     I<-length(rc)
                                  } else {
                                       h<- min( rc_2 )
                                       I<- (h-1):length(rc)
                                 }
                     Nat[i]<- sum(uN[AT][I])
                     Yat[i]<- sum(uD[AT][I])
                     }
    }
    # right end => upper limit
    TB<- uTimes>uTimes[i] # greater than t
    if (sum(uN[TB])<=m){
      Ntb[i]<- sum(uN[TB])
      Ytb[i]<- sum(uD[TB])
    }  else {
             cc<- cumsum(uN[TB])
             # sometimes we can get exactly m values or greater than m values greater than or equal to t
            if (any(cc==m)){
                           I<- cc<=m
                           Ntb[i]<- sum(uN[TB][I])
                           Ytb[i]<- sum(uD[TB][I])
                           } else {
                           # Let  cc[h] be the largest cc value < m
                           # then we pick values from 1:(h+1)
                           # c(1:length(cc))[cc<m] can be interger(0). We set cc_2.
                           cc_2<-c(1:length(cc))[cc<m]
                                  if(length(cc_2)==0){cc_2=0}        
                                    h<- max( cc_2 )
                                    I<-   1:(h+1)
                       Ntb[i]<- sum(uN[TB][I])
                       Ytb[i]<- sum(uD[TB][I])
                                  }
              }
    
     } 
    qcl=rep(0, (k-1))
    qcu=rep(0, (k-1))
    qcl=qbeta((1-Confidence.level)/2, Yat, Nat-Yat+1)
    qcu=qbeta(1-((1-Confidence.level)/2), Ytb+1, Ntb-Ytb)
  
# If qcl is greater than qcu,
  NPzvalue=isoreg(C,D)$yf
   NPT<-rep(NA, (k-1))
   for(i in 1:(k-1)){
   NPT[i]<-NPT_function(C, NPzvalue, uTimes[i])
   }
     
  II<-qcl>qcu #If the Lower limit is greater than the Upper limit, then the Lower and Upper are NPMLE.
  
 qcl[II]=NPT[II]
 qcu[II]=NPT[II]
 


#Returning all CIs at assessment times (Cis)
l_uTimes<-uTimes[-(length(uTimes))]
u_uTimes<-uTimes[-1]
C_intervals=noquote(paste0("[",round(l_uTimes,4),",", round(u_uTimes,4), ")"))

#MPF edit: add left value of times interval, and NPMLE
#out_lower_upper<-data.frame(C_intervals, qcl,qcu)
#colnames(out_lower_upper)<-c("Intervals","Lower CL", "Upper CL")
out_lower_upper<-data.frame(C_intervals, l_uTimes, NPT, qcl,qcu)
colnames(out_lower_upper)<-c("Intervals","times","NPMLE","Lower CL", "Upper CL")

#Returning CIs at times (tis)
low_ci<-rep(0,length(times))
upp_ci<-rep(0,length(times))
# MPF edit: add NPT_times
NPT_times<- rep(0,length(times))
for(i in 1:length(times)){
	loc_t<-0
    loc_t<-which(l_uTimes<=times[i] & times[i]<u_uTimes)
    # MPF edit: add NPT_times calculation
    NPT_times[i]<- NPT_function(C, NPzvalue, times[i])
    low_ci[i]<-qcl[loc_t]
    upp_ci[i]<-qcu[loc_t]
}
out_times_lower_upper<-data.frame(times, NPT_times, low_ci,upp_ci)
colnames(out_times_lower_upper)<-c("times","NPMLE","Lower CL", "Upper CL")
#MPF edit: name elements of list
out_cis<-list(ciTable_all=out_lower_upper, ciTable_times=out_times_lower_upper) #Show "out_lower_upper" and "out_times_lower_upper".
#out_cis<-out_times_lower_upper #Show only "out_times_lower_upper"
return(out_cis)


}
	
	

#2-a) kernel_function_density ##################################################### 
 kernel_function_density <-
function(type_kernel,u)
# INPUTS:
#   "type_kernel" kernel function: 	"n" Normal, 
#                                   "t" Triweight         
#   "u" array or single value where the kernel is evaluated

{
	
		if(type_kernel == "n")	
		{
		result <- dnorm(u)
		return(result)
		}
		 
		 
	else 
		if(type_kernel == "t")	
		{
		  result <- u
         	Logic0 <- (u <= -1)
          Logic1 <- (u >= 1)
       	  Logic2 <- (u > -1 & u < 1)
       	  Logic3 <- (u > -1 & u < 1)
          result[Logic0] <- 0
	        result[Logic1] <- 0
		  Uval <- result[Logic2]
      result[Logic2] <- (35/32)*((1-(Uval^(2)))^(3))
      return(result)
	  }
}

#2-b) kernel_function_derivative ####################################################
kernel_function_derivative <-
function(type_kernel,u)
# INPUTS:
#   "type_kernel" kernel function: 	"n" Normal, 
#                                   "t" Triweight         
#   "u" array or single value where the kernel is evaluated

{
	
		if(type_kernel == "n")	
		{
		result <- (-u)*dnorm(u)
		return(result)
		}
		 
		 
	else 
		if(type_kernel == "t")	
		{
		  result <- u
         	Logic0 <- (u <= -1)
          Logic1 <- (u >= 1)
       	  Logic2 <- (u > -1 & u < 1)
       	  Logic3 <- (u > -1 & u < 1)
          result[Logic0] <- 0
	        result[Logic1] <- 0
		  Uval <- result[Logic2]
      result[Logic2] <- (35/32)*(-6)*(Uval)*((1-(Uval^(2)))^(2))
      return(result)
	  }
}

#2-c) kernel_function_distribution #####################################################
kernel_function_distribution <-
function(type_kernel,u)
# INPUTS:
#   "type_kernel" kernel function: 	"n" Normal, 
#                                   "t" Triweight         
#   "u" array or single value where the kernel is evaluated

{
	
		if(type_kernel == "n")	
		{
		result <- pnorm(u)
		return(result)
		}
		 
		 
	else 
		if(type_kernel == "t")	
		{
		  result <- u
         	Logic0 <- (u <= -1)
          Logic1 <- (u >= 1)
       	  Logic2 <- (u > -1 & u < 1)
       	  Logic3 <- (u > -1 & u < 1)
          result[Logic0] <- 0
	        result[Logic1] <- 1
		  Uval <- result[Logic2]
      result[Logic2] <- (((35/32) * Uval) -((35/32)* (Uval^3))+((21/32)* (Uval^5)))-((5/32)*(Uval^7)) + 0.5
      return(result)
	  }
}



#2-d) ghat_function #####################################################
ghat_function<-function(x_v, y_v ,TV){
                 
                  kkt=0    
                  R_x_v<-c(-Inf, x_v, Inf)
                  kkt<-(0:length(x_v))[R_x_v[-(length(x_v)+2)]<=TV & TV<R_x_v[-1]] 

                  R_y_v<-c(y_v[1], y_v, y_v[(length(y_v))])
                  R_x_v<-c(x_v[1], x_v, x_v[(length(x_v))])

g_hat<-0
    if(kkt==(length(x_v))){
           g_hat=y_v[(length(y_v))]
       }else if(kkt==0){
           g_hat=y_v[1]
       }else if(kkt!=(length(x_v)) & kkt!=0){
           g_hat=R_y_v[(kkt+1)]+(R_y_v[(kkt+2)]-R_y_v[(kkt+1)])*((TV-R_x_v[(kkt+1)])/(R_x_v[(kkt+2)]-R_x_v[(kkt+1)]))
       }

return(g_hat)
}

#2-e) NPT_function #####################################################
NPT_function<-function(data, npv ,TV){
                  npv<-sort(npv)
                  data<-sort(data)

                  kkt=0    
                  R_data<-c(0, data, Inf)
                  kkt<-(0:length(data))[R_data[-(length(data)+2)]<=TV & TV<R_data[-1]] 

                  R_npv<-c(min(npv), npv, max(npv))
                  R_data<-c(min(data), data, max(data))

NP_hat<-0
    if(kkt==(length(data))){
           NP_hat=max(npv)
       }else if(kkt==0){
           NP_hat=0 # Set 0 
       }else if(kkt!=(length(data)) & kkt!=0){
           NP_hat=R_npv[(kkt+1)]+(R_npv[(kkt+2)]-R_npv[(kkt+1)])*((TV-R_data[(kkt+1)])/(R_data[(kkt+2)]-R_data[(kkt+1)]))
       }

return(NP_hat)
}

#2-f) CSCI_ABA function ##################################################### 
CSCI_ABA<-function(C,D,times,Confidence.level){
  alp<- 1-Confidence.level
  o<-order(C)
  C<- C[o]
  D<- D[o] 
  
  uTimes<- sort(unique(c(0,Inf,C)))
  k<- length(uTimes)
  n<- length(C)
  uN<-uD<- rep(0,k)
 # may be a more efficient way to do this with table?
  for (i in 1:k){
    I<- C==uTimes[i]
    uN[i]<- length(C[I])
    uD[i]<- sum(D[I])
  }
  
  
  
#NPT_function
   NPzvalue=isoreg(C,D)$yf
   NPT<-NULL
   for(i in 1:(k-1)){
   NPT[i]<-NPT_function(C, NPzvalue, uTimes[i])
   }

#Bandwidth
ap=approx(C, NPzvalue, n=512, ties=mean)
dapy=c(ap$y[1], diff(ap$y))
aux<-outer(uTimes[-length(uTimes)],ap$x,"-")

bw=.9*(min(sd(C), (IQR(C)/1.34))*n^(-1/5))#Silverman (1986), Scott (1992), the intitial bandwidth

#Kernel function derivative
aux_dv <-(1/(bw^2))* kernel_function_derivative("n", aux/bw)
fdhat<-NULL
fdhat=aux_dv%*%dapy
fdhat2=fdhat^2
fdhat2[fdhat2<=.001]<-.001 #adjustment for the squard of the first derivative=0

#Kernel function distribution
aux_D <- kernel_function_distribution("n", aux/bw)
Fhat<-NULL
Fhat=aux_D%*%dapy
Fhat_m<-Fhat
Fhat_m[Fhat_m>=.99]<-.99
Fhat_m[Fhat_m<=.01]<-.01 #adjustment for F(t)=zero or 1    
   
#ghat function
    ghat<-NULL
      dTx<-density(C, from=min(C), to=max(C))$x
      dTy<-density(C, from=min(C), to=max(C))$y
      for(i in 1: (k-1)){
                             ghat[i]=ghat_function(dTx, dTy, uTimes[i])
                            }
      ghat[ghat<=.0001]<-.0001 #adjustment g(t)=0
      
#fhat
c_hat=((((Fhat_m*(1-Fhat_m))/ghat)*(1/(2*sqrt(pi))))^(1/5))*((fdhat2)^(-1/5))
h_hat=(c_hat)*n^(-1/5)

Naux<-matrix(NA, nrow=(k-1), ncol=512)
for(i in 1:(k-1)){
Naux[i,]=aux[i,]/(h_hat[i,1])
}

Naux <- kernel_function_distribution("n", Naux)
NFhat=Naux%*%dapy # New F(t)

if((k-1)>=2){
for(i in 1:((k-1)-1)){
if(NFhat[i+1]<=NFhat[i]){
NFhat[i+1]<-NFhat[i]}
}
}

NFhat[NFhat<=.01]<-.01
NFhat[NFhat>=.99]<-.99  #adjustment for F(t)=zero or 1 and monotonicity

Naux<-matrix(NA, nrow=(k-1), ncol=512)
for(i in 1:(k-1)){
Naux[i,]=aux[i,]/(h_hat[i,1])
}
h_hat7=(c_hat)*n^(-1/5)
aux_density_p <-kernel_function_density("n", Naux)
aux_density<-matrix(NA, nrow=(k-1), ncol=512)
for(i in 1:(k-1)){
aux_density[i,]=aux_density_p[i,]/h_hat7[i,1]
}

fhat<-NULL
fhat=aux_density%*%dapy
fhat[fhat<=.0001]<-.0001 #adjustment f(t)=0



#Ratio of f over g.
    ratiofg=fhat/ghat

#dfnew function  
  dfnew<-function(n, ratio.fg, N.Fhat){
   aa=(1/((4*n)^2))*(ratio.fg^2)
   bb=-(1/((4*n)^2))*(ratio.fg^2)
   dd=((N.Fhat)^2-(N.Fhat))
   zz=matrix(c(dd,0,bb,aa), ncol=1)
   as=polyroot(zz)
   Re(as[1])}


qcl=rep(0, (k-1))
qcu=rep(0, (k-1))
qcl_midp=rep(0, (k-1))
qcu_midp=rep(0, (k-1))

#The beginning of q loop ##############################################################

for (q in 1: (k-1)){
   	
#Finding m
m=0
m=dfnew(n, ratiofg[q], NFhat[q])

if(m<=2){m=2}
m=min(m,  (NFhat[q]*(4*n)/ratiofg[q]), ((1-NFhat[q])*(4*n)/ratiofg[q]), n)
if(m<=2){m=2}

N_at_T<-0
N_at_T=uN[uTimes==uTimes[q]]
D_at_T<-0
D_at_T=uD[uTimes==uTimes[q]]
if(length(N_at_T)==0){N_at_T=0}
if(length(D_at_T)==0){D_at_T=0}


AT<- uTimes<=uTimes[q] #less than or equal to
S_Left=sum(uN[AT])
TB<- uTimes>uTimes[q] #greater than
S_Right=sum(uN[TB])


mm<-0

if(m-N_at_T<=0){
	Nab=N_at_T
	Dab=D_at_T
	}else{
		mm<-ceiling(m/2)
        mm_new=min(S_Left, S_Right, mm)
        #mm_new


#kim_add<-function(x,y){
#	if(length(x)!=length(y)){
#		nl<-max(length(x),length(y))
#		length(x)<-nl
#		length(y)<-nl
#		x[is.na(x)]<-0
#		y[is.na(y)]<-0
#	}
#	x+y
#}
#kim_add(rev(uN[AT]),uN[TB])

if(length(uN[AT])>length(uN[TB])){uNTB<-c(uN[TB], rep(0, length(uN[AT])-length(uN[TB])))
	} else {uNTB<-uN[TB]
				}
if(length(uN[AT])<length(uN[TB])){ruNAT<-c(rev(uN[AT]), rep(0, length(uN[TB])-length(uN[AT])))
	} else {ruNAT<-rev(uN[AT])
		}
        		
	    
if(length(uD[AT])>length(uD[TB])){uDTB<-c(uD[TB], rep(0, length(uD[AT])-length(uD[TB])))
	} else {uDTB<-uD[TB]
				}
if(length(uD[AT])<length(uD[TB])){ruDAT<-c(rev(uD[AT]), rep(0, length(uD[TB])-length(uD[AT])))
	} else {ruDAT<-rev(uD[AT])
		}
		
    l_m<-min(which(cumsum(ruNAT)>=mm_new))
	u_m<-min(which(cumsum(uNTB)>=mm_new))
	Nab<-0
    Nab<-cumsum(ruNAT)[l_m]+cumsum(uNTB)[u_m]
    #Nab
    Dab<-0
    Dab<-cumsum(ruDAT)[l_m]+cumsum(uDTB)[u_m]
 }
    
qcl[q]=qbeta((1-Confidence.level)/2, Dab, Nab-Dab+1)
qcu[q]=qbeta(1-((1-Confidence.level)/2), Dab+1, Nab-Dab)
qcl_midp[q]=binom.exact(Dab,Nab,  conf.level=Confidence.level, midp=TRUE)$conf.int[1]
qcu_midp[q]=binom.exact(Dab,Nab,  conf.level=Confidence.level, midp=TRUE)$conf.int[2]

}#end of q loop ###############################################################################

#Edge adjustment ##################################################################################
qcl[1]<-0
qcu[(k-1)]<-1
qcl_midp[1]<-0
qcu_midp[(k-1)]<-1

#If NPMLE is outside of CI, then replace Lower or Upper with NPMLE.#############################
Iqcu<-qcu<NPT
qcu[Iqcu]<-NPT[Iqcu]
Iqcl<-qcl>NPT
qcl[Iqcl]<-NPT[Iqcl]
IqcuM<-qcu_midp<NPT
qcu_midp[IqcuM]<-NPT[IqcuM]
IqclM<-qcl_midp>NPT
qcl_midp[IqclM]<-NPT[IqclM]

#adjustment qcl and qcu ##################################################################
for(i in 1:(k-1)){
	qcl[i]<-max(qcl[1:i])
}

for(i in 1:(k-1)){
	qcu[i]<-min(qcu[i:(length(qcu))])
}

for(i in 1:(k-1)){
	qcl_midp[i]<-max(qcl_midp[1:i])
}

for(i in 1:(k-1)){
	qcu_midp[i]<-min(qcu_midp[i:(length(qcu_midp))])
}

#Returning all CIs at assessment times (Cis) ##############################################
l_uTimes<-uTimes[-(length(uTimes))]
u_uTimes<-uTimes[-1]
C_intervals=noquote(paste0("[",round(l_uTimes,4),",", round(u_uTimes,4), ")"))
#MPF edit: add left end of interval and NPMLE to output
#out_lower_upper<-data.frame(C_intervals, qcl,qcu, qcl_midp,qcu_midp)
#colnames(out_lower_upper)<-c("Intervals","Lower CL", "Upper CL", "midP-Lower CL", "midP-Upper CL")
out_lower_upper<-data.frame(C_intervals, l_uTimes, NPT, qcl,qcu, qcl_midp,qcu_midp)
colnames(out_lower_upper)<-c("Intervals","times", "NPMLE", "Lower CL", "Upper CL", "midP-Lower CL", "midP-Upper CL")



#Returning CIs at times (tis) ##############################################################
low_ci<-rep(0,length(times))
upp_ci<-rep(0,length(times))
low_ci_midp<-rep(0,length(times))
upp_ci_midp<-rep(0,length(times))
# MPF edit: add NPT values at times
NPT_times<- rep(0,length(times))
for(i in 1:length(times)){
	loc_t<-0
    loc_t<-which(l_uTimes<=times[i] & times[i]<u_uTimes)
    low_ci[i]<-qcl[loc_t]
    upp_ci[i]<-qcu[loc_t]
    low_ci_midp[i]<-qcl_midp[loc_t]
    upp_ci_midp[i]<-qcu_midp[loc_t]
    #MPF edit: calculate NPT values at times
    NPT_times[i]<-NPT_function(C, NPzvalue, times[i])
}
# MPF edit: add NPMLE
#out_times_lower_upper<-data.frame(times, low_ci,upp_ci,low_ci_midp,upp_ci_midp)
#colnames(out_times_lower_upper)<-c("times","Lower CL", "Upper CL", "midP-Lower CL", "midP-Upper CL")
out_times_lower_upper<-data.frame(times, NPT_times, low_ci,upp_ci,low_ci_midp,upp_ci_midp)
colnames(out_times_lower_upper)<-c("times","NPMLE","Lower CL", "Upper CL", "midP-Lower CL", "midP-Upper CL")

#MPF edit: name list elements
out_cis<-list(ciTable_all=out_lower_upper, ciTable_times=out_times_lower_upper) #Show "out_lower_upper" and "out_times_lower_upper".
#out_cis<-out_times_lower_upper #Show only "out_times_lower_upper"
return(out_cis)

}


#3-a) CSCI_Likelihood function #####################################################
CSCI_Likelihood<-function(C,D,times,Confidence.level){

  o<-order(C)
  C<- C[o]
  D<- D[o] 
  
  #MPF Question: why do you define uTimes=sort(unique(c(0,Inf,C))) for the other methods
  #            but not for type="LIKELIHOOD" ?  
  uTimes<- sort(unique(C))
  k<- length(uTimes)
  n<- length(C)
  uN<-uD<- rep(0,k)
 # may be a more efficient way to do this with table?
  for (i in 1:k){
    I<- C==uTimes[i]
    uN[i]<- length(C[I])
    uD[i]<- sum(D[I])
  }
  

# MPF edit: allow quan_p and xp_hat to have different input values
#   through the control function
#quan_p<-c(.25,.50,.75,.80,.85,.90,.95,.99)
#xp_hat<-c(.06402, .28506, .80694, .98729, 1.22756, 1.60246, 2.26916, 3.83630)
quan_p<- control$quan_p
xp_hat<- control$xp_hat
  
  
  
d_alpha<-xp_hat[which(quan_p==Confidence.level)]

if (is.na(match(Confidence.level, quan_p))) 
  stop("Choose the Confidence level= .25,.50,.75,.80,.85,.90,.95, or.99. See Table 2 in Banerjee and Wellner (2001).")


#MPF edit: allow different values of intF
#intF=1000
intF<-control$intF
F=c(1:(intF-1)/intF)



isoci<-function(x,y){
    out<-isoreg(x,y)
    cbind(out$x,out$yf)
}


    Ures=0   
        
    NP=isoci(C,D)
    NPF=NP[,2]

    I1<-NPF!=0 & NPF!=1
    mBB1<-D[I1]
    MNPF<-NPF[I1]
    
    Ures=sum(mBB1*log(MNPF)+(1-mBB1)*log(1-MNPF))
  

Lj<-rep(NA, length(times))
Uj<-rep(NA, length(times))


for(j in 1:length(times)){
	
	Res=rep(0, length(F))  
	I_ci<-C<=times[j] 
	kk<-length(C[I_ci])
	
	for(l in 1:length(F)){
 
        if(kk<(n-1)& kk>=2){
            NP1=isoci(C[1:kk], D[1:kk])
            NPF1=NP1[,2]
            NPF1[1:kk][NPF1[1:kk]>F[l]]<-F[l]
            
            NP2=isoci(C[(kk+1):n], D[(kk+1):n])
            NPF2=NP2[,2]
            NPF2[1:(n-kk)][NPF2[1:(n-kk)]<F[l]]<-F[l]
            
            NNPF=c(NPF1,NPF2)

            I<-NNPF!=0 & NNPF!=1
            mBB2<-D[I]
            MNNPF<-NNPF[I]

            Res[l]=sum(mBB2*log(MNNPF)+(1-mBB2)*log(1-MNNPF))

        }else if (kk>=(n-1)){
            NP1=isoci(C,D)
            NPF1=NP1[,2]
            NPF1[NPF1>F[l]]<-F[l]
           
            I<-NPF1!=0 & NPF1!=1
            mBB2<-D[I]
            MNNPF<-NPF1[I]

            Res[l]=sum(mBB2*log(MNNPF)+(1-mBB2)*log(1-MNNPF))

        }else if (kk<2){
          
            NP2=isoci(C,D)
            NPF2=NP2[,2]
            NPF2[NPF2<F[l]]<-F[l]
            
            I<-NPF2!=0 & NPF2!=1
            mBB2<-D[I]
            MNNPF<-NPF2[I]

            Res[l]=sum(mBB2*log(MNNPF)+(1-mBB2)*log(1-MNNPF))
        }
    }
	
	R2R=0
    R2R=2*(Ures-Res)
    RLow=0
    RLow=min(which(R2R<d_alpha))/intF
    RUpp=0
    RUpp=max(which(R2R<d_alpha))/intF
    RLength=RUpp-RLow


Lj[j]<-RLow
Uj[j]<-RUpp

}

# MPF edit: may be very inefficient, but I just copied the NPMLE function
#    and the way the NPMLE is calculated for the times values 
#    from the CSCI_VALID function
NPzvalue=isoreg(C,D)$yf
NPT<-NULL
for(i in 1:(k-1)){
  NPT[i]<-NPT_function(C, NPzvalue, uTimes[i])
}
NPT_times<- rep(0,length(times))
for(i in 1:length(times)){
  NPT_times[i]<- NPT_function(C, NPzvalue, times[i])
}
# MPF edit: added NPMLE to output
#   and made output into the same format as for the other two types
#out_lower_upper<-cbind(times,  Lj,Uj)
#colnames(out_lower_upper)<-c("times","Lower CL", "Upper CL")
#return(out_lower_upper)
out_lower_upper<- data.frame(times, NPT_times,  Lj,Uj)
colnames(out_lower_upper)<-c("times","NPMLE","Lower CL", "Upper CL")
out_cis<-list(ciTable_all=NULL, ciTable_times=out_lower_upper) 
return(out_cis)

}

#Select one of VALID, ABA, or LIKELIHOOD.#####################################################
	   
    if (toupper(type_CSCI)=="VALID"){   	
    	return(CSCI_Valid(C,D,times,Confidence.level))  
    	} else if (toupper(type_CSCI)=="ABA"){ 
    		return(CSCI_ABA(C,D,times,Confidence.level))
    		} else if (toupper(type_CSCI)=="LIKELIHOOD"){ 
    		return(CSCI_Likelihood(C,D,times,Confidence.level))
    		}    
  

}
############################################################################################### 
#the end of the CSCI function.#################################################################
###############################################################################################
# MPF edits: add controlCSCI function
# controlCSCI=function that gives control parameters for algorithms for the 
# MPF Question: Do we want to add any control parameters for the kernel functions? 
#               Perhaps do that at a latter version of the software
controlCSCI<-function(power=2/3,
                      quan_p=c(.25,.50,.75,.80,.85,.90,.95,.99),
                      xp_hat=c(.06402, .28506, .80694, .98729, 1.22756, 1.60246, 2.26916, 3.83630),
                      intF=1000){
  if (power<0 | power>1) stop("power must be in (0,1)")
  # get quan_p and xp_hat from Table 2 in Banerjee and Wellner (2001)
  if (any(quan_p<0 | quan_p>1)) stop("quan_p must be in (0,1)")
  if (intF<3) stop("intF<3")
  list(power=power, quan_p=quan_p, xp_hat=xp_hat,intF=intF)
}
#data(hepABulg)
#library(exactci)
#v<-CSCI(C=hepABulg$age,D=hepABulg$testPos,type="VALID")
#v<-CSCI(C=hepABulg$age,D=hepABulg$testPos,times=c(10,20,50),type="VALID")
#v<-CSCI(C=hepABulg$age,D=hepABulg$testPos,times=c(10,20,50),type="ABA")
#v<-CSCI(C=hepABulg$age,D=hepABulg$testPos,times=c(10,20,50),type="LIKELIHOOD")
#str(v)
#v
