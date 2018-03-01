cat("\014")
rm(list=ls())
dev.off()
library("optimx")
library("GenSA")
library("survival")
library("Hmisc")
library("rootSolve")
library("stats4")
library("VGAM")
library("lattice")
library("splines")
library("zoo")
library("forecast")
library("pracma")
library("base")


Shoreline <- read.csv("D:\\YEN\\NhaTrang data\\Data_Rafael_breakwave_SWAN\\NT_Shoreline_Rafael.csv")
Forcing <- read.csv("D:\\YEN\\NhaTrang data\\Data_Rafael_breakwave_SWAN\\Hs_Tp_Dir_breakwave_Jan2013_July2014_delT_1h.csv")

gamma=0.78
g=9.81
rho=1000

date_start <- "01/08/2013 00:00:00"
date_end <- "08/07/2014 00:00:00"

##SHORELINE##
##################################################################################################################
time_S <- ISOdate(paste( Shoreline$year), paste(Shoreline$month), paste(Shoreline$day), paste(Shoreline$hour),0,0)
timeS <- as.factor(format(time_S,'%d/%m/%Y %H:%M:%S'))
i0 <- which(timeS==date_start)
i0f <- which(timeS==date_end)
# dt_timeS_hour <- as.numeric(diff(time_S[i0:i0f]))   #hours
# dt_timeS <- dt_timeS_hour/24   #days
dt_timeS <- as.numeric(diff(time_S[i0:i0f]))

Shoreline_P1 <- Shoreline[,5]

XS_S <- Shoreline_P1[i0:i0f]    ## Shoreline of interest (shorline of considered time-series)
XS_S <- -XS_S 

# plot(time_S[i0:i0f],XS_S,type="l", ylim=range(15:30), xlab="years",ylab="S (m)", col="red", cex.lab = 1.7, cex.axis = 1.7 )


# ##    POINT 5 _ ECMWF   ##
# ##################################################################################################################
# timeF_P5 <- ISOdate(paste(P5_ECMWF$year), paste(P5_ECMWF$month), paste(P5_ECMWF$day), paste(P5_ECMWF$hour),0,0)
# time_P5 <- as.factor(format(timeF_P5,'%d/%m/%Y %H:%M:%S'))
# i0_P5 <- which(time_P5==date_start)
# i0f_P5 <- which(time_P5==date_end)
# vec_ind_P5 <- i0_P5:i0f_P5
# 
# Hs_P5 <- P5_ECMWF[,5]
# Tp_P5 <- P5_ECMWF[,6]
# Dir_P5 <- P5_ECMWF[,7]


##FORCING##
##################################################################################################################

time_F <- ISOdate(paste(Forcing$year), paste(Forcing$month), paste(Forcing$day), paste(Forcing$hour),0,0)
timeF <- as.factor(format(time_F,'%d/%m/%Y %H:%M:%S'))
ii0 <- which(timeF==date_start)
iif <- which(timeF==date_end)
vec_ind <- ii0:iif

Hs_F <- Forcing[,5]
Hs_b <- Forcing[,8]
Hs_b_Snell <- Forcing[,11]
Tp_F <- Forcing[,7]
Dir_10 <- Forcing[,6]
dir_10 <- Dir_10 - 90
Dir_F <- Forcing[,9]
Dir_P5 <- Forcing[,10]
Dirb_Snell <- Forcing[,12]


dir_F <- Dir_F - 90

Hs <- Hs_F[vec_ind]
Hsb <- Hs_b[vec_ind]
Hsb_Snell <- Hs_b_Snell[vec_ind]
Tp <- Tp_F[vec_ind]
Dir <- Dir_F[vec_ind]
dir <- dir_F[vec_ind]

mean_Dir_F <- mean(Dir_F[ Dir_F != 0 ])
mean_dir_F <- mean(dir_F[ Dir_F != 0 ])

mean_Dir_10 <- mean(Dir_10[ Dir_F != 0 ])
mean_dir_10 <- mean(dir_10[ Dir_F != 0 ])

mean_Dirb_Snell <- mean(Dirb_Snell[ Dir_F != 0 ])

# ##################################################################################################################
# Dir_F <- Dir_10
# mean_Dir_F <- mean_Dir_10
# ##################################################################################################################

##################################################################################################################
# Dir_F <- Dirb_Snell
# mean_Dir_F <- mean_Dirb_Snell
# Hs_b <- Hs_b_Snell
##################################################################################################################


plot(time_F[vec_ind],Hs_F[vec_ind],type="l", xlab="years",ylab="Hs (m)",ylim=range(0:2), cex.lab = 1.7, cex.axis = 1.7 )
par(new=TRUE)
plot(time_F[vec_ind],Hs_b[vec_ind],type="l", xlab="",ylab="", col="blue",ylim=range(0:2),cex.lab = 1.7, cex.axis = 1.7 )
legend("topleft", legend = c("Hs_10m","Hs_b"),lty=c(1,1),col=c("black","blue"),bty = "n",cex = 1)

plot(time_F[vec_ind],Tp_F[vec_ind],type="l", xlab="years",ylab="Tp (s)", cex.lab = 1.7, cex.axis = 1.7 )
# par(new=TRUE)
# plot(timeF_P5[vec_ind_P5],Tp_P5[vec_ind_P5],type="l", xlab="",ylab="", col="blue",ylim=range(0:3.5),cex.lab = 1.7, cex.axis = 1.7 )


plot(time_F[vec_ind],Dir_F[vec_ind],type="l", xlab="years",ylab="Dir",ylim=range(0:360), cex.lab = 1.7, cex.axis = 1.7 )
par(new=TRUE)
plot(time_F[vec_ind],Dir_10[vec_ind],type="l", xlab="",ylab="", col="blue",ylim=range(0:360),cex.lab = 1.7, cex.axis = 1.7 )
par(new=TRUE)
plot(time_F[vec_ind],Dir_P5[vec_ind],type="l", xlab="",ylab="", col="red",ylim=range(0:360),cex.lab = 1.7, cex.axis = 1.7 )
legend("topleft", legend = c("Dir_b","Dir_10","Dir_P5"),lty=c(1,1),col=c("black","blue","red"),bty = "n",cex = 1)


ia <- match(timeS,timeF)    # vector of indices element index of forcing data corresponding w/ shoreline data

dt_timeS_n <- dt_timeS

XS <- XS_S
XSb <- mean(XS)
XS_0 <- XS - XSb
XS <- XS_0
dXS <- diff(XS)  
dXS_dt <- dXS/dt_timeS_n    ## Shoreline change rate (m/d)

index_not0 <- which(!Dir_F == 0)
fluc <- Dir_F - mean_Dir_F
mean_fluc <- mean(fluc[index_not0])
fluc_10 <- Dir_10 - mean_Dir_10
mean_fluc_10 <- mean(fluc_10)
Q <- 0.1*(Hs_b**2.5)*sin(2*fluc*pi/180)
sum_Q <- sum(Q)
corr1 <- mean((Hs_b**2.5)*cos(2*fluc*pi/180))
corr2 <- mean((Hs_b**2.5)*sin(2*fluc*pi/180))
tan <- corr2/corr1
psi <- atan(corr2/corr1)    ## radian
psi_deg <- psi*180/pi    ## degree
print(tan)
print(psi_deg)


plot(time_F[index_not0],fluc[index_not0],type="l", xlab="years",ylab=expression(alpha~fluc) ,xlim=range(time_F),ylim=range(-100:300), cex.lab = 1.7, cex.axis = 1.7 )
par(new=TRUE)
plot(time_F,fluc_10,type="l", xlab="",ylab="",xlim=range(time_F),ylim=range(-100:300),col='red', cex.lab = 1.7, cex.axis = 1.7 )

plot(time_F[vec_ind],fluc_10[vec_ind],type="l",  xlab="years",ylab=expression(alpha~fluc),xlim=range(time_F[vec_ind]),ylim=range(-100:300),col='red', cex.lab = 1.7, cex.axis = 1.7 )

plot(time_F[vec_ind],fluc[vec_ind],pch=1,  xlab="years",ylab=expression(alpha~fluc),xlim=range(time_F[vec_ind]),ylim=range(-5:10),col='red', cex.lab = 1.7, cex.axis = 1.7 )


k <- matrix(0,nrow=length(timeF), ncol=1)
h <- 16

Airy <- function(h,Tp){
  wg<-2*pi/Tp
  Airy0<-wg*wg/g
  #print(Airy0)
  kkit<-Airy0
  for(it in 1:500){
    kkit<-Airy0/tanh(kkit*h)
  }
  return(kkit)
}

for(ik in 1:length(timeF)){
  k[ik]<-Airy(16,Tp_F[ik])
}

c <- sqrt(g*tanh(h*k)/(k))
Cg <- (c/2)*(1+(2*h*k)/sinh(2*h*k))


Pow <- Hs_F*Hs_F/16*rho*g*Cg
Hs_F_3_2 <- Hs_b**1.5
Hs_F_3_2_cos <- (Hs_b**1.5)*cos(2*fluc*pi/180)
Hs_F_3_2_sin <- (Hs_b**1.5)*sin(2*fluc*pi/180)

total <- Hs_F_3_2_cos + Hs_F_3_2_sin

plot(time_F[vec_ind],Hs_F_3_2_cos[vec_ind],type="l", xlab="years",ylab="Hs_F_3_2_cos", cex.lab = 1.7, cex.axis = 1.7 )

plot(time_F[vec_ind],Hs_F_3_2_sin[vec_ind],type="l", xlab="years",ylab="Hs_F_3_2_sin", cex.lab = 1.7, cex.axis = 1.7 )

plot(time_F[vec_ind],total[vec_ind],type="l", xlab="years",ylab="", cex.lab = 1.7, cex.axis = 1.7 )

## Calculate sediment fall velocity
del <- 1.65
g <- 9.81
nu <- 10^(-6)
d50 <- 0.0004
d_st <- ((del*g/nu^2)^(1/3))*d50
w_s <- nu/d50*(sqrt(25+1.2*d_st^2)-5)^1.5
# w_s=0.037

## Calculate Dean parameter
omega_F <- Hs_F/(Tp_F*w_s)    # Dean parameter of all forcing time series 
omega <- Hs/(Tp*w_s)          # Dean parameter of considered time series 



################################################################################################################
## Calculate Equilibrium Dean parameter time series for every hour
phi_all <- c(5,10,15,20,25,40,50)          
phi_all <- 20

for (i_phi in 1:length(phi_all))
{
  phi <- phi_all[i_phi]    ## memory decay in days
  
  n_D <- 24
  phi_h=phi*n_D
  DD=2*phi_h  # number of values in the computation of omega_eq
  qq=10^(-1/phi_h)
  deno=(1-qq^(DD+1))/(1-qq)
  
  omega_eqlbrm <- matrix(0,nrow=length(timeF), ncol=1)
  
  Vqq <- matrix(1,nrow=(DD+1), ncol=1)
  for(ii in 0:(DD)){ Vqq[ii+1]=qq^ii}
  Vqq=rev(Vqq)
  
  for (kk in (DD+1):length(timeF))
  {range=(kk-DD):kk
  omega_eqlbrm[kk] = omega_F[range]%*%Vqq}
  omega_eqlbrm=omega_eqlbrm/deno
  
  Domega=omega_eqlbrm-omega_F
  sigma_omega <- sqrt(mean((Domega-mean(Domega))^2))
  
  plot(time_F[vec_ind],omega,type="l", xlab="years",ylab="omega",ylim=range(0,4), xaxt='n', cex.lab = 1.7, cex.axis = 1.7 )
  par(new=TRUE)
  plot(time_F[vec_ind],omega_eqlbrm[vec_ind],type="l", xlab="",ylab="",ylim=range(0,4), xaxt='n',col="red", cex.lab = 1.7, cex.axis = 1.7 )
  axis.POSIXct(1, at = as.Date(time_F[vec_ind]), format= "%m-%Y", labels = TRUE)
  legend("topleft", legend = paste("phi =", round(phi, 1)), bty = "n",cex = 1)
  
  plot(time_F,omega_F,type="l", xlab="years",ylab="omega", xaxt='n',ylim=range(0,4), cex.lab = 1.7, cex.axis = 1.7 )
  par(new=TRUE)
  plot(time_F,omega_eqlbrm,type="l", xlab="years",ylab="omega", xaxt='n',ylim=range(0,4),col="red", cex.lab = 1.7, cex.axis = 1.7 )
  axis.POSIXct(1, at = as.Date(time_F), format= "%m-%Y", labels = TRUE)
  legend("topleft", legend = paste("phi =", round(phi, 1)), bty = "n",cex = 1)
  
  ## Calculate the 2nd term of the formulation for every hour [ FF = (omega_eqlbrm - omega)P^0.5 ]
  
  FF <- matrix(0,nrow=length(timeF), ncol=1)
  FF_pl <- matrix(0,nrow=length(timeF), ncol=1)
  FF_mi <- matrix(0,nrow=length(timeF), ncol=1) 
  
  ######################################################################################################################
  ## Calculate Mean value of FF_plus and FF_minus for every survey time interval
  
  FF <- (Domega)*(sqrt(Pow))/sigma_omega
  # mean_FF <- mean(FF)
  FF_pl <- FF
  FF_mi <- FF
  for (kk in 1:length(timeF))
  {
    if (FF_pl[kk] < 0)
    {FF_pl[kk] = 0}
    if (FF_mi[kk] > 0)
    {FF_mi[kk] = 0}
  }
  
  indr=ia[1]:ia[length(ia)]
  r <- abs(sum(FF_pl[indr])/sum(FF_mi[indr]))
  

  m_Hs_F_3_2_cos <- matrix(0,nrow=length(dt_timeS_n), ncol=1)
  m_Hs_F_3_2_sin <- matrix(0,nrow=length(dt_timeS_n), ncol=1)
  mFF_pl <- matrix(0,nrow=length(dt_timeS_n), ncol=1)
  mFF_mi <- matrix(0,nrow=length(dt_timeS_n), ncol=1)
  
  
  for (k in 1:length(dt_timeS_n))
  {
    m <- ia[k]
    n <- ia[(k+1)]
    m_Hs_F_3_2_cos[k] <- sum(Hs_F_3_2_cos[(m+1):n])/(n-m)*dt_timeS_n[k]
    m_Hs_F_3_2_sin[k] <- sum(Hs_F_3_2_sin[(m+1):n])/(n-m)*dt_timeS_n[k]
    mFF_pl[k] <- sum(FF_pl[(m+1):n])/(n-m)*dt_timeS_n[k]
    mFF_mi[k] <- sum(FF_mi[(m+1):n])/(n-m)*dt_timeS_n[k]
    
  }
  
  
  # r <- abs(sum(mFF_pl)/sum(mFF_mi))
  # r <- 0.55
  
  Ds <- matrix(nrow=length(dt_timeS_n), ncol=1)
  S_calib <- matrix(0,nrow=length(dt_timeS_n)-1, ncol=1)
  S_calib <- rbind(XS[1],S_calib)
  
  rr <- function(a){
    
    for (k in 1:length(dt_timeS_n))
    {
      
      Ds[k] <- (a[1]*(r*mFF_mi[k] + mFF_pl[k]) + a[2]*m_Hs_F_3_2_cos[k] + a[3]*m_Hs_F_3_2_sin[k])
      
      S_calib[k+1] <- S_calib[k] + Ds[k] 
    }
    # avrg_S_calib <- mean(S_calib)
    # del_meanS <- XSb - avrg_S_calib
    # S_calib <- S_calib + del_meanS
    
    rr <- sqrt(mean((XS-S_calib)^2))}
  
  op <- GenSA(fn=rr, lower=c(-10,-10,-10), upper=c(10,10,10))
  
  aa <-op$par
  
  print(aa)
  
  for (k in 1:length(dt_timeS_n))
  {
    Ds[k] <-  (aa[1]*(r*mFF_mi[k] + mFF_pl[k]) + aa[2]*m_Hs_F_3_2_cos[k] + aa[3]*m_Hs_F_3_2_sin[k])
    
    S_calib[k+1] <- S_calib[k] + Ds[k] 
  }
  
  # avrg_S_calib <- mean(S_calib)
  # del_meanS <- XSb - avrg_S_calib
  # S_calib <- S_calib + del_meanS
  
  mean_S_calib <- mean(S_calib)
  d_meanS <- XSb - mean_S_calib
  RMSE <- sqrt(mean((XS-S_calib)^2))
  NMSE <- sum((XS-S_calib)^2)/sum((XS-XSb)^2)
  
  plot(time_S[i0:i0f],XS,type="l",ylim=rev(range(-10:6)), xlab="years",ylab="Shoreline position (m)", col="red", cex.lab = 1.7, cex.axis = 1.7 )
  par(new=TRUE)
  lines(time_S[i0:i0f],S_calib,ylim=rev(range(-10:6)), xlab="",ylab="",lwd=3, col="black", cex.lab = 1.7, cex.axis = 1.7 )
  mtext("P1",cex = 1.7)
  legend("topleft", legend = c(
    paste("phi =", round(phi, 1)),
    paste("c_cs =",round(aa[1],4)),
    paste("a =",round(aa[2],4)), 
    paste("b =",round(aa[3],4)), 
    paste("r =",round(r,4)),
    paste("RMSE =",round(RMSE,4)),
    paste("NMSE =",round(NMSE,4)),
    paste("d_meanS_m =",round(d_meanS,4))),
    bty = "n",cex = 1)
  
  plot(time_S[i0:i0f],XS,pch=1,ylim=rev(range(-10:6)), xlab="2013-2014",ylab="S (m)", col="red", cex.lab = 1.2, cex.axis = 1.2 )
  par(new=TRUE)
  lines(time_S[i0:i0f],S_calib,ylim=rev(range(-10:6)), xlab="",ylab="",lwd=3, col="black", cex.lab = 1.2, cex.axis = 1.2 )
  
  
  ##CROSS-SHORE CONTRIBUTION
  ######################################################################################################################
  Ds_cr <- matrix(nrow=length(dt_timeS_n), ncol=1)
  S_calib_cr <- matrix(0,nrow=length(dt_timeS_n)-1, ncol=1)
  S_calib_cr <- rbind(XS[1],S_calib_cr)
  
  cr <- c(aa[1],0,0)
  
  for (k_cr in 1:length(dt_timeS_n))
  {
    Ds_cr[k_cr] <-  aa[1]*(r*mFF_mi[k_cr] + mFF_pl[k_cr]) 
    
    S_calib_cr[k_cr+1] <- S_calib_cr[k_cr] + Ds_cr[k_cr] 
  }
  plot(time_S[i0:i0f],XS,type="l",ylim=rev(range(-10:6)), xlab="years",ylab="Shoreline position (m)", col="red", cex.lab = 1.7, cex.axis = 1.7 )
  par(new=TRUE)
  lines(time_S[i0:i0f],S_calib_cr,ylim=rev(range(-10:6)), xlab="",ylab="",lwd=3, col="black", cex.lab = 1.7, cex.axis = 1.7 )
  mtext("P1_cross_shore",cex = 1.7)
  legend("topleft", legend = c(
    paste("phi =", round(phi, 1)),
    paste("c_cs =",round(aa[1],4)), 
    paste("r =",round(r,4))),
    bty = "n",cex = 1)
  
  plot(time_S[i0:i0f],XS,pch=1,ylim=rev(range(-10:6)), xlab="2013-2014",ylab="S (m)", col="red", cex.lab = 1.2, cex.axis = 1.2 )
  par(new=TRUE)
  lines(time_S[i0:i0f],S_calib_cr,ylim=rev(range(-10:6)), xlab="",ylab="",lwd=3, col="black", cex.lab = 1.2, cex.axis = 1.2 )
  
  
  
  ##LONGSHORE CONTRIBUTION
  ######################################################################################################################
  Ds_lo <- matrix(nrow=length(dt_timeS_n), ncol=1)
  S_calib_lo <- matrix(0,nrow=length(dt_timeS_n)-1, ncol=1)
  S_calib_lo <- rbind(XS[1],S_calib_lo)
  
  lo <- c(0,aa[2],aa[3])
  
  for (k_lo in 1:length(dt_timeS_n))
  {
    Ds_lo[k_lo] <-  aa[2]*m_Hs_F_3_2_cos[k_lo] + aa[3]*m_Hs_F_3_2_sin[k_lo]
    
    S_calib_lo[k_lo+1] <- S_calib_lo[k_lo] + Ds_lo[k_lo] 
  }
  plot(time_S[i0:i0f],XS,type="l",ylim=rev(range(-10:6)), xlab="years",ylab="Shoreline position (m)", col="red", cex.lab = 1.7, cex.axis = 1.7 )
  par(new=TRUE)
  lines(time_S[i0:i0f],S_calib_lo,ylim=rev(range(-10:6)), xlab="",ylab="",lwd=3, col="black", cex.lab = 1.7, cex.axis = 1.7 )
  mtext("P1_long_shore",cex = 1.7)
  legend("topleft", legend = c(
    paste("phi =", round(phi, 1)),
    paste("a =",round(aa[2],4)), 
    paste("b =",round(aa[3],4))),
    bty = "n",cex = 1)
  
  
  plot(time_S[i0:i0f],XS,pch=1,ylim=rev(range(-10:6)), xlab="2013-2014",ylab="S (m)", col="red", cex.lab = 1.2, cex.axis = 1.2 )
  par(new=TRUE)
  lines(time_S[i0:i0f],S_calib_lo,ylim=rev(range(-10:6)), xlab="",ylab="",lwd=3, col="black", cex.lab = 1.2, cex.axis = 1.2 )
  
  
  
  ####################################################################################################################
  ####################################################################################################################
  ## Calculate dsdt with constant time interval for entire series
  
  
  tm <- 1/24          # chosen time interval (days)
  dt_tm <- rep.int(tm,sum(dt_timeS)/tm)
  im <- seq(ia[1],ia[length(ia)],tm*n_D)
  
  
  m_Hs_F_3_2_cos_m <- matrix(0,nrow=length(dt_tm), ncol=1)
  m_Hs_F_3_2_sin_m <- matrix(0,nrow=length(dt_tm), ncol=1)
  mFF_pl_m <- matrix(0,nrow=length(dt_tm), ncol=1)
  mFF_mi_m <- matrix(0,nrow=length(dt_tm), ncol=1)
  
  
  for (k_m in 1:length(dt_tm))
  {
    m_m <- im[k_m]
    n_m <- im[(k_m+1)]
    m_Hs_F_3_2_cos_m[k_m] <- sum(Hs_F_3_2_cos[(m_m+1):n_m])/(n_m-m_m)*dt_tm[k_m]
    m_Hs_F_3_2_sin_m[k_m] <- sum(Hs_F_3_2_sin[(m_m+1):n_m])/(n_m-m_m)*dt_tm[k_m]
    mFF_pl_m[k_m] <- sum(FF_pl[(m_m+1):n_m])/(n_m-m_m)*dt_tm[k_m]
    mFF_mi_m[k_m] <- sum(FF_mi[(m_m+1):n_m])/(n_m-m_m)*dt_tm[k_m]
    
  }
  
  
  Ds_m <- matrix(nrow=length(dt_tm), ncol=1)
  S_calib_m <- matrix(0,nrow=length(dt_tm)-1, ncol=1)
  S_calib_m <- rbind(XS[1],S_calib_m)
  
  for (k_m in 1:length(dt_tm))
  {
    Ds_m[k_m] <-  (aa[1]*(r*mFF_mi_m[k_m] + mFF_pl_m[k_m]) + aa[2]*m_Hs_F_3_2_cos_m[k_m] + aa[3]*m_Hs_F_3_2_sin_m[k_m])
    
    S_calib_m[k_m+1] <- S_calib_m[k_m] + Ds_m[k_m]
  }
  
  # time_n_i <- seq(ISOdate(2013,5,25,12), by = "6 hour", length.out = length(S_calib_m))
  time_n_i <- seq(ISOdate(2013,8,1,0), by = "1 hour", length.out = length(S_calib_m))
  
  ## Calculate RMSE
  ttS <- c(0,cumsum(dt_timeS)*24)
  XS_intp <- approx(ttS,XS,method="linear",xout = seq(0,length(time_n_i)-1, by = 1))$y
  RMSE_tm <- sqrt(mean((XS_intp-S_calib_m)^2))
  # NMSE_tm <- sum((XS[1:length(S_m_tS)]-S_m_tS)^2)/sum((XS_m-XSb_m)^2)
  # print(RMSE_tm)
  
  plot(time_S[i0:i0f],XS,type="l",ylim=rev(range(-10:6)), xlab="years",ylab="Shoreline position (m)", col="red", cex.lab = 1.7, cex.axis = 1.7 )
  par(new=TRUE)
  lines(time_n_i,S_calib_m,ylim=rev(range(-10:6)), xlab="",ylab="",lwd=3, col="black", cex.lab = 1.7, cex.axis = 1.7 )
  mtext("P1",cex = 1.7)
  legend("topleft", legend = c(
    paste("phi =", round(phi, 1)),
    paste("c_cs =",round(aa[1],4)),
    paste("a =",round(aa[2],4)),
    paste("b =",round(aa[3],4)),
    paste("r =",round(r,4)),
    paste("RMSE_tm =",round(RMSE_tm,4)),
    # paste("NMSE =",round(NMSE,4)),
    paste("d_meanS_m =",round(d_meanS,4))),
    bty = "n",cex = 1)
  
  
}
