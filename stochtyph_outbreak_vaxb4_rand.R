
stochtyph_identify = function(tmax=4173, dt=.1, vax=FALSE, catchup=FALSE, threshIR=1000, threshSD=9, t_to_vax=NA,t_outbreak=t_outbreak){
  
  tf= tmax/dt ### number of timesteps to simulate
  
  t_diff <- t2-t1
  
  t1=t_outbreak
  t2=t_outbreak+t_diff
  
  veff    <- sample(veff.samps,1) #.875
  vwane   <- sample(omega_v.samps,1)
  vcov_C  <- runif(1,.6,.9) #.6
  rep.adj <- sample(adj_fact.samps,1)
  
  tvacc       <- round(2400/dt + t_to_vax)  #round(t1 - t0 + t_to_vax) 
  v           <- matrix(0,nrow=tspan, ncol=al) 
  massvacc    <- matrix(0,nrow=tspan, ncol=al)
  
  if (vax==TRUE) {
    v[tvacc:tf,agevacc] <- vcov_R.t[1:length(tvacc:tf)]
    
    if (catchup==TRUE) {massvacc[(tvacc):(tvacc+(4/dt)),agevacc:4] <- min((1-(1-vcov_C)^.25),.625)}
    
  }

  #intitialize empty vectors/matrices for each compartment
  S1 = I1 = R = S2 = I2 = C = N = V1 = V2 = matrix(0, nrow=tf, ncol=al)
  W = D = D.adj = outbreak_tstep = rep(0,tf)
  
  #initialize empty vectors/matrices for transitions
  W_enter  = W_leave = rep(0,tf)
  births = deaths_S1 = deaths_I1 = deaths_R = deaths_S2 = deaths_I2 = deaths_C = deaths_V1 = deaths_V2 = 
    age_out_S1 = age_out_I1 = age_out_R = age_out_S2 = age_out_I2 = age_out_C = age_out_V1 = age_out_V2 =
    age_in_S1 = age_in_I1 = age_in_R = age_in_S2 = age_in_I2 = age_in_C = age_in_V1 = age_in_V2 =
    S1_I1 = I1_R = I1_C = R_S2 = S2_S1 = S2_I2 = I2_R = I2_C = S1_V1 = V1_S1 = S2_V2 = V2_S2 = R_V2 = 
    catchupvax_S1 = catchupvax_S2 = catchupvax_R = routinevax_S1 = routinevax_S2 = routinevax_R = 
    Ncatchupvax_S1 = Ncatchupvax_S2 = Ncatchupvax_R = Ncatchupvax_I2 = Ncatchupvax_C = 
    Nroutinevax_S1 = Nroutinevax_S2 = Nroutinevax_R = Nroutinevax_I2 = Nroutinevax_C = 
    NVaxrout = NVaxcamp = matrix(0, nrow=tf, ncol=al)
  
   ### initial number of people
  S1[1,] = round(sapply(sapply(c(.9*N0[1:7],.89*N0[8:al]), FUN=function(x) x-17), FUN=function(x) x*dt))
  I1[1,] = sapply(rep(17,al), FUN=function(x) x*dt)
  R[1,]  = rep(0,al)*dt
  S2[1,] = round(sapply(sapply(.1*N0, FUN=function(x) x-3), FUN=function(x) x*dt))
  I2[1,] = ceiling(sapply(rep(3,al), FUN=function(x) x*dt))
  C[1,]  = round(sapply(c(rep(0,7),round(.01*N0[8:al])), FUN=function(x) x*dt))
  
  V1[1,] = V2[1,] = rep(0,al)*dt
  
  W[1] = 0
  
  N[1,] = S1[1,] + I1[1,] + R[1,] + S2[1,] + I2[1,] + C[1,]
  
  D[1]     = rbinom(1, round(sum(I1[1,])), rep*.25) #adjust for underreporting
  D.adj[1] = D[1]*rep.adj
  # D[1] = rbinom(1, sum(I1[1,]), 0.0517*.25)
  
  ### Start loop
  
  t = 2
  
  theta2 <- rep(0,al)
  
  trans <- 1
  dur   <- 3
  
  newcases_agg <- 0
  outbreak     <- 0
  vax_t        <- 0
  catchup_camp <- 0
  D.wk         <- NULL
  D.wk.adj     <- NULL
  D.mth        <- NULL
  outbreak_mth <- NULL
  pop.wk       <- NULL
  pop.mth      <- NULL
  IR.mth       <- NULL
  SD.mth       <- NULL
  outbreak_t   <- 9999999999
  
  while ((t<=tf)){
    # while (t<=38240){
    print(c(t))
    
    if (t>=(t1)) {
      if (t<(t2)) {
        trans <- 1+(rt-1)*(t-(t1))/((t2-t1))
        delta  <- 1/(dur*trans)
      } else if (t>=(t2)) {
        trans <- rt
        delta  <- 1/(dur*rt)
      }
    }
    
    ### Force of infection
    lambdap <- matrix(sapply(betap, FUN=function(x) x*(1+q*cos(2*pi*(t*dt-lag)/(52.18)))),ncol=18) %*% (I1[(t-1),]+I2[(t-1),]+rC*C[(t-1),])/sum(N[(t-1),])*trans
    lambdaw <- 0
    
    
    ### Transitions
    births[t,]    <- rpois(al, B[(t-1),]*sum(N[(t-1),])*dt)
    
    S1_I1[t,]     <- rpois(al,((lambdap)* S1[(t-1),])*dt)
    
    deaths_S1[t,] <- sapply(mu[(t-1)]*S1[(t-1),]*dt, FUN = function(x) ifelse(x<0, -rpois(1,-x), rpois(1,x)))
    
    I1_R[t,]      <- rpois(al,delta*(1-alpha-theta)*I1[(t-1),]*dt)
    
    I1_C[t,]      <- rpois(al,delta*(theta)*I1[(t-1),]*dt)
    
    deaths_I1[t,] <- sapply(((delta*alpha)+mu[(t-1)])*I1[(t-1),]*dt, FUN = function(x) ifelse(x<0, -rpois(1,-x), rpois(1,x)))
   
    R_S2[t,]      <- rpois(al,omega*R[(t-1),]*dt)
    
    deaths_R[t,]  <- sapply(mu[(t-1)]*R[(t-1),]*dt, FUN = function(x) ifelse(x<0, -rpois(1,-x), rpois(1,x)))
    
    S2_S1[t,]     <- rpois(al,epsilon*S2[(t-1),]*dt)
    
    S2_I2[t,]     <- rpois(al,((lambdap)* S2[(t-1),])*dt)
    
    deaths_S2[t,] <- sapply(mu[(t-1)]*S2[(t-1),]*dt, FUN = function(x) ifelse(x<0, -rpois(1,-x), rpois(1,x)))
   
    I2_R[t,]      <- rpois(al,delta*(1-theta2)*I2[(t-1),]*dt)
    
    I2_C[t,]      <- rpois(al,delta*(theta2)*I2[(t-1),]*dt)
    
    deaths_I2[t,] <- sapply(mu[(t-1)]*I2[(t-1),]*dt, FUN = function(x) ifelse(x<0, -rpois(1,-x), rpois(1,x)))

    deaths_C[t,]  <- sapply(mu[(t-1)]*C[(t-1),]*dt, FUN = function(x) ifelse(x<0, -rpois(1,-x), rpois(1,x)))
    
    # aggregate cases at the week level
    if(t %% (1/dt) == 0) {
      # rep <- sample(rep.samps,1)
      newcases_agg_wk     <- rbinom(n = 1, size = round(sum(I1[(t-1/dt):(t-1),])), prob = rep*.25)
      newcases_agg_wk.adj <- newcases_agg_wk*rep.adj
      D.wk                <- c(D.wk,newcases_agg_wk)
      D.wk.adj            <- c(D.wk.adj,newcases_agg_wk.adj)
      pop_agg_wk          <- sum(N[(t-1/dt):(t-1),])
      pop.wk              <- c(pop.wk,pop_agg_wk)
    } 
   
    # aggregate cases at the month level 
    if(t %% (4/dt) == 0) {
      newcases_agg_mth <- sum(D.wk[(length(D.wk)-3):(length(D.wk))])
      D.mth           <- c(D.mth,newcases_agg_mth)
      pop_agg_mth     <- sum(pop.wk[(length(pop.wk)-3):(length(pop.wk))])
      pop.mth         <- c(pop.mth,pop_agg_mth)
      new.IR.mth      <- (newcases_agg_mth/pop_agg_mth)*100000
      IR.mth          <- c(IR.mth,new.IR.mth)
      new.SD.mth      <- (mean(D.mth) + (threshSD*(sd(D.mth))))
      if ((t>(2400/dt)) ) {
        if ((new.IR.mth >= threshIR)) {
          outbreak_t    <- c(outbreak_t,t)
          outbreak      <- 1
        } else if ( newcases_agg_mth >= new.SD.mth) {
          outbreak_t    <- c(outbreak_t,t)
          outbreak      <- 1
        } 
      }
    } 
    
    #Vaccination transitions
    V1_S1[t,]         <- rpois(al,vwane*V1[(t-1),]*dt)
    
    V2_S2[t,]         <- rpois(al,vwane*V2[(t-1),]*dt)
    
    deaths_V1[t,]     <- sapply(mu[(t-1)]*V1[(t-1),]*dt, FUN = function(x) ifelse(x<0, -rpois(1,-x), rpois(1,x)))

    deaths_V2[t,]     <- sapply(mu[(t-1)]*V2[(t-1),]*dt, FUN = function(x) ifelse(x<0, -rpois(1,-x), rpois(1,x)))

    Ncatchupvax_S1[t,] <- rpois(al, (massvacc[t,])*S1[(t-1),]*dt)
    Ncatchupvax_S2[t,] <- rpois(al, (massvacc[t,])*S2[(t-1),]*dt)
    Ncatchupvax_R[t,]  <- rpois(al, (massvacc[t,])*R[(t-1),]*dt)
    Ncatchupvax_I2[t,] <- rpois(al, (massvacc[t,])*I2[(t-1),]*dt) # "wasted" vaccine doses
    Ncatchupvax_C[t,]  <- rpois(al, (massvacc[t,])*C[(t-1),]*dt)  # "wasted" vaccine doses
    
    catchupvax_S1[t,] <- rpois(al, veff*Ncatchupvax_S1[t,])
    catchupvax_S2[t,] <- rpois(al, veff*Ncatchupvax_S2[t,])
    catchupvax_R[t,]  <- rpois(al, veff*Ncatchupvax_R[t,])
    
   
    #aging out
    age_out_S1[t,] <- rpois(al,u_out*S1[(t-1),]*dt)
    age_out_I1[t,] <- rpois(al,u_out*I1[(t-1),]*dt)
    age_out_R[t,]  <- rpois(al,u_out* R[(t-1),]*dt)
    age_out_S2[t,] <- rpois(al,u_out*S2[(t-1),]*dt)
    age_out_I2[t,] <- rpois(al,u_out*I2[(t-1),]*dt)
    age_out_C[t,]  <- rpois(al,u_out* C[(t-1),]*dt)
    age_out_V1[t,] <- rpois(al,u_out* V1[(t-1),]*dt)
    age_out_V2[t,] <- rpois(al,u_out* V2[(t-1),]*dt)

    #aging in 
    age_in_S1[t,] <- c(0, age_out_S1[t,1:(al-1)]) 
    age_in_I1[t,] <- c(0, age_out_I1[t,1:(al-1)]) 
    age_in_R[t,]  <- c(0, age_out_R[t,1:(al-1)]) 
    age_in_S2[t,] <- c(0, age_out_S2[t,1:(al-1)]) 
    age_in_I2[t,] <- c(0, age_out_I2[t,1:(al-1)]) 
    age_in_C[t,]  <- c(0, age_out_C[t,1:(al-1)]) 
    age_in_V1[t,] <- c(0, age_out_V1[t,1:(al-1)]) 
    age_in_V2[t,] <- c(0, age_out_V2[t,1:(al-1)]) 
    
    #people vaccinated (children 9mo)
    Nroutinevax_S1[t,] <- rpois(al, v[t,]*age_in_S1[t,])
    Nroutinevax_S2[t,] <- rpois(al, v[t,]*age_in_S2[t,])
    Nroutinevax_R[t,]  <- rpois(al, v[t,]*age_in_R[t,] )
    Nroutinevax_I2[t,] <- rpois(al, v[t,]*age_in_I2[t,]) # these are "wasted" vaccine doses
    Nroutinevax_C[t,]  <- rpois(al, v[t,]*age_in_C[t,] ) # these are "wasted" vaccine doses
    
    #people who successfully mount vaccine response
    routinevax_S1[t,] <- rpois(al, veff*Nroutinevax_S1[t,])
    routinevax_S2[t,] <- rpois(al, veff*Nroutinevax_S2[t,])
    routinevax_R[t,]  <- rpois(al, veff*Nroutinevax_R[t,])
    
    W_enter[t]       <- 0
    W_leave[t]       <- 0
    
    ### Transmission dynamics
    S1[t,] <- sapply(S1[(t-1),] - deaths_S1[t,] - age_out_S1[t,] + age_in_S1[t,] + births[t,]  - routinevax_S1[t,] - catchupvax_S1[t,] - S1_I1[t,]  + S2_S1[t,] + V1_S1[t,] , FUN = function(x) max(x,0))
    I1[t,] <- sapply(I1[(t-1),] - deaths_I1[t,] - age_out_I1[t,] + age_in_I1[t,] + S1_I1[t,] - I1_R[t,]  - I1_C[t,] , FUN = function(x) max(x,0))
    R[t,]  <- sapply(R[(t-1),]  - deaths_R[t,]  - age_out_R[t,]  + age_in_R[t,]  - routinevax_R[t,] - catchupvax_R[t,] + I1_R[t,] + I2_R[t,] - R_S2[t,] , FUN = function(x) max(x,0))
    S2[t,] <- sapply(S2[(t-1),] - deaths_S2[t,] - age_out_S2[t,] + age_in_S2[t,] - routinevax_S2[t,] - catchupvax_S2[t,] - S2_I2[t,] + R_S2[t,] - S2_S1[t,] + V2_S2[t,] , FUN = function(x) max(x,0))
    I2[t,] <- sapply(I2[(t-1),] - deaths_I2[t,] - age_out_I2[t,] + age_in_I2[t,] + S2_I2[t,] - I2_C[t,] - I2_R[t,] , FUN = function(x) max(x,0))
    C[t,]  <- sapply(C[(t-1),]  - deaths_C[t,]  - age_out_C[t,]  + age_in_C[t,]  + I2_C[t,] + I1_C[t,] , FUN = function(x) max(x,0))
    V1[t,] <- sapply(V1[(t-1),] - deaths_V1[t,] - age_out_V1[t,] + age_in_V1[t,] + routinevax_S1[t,] + catchupvax_S1[t,] - V1_S1[t,] , FUN = function(x) max(x,0))
    V2[t,] <- sapply(V2[(t-1),] - deaths_V2[t,] - age_out_V2[t,] + age_in_V2[t,] + routinevax_S2[t,] + routinevax_R[t,] + catchupvax_S2[t,] + catchupvax_R[t,] - V2_S2[t,] , FUN = function(x) max(x,0))
    
    W[t]   <- 0
    
    N[t,]  <- (S1+I1+R+S2+I2+C+V1+V2)[t,]
    
    NVaxrout[t,] <- Nroutinevax_S1[t,] + Nroutinevax_S2[t,] + Nroutinevax_R[t,] + Nroutinevax_I2[t,] + Nroutinevax_C[t,] 
    NVaxcamp[t,] <- Ncatchupvax_S1[t,] + Ncatchupvax_S2[t,] + Ncatchupvax_R[t,] + Ncatchupvax_I2[t,] + Ncatchupvax_C[t,]
    
    t = t+1
    
  }
  
  
  #aggregate to single time unit
  NewVaxrout  = rep(0,tmax)
  NewVaxcamp  = rep(0,tmax)
  NewPop  = rep(0,tmax)
  for (i in 1:tmax){
    # NewCase[i] = sum(I1[((i-1)*(1/dt)+1):(i/dt),])
    NewVaxrout[i]  = sum(NVaxrout[((i-1)*(1/dt)+1):(i/dt),])
    NewVaxcamp[i]  = sum(NVaxcamp[((i-1)*(1/dt)+1):(i/dt),])
    NewPop[i] = median(N[((i-1)*(1/dt)+1):(i/dt),])
  }

  FirstInf     = matrix(0, ncol=al, nrow=tmax)
  FirstInfNorm = matrix(0, ncol=al, nrow=tmax)
  for (i in 1:tmax){
    for (j in 1:al){
      FirstInf[i,j] = sum(I1[((i-1)*(1/dt)+1):(i/dt),j])
    }
    FirstInfNorm[i,] = FirstInf[i,]/sum(FirstInf[i,])
  }
  
  MeanAge       <- FirstInfNorm %*% c(.375,2.875,seq(7.5,82.5,5))
  BaselineCases <- FirstInf[2401,]

  # vax_start      <- min(which(v1[,2]>0))
  outbreak_start <- min(outbreak_t)

  return(list(D.wk,NewPop,NewVaxrout,NewVaxcamp,MeanAge,BaselineCases, 
              outbreak_start, D.mth, IR.mth,D.wk.adj,
              veff,vwane,vcov_C,rep.adj))
}


