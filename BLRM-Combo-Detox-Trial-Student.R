# This code implements the two-dimensional BLRM design
# for the dual-agent combination trial with one compound (methadone) being defined externally
# and having a
# t Student's distribution
# as proposed by P Mozgunov et.al (2020) in
# ``A dose-finding design for dual-agent trials with patient-specific
# doses for one agent with application to an opiate detoxification trial''


########    Required packages
library("rjags")


#######     Specifying the MCMC Model for JAGS 
model1.string <-"
model {

for (i in 1:N){
logit(p[i]) <- alpha0[1] + alpha1[1] * log(b[i]/B.ref)
odds1[i]<- p[i]/(1-p[i])
logit(q[i])<- alpha0[2] + alpha1[2] * log(m[i]/M.ref)
odds2[i]<-q[i]/(1-q[i])
odds0[i]<-odds1[i] + odds2[i] + odds1[i] * odds2[i]
odds[i]<-odds0[i]*exp(eta * (b[i]/B.ref) * (m[i]/M.ref) ) 
tox[i]<-odds[i]/(1+odds[i]) # transforming odds to probabilities 
s[i] ~ dbern(tox[i])		 # Binomial model
}

eta~dnorm(0,priorEta)             # Normal prior for the interaction parameter eta
for(t in 1:2){
theta[1:2, t] ~ dmnorm(priorMean[1:2, t], # setting up bivariate normal prior for each single-agent model
priorPrec[1:2,
1:2,
t])
## extract actual coefficients
alpha0[t] <- theta[1, t]         # Tranforming variables so that alpha1 is always positive
alpha1[t] <- exp(theta[2, t])
}

}
"
model1.spec<-textConnection(model1.string)


### Scenario Function
true.tox<-function(a01,a11,a02,a12,e,b,m,B.ref=60,M.ref=60){
  term1 <- a01 + a11 * log(b/B.ref)
  p1<-exp(term1)/(1+exp(term1))
  term2 <- a02 + a12 * log(m/M.ref)
  p2<-exp(term2)/(1+exp(term2))
  odds1<-p1/(1-p1)
  odds2<-p2/(1-p2)
  odds0<-odds1 + odds2 + odds1*odds2
  odds<-odds0 * exp(e * (b/B.ref) * (m/M.ref)) 
  distr<-odds/(1+odds) 
  return(distr)
}

# End of Scenario Function

iter<-5000                              # Number of Samples Used by MCMC (increase to 40000 in individual trials)

testing.set<-112

baclofen<-c(10,30,60,90)
B.ref<-60
M.ref<-60

# Sample Size
cohort<-3                                 # Cohort Size             
ncohort<-16                              # Number of Cohorts
npatients<-cohort * ncohort               # Number of Patients


# Target Intervals
low.target<-  0.15                        # Lower bound of the target interval
upp.target<-  0.25                        # Upper bound of the target interval
target<-      0.20                        # Target probability

# Overdose Control
combo.safety<-0.25                        # Probability controlling overdosing
safety<-0.25                          # Probability controlling stopping for safety (for the first combination only)
futility<-0.925


safety.stopping<-T
futility.stopping<-T

mean.methadone<-51.4
sd.methadone<-23.3
min.methadone<-10
max.methadone.safety<-60
min.baclofen.safety<-30
min.methadone.test<-1
max.methadone.futility<-120
max.baclofen.futility<-90


# Degrees of Freedom for the Predictive Sample
DF<-10
####

##### Prior Belief Scenario
A01.true<-c(-4.20)
A02.true<-c(-5.45)
A11.true<-c(0.90)
A12.true<-c(0.05)
E.true<-c(0.85)


##### Low Toxicity 1 Scenario
# A01.true<-c(-4.20)
# A02.true<-c(-5.3)
# A11.true<-c(1.00)
# A12.true<-c(0.20)
# E.true<-c(1.00)


##### Low Toxicity 2 Toxicity Scenario
# A01.true<-c(-3.4)
# A02.true<-c(-5.6)
# A11.true<-c(0.50)
# A12.true<-c(0.50)
# E.true<-c(1.00)


##### Medium Toxicity 1 Toxicity Scenario
# A01.true<-c(-2.5)
# A02.true<-c(-3.6)
# A11.true<-c(1.5)
# A12.true<-c(1.0)
# E.true<-c(0.40)

##### Medium Toxicity 2 Toxicity Scenario
# A01.true<-c(-2.5)
# A02.true<-c(-2.0)
# A11.true<-c(2.5)
# A12.true<-c(1.5)
# E.true<-c(0.10)

##### High Toxicity 1 Toxicity Scenario 
# A01.true<-c(-2.5)
# A02.true<-c(-2.5)
# A11.true<-c(1)
# A12.true<-c(0.2)
# E.true<-c(0.5)


##### High Toxicity 2 Toxicity Scenario
# A01.true<-c(-2.0)
# A02.true<-c(-2.0)
# A11.true<-c(1.0)
# A12.true<-c(1.5)
# E.true<-c(0.40)


# ##### Unsafe Toxic Scenario
# A01.true<-c(-1.35)
# A02.true<-c(-1.35)
# A11.true<-c(0.90)
# A12.true<-c(0.0)
# E.true<-c(0.85)


#######     Specifying Prior Distribution

# This is calibrated operational prior
priorMean<-mat.or.vec(2,2)
priorMean[,1]<-c(-3.75,0.40)  # Prior Means for alpha0 and alpha1 for the first agent
priorMean[,2]<-c(-3.25,0.05)  # Prior Means for alpha0 and alpha1 for the second agent

priorPrec<-priorVar<-array(0,dim=c(2,2,2))
# priorVar[,,1]<-matrix(c(1.60^2,0.0,0.0,0.15^2),2,2)   # Prior Covariance Matrix for the first agent
# priorVar[,,2]<-matrix(c(1.60^2,0.0,0.0,0.15^2),2,2)   # Prior Covariance Matrix for the second agent
priorVar[,,1]<-matrix(c(0.50^2,0.0,0.0,0.35^2),2,2)   # Prior Covariance Matrix for the first agent
priorVar[,,2]<-matrix(c(0.50^2,0.0,0.0,0.35^2),2,2)   # Prior Covariance Matrix for the second agent
priorPrec[,,1]<-solve(priorVar[,,1])                  # Transforming to Precision for MCMC Model
priorPrec[,,2]<-solve(priorVar[,,2])

priorEta<-1.25                                         # Prior Precision (!) for the interaction term

#######     End of Specifying Prior Distribution

nsims<-2000

all.testing.accept<-mat.or.vec(nsims,1)
all.testing.best<-mat.or.vec(nsims,1)
all.testing.non.tox<-mat.or.vec(nsims,1)
all.ss<-all.dlts<-mat.or.vec(nsims,1)
all.b<-all.s<-all.m<-mat.or.vec(nsims,npatients)

counter<-counter.futility<-0

  for (z in 1:nsims){
    stop<-stop.safety<-stop.futility<-0
    
    b<-c(10,10,10)
    m<-round(rnorm(cohort,mean=mean.methadone,sd=sd.methadone),0)
    m[which(m<=min.methadone)]<-min.methadone
    
    s<-c()
    for (t in 1:3){
      s[t]<-rbinom(1,1,prob=true.tox(a01=A01.true,a02=A02.true,a11=A11.true,a12=A12.true,e=E.true,b=b[t],m=m[t]))
    }
    
    i<-1
    
    while(i < ncohort){
      N<-length(b)
      # Running MCMC model
      model1.spec<-textConnection(model1.string) # Uploading Model
      # Uploading the data
      mydata <- list(N=N,b=b,s=s,m=m,B.ref=B.ref,M.ref=M.ref,priorMean=priorMean,priorPrec=priorPrec,priorEta=priorEta)
      # Initialising the model
      jags <- jags.model(model1.spec,data =mydata,n.chains=1,n.adapt=iter,quiet=TRUE)
      # Updating the model
      update(jags, iter,progress.bar="none")
      # Storing the Samples
      tt<-jags.samples(jags,c('alpha0','alpha1','eta'),iter,progress.bar="none")
      # Enf of Running of MCMC Model
      
      
      # Storing Sample from the MCMC output as separate variables
      a01<-tt$alpha0[1,,] # Samples of the intercept parameter for the first agent
      a02<-tt$alpha0[2,,] # Samples of the intercept parameter for the second agent
      a11<-tt$alpha1[1,,] # Samples of the slope parameter for the first agent
      a12<-tt$alpha1[2,,] # Samples of the slope parameter for the second agent
      e<-tt$eta[1,,]      # Samples of the interation parameters
      
      m.new<-round(rnorm(cohort,mean=mean.methadone,sd=sd.methadone),0)
      m.new[which(m.new<=min.methadone)]<-min.methadone
      
      b.new<-s.new<-c()
      
      prob.mtd<-prob.above<-c()
      for (t in 1:cohort){
        for (j in 1:length(baclofen)){
          term1 <- a01 + a11 * log(baclofen[j]/B.ref)
          p1<-exp(term1)/(1+exp(term1))
          term2 <- a02 + a12 * log(m.new[t]/M.ref)
          p2<-exp(term2)/(1+exp(term2))
          odds1<-p1/(1-p1)
          odds2<-p2/(1-p2)
          # Distribution of odds under no interaction for the combination (j,k)
          odds0<-odds1 + odds2 + odds1*odds2
          # Distribution of odds for the combination (j,k)
          odds<-odds0 * exp(e * (baclofen[j]/B.ref) * (m.new[t]/M.ref))
          # Distribution of the toxicity probability for the combination (j,k)
          distr<-odds/(1+odds)
          # Probability that the toxicity probability for combo (j,k) is in the target interval
          prob.mtd[j]<-mean(distr<upp.target)-mean(distr<low.target)
          # Probability that the toxicity probability for combo (j,k) is above the target interval
          prob.above[j]<-mean(distr>(upp.target))
        }
        
        prob.mtd[which(prob.above>combo.safety)]<-0
        
        
        # Do not skip doses restriction
        max.curr.dose<-which(baclofen==max(b))
        if(!(max.curr.dose==length(baclofen) | max.curr.dose==(length(baclofen)-1) ) ){
          prob.mtd[(max.curr.dose+2):length(baclofen)]<-0
        }

        
        # Coherence Restriction
        last.dlts<-s[(length(s)-cohort+1):length(s)]
        last.doses<-b[(length(b)-cohort+1):length(b)]
        
        if(sum(last.dlts)!=0){
          if(which(baclofen==min(last.doses[which(last.dlts==1)]))!=length(baclofen)){
            non.coherent<-(which(baclofen==min(last.doses[which(last.dlts==1)]))+1):length(baclofen)
            prob.mtd[non.coherent]<-0
          }
        }
        
  
        
        b.new[t]<-baclofen[which.max(prob.mtd)]
        s.new[t]<-rbinom(1,1,prob=true.tox(a01=A01.true,a02=A02.true,a11=A11.true,a12=A12.true,e=E.true,b=b.new[t],m=m.new[t]))
      }
      
      # Safety Stopping
      if(safety.stopping==T){
        
        term1 <- a01 + a11 * log(min.baclofen.safety/B.ref)
        p1<-exp(term1)/(1+exp(term1))
        term2 <- a02 + a12 * log(max.methadone.safety/M.ref)
        p2<-exp(term2)/(1+exp(term2))
        odds1<-p1/(1-p1)
        odds2<-p2/(1-p2)
        odds0<-odds1 + odds2 + odds1*odds2
        odds<-odds0 * exp(e * (min.baclofen.safety/B.ref) * (max.methadone.safety/M.ref))
        distr<-odds/(1+odds)
        prob.above.safety<-mean(distr>(upp.target),na.rm=T)
        
        

        if(prob.above.safety>safety){
          stop<-stop.safety<-1
          break()
        }
        
      }  

      
      
      # Futility Stopping
      if(futility.stopping==T){
        
        term1 <- a01 + a11 * log(max.baclofen.futility/B.ref)
        p1<-exp(term1)/(1+exp(term1))
        term2 <- a02 + a12 * log(max.methadone.futility/M.ref)
        p2<-exp(term2)/(1+exp(term2))
        odds1<-p1/(1-p1)
        odds2<-p2/(1-p2)
        odds0<-odds1 + odds2 + odds1*odds2
        odds<-odds0 * exp(e * (max.baclofen.futility/B.ref) * (max.methadone.futility/M.ref))
        distr<-odds/(1+odds)
        prob.below.futility<-mean(distr<(low.target),na.rm=T)
        
        
        if(prob.below.futility>futility & max(b)==max(baclofen)){
          counter.futility<-counter.futility+1
          stop.futility<-1
          break()
        }
        
      }  

      
      
      b<-c(b,b.new)
      m<-c(m,m.new)
      s<-c(s,s.new)
      i<-i+1
    }
    
    
    
    all.ss[z]<-length(m)
    all.dlts[z]<-sum(s)
    all.m[z,1:length(m)]<-m
    all.s[z,1:length(s)]<-s
    all.b[z,1:length(b)]<-b
    
    
    if(stop==0){
      
      
      N<-length(b)
      model1.spec<-textConnection(model1.string) # Uploading Model
      mydata <- list(N=N,b=b,s=s,m=m,B.ref=B.ref,M.ref=M.ref,priorMean=priorMean,priorPrec=priorPrec,priorEta=priorEta)
      jags <- jags.model(model1.spec,data =mydata,n.chains=1,n.adapt=iter,quiet=TRUE)
      update(jags, iter,progress.bar="none")
      tt<-jags.samples(jags,c('alpha0','alpha1','eta'),iter,progress.bar="none")
      a01<-tt$alpha0[1,,] # Samples of the intercept parameter for the first agent
      a02<-tt$alpha0[2,,] # Samples of the intercept parameter for the second agent
      a11<-tt$alpha1[1,,] # Samples of the slope parameter for the first agent
      a12<-tt$alpha1[2,,] # Samples of the slope parameter for the second agent
      e<-tt$eta[1,,]      # Samples of the interation parameters
      
      m.testing<-round(sqrt(sd.methadone)*qt(runif(testing.set,0,1),df=DF)+mean.methadone,0)
      m.testing[which(m.testing<=min.methadone.test)]<-min.methadone.test
      b.testing.true<-b.testing<-c()
      prob.testing.true<-mat.or.vec(testing.set,1)
      
      for (u in 1:length(m.testing)){
        b.testing.true[u]<-  baclofen[which.min(abs(true.tox(a01=A01.true,a02=A02.true,a11=A11.true,a12=A12.true,e=E.true,b=baclofen,m=m.testing[u])-target))]
      }
      
      prob.mtd<-prob.above<-c()
      for (t in 1:length(m.testing)){
        for (j in 1:length(baclofen)){
          term1 <- a01 + a11 * log(baclofen[j]/B.ref)
          p1<-exp(term1)/(1+exp(term1))
          term2 <- a02 + a12 * log(m.testing[t]/M.ref)
          p2<-exp(term2)/(1+exp(term2))
          odds1<-p1/(1-p1)
          odds2<-p2/(1-p2)
          # Distribution of odds under no interaction for the combination (j,k)
          odds0<-odds1 + odds2 + odds1*odds2
          # Distribution of odds for the combination (j,k)
          odds<-odds0 * exp(e * (baclofen[j]/B.ref) * (m.testing[t]/M.ref))
          # Distribution of the toxicity probability for the combination (j,k)
          distr<-odds/(1+odds)
          # Probability that the toxicity probability for combo (j,k) is in the target interval
          prob.mtd[j]<-mean(distr<upp.target)-mean(distr<low.target)
          # Probability that the toxicity probability for combo (j,k) is above the target interval
          prob.above[j]<-mean(distr>(upp.target))
        }
        
        prob.mtd[which(!prob.above<=combo.safety)]<-0
        b.testing[t]<-baclofen[which.max(prob.mtd)]
      }
      
      
      
      true.testing.at.estimated<-c()
      for (u in 1:length(m.testing)){
        true.testing.at.estimated[u]<-true.tox(a01=A01.true,a02=A02.true,a11=A11.true,a12=A12.true,e=E.true,b=b.testing[u],m=m.testing[u])
        true.testing<-true.tox(a01=A01.true,a02=A02.true,a11=A11.true,a12=A12.true,e=E.true,b=baclofen,m=m.testing[u])
        if(any(true.testing>=low.target & true.testing <=upp.target)){
          true.set<-baclofen[which(true.testing>=low.target & true.testing <=upp.target)]
          if(b.testing[u] %in% true.set){
            prob.testing.true[u]<-1
          } 
        }else{
          if(b.testing[u]==baclofen[which.min(abs(true.testing-target))]){
            prob.testing.true[u]<-1
          }
        }
      }
      
      all.testing.accept[z]<-mean(prob.testing.true)
      all.testing.non.tox[z]<-length(which(true.testing.at.estimated<=upp.target ))/testing.set
      all.testing.best[z]<-sum(b.testing.true==b.testing)/testing.set
    }else{
      counter<-counter+1
      cat("stopped","\n")
      all.testing.accept[z]<-0
      all.testing.best[z]<-0
    }
    cat(z, "out of",nsims,"\n")
  }
  
  # Proportion of Correct Selections in the Predictive Sample
  mean(all.testing.accept[all.testing.accept!=0])
  # Proprtion of Safe Selections in the Predictive Sample
  mean(all.testing.non.tox)
  # Proportion of DLTs
  mean(all.dlts)/mean(all.ss)
  # Average Sample Size
  mean(all.ss)
  #Proportion of Trials stopped for unsafe
  counter/nsims
  # Proportion of Trials stopped for `futility`
  counter.futility/nsims
