require("ggplot2")
require("survminer")
require("survival")
require("adjustedCurves")
require("riskRegression")

m2a=read.csv("Multivariable_MetaData.csv")

# Candidate Variable from univariable analysis based on p<0.2

Candidate_Var=c("EMM", "emm.pattern", "TET_NS", "EXOTOXIN_A",
                "EXOTOXIN_C","EXOTOXIN_G","EXOTOXIN_J",
                "EXOTOXIN_S","SURFACE_FBAA","SURFACE_PRTF2",
                "SURFACE_R28", "SURFACE_SFB1","SURFACE_SOF",
                "CAPSULE", "ROCA", "SDA1", "SIC","AGE_GROUP6",
                "SYN6","RESI4","ALCOHOL","CIRR","PWID",
                "SEASON2", "REGION4")


## Correlation with emm
m1=m2a
dat=m1[, c(Candidate_Var, "FINAL_EMMTYPE")]
EMM=factor(m1$FINAL_EMMTYPE)
dat$FINAL_EMMTYPE=NULL
dat$EMM=NULL

m5=NULL
for (j1 in colnames(dat)){
  t1=(assocstats(table(factor(dat[, j1]), EMM))) 
  m5=rbind(m5, data.frame(VAR1="EMM", VAR2=j1, CramerV=t1$cramer))
}
m5[order(m5$CramerV),]
EMM_CORT=(m5$VAR2)[m5$CramerV<0.8]

dat2=m1[, EMM_CORT]
NC=dim(dat2)[2]-1
COL=colnames(dat2)
m6=NULL
for (j1 in 1:NC){
  for (j2 in (j1+1):(NC+1)){
    t1=(assocstats(table(factor(dat2[, j1]), factor(dat2[, j2]))))
    m6=rbind(m6, data.frame(VAR1=COL[j1], VAR2=COL[j2], CramerV=t1$cramer))
  }
}
m6[order(m6$CramerV),]
m6a=m6[m6$CramerV>0.8,]
ADDCORT=unique(c(m6a$VAR1, m6a$VAR2))

UNVAR=setdiff(EMM_CORT, ADDCORT) 
# Remaining variables after removing highly correlated variables

# data set for Cox regression

m2=read.csv("Multivariable_Event_MetaData.csv")
m1=cbind(m2[, c(1:3)],
         as.data.frame(lapply(m2[,4:ncol(m2)], as.factor)))

m1$AGE_GROUP6=relevel(m1$AGE_GROUP6, ref = "18-34")
m1$RACE3=relevel(m1$RACE3, ref = "White")
m1$SYN6=relevel(m1$SYN6, ref = "CELLSYN")
m1$RESI4=relevel(m1$RESI4, ref = "Private")
m1$SEASON2=relevel(m1$SEASON2, ref = "Q2Q3")
m1$REGION4=relevel(m1$REGION4, ref = "NE(CT, NY, MD)")

EMM_Levels=as.character(c(1, 89 ,12 ,28 ,4,11,77 ,3,6, 92, "OTHER"))
m1$EMM=factor(m1$EMM, levels=EMM_Levels)
m1$EMM=relevel(m1$EMM, ref = "89")
m1$emm.pattern=relevel(m1$emm.pattern, ref = "A-C")

cox_dat=m1[, c("time", "status", "EMM", UNVAR)] 

##Initial Full Model
fit = coxph(Surv(time, status)~., data=cox_dat)
summary(fit) 

##Model selection
selectedMod <- step(fit) #Choose a model by AIC in a Stepwise Algorithm
summary(selectedMod) #Selected model

## Test proportional hazards (PH) assumption 
test.ph = cox.zph(selectedMod)
ggcoxzph(test.ph) # Plot scaled Schoenfeld residuals


##Plot adjusted probability based on the fitted Cox model

fit2 = coxph(Surv(time, status) ~ EMM + EXOTOXIN_C + 
                         RESI4 + ALCOHOL + CIRR + PWID + 
                         SEASON2 + REGION4, data=cox_dat, x=TRUE)
#EMM
adjsurv <- adjustedsurv(data=cox_dat,
                        variable="EMM",
                        ev_time="time",
                        event="status",
                        method="direct",
                        outcome_model=fit2,
                        conf_int=TRUE)

plot(adjsurv, conf_int=TRUE)

#Residence
adjsurv <- adjustedsurv(data=cox_dat,
                        variable="RESI4",
                        ev_time="time",
                        event="status",
                        method="direct",
                        outcome_model=fit2,
                        conf_int=TRUE)

plot(adjsurv, conf_int=TRUE)

#Injection drug use
adjsurv <- adjustedsurv(data=cox_dat,
                        variable="PWID",
                        ev_time="time",
                        event="status",
                        method="direct",
                        outcome_model=fit2,
                        conf_int=TRUE)
plot(adjsurv, conf_int=TRUE)

#Season
adjsurv <- adjustedsurv(data=cox_dat,
                        variable="SEASON2",
                        ev_time="time",
                        event="status",
                        method="direct",
                        outcome_model=fit2,
                        conf_int=TRUE)

plot(adjsurv, conf_int=TRUE)
