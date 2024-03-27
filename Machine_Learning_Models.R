require("ggplot2")
require("survival")
require("mlr3")
require("mlr3proba")
require("mlr3verse")
require("mlr3viz")
require("mlr3learners")
require("randomForestSRC")
require("SurvMetrics")

m1=read.csv("ML_data.csv")

##Prepare train and test data 
set.seed(2024) 
rows_train=sample(1:nrow(m1), round(nrow(m1)/2), replace = F)
rows_test=(1:nrow(m1))[-rows_train]
dat_train=m1[rows_train, ]
dat_test=m1[rows_test, ]

##Compare Models 
measure = msr("surv.cindex")

Lnum0=c("surv.rfsrc", "surv.coxph", "surv.ranger",
        "surv.aorsf", 
        "surv.cforest","surv.coxboost", "surv.ctree",
        "surv.cv_coxboost", "surv.gbm", 
        "surv.glmboost", "surv.glmnet", "surv.rpart",
        "surv.xgboost")

task2 = as_task_surv(dat_train, time="time", event = "status", 
                     id="cluster")
set.seed(100)
bmr = benchmark(benchmark_grid(task2, lrns(Lnum0),
                                rsmp("cv", folds = 10)))
bmr$aggregate(measure) 
# performance metrics (C-index) for all models 
# Best model: Random Survival Forest

#Hyper parameter optimization 
set.seed(100)
o.tune <- randomForestSRC::tune(Surv(time, status) ~ ., 
                                data = dat_train[, -1])
o.tune$optimal
#nodesize     mtry 
#20       20 

#Model train
set.seed(100)
RFSRC.train <- rfsrc(Surv(time, status) ~ ., data = dat_train[, -1],
                     nodesize=20, mtry=20, importance = "TRUE")

##Plot variable importance 
plot(RFSRC.train, cex=0.5)

##Predicted Risk score distribution
predicted_score=RFSRC.train$predicted
score_Q1Q4=quantile(predicted_score, 
                    probs = c(0.25,0.75))

# use 25 percentile and 75 percentile of the predicted Risk score
# as risk group cutoff

##Plot Risk score distribution
df1=dat_train[, c("time", "status")]
df1$score=predicted_score
df1$RISK_GROUP="MODERATE"
df1$RISK_GROUP[df1$score<=score_Q1Q4[1]]="LOW"
df1$RISK_GROUP[df1$score>score_Q1Q4[2]]="HIGH"

ggplot(df1, aes(x=score)) + 
  geom_histogram(binwidth=1,fill="white",color="black")+
  annotate("text", x=score_Q1Q4[1], y=60, 
           label=paste0("Q1:", round(score_Q1Q4[1],1)))+
  annotate("text", x=score_Q1Q4[2], y=60, 
           label=paste0("Q3:", round(score_Q1Q4[2],1)))+
  geom_vline(xintercept =score_Q1Q4, linetype="dotted", 
             color = "black", linewidth=0.5)+
  xlab("Predicted Risk Score")+
  ylab("Count")

# Predictive performance in validation data set

o.pred <- predict(object = RFSRC.train, dat_test[,-1])

predicted_score=o.pred$predicted

df1=dat_test[, c("time", "status")]
df1$score=predicted_score
df1$RISK_GROUP="MODERATE"
df1$RISK_GROUP[df1$score<=score_Q1Q4[1]]="LOW"
df1$RISK_GROUP[df1$score>score_Q1Q4[2]]="HIGH"

#Plot C-index

time <- df1$time
status <- df1$status
predicted <- df1$score
C_index=Cindex(Surv(time, status), 1-predicted)

C_index_bootstrap=NULL
set.seed(100)
for (j1 in 1:1000){
  resample=sample(nrow(df1), replace=T)
  time <- df1$time[resample]
  status <- df1$status[resample]
  predicted <- df1$score[resample]
  C_index_bootstrap=c(C_index_bootstrap, 
                      Cindex(Surv(time, status), 1-predicted))
  print(j1)
}

df2=data.frame(DAT="TEST", CINDEX=C_index, CINDEX_BS=C_index_bootstrap)
C_index_CI=quantile(df2$CINDEX_BS, 
                    probs = c(0.025,0.975))

ggplot(df2, aes(DAT, CINDEX_BS))+
  geom_boxplot()+
  xlab("Validation")+
  ylab("C-index")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())

#Risk Score Distribution

ggplot(df1, aes(x=score, fill=RISK_GROUP, color=RISK_GROUP)) + 
  geom_histogram(binwidth=1)+
  annotate("text", x=0, y=60, 
           label=paste0("<=", round(score_Q1Q4[1],1)), 
           hjust=0,color=COL3[2])+
  annotate("text", x=40, y=60, 
           label=paste0(">", round(score_Q1Q4[2],1)),
           color=COL3[1])+
  annotate("text", x=25, y=60, 
           label=paste0(round(score_Q1Q4[1],1), "-",
                        round(score_Q1Q4[2],1)),
           color=COL3[3])+
  geom_vline(xintercept =score_Q1Q4, linetype="dotted", 
             color = "black", linewidth=0.5)+
  xlab("Predicted Risk Score")+
  ylab("Count")

##Plot risk groups
df1a=df1
df1$RISK_GROUP=factor(df1$RISK_GROUP, 
                      levels=c("LOW", "MODERATE", "HIGH"))
df1$status=factor(df1$status, levels=c(T,F))

t1=table(df1$RISK_GROUP, df1$status)
TEXTLABEL=paste0(round(100*t1[,1]/(t1[,1]+t1[,2])), "%")
TEXTY=t1[,1]+t1[,2]

ggplot(df1, aes(RISK_GROUP))+
  geom_bar(aes(fill = status))+
  scale_fill_manual(name="Cluster\nFormation",values=COL2[c(1,2)])+
  annotate("text", x=1:3, y=TEXTY, label=TEXTLABEL, vjust=-0.2, 
           color=COL2[1], size=5)+
  ylab("Count")

##Plot survival curve by RISK_GROUP

df1=df1a

df1$RISK_GROUP=factor(df1$RISK_GROUP, 
                      levels=c("HIGH", "MODERATE", "LOW"))

fit <- survfit(Surv(time, status) ~ RISK_GROUP, data = df1)

SUR150=summary(fit, times = c(90))
YTEXT3=paste0(
  sprintf("%.2f", round(1-SUR150$surv,2)), "(", 
  sprintf("%.2f", round(1-SUR150$upper,2)), "-", 
  sprintf("%.2f", round(1-SUR150$lower,2)), ")")

F1A=
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = df1,             # data used to fit survival curves.
    risk.table = TRUE,       # show risk table.
    pval = F,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for 
    # point estimates of survival curves.
    xlim = c(0,750),        # present narrower X axis, but not affect
    # survival estimates.
    ylim=c(0,0.5),
    cumevents=T,            # show the cumulative number of events
    xlab = "Days after detection",   # customize X axis label.
    ylab = "Cumulative incidence",   # customize Y axis label.
    break.time.by = 150,     # break X axis in time intervals by 50.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    tables.y.text.col = T, # colour risk table text annotations.
    tables.y.text = F, # show bars instead of names in text annotations
    # in legend of risk table
    fun = "event",
    legend = "top",
    legend.title="",
    censor=F
  )

F1=F1A
F1$plot <- F1$plot+
  ggplot2::annotate("text", 
                    x = 0, y = c(0.3, 0.35, 0.4), colour=COL3[3:1],
                    label = rev(YTEXT3), size = 4, hjust=0, vjust=0)+
  ggplot2::annotate("text", 
                    x = 0, y = 0.45,
                    label = "90-day cumulative incidence (95% CI)", 
                    size = 4, hjust=0, vjust=0)
F1