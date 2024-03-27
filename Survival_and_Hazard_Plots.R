require("ggplot2")
require("survival")
require("survminer")
require("bshazard")

dat=read.csv("Event_Data.csv")

#Survival Curve
fit <- survfit(Surv(time, status) ~ 1, data = dat)
F1A=
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = dat,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,750),        # present narrower X axis, but not affect
                          # survival estimates.
  cumevents=T,            # show the cumulative number of events
  xlab = "Days after detection",   # customize X axis label.
  ylab = "Cumulative incidence",   # customize Y axis label.
  break.time.by = 150,     # break X axis in time intervals by 50.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  tables.y.text.col = T, # colour risk table text annotations.
  tables.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  fun = "event",
  legend = "none",
  legend.title="",
  censor=F
)

F1A

## Hazard Function
dat$status=as.integer(dat$status)

haz_est <- bshazard(Surv(time, status) ~ 1, data = dat)
plot(haz_est)

m5=data.frame(time=haz_est$time, 
              est = haz_est$hazard, 
              LCI=haz_est$lower.ci, UCI=haz_est$upper.ci)
m5a=m5[m5$time <=750, ]

F1B=
ggplot(m5a, aes(x = time, y = est)) +
  geom_line(color="red", linewidth=1.5, alpha=0.7)+
  xlab("Days after detection") + ylab("Hazard")+
  geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.2)+
  xlim(c(0,750))+
  scale_x_continuous(breaks = 0:5*150)

F1B  