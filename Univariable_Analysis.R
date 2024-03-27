require("survival")
require("DescTools")

m3b=read.csv("Univariable_Analysis_Data.csv")

m6=NULL
for (j1 in 4:ncol(m3b)){
  print(j1)
  dat3=data.frame(time=m3b$time, status=m3b$status, VAR=factor(m3b[, j1]))
  t1=table(dat3$VAR==1, dat3$status)
  N1=t1[2, 1]+t1[2, 2]
  N2=t1[2, 2]
  P2=round(N2/N1*100, 1)
  P3=round(N2/(N1*mean(dat3$time/365))*100, 1)
  
  x1=round(as.numeric(BinomCI(N2, N1))*100,1)
  
  fit=survdiff(Surv(time, status) ~ VAR, data = dat3)
  P_Global=pchisq(fit$chisq, length(fit$n)-1, lower.tail = FALSE)
  
  m6=rbind(m6, 
           data.frame(VAR=colnames(m3b)[j1],
                      N0=N1,
                      N2=paste0(N2, "(", P2, "%)"),
                      N3=paste0( x1[2], "% to ", x1[3], "%"),
                      P=P_Global,
                      )
  )
}

m6$P2=p.adjust(p=m6$P, method = "BH")
#Adjust P-values for Multiple Comparisons 
