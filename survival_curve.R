#Script two
#Survival curve prediction

library(survival)
library(survminer)
kmsurvo<-Surv(combined$Rectime,combined$Censor.Status)
sfit <- survfit(Surv(combined$Rectime,combined$Censor.Status)~SPIB, data=combined)
ggsurvplot(sfit, pval=TRUE, 
           legend.labs=c("Normal", "Altered"), legend=c(.70,.70),  
           title="SPIB")
ggsurvplot(sfit, legend = c(0.2, 0.2))
plot(sfit)

ggsurv<-ggsurvplot(sfit,  pval=TRUE, 
           legend.labs=c("Normal", "Altered"), 
           font.x = c(16, "bold.italic", "red"),
           font.y = c(16, "bold.italic", "darkred"),
           legend=c(.70,.70),  
           title="SPIB", main = "Survival curve",
           font.main = c(16, "bold", "darkblue"),
           font.tickslab = c(14, "plain", "darkgreen"))

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, color = "black", face = "bold"))
ggsurv
