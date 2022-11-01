library(tidyverse)
library(marked)
library(patchwork)

capture_history <- readRDS("./survival_data.rds")

## 
redstart.proc=process.data(capture_history[,c(-1)], begin.time=2010)

# Create design data 
design.Phi=list(static=c("sex","HABITAT", "delay"))

design.parameters=list(Phi=design.Phi)

ddl=make.design.data(redstart.proc,parameters=design.parameters)


fit.models=function()
   {
   Phi.2=list(formula=~ delay + HABITAT + sex )
   p.2=list(formula=~ HABITAT)
   cml=create.model.list(c("Phi","p"))
   results=crm.wrapper(cml,data=redstart.proc, ddl=ddl,
                       external=FALSE,accumulate=FALSE, hessian=T)
   return(results)
}


model <- fit.models()

model[[1]]


p_estimates <- model[[1]]$results$reals$Phi

## Plotting Model Results

## plotting detection probability by year

ggplot(data = p_estimates %>% 
          group_by(delay, sex, HABITAT) %>% 
          summarize(estimate = mean(estimate), lcl = mean(lcl), ucl=mean(ucl)) %>%
          filter(), 
       aes(x=delay, y=estimate, ymin=lcl,ymax=ucl, group=sex)) + 
   geom_ribbon(aes(x=delay, ymin = lcl, ymax = ucl), alpha=0.4) +
   geom_point(aes(color=factor(sex)),size=3, position=position_dodge(width=0.3)) + 
   geom_smooth(aes(color=factor(sex)),method="lm", se=F, linetype=2) + scale_colour_manual(values = c("goldenrod1","black")) +
   theme_classic(base_size = 20) + xlab("Relative Migratory Timing (days)") + ylab("Apparent Annual Survival") +
   theme(legend.position = "top", legend.title = element_blank()) + 
   ylim(0,1) + facet_wrap(~HABITAT) #+  geom_errorbar(aes(color=factor(sex)),width=0.2, position=position_dodge(width=0.3))  

   