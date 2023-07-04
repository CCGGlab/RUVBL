plot_CPH_forest<- function(clin_merged,surv_object, vars){
      # Df
      HR_multi_df<- data.frame(
            cond=vars,
            HR=NA,
            ci_l=NA,
            ci_h=NA,
            p=NA,
            n1=NA,
            n=NA
      )
      # CoxPH
      fit.coxph <- coxph(as.formula(paste0("surv_object ~ ", paste0(vars, collapse = "+"))), data = clin_merged)
      for(var in vars){
            HR_multi_df[HR_multi_df$cond==var,c("HR","ci_l","ci_h")]<- summary(fit.coxph)$conf.int[paste0(var), c(1,3,4)]
            HR_multi_df[HR_multi_df$cond==var,"p"]<- summary(fit.coxph)$coefficients[paste0(var),5]
            HR_multi_df[HR_multi_df$cond==var,c("n1","n")]<- c(summary(fit.coxph)$nevent, summary(fit.coxph)$n)
      }
      # Plot
      HR_multi_df$cond<- factor(HR_multi_df$cond, levels=rev(unique(HR_multi_df$cond)))
      fp_multi <- ggplot(data=HR_multi_df, aes(x=cond, y=HR, ymin=ci_l, ymax=ci_h, label=paste("P =",signif(p,2)))) +
            geom_pointrange(size=.3, fatten = 5) + 
            geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
            geom_text(aes(y=HR, label=signif(HR,2)), vjust=-1, size=7/3) +
            geom_text(aes(y=6.5, label=paste0("P=",signif(p,2))), fontface="italic", hjust=0, size=6/3) +
            coord_flip() +  # flip coordinates (puts labels on y axis)
            xlab("") +
            ylab("Hazard Ratio (Mean +/- 95% CI)") +
            scale_y_continuous(limits = c(0.15,10), trans='log10') +
            theme_bw() +  # use a white background
            theme(
                  axis.text = element_text(size=8),
                  axis.title = element_text(size=8)
            )
      return(fp_multi)
}

