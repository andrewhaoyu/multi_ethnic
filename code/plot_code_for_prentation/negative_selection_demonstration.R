MAF = seq(0.01,0.5,by=0.01)
N = length(MAF)
h2 = 0.4

beta_sn = h2/(N*sqrt(2*MAF*(1-MAF)))

beta_mn = (2*MAF*(1-MAF))^-0.125
u_beta_mn = beta_mn*sqrt(2*MAF*(1-MAF))
scale_factor = h2/sum(u_beta_mn^2)
beta_mn = scale_factor*beta_mn

beta_nn = (2*MAF*(1-MAF))^0
u_beta_nn = beta_nn*sqrt(2*MAF*(1-MAF))
scale_factor = h2/sum(u_beta_nn^2)
beta_nn = scale_factor*beta_nn

beta = c(beta_sn,beta_mn,beta_nn)

Negative_selection = factor(c(rep("Strong negative selection",N),
                       rep("Mild negative selection",N),
                       rep("No negative selection",N)),
                       levels = c("Strong negative selection",
                                  "Mild negative selection",
                                  "No negative selection"))
plot.data = data.frame(MAF=rep(MAF,3),beta,Negative_selection)



p = ggplot(plot.data,aes(MAF,beta,col=Negative_selection))+
  geom_point()+
  geom_line()+theme_Publication()+
  scale_colour_Publication()+
  #theme(axis.text.y=element_blank())+
  ylab("Beta")+
  xlab("MAF")+
  ggtitle("Negative selection demonstration")+
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(face = "bold",size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(face = "bold",size = rel(1.2)))+
  labs(col='Negative selection') 
png(file = paste0("./result/presentation_plot/Negative_selection.png"),width = 8, height = 6, units = "in",res = 300)
print(p)
dev.off()

