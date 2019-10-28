#goal: plot prediction results
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/')
library(ggplot2)
p.thr <- c(1E-8,5E-8,1E-7,5E-7,1E-6,5E-6,1E-5,5E-5,1E-4,5E-4,1E-3,5E-3,1E-2,0.1,0.3,0.5)
log10P <- -log10(p.thr)
log10P <- round(log10P,1)
n.cut <- length(p.thr)
prop.result <- matrix(0,3,n.cut)
n.snp.result <- matrix(0,3,n.cut)
vad.r2 <- rep(0,3)
pop <- c("European",
         "African",
         "Latino")
gr = 2
for(i in 1:3){
  load(paste0("./multi_ethnic/result/LDP_summary_",i,"_",gr)) 
  n.snp.result[i,] <- round(LDP.result[[1]],0)
  prop.result[i,] <- round(LDP.result[[2]],2)
  vad.r2[i] <- mean(LDP.result[[5]])
  r2.test <- LDP.result[[3]]
  data <- data.frame(p.thr,r2.test,prop.result[i,])
  colnames(data) <- c("Pthr","R2","Prop")
  
  
  png(paste0("./multi_ethnic/result/prediction_",pop[i],"_",gr,".png"),width = 24, height=12,unit = "cm",res = 300)
 print({
   ggplot(data)+geom_line(aes(-log10(Pthr),R2))+
     fte_theme()+
     xlab("-log10(P-value)")+
     ylab("Adjusted R2")+
     ggtitle(paste0("Prediction result for ",pop[i]," population"))+
     scale_x_discrete(limit= log10P)
 }) 
 dev.off()
 png(paste0("./multi_ethnic/result/proportion_",pop[i],"_",gr,".png"),width = 24, height=12,unit = "cm",res = 300)
 print({
   ggplot(data)+geom_line(aes(-log10(Pthr),Prop))+
     fte_theme()+
     xlab("-log10(P-value)")+
     ylab("Proportion of causal SNPs in the selected SNPs")+
     ggtitle(paste0("Prediction result for ",pop[i]," population"))+
     scale_x_discrete(limit= log10P)
 }) 
 dev.off() 
  
}

write.csv(n.snp.result,"./multi_ethnic/result/number_of_selected_SNP.csv")
write.csv(prop.result,"./multi_ethnic/result/causal_proportion.csv")

##plot the results for EB and original regression coef
load("./multi_ethnic/result/test_result_ori_EB.Rdata")

p.thr <- c(10^-8,5E-8,10^-7,5E-7,10^-6,5E-6,10^-5,5E-5,10^-4,5E-4,10^-3,5E-3,10^-2,0.1,0.3,0.5)

p.thr2 <- c(10^-8,5E-8,10^-7,5E-7,10^-6,5E-6,10^-5,5E-5,10^-4,5E-4,10^-3,5E-3,10^-2,0.1,0.3,0.5)


library(corrplot)
library(colortools)
ind <- 3
pop <- c("African","Latino")
for(i in 1:2){
for(j in 1:2){
  #column is the African
  #row is Europe
    # result <- matrix(data[,temp],length(p.thr),length(p.thr))
    # colnames(result) <- round(-log10(p.thr),1)
    # rownames(result) <- round(-log10(p.thr),1)
    # 
  temp <- data[,c(1:2,ind)] 
  colnames(temp)[3] <- "r2"
  temp = temp %>% 
    mutate(logp_train = as.factor(round(-log10(pthr_train),1))) %>% 
    mutate(logp_ref = as.factor(round(-log10(pthr_ref),1))) 
  png(paste0("./multi_ethnic/result/dynamic_",pop[i],"_method",j,".png"),width = 16, height=12,unit = "cm",res = 300)
      print({
        base_size <- 9
        ggplot(data = temp, aes(x=logp_train, y=logp_ref)) + 
           geom_tile(aes(fill=r2),colour = "white")+ 
          scale_fill_gradient2(low = "grey94",
                              high = "dodgerblue4") +
          fte_theme()+
          labs(x = paste0(pop[i]," -log10(P-value)"),
               y = paste0("European -log10(P-value)")) +
          ggtitle(paste0("Adjusted r2 in ",pop[i]," population"))+
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) + 
          labs(fill = "Adjusted r2")
          
      }
      ) 
      dev.off() 
      ind <- ind+1
      }
    
}

i <-2 
j <- 2
idx <- which.max(data[,5])
data[idx,]


