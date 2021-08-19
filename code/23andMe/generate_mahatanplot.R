com.args = commandArgs(trailingOnly = T)
#i1 ethnic group
#i2 trait
#i3 =1 means mega+hapmap i3 = 2 means all SNPs
i1 = as.numeric(com.args[[1]])
i2 = as.numeric(com.args[[2]])
i3 = as.numeric(com.args[[3]])
eth <- c("EUR","AFR","AMR","EAS","SAS")
eth_name = c("European","African American",
             "Latino","East Asian",
             "South Asian")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")

trait_name = c("Any CVD","Depression",
               "Heart metabolic disease burden",
               "Height",
               "Sing back musical note",
               "Migraine Diagnosis",
               "Morning Person")
setwd("/data/zhangh24/multi_ethnic/data/cleaned/")
source("/data/zhangh24/multi_ethnic/code/stratch/theme_publication.R")
library(data.table)
library(qqman)
library(dplyr)
if(i3==1){
  data <- fread(paste0(eth[i1],"/sumdat/",trait[i2],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T)  
}else{
  data <- fread(paste0(eth[i1],"/sumdat/",trait[i2],"_passQC_noNA_matchinfo_matchMAFnoNA_common.txt"),header=T)  
  data = data %>% 
    mutate(CHR=gsub("chr","",scaffold)) %>% 
    rename(rsid = assay.name,
           BP = position,
           FREQ_A1 = freq.a,
           P = pvalue) %>% 
    mutate(CHR= ifelse(CHR=="X",23,CHR)) %>% 
    select(rsid,CHR,BP,FREQ_A1,P)
}


dat = data %>% 
  mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>% 
  select(rsid,CHR,BP,P,MAF) %>% 
  rename(SNP = rsid)

sample_size <- as.data.frame(fread("/data/zhangh24/multi_ethnic/data/23_sample_size.csv",header=T))
library(readr)
library(dplyr)
x = dat$P
z = qnorm(x / 2)
lambda = round(median(z^2) / qchisq(0.5,1), 3)

idx <- which(sample_size$eth==eth[i1]&
               sample_size$Disease==trait[i2])
N.effect  <- sample_size[idx,"N_effect"]
#rescale lambda to 1000 subjects
if(i2%in%c(1:2,5:7)){
  lambda_1000 = round(1+500*(lambda-1)/N.effect ,3)
}else{
  lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)
}


convert.qval.pval = function(qvalues) {
  # you need to know the estimate of pi0 used to create the q-value
  # that's the maximum q-value (or very, very close to it)
  pi0 = max(qvalues)
  # compute m0, the estimated number of true nulls
  m0 = length(qvalues) * pi0
  # then you multiply each q-value by the proportion of true nulls
  # expected to be under it (the inverse of how you get there from
  # the p-value):
  return(qvalues * rank(qvalues) / m0)
}
p.pwas <- 5E-08
#q.pwas <- convert.qval.pval(c(tmp, 0.05))[(nrow(dat)+1)]

nCHR <- length(unique(dat$CHR))
dat$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(dat$CHR)){
  nbp[i] <- max(dat[dat$CHR == i,]$BP)
  dat$BPcum[dat$CHR == i] <- dat$BP[dat$CHR == i] + s
  s <- s + nbp[i]
}
library(dplyr)
axis.set <- dat %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(dat$P)))) + 2 
sig1 <- p.pwas
#sig2 <- q.pwas

#sigline <- data.frame(sig=c(-log10(sig1),-log10(sig2)),val=c(paste0("P=",signif(sig1,2)),"FDR=0.05"))
sigline <- data.frame(sig=c(-log10(sig1)),val=c(paste0("P=",signif(sig1,2))))
library(ggplot2)
manhplot <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
                            color = as.factor(CHR), size = -log10(P))) +
  geom_point(alpha = 0.8, size=0.8) + 
  # annotate("text", x=1649636522, y=9, label= "ABO") +
  # annotate("text", x=1749636522, y=8, size=3, color="grey40",
  #          label= "C2")+
  # annotate("text", x=2026257501, y=10.5, label= "OAS1") +
  # annotate("text", x=2126257501, y=9.6, size=3, color="grey40",
  #          label= "A2") +
  # annotate("text", x=2126307501, y=8.2, size=3, color="grey40",
  #          label= "C2") +
  # annotate("text", x=2126307501, y=7.6, size=3, color="grey40",
  #          label= "B2") +
  # annotate("text", x=2655970388, y=6, label= "BCAT2") +
# annotate("text", x=2755970388, y=4.7, size=3, color="grey40",
#          label= "C2") +
scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(data = sigline, aes(yintercept = sig), color= "red", linetype="dashed") +
  guides(color = F) + 
  labs(x = NULL, 
       y = "-log10(p)", 
       linetype = "",
       title = paste0(trait_name[i2]," for ",eth_name[i1]))+
  #subtitle = "A2: Critically ill COVID19+ vs. population controls;\nB1: Hospitalized COVID19+ vs non-hospitalized COVID19+;\nB2: Hospitalized COVID19+ vs. population controls;\nC2: Reported SARS-CoV-2 infection vs. population controls") + 
  theme_Publication()+
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8)
  )
if(i3==1){
  outpath <-"/data/zhangh24/multi_ethnic/result/cleaned/summary_plot/mega_hap"  
}else{
  outpath <-"/data/zhangh24/multi_ethnic/result/cleaned/summary_plot/im"  
}


ggsave(filename=paste0("man_",eth[i1],"_",trait[i2],".png"),
       plot=manhplot, device="png",
       path=outpath,
       width=9, height=4, units="in", dpi=300)

#qq(dat$P)










library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")



# option_list <- list(
#   make_option("--input", type="character", default="",
#               help="Input file, tab delimited; required columns: 'MAF' and 'PVALUE'"),
#   make_option("--prefix", type="character", default="",
#               help="Prefix of output files"),
#   make_option("--top.size", type="numeric", default=0.125,
#               help="top size = proportion of total length y axis [default=0.125]"),
#   make_option("--break.top", type="numeric", default=15,
#               help="set axis break at -log10(P) [default=15]"),
#   make_option("--width", type="numeric", default=900,
#               help="Width QQ plot in pixel [default=900]"),
#   make_option("--height", type="numeric", default=900,
#               help="Height QQ plot in pixel [default=900]"),
#   make_option("--pointsize", type="numeric", default=16,
#               help="Point size of plots [default=16]"),
#   make_option("--maf", type="character", default="",
#               help="name of column with MAF [default='']"),
#   make_option("--af", type="character", default="",
#               help="name of column with AF [default='']"),
#   make_option("--pvalue", type="character", default="PVALUE",
#               help="name of column with p.value [default='PVALUE']"),
#   make_option("--log10p", type="logical", default=F,
#               help="Input p.value column with -log10(p.value) [default=F]"),
#   make_option("--maintitle", type="character", default="",
#               help="Plot title")
# )
# 
# parser <- OptionParser(usage="%prog [options]", option_list=option_list)
# # 
# args <- parse_args(parser, positional_arguments = 0)
# opt <- args$options


qqplotdata <- function(logpvector){
  o = sort(logpvector,decreasing=T)
  e = -log10(ppoints(length(o)))       
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]
  
  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}
yLine <- c(-log10(5E-8))
colLine <- c("red")
dat$log10P = -log10(dat$P)
gwas = as.data.frame(dat)
# Determine frequency bins and create variable for binned QQ plot

minMAF <- min(gwas$MAF)

freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
gwas$freqbin <- cut(gwas$MAF, freqbins,include.lowest=T)
freqtable <- table(gwas$freqbin)
freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
freqtable <- freqtable[freqtable > 0]

## Generate QQ plot data by frequency bin
fbin <- character(0)
fN <- integer(0)
fx <- numeric(0)
fy <- numeric(0)
fcol <- character(0)
legendcol <- character(0)
conf <- list()
allcols <- brewer.pal(4,"Set1")
ycol <- "log10P"
for(f in 1:length(freqtable)){
  fbin <- c(fbin,names(freqtable)[f])
  fsnps <- which(gwas$freqbin ==names(freqtable)[f])
  plotdata <- qqplotdata(gwas[[ycol]][fsnps])
  fN <- c(fN,freqtable[f])
  fx <- c(fx,plotdata$e)
  fy <- c(fy,plotdata$o)
  fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
  conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                          'y'=c(plotdata$c975,rev(plotdata$c025)))
  legendcol <- c(legendcol,allcols[f])
}
legendtext <- paste0("MAF=",fbin,"; N SNPs=",format(fN,big.mark=",",scientific=FALSE))
opt =  list(break.top = 15,
            top.size = 0.125)


png(filename = paste0(outpath,"/QQ_",eth[i1],"_",trait[i2],".png"), width = 8, height = 8, units = "in",res=300)
xlim <- c(0,max(fx,na.rm=T))
ylim <- c(0,max(fy,na.rm=T))
maxY <- max(fy,na.rm=T)
print("okkkk2")
par(mar=c(5.1,5.1,4.1,1.1))
print("okkkk3")

lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
lab2 <- lab2[lab2 > max(lab1)]

# resulting range of top scale in bottom scale units
top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
top.data = max(lab2)-opt$break.top

# function to rescale the top part
rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
rescaled.y = rescale(fy[fy>opt$break.top])
plot(0,0,
     ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
     xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
     ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
     cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
     main=opt$maintitle,pch=19)

# Plot confidence intervals	
for(p in 1:length(conf)){
  polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
          col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
          border = NA)
}

# add points
points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)

# identify line & add axis break
lines(xlim,xlim,col="black",lty = 2)
axis(1,cex.axis=1.5,cex.lab=1.5)
par(las=1)
axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
par(las=0)
box()
par(las=0)
points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
par(las=1)
axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
abline(h=ifelse(yLine<opt$break.top,
                yLine,
                rescale(yLine)),
       col=colLine,lwd=1.5,lty=2)
legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
text(5,1,expression(paste(lambda[1000]," = ")),cex = 1.5)
text(5.7,1,paste(lambda_1000),cex = 1.5)

title(paste0(trait_name[i2]," for ",eth_name[i1]))
dev.off()

# plot(1,1)
# text(1,0.8,expression(paste(lambda[1000]," = ",buquote(.(lambda_1000)))),cex = 1.5)

lambda_vec = c(lambda,lambda_1000)
save(lambda_vec,file = paste0("/data/zhangh24/multi_ethnic/result/cleaned/lambda_value/lambda_vec_",i1,"_",i2,"_",i3,".rdata"))








head(sample_size)
sample_size = sample_size %>% 
  mutate(N_case = ifelse(is.na(N_case),0,N_case))
sample_size = sample_size %>% 
  mutate(N_total = N_control + N_case)
avg_sample = sample_size %>% 
  group_by(eth) %>% 
  summarize(avg_N = mean(N_total))
