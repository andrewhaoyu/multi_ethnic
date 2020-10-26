library(ggplot2)
myleg<-read.csv(text="lett,num
a,1
a,2
b,1
b,2
h,1
h,2
h,3
h,4")
ggplot(myleg,aes(lett,alpha=factor(num),fill=lett)) +geom_bar(position=position_stack(reverse=T)) +scale_alpha_discrete(range=c(1,.1), name="alpha legend",labels=c("alpha lab 4","alpha lab 3","alpha lab 2", "alpha lab 1")) +labs(title="initial bar plot for data")

library(ggplot2)
library(grid)
library(gridExtra)
myleg <- structure(list(lett = structure(c(1L, 1L, 2L, 2L, 3L, 3L, 3L,
                                           3L), .Label = c("a", "b", "h"), class = "factor"), num = c(1L, 
                                                                                                      2L, 1L, 2L, 1L, 2L, 3L, 4L)), .Names = c("lett", "num"), 
                   class = "data.frame", row.names = c(NA, -8L))

getLegend <- function(p) {
  g <- ggplotGrob(p)
  k <- which(g$layout$name=="guide-box")
  g$grobs[[k]]
}

p1 <- ggplot(myleg,aes(lett,alpha=factor(num),fill=lett)) +geom_bar(position="stack",fill="#f8766d") +scale_alpha_discrete(name="red legend",labels=c("red lab 2","red lab 1"),breaks=c("3","4"))
p2 <- ggplot(myleg,aes(lett,alpha=factor(num),fill=lett)) +geom_bar(position="stack",fill="#00ba38") +scale_alpha_discrete(name="green legend",labels=c("green lab 2","green lab 1"),breaks=c("3","4"))
p3 <- ggplot(myleg,aes(lett,alpha=factor(num),fill=lett)) +geom_bar(position="stack",fill="#619cff") +scale_alpha_discrete(name="blue legend",labels=c("blue lab 4","blue lab 3","blue lab 2", "blue lab 1"))

p <- ggplot(myleg,aes(lett,alpha=factor(num),fill=lett)) +
  geom_bar(position=position_stack(reverse=T)) +
  scale_alpha_discrete(range=c(1,.1), name="alpha legend",
                       labels=c("alpha lab 4","alpha lab 3","alpha lab 2", "alpha lab 1")) +
  labs(title="initial bar plot for data")
g <- ggplotGrob(p)

k <- which(g$layout$name=="guide-box")
g$grobs[[k]] <- grid.arrange(getLegend(p1),getLegend(p2),getLegend(p3),ncol=1)
grid.draw(g)
