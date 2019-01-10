library(ggplot2)
library(plyr)
library(reshape2)
p65.tt <- read.table("p65_distr_across_labeled.bed", header = T)
STAT6.tt <- read.table("STAT6_distr_across_labeled.bed", header = T)

rownames(p65.tt) <- p65.tt$Factor
rownames(STAT6.tt) <- STAT6.tt$Factor

ord <- order(STAT6.tt$Count, decreasing = F)
STAT6.tt <- STAT6.tt[ord,]
p65.tt <- p65.tt[rownames(STAT6.tt),]

tt <- rbind(cbind(STAT6.tt, TF = "STAT6"), cbind(p65.tt, TF = "p65"))
tt$Factor <- factor(tt$Factor,levels = STAT6.tt$Factor)
tt$TF <- factor(tt$TF,levels = c("STAT6", "p65"), ordered = T)

gg <- ggplot(data=tt) +
geom_bar(aes(Factor,Count,group=TF,fill=TF), stat = "identity",subset(tt,tt$TF=="p65")) +
geom_bar(aes(Factor,-Count,group=TF,fill=TF), stat = "identity",subset(tt,tt$TF=="STAT6")) +
ylim(-1300,1300) +
coord_flip() + theme_bw() + scale_fill_manual(values = c("red", "blue"))

ggsave(gg, filename = "balance.pdf")

percent.tt <- tt
percent.tt[percent.tt$TF == "STAT6","Count"] = subset(percent.tt, TF == "STAT6")$Count/sum(subset(percent.tt, TF == "STAT6")$Count)
percent.tt[percent.tt$TF == "p65","Count"] = subset(percent.tt, TF == "p65")$Count/sum(subset(percent.tt, TF == "p65")$Count)
percent.gg <- ggplot(data=percent.tt) +
geom_bar(aes(Factor,-Count,group=TF,fill=TF), stat = "identity",subset(percent.tt,percent.tt$TF=="STAT6")) +
geom_bar(aes(Factor,Count,group=TF,fill=TF), stat = "identity",subset(percent.tt,percent.tt$TF=="p65")) +
ylim(-0.23,0.23) +
coord_flip() + theme_bw() + scale_fill_manual(values = c("red", "blue"))

ggsave(percent.gg, filename = "balance_percent.pdf")
