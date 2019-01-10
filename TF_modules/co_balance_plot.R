library(ggplot2)
library(plyr)
library(reshape2)
tt <- read.table("co_labeled.bed", header = T)
rownames(tt) <- tt$Factor

ord <- order(tt$Count, decreasing = F)
tt <- tt[ord,]

tt$Factor <- factor(tt$Factor,levels = tt$Factor)

gg <- ggplot(data=tt) +
geom_bar(aes(Factor,Count,fill="red"), stat = "identity") +
ylim(0,35000) +
coord_flip() + theme_bw() + scale_fill_manual(values = c("red"))

ggsave(gg, filename = "co_balance.pdf")


gg2 <- ggplot(data=tt[-c(11:13,15),]) +
geom_bar(aes(Factor,Count,fill="red"), stat = "identity") +
ylim(0,6000) +
#coord_flip() +
 theme_bw() + scale_fill_manual(values = c("red"))

ggsave(gg2, filename = "co_balance_2.pdf")

gg <- ggplot(data=tt[c(11:13,15),]) +
geom_bar(aes(Factor,Count,fill="red"), stat = "identity") +
ylim(0,35000) +
coord_flip() + theme_bw() + scale_fill_manual(values = c("red"))

ggsave(gg, filename = "co_balance_onlies.pdf")
