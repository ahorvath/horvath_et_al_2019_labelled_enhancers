library(e1071)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(randomForest)
library(rfUtilities)

ttt <- NULL
samples <- c("mm_BMDM_IL4_0h_PU1_CS286","mm_BMDM_Irf8_veh_CS1055","mm_BMDM_JunB_veh_CS1054","mm_BMDM_RUNX1_veh", "mm_BMDM_CEBPa_veh")
for (sample.name in samples) {
	ATAC_0h_not_PU1_0h.tt <- read.table(paste0("../boxplots/",sample.name, "_on_ATAC_0h_not_PU1_0h_distal.tsv"), header = F)
	ATAC_0h_not_PU1_0h.tt <- ATAC_0h_not_PU1_0h.tt[sample(1:nrow(ATAC_0h_not_PU1_0h.tt), 1000),]

	PU1_0h_and_ATAC_0h.tt <- read.table(paste0("../boxplots/",sample.name, "_on_PU1_0h_and_ATAC_0h_distal.tsv"), header = F)
	PU1_0h_and_ATAC_0h.tt <- PU1_0h_and_ATAC_0h.tt[sample(1:nrow(PU1_0h_and_ATAC_0h.tt), 1000),]

	PU1_0h_not_ATAC_0h.tt <- read.table(paste0("../boxplots/",sample.name, "_on_PU1_0h_not_ATAC_0h_distal.tsv"), header = F)
	PU1_0h_not_ATAC_0h.tt <- PU1_0h_not_ATAC_0h.tt[sample(1:nrow(PU1_0h_not_ATAC_0h.tt), 1000),]

	tt <- c(ATAC_0h_not_PU1_0h.tt[,4], PU1_0h_and_ATAC_0h.tt[,4], PU1_0h_not_ATAC_0h.tt[,4])
	print(head(tt))

	ttt <- cbind(ttt, tt)
}
#colnames(ttt) <- samples
colnames(ttt) <- c("PU1", "IRF8", "JUNB", "RUNX1", "CEBPA")
ttt <- as.data.frame(cbind(ttt, Group = rep(c("ATAC_0h_not_PU1_0h", "PU1_0h_and_ATAC_0h", "PU1_0h_not_ATAC_0h"), c(1000,1000,1000))))

ttt[,1:5] <- apply(ttt[,1:5], 1:2, as.numeric)

tttt <- ttt
tttt$Group <- factor(tttt$Group, levels = c("PU1_0h_and_ATAC_0h", "ATAC_0h_not_PU1_0h", "PU1_0h_not_ATAC_0h"))

ttt.PU1.cistrome <- subset(tttt, Group != "ATAC_0h_not_PU1_0h")
ttt.PU1.cistrome$Group <- factor(ttt.PU1.cistrome$Group)

null.ttt.PU1.cistrome <- ttt.PU1.cistrome[,c(5,6)]
rf.null.ttt.PU1.cistrome <- randomForest(Group ~  ., data = null.ttt.PU1.cistrome, importance = T, ntree = 1000, cross = 10)
null.model.acc <- sum(diag(table(rf.null.ttt.PU1.cistrome$predicted, null.ttt.PU1.cistrome$Group)))/sum(table(rf.null.ttt.PU1.cistrome$predicted, null.ttt.PU1.cistrome$Group))
result.table <- c(1, null.model.acc, 0)

set.seed(8888)
for (number.of.partners in 1:4) {
	cc <- combn(4, number.of.partners)
	for (col.cc in 1:ncol(cc)) {
		partner.index <- cc[,col.cc] + 1
		current.ttt.PU1.cistrome <- ttt.PU1.cistrome[,c(1, partner.index, 6)]
		print(partner.index)
		rf.current.ttt.PU1.cistrome <- randomForest(Group ~  ., data = current.ttt.PU1.cistrome, importance = T, ntree = 1000, cross = 10)
		current.model.acc <- sum(diag(table(rf.current.ttt.PU1.cistrome$predicted, current.ttt.PU1.cistrome$Group)))/sum(table(rf.current.ttt.PU1.cistrome$predicted, current.ttt.PU1.cistrome$Group))
		result.table <- rbind(result.table,c(paste(partner.index, collapse = "_"), current.model.acc, nrow(cc)))
	}
}

df.result.table <- data.frame(result.table)
df.result.table[,2:3] <- apply(df.result.table[,2:3],1:2,as.numeric)
colnames(df.result.table) <- c("Partners", "Accuracy", "Size")
df.result.table$Size <- factor(df.result.table$Size)
gg <- ggplot(df.result.table, aes(Size, Accuracy, fill = Size)) + geom_boxplot(outlier.size = 0) + geom_point(aes(y = Accuracy), shape = 21, size = 6, position = position_dodge(width = 0.75)) + theme_bw() + scale_fill_brewer(palette = "Reds") + scale_colour_brewer(palette = "Reds") + theme(axis.ticks.x = element_blank())

ggsave(gg, filename = "additive_boxplots.jpg")

rf.PU1.cistrome.tttt <- rf.current.ttt.PU1.cistrome

PU1.cistrome.importance.tt <- rf.PU1.cistrome.tttt$importance[,4, drop = F]
PU1.cistrome.importance.tt <- PU1.cistrome.importance.tt[order(PU1.cistrome.importance.tt[,1], decreasing = F),, drop = F]
PU1.cistrome.importance.tt <- as.data.frame(cbind(Factor = rownames(PU1.cistrome.importance.tt), PU1.cistrome.importance.tt))
PU1.cistrome.importance.tt[,2] <- as.numeric(as.character(PU1.cistrome.importance.tt[,2]))

PU1.cistrome.importance.tt$Factor <- factor(PU1.cistrome.importance.tt$Factor, levels = PU1.cistrome.importance.tt$Factor)

PU1.cistrome.gg <- ggplot(PU1.cistrome.importance.tt, aes(Factor, MeanDecreaseGini)) + geom_bar(stat = "identity", fill = "azure4", width = 0.6) + theme_bw() +  coord_flip()

ggsave(PU1.cistrome.gg, filename = "PU1_cistrome_openness.jpg")

PU1.cistrome.importance.tt <- rf.PU1.cistrome.tttt$importance[,3, drop = F]
PU1.cistrome.importance.tt <- PU1.cistrome.importance.tt[order(PU1.cistrome.importance.tt[,1], decreasing = F),, drop = F]
PU1.cistrome.importance.tt <- as.data.frame(cbind(Factor = rownames(PU1.cistrome.importance.tt), PU1.cistrome.importance.tt))
PU1.cistrome.importance.tt[,2] <- as.numeric(as.character(PU1.cistrome.importance.tt[,2]))

PU1.cistrome.importance.tt$Factor <- factor(PU1.cistrome.importance.tt$Factor, levels = PU1.cistrome.importance.tt$Factor)

PU1.cistrome.gg <- ggplot(PU1.cistrome.importance.tt, aes(Factor, MeanDecreaseAccuracy)) + geom_bar(stat = "identity", fill = "#E41A1C") + theme_bw() +  coord_flip()

ggsave(PU1.cistrome.gg, filename = "PU1_cistrome_openness_meanDec.jpg")

