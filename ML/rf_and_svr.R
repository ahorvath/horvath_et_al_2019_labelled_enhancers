library(e1071)
library(ROCR)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(randomForest)
library(rfUtilities)
library(penalizedSVM)
library(MASS)
library(ggplot2)
library(pheatmap)
library(caret)

rmse <- function(error) {
  sqrt(mean(error^2))
}

ATAC.high.acc.signal.table <- read.table("../../50M_PU1_v10/NFR_vs_LDTFs/ML_close_vs_open_openness/boxplots/mm_BMDM_ATACseq_veh_1_on_five_factors_union_and_ATAC_0h_distal.tsv")
ATAC.low.acc.signal.table <- read.table("../../50M_PU1_v10/NFR_vs_LDTFs/ML_close_vs_open_openness/boxplots/mm_BMDM_ATACseq_veh_1_on_five_factors_union_not_ATAC_0h_distal.tsv")

rownames(ATAC.low.acc.signal.table) <- c(paste0(ATAC.low.acc.signal.table[,1],":",ATAC.low.acc.signal.table[,2] ,"-", ATAC.low.acc.signal.table[,3]))
rownames(ATAC.high.acc.signal.table) <- c(paste0(ATAC.high.acc.signal.table[,1],":",ATAC.high.acc.signal.table[,2] ,"-", ATAC.high.acc.signal.table[,3]))
colnames(ATAC.low.acc.signal.table)[4] <- "ATAC"
colnames(ATAC.high.acc.signal.table)[4] <- "ATAC"

high.acc <- read.table("../../50M_PU1_v10/NFR_vs_LDTFs/ML_close_vs_open_openness/high_distals_wo_repeats.bed")
labeled <- read.table("../../50M_PU1_v10/NFR_vs_LDTFs/ML_close_vs_open_openness/low_distals_wo_repeats.bed")

rownames(labeled) <- c(paste0(labeled[,1],":",labeled[,2] ,"-", labeled[,3]))
rownames(high.acc) <- c(paste0(high.acc[,1],":",high.acc[,2] ,"-", high.acc[,3]))
high.acc <- high.acc[,-c(1:3)]
labeled <- labeled[,-c(1:3)]
colnames(high.acc) <- c("PU1", "IRF8", "JunB", "CEBPA", "RUNX1")
colnames(labeled) <- c("PU1", "IRF8", "JunB", "CEBPA", "RUNX1")

train.high.acc.sites <- sample(rownames(high.acc), 1000)
train.labeled.sites <- sample(rownames(labeled), 1000)

train.labeled.vs.high.acc.cistrome <- rbind(cbind(labeled[train.labeled.sites,  1:5], ATAC = ATAC.low.acc.signal.table[train.labeled.sites,4], Group = rep("labeled", 1000)), 
				            cbind(high.acc[train.high.acc.sites, 1:5], ATAC = ATAC.high.acc.signal.table[train.high.acc.sites,4], Group = rep("HighAcc",1000)))

train.labeled.vs.high.acc.cistrome$Group <- factor(train.labeled.vs.high.acc.cistrome$Group,  levels = c("HighAcc", "labeled")) 

valid.labeled.sites <- sample(setdiff(rownames(labeled), train.high.acc.sites), 1000)
valid.high.acc.sites <- sample(setdiff(rownames(high.acc), train.labeled.sites), 1000)

valid.labeled.vs.high.acc.cistrome <- rbind(cbind(labeled[valid.labeled.sites, 1:5], ATAC = ATAC.low.acc.signal.table[valid.labeled.sites,4], Group = rep("labeled", 1000)),

valid.labeled.vs.high.acc.cistrome$Group <- factor(valid.labeled.vs.high.acc.cistrome$Group,  levels = c("HighAcc", "labeled")) 


class.rf.labeled.vs.high <- randomForest(log2(train.labeled.vs.high.acc.cistrome[,1:5]+1), train.labeled.vs.high.acc.cistrome$Group, cross = 10)

save.image()



pheatmap(table(class.rf.labeled.vs.high$predicted,train.labeled.vs.high.acc.cistrome$Group)/2000)

importance.tt <- class.rf.labeled.vs.high$importance[,1, drop = F]
importance.tt <- importance.tt[order(importance.tt[,1], decreasing = F),, drop = F]
importance.tt <- as.data.frame(cbind(Factor = rownames(importance.tt), importance.tt))
importance.tt[,2] <- as.numeric(as.character(importance.tt[,2]))

importance.tt$Factor <- factor(importance.tt$Factor, levels = importance.tt$Factor)

gg.contr <- ggplot(importance.tt, aes(Factor, MeanDecreaseGini)) + geom_bar(stat = "identity", fill = "azure4", width = 0.6) + theme_bw()

ggsave(gg.contr, filename = "union_importance.pdf")

train.predictions = predict(class.rf.labeled.vs.high,type="class",newdata=log2(train.labeled.vs.high.acc.cistrome[,1:5]+1))
valid.predictions = predict(class.rf.labeled.vs.high,type="class",newdata=log2(valid.labeled.vs.high.acc.cistrome[,1:5]+1))

options(scipen = 999, digits=10)

cc.train <- confusionMatrix(class.rf.labeled.vs.high$predicted, train.labeled.vs.high.acc.cistrome$Group)
cc.valid <- confusionMatrix(valid.predictions, valid.labeled.vs.high.acc.cistrome$Group)

predictions=predict(class.rf.labeled.vs.high,type="prob",newdata=log2(valid.labeled.vs.high.acc.cistrome[,1:5]+1))[,2]    
pred=prediction(predictions,valid.labeled.vs.high.acc.cistrome$Group)

perf_AUC=performance(pred,"auc") #Calculate the AUC value

AUC=perf_AUC@y.values[[1]]

perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve

df <- data.frame(TPR = perf_ROC@y.values[[1]], FPR = perf_ROC@x.values[[1]])

gg <- ggplot(df, aes(FPR, TPR)) + geom_line(size = 1, alpha = 1, color = "blue") + labs(x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)") + theme_bw() 
ggsave(gg, filename = "ROC_valid.pdf")

#class.rf.nulabeled.ttt.PU1.cistrome <- randomForest(ttt.close_vs_open.cistrome[,7:11], ttt.close_vs_open.cistrome$Group, cross = 10)
#svm.nulabeled.ttt.PU1.cistrome <- svm(Group ~  ., data = ttt.close_vs_open.cistrome[,7:12], cross = 10)
svr.labeled.vs.high <- svm(ATAC ~ ., log2(train.labeled.vs.high.acc.cistrome[,c(1:6)]+1), cross = 10)


cor(svr.labeled.vs.high$fitted, log2(train.labeled.vs.high.acc.cistrome[,6]+1))

predicted.values <- predict(svr.labeled.vs.high, log2(valid.labeled.vs.high.acc.cistrome[,c(1:5)]+1))

valid.cor <- cor(predicted.values, log10(valid.labeled.vs.high.acc.cistrome$ATAC+1)) 
valid.error <- predicted.values  - log10(valid.labeled.vs.high.acc.cistrome$ATAC+1)

#dat <- as.data.frame(cbind(measured = log10(train.labeled.vs.high.acc.cistrome$ATAC+1), fitted = svr.labeled.vs.high$fitted))
dat <- as.data.frame(cbind(measured = log2(valid.labeled.vs.high.acc.cistrome$ATAC+1), fitted = predicted.values))
gg <- ggplot(data = dat, aes(measured, fitted)) + stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) + scale_fill_continuous(low = "white", high = "dodgerblue4") + xlim(0,6) + ylim(0,6) + theme_bw() + geom_point(alpha = 0.1, shape = 20) 

ggsave(gg, filename = "contour.pdf")
ggsave(gg, filename = "contour.png")

##############

predictions=as.vector(rf.current.ttt.PU1.cistrome$votes[,2])
pred=prediction(predictions,current.ttt.PU1.cistrome$Group)

perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUC=perf_AUC@y.values[[1]]

perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve

df <- data.frame(TPR = perf_ROC@y.values[[1]], FPR = perf_ROC@x.values[[1]])

gg <- ggplot(df, aes(FPR, TPR)) + geom_line(size = 1, alpha = 1, color = "blue") + labs(x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)") + theme_bw()



set.seed(8888)
max.number.of.partners = 5
result.table = NULL
for (number.of.partners in 1:max.number.of.partners) {
#for (number.of.partners in 1) {
        cc <- combn(max.number.of.partners, number.of.partners)
        for (col.cc in 1:ncol(cc)) {
                partner.index <- cc[,col.cc]
                print(partner.index)
                #rf.current.ttt.PU1.cistrome <- randomForest(Group ~  ., data = current.ttt.PU1.cistrome, importance = T, ntree = 1000, cross = 10)
		current.model <- randomForest(log2(train.labeled.vs.high.acc.cistrome[,partner.index, drop = F]+1), train.labeled.vs.high.acc.cistrome$Group, cross = 10)
		current.model.train <- confusionMatrix(current.model$predicted, train.labeled.vs.high.acc.cistrome$Group)

		valid.predictions = predict(current.model,type="class", newdata=log2(valid.labeled.vs.high.acc.cistrome[, partner.index, drop = F]+1))	
		current.model.valid <- confusionMatrix(valid.predictions, valid.labeled.vs.high.acc.cistrome$Group)

		current.predictions = predict(current.model,type="prob",newdata=log2(valid.labeled.vs.high.acc.cistrome[, partner.index, drop = F]+1))[,2]
		current.pred = prediction(current.predictions,valid.labeled.vs.high.acc.cistrome$Group)

		current.perf_AUC = performance(current.pred,"auc") #Calculate the AUC value
		current.AUC = current.perf_AUC@y.values[[1]]		

                result.table <- rbind(result.table,c(paste(partner.index, collapse = "_"), current.model.valid$overall[1], current.AUC, nrow(cc)))
        }
}

df.result.table <- data.frame(result.table)
df.result.table[,2:4] <- apply(df.result.table[,2:4],1:2,as.numeric)
colnames(df.result.table) <- c("Partners", "Accuracy", "AUC", "Size")
df.result.table$Size <- factor(df.result.table$Size)
gg <- ggplot(df.result.table, aes(Size, Accuracy, fill = Size)) + geom_boxplot(outlier.size = 0) + geom_point(aes(y = Accuracy), shape = 21, size = 6, position = position_dodge(width = 0.75)) + theme_bw() + scale_fill_brewer(palette = "Reds") + ylim(0.6,0.83) + scale_colour_brewer(palette = "Reds") + theme(axis.ticks.x = element_blank())

ggsave(gg, filename = "additive_boxplots.pdf")





