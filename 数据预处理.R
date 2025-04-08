rm(list=ls())

setwd("D:\\KY\\SZ")

library(openxlsx)

## 蛋白数据读取
pro <- read.xlsx("SZ_PRO.xlsx")
rownames(pro) <- pro[,1]
pro <- pro[,-1]
pro_info <- pro[,c(1,107:117)]
data_pro <- data.frame(t(pro[,-c(1,107:117)]))

#write.csv(pro_info,"pro_info.csv")
#write.csv(data_pro,"pro_raw_expr.csv")

## PTM数据读取
ptm <- read.xlsx("SZ_PTM.xlsx")
rownames(ptm) <- paste0("ptm", 1:2664)
ptm_info <- ptm[,c(1:2,108:116)]
data_ptm <- data.frame(t(ptm[,-c(1:2,108:116)]))

#write.csv(ptm_info,"ptm_info.csv")
#write.csv(data_ptm,"ptm_raw_expr.csv")

## META数据读取
meta <- read.xlsx("SZ_META.xlsx")
rownames(meta) <- paste0("meta", 1:1535)
meta_info <- meta[,c(1:10,115:122)]
data_meta <- data.frame(t(meta[,-c(1:10,115:122)]))

#write.csv(meta_info,"meta_info.csv")
#write.csv(data_meta,"meta_raw_expr.csv")

library(missForest)

## 补缺失值
set.seed(123)  
imputed_pro <- missForest(data_pro)
imputed_ptm <- missForest(data_ptm)
imputed_meta <- missForest(data_meta)

# 获取插补后的数据
completed_pro <- imputed_pro$ximp
completed_ptm <- imputed_ptm$ximp
completed_meta <- imputed_meta$ximp

#write.csv(completed_pro,"pro_nona.csv")
#write.csv(completed_ptm,"ptm_nona.csv")
#write.csv(completed_meta,"ptm_nona.csv")

## 标准化
pro_scale <- data.frame(scale(completed_pro))
ptm_scale <- data.frame(scale(completed_ptm))
meta_scale <- data.frame(scale(completed_meta))

#write.csv(pro_scale,"pro_nona_scale.csv")
#write.csv(ptm_scale,"ptm_nona_scale.csv")
#write.csv(meta_scale,"ptm_nona_scale.csv")

## 合并三个组学数据
m0 <- merge(pro_scale, ptm_scale, by = "row.names", all = TRUE) 
rownames(m0)<-m0[,1]
m0<-m0[,-1]
m0<-m0[-which(rownames(m0)=="SZ_F_66.07255"),]
m <- merge(m0, meta_scale, by = "row.names", all = TRUE) 
rownames(m)<-m[,1]
m<-m[,-1]

#write.csv(m,"all_nona_scale.csv")

library(caret)

## 去冗余
cor_matrix <- cor(m) # 计算相关性矩阵
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.75) # 选择与目标变量相关性较强的特征
selected_features <- m[, -highly_correlated]

#write.csv(selected_features,"all_nona_scale_cor.csv")

## 分组
group <- as.data.frame(c(rep(1, 51), rep(0, 53))) #1 NC 0 SZ
colnames(group) <- "label"

data <- cbind(group,m) #4566
data_cor <- cbind(group,selected_features) #2053

group0 <- as.data.frame(c(rep(1, 51), rep(0, 54))) #1 NC 0 SZ
colnames(group0) <- "label"

pro_label <- cbind(group0,pro_scale)
ptm_label <- cbind(group0,ptm_scale)
meta_label <- cbind(group,meta_scale)

write.csv(data,"all_nona_scale_label.csv")
write.csv(data_cor,"all_nona_scale_cor_label.csv")

write.csv(pro_label,"pro_nona_scale_label.csv")
write.csv(ptm_label,"ptm_nona_scale_label.csv")
write.csv(meta_label,"meta_nona_scale_label.csv")























