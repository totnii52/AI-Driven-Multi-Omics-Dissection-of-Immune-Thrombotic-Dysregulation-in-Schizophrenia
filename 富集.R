library(clusterProfiler)
library(org.Hs.eg.db)  # 如果分析人类数据

gene = c("IGKC","F10","C9","CD5L","IGHG1","VWF","DPP4","CFI","F2","PLG","A1BG","KRT85","FETUB","TF","ORM2","OAF")

# 假设 gene_list 是基因向量
enrich_BP <- enrichGO(gene          = gene,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "SYMBOL",  # 基因名类型
                      ont           = "BP",     # 分析 GO 的类型：BP, MF, CC
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

#Net <- enrichmentNetwork(enrich_BP@result, 
#                         repelLabels = TRUE,        # 尝试避免标签重叠
#                         drawEllipses = TRUE,  
#                         minClusterSize = 3,
#                         fontSize = 3)         

#ggsave("D:\\KY\\multi_omics_SZ\\enrichment_network_BP_plot.pdf", plot = Net, width = 8, height = 6)
# 转换为数据框
enrich_BP_df <- as.data.frame(enrich_BP)

# 按 Gene Count 从大到小排序
enrich_BP_df <- enrich_BP_df[order(-enrich_BP_df$Count),]


# 设置保存路径
png(file = "D:\\KY\\multi_omics_SZ\\BP.png", width = 800, height = 400)

# 创建ggplot
ggplot(enrich_BP_df[1:10, ], aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "#9B59B6", high = "#D55E00") +  # 蓝色到橙色的渐变
  labs(x = "GO Term", y = "Gene Count", fill = "p-value Adjusted") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18),  # 放大文字
        axis.title = element_text(size = 20),  # 放大标题
        plot.title = element_text(hjust = 0.5, size = 20),  # 放大标题
        legend.title = element_text(size = 14),  # 放大图例标题
        legend.text = element_text(size = 14)) +  # 放大图例文字
  ggtitle("BP")

# 关闭PDF设备
dev.off()


# 设置保存路径
pdf(file = "D:\\KY\\multi_omics_SZ\\enrich_go_BP.pdf", width = 8, height = 8)
# 可视化
barplot(enrich_BP, showCategory=15)
dev.off()

enrich_MF <- enrichGO(gene          = gene,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "SYMBOL",  # 基因名类型
                      ont           = "MF",     # 分析 GO 的类型：BP, MF, CC
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

# 设置保存路径
pdf(file = "D:\\KY\\multi_omics_SZ\\enrich_go_MF.pdf", width = 8, height = 8)
# 可视化
barplot(enrich_MF, showCategory=15)
dev.off()

# 转换为数据框
enrich_MF_df <- as.data.frame(enrich_MF)

# 按 Gene Count 从大到小排序
enrich_MF_df <- enrich_MF_df[order(-enrich_MF_df$Count),]

# 设置保存路径
png(file = "D:\\KY\\multi_omics_SZ\\MF.png", width = 700, height = 400)

# 创建ggplot
ggplot(enrich_MF_df[1:10, ], aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "#56B4E9", high = "#E69F00") +  # 蓝色到橙色的渐变
  labs(x = "GO Term", y = "Gene Count", fill = "p-value Adjusted") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18),  # 放大文字
        axis.title = element_text(size = 20),  # 放大标题
        plot.title = element_text(hjust = 0.5, size = 20),  # 放大标题
        legend.title = element_text(size = 14),  # 放大图例标题
        legend.text = element_text(size = 14)) +  # 放大图例文字
  ggtitle("MF")

# 关闭PDF设备
dev.off()

enrich_CC <- enrichGO(gene          = gene,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "SYMBOL",  # 基因名类型
                      ont           = "CC",     # 分析 GO 的类型：BP, MF, CC
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

# 设置保存路径
pdf(file = "D:\\KY\\R\\enrich_go_CC.pdf", width = 6, height = 6)
# 可视化
barplot(enrich_CC, showCategory=10)
dev.off()

# 转换为数据框
enrich_CC_df <- as.data.frame(enrich_CC)

# 按 Gene Count 从大到小排序
enrich_CC_df <- enrich_CC_df[order(-enrich_CC_df$Count),]

# 设置保存路径
png(file = "D:\\KY\\multi_omics_SZ\\CC.png", width = 750, height = 400)

# 创建ggplot
ggplot(enrich_CC_df[1:10, ], aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "#1ABC9C", high = "#F1C40F") +  # 蓝色到橙色的渐变
  labs(x = "GO Term", y = "Gene Count", fill = "p-value Adjusted") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18),  # 放大文字
        axis.title = element_text(size = 20),  # 放大标题
        plot.title = element_text(hjust = 0.5, size = 20),  # 放大标题
        legend.title = element_text(size = 14),  # 放大图例标题
        legend.text = element_text(size = 14)) +  # 放大图例文字
  ggtitle("CC")

# 关闭PDF设备
dev.off()