
# library
library(dplyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(scales)

# load
load("Plasmids_carrying_Tn_IS_Integron_number_Length.RData")
load("sample_name_sampleName_df.RData")
load("sample_name_number_df.RData")
sample_name_sampleName_df$sampleType <- sapply(strsplit(sample_name_sampleName_df$sample_name, "_", fixed = T), "[[", 2)

# 
plotDat_ARG <- plotDat  # 1637 (352*3*3 = 3168) 

# 
plotDat_ARG_width <- dcast(plotDat_ARG, sampleName + MGE~Plasmid_type,
                           value.var = "normalization_number")
plotDat_ARG_width[is.na(plotDat_ARG_width)] <- 0

plotDat_ARG_S <- melt(plotDat_ARG_width, variable.name = "Plasmid_type", value.name = "normalization_number")  # 2913
plotDat_ARG_S$Plasmid_type <- as.character(plotDat_ARG_S$Plasmid_type) 

# 
plotDat_ARG_width <- dcast(plotDat_ARG_S, sampleName + Plasmid_type ~ MGE,
                           value.var = "normalization_number")
plotDat_ARG_width[is.na(plotDat_ARG_width)] <- 0 # 352*3=1056

plotDat_ARG_S <- melt(plotDat_ARG_width, variable.name = "MGE", 
                      value.name = "normalization_number")  # 3168
plotDat_ARG_S$MGE <- as.character(plotDat_ARG_S$MGE)

#### ---------------
plotDat_ARG_S$sampleType <- sapply(plotDat_ARG_S$sampleName,
                                   function(x) sample_name_sampleName_df$sampleType[which(sample_name_sampleName_df$sampleName == x)])

# factor
plotDat_ARG_S$MGE <- factor(plotDat_ARG_S$MGE, levels = c("Plasmid_IS","Plasmid_Transposon",  "Plasmid_Integron"))
plotDat_ARG_S$sampleType <- factor(plotDat_ARG_S$sampleType, 
                                   levels = c("tailings", "Sludge", "farmland", "forest", "grass", "gobi"))


# 构建一个函数，以一位小数点的科学计数法显示
format_science <- function(x) {
  formatC(x, format = "e", digits = 1)
}

# A.Conjugative plasmid --------------------------------------------------------
plotDat_ARG_conj <- plotDat_ARG_S %>% filter(Plasmid_type == "Conj") # 

P_seperate_conj <- ggboxplot(data = plotDat_ARG_conj, x = "MGE", y = "normalization_number", alpha=0.7, width = 0.8, size = 0.45, alpha=0.4, color="black",outlier.shape = NA) +
  geom_jitter(aes(x = MGE, y = normalization_number,color = MGE),  width = 0.2, size=2.5, alpha=0.45) + 
  #geom_jitter(width = 0.2, shape = 21, aes(fill = MGE), alpha = 0.7, size =3.5, color = "white") +
  #geom_boxplot(outlier.shape = NA, aes(color = MGE), size = 0.8,position = position_nudge(x= 0.1, y= 0),side = "R", adjust = 1.2, fill = NA) +
  facet_wrap(.~sampleType, nrow = 1, scale = "free_y") + 
  scale_fill_manual(values = c('#fb8c00','#689f38', "#c0ca33")) +
  scale_color_manual(values = c('#fb8c00','#689f38', "#c0ca33")) +
  theme_pubr(base_size = 14, border = F) +
  xlab("") +ylab("") + 
  #ylim(0, 20000) +
  theme(legend.position = "right")   + 
  theme(axis.ticks = element_line(colour = "gray0"),
        axis.text = element_text(colour = "gray0"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  #stat_compare_means(label.y = c(7), label.x = c(2.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.x = element_blank(), #去掉X轴坐标
        axis.ticks.x = element_blank()) + 
  scale_y_continuous(labels = format_science) 
#coord_cartesian(ylim = c(1.0e-06, 3.3e-05)) 

P_seperate_conj


# B.Mobilizable plasmid --------------------------------------------------------
plotDat_ARG_mob <- plotDat_ARG_S %>% filter(Plasmid_type == "mob_unconj") 

P_seperate_mob <- 
  ggboxplot(data = plotDat_ARG_mob, x = "MGE", y = "normalization_number", alpha=0.7, width = 0.8, size = 0.45, alpha=0.4, color="black",outlier.shape = NA) +
  geom_jitter(aes(x = MGE, y = normalization_number,color = MGE),  width = 0.2, size=2.5, alpha=0.45) + 
  #geom_jitter(width = 0.2, shape = 21, aes(fill = MGE), alpha = 0.7, size =3.5, color = "white") +
  #geom_boxplot(outlier.shape = NA, aes(color = MGE), size = 0.8,position = position_nudge(x= 0.1, y= 0),side = "R", adjust = 1.2, fill = NA) +
  facet_wrap(.~sampleType, nrow = 1, scale = "free_y") + 
  scale_fill_manual(values = c('#fb8c00','#689f38', "#c0ca33")) +
  scale_color_manual(values = c('#fb8c00','#689f38', "#c0ca33")) +
  theme_pubr(base_size = 14, border = F) +
  xlab("") +ylab("") + 
  #ylim(0, 20000) +
  theme(legend.position = "right")   + 
  theme(axis.ticks = element_line(colour = "gray0"),
        axis.text = element_text(colour = "gray0"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  #stat_compare_means(label.y = c(7), label.x = c(2.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.x = element_blank(), #去掉X轴坐标
        axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = format_science) 
#coord_cartesian(ylim = c(1.0e-06, 3.3e-05)) 

P_seperate_mob

# c.nonMobilizable plasmid -----------------------------------------------------
plotDat_ARG_unmob <- plotDat_ARG_S %>% filter(Plasmid_type == "unmob") 

# seperate
P_seperate_unmob <- ggboxplot(data = plotDat_ARG_unmob, x = "MGE", y = "normalization_number", alpha=0.7, width = 0.8, size = 0.45, alpha=0.4, color="black",outlier.shape = NA) +
  geom_jitter(aes(x = MGE, y = normalization_number,color = MGE),  width = 0.2, size=2.5, alpha=0.45) + 
  #geom_jitter(width = 0.2, shape = 21, aes(fill = MGE), alpha = 0.7, size =3.5, color = "white") +
  #geom_boxplot(outlier.shape = NA, aes(color = MGE), size = 0.8,position = position_nudge(x= 0.1, y= 0),side = "R", adjust = 1.2, fill = NA) +
  facet_wrap(.~sampleType, nrow = 1, scale = "free_y") + 
  scale_fill_manual(values = c('#fb8c00','#689f38', "#c0ca33")) +
  scale_color_manual(values = c('#fb8c00','#689f38', "#c0ca33")) +
  theme_pubr(base_size = 14, border = F) +
  xlab("") +ylab("") + 
  #ylim(0, 20000) +
  theme(legend.position = "right")   + 
  theme(axis.ticks = element_line(colour = "gray0"),
        axis.text = element_text(colour = "gray0"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  #stat_compare_means(label.y = c(7), label.x = c(2.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.x = element_blank(), #去掉X轴坐标
        axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = format_science) 

#coord_cartesian(ylim = c(1.0e-06, 3.3e-05))  
P_seperate_unmob

