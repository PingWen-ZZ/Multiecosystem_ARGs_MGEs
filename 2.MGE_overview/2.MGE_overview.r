
# library -------------------------------------------
library(dplyr)
library(ggplot2) 
library(ggpubr) 
library(ggprism)
library(gg.gap)
library(ggbreak)
library(Hmisc)
library(data.table)
library(reshape2)
library(stringr)
library(cowplot) 
library(ggpubr) 
library(patchwork)
library(eulerr)

# 1. MGE composition ----------------------------------------
MGEs_diversity <- fread("MGEs_diversity.txt")

p1 <- ggplot(MGEs_diversity_basedSequence) +
  aes(x = MGEs, y = MGEnumber, fill = MGEs ,color = MGEs) +
  geom_bar(alpha = 0.95 ,width = 0.8, stat = "identity") +
  scale_fill_manual(values = c("#e0e0e0", "#1976d2","#dce775","#f1e54c",  "#2e7d32" ,"#d0c0a5","#d7660d")) +
  scale_color_manual(values = c("#e0e0e0", "#1976d2","#dce775","#f1e54c",  "#2e7d32" ,"#d0c0a5","#d7660d")) +
  #scale_y_discrete(expand = c(0, 0))+
  expand_limits(y=c(0,30000000))+
  theme_test(base_size = 14)+ 
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = "white"),
        legend.position = "NULL") +
  labs(x = "", y = "Number of MGEs")+
  #labs(subtitle = "Method Genizi") +
  #labs(title = "Gobi") +
  theme(axis.title = element_text(size = 14,
                                  vjust = 1)) + theme(axis.text.x = element_text(vjust = 1,hjust = 1)) + coord_flip()

p1

p2 <- p1 + coord_flip() +  scale_y_break(c(100000,2000000), scales = 0.3, space = 0.1,)

p2

# 
MGEs_diversity_basedSequence$AllMGEnumber <- rep("29629712")
MGEs_diversity_basedSequence$AllMGEnumber <- as.numeric(MGEs_diversity_basedSequence$AllMGEnumber)

MGEs_diversity_basedSequence <- MGEs_diversity_basedSequence %>% mutate(Freq=MGEnumber/AllMGEnumber*100)


MGEs_diversity_percentage <- ggdonutchart(MGEs_diversity_basedSequence[,c(1,2)], "MGEnumber",
                                          label = "MGEs",                               
                                          fill = "MGEs",                            
                                          color = "white",   size = 0,                           
                                          palette = c("#1976d2","#e0e0e0","#dce775", "#f1e54c",  "#2e7d32" ,"#d0c0a5","#d7660d"),
                                          ggtheme = theme_pubr()) + # #bbded6
  theme(axis.text.x = element_blank(), legend.position = "NULL")

MGEs_diversity_percentage

# 2.classifed MGEs -----------------------------------
# load
MGE_ratio_information <- fread("diamond_ratio.csv")
MGE_ratio_information <- MGE_ratio_information %>% select(MGE, Database_Known_percent, Database_Unknown_percent)
MGE_ratio_information <- melt(MGE_ratio_information, variable.name = "Category")

# plot
MGE_ratio_information$MGE <- factor(MGE_ratio_information$MGE, levels = c("Plasmid","Phage","ICE", "Transposon", "IS","Integron"))

MGE_ratio_1 <- ggplot(MGE_ratio_information, aes(x = MGE, y = value, fill=Category))  + geom_col(width = 0.6) +  
  scale_fill_manual(values = c("#9e9d24", "#e0e0e0" )) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
  theme_pubr(border = T, legend = "top", base_size = 14) + 
  scale_y_continuous(expand = c(0.0001, 0.015), limits = c(0,1)) +
  scale_x_discrete(expand = c(0.12, 0.02)) +
  xlab("") + ylab("") + #+ guides(fill = guide_legend(nrow=2)) 
  theme(legend.position="right") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  coord_cartesian(ylim = c(0, 0.1))

MGE_ratio_1

MGE_ratio_2 <- ggplot(MGE_ratio_information, aes(x = MGE, y = value, fill=Category))  + geom_col(width = 0.6) +  
  scale_fill_manual(values = c("#9e9d24", "#e0e0e0" )) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
  theme_pubr(border = T, legend = "top", base_size = 14) + 
  scale_y_continuous(expand = c(0.0001, 0.015), limits = c(0,1)) +
  scale_x_discrete(expand = c(0.12, 0.02)) +
  xlab("") + ylab("Composition of MGEs") + #+ guides(fill = guide_legend(nrow=2)) 
  theme(legend.position="right") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(), #去掉X轴坐标
        axis.ticks.x = element_blank()) + 
  theme(panel.border = element_rect(fill=NA,color="black", size=0.3, linetype="solid")) +
  
  coord_cartesian(ylim = c(0.84, 1))

MGE_ratio_2

MGE_ratio <- ggarrange(MGE_ratio_2,MGE_ratio_1,heights=c(4/7, 3/7),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v") 

MGE_ratio

# 3. MGE abundance  --------------------------------
load("1_TotalMGE_TotalMGEtypes_basedContigsUnique.RData")

totalMGE_Depth_df$sampleType <- factor(totalMGE_Depth_df$sampleType, levels= c("tailings", "Sludge", "farmland","forest","grass","gobi"))

p1<- ggboxplot(data = totalMGE_Depth_df, 
               
               x="sampleType", y="total.MGE.DepthPG", fill = "sampleType", alpha=0.7, width = 0.6, size = 0.5, 
               
               color="lightslategray",outlier.shape = NA) +
  
  geom_jitter(aes(x=sampleType, y=total.MGE.DepthPG, color=sampleType), width = 0.2, size=1.5, alpha=0.5) + 
  
  scale_fill_manual(values = c('#EA5C15','#7986cb',"#addd8e" , "#7fcdbb","#2c7fb8","#bdbdbd")) +
  
  scale_color_manual(values = c('#EA5C15','#7986cb',"#addd8e" , "#7fcdbb","#2c7fb8","#bdbdbd")) +
  
  ylab("Total MGE abundance (coverage, /Gb)") + xlab("") + 
  
  theme_test(base_size = 14) + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  
  ylim(0, 150000) +
  
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) + 
  
  scale_x_discrete(labels=c("Tailings", "Sewage","Farmland","Forest", "Grassland","Desert")) +
  
  theme(legend.position = "none", panel.grid = element_blank()) 

p1

# 4. MGE composition -------------------------------

# load
load("/data/wenp/R/Plasmid/1_All_MGE_mapping_corrected_by_JXforest_1kb.RData")

# 
allMGE_df.l <- All_MGE_mapping_final_df_1kb
#allMGE_df.l <- allMGE_df.l %>% filter(MGE_types != "Others") 

# 
allMGE_df.l <- allMGE_df.l %>% 
  select(Contig, DepthPG, sampleName, MGE_types, sampleType) %>% 
  unique() 

allMGEtype_df.1 <- allMGE_df.l %>% group_by(sampleName, MGE_types) %>%  # 
  summarise(sampleType=unique(sampleType),DepthPG=sum(DepthPG))

unique(allMGEtype_df.1$MGE_types)

# 
total.MGEDepth_composition_df <- allMGE_df.l %>% 
  group_by(sampleName) %>% summarise(total.MGE.DepthPG=sum(DepthPG),numMGE=sum(DepthPG >0), sampleType=unique(sampleType))

allMGEtype_df.1$totalMGE_depth <- sapply(allMGEtype_df.1$sampleName, 
                                         function(x) total.MGEDepth_composition_df$total.MGE.DepthPG[which(total.MGEDepth_composition_df$sampleName == x)])

allMGEtype_df.1$percent <- allMGEtype_df.1$DepthPG/allMGEtype_df.1$totalMGE_depth
sum(allMGEtype_df.1$percent)

#
allMGEtype_df.2 <- allMGEtype_df.1 %>% group_by(sampleType, MGE_types) %>% 
  summarise(sum=sum(percent)) 
unique(allMGEtype_df.2$sampleType)

allMGEtype_df.2$num <- rep(c(30, 87, 81,12,27, 115), each = 7)

allMGEtype_df.2 <- allMGEtype_df.2 %>% mutate(percent = sum/num)

allMGEtype_df.2 <- allMGEtype_df.2 %>% group_by(sampleType, MGE_types) %>% 
  summarise(percent=mean(percent))

# color
MGEType_rankABC <- unique(allMGEtype_df.2$MGE_types)[order(unique(allMGEtype_df.2$MGE_types))]
MGEType_fctLevel <- c(MGEType_rankABC[MGEType_rankABC != "Others"], "Others") 

# plot
allMGEtype_df.2$MGE_type_fct <- factor(allMGEtype_df.2$MGE_types, levels = MGEType_fctLevel)
allMGEtype_df.2$sampleType <- factor(allMGEtype_df.2$sampleType, levels = c("gobi","grass","forest", "farmland", "Sludge","tailings"))

p3 <- ggplot(allMGEtype_df.2, aes(x = sampleType, y = percent, fill=MGE_type_fct))  + geom_col(width = 0.6) +  
  scale_fill_manual(values = c("#1976d2","#dce775","#d0c0a5","#f0e442",  "#d55e00", "#2e7d32", "#e0e0e0" )) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
  theme_pubr(border = T, legend = "top", base_size = 14) + 
  scale_y_continuous(expand = c(0.0001, 0.015), limits = c(0,1)) +
  scale_x_discrete(expand = c(0.12, 0.02)) +
  xlab("") + ylab("Composition of MGE types") + #+ guides(fill = guide_legend(nrow=2)) 
  theme(legend.position="right") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + coord_flip() +
  scale_x_discrete(labels=c("Desert", "Grassland","Forest","Farmland","Sewage","Tailings")) 

p3 +labs(fill="MGE Types")

# 5. MGE abundance per sample site --------------------------------------------

# load
load("1_TotalMGE_TotalMGEtypes_basedContigsUnique.RData")
Farmland_location <- fread("Farmland_location.txt")
Forest_location <- fread("Forest_location.txt")
Grass_location <- fread("Grass_location.txt")
Gobi_location <- fread("Gobi_location.txt")
Sewage_location <- fread("Sewage_location.txt")
Tailings_location <- fread("Tailings_location.txt")

# color for ecosystem ----------------
spType_color_df <- cbind.data.frame(spType=c("tailings","Sludge","farmland", "forest","grass","gobi"),
                                    color=c( '#EA5C15','#7986cb',"#addd8e" , "#7fcdbb","#2c7fb8","#bdbdbd"), 
                                    stringsAsFactors=F)

colnames(totalMGE_Depth_df)[which(colnames(totalMGE_Depth_df) == "total.MGE.DepthPG")] <- "Abundance" 

# (1) plot for Tailings data ----
dat_Tailings <- totalMGE_Depth_df %>% 
  dplyr::filter(grepl("tailings",sampleType,perl = T)) %>%
  as.data.frame()

dat_Tailings$location <- sapply(dat_Tailings$sampleName, function(x) Tailings_location$location[which(Tailings_location$sampleID == x)])
dat_Tailings$Habitat <- rep("tailings", nrow(dat_Tailings))

dat_Tailings_bysite <- dat_Tailings %>% dplyr::group_by(location, Habitat) %>%
  dplyr::summarise(MGE_Abund_avg=mean(Abundance),
                   numMGE_avg=mean(numMGE),
                   MGE_Abund_sd=sd(Abundance),
                   numMGE_sd=sd(numMGE)) %>% as.data.frame() %>% 
  base::merge(All_location_df %>% select(-Habitat), by=c("location")) %>% 
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude)) 


dat_Tailings_bysite$location <- factor(dat_Tailings_bysite$location, levels=dat_Tailings_bysite$location)

colors <- sapply(unique(dat_Tailings_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)]) #两个表格匹配颜色

# 
p_Tailings <- ggplot(dat_Tailings_bysite) +
  geom_errorbar(aes(ymin=MGE_Abund_avg-MGE_Abund_sd, ymax=MGE_Abund_avg+MGE_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=MGE_Abund_avg), width=0.6, color="black", fill=colors, alpha = 0.9) + #geom_col与geom_bar一样也是柱状图，两者差别是，它能映射到Y值
  #sacle_y_continuous(limits = c(4000, 50000)) +
  
  theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,150000)) +
  theme(axis.ticks = element_line(colour = "gray5"), axis.text = element_text(colour = "gray5")) + 
  #scale_x_discrete(labels= sub("_Tail","",unique(dat_Tailings_bysite$location))) + ylab("Total MGE abundance (coverage, x/Gb") + xlab("")
  scale_x_discrete(labels= substr(unique(dat_Tailings_bysite$location), 1, 7)) + 
  ylab("Total MGE abundance (coverage, x/Gb") + xlab("") # 要切割下字符串，不然组图的时候会图形展示很不对称

p_Tailings

# (2) plot for Sewage data ----
dat_Sewage <- totalMGE_Depth_df %>% 
  dplyr::filter(grepl("Sludge",sampleType,perl = T)) %>%
  as.data.frame()

dat_Sewage$location <- sapply(dat_Sewage$sampleName, function(x) Sewage_location$location[which(Sewage_location$sampleID == x)])
dat_Sewage$Habitat <- rep("Sludge", nrow(dat_Sewage))

dat_Sewage_bysite <- dat_Sewage %>% dplyr::group_by(location, Habitat) %>%
  dplyr::summarise(MGE_Abund_avg=mean(Abundance),
                   numMGE_avg=mean(numMGE),
                   MGE_Abund_sd=sd(Abundance),
                   numMGE_sd=sd(numMGE)) %>% as.data.frame() 

dat_Sewage_bysite$location <- unlist(dat_Sewage_bysite$location)

dat_Sewage_bysite <- merge(dat_Sewage_bysite, All_location_df %>% select(-Habitat), 
                           by.x =c("location"), by.y = c("location_raw")) %>%
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude)) 

dat_Sewage_bysite$location.y <- factor(dat_Sewage_bysite$location.y, levels=dat_Sewage_bysite$location.y)

colors <- sapply(unique(dat_Sewage_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)]) #两个表格匹配颜色

# 
p_Sewage <- ggplot(dat_Sewage_bysite) +
  geom_errorbar(aes(ymin=MGE_Abund_avg-MGE_Abund_sd, ymax=MGE_Abund_avg+MGE_Abund_sd, x=location.y), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location.y, y=MGE_Abund_avg), width=0.6, color="black", fill=colors, alpha = 0.9) + #
  theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,150000)) +
  theme(axis.ticks = element_line(colour = "gray5"), axis.text = element_text(colour = "gray5")) + 
  scale_x_discrete(labels= substr(unique(dat_Sewage_bysite$location.y), 1, 8)) + 
  ylab("Total MGE abundance (coverage, x/Gb") + xlab("") # 要切割下字符串，不然组图的时候会图形展示很不对称

p_Sewage

# (3) plot for Forest data ----
dat_Forest <- totalMGE_Depth_df %>% 
  dplyr::filter(grepl("forest",sampleType,perl = T)) %>%
  as.data.frame()

dat_Forest$location <- sapply(dat_Forest$sampleName, function(x) Forest_location$location[which(Forest_location$sampleTyple == x)])
dat_Forest$Habitat <- rep("forest", nrow(dat_Forest))

dat_Forest_bysite <- dat_Forest %>% dplyr::group_by(location, Habitat) %>%
  dplyr::summarise(MGE_Abund_avg=mean(Abundance),
                   numMGE_avg=mean(numMGE),
                   MGE_Abund_sd=sd(Abundance),
                   numMGE_sd=sd(numMGE)) %>% as.data.frame() %>% 
  base::merge(all_location %>% select(-Habitat), by=c("location")) %>% 
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude)) 

dat_Forest_bysite$location <- factor(dat_Forest_bysite$location, levels=dat_Forest_bysite$location)

colors <- sapply(unique(dat_Forest_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)]) 

# 
p_Forest <- ggplot(dat_Forest_bysite) +
  geom_errorbar(aes(ymin=MGE_Abund_avg-MGE_Abund_sd, ymax=MGE_Abund_avg+MGE_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=MGE_Abund_avg), width=0.6, color="black", fill=colors, alpha = 0.9) + 
  theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,150000)) +
  theme(axis.ticks = element_line(colour = "gray5"), axis.text = element_text(colour = "gray5")) + 
  scale_x_discrete(labels= sub("_Forest","",unique(dat_Forest_bysite$location))) + ylab("Total MGE abundance (coverage, x/Gb") + xlab("")

p_Forest

# (4) plot for Farmland data ----
dat_Farmland <- totalMGE_Depth_df %>% 
  dplyr::filter(grepl("farmland",sampleType,perl = T)) %>%
  as.data.frame()

dat_Farmland$location <- sapply(dat_Farmland$sampleName, function(x) Farmland_location$location[which(Farmland_location$sampleTyple == x)])
dat_Farmland$Habitat <- rep("farmland", nrow(dat_Farmland))

dat_Farmland_bysite <- dat_Farmland %>% dplyr::group_by(location, Habitat) %>%
  dplyr::summarise(MGE_Abund_avg=mean(Abundance),
                   numMGE_avg=mean(numMGE),
                   MGE_Abund_sd=sd(Abundance),
                   numMGE_sd=sd(numMGE)) %>% as.data.frame() %>% 
  base::merge(all_location %>% select(-Habitat), by=c("location")) %>% 
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude)) 

dat_Farmland_bysite$location <- factor(dat_Farmland_bysite$location, levels=dat_Farmland_bysite$location)

colors <- sapply(unique(dat_Farmland_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)]) 

# 
p_Farmland <- ggplot(dat_Farmland_bysite) +
  geom_errorbar(aes(ymin=MGE_Abund_avg-MGE_Abund_sd, ymax=MGE_Abund_avg+MGE_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=MGE_Abund_avg), width=0.6, color="black", fill=colors,alpha = 0.9) +
  theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,150000)) +
  theme(axis.ticks = element_line(colour = "gray5"), axis.text = element_text(colour = "gray5")) + 
  scale_x_discrete(labels= sub("_Farmland","",unique(dat_Farmland_bysite$location))) + ylab("Total MGE abundance (coverage, x/Gb") + xlab("")

p_Farmland

# (5) plot for Gobi data ----
dat_Gobi <- totalMGE_Depth_df %>% 
  dplyr::filter(grepl("gobi",sampleType,perl = T)) %>%
  as.data.frame()

dat_Gobi$location <- sapply(dat_Gobi$sampleName, function(x) Gobi_location$location[which(Gobi_location$sampleTyple == x)])
dat_Gobi$Habitat <- rep("gobi", nrow(dat_Gobi))

dat_Gobi_bysite <- dat_Gobi %>% dplyr::group_by(location, Habitat) %>%
  dplyr::summarise(MGE_Abund_avg=mean(Abundance),
                   numMGE_avg=mean(numMGE),
                   MGE_Abund_sd=sd(Abundance),
                   numMGE_sd=sd(numMGE)) %>% as.data.frame() %>% 
  base::merge(all_location %>% select(-Habitat), by=c("location")) %>% 
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude)) 


dat_Gobi_bysite$location <- factor(dat_Gobi_bysite$location, levels=dat_Gobi_bysite$location)

colors <- sapply(unique(dat_Gobi_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)]) 

#
p_Gobi<- ggplot(dat_Gobi_bysite) +
  geom_errorbar(aes(ymin=MGE_Abund_avg-MGE_Abund_sd, ymax=MGE_Abund_avg+MGE_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=MGE_Abund_avg), width=0.6, color="black", fill=colors,alpha = 0.9) 
theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,150000)) +
  theme(axis.ticks = element_line(colour = "gray5"), axis.text = element_text(colour = "gray5")) + 
  scale_x_discrete(labels= sub("_Gobi","",unique(dat_Gobi_bysite$location))) + ylab("Total MGE abundance (coverage, x/Gb") + xlab("")

p_Gobi

# (6) plot for Grass data ----
dat_Grass <- totalMGE_Depth_df %>% 
  dplyr::filter(grepl("grass",sampleType,perl = T)) %>%
  as.data.frame()

dat_Grass$location <- sapply(dat_Grass$sampleName, function(x) Grass_location$location[which(Grass_location$sampleTyple == x)])
dat_Grass$Habitat <- rep("grass", nrow(dat_Grass))

dat_Grass_bysite <- dat_Grass %>% dplyr::group_by(location, Habitat) %>%
  dplyr::summarise(MGE_Abund_avg=mean(Abundance),
                   numMGE_avg=mean(numMGE),
                   MGE_Abund_sd=sd(Abundance),
                   numMGE_sd=sd(numMGE)) %>% as.data.frame() %>% 
  base::merge(all_location %>% select(-Habitat), by=c("location")) %>% 
  dplyr::arrange(latitude) %>% dplyr::filter(!is.na(latitude)) 


dat_Grass_bysite$location <- factor(dat_Grass_bysite$location, levels=dat_Grass_bysite$location)

colors <- sapply(unique(dat_Grass_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)]) 

# 
p_Grass <- ggplot(dat_Grass_bysite) +
  geom_errorbar(aes(ymin=MGE_Abund_avg-MGE_Abund_sd, ymax=MGE_Abund_avg+MGE_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=MGE_Abund_avg), width=0.6, color="black", fill=colors, alpha = 0.9) + 
  theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,150000)) +
  theme(axis.ticks = element_line(colour = "gray5"), axis.text = element_text(colour = "gray5")) + 
  scale_x_discrete(labels= sub("_Grass","",unique(dat_Grass_bysite$location))) + ylab("Total MGE abundance (coverage, x/Gb") + xlab("")

p_Grass

layout <- "
AAAAAAAABBBBBBBCCCCCCCDDDDDDDEEEFF
"
p_Tailings + p_Sewage + p_Farmland + p_Forest + p_Grass + p_Gobi + plot_layout(design = layout)


# 6. MGE composition per sample site -------------------------------------------
load("2_Total_MGE_composition_bysampleType.RData")
All_location_df <- fread("All_location_df.txt")

# tailings -----------------------------
dat_Tailings_df <- dat_Tailings_df %>% group_by(sampleName, MGE_types) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), percent=sum(percent))

plotDat <- dat_Tailings_df %>% group_by(location, MGE_types) %>% summarise(percent=mean(percent)) %>% as.data.frame() 
sum(plotDat$percent)

plotDat$latitute <- sapply(plotDat$location, 
                           function(x) All_location_df$latitude[which(All_location_df$location == x)]) 

# color
MGEType_rankABC <- unique(plotDat$MGE_types)[order(unique(plotDat$MGE_types))]
MGEType_fctLevel <- c(MGEType_rankABC[MGEType_rankABC != "Others"], "Others") 

MGEType_fctLevel <- c("Others", "ICE", "Integron", "IS", "Phage", "Plasmids", "Transposon")

plotDat_Tailings <- plotDat %>%  arrange(latitute)
plotDat_Tailings$location <- factor(plotDat_Tailings$location,levels = unique(plotDat_Tailings$location))
plotDat_Tailings$MGE_type_fct <- factor(plotDat_Tailings$MGE_types,
                                        
                                        levels = MGEType_fctLevel)

p_tailings_df <- ggplot(plotDat_Tailings, aes(x = location, y = percent, fill=MGE_type_fct))  + geom_col() +  
  
  # scale_fill_manual(values = c("#b2dfdb","#f0e442",   "#aed581" ,"#ffe0b2","#d55e00", "#e0e0e0")) +
  
  scale_fill_manual(values = c("#e0e0e0","#1976d2","#dce775","#d0c0a5","#f0e442",  "#d55e00", "#2e7d32")) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("Composition of MGE types") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Tail","",unique(plotDat_Tailings$location))) 

p_tailings_df


# Sewage -----------------------------

dat_Sludge_df <- dat_Sludge_df %>% group_by(sampleName, MGE_types) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), percent=sum(percent))

plotDat <- dat_Sludge_df %>% group_by(location, MGE_types) %>% summarise(percent=mean(percent)) %>% as.data.frame() 
sum(plotDat$percent)

str(plotDat)

plotDat$location <- as.vector(plotDat$location) 

plotDat$latitute <- sapply(plotDat$location, 
                           function(x) All_location_df$latitude[which(All_location_df$location == x)]) 

plotDat$location <- sapply(strsplit(plotDat$location, "_", fixed = T), "[[", 2) 
# color
MGEType_rankABC <- unique(plotDat$MGE_types)[order(unique(plotDat$MGE_types))]
MGEType_fctLevel <- c(MGEType_rankABC[MGEType_rankABC != "Others"], "Others") 

MGEType_fctLevel <- c("Others", "ICE", "Integron", "IS", "Phage", "Plasmids", "Transposon")

plotDat_Sewage <- plotDat %>%  arrange(latitute)
plotDat_Sewage$location <- factor(plotDat_Sewage$location,levels = unique(plotDat_Sewage$location))
plotDat_Sewage$MGE_type_fct <- factor(plotDat_Sewage$MGE_types,
                                      
                                      levels = MGEType_fctLevel)

p_sewage_df <- ggplot(plotDat_Sewage, aes(x = location, y = percent, fill=MGE_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = c("#e0e0e0","#1976d2","#dce775","#d0c0a5","#f0e442",  "#d55e00", "#2e7d32")) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("Composition of MGE types") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) 

p_sewage_df

# Farmland -----------------------------

dat_Farmland_df <- dat_Farmland_df %>% group_by(sampleName, MGE_types) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), percent=sum(percent))

plotDat <- dat_Farmland_df %>% group_by(location, MGE_types) %>% summarise(percent=mean(percent)) %>% as.data.frame() 
sum(plotDat$percent)

#
plotDat$latitute <- sapply(plotDat$location, 
                           function(x) All_location_df$latitude[which(All_location_df$location == x)]) 

# color
MGEType_rankABC <- unique(plotDat$MGE_types)[order(unique(plotDat$MGE_types))]
MGEType_fctLevel <- c(MGEType_rankABC[MGEType_rankABC != "Others"], "Others") 

MGEType_fctLevel <- c("Others", "ICE", "Integron", "IS", "Phage", "Plasmids", "Transposon")

plotDat_Farmland <- plotDat %>%  arrange(latitute)
plotDat_Farmland$location <- factor(plotDat_Farmland$location,levels = unique(plotDat_Farmland$location))
plotDat_Farmland$MGE_type_fct <- factor(plotDat_Farmland$MGE_types,
                                        
                                        levels = MGEType_fctLevel)

p_farmland_df <- ggplot(plotDat_Farmland, aes(x = location, y = percent, fill=MGE_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = c("#e0e0e0","#1976d2","#dce775","#d0c0a5","#f0e442",  "#d55e00", "#2e7d32")) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("Composition of MGE types") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Farmland","",unique(plotDat_Farmland$location))) 

p_farmland_df

# Forest -----------------------------

dat_Forest_df <- dat_Forest_df %>% group_by(sampleName, MGE_types) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), percent=sum(percent))

plotDat <- dat_Forest_df %>% group_by(location, MGE_types) %>% summarise(percent=mean(percent)) %>% as.data.frame() 
sum(plotDat$percent)

# 
plotDat$latitute <- sapply(plotDat$location, 
                           function(x) All_location_df$latitude[which(All_location_df$location == x)]) 

# color
MGEType_rankABC <- unique(plotDat$MGE_types)[order(unique(plotDat$MGE_types))]
MGEType_fctLevel <- c(MGEType_rankABC[MGEType_rankABC != "Others"], "Others") 

MGEType_fctLevel <- c("Others", "ICE", "Integron", "IS", "Phage", "Plasmids", "Transposon")

plotDat_Forest <- plotDat %>%  arrange(latitute)
plotDat_Forest$location <- factor(plotDat_Forest$location,levels = unique(plotDat_Forest$location))
plotDat_Forest$MGE_type_fct <- factor(plotDat_Forest$MGE_types,
                                      
                                      levels = MGEType_fctLevel)

p_forest_df <- ggplot(plotDat_Forest, aes(x = location, y = percent, fill=MGE_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = c("#e0e0e0","#1976d2","#dce775","#d0c0a5","#f0e442",  "#d55e00", "#2e7d32")) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("Composition of MGE types") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Forest","",unique(plotDat_Forest$location))) 

p_forest_df

# Grass -----------------------------

dat_Grass_df <- dat_Grass_df %>% group_by(sampleName, MGE_types) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), percent=sum(percent))

plotDat <- dat_Grass_df %>% group_by(location, MGE_types) %>% summarise(percent=mean(percent)) %>% as.data.frame() 
sum(plotDat$percent)

# 
plotDat$latitute <- sapply(plotDat$location, 
                           function(x) All_location_df$latitude[which(All_location_df$location == x)]) 

# color
MGEType_rankABC <- unique(plotDat$MGE_types)[order(unique(plotDat$MGE_types))]
MGEType_fctLevel <- c(MGEType_rankABC[MGEType_rankABC != "Others"], "Others") 

MGEType_fctLevel <- c("Others", "ICE", "Integron", "IS", "Phage", "Plasmids", "Transposon")

plotDat_Grass <- plotDat %>%  arrange(latitute)
plotDat_Grass$location <- factor(plotDat_Grass$location,levels = unique(plotDat_Grass$location))
plotDat_Grass$MGE_type_fct <- factor(plotDat_Grass$MGE_types,
                                     
                                     levels = MGEType_fctLevel)

p_grass_df <- ggplot(plotDat_Grass, aes(x = location, y = percent, fill=MGE_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = c("#e0e0e0","#1976d2","#dce775","#d0c0a5","#f0e442",  "#d55e00", "#2e7d32")) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("Composition of MGE types") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Grass","",unique(plotDat_Grass$location))) 

p_grass_df

# gobi -----------------------------

dat_Gobi_df <- dat_Gobi_df %>% group_by(sampleName, MGE_types) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), percent=sum(percent))

plotDat <- dat_Gobi_df %>% group_by(location, MGE_types) %>% summarise(percent=mean(percent)) %>% as.data.frame() 
sum(plotDat$percent)

# 
plotDat$latitute <- sapply(plotDat$location, 
                           function(x) All_location_df$latitude[which(All_location_df$location == x)]) 

# color
MGEType_rankABC <- unique(plotDat$MGE_types)[order(unique(plotDat$MGE_types))]
MGEType_fctLevel <- c(MGEType_rankABC[MGEType_rankABC != "Others"], "Others") 

MGEType_fctLevel <- c("Others", "ICE", "Integron", "IS", "Phage", "Plasmids", "Transposon")

plotDat_Gobi <- plotDat %>%  arrange(latitute)
plotDat_Gobi$location <- factor(plotDat_Gobi$location,levels = unique(plotDat_Gobi$location))
plotDat_Gobi$MGE_type_fct <- factor(plotDat_Gobi$MGE_types,
                                    
                                    levels = MGEType_fctLevel)

p_gobi_df <- ggplot(plotDat_Gobi, aes(x = location , y = percent, fill=MGE_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = c("#e0e0e0","#1976d2","#dce775","#d0c0a5","#f0e442",  "#d55e00", "#2e7d32")) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  xlab("") + ylab("Composition of MGE types") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Gobi","",unique(plotDat_Gobi$location))) 

p_gobi_df

layout <- "
AAAAAAAABBBBBBBCCCCCCCDDDDDDDEEF
"
p_tailings_df + p_sewage_df + p_farmland_df + p_forest_df + p_grass_df + p_gobi_df + plot_layout(design = layout)
