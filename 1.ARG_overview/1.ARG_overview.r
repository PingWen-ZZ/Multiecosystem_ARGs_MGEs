
# library ------------------
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot) 
library(ggpubr) 
library(patchwork)
library(eulerr)
library(data.table) 
library(Hmisc)
library(scales)
library(ggsci)
library(reshape2)
library(eulerr)

# 1. Total abundance diversity and composition ------------------------------------

# load
total.ARGDepth_df <- fread("total.abundance.diversity.txt")

# plot -----
# ARG abundance  -----------
total.ARGDepth_df$sampleType <- factor(total.ARGDepth_df$sampleType, levels= c("T", "S", "Fa","Fo","Gr","Go"))

p1<- ggboxplot(data = total.ARGDepth_df, 
               
               x="sampleType", y="total.ARG.DepthPG", fill = "sampleType", alpha=0.7, width = 0.6, size = 0.5, 
               
               color="lightslategray",outlier.shape = NA) +
  
  geom_jitter(aes(x=sampleType, y=total.ARG.DepthPG, color=sampleType), width = 0.2, size=1.5, alpha=0.5) + 
  
  scale_fill_manual(values = c('#EA5C15','#7986cb',"#addd8e" , "#7fcdbb","#2c7fb8","#bdbdbd")) +
  
  scale_color_manual(values = c('#EA5C15','#7986cb',"#addd8e" , "#7fcdbb","#2c7fb8","#bdbdbd")) + 
  
  theme_test() + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  
  ylab("Total ARG abundance (coverage /Gb)") + xlab("") + 
  
  ylim(50, 4100) +
  
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) + 
  
  scale_x_discrete(labels=c("Tailings", "Sewage","Farmland","Forest", "Grassland","Gobi desert")) +
  
  theme(legend.position = "none", panel.grid = element_blank())

p1

# ARG number ---------
p2<- ggboxplot(data = total.ARGDepth_df, 
               
               x="sampleType", y="numARG", fill = "sampleType", alpha=0.6, width = 0.6, size = 0.5, 
               
               color="lightslategray",outlier.shape = NA) +
  
  geom_jitter(aes(x=sampleType, y=numARG, color=sampleType), width = 0.2, size=1.5, alpha=0.5) + 
  
  scale_fill_manual(values = c('#EA5C15','#7986cb',"#addd8e" , "#7fcdbb","#2c7fb8","#bdbdbd")) +
  
  scale_color_manual(values = c('#EA5C15','#7986cb',"#addd8e" , "#7fcdbb","#2c7fb8","#bdbdbd")) +
  
  theme_test() + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  
  
  ylab("Number of ARG subtypes") + xlab("") + 
  
  ylim(50, 550) +
  
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) + 
  
  scale_x_discrete(labels=c("Tailings", "Sewage","Farmland","Forest", "Grassland","Gobi desert")) +
  
  theme(legend.position = "none", panel.grid = element_blank())

p2

# ARG composition --------------------------------------------------

drug_depth_df <- fread("drug.composition.txt")

drugType_colorsDf <- cbind.data.frame(
  
  Colors=c(pal_npg("nrc")(10), "gray","gold3","orchid"),
  
  DrugTypes=c("aminoglycoside", "bacitracin", "beta-lactam","glycopeptide",
              
              "MLS","tetracycline","phenicol","multidrug",
              
              "rifamycin","sulfonamide","Others","mupirocin","fosmidomycin"), 
  
  stringsAsFactors=F
  
)

drugType_rankABC <- unique(drug_depth_df$drug_type)[order(unique(drug_depth_df$drug_type))]
drugType_fctLevel <- c(drugType_rankABC[drugType_rankABC != "Others"], "Others") #将others放到末尾

drugType_pallette <-sapply(drugType_fctLevel, 
                           
                           function(x) drugType_colorsDf$Colors[which(drugType_colorsDf$DrugTypes == x)])

# plot
drug_depth_df$drug_type_fct <- factor(drug_depth_df$drug_type, levels = drugType_fctLevel)
drug_depth_df$sampleType <- factor(drug_depth_df$sampleType, levels = c("Gr", "Fo","Fa","Go", "S", "T"))

p3 <- ggplot(drug_depth_df, aes(x = sampleType, y = drug_percent, fill=drug_type_fct))  + geom_col(width = 0.6) +  
  scale_fill_manual(values = c("#9FCD99", "#FFDD71","#DFC27D", "#ce3d36", "#e9a91f","#0b8e36", "#D9F0D3", "#D1E5F0", "#FFB977","#B2AD8F","gray")) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + 
  theme_test(base_size = 14) +
  scale_y_continuous(expand = c(0.0001, 0.015), limits = c(0,1)) +
  scale_x_discrete(expand = c(0.12, 0.02)) +
  xlab("") + ylab("Composition of ARG types") + #+ guides(fill = guide_legend(nrow=2)) 
  theme(legend.position="right") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + coord_flip()

p3 + labs(fill="ARG Types")

layout <- "
AAAABBBBCCCCCCCCC
"
p1 + p2 + p3 + plot_layout(design = layout)

# 2. Abundance and diversity per sample sites ----------------------------
# load 
dat_Farmland_bysite <- fread("dat_Farmland_bysite.txt")
dat_Forest_bysite <- fread("dat_Forest_bysite.txt")
dat_Grass_bysite <- fread("dat_Grass_bysite.txt")
dat_Gobi_bysite <- fread("dat_Gobi_bysite.txt")
dat_Sewage_bysite <- fread("dat_Sewage_bysite.txt")

# color for ecosystem 
spType_color_df <- cbind.data.frame(spType=c("S","Fa","Fo","Gr","Go"),
                                    color=c('#7986cb', "#addd8e" , "#7fcdbb","#2c7fb8","#deebf7"), 
                                    stringsAsFactors=F)

# (1) plot for Farmland data ----
dat_Farmland_bysite$location <- factor(dat_Farmland_bysite$location, levels=dat_Farmland_bysite$location)

colors <- sapply(unique(dat_Farmland_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)]) 

# 
p1 <- ggplot(dat_Farmland_bysite) +
  geom_errorbar(aes(ymin=ARG_Abund_avg-ARG_Abund_sd, ymax=ARG_Abund_avg+ARG_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=ARG_Abund_avg), width=0.6, color="black", fill=colors) + #geom_col与geom_bar一样也是柱状图，两者差别是，它能映射到Y值
  #sacle_y_continuous(limits = c(0,6000)) +
  
  theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,3500)) +
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) + 
  scale_x_discrete(labels= sub("_Farmland","",unique(dat_Farmland_bysite$location))) + ylab("") + xlab("")

p1

# 
p2 <- ggplot(dat_Farmland_bysite) +
  geom_point(aes(x=location, y=numsubtype_avg), size=3, shape=21, fill="#FFF6F5") +
  scale_y_continuous(limits = c(200,600), position = "right") + 
  #scale_y_continuous(position = "right") +
  theme_half_open(font_size=9, rel_small = 1) + #theme_half_open这个主题只有一半的坐标轴显示，不是四周全封闭的,属于cowplot包里面的
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") + ylab("")
p2
# 
#aligned_plots <- align_plots(p1, p2, align = "hv", axis = "tblr") 
aligned_plots <- align_patches(p1, p2) 
p_Farmland <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p_Farmland

# (2) plot for Forest data -----------
dat_Forest_bysite$location <- factor(dat_Forest_bysite$location, levels=dat_Forest_bysite$location)

colors <- sapply(unique(dat_Forest_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)])

# paint
p1 <- ggplot(dat_Forest_bysite) +
  geom_errorbar(aes(ymin=ARG_Abund_avg-ARG_Abund_sd, ymax=ARG_Abund_avg+ARG_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=ARG_Abund_avg), width=0.6, color="black", fill=colors) + 
  #sacle_y_continuous(limits = c(0,6000)) +
  
  theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,3500)) +
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) + 
  scale_x_discrete(labels= sub("_Forest","",unique(dat_Forest_bysite$location))) + ylab("") + xlab("")
p1

p2 <- ggplot(dat_Forest_bysite) +
  geom_point(aes(x=location, y=numsubtype_avg), size=3, shape=21, fill="#F3FBEF") +
  scale_y_continuous(limits = c(200,600), position = "right") + 
  #scale_y_continuous(position = "right") +
  theme_half_open(font_size=9, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") + ylab("")
p2
# 
aligned_plots <- align_patches(p1, p2) 
p_Forest <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p_Forest

# (3) plot for Grass data ----------------
dat_Grass_bysite$location <- factor(dat_Grass_bysite$location, levels=dat_Grass_bysite$location)

colors <- sapply(unique(dat_Grass_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)])

# paint
p1 <- ggplot(dat_Grass_bysite) +
  geom_errorbar(aes(ymin=ARG_Abund_avg-ARG_Abund_sd, ymax=ARG_Abund_avg+ARG_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=ARG_Abund_avg), width=0.6, color="black", fill=colors) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,3500)) +
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) + 
  scale_x_discrete(labels= sub("_Grass","",unique(dat_Grass_bysite$location))) + ylab("") + xlab("")
p1

p2 <- ggplot(dat_Grass_bysite) +
  geom_point(aes(x=location, y=numsubtype_avg), size=3, shape=21, fill="#FFFFED") +
  scale_y_continuous(limits = c(200,600), position = "right") + 
  #scale_y_continuous(position = "right") +
  theme_half_open(font_size=9, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") + ylab("")
p2

aligned_plots <- align_patches(p1, p2) 
p_Grass <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p_Grass

# (4) plot for Gobi data -------------------
dat_Gobi_bysite$location <- factor(dat_Gobi_bysite$location, levels=dat_Gobi_bysite$location)

colors <- sapply(unique(dat_Gobi_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)])

# paint
colnames(allARG_df.l)[which(colnames(allARG_df.l) == "DepthPG")] <- "Abundance" 

p1 <- ggplot(dat_Gobi_bysite) +
  geom_errorbar(aes(ymin=ARG_Abund_avg-ARG_Abund_sd, ymax=ARG_Abund_avg+ARG_Abund_sd, x=location), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location, y=ARG_Abund_avg), width=0.6, color="black", fill=colors) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,3500)) +
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) + 
  scale_x_discrete(labels= sub("_Gobi","",unique(dat_Gobi_bysite$location))) + ylab("") + xlab("")
p1

p2 <- ggplot(dat_Gobi_bysite) +
  geom_point(aes(x=location, y=numsubtype_avg), size=3, shape=21, fill="white") +
  scale_y_continuous(limits = c(200,600), position = "right") + 
  #scale_y_continuous(position = "right") +
  theme_half_open(font_size=9, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") + ylab("No.of ARG subtypes")
p2

aligned_plots <- align_patches(p1, p2) 
p_Gobi <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]]) 
p_Gobi

# (5) plot for Sewage data ----
dat_Sewage_bysite$location.y <- factor(dat_Sewage_bysite$location.y, levels=dat_Sewage_bysite$location.y)

colors <- sapply(unique(dat_Sewage_bysite$Habitat),
                 function(x) spType_color_df$color[which(spType_color_df$spType == x)]) 

# paint

p1 <- ggplot(dat_Sewage_bysite) +
  geom_errorbar(aes(ymin=ARG_Abund_avg-ARG_Abund_sd, ymax=ARG_Abund_avg+ARG_Abund_sd, x=location.y), width=0.2, position=position_dodge(0.9)) + 
  geom_col(aes(x=location.y, y=ARG_Abund_avg), width=0.6, color="black", fill=colors) + 
  theme_test() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0.00006, 50), limits = c(0,3500)) +
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) + 
  ylab("Total ARG abundance (coverage, x/Gb") + xlab("")

p1

p2 <- ggplot(dat_Sewage_bysite) +
  geom_point(aes(x=location.y, y=numsubtype_avg), size=3, shape=21, fill="#FFF6F5") +
  scale_y_continuous(limits = c(200,600), position = "right") + 
  #scale_y_continuous(position = "right") +
  theme_half_open(font_size=9, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") + ylab("")
p2
p_Sewage <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p_Sewage


library(patchwork)
# patchwork自定义布局

layout <- "
AAAAAAABBBBBBBCCCCCCDDDEE
"
p_Sewage + p_Farmland + p_Forest + p_Grass + p_Gobi + plot_layout(design = layout)

# 3. ARG composition per site ---------------------------------------------------

# load
dat_Farmland_df <- fread("dat_Farmland.txt")
dat_Forest_df <- fread("dat_Forest.txt")
dat_Gobi_df <- fread("dat_Gob.txt")
dat_Grass_df <- fread("dat_Grass.txt")
dat_Sewage_df <- fread("dat_Sewage.txt")
all_location <- fread("all_location.txt")

# Farmland ------------------
dat_Farmland_df <- dat_Farmland_df %>% group_by(sampleID, drug_type) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), drug_percent=sum(drug_percent))

plotDat <- dat_Farmland_df %>% group_by(location, drug_type) %>% summarise(drug_percent=mean(drug_percent)) %>% as.data.frame() 

# 
plotDat$latitute <- sapply(plotDat$location, 
                           function(x) all_location$latitude[which(all_location$location == x)]) 
# color
drugType_colorsDf <- cbind.data.frame(
  
  Colors=c("#9FCD99", "#FFDD71","#DFC27D", "#ce3d36", "#e9a91f","#0b8e36", "#D9F0D3", "#D1E5F0", "#FFB977","#B2AD8F","gray"),
  
  DrugTypes=c("aminoglycoside", "bacitracin", "beta-lactam","glycopeptide",
              
              "MLS","multidrug","peptide", "phenicol","sulfonamide","tetracycline", "Others"), 
  
  stringsAsFactors=F
  
)

drugType_rankABC <- unique(plotDat$drug_type)[order(unique(plotDat$drug_type))] 
drugType_fctLevel <- c(drugType_rankABC[drugType_rankABC != "Others"], "Others") #Others放最后

drugType_pallette <-sapply(drugType_fctLevel, 
                           
                           function(x) drugType_colorsDf$Colors[which(drugType_colorsDf$DrugTypes == x)])

# plot ARG composition
plotDat_Farmland <- plotDat %>%  arrange(latitute) 

plotDat_Farmland$location <- factor(plotDat_Farmland$location,levels = unique(plotDat_Farmland$location))

plotDat_Farmland$drug_type_fct <- factor(plotDat_Farmland$drug_type,
                                         
                                         levels = drugType_fctLevel)

p_Farmland <- ggplot(plotDat_Farmland, aes(x = location, y = drug_percent, fill=drug_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = drugType_pallette) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Farmland","",unique(plotDat_Farmland$location))) + 
  
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) 

p_Farmland

# Forest ------------------------
dat_Forest_df <- dat_Forest_df %>% group_by(sampleID, drug_type) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), drug_percent=sum(drug_percent))

plotDat <- dat_Forest_df %>% group_by(location, drug_type) %>% summarise(drug_percent=mean(drug_percent)) %>% as.data.frame() 

plotDat$latitute <- sapply(plotDat$location, 
                           function(x) all_location$latitude[which(all_location$location == x)]) 

# plot ARG composition
plotDat_Forest <- plotDat %>%  arrange(latitute) 

plotDat_Forest$location <- factor(plotDat_Forest$location,levels = unique(plotDat_Forest$location))

plotDat_Forest$drug_type_fct <- factor(plotDat_Forest$drug_type,
                                       
                                       levels = drugType_fctLevel)

p_Forest <- ggplot(plotDat_Forest, aes(x = location, y = drug_percent, fill=drug_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = drugType_pallette) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Forest","",unique(plotDat_Forest$location))) + 
  
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) 

p_Forest

# Gobi ----------------------
dat_Gobi_df <- dat_Gobi_df %>% group_by(sampleID, drug_type) %>% 
  summarise(DepthPGb=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), drug_percent=sum(drug_percent))

plotDat <- dat_Gobi_df %>% group_by(location, drug_type) %>% summarise(drug_percent=mean(drug_percent)) %>% as.data.frame() 

plotDat$latitute <- sapply(plotDat$location, 
                           function(x) all_location$latitude[which(all_location$location == x)]) 

# plot ARG composition
plotDat_Gobi <- plotDat %>%  arrange(latitute) 

plotDat_Gobi$location <- factor(plotDat_Gobi$location,levels = unique(plotDat_Gobi$location))

plotDat_Gobi$drug_type_fct <- factor(plotDat_Gobi$drug_type,
                                     
                                     levels = drugType_fctLevel)

p_Gobi <- ggplot(plotDat_Gobi, aes(x = location, y = drug_percent, fill=drug_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = drugType_pallette) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Gobi","",unique(plotDat_Gobi$location))) + 
  
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) 

p_Gobi

# Sewage ------------------------------------------------

dat_Sewage_df <- dat_Sewage_df %>% group_by(sampleID, drug_type) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), drug_percent=sum(drug_percent)) 
plotDat <- dat_Sewage_df %>% group_by(location, drug_type) %>% summarise(drug_percent=mean(drug_percent)) %>% as.data.frame() 

plotDat$latitute <- sapply(plotDat$location, 
                           function(x) All_location_df$latitude[which(All_location_df$location_raw == x)]) 

plotDat$location_simplify <- sapply(plotDat$location, 
                                    function(x) All_location_df$location[which(All_location_df$location_raw == x)]) 

# plot ARG composition
plotDat_Sewage <- plotDat %>%  arrange(latitute) 

plotDat_Sewage$location_simplify <- factor(plotDat_Sewage$location_simplify,levels = unique(plotDat_Sewage$location_simplify))

plotDat_Sewage$drug_type_fct <- factor(plotDat_Sewage$drug_type,
                                       
                                       levels = drugType_fctLevel)


p_Sewage <- ggplot(plotDat_Sewage, aes(x = location_simplify, y = drug_percent, fill=drug_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = drugType_pallette) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("Composition of ARG types") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) 

p_Sewage


# Grass ------------------------

dat_Grass_df <- dat_Grass_df %>% group_by(sampleID, drug_type) %>% 
  summarise(DepthPG=sum(DepthPG), sampleType=unique(sampleType), location=unique(location), drug_percent=sum(drug_percent)) 
plotDat <- dat_Grass_df %>% group_by(location, drug_type) %>% summarise(drug_percent=mean(drug_percent)) %>% as.data.frame() 

plotDat$latitute <- sapply(plotDat$location, 
                           function(x) all_location$latitude[which(all_location$location == x)]) 

# plot ARG composition
plotDat_Grass <- plotDat %>%  arrange(latitute) 

plotDat_Grass$location <- factor(plotDat_Grass$location,levels = unique(plotDat_Grass$location))

plotDat_Grass$drug_type_fct <- factor(plotDat_Grass$drug_type,
                                      
                                      levels = drugType_fctLevel)


p_Grass <- ggplot(plotDat_Grass, aes(x = location, y = drug_percent, fill=drug_type_fct))  + geom_col() +  
  
  scale_fill_manual(values = drugType_pallette) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  
  xlab("") + ylab("") + guides(fill = guide_legend(nrow = 2)) +
  
  theme(legend.position="bottom") +
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank()) +
  
  scale_x_discrete(labels= sub("_Grass","",unique(plotDat_Grass$location))) + 
  
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black")) 

p_Grass

#
layout <- "
AAAAAAABBBBBBBCCCCCCDDE
"
p_Sewage + p_Farmland + p_Forest + p_Grass + p_Gobi + plot_layout(design = layout)
