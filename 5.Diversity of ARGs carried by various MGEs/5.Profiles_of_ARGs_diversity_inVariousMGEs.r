
# library
library(ggplot2)
library(scales)
library(ggpubr)
library(ggprism)
library(tidyverse)
library(hrbrthemes)
library(ggthemes)
library(dplyr)
library(reshape2)
library(venn)
library(ggvenn)
library(stringr)
library(ragg)
library(ggpattern)

# 1. transfer MGEs --------------------------------------------------

# 
load("1_Combined_results_diversity_total.RData")
load("1_Ranking_total.RData")

# 过滤出overall
Overall_df <- combined_df %>% filter(drug_type == "Overall")

# 
Overall_df$MGE_type_final <- factor(Overall_df$MGE_type_final, levels = c( 'Plasmid',"IS","Transposon","Phage", "Integron",'ICE'))
Overall_df$sampleType <-  factor(Overall_df$sampleType, levels = c('Tailings','Sewage','Farmland','Forest','Grass','Gobi'))

p_total <- Overall_df %>%
  ggplot() +
  geom_errorbar(aes(x = MGE_type_final, fill = MGE_type_final,
                    ymin = num_percent_mean-sd, ymax = num_percent_mean+sd),
                width = 0.5,
                colour="black", alpha=1, size=0.3) + 
  geom_bar( aes(x = MGE_type_final, fill = MGE_type_final, weight = num_percent_mean),
            position = 'dodge',   color = NA,  width = 0.85) + 
  scale_fill_manual(values = c("#d7660d","#d0c0a5","#2e7d32" ,"#dce775", "#f1e54c","#1976d2")) + 
  scale_y_continuous(expand = c(0,0)) + # y = 0
  facet_wrap(sampleType~drug_type, nrow =2, scales = "fixed") + 
  theme_prism(base_line_size = 0.3, border = T)  + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(legend.position = 'none') + xlab("") +ylim(0, 1) + 
  ylab("Ratio of MGE-carrying ARG subtypes to the total ARG subtypes") 

p_total


# 2. different Plasmids --------------------------------------------

load("1_Combined_results_diversity_plasmid.RData")
load("1_Ranking_plasmid.RData")

# 过滤出overall
Overall_df <- combined_df %>% filter(drug_type == "Overall")

# 
Overall_df$Plasmid_type <- factor(Overall_df$Plasmid_type, levels = c( 'Unmobilizable',"Mobilizable",'Conjugative'))
Overall_df$sampleType <-  factor(Overall_df$sampleType, levels = c('Tailings','Sewage','Farmland','Forest','Grass','Gobi'))

p_total <- Overall_df %>% 
  ggplot() +
  geom_errorbar(aes(x = Plasmid_type, fill = Plasmid_type,
                    ymin = num_percent_mean-sd, ymax = num_percent_mean+sd),
                width = 0.5,
                colour="black", alpha=1, size=0.3) + 
  geom_bar( aes(x = Plasmid_type, fill = Plasmid_type, weight = num_percent_mean),
            position = 'dodge',   color = NA,  width = 0.85) + 
  scale_fill_manual(values = c("#497549", "#e9a91f", "#ce3d36")) + 
  scale_y_continuous(expand = c(0,0)) + # y = 0
  facet_wrap(sampleType~drug_type, nrow =2, scales = "fixed") + 
  theme_prism(base_line_size = 0.3, border = T)  + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(legend.position = 'none') + xlab("") +ylim(0, 1) + 
  ylab("Ratio of plasmid-carrying ARG subtypes to the total ARG subtypes") 

p_total

# 3.venn ---------------------------------------------------
load("1_Venn_combined_df.RData")

Venn_combined_df <- Venn_combined_df %>% filter(sampleType == "Tailings")

# 
unique(Venn_combined_df$sampleType) # "farmland" "forest"   "tailings" "Sludge"   "gobi"     "grass"

# Tailings
Venn_combined_df_Tailings_plasmids <- Venn_combined_df %>% filter(sampleType == "Tailings") %>% filter(Plamids == "Plasmids") # 588
aa <- Venn_combined_df_Tailings_plasmids$ARG %>% unique()

Venn_combined_df_Tailings_phage <- Venn_combined_df %>% filter(sampleType == "Tailings") %>% filter(Phage == "Phage")
bb <- Venn_combined_df_Tailings_phage$ARG

Venn_combined_df_Tailings_ICE <- Venn_combined_df %>% filter(sampleType == "Tailings") %>% filter(ICE == "ICE")
cc <- Venn_combined_df_Tailings_ICE$ARG

x <- list(Plasmids=aa, Phage=bb, ICE=cc)
x

p_tailings <- ggvenn::ggvenn(x,c("Plasmids","ICE","Phage"), fill_color = c("#d7660d", "#1976d2","#fdd835"), 
                             fill_alpha = 0.6,#改变图形透明度
                             stroke_color="#616161",#边界线的颜色
                             stroke_alpha = 0.5,#边界线的透明度
                             text_size = 6)#调节字体大小
p_tailings <- p_tailings + 
  labs(subtitle = "Tailings") + 
  theme(plot.subtitle = element_text(size = 16, colour = "black"))

p_tailings

# sewage
Venn_combined_df_Sludge_plasmids <- Venn_combined_df %>% filter(sampleType == "Sludge") %>% filter(Plamids == "Plasmids")
aa <- Venn_combined_df_Sludge_plasmids$ARG

Venn_combined_df_Sludge_phage <- Venn_combined_df %>% filter(sampleType == "Sludge") %>% filter(Phage == "Phage")
bb <- Venn_combined_df_Sludge_phage$ARG

Venn_combined_df_Sludge_ICE <- Venn_combined_df %>% filter(sampleType == "Sludge") %>% filter(ICE == "ICE")
cc <- Venn_combined_df_Sludge_ICE$ARG

x <- list(Plasmids=aa, Phage=bb, ICE=cc)
x

p_Sludge <-ggvenn::ggvenn(x,c("Plasmids","ICE","Phage"), fill_color = c("#d7660d", "#1976d2","#fdd835"), 
                          fill_alpha = 0.6,#改变图形透明度
                          stroke_color="#616161",#边界线的颜色
                          stroke_alpha = 0.5,#边界线的透明度
                          text_size = 6)#调节字体大小
p_Sludge <- p_Sludge + 
  labs(subtitle = "Sewage") + 
  theme(plot.subtitle = element_text(size = 16, face = "bold", 
                                     colour = "black"))
p_Sludge

# farmland
Venn_combined_df_farmland_plasmids <- Venn_combined_df %>% filter(sampleType == "Farmland") %>% filter(Plamids == "Plasmids")
aa <- Venn_combined_df_farmland_plasmids$ARG

Venn_combined_df_farmland_phage <- Venn_combined_df %>% filter(sampleType == "Farmland") %>% filter(Phage == "Phage")
bb <- Venn_combined_df_farmland_phage$ARG

Venn_combined_df_farmland_ICE <- Venn_combined_df %>% filter(sampleType == "Farmland") %>% filter(ICE == "ICE")
cc <- Venn_combined_df_farmland_ICE$ARG

x <- list(Plasmids=aa, Phage=bb, ICE=cc)
x

p_farmland <-ggvenn::ggvenn(x,c("Plasmids","ICE","Phage"), fill_color = c("#d7660d", "#1976d2","#fdd835"), 
                            fill_alpha = 0.6,#改变图形透明度
                            stroke_color="#616161",#边界线的颜色
                            stroke_alpha = 0.5,#边界线的透明度
                            text_size = 6)#调节字体大小
p_farmland <- p_farmland + 
  labs(subtitle = "Farmland") + 
  theme(plot.subtitle = element_text(size = 16, face = "bold", 
                                     colour = "black"))
p_farmland

# forest
Venn_combined_df_forest_plasmids <- Venn_combined_df %>% filter(sampleType == "Forest") %>% filter(Plamids == "Plasmids")
aa <- Venn_combined_df_forest_plasmids$ARG

Venn_combined_df_forest_phage <- Venn_combined_df %>% filter(sampleType == "Forest") %>% filter(Phage == "Phage")
bb <- Venn_combined_df_forest_phage$ARG

Venn_combined_df_forest_ICE <- Venn_combined_df %>% filter(sampleType == "Forest") %>% filter(ICE == "ICE")
cc <- Venn_combined_df_forest_ICE$ARG

x <- list(Plasmids=aa, Phage=bb, ICE=cc)
x

p_forest <-ggvenn::ggvenn(x,c("Plasmids","ICE","Phage"), fill_color = c("#d7660d", "#1976d2","#fdd835"), 
                          fill_alpha = 0.6,#改变图形透明度
                          stroke_color="#616161",#边界线的颜色
                          stroke_alpha = 0.5,#边界线的透明度
                          text_size = 6)#调节字体大小
p_forest <- p_forest + 
  labs(subtitle = "Forest") + 
  theme(plot.subtitle = element_text(size = 16, face = "bold", 
                                     colour = "black"))
p_forest

# grass
Venn_combined_df_grass_plasmids <- Venn_combined_df %>% filter(sampleType == "Grass") %>% filter(Plamids == "Plasmids")
aa <- Venn_combined_df_grass_plasmids$ARG

Venn_combined_df_grass_phage <- Venn_combined_df %>% filter(sampleType == "Grass") %>% filter(Phage == "Phage")
bb <- Venn_combined_df_grass_phage$ARG

Venn_combined_df_grass_ICE <- Venn_combined_df %>% filter(sampleType == "Grass") %>% filter(ICE == "ICE")
cc <- Venn_combined_df_grass_ICE$ARG

x <- list(Plasmids=aa, Phage=bb, ICE=cc)
x

p_grass <-ggvenn::ggvenn(x,c("Plasmids","ICE","Phage"), fill_color = c("#d7660d", "#1976d2","#fdd835"), 
                         fill_alpha = 0.6,#改变图形透明度
                         stroke_color="#616161",#边界线的颜色
                         stroke_alpha = 0.5,#边界线的透明度
                         text_size = 6)#调节字体大小
p_grass <- p_grass + 
  labs(subtitle = "Grass") + 
  theme(plot.subtitle = element_text(size = 16, face = "bold", 
                                     colour = "black"))
p_grass


# gobi
Venn_combined_df_gobi_plasmids <- Venn_combined_df %>% filter(sampleType == "Gobi") %>% filter(Plamids == "Plasmids")
aa <- Venn_combined_df_gobi_plasmids$ARG

Venn_combined_df_gobi_phage <- Venn_combined_df %>% filter(sampleType == "Gobi") %>% filter(Phage == "Phage")
bb <- Venn_combined_df_gobi_phage$ARG

Venn_combined_df_gobi_ICE <- Venn_combined_df %>% filter(sampleType == "Gobi") %>% filter(ICE == "ICE")
cc <- Venn_combined_df_gobi_ICE$ARG

x <- list(Plasmids=aa, Phage=bb, ICE=cc)
x

p_gobi <-ggvenn::ggvenn(x,c("Plasmids","ICE","Phage"), fill_color = c("#d7660d", "#1976d2","#fdd835"), 
                        fill_alpha = 0.6,#改变图形透明度
                        stroke_color="#616161",#边界线的颜色
                        stroke_alpha = 0.5,#边界线的透明度
                        text_size = 6)#调节字体大小
p_gobi <- p_gobi + 
  labs(subtitle = "Gobi") + 
  theme(plot.subtitle = element_text(size = 16,  
                                     colour = "black"))
p_gobi

library(patchwork)

p_tailings + p_Sludge + p_farmland + p_forest + p_grass + p_gobi

# 4. bar -----------------------------------------

# （1）Tailings --------

PlotData <- fread("PlotData.tailings.txt")

# 
datagroup <- PlotData$Plasmid_type %>% unique() # "Conj"  "ICE" "Phage"  "mob_unconj" "unmob" 

allplotdata <- tibble('Plasmid_type' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(PlotData) %>% arrange(Plasmid_type) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1 6  20

segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))

# 
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$mean), to = max(allplotdata$mean), 0.5),
                 'coordytext' = as.character(round(coordylocation, 4)), 
                 'x' = 1) # Y轴的范围

# 
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

# 绘图
allplotdata$Plasmid_type <- factor(allplotdata$Plasmid_type, levels = c("unmob","unconj","Conj", "Phage", "ICE"))

p_Tailings <- ggplot()+
  geom_bar_pattern(data = allplotdata, aes(x=xid , y=mean,pattern_type = Plasmid_type, pattern_fill = Plasmid_type, fill=Plasmid_type),
                   position = position_dodge(preserve = "single"),color = "white", stat = "identity",
                   #pattern_fill = "black",
                   pattern ='magick',
                   #pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   pattern_fill = "black") +
  scale_pattern_type_discrete(choices = c( 'right30', 'gray90','gray70', 'gray100','gray100'))+#values = c('stripe','circle','placeholder','none','crosshatch')
  scale_fill_manual(values=c("#d55e00",'#d55e00' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  coord_polar() +
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=3) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.11) +  # 
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=3, angle=0) + # 加坐标轴信息
  geom_segment(data = griddata, 
               aes(x = locationx - 0.5, xend = locationx + 0.5, y = locationy, yend = locationy), # 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-4, 5.5)) +  # 前面的控制的是圈的位置
  #scale_fill_manual(values=c("pink",'orange' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  theme_void() +
  theme(legend.position = 'none')

p_Tailings

# （2）Sewage ------------------------------------------------------------------

PlotData <- fread("PlotData.sewage.txt")
# 
datagroup <- PlotData$Plasmid_type %>% unique() # "Conj"  "ICE" "Phage"  "mob_unconj" "unmob" 

allplotdata <- tibble('Plasmid_type' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(PlotData) %>% arrange(Plasmid_type) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1 6  20

segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))

# 
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$mean), to = max(allplotdata$mean), 0.5),
                 'coordytext' = as.character(round(coordylocation, 4)), 
                 'x' = 1) 

# 自定义坐标轴的网格
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

# 绘图
allplotdata$Plasmid_type <- factor(allplotdata$Plasmid_type, levels = c("unmob","unconj","Conj", "Phage", "ICE"))

p_Sewage <- ggplot()+
  geom_bar_pattern(data = allplotdata, aes(x=xid , y=mean,pattern_type = Plasmid_type, pattern_fill = Plasmid_type, fill=Plasmid_type),
                   position = position_dodge(preserve = "single"),color = "white", stat = "identity",
                   #pattern_fill = "black",
                   pattern ='magick',
                   #pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   pattern_fill = "black") +
  scale_pattern_type_discrete(choices = c( 'right30', 'gray90','gray70', 'gray100','gray100'))+#values = c('stripe','circle','placeholder','none','crosshatch')
  scale_fill_manual(values=c("#d55e00",'#d55e00' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  coord_polar() +
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=3) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.11) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=3, angle=0) + 
  geom_segment(data = griddata, 
               aes(x = locationx - 0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-4, 5.5)) +  
  #scale_fill_manual(values=c("pink",'orange' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  theme_void() +
  theme(legend.position = 'none')

p_Sewage

# （3）Farmland ----------------------------------------------------------------
PlotData <- fread("PlotData.farmland.txt")

datagroup <- PlotData$Plasmid_type %>% unique() # "Conj"  "ICE" "Phage"  "mob_unconj" "unmob" 

allplotdata <- tibble('Plasmid_type' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(PlotData) %>% arrange(Plasmid_type) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 


firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1 6  20

segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))

# 自定坐标轴 
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$mean), to = max(allplotdata$mean), 0.5),
                 'coordytext' = as.character(round(coordylocation, 4)), 
                 'x' = 1) 

# 自定义坐标轴的网格
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

# 绘图
allplotdata$Plasmid_type <- factor(allplotdata$Plasmid_type, levels = c("unmob","unconj","Conj", "Phage", "ICE"))

p_Farmland <- ggplot()+
  geom_bar_pattern(data = allplotdata, aes(x=xid , y=mean,pattern_type = Plasmid_type, pattern_fill = Plasmid_type, fill=Plasmid_type),
                   position = position_dodge(preserve = "single"),color = "white", stat = "identity",
                   #pattern_fill = "black",
                   pattern ='magick',
                   #pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   pattern_fill = "black") +
  scale_pattern_type_discrete(choices = c( 'right30', 'gray90','gray70', 'gray100','gray100'))+#values = c('stripe','circle','placeholder','none','crosshatch')
  scale_fill_manual(values=c("#d55e00",'#d55e00' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  coord_polar() +
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=3) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.11) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + # unmob, conj等标签的位
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=3, angle=0) +
  geom_segment(data = griddata, 
               aes(x = locationx - 0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-4, 5.5)) + 
  #scale_fill_manual(values=c("pink",'orange' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  theme_void() +
  theme(legend.position = 'none')

p_Farmland

# （4）Forest ------------------------------------------------------------------
PlotData <- fread("PlotData.forest.txt")
# 
datagroup <- PlotData$Plasmid_type %>% unique() # "Conj"  "ICE" "Phage"  "mob_unconj" "unmob" 

allplotdata <- tibble('Plasmid_type' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(PlotData) %>% arrange(Plasmid_type) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1 6  20

segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))

# 
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$mean), to = max(allplotdata$mean), 0.5),
                 'coordytext' = as.character(round(coordylocation, 4)), 
                 'x' = 1) 

# 
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

# 
allplotdata$Plasmid_type <- factor(allplotdata$Plasmid_type, levels = c("unmob","unconj","Conj", "Phage", "ICE"))

p_Forest <- ggplot()+
  geom_bar_pattern(data = allplotdata, aes(x=xid , y=mean,pattern_type = Plasmid_type, pattern_fill = Plasmid_type, fill=Plasmid_type),
                   position = position_dodge(preserve = "single"),color = "white", stat = "identity",
                   #pattern_fill = "black",
                   pattern ='magick',
                   #pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   pattern_fill = "black") +
  scale_pattern_type_discrete(choices = c( 'right30', 'gray90','gray70', 'gray100','gray100'))+#values = c('stripe','circle','placeholder','none','crosshatch')
  scale_fill_manual(values=c("#d55e00",'#d55e00' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  coord_polar() +
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=3) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.11) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=3, angle=0) + 
  geom_segment(data = griddata, 
               aes(x = locationx - 0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-4, 5.5)) +  
  #scale_fill_manual(values=c("pink",'orange' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  theme_void() +
  theme(legend.position = 'none')

p_Forest

# （5）Grass -------------------------------------------------------------------
PlotData <- fread("PlotData.grass.txt")
# 
datagroup <- PlotData$Plasmid_type %>% unique() # "Conj"  "ICE" "Phage"  "mob_unconj" "unmob" 

allplotdata <- tibble('Plasmid_type' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(PlotData) %>% arrange(Plasmid_type) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1 6  20

segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))

#
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$mean), to = max(allplotdata$mean), 0.5),
                 'coordytext' = as.character(round(coordylocation, 4)), 
                 'x' = 1) 

# 自定义坐标轴的网格
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

# 绘图
allplotdata$Plasmid_type <- factor(allplotdata$Plasmid_type, levels = c("unmob","unconj","Conj", "Phage", "ICE"))

p_Grass <- ggplot()+
  geom_bar_pattern(data = allplotdata, aes(x=xid , y=mean,pattern_type = Plasmid_type, pattern_fill = Plasmid_type, fill=Plasmid_type),
                   position = position_dodge(preserve = "single"),color = "white", stat = "identity",
                   #pattern_fill = "black",
                   pattern ='magick',
                   #pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   pattern_fill = "black") +
  scale_pattern_type_discrete(choices = c( 'right30', 'gray90','gray70', 'gray100','gray100'))+#values = c('stripe','circle','placeholder','none','crosshatch')
  scale_fill_manual(values=c("#d55e00",'#d55e00' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  coord_polar() +
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=3) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.11) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=3, angle=0) + 
  geom_segment(data = griddata, 
               aes(x = locationx - 0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-4, 5.5)) +  
  #scale_fill_manual(values=c("pink",'orange' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  theme_void() +
  theme(legend.position = 'none')

p_Grass

# （6）Gobi -------------------------------------------------------------------
PlotData <- fread("PlotData.gobi.txt")

# 
datagroup <- PlotData$Plasmid_type %>% unique() # "Conj"  "ICE" "Phage"  "mob_unconj" "unmob" 

allplotdata <- tibble('Plasmid_type' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(PlotData) %>% arrange(Plasmid_type) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1 6  20

segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))

#
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$mean), to = max(allplotdata$mean), 0.5),
                 'coordytext' = as.character(round(coordylocation, 4)), 
                 'x' = 1) 

# 
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

# 
allplotdata$Plasmid_type <- factor(allplotdata$Plasmid_type, levels = c("unmob","unconj","Conj", "Phage", "ICE"))

p_Gobi <- ggplot()+
  geom_bar_pattern(data = allplotdata, aes(x=xid , y=mean,pattern_type = Plasmid_type, pattern_fill = Plasmid_type, fill=Plasmid_type),
                   position = position_dodge(preserve = "single"),color = "white", stat = "identity",
                   #pattern_fill = "black",
                   pattern ='magick',
                   #pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   pattern_fill = "black") +
  scale_pattern_type_discrete(choices = c( 'right30', 'gray90','gray70', 'gray100','gray100'))+#values = c('stripe','circle','placeholder','none','crosshatch')
  scale_fill_manual(values=c("#d55e00",'#d55e00' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  coord_polar() +
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=3) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.11) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=3, angle=0) + 
  geom_segment(data = griddata, 
               aes(x = locationx - 0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-4, 5.5)) +  
  #scale_fill_manual(values=c("pink",'orange' ,"#d55e00",'#f0e442' ,"#1976d2")) +
  theme_void() +
  theme(legend.position = 'none')

p_Gobi

library(patchwork)

p_Tailings + p_Sewage + p_Farmland + p_Forest + p_Grass + p_Gobi




