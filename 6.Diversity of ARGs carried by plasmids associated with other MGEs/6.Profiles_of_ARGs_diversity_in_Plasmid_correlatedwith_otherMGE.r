

# library ----------------------------
library(ggplot2)
library(scales)
library(ggpubr)
library(ggprism)
library(tidyverse)
library(hrbrthemes)
library(ggthemes)
library(dplyr)
library(ggpattern)
library(Hmisc)
library(reshape2)
library(venn)
library(ggvenn)
library(stringr)
library(ragg)

# 1. bar ------------------------------

# 
load("1_Combined_results_diversity.RData")
load("1_Ranking.RData")

# 过滤出overall
Overall_df <- combined_df %>% filter(drug_type == "Overall")

# 
Overall_df$MGE_type_final <- factor(Overall_df$MGE_type_final, levels = c( 'IS',"Transposon",'Integron'))
Overall_df$sampleType <-  factor(Overall_df$sampleType, levels = c('Tailings','Sewage','Farmland','Forest','Grass','Gobi'))

p_total <- Overall_df  %>%
  ggplot() +
  geom_errorbar(aes(x = MGE_type_final, fill = MGE_type_final,
                    ymin = num_percent_mean-sd, ymax = num_percent_mean+sd),
                width = 0.5,
                colour="black", alpha=1, size=0.3) +
  geom_bar( aes(x = MGE_type_final, fill = MGE_type_final, weight = num_percent_mean),
            position = 'dodge',   color = NA,  width = 0.85) + 
  scale_fill_manual(values = c('#d0c0a5','#2e7d32' , "#dce775")) +
  scale_y_continuous(expand = c(0,0)) + # y = 0
  facet_wrap(sampleType~drug_type, nrow =3, scales = "fixed")+
  theme_prism(base_line_size = 0.3, border = T)  + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  theme(legend.position = 'none') + xlab("") +ylim(0, 0.51) + 
  ylab("Ratio of plasmid-transposon/IS/integron-carrying \n ARG subtypes to plasmid-borne ARG subtypes") 

p_total


# 2. venn -------------------------

# load
load("1_plasmids_Tn_IS_Integron_commonDf.RData")
# check data
unique(plasmids_Integron_Tn_IS_df$MGE_type_final) # "Integron"   "Transposon" "IS"

#
Venn_df <- plasmids_Integron_Tn_IS_df %>% 
  filter(DepthPG > 0) %>% 
  select(ARG, MGE_type_final,sampleType) %>% unique()

unique(Venn_df$MGE_type_final) # "Integron"   "Transposon" "IS"

# 
Venn_df_Integron <- Venn_df %>% filter(MGE_type_final == "Integron") 
colnames(Venn_df_Integron)[colnames(Venn_df_Integron) == "MGE_type_final"] <- "Integron"

Venn_df_Transposon <- Venn_df %>% filter(MGE_type_final == "Transposon")
colnames(Venn_df_Transposon)[colnames(Venn_df_Transposon) == "MGE_type_final"] <- "Transposon"

Venn_df_IS <- Venn_df %>% filter(MGE_type_final == "IS")
colnames(Venn_df_IS)[colnames(Venn_df_IS) == "MGE_type_final"] <- "IS"

Venn_combined_df <- merge(merge(Venn_df_Integron, Venn_df_Transposon, by = c("sampleType", "ARG"), all = TRUE), 
                          Venn_df_IS,  by = c("sampleType", "ARG"), all = TRUE)

# save(Venn_combined_df, file = "1_Venn_combined_df.RData")

load("1_Venn_combined_df.RData")

# 
unique(Venn_combined_df$sampleType) # "Farmland" "Forest"   "Gobi"     "Grass"    "Sludge"   "Tailings"

# Tailings
Venn_combined_df_Tailings_Integron <- Venn_combined_df %>% filter(sampleType == "Tailings") %>% 
  filter(Integron == "Integron")
aa <- Venn_combined_df_Tailings_Integron$ARG

Venn_combined_df_Tailings_Transposon <- Venn_combined_df %>% filter(sampleType == "Tailings") %>% 
  filter(Transposon == "Transposon")
bb <- Venn_combined_df_Tailings_Transposon$ARG

Venn_combined_df_Tailings_IS <- Venn_combined_df %>% filter(sampleType == "Tailings") %>% 
  filter(IS == "IS")
cc <- Venn_combined_df_Tailings_IS$ARG

x <- list(Integron=aa, Transposon=bb, IS=cc)
x

p_Tailings <-ggvenn::ggvenn(x,c("IS","Transposon","Integron"), fill_color = c("#d0c0a5","#2e7d32", "#dce775"), 
                            fill_alpha = 0.6,#改变图形透明度
                            stroke_color="#616161",#边界线的颜色
                            stroke_alpha = 0.5,#边界线的透明度
                            text_size = 6)#调节字体大小
p_Tailings <- p_Tailings + 
  labs(subtitle = "Tailings") + 
  theme(plot.subtitle = element_text(size = 16, colour = "black"))

p_Tailings

# sewage
Venn_combined_df_Sludge_Integron <- Venn_combined_df %>% filter(sampleType == "Sludge") %>% filter(Integron == "Integron")
aa <- Venn_combined_df_Sludge_Integron$ARG

Venn_combined_df_Sludge_Transposon <- Venn_combined_df %>% filter(sampleType == "Sludge") %>% filter(Transposon == "Transposon")
bb <- Venn_combined_df_Sludge_Transposon$ARG

Venn_combined_df_Sludge_IS <- Venn_combined_df %>% filter(sampleType == "Sludge") %>% filter(IS == "IS")
cc <- Venn_combined_df_Sludge_IS$ARG

x <- list(Integron=aa, Transposon=bb, IS=cc)
x

p_Sludge <-ggvenn::ggvenn(x,c("IS","Transposon","Integron"), fill_color = c("#d0c0a5","#2e7d32", "#dce775"),
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
Venn_combined_df_farmland_Integron <- Venn_combined_df %>% filter(sampleType == "Farmland") %>% filter(Integron == "Integron")
aa <- Venn_combined_df_farmland_Integron$ARG

Venn_combined_df_farmland_Transposon <- Venn_combined_df %>% filter(sampleType == "Farmland") %>% filter(Transposon == "Transposon")
bb <- Venn_combined_df_farmland_Transposon$ARG

Venn_combined_df_farmland_IS <- Venn_combined_df %>% filter(sampleType == "Farmland") %>% filter(IS == "IS")
cc <- Venn_combined_df_farmland_IS$ARG

x <- list(Integron=aa, Transposon=bb, IS=cc)
x

p_farmland <-ggvenn::ggvenn(x,c("IS","Transposon","Integron"), fill_color = c("#d0c0a5","#2e7d32", "#dce775"),
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
Venn_combined_df_forest_Integron <- Venn_combined_df %>% filter(sampleType == "Forest") %>% filter(Integron == "Integron")
aa <- Venn_combined_df_forest_Integron$ARG

Venn_combined_df_forest_Transposon <- Venn_combined_df %>% filter(sampleType == "Forest") %>% filter(Transposon == "Transposon")
bb <- Venn_combined_df_forest_Transposon$ARG

Venn_combined_df_forest_IS <- Venn_combined_df %>% filter(sampleType == "Forest") %>% filter(IS == "IS")
cc <- Venn_combined_df_forest_IS$ARG

x <- list(Integron=aa, Transposon=bb, IS=cc)
x

p_forest <-ggvenn::ggvenn(x,c("IS","Transposon","Integron"), fill_color = c("#d0c0a5","#2e7d32", "#dce775"), 
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
Venn_combined_df_grass_Integron <- Venn_combined_df %>% filter(sampleType == "Grass") %>% filter(Integron == "Integron")
aa <- Venn_combined_df_grass_Integron$ARG

Venn_combined_df_grass_Transposon <- Venn_combined_df %>% filter(sampleType == "Grass") %>% filter(Transposon == "Transposon")
bb <- Venn_combined_df_grass_Transposon$ARG

Venn_combined_df_grass_IS <- Venn_combined_df %>% filter(sampleType == "Grass") %>% filter(IS == "IS")
cc <- Venn_combined_df_grass_IS$ARG

x <- list(Integron=aa, Transposon=bb, IS=cc)
x

p_grass <-ggvenn::ggvenn(x,c("IS","Transposon","Integron"), fill_color = c("#d0c0a5","#2e7d32", "#dce775"),
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
Venn_combined_df_gobi_Integron <- Venn_combined_df %>% filter(sampleType == "Gobi") %>% filter(Integron == "Integron")
aa <- Venn_combined_df_gobi_Integron$ARG

Venn_combined_df_gobi_Transposon <- Venn_combined_df %>% filter(sampleType == "Gobi") %>% filter(Transposon == "Transposon")
bb <- Venn_combined_df_gobi_Transposon$ARG

Venn_combined_df_gobi_IS <- Venn_combined_df %>% filter(sampleType == "Gobi") %>% filter(IS == "IS")
cc <- Venn_combined_df_gobi_IS$ARG

x <- list(Integron=aa, Transposon=bb, IS=cc)
x

p_gobi <-ggvenn::ggvenn(x,c("IS","Transposon","Integron"), fill_color = c("#d0c0a5","#2e7d32", "#dce775"),
                        fill_alpha = 0.6,#改变图形透明度
                        stroke_color="#616161",#边界线的颜色
                        stroke_alpha = 0.5,#边界线的透明度
                        text_size = 6)#调节字体大小
p_gobi <- p_gobi + 
  labs(subtitle = "Gobi") + 
  theme(plot.subtitle = element_text(size = 16,  
                                     colour = "black"))
p_gobi

# 3. circle bar ---------------------

# Tailings ----------------------
combined_df_T <- fread("PlotData.tailings.txt")

# 
datagroup <- combined_df_T$MGE_type_final %>% unique() # "IS" "Transposon"

allplotdata <- tibble('MGE_type_final' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(combined_df_T) %>% arrange(MGE_type_final) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1  10

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
allplotdata$MGE_type_final <- factor(allplotdata$MGE_type_final, levels = c("IS", "Transposon",  "Integron"))

p_Tailings <- ggplot() + 
  geom_bar(data = allplotdata, aes(x = xid, y = mean, fill = MGE_type_final), stat = 'identity') + 
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=4) + 
  coord_polar() + 
  #ylim(-0.0005,0.001) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.1) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=4 , angle=0) + 
  geom_segment(data = griddata, 
               aes(x = locationx-0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-2, 4)) +  
  scale_fill_manual(values=c('#d0c0a5','#2e7d32', "#dce775")) +
  theme_void() +
  theme(legend.position = 'none')

p_Tailings

# （2）Sewage -------------------------------------
combined_df_S <- fread("PlotData.sewage.txt")

# 对数据进行转换
datagroup <- combined_df_S$MGE_type_final %>% unique() # "IS" "Transposon"

allplotdata <- tibble('MGE_type_final' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(combined_df_T) %>% arrange(MGE_type_final) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1  10

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

# 绘图

allplotdata$MGE_type_final <- factor(allplotdata$MGE_type_final, levels = c("IS", "Transposon",  "Integron"))

p_Sewage <- ggplot() + 
  geom_bar(data = allplotdata, aes(x = xid, y = mean, fill = MGE_type_final), stat = 'identity') + 
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=4) + 
  coord_polar() + 
  #ylim(-0.0005,0.001) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.1) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=4 , angle=0) +
  geom_segment(data = griddata, 
               aes(x = locationx-0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-2, 4)) +  
  scale_fill_manual(values=c('#d0c0a5','#2e7d32', "#dce775")) +
  theme_void() +
  theme(legend.position = 'none')

p_Sewage

# （3）Farmland -----------------------------------
combined_df_Fa <- fread("PlotData.farmland.txt")
# 
datagroup <- combined_df_Fa$MGE_type_final %>% unique() # "IS" "Transposon"

allplotdata <- tibble('MGE_type_final' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(combined_df_T) %>% arrange(MGE_type_final) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1  10

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
allplotdata$MGE_type_final <- factor(allplotdata$MGE_type_final, levels = c("IS", "Transposon",  "Integron"))

p_Farmland <- ggplot() + 
  geom_bar(data = allplotdata, aes(x = xid, y = mean, fill = MGE_type_final), stat = 'identity') + 
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=4) + 
  coord_polar() + 
  #ylim(-0.0005,0.001) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.1) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=4 , angle=0) +
  geom_segment(data = griddata, 
               aes(x = locationx-0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-2, 4)) +  
  scale_fill_manual(values=c('#d0c0a5','#2e7d32', "#dce775")) +
  theme_void() +
  theme(legend.position = 'none')

p_Farmland

# （4）Forest -------------------------------------

combined_df_Fo <- fread("PlotData.forest.txt")
# 
datagroup <- combined_df_Fo$MGE_type_final %>% unique() # "IS" "Transposon"

allplotdata <- tibble('MGE_type_final' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(combined_df_T) %>% arrange(MGE_type_final) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1  10

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
allplotdata$MGE_type_final <- factor(allplotdata$MGE_type_final, levels = c("IS", "Transposon",  "Integron"))

p_Forest <- ggplot() + 
  geom_bar(data = allplotdata, aes(x = xid, y = mean, fill = MGE_type_final), stat = 'identity') + 
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=4) + 
  coord_polar() + 
  #ylim(-0.0005,0.001) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.1) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=4 , angle=0) + 
  geom_segment(data = griddata, 
               aes(x = locationx-0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-2, 4)) +  
  scale_fill_manual(values=c('#d0c0a5','#2e7d32', "#dce775")) +
  theme_void() +
  theme(legend.position = 'none')

p_Forest

# （5）Grass -----------------------------------
combined_df_Gr <- fread("PlotData.grass.txt")

# 
datagroup <- combined_df_Gr $MGE_type_final %>% unique() # "IS" "Transposon"

allplotdata <- tibble('MGE_type_final' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(combined_df_T) %>% arrange(MGE_type_final) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

# 
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1  10

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
allplotdata$MGE_type_final <- factor(allplotdata$MGE_type_final, levels = c("IS", "Transposon",  "Integron"))

p_Grass <- ggplot() + 
  geom_bar(data = allplotdata, aes(x = xid, y = mean, fill = MGE_type_final), stat = 'identity') + 
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=4) + 
  coord_polar() + 
  #ylim(-0.0005,0.001) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.1) +  
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) + 
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=4 , angle=0) + 
  geom_segment(data = griddata, 
               aes(x = locationx-0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-2, 4)) +  
  scale_fill_manual(values=c('#d0c0a5','#2e7d32', "#dce775")) +
  theme_void() +
  theme(legend.position = 'none')

p_Grass

# （6）Gobi -----------------------------------
combined_df_Go <- fread("PlotData.grass.txt")
datagroup <- combined_df_Go$MGE_type_final %>% unique() # "IS" "Transposon"

allplotdata <- tibble('MGE_type_final' = datagroup,
                      'ARG' = paste0('empty_ARG_', seq_along(datagroup)),
                      'mean' = 0) %>% 
  bind_rows(combined_df_T) %>% arrange(MGE_type_final) %>% dplyr::mutate(xid = 1:n()) %>% 
  dplyr::mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  dplyr::mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  dplyr::mutate(angle = ifelse(angle < -90, angle+180, angle)) 

#
firstxid <- which(str_detect(allplotdata$ARG, pattern = "empty_ARG")) #  1  10

segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))

#
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$mean), to = max(allplotdata$mean), 0.5),
                 'coordytext' = as.character(round(coordylocation, 4)), #
                 'x' = 1) #

# 
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

#

allplotdata$MGE_type_final <- factor(allplotdata$MGE_type_final, levels = c("IS", "Transposon",  "Integron"))

p_Gobi <- ggplot() + 
  geom_bar(data = allplotdata, aes(x = xid, y = mean, fill = MGE_type_final), stat = 'identity') + 
  geom_text(data = allplotdata %>% filter(!str_detect(ARG, pattern = "empty_ARG")), 
            aes(x = xid, label = ARG, y = mean + 0.001, angle = angle, hjust = hjust),
            color="black", size=4) + 
  coord_polar() + 
  #ylim(-0.0005,0.001) +
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -0.1, yend=-0.1) + 
  #geom_text(data = segment_data, aes(x = labelx, label = label), y = -0.0015) +
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=4 , angle=0) + 
  geom_segment(data = griddata, 
               aes(x = locationx-0.5, xend = locationx + 0.5, y = locationy, yend = locationy), 
               colour = "black", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-2, 4)) +  
  scale_fill_manual(values=c('#d0c0a5','#2e7d32', "#dce775")) +
  theme_void() +
  theme(legend.position = 'none')

p_Gobi

# 
library(patchwork)

p_Tailings + p_Sewage + p_Farmland + p_Forest + p_Grass + p_Gobi






