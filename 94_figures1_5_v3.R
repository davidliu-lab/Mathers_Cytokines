library(tidyr)
library(rpart)
library(rpart.plot)
library(epitools)
library(sjlabelled)
library(genetics)
library(caret)
library(ROCit)
library(pROC)
library(tibble)
library(readxl)
library(readr)
library(effects)
library(sjPlot)
library(DT)
library(gtsummary)
library(gt)
library(ggthemr)
library(ggcorrplot)
library(ggpmisc)
library(fastcluster)
library(purrr)
library(grid)
library(ComplexHeatmap)
library(uwot)
library(stringr)
library(readxl)
library(scater)
library(dplyr)
library(S4Vectors)
library(ComplexHeatmap)
library(circlize)
library(flowCore)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggpubr)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(diffcyt)
require(data.table)
library(xlsx)
library(factoextra)
library(seriation)
library(survminer)
library(survival)
library(ggrepel)
library(ggcorrplot)
library(corrplot)
library(here)
library(gridExtra)
library(grid)
library(scales)
library(npsm)
library(plyr)
library(ggbreak)
library(gridtext)
library(dict)
library(VennDiagram)
library(RColorBrewer)
library(ggrisktable)

#color coding -------------------------------------------------------------------------------
dec.color <- "#B5CDA3"
inc.color <- "#ECDBBA"
DCB.color <- "#3C548899"
NCB.color <- "#B09C8599"
priorIO.N.color <- "#CEB5B9"
priorIO.Y.color <- "#98A1B1"

#Figure1A------------------------------------------------------------------------------------



cytokines.3t.hm.t <- Heatmap(
  matrix = tmp.pval.3t.reci.mtx.t,
  name= "ad_cd3_func_pre_YvsN",
  cluster_columns = FALSE,cluster_rows = FALSE,
  
  col = colorRamp2(c(-20,0,20),
                   c(NCB.color,"white",DCB.color)),
  column_title_side = ifelse(length(by) ==2,"top", "bottom"),
  show_row_dend = FALSE,
  show_column_dend = FALSE, 
  row_names_side = "left",
  column_names_rot =  30,column_names_side = "top",
  row_dend_width = unit(2, "cm"),
  rect_gp = gpar(col = "black"),
  width = unit(12,"cm"),
  height=unit(3,"cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(between(abs(tmp.pval.3t.reci.mtx.t[i, j]),1/0.05,1/0.01)==TRUE)
      grid.text(sprintf("*",tmp.pval.3t.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (between(abs(tmp.pval.3t.reci.mtx.t[i, j]),1/0.01,1/0.001)==TRUE)
      grid.text(sprintf("**",tmp.pval.3t.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (between(abs(tmp.pval.3t.reci.mtx.t[i, j]),1/0.001,1/0.0001)==TRUE)
      grid.text(sprintf("***",tmp.pval.3t.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (abs(tmp.pval.3t.reci.mtx.t[i, j]) >= 1/0.0001)
      grid.text(sprintf("****",tmp.pval.3t.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    grid_height = unit(1, "cm"), grid_width = unit(1, "mm"),direction = "vertical",
    at = c(-20,-10,10,20),
    labels = c("p <= 0.05","p <= 0.1",
               "p <=0.1","p <=0.05"),
    title = "p value"
  )
)

#Figure1B -------------------------------------------------------------------------
luminex_c1c2_df_v0117_v2 %>%
  dplyr::mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  dplyr::mutate(CXCL5 = log2(CXCL5)) %>%
  dplyr::select(CXCL5,DCB,timepoint) %>%
  filter(timepoint == "pre") %>%
  ggplot(aes(x=factor(DCB,levels = c("Y","N")),
             y=CXCL5, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="", y = "Log2(CXCL5)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",hide.ns = TRUE,size = 20,
                     label.y.npc = 0.85,
                     label.x.npc = 0.5,
                     show.legend = FALSE
  )+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))

luminex_c1c2_df_v0117_v2 %>%
  mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  mutate(IL.23 = log2(IL.23)) %>%
  dplyr::select(IL.23,DCB,timepoint) %>%
  filter(timepoint == "pre") %>%
  ggplot(aes(x=factor(DCB,levels = c("Y","N")),
             y=IL.23, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="", y = "Log2(IL.23)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",hide.ns = TRUE,size = 20,
                     label.y.npc = 0.85,
                     label.x.npc = 0.5,
                     show.legend = FALSE
  )+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))
#p = 0.044

luminex_c1c2_df_v0117_v2 %>%
  dplyr::mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  dplyr::mutate(CXCL6 = log2(CXCL6)) %>%
  dplyr::select(CXCL6,DCB,timepoint) %>%
  dplyr::filter(timepoint == "pre") %>%
  ggplot(aes(x=factor(DCB,levels = c("Y","N")),
             y=CXCL6, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="", y = "Log2(CXCL6)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",hide.ns = TRUE,size = 20,
                     label.y.npc = 0.85,
                     label.x.npc = 0.5,
                     show.legend = FALSE
  )+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))

luminex_c1c2_df_v0117_v2 %>%
  mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  mutate(IL10 = log2(IL10)) %>%
  dplyr::select(IL10,DCB,timepoint) %>%
  dplyr::filter(timepoint == "pre") %>%
  ggplot(aes(x=factor(DCB,levels = c("Y","N")),
             y=IL10, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="", y = "Log2(IL10)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",hide.ns = TRUE,size = 20,
                     label.y.npc = 0.85,
                     label.x.npc = 0.5,
                     show.legend = FALSE
  )+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))
#p = 0.019



#Figure1C -------------------------------------------------------------------------
cytokines.longitudinal.hm.t <- Heatmap(
  matrix = lcr.mtx.reci.mtx.t,
  name= "ad_cd3_func_pre_YvsN",
  cluster_columns = FALSE,cluster_rows = FALSE,
  col = colorRamp2(c(-20,0,20),
                   c(dec.color,"white",inc.color)),
  column_title_side = ifelse(length(by) ==2,"top", "bottom"),
  show_row_dend = FALSE,
  show_column_dend = FALSE, 
  row_names_side = "left",
  column_names_rot =  30,column_names_side = "top",
  row_dend_width = unit(2, "cm"),
  rect_gp = gpar(col = "black"),
  width = unit(12,"cm"),
  height=unit(3,"cm"),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(between(abs(lcr.mtx.reci.mtx.t[i, j]),1/0.05,1/0.01)==TRUE)
      grid.text(sprintf("*",lcr.mtx.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (between(abs(lcr.mtx.reci.mtx.t[i, j]),1/0.01,1/0.001)==TRUE)
      grid.text(sprintf("**",lcr.mtx.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (between(abs(lcr.mtx.reci.mtx.t[i, j]),1/0.001,1/0.0001)==TRUE)
      grid.text(sprintf("***",lcr.mtx.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (abs(lcr.mtx.reci.mtx.t[i, j]) >= 1/0.0001)
      grid.text(sprintf("****",lcr.mtx.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))},
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    grid_height = unit(1, "cm"), grid_width = unit(1, "mm"),direction = "vertical",
    at = c(-20,-10,10,20),
    labels = c("p <= 0.001","p <= 0.1", 
               "p <=0.1","p <=0.001"),
    title = "p value"
  )
)



#Figure2A ---------------------------------------------------------------------------------------
c1c2_prop_df_v1_labs_updated %>%
  mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  dplyr::select(Abs_Lymphs,DCB,timepoint) %>%
  mutate(Abs_Lymphs = log2(Abs_Lymphs)) %>%
  ggplot(aes(x=factor(timepoint,levels = c("pre","post1","post2")),
             y=Abs_Lymphs, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="Time point", y = "Log2(Absolute Lymphocytes)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",hide.ns = TRUE,size = 20,
                     label.y.npc = 0.85,show.legend = FALSE)+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))+ 
  geom_segment(aes(x = 1.9, y = 1.8, xend = 2.1, yend = 1.8, colour = "black"))

c1c2_prop_df_v1_labs_updated %>%
  mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  mutate(Abs_Neuts = log2(Abs_Neuts)) %>%
  dplyr::select(Abs_Neuts,DCB,timepoint) %>%
  ggplot(aes(x=factor(timepoint,levels = c("pre","post1","post2")),
             y=Abs_Neuts, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="Time point", y = "Log2(Absolute Neutrophils)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",hide.ns = TRUE,
                     size=8,label.y.npc = 0.85,show.legend = FALSE)+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))+ geom_text(
    label="p = 0.13", x= 2,y=4,color = "black", size = 10
  )

c1c2_prop_df_v1_labs_updated %>%
  mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  dplyr::select(NLR,DCB,timepoint) %>%
  mutate(NLR = log2(NLR)) %>%
  ggplot(aes(x=factor(timepoint,levels = c("pre","post1","post2")),
             y=NLR, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="Time point", y = "log2(NLR)") +
  stat_compare_means(method = "wilcox.test",size =20,
                     label = "p.signif",hide.ns = TRUE,
                     label.y.npc = 0.85,show.legend = FALSE)+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))+ 
  geom_segment(aes(x = 1.9, y = 5, xend = 2.1, yend = 5, colour = "black"))





#Figure2B ---------------------------------------------------------------------------------------
c1c2_prop_df_v1 %>%
  drop_na(Abs_Lymphs) %>%
  dplyr::mutate(Abs_Lymphs = log2(Abs_Lymphs)) %>%
  ggline(x = "timepoint",
         y = "Abs_Lymphs", add = "mean_se",position = position_dodge(0.05),
         color = "DCB",
         palette = c(NCB.color,DCB.color),
         stroke = 2,size=1,
         scales="free",
         order = c("pre", "post1", "post2"))+
  ylab("Log2(Absolute Lymphocytes)")+
  theme(text = element_text(size = 30))





#Figure2C ---------------------------------------------------------------------------------------
#PFS
c1c2_prop_df_v1_labs_fc_survival_PFS <- survfit(Surv(PFS/30,progress) ~ by_Abs_Lymphs_fc_median,
                                                data = c1c2_prop_df_v1_labs_fc_labs)
p <- ggsurvplot(c1c2_prop_df_v1_labs_fc_survival_PFS, 
                data = c1c2_prop_df_v1_labs_fc_labs,
                size =3, #size of the line (line type can be changed too)
                surv.median.line = "hv",
                title = "Progression Free Survival",
                xlab = "Time(Months)",
                legend.labs =
                  c(paste("High","FC of Abs Lymph",sep=""), paste("Low","FC of Abs Lymph",sep=" ")),   
                pval.method = TRUE,
                pval = TRUE,risk.table = TRUE, risk.table.y.text.col = TRUE,
                palette = c("#B7657B","#3A3845"),
                fontsize = 8,
                risk.table.y.text = TRUE,
                font.title = c(30, "bold", "darkblue"),
                font.subtitle = c(30, "bold.italic", "purple"),
                font.caption = c(30, "plain", "orange"),
                font.x = c(30, "bold", "black"),
                font.y = c(30, "bold", "black"),
                font.main = c(30, "bold", "black"),
                font.tickslab = c(30, "plain", "black"),
                pval.size = 10,
                pval.coord = c(30,0.9),
                pval.method.coord = c(30,1)
)

p$plot <- p$plot +
  theme(legend.text = element_text(size = 20, color = "black", face = "bold"))+
  annotate("text", label = "High FC of Abs Lymph", x =30, y = 0.6, size = 10, colour = color.high )+
  annotate("text", label = "Low FC of Abs Lymph", x =20, y = 0.05, size = 10, colour = color.low)

p$table <- ggrisktable(c1c2_prop_df_v1_labs_fc_survival_PFS,
                       data = c1c2_prop_df_v1_labs_fc_labs,
                       color = "by_Abs_Lymphs_fc_median",
                       palette = c(color.high,color.low),
                       y.text = T,
                       fontsize=10,
                       y.text.col = TRUE, xlab = "",font.tickslab = c(20, "bold"),
                       tables.theme = theme_cleantable()
)+scale_y_discrete(labels=c("Low FC of Abs Lymph","High FC of Abs Lymph"))





#Figure2D ---------------------------------------------------------------------------------------
c1c2_prop_df_v1_labs_fc_survival_OS <- survfit(Surv(OS/30,death) ~ by_Abs_Lymphs_fc_median,
                                               data = c1c2_prop_df_v1_labs_fc_labs)
p <- ggsurvplot(c1c2_prop_df_v1_labs_fc_survival_OS, 
                data = c1c2_prop_df_v1_labs_fc_labs,
                size =3, 
                surv.median.line = "hv",
                title = "Overall Survival",
                xlab = "Time(Months)",
                legend.labs =
                  c(paste("High","FC of Abs Lymph",sep=" "), paste("Low","FC of Abs Lymph",sep=" ")),
                pval.method = TRUE,
                pval = TRUE,risk.table = TRUE, risk.table.y.text.col = TRUE,
                palette = c("#B7657B","#3A3845"),
                fontsize = 8,
                risk.table.y.text = TRUE,
                font.title = c(30, "bold", "darkblue"),
                font.subtitle = c(30, "bold.italic", "purple"),
                font.caption = c(30, "plain", "orange"),
                font.x = c(30, "bold", "black"),
                font.y = c(30, "bold", "black"),
                font.main = c(30, "bold", "black"),
                font.tickslab = c(30, "plain", "black"),
                pval.size = 10,
                pval.coord = c(30,0.9),
                pval.method.coord = c(30,1)
)

p$plot <- p$plot +
  theme(legend.text = element_text(size = 20, color = "black", face = "bold"))+
  annotate("text", label = "High FC of Abs Lymph", x =35, y = 0.75, size = 10, colour = color.high )+
  annotate("text", label = "Low FC of Abs Lymph", x =20, y = 0.2, size = 10, colour = color.low)
class(c1c2_prop_df_v1_labs_fc_labs$by_Abs_Lymphs_fc_median)
levels(c1c2_prop_df_v1_labs_fc_labs$by_Abs_Lymphs_fc_median)
p$table <- ggrisktable(c1c2_prop_df_v1_labs_fc_survival_OS,
                       data = c1c2_prop_df_v1_labs_fc_labs,
                       color = "by_Abs_Lymphs_fc_median",
                       palette = c(color.high,color.low),
                       y.text = T,
                       fontsize=10,
                       y.text.col = TRUE, xlab = "",font.tickslab = c(20, "bold"),
                       tables.theme = theme_cleantable()
)+scale_y_discrete(labels=c("Low FC of Abs Lymph","High FC of Abs Lymph"))




#Figure2E ---------------------------------------------------------------------------------------
##PFS
colnames(c1c2_prop_df_v1_labs_survival_post1)
c1c2_prop_df_v1_labs_survival_PFS <- survfit(Surv(PFS/30,progress) ~ by.NLR.median,
                                             data = c1c2_prop_df_v1_labs_survival_post1)
p <- ggsurvplot(c1c2_prop_df_v1_labs_survival_PFS, 
                size =3, 
                surv.median.line = "hv",
                title = "Progression Free Survival",
                xlab = "Time(Months)",
                legend.labs =
                  c(paste("High","NLR",sep=" "), paste("Low","NLR",sep=" ")),   
                pval.method = TRUE,
                data = c1c2_prop_df_v1_labs_survival_post1,
                pval = TRUE,risk.table = TRUE, risk.table.y.text.col = TRUE,
                palette = c("#B7657B","#3A3845"),
                fontsize = 8,
                risk.table.y.text = TRUE,
                font.title = c(30, "bold", "darkblue"),
                font.subtitle = c(30, "bold.italic", "purple"),
                font.caption = c(30, "plain", "orange"),
                font.x = c(30, "bold", "black"),
                font.y = c(30, "bold", "black"),
                font.main = c(30, "bold", "black"),
                font.tickslab = c(30, "plain", "black"),
                pval.size = 10,
                pval.coord = c(30,0.9),
                pval.method.coord = c(30,1)
)

p$plot <- p$plot +
  theme(legend.text = element_text(size = 20, color = "black", face = "bold"))+
  annotate("text", label = "High NLR", x =10, y = 0.1, size = 10, colour = color.high )+
  annotate("text", label = "Low NLR", x =20, y = 0.57, size = 10, colour = color.low)

p$table <- ggrisktable(c1c2_prop_df_v1_labs_survival_PFS,
                       data = c1c2_prop_df_v1_labs_survival_post1,
                       color = "by.NLR.median",
                       palette = c(color.high,color.low),
                       y.text = T,
                       fontsize=10,
                       y.text.col = TRUE, xlab = "",font.tickslab = c(20, "bold"),
                       tables.theme = theme_cleantable()
)+scale_y_discrete(labels=c("Low NLR","High NLR"))





#Figure2F ---------------------------------------------------------------------------------------
##OS
c1c2_prop_df_v1_labs_survival_OS <- survfit(Surv(OS/30,death) ~ by.NLR.median,
                                            data = c1c2_prop_df_v1_labs_survival_post1)

p <- ggsurvplot(c1c2_prop_df_v1_labs_survival_OS, 
                size =3, #size of the line (line type can be changed too)
                surv.median.line = "hv",
                title = "Overall Survival",
                xlab = "Time(Months)",
                legend.labs =
                  c(paste("High","NLR",sep=" "), paste("Low","NLR",sep=" ")),   
                pval.method = TRUE,
                data = c1c2_prop_df_v1_labs_survival_post1,
                pval = TRUE,risk.table = TRUE, risk.table.y.text.col = TRUE,
                palette = c("#B7657B","#3A3845"),
                fontsize = 8,
                risk.table.y.text = TRUE,
                font.title = c(30, "bold", "darkblue"),
                font.subtitle = c(30, "bold.italic", "purple"),
                font.caption = c(30, "plain", "orange"),
                font.x = c(30, "bold", "black"),
                font.y = c(30, "bold", "black"),
                font.main = c(30, "bold", "black"),
                font.tickslab = c(30, "plain", "black"),
                pval.size = 10,
                pval.coord = c(30,0.9),
                pval.method.coord = c(30,1)
)

p$plot <- p$plot +
  theme(legend.text = element_text(size = 20, color = "black", face = "bold"))+
  annotate("text", label = "High NLR", x =15, y = 0.4, size = 10, colour = color.high )+
  annotate("text", label = "Low NLR", x =22, y = 0.95, size = 10, colour = color.low)

p$table <- ggrisktable(c1c2_prop_df_v1_labs_survival_OS,
                       data = c1c2_prop_df_v1_labs_survival_post1,
                       color = "by.NLR.median",
                       palette = c(color.high,color.low),
                       y.text = T,
                       fontsize=10,
                       y.text.col = TRUE, xlab = "",font.tickslab = c(20, "bold"),
                       tables.theme = theme_cleantable()
)+scale_y_discrete(labels=c("Low NLR","High NLR"))



#Figure3A---------------------------------------------------------------------------------
triangle.hm.all <- Heatmap(tmp.scs.corr, rect_gp = gpar(type = "none"), 
                           name = "Correlation Coefficient (R)",
                           col = colorRamp2(c(-1,0,1),
                                            c("#D1D1D1", "white","#94B49F")),
                           row_names_side = "left",show_row_names = TRUE,
                           column_names_rot = 30,
                           show_column_names = TRUE,column_names_side = "bottom",
                           column_dend_side = "bottom",
                           height = unit(30,"cm"),
                           width = unit(30,"cm"),
                           cell_fun = function(j, i, x, y, w, h, fill) {
                             if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                               grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                               if(tmp.scs.corr.pval[i, j] < 0.05){
                                 grid.text(
                                   paste0(sprintf("%.2f",tmp.scs.corr[i, j])),
                                   x, y, gp = gpar(fontsize = 12))
                               }
                             }
                           }
                           
)




#Figure4A---------------------------------------------------------------------------------
luminex_c1c2_df_v0117_v2 %>%
  filter(grepl("pre",sample_id)) %>% 
  ggplot(
    aes(x=log2(IL10),y=log2(CXCL6), color=DCB, size=log2(`IL.23`))) +
  geom_point(alpha=0.7,label = FALSE) +
  labs(x="Log2(IL10)", y = "Log2(CXCL6)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  scale_colour_manual(values=c("#B09C8599","#3C548899"))+ geom_hline(yintercept=7.3, linetype="dashed", color = "black")+
  geom_vline(xintercept= round(median(log2(luminex_c1c2_df_v0117_v2$IL10)),2),linetype = "longdash")+
  annotate("text", label = paste0("median(log2(IL10)) \n =",
                                  round(median(log2(luminex_c1c2_df_v0117_v2$IL10)),2)), 
           x =round(median(log2(luminex_c1c2_df_v0117_v2$IL10)),2), y = 10, size = 10, colour = "black")+
  theme(text = element_text(size = 30),
        legend.key=element_blank())+ #this gets rid of the grey background in the legend 
  guides(color = guide_legend(override.aes = list(size = 10)))






#Figure4B---------------------------------------------------------------------------------
rocs_v0126 <- list()
rocs_v0126[["CXCL6"]] <- roc(df_v2_0126$DCB, predict(glm(as.factor(DCB) ~ CXCL6_cat, data = df_v2_0126, family =binomial),newdata=df_v2_0126, type = "response") ,smoothed = TRUE)
rocs_v0126[["CXCL6+IL10"]] <- roc(df_v2_0126$DCB,predict(glm(as.factor(DCB) ~ CXCL6_cat+IL10_cat, data = df_v2_0126, family =binomial),newdata=df_v2_0126, type = "response") ,smoothed = TRUE)
rocs_v0126[["CXCL6+IL10+IL23"]] <- roc(df_v2_0126$DCB, predict(glm(as.factor(DCB) ~ CXCL6_cat+IL10_cat+IL.23_cat, data = df_v2_0126, family =binomial),newdata=df_v2_0126, type = "response") ,smoothed = TRUE)

ggroc(rocs_v0126,legacy.axes = TRUE, 
      size = 1.5)+
  theme_bw()+
  geom_abline()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         title=element_blank(),
         axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("#6C4A4A","#C89595","#C7BEA2"))+
  theme(text = element_text(size = 30))+
  annotate("text", label = paste0("AUC =",round(rocs_v0126$CXCL6$auc,2)), x =0.75, y = 0.5, size = 10, colour = "#6C4A4A")+
  annotate("text", label = paste0("AUC =",round(rocs_v0126$`CXCL6+IL10`$auc,2)), x =0.75, y = 0.4, size = 10, colour = "#C89595")+
  annotate("text", label = paste0("AUC =",round(rocs_v0126$`CXCL6+IL10+IL23`$auc,2)), x =0.75, y = 0.3, size = 10, colour = "#C7BEA2")


glm.fit.v2 <- glm(as.factor(DCB) ~ CXCL6_cat+IL10_cat+IL.23_cat, 
                  data = df_v2, family =binomial)
od_ratio_table_v2 <- as.data.frame(exp(cbind("Odds ratio" = coef(glm.fit.v2), 
                                             confint.default(glm.fit.v2, level = 0.95)))) %>% 
  cbind("P values" = coef(summary(glm.fit.v2))[,'Pr(>|z|)']) %>%
  rownames_to_column("Covariate")
colnames(od_ratio_table_v2) <- c("Covariate","Odds_ratio","ciLow","ciHigh","P values")


od_ratio_table_gt_v2_v2 <- od_ratio_table_v2 %>%
  mutate(`P values` = formatC(`P values`, format = "e", digits = 2)) %>%
  dplyr::mutate(across(2:4, round, 2)) %>%
  filter(Covariate != "(Intercept)") %>%
  mutate(Covariate = paste("High",sub("_[^_]+$", "",Covariate))) %>%
  gt() %>% tab_options(table.font.size=25) 







#Figure4C---------------------------------------------------------------------------------
surv_object_ci <- Surv(surv_mat_v230309$OS/30, event = surv_mat_v230309$death)
surv_object_pfs_ci <- Surv(surv_mat_v230309$PFS/30, event = surv_mat_v230309$progress)

fit_ci <- survfit(surv_object_ci ~ strata, data = surv_mat_v230309)
surv_mat_v230309$strata <- as.factor(surv_mat_v230309$strata)

p <- ggsurvplot(fit_ci, 
                size =3, 
                legend.labs = levels(surv_mat_v230309$strata),
                title = "Overall Survival",
                xlab = "Time(Months)", 
                pval.method = TRUE,
                data = surv_mat_v230309, pval = TRUE,risk.table = TRUE, risk.table.y.text.col = TRUE,
                palette = rev(c("#181D31","#678983","#898121","#B2B1B9")),
                fontsize = 8,
                risk.table.y.text = TRUE,
                font.title = c(30, "bold", "darkblue"),
                font.subtitle = c(30, "bold.italic", "purple"),
                font.caption = c(30, "plain", "orange"),
                font.x = c(30, "bold", "black"),
                font.y = c(30, "bold", "black"),
                font.main = c(30, "bold", "black"),
                font.tickslab = c(30, "plain", "black"),
                pval.size = 10,
                pval.coord = c(36,0.9),
                pval.method.coord = c(36,1)
)

p$plot <- p$plot +
  theme(legend.text = element_text(size = 30, color = "black", face = "bold"))+
  annotate("text", label = "score3", 
           x =25, y = 0.9, size = 10, colour = "#181D31")+
  annotate("text", label = "score2", x =30, y = 0.75, size = 10, colour = "#678983")+
  annotate("text", label = "score1", x =18, y = 0.5, size = 10, colour = "#898121")+
  annotate("text", label = "score0", x =10, y = 0.4, size = 10, colour = "#B2B1B9")

class(levels(p$data.survtable$group))
p$table <- ggrisktable(fit_ci,
                       data = surv_mat_v230309,
                       color = "strata",
                       palette = rev(c("#181D31","#678983","#898121","#B2B1B9")), #the table annotation goes from top to bottom
                       y.text = T,
                       fontsize=10,
                       y.text.col = TRUE, xlab = "",font.tickslab = c(20, "bold"),
                       tables.theme = theme_cleantable())+
  scale_y_discrete(labels=rev(levels(surv_mat_v230309$strata))) #the labeling in the table goes from bottom to top

fit_pfs_ci <- survfit(surv_object_pfs_ci ~ strata, data = surv_mat_v230309)
surv_mat_v230309$strata <- as.factor(surv_mat_v230309$strata)
p <- ggsurvplot(fit_pfs_ci, 
                size =3, 
                legend.labs = levels(surv_mat_v230309$strata),
                title = "Progression Free Survival",
                xlab = "Time(Months)", 
                pval.method = TRUE,
                data = surv_mat_v230309, pval = TRUE,risk.table = TRUE, risk.table.y.text.col = TRUE,
                palette = rev(c("#181D31","#678983","#898121","#B2B1B9")),
                fontsize = 8,
                risk.table.y.text = TRUE,
                font.title = c(30, "bold", "darkblue"),
                font.subtitle = c(30, "bold.italic", "purple"),
                font.caption = c(30, "plain", "orange"),
                font.x = c(30, "bold", "black"),
                font.y = c(30, "bold", "black"),
                font.main = c(30, "bold", "black"),
                font.tickslab = c(30, "plain", "black"),
                pval.size = 10,
                pval.coord = c(36,0.9),
                pval.method.coord = c(36,1)
)


p$plot <- p$plot +
  theme(legend.text = element_text(size = 30, color = "black", face = "bold"))+
  annotate("text", label = "score3", 
           x =25, y = 0.7, size = 10, colour = "#181D31")+
  annotate("text", label = "score2", x =20, y = 0.55, size = 10, colour = "#678983")+
  annotate("text", label = "score1", x =10, y = 0.35, size = 10, colour = "#898121")+
  annotate("text", label = "score0", x =5, y = 0.2, size = 10, colour = "#B2B1B9")

class(levels(p$data.survtable$group))
p$table <- ggrisktable(fit_ci,
                       data = surv_mat_v230309,
                       color = "strata",
                       palette = rev(c("#181D31","#678983","#898121","#B2B1B9")),
                       fontsize=10,
                       y.text.col = TRUE, xlab = "",font.tickslab = c(20, "bold"),
                       tables.theme = theme_cleantable())+
  scale_y_discrete(labels=rev(levels(surv_mat_v230309$strata))) 




#Figure4D---------------------------------------------------------------------------------
ggforest3(surv_mat_v230309_model_pfs, font.x.size = 20,fontsize = 2,
          main = "Hazard ratio (PFS)")

ggforest3(surv_mat_v230309_model_pfs,
          font.x.size = 20,fontsize = 2,
          title.size = 30, 
          main = "Hazard ratio (PFS)")

surv_mat_score1to2$`feature score` <- as.character(surv_mat_score1to2$`feature score`)
table(surv_mat_score1to2$`feature score`)
surv_mat_score1to2_model_os <- coxph(Surv(surv_mat_score1to2$OS/30, 
                                          event = surv_mat_score1to2$death) ~ `feature score`,
                                     data = surv_mat_score1to2)
ggforest3(surv_mat_score1to2_model_os,
          font.x.size = 20,fontsize = 2,
          title.size = 30, 
          main = "Hazard ratio (OS)")


#Figure5A---------------------------------------------------------------------------------
ggplot() + geom_bar(aes(y = Proportion,
                        x=priorIO,
                        fill = DCB),
                    data = prop.m,
                    stat = "identity",
                    width = 0.5)+
  xlab("priorIO")+
  ylab("DCB")+
  labs(title="priorIO vs. DCB") +
  theme_classic()+
  scale_fill_manual(values = c("DCB_Y" = DCB.color, "DCB_N" = NCB.color))+
  theme(text = element_text(size = 30),
        legend.title = element_blank())+
  geom_text(aes(x=1,y=0.2),size = 12,data = prop.m,label = addmargins(c1.c2.DCB.prior)[2,1])+
  geom_text(aes(x=1,y=0.8),size = 12,data = prop.m,label = addmargins(c1.c2.DCB.prior)[1,1])+
  geom_text(aes(x=2,y=0.2),size = 12,data = prop.m,label = addmargins(c1.c2.DCB.prior)[2,2])+
  geom_text(aes(x=2,y=0.8),size = 12,data = prop.m,label = addmargins(c1.c2.DCB.prior)[1,2])


#Figure5B---------------------------------------------------------------------------------
cytokines.pre.for.priorIO.hm.t <- Heatmap(
  matrix = tmp.pre.pval.for.priorIO.reci.mtx.t,
  name= "ad_cd3_func_pre_YvsN",
  cluster_columns = FALSE,cluster_rows = FALSE,
  col = colorRamp2(c(-20,0,20),
                   c("#CEB5B9","white","#98A1B1")),
  column_title_side = ifelse(length(by) ==2,"top", "bottom"),
  show_row_dend = FALSE,
  show_column_dend = FALSE, 
  row_names_side = "left",
  column_names_rot =  60,column_names_side = "top",
  row_dend_width = unit(2, "cm"),
  rect_gp = gpar(col = "black"),
  width = unit(12,"cm"),
  height=unit(1,"cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(between(abs(tmp.pre.pval.for.priorIO.reci.mtx.t[i, j]),1/0.05,1/0.01)==TRUE)
      grid.text(sprintf("*",tmp.pre.pval.for.priorIO.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (between(abs(tmp.pre.pval.for.priorIO.reci.mtx.t[i, j]),1/0.01,1/0.001)==TRUE)
      grid.text(sprintf("**",tmp.pre.pval.for.priorIO.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (between(abs(tmp.pre.pval.for.priorIO.reci.mtx.t[i, j]),1/0.001,1/0.0001)==TRUE)
      grid.text(sprintf("***",tmp.pre.pval.for.priorIO.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
    else if (abs(tmp.pre.pval.for.priorIO.reci.mtx.t[i, j]) >= 1/0.0001)
      grid.text(sprintf("****",tmp.pre.pval.for.priorIO.reci.mtx.t[i, j]),
                x, y, gp = gpar(fontsize = 20))
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    grid_height = unit(1, "cm"), grid_width = unit(1, "mm"),direction = "vertical",
    at = c(-20,-10,10,20),
    labels = c("p <= 0.05","p <= 0.1",
               "p <=0.1","p <=0.05"),
    title = "p value"
  )
)




#Figure5C---------------------------------------------------------------------------------
roc_list_io_exp_naive_v2 <- list()
roc_list_io_exp_naive_v2[["IO naive"]] <- roc(surv_mat_v230323_io_naive$DCB, 
                                              predict(glm(as.factor(DCB) ~ `CXCL6_cat`+`IL10_cat`+`IL.23_cat`, 
                                                          data = df_v2_0126_cytokines_only_na_omit, 
                                                          family =binomial),
                                                      newdata=surv_mat_v230323_io_naive, type = "response"),
                                              smoothed = TRUE)
roc_list_io_exp_naive_v2[["IO experienced"]] <- roc(surv_mat_v230323_io_exp$DCB, 
                                                    predict(glm(as.factor(DCB) ~ `CXCL6_cat`+`IL10_cat`+`IL.23_cat`, 
                                                                data = df_v2_0126_cytokines_only_na_omit, 
                                                                family =binomial),
                                                            newdata=surv_mat_v230323_io_exp, type = "response"),
                                                    smoothed = TRUE)

ggroc(roc_list_io_exp_naive_v2,legacy.axes = TRUE, 
      size = 1.5)+
  xlab("False positive rate") + ylab("True positive rate")+
  theme_bw()+
  geom_abline()+
  scale_color_manual(values = c(priorIO.N.color,priorIO.Y.color))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         legend.title= element_blank(),
         axis.line = element_line(colour = "black"),
         text = element_text(size = 30))+
  ggtitle("AUC of Three-feature models \n (CXCL6, IL-23, IL-10)")+
  annotate("text", 
           label = paste0("AUC =",round(roc_list_io_exp_naive_v2$`IO naive`$auc,2)), 
           x =0.75, y = 0.6, size = 10, colour = priorIO.N.color)+
  annotate("text", 
           label = paste0("AUC =",round(roc_list_io_exp_naive_v2$`IO experienced`$auc,2)),
           x =0.75, y = 0.5, size = 10, colour = priorIO.Y.color)


#Figure5D---------------------------------------------------------------------------------

#CXCL6
luminex_c1c2_df_v0117_v2 %>%
  mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  mutate(CXCL6 = log2(CXCL6)) %>%
  dplyr::select(CXCL6,DCB,timepoint,priorIO) %>%
  mutate(priorIO = if_else(priorIO == "Y","with priorIO","no priorIO")) %>%
  filter(timepoint == "pre") %>%
  ggplot(aes(x=factor(priorIO,levels = c("with priorIO","no priorIO")),
             y=CXCL6, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="", y = "Log2(CXCL6)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",hide.ns = TRUE,size = 20,
                     #size=8,
                     label.y.npc = 0.85,
                     label.x.npc = 0.5,
                     show.legend = FALSE
  )+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))



#Figure5E---------------------------------------------------------------------------------
#IL10
luminex_c1c2_df_v0117_v2 %>%
  mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  mutate(IL10 = log2(IL10)) %>%
  dplyr::select(IL10,DCB,timepoint,priorIO) %>%
  mutate(priorIO = if_else(priorIO == "Y","with priorIO","no priorIO")) %>%
  filter(timepoint == "pre") %>%
  ggplot(aes(x=factor(priorIO,levels = c("with priorIO","no priorIO")),
             y=IL10, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="", y = "Log2(IL10)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",hide.ns = TRUE,size = 20,
                     label.y.npc = 0.85,
                     label.x.npc = 0.5,
                     show.legend = FALSE
  )+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))

#Figure5F---------------------------------------------------------------------------------
luminex_c1c2_df_v0117_v2 %>%
  mutate(DCB = factor(DCB, levels = c("Y","N"), ordered = TRUE)) %>%
  mutate(IL23 = log2(IL.23)) %>%
  dplyr::select(IL23,DCB,timepoint,priorIO) %>%
  mutate(priorIO = if_else(priorIO == "Y","with priorIO","no priorIO")) %>%
  filter(timepoint == "pre") %>%
  ggplot(aes(x=factor(priorIO,levels = c("with priorIO","no priorIO")),
             y=IL23, color=DCB)) +
  geom_jitter( position=position_jitterdodge(), size=3, alpha=0.9) +
  labs(x="", y = "Log2(IL23)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",hide.ns = TRUE,size = 20,
                     label.y.npc = 0.85,
                     label.x.npc = 0.5,
                     show.legend = FALSE
  )+
  stat_summary(fun.data =median_hilow,  conf.int=.5, position=position_dodge(.9), geom = "pointrange", aes(group=(DCB)), color="black", width = 0.5)+
  stat_summary(fun.y = median, position=position_dodge(.9), geom = "crossbar", aes(group=DCB), color="black", width = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual("DCB",values = c("Y" = DCB.color,"N" = NCB.color))+
  theme(text = element_text(size = 30))
