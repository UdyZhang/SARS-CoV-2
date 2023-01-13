library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggsignif)
library(dplyr)
ibrary(ComplexHeatmap)
library(ggridges)
library(readxl)
library(Biostrings)
library(seqinr)

#read data
fn <- list.files("./", pattern = "*.csv", full.names = T)
file_tab2 <- read.table("data_sec_filetab.txt", header = T, sep = "\t", stringsAsFactors = F)
vgenelen_tab <- read.table("Vgene_length.txt", sep = "\t", header = T, stringsAsFactors = F)

data_vdj_raw2 <- list()
data_vdj_single2 <- list()
for (i in 1:11) {
  data_vdj_raw2[[i]] <- read.csv(fn[i])
  h <- data_vdj_raw2[[i]][, -c(2:5,11,12,29:31)]
  n <- table(h$barcode) %>% data.frame()
  h <- h[which(h$barcode %in% n$Var1[which(n$Freq == 2)]), ]
  n <- h[which(h$chain == "IGH"), ]
  h <- h[-which(h$chain == "IGH"), ]
  names(n)[-1] <- paste0(names(n)[-1], "_H")
  names(h)[-1] <- paste0(names(h)[-1], "_L")
  data_vdj_single2[[i]] <- merge(n, h, by = c("barcode"), all = T)
  data_vdj_single2[[i]]$cell_id <- paste0(file_tab2$sample[i], "_", data_vdj_single2[[i]]$barcode)
}
names(data_vdj_single2) <- file_tab2$sample

#filter functional ab
pos <- which(data_vdj_single_tab$v_gene_H %in% Hvgene_f$vgene & data_vdj_single_tab$v_gene_L %in% Lvgene_f$vgene)
data_vdj_single_tab <- data_vdj_single_tab[pos,]
pos <- intersect(grep(pattern = "^C([A-Z]+)W$", data_vdj_single_tab$cdr3_H),
                 grep(pattern = "^C([A-Z]+)F$", data_vdj_single_tab$cdr3_L)) #CDR3
data_vdj_single_tab <- data_vdj_single_tab[pos,]


#SHM
##all
n <- data_vdj_shm_aa[which(data_vdj_shm_aa$group == "Delta-\nInfected"), ]
m <- data_vdj_shm_aa[which(data_vdj_shm_aa$group %in% c("Non-\nvaccinated", "Known mAb")), ]
n <- rbind(n, m)
t <- unique(n$group) %>% as.character()
n$group <- factor(n$group, levels = t[c(3:1,4)])

par(mai = c(1,2,.1,.1))
ggplot(n, aes(x = SHM_aa_H, y = group, fill = group)) +
  geom_density_ridges(alpha = 0.6, from = 0, to = 40, quantile_lines = T, 
                      quantile_fun = function(SHM_aa_H,...)mean(SHM_aa_H),
                      color = "white", vline_color = "black", vline_linetype = "dashed",
                      scale = 1.7) +
  theme_ridges(grid = F, center = T) +
  scale_y_discrete(labels = rev(c("Delta-\nInfected", 
                                  "Non-\nvaccinated", "Known\nmAb")),
                   expand = c(0.02, 0.02)) +
  scale_fill_manual(values = c(rev(col_group), "#FF7F00")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.spacing.x = unit(8, "mm"),
        axis.line.x.bottom = element_line(colour = "black", size = 0.5),
        axis.text.y.left = element_text(size = 12, family = "Arial", hjust = 0.5),
        axis.text.x.bottom = element_text(size = 10, family = "Arial"),
        plot.title = element_text(hjust = 0.5, size = 14, family = "Arial", face = "bold")) +
  labs(title = "Number of somatic hypermutation\nin VH (amino acid)", y = NULL, x = NULL) +
  guides(fill = "none")
par(opar)

##low/high SHM proportion
m <- mean(data_vdj_shm_aa$SHM_aa_H[which(data_vdj_shm_aa$group == "Known mAb")])
n <- data_vdj_shm_aa[, c(1,2,4)]
n <- n[which(n$group %in% unique(n$group)[2:3]),]
n$group <- as.character(n$group)
n <- merge(n, data.frame(mAb = data_vdj_single_tab$cell_id,
                         sample = data_vdj_single_tab$ind,
                         stringsAsFactors = F), by = "mAb", all = F, sort = F)
t <- unique(n$sample)[c(7,10,11,12,15:20)]#23
nt <- list()
for (i in 1:length(t)) {
  h <- n[which(n$sample == t[i]),]
  h$shm_level <- "High"
  h$shm_level[which(h$SHM_aa_H == 0)] <- "Low"
  m <- table(h[, -c(1,2)]) %>% prop.table() %>% data.frame()
  m$Freq <- m$Freq * 100
  nt[[i]] <- m
}

h <- do.call(rbind, nt)

ggplot(h[which(h$shm_level == "Low"),], aes(x = group, y = Freq)) +
  stat_boxplot(geom = "errorbar", aes(color = group), width = 0.2) +
  geom_boxplot(aes(color = group), fill = "white", outlier.shape = NA, width = 0.5) +
  geom_jitter(size = 2, aes(fill = group), shape = 21, alpha = 0.6, width = 0.2) +
  geom_signif(comparisons = list(c("Delta Breakthrough", "Non-vaccinated")), 
              test = "t.test", step_increase = .05, map_signif_level = T, 
              y_position = 16) +
  scale_fill_manual(values = col_group[1:2]) +
  scale_color_manual(values = col_group[1:2]) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 6)) +
  scale_x_discrete(labels = c("Delta\nBreakthrough", "Non-\nvaccinated")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "grey90"),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.left = element_text(vjust = 2, size = 12, family = "Arial", face = "bold"), 
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = 10, family = "Arial", color = "black")) +
  labs(x = NULL, y = "Proportion of unmutated VH (%)") +
  guides(fill = "none", color = "none")
 
##selected & validated
n <- n[which(n$group != "FALSE"), c(2,5,15)]
n$single <- n$SHM_aa_H
m <- table(n$vgene_h)
n$single[which(n$vgene_h %in% names(m)[which(m > 2)])] <- NA
t <- data_vdj_shm_aa[which(data_vdj_shm_aa$group == "Known mAb"),]
t$vgene_h <- multi_gsub(t$mAb, ab_db$Name, ab_db$`Heavy V Gene`)
t$single <- NA
t <- t[which(t$vgene_h %in% unique(n$vgene_h)), ]
n <- rbind(n, t[,c(5,2,4,6)])
n$source <- "Binding mAb"
n$source[which(n$group == "Known mAb")] <- "Known mAb"
n$source <- factor(n$source, levels = unique(n$source)[2:1])
n$vgene_h <- multi_gsub(n$vgene_h, c("IGHV1-69D", "IGHV3-23D", "IGHV3-30-3"),
                        c("IGHV1-69", "IGHV3-23", "IGHV3-30"))
n$vgene_h <- factor(n$vgene_h, levels = rev(Hvgene_f$vgene))
n$group <- factor(n$group, levels = unique(n$group)[c(2,1,3)])

ggplot(n, aes(y = vgene_h)) +
  geom_density_ridges(aes(x = SHM_aa_H, fill = source),
                      alpha = 0.6, from = 0, to = 40, quantile_lines = T, quantiles = 2,
                      color = "white", vline_color = "black", vline_linetype = "dashed",
                      scale = 1.7) +
  theme_ridges(grid = F, center = T) +
  geom_point(aes(x = single), shape = 21, fill = "#E41A1C", alpha = 0.6, size = 1, color = "white") +
  scale_fill_manual(values = c("#4DAF4A", "#E41A1C")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.spacing.x = unit(8, "mm"),
        axis.line.x.bottom = element_line(colour = "black", size = 0.5),
        axis.text.y.left = element_text(size = 10, family = "Arial", hjust = 1),
        axis.text.x.bottom = element_text(size = 10, family = "Arial"),
        axis.title.x.bottom = element_text(size = 14, family = "Arial", face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = "Number of somatic hypermutation\n(amino acid on VH)", y = NULL)

n$high <- "T"
n$high[which(n$SHM_aa_H <= 5)] <- "F"
df <- unique(n[, c(1,5)])
df$high_prop <- NA
for (i in 1:nrow(df)) {
  m <- n[which(n$vgene_h == df$vgene_h[i] & n$source == df$source[i]),]
  df$high_prop[i] <- length(which(m$high == "T")) / nrow(m) * 100
}

ggplot(df, aes(x = high_prop, y = vgene_h, group = vgene_h)) +
  geom_line() +
  geom_point(aes(fill = source), shape = 21, size = 3) +
  scale_fill_manual(values = c("#4DAF4A", "#E41A1C")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", fill = NA, size = .5),
        axis.line.x.bottom = element_line(colour = "black", size = 0.5),
        axis.line.y.left = element_line(colour = "black", size = 0.5),
        axis.text.y.left = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x.bottom = element_text(size = 10, family = "Arial", color = "black"),
        axis.title.x.bottom = element_text(vjust = 0, size = 12, family = "Arial", face = "bold"),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white")) +
  labs(y = NULL, x = "Proportion of mAbs with\nmore than 5 SHM on VH (%)") +
  guides(fill = guide_legend(title = ""))

##cross-variant ~ known
h <- data_vdj_shm_aa[which(data_vdj_shm_aa$mAb %in% ab_bind_tab$mAb[2:23]), ]
h$group <- "Cross-variant binding"
n <- data_vdj_shm_aa[which(data_vdj_shm_aa$mAb %in% ab_selc$mAb[which(ab_selc$group == "TRUE")]),]
n$group <- "Binding mAb"
h <- rbind(h, n)
m <- data_vdj_shm_aa[which(data_vdj_shm_aa$group %in% c("Known mAb")), ]
n <- rbind(m, h)
t <- unique(n$group) %>% as.character()
n$group <- factor(n$group, levels = t[c(2,3,1)])

signif_df <- data.frame(start = t[c(2,3,2)],
                        end = t[c(3,1,1)],
                        p = NA, y = NA, stringsAsFactors = F)

for (i in 1:3) {
  x <- n$SHM_aa_H[which(n$group == signif_df[i, 1])]
  y <- n$SHM_aa_H[which(n$group == signif_df[i, 2])]
  signif_df$p[i] <- wilcox.test(x,y)$p.value
}
signif_df$y <- c(37,35)

ggplot(n[which(n$SHM_aa_H < 40), ], aes(y = SHM_aa_H, x = group)) +
  geom_violin(color = "black", alpha = 0.6, aes(fill = group)) +
  stat_summary(fun = mean, geom="point") +
  stat_summary(fun.data = mean_sdl,fun.args = list(mult=1), geom='errorbar',width=.2, color = "black") + 
  geom_signif(data = signif_df,
              aes(xmin = start, xmax = end, annotations = label, y_position = y),
              manual = T, vjust = 0.5, tip_length = 0.01) +
  scale_fill_manual(values = c("red", "grey60", col_group[3])) +
  scale_y_continuous(limits = c(0, 40), breaks = c(0,2,5,10,20,30,40), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Cross-variant\nbinding", "Binding\nmAb", "Known\nmAb")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", fill = NA, size = .5),
        axis.line.x.bottom = element_line(colour = "black", size = 0.5),
        axis.line.y.left = element_line(colour = "black", size = 0.5),
        axis.text.y.left = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x.bottom = element_text(size = 10, family = "Arial", color = "black"),
        axis.title.y.left = element_text(vjust = 2, size = 12, family = "Arial", face = "bold"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, vjust = 0.5, family = "Arial", face = "bold")) +
  labs(y = "Number of somatic hypermutation\nin VH (amino acid)", x = NULL) +
  guides(fill = "none")

#V gene
##selected & validated
n <- ab_selc[, c(2,15)]
n$vgene_h <- multi_gsub(n$vgene_h, c("IGHV1-69D", "IGHV3-23D", "IGHV3-30-3"),
                        c("IGHV1-69", "IGHV3-23", "IGHV3-30"))
n$vgene_h <- factor(n$vgene_h, levels = rev(Hvgene_f$vgene))
n$group <- factor(n$group, levels = unique(n$group)[2:1])

ggplot(n, aes(x = vgene_h, fill = group)) +
  geom_bar(color = "black", alpha = 0.6) +
  scale_fill_manual(values = c("white", "#E41A1C"), 
                    labels = c("No binding (54)", "Binding (63)")) +
  scale_y_continuous(limits = c(0, 25), expand = c(0,0), breaks = seq(0, 25, 5)) +
  coord_flip() +
  theme(panel.background = element_rect(fill = "white"),
        panel.spacing = unit(.3, "lines"),
        panel.border = element_blank(),
        axis.line.x.bottom = element_line(colour = "black", size = 0.5),
        axis.line.y.left = element_line(colour = "black", size = 0.5),
        axis.text.y.left = element_text(size = 10, family = "Arial", vjust = 0.5, color = "black"),
        axis.text.x.bottom = element_text(size = 10, family = "Arial", vjust = 0.5, 
                                          hjust = 0.5, color = "black"),
        axis.title.x.bottom = element_text(size = 14, family = "Arial", face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold"),
        plot.margin = margin(.3, .3, .3, .3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizon",
        legend.text = element_text(size = 10, family = "Arial", color = "black")) +
  labs(x = NULL, y = "Number") +
  guides(fill = guide_legend(title = "117 selected mAbs", title.position = "top",
                             title.hjust = 0.5))

#isotype
##boxplot
n <- data_vdj_shm_aa[, c(1,4)]
n <- n[which(n$group %in% unique(n$group)[2:3]),]
n$ind <- apply(n, 1, function(x){strsplit(x[1], "_")[[1]][1]})
n <- n[which(n$ind %in% unique(n$ind)[c(7,10:12,15:20)]),]
n$group <- as.character(n$group)
n <- merge(n, data.frame(mAb = data_vdj_single_tab$cell_id, 
                         isotype = data_vdj_single_tab$c_gene_H,
                         sample = data_vdj_single_tab$ind,
                         stringsAsFactors = F), by = "mAb", all = F, sort = F)
n <- n[which(n$isotype %in% isotype),]
t <- unique(n$sample)
nt <- list()
for (i in 1:length(t)) {
  m <- table(n[which(n$sample == t[i]),-1]) %>% prop.table() %>% data.frame()
  m$Freq <- m$Freq * 100
  m <- arrange(m, isotype)
  #m$label <- paste0(m$isotype, "\n", round(m$Freq, 2), "%")
  nt[[i]] <- m
}

h <- do.call(rbind, nt)
h$isotype <- factor(h$isotype, levels = isotype[c(2,1,3:8)])
h$group <- gsub("Infected Only", "Non-vaccinated", h$group)

signif_df <- data.frame(isotype = rep(isotype),
                        start = rep("Delta Breakthrough", 8),
                        end = rep("Non-vaccinated", 8),
                        p = NA, y = 0, stringsAsFactors = F)
signif_df$isotype <- factor(signif_df$isotype, levels = isotype)
for (j in 1:8) {
  x <- h$Freq[which(h$group == signif_df[j, 2] & h$isotype == signif_df[j, 1])]
  y <- h$Freq[which(h$group == signif_df[j, 3] & h$isotype == signif_df[j, 1])]
  signif_df$p[j] <- t.test(x,y)$p.value
  signif_df$y[j] <- max(c(x,y)) + 5
}

ggplot(h, aes(x = group, y = Freq)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  geom_jitter(size = 1.5, aes(fill = group), shape = 21, alpha = 0.6, width = 0.2) +
  geom_signif(data = signif_df, aes(xmin = start, xmax = end, annotations = label, y_position = y), manual = T, vjust = 0.5, tip_length = 0.01) +
  scale_fill_manual(values = col_group[1:2], labels = c("Delta-Infected", "Non-vaccinated")) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, 25)) +
  theme_bw() +
  labs(x = NULL, y = "Isotype Proportion (%)") +
  guides() +
  facet_grid(. ~ isotype)

#polar bar
n <- data_vdj_shm_aa[, c(1,4)]
n <- n[which(n$group %in% unique(n$group)[2:3]),]
n$ind <- apply(n, 1, function(x){strsplit(x[1], "_")[[1]][1]})
n <- n[which(n$ind %in% unique(n$ind)[c(7,10:12,15:20)]),]
n$group <- as.character(n$group)
n <- merge(n, data.frame(mAb = data_vdj_single_tab$cell_id, 
                         isotype = data_vdj_single_tab$c_gene_H,
                         stringsAsFactors = F), by = "mAb", all = F, sort = F)
n <- n[which(n$isotype %in% isotype),]
t <- unique(n$group)
m <- rbind(table(n[which(n$group == t[1]),-c(1,3)]) %>% prop.table() %>% data.frame(),
           table(n[which(n$group == t[2]),-c(1,3)]) %>% prop.table() %>% data.frame())
m$Freq <- m$Freq * 100
m$fold <- NA
m$fold[1:8] <- m$Freq[1:8]/m$Freq[9:16]
m$fold[9:16] <- m$Freq[9:16]/m$Freq[1:8]
m$fold <- log(m$fold)
m$isotype <- factor(m$isotype, levels = isotype)
m <- arrange(m, isotype)
m$label <- paste0(m$isotype, "\n", round(m$Freq, 2), "%")

p <- list()
for (i in 1:2) {
  h <- m[which(m$group == t[i]), ]
  h$label <- factor(h$label, levels = h$label[c(3:8,1:2)])
  myAngle <- seq(-20,-340,length.out = 8)
  p[[i]] <- ggplot(h) +
    geom_bar(aes(x = label, y = Freq, fill = isotype),
             stat = "identity", width = 1, color = "black") +
    coord_polar(theta = "x",start=0) +
    scale_fill_manual(values = brewer.pal(8,"Paired")[c(2,1,6,8,5,7,4,3)]) +
    ylim(c(0,55))+
    theme_light()+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_line(colour = "grey80",size=.25),
          axis.text.x = element_text(size = 7,colour="black",angle = myAngle,
                                     hjust = .5, family = "Arial",
                                     color = c( "red", rep("black",5), "blue", "black")),
          axis.text.y.left = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 10, family = "Arial", face = "bold", 
                                    hjust = .5, vjust = 2),
          legend.position = "right",
          legend.title = element_text(family = "Arial", hjust = 0.5, size = 8),
          legend.text = element_text(family = "Arial", size = 7)) +
    guides(fill = guide_colourbar(title = "logFC", barwidth = 0.6, barheight = 4)) +
    labs(title = c("Delta-Infected", "Non-vaccinated")[i])
}
ggarrange(plotlist = p, ncol = 2)

#22 binding mAbs
ab_bind_tab <- read.table("mAb_binding_tab.txt", sep = "\t", header = T, stringsAsFactors = F)
ab_bind_tab <- merge(ab_bind_tab, ab_validated[,1:2], by = "mAb", all = F, sort = F)
m <- ab_bind_tab[-c(1,24),]
mat <- as.matrix(m[,-c(1,13)])
rownames(mat) <- m$mAb
mat <- mat * 100
mat <- t(mat)
col_fun <- circlize::colorRamp2(c(0, 50, 100,200), c("white", "#FFD92F", "#FB8072", "#984EA3"))
col_row <- rep("black", 22)
col_row[c(3,6,9,14,20,22)] <- "red"
ht <- Heatmap(mat, show_row_dend = F, show_column_dend = F, show_heatmap_legend = T, 
              col = col_fun, column_order = colnames(mat), row_order = rownames(mat),
              row_names_side = "right", column_names_side = "top", 
              column_names_rot = 45, 
              border = "black", 
              rect_gp = gpar(col = "white", lwd = 1.5),
              row_names_gp = gpar(family = "Arial"),
              row_labels = c(gsub("[.]", " ", rownames(mat)[1:8]),
                             "Omicron BA.1 RBD", "Omicron BA.1 S1+S2", "SARS-CoV-1 S1"),
              column_names_gp = gpar(family = "Arial", col = "red"),
              heatmap_legend_param = list(title = "Relative binding (% of WT S1)",
                                          family = "Arial", border = "black",
                                          #title_position = "leftcenter-rot",
                                          legend_width = unit(5, "cm"),
                                          direction = "horizontal",
                                          title_gp = gpar(fontsize = 12, face = "bold")))
draw(ht, heatmap_legend_side = "bottom")

