#### Fig1 I-M

### Fig1 I-J

library(ggplot2)
library(rstatix)
library(ggh4x)
library(dplyr)
library(ggpubr)

##F/B
metadata = read.table("metadata.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
data = read.table("sum_p.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)[c("Firmicutes","Bacteroidetes"),] %>%
       t() %>% as.data.frame() %>% mutate(Group=row.names(.))
dat = inner_join(metadata,data,by = c("Description"="Group")) %>% mutate(group=paste0("Mouth_",gsub(".*_","",Group)))

dat$group <- factor(dat$group,levels = c(paste0("Mouth_",c(1,3,6,9))))
colnames(dat)[c(4,5)] = c("Firmicutes","Bacteroidetes")
head(dat)

lower_quantile <- quantile(dat$Firmicutes, 0.25)
upper_quantile <- quantile(dat$Firmicutes, 0.75)
lower_s <- min(dat$Firmicutes)
upper_s <- max(dat$Firmicutes)
df_sig <- dat %>% group_by(group) %>% wilcox_test(Firmicutes ~ Group2) %>% rstatix::adjust_pvalue() %>% add_significance("p.adj") %>% add_xy_position(x = "group")

p=ggplot(dat,aes(group, Firmicutes))+
  geom_rect(xmin =0.5,xmax=1.5,ymin=lower_s,ymax=upper_s,fill="#3B9AB2")+
  geom_rect(xmin =1.5,xmax=2.5,ymin =lower_s,ymax= upper_s,fill="#7294D4")+
  geom_rect(xmin =2.5,xmax=3.5,ymin=lower_s,ymax= upper_s, fill="#9E8CC6")+
  geom_rect(xmin =3.5,xmax=4.5,ymin =lower_s,ymax= upper_s,fill="#E6A0C4")+
  geom_linerange(xmin=0.5,xmax=4.5,y=upper_quantile,lty=2,color ="grey70",linewidth=0.4)+
  geom_linerange(xmin=0.5,xmax=4.5,y= lower_quantile,lty=2,color ="grey70",linewidth=0.4)+
  geom_linerange(x=0.5,ymin =lower_s,ymax=upper_s,lty=2,color = "grey70",linewidth=0.4)+
  geom_linerange(x=1.5,ymin =lower_s,ymax=upper_s,lty=2,color = "grey70",linewidth=0.4)+
  geom_linerange(x=2.5,ymin =lower_s,ymax=upper_s,lty=2,color = "grey70",linewidth=0.4)+
  geom_linerange(x=3.5,ymin =lower_s,ymax=upper_s,lty=2,color = "grey70",linewidth=0.4)+
  geom_linerange(x=4.5,ymin =lower_s,ymax=upper_s,lty=2,color ="grey70",linewidth=0.4)+
  geom_boxplot(aes(fill=Group2,color =Group2),position = position_dodge(0.8),outlier.color =NA,linewidth =1,width =0.5)+
  stat_summary(aes(fill=Group2),fun = mean, geom ="point", size=2,color ="white",show.legend =F,position = position_dodge(0.8))+
  geom_text_aimed(data=df_sig,aes(x=group,y= (upper_s-10), label = p.adj.signif, group = group),angle=90,size=5)+
  theme_void()+coord_polar()+theme(legend.position=c(0.9,0.15))+
  geom_rect(xmin =0.5,xmax=1.5,ymin=(upper_s+1),ymax=(upper_s+10),fill="#3B9AB2")+
  geom_rect(xmin =1.5,xmax=2.5,ymin=(upper_s+1),ymax=(upper_s+10),fill="#7294D4")+
  geom_rect(xmin =2.5,xmax=3.5,ymin=(upper_s+1),ymax=(upper_s+10), fill="#9E8CC6")+
  geom_rect(xmin =3.5,xmax=4.5,ymin=(upper_s+1),ymax=(upper_s+10),fill="#E6A0C4")+
  geom_text(data=df_sig, aes(x=x,y=(upper_s+6),label=group),size=4, angle=c(315,45,-45,45),fontface="bold",inherit.aes =F)+
  scale_y_continuous(limits=c(-20,(upper_s+10))) + labs(title = "Firmicutes")+
  theme(plot.title = element_text(hjust = 0.5 ,vjust=-5,face="bold",size = 16),legend.title = element_blank())
ggsave("./fig1_I_1.pdf",p,width = 6,height=6)


lower_quantile <- quantile(dat$Bacteroidetes, 0.25)
upper_quantile <- quantile(dat$Bacteroidetes, 0.75)
lower_s <- min(dat$Bacteroidetes)
upper_s <- max(dat$Bacteroidetes)
df_sig <- dat %>% group_by(group) %>% wilcox_test(Bacteroidetes ~ Group2) %>% rstatix::adjust_pvalue() %>% add_significance("p.adj") %>% add_xy_position(x = "group")

p=ggplot(dat,aes(group, Bacteroidetes))+
  geom_rect(xmin =0.5,xmax=1.5,ymin=lower_s,ymax=upper_s,fill="#3B9AB2")+
  geom_rect(xmin =1.5,xmax=2.5,ymin =lower_s,ymax= upper_s,fill="#7294D4")+
  geom_rect(xmin =2.5,xmax=3.5,ymin=lower_s,ymax= upper_s, fill="#9E8CC6")+
  geom_rect(xmin =3.5,xmax=4.5,ymin =lower_s,ymax= upper_s,fill="#E6A0C4")+
  geom_linerange(xmin=0.5,xmax=4.5,y=upper_quantile,lty=2,color ="grey70",linewidth=0.4)+
  geom_linerange(xmin=0.5,xmax=4.5,y= lower_quantile,lty=2,color ="grey70",linewidth=0.4)+
  geom_linerange(x=0.5,ymin =lower_s,ymax=upper_s,lty=2,color = "grey70",linewidth=0.4)+
  geom_linerange(x=1.5,ymin =lower_s,ymax=upper_s,lty=2,color = "grey70",linewidth=0.4)+
  geom_linerange(x=2.5,ymin =lower_s,ymax=upper_s,lty=2,color = "grey70",linewidth=0.4)+
  geom_linerange(x=3.5,ymin =lower_s,ymax=upper_s,lty=2,color = "grey70",linewidth=0.4)+
  geom_linerange(x=4.5,ymin =lower_s,ymax=upper_s,lty=2,color ="grey70",linewidth=0.4)+
  geom_boxplot(aes(fill=Group2,color =Group2),position = position_dodge(0.8),outlier.color =NA,linewidth =1,width =0.5)+
  stat_summary(aes(fill=Group2),fun = mean, geom ="point", size=2,color ="white",show.legend =F,position = position_dodge(0.8))+
  geom_text_aimed(data=df_sig,aes(x=group,y= (upper_s-10), label = p.adj.signif, group = group),angle=90,size=5)+
  theme_void()+coord_polar()+theme(legend.position=c(0.9,0.15))+
  geom_rect(xmin =0.5,xmax=1.5,ymin=(upper_s+1),ymax=(upper_s+10),fill="#3B9AB2")+
  geom_rect(xmin =1.5,xmax=2.5,ymin=(upper_s+1),ymax=(upper_s+10),fill="#7294D4")+
  geom_rect(xmin =2.5,xmax=3.5,ymin=(upper_s+1),ymax=(upper_s+10), fill="#9E8CC6")+
  geom_rect(xmin =3.5,xmax=4.5,ymin=(upper_s+1),ymax=(upper_s+10),fill="#E6A0C4")+
  geom_text(data=df_sig, aes(x=x,y=(upper_s+6),label=group),size=4, angle=c(315,45,-45,45),fontface="bold",inherit.aes =F)+
  scale_y_continuous(limits=c(-20,(upper_s+10))) + labs(title = "Bacteroidetes")+
  theme(plot.title = element_text(hjust = 0.5 ,vjust=-5,face="bold",size = 16),legend.title = element_blank())
ggsave("./fig1_I_2.pdf",p,width = 6,height=6)


##F/B ratio
dat$FB_ratio <- dat$Firmicutes/dat$Bacteroidetes
dat = dat[is.finite(dat$FB_ratio),]
mean_FB = c()
for(i in unique(dat$Group)){
  tps = subset(dat,Group%in%i)
  tempdata = data.frame(Group2 = tps[1,1],Group = tps[1,2],group = tps[1,6],mean_FB = log2(mean(tps$FB_ratio)))
  mean_FB = mean_FB %>% rbind(tempdata)
}

mean_FB%>% mutate(new_x = 1:8) -> mean_FB
number <- nrow(mean_FB)
angle <- 90 - 360*(mean_FB$new_x - 0.5)/number
mean_FB$hjust <- ifelse(angle < -90, 1, 0)
mean_FB$angle <- ifelse(angle < -90, angle + 180, angle)
head(mean_FB)

p = ggplot(data = mean_FB, aes(x = as.factor(new_x), y = mean_FB, fill = Group2)) +
  geom_bar(stat = "identity")+ ylim(-3, 9)+
  coord_polar(start = 0) +scale_fill_manual(values = c('#ee3124','#50bcb6')) +
  theme_minimal() +theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank())+
  geom_text(aes(label = "DEG"),x = -0, y = -750, size = 4, fontface = 'bold', family = "serif")+
  theme(legend.position=c(0.5,0.5),legend.title = element_blank())+
  geom_text(data = mean_FB,aes(x = new_x,y = mean_FB+0.2,
                               label = Group,hjust = hjust),color = 'black', size = 4, fontface = 'bold', 
            angle = mean_FB$angle,inherit.aes = FALSE )+ 
  labs(title = "log[2](F/B ratio)")+theme(plot.title = element_text(hjust = 0.5 ,vjust=-1,face="bold",size = 16))
ggsave("./fig1_J.pdf",p,width = 6,height=6)



### Fig1 K-M

library(tidyverse)
library(glmnet)
library(ggplot2)

# =========================================================
# 0. paths
# =========================================================
otu_file  <- file.path("./otutab.txt")
meta_file <- file.path("./metadata.txt")

# =========================================================
# 1. read data
# =========================================================
otu <- read.table(
  otu_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- read.table(
  meta_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(otu) <- otu[[1]]
otu <- otu[, -1, drop = FALSE]

# =========================================================
# 2. metadata parsing
# =========================================================
meta <- meta %>%
  dplyr::rename(
    sample_id = SampleID,
    group = Group2,
    group_age = Group
  ) %>%
  mutate(
    group = factor(group, levels = c("WT", "RTT")),
    age_month = as.numeric(sub(".*_(\\d+)$", "\\1", group_age)),
    animal_id = sub("^([RC]_\\d+)_.*$", "\\1", sample_id)
  )

common_samples <- intersect(colnames(otu), meta$sample_id)

otu <- otu[, common_samples, drop = FALSE]
meta <- meta %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, colnames(otu)))

otu <- otu[, meta$sample_id, drop = FALSE]
stopifnot(identical(colnames(otu), meta$sample_id))

# =========================================================
# 3. sample x ASV count matrix
# =========================================================
X_counts <- t(as.matrix(otu))
mode(X_counts) <- "numeric"
rownames(X_counts) <- meta$sample_id

# =========================================================
# 4. filter ASVs based on WT prevalence
# =========================================================
wt_idx <- meta$group == "WT"
prev_wt <- colMeans(X_counts[wt_idx, , drop = FALSE] > 0)
keep_asv <- prev_wt >= 0.10
X_counts <- X_counts[, keep_asv, drop = FALSE]

# =========================================================
# 5. relative abundance + CLR
# =========================================================
pseudo <- 0.5
X_rel <- sweep(X_counts + pseudo, 1, rowSums(X_counts + pseudo), FUN = "/")

clr_transform <- function(v) {
  lv <- log(v)
  lv - mean(lv)
}

X_clr <- t(apply(X_rel, 1, clr_transform))
X_clr <- as.matrix(X_clr)

# =========================================================
# 6. WT-only microbiome age model
#    leave-one-animal-out CV
# =========================================================
X_wt <- X_clr[meta$group == "WT", , drop = FALSE]
y_wt <- meta$age_month[meta$group == "WT"]
id_wt <- meta$animal_id[meta$group == "WT"]
foldid_wt <- as.integer(factor(id_wt))

alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)
cv_list <- list()
cv_summary <- data.frame()

for (a in alpha_grid) {
  cvfit <- cv.glmnet(
    x = X_wt,
    y = y_wt,
    family = "gaussian",
    alpha = a,
    foldid = foldid_wt,
    standardize = TRUE,
    type.measure = "mse"
  )
  
  cv_list[[as.character(a)]] <- cvfit
  cv_summary <- rbind(
    cv_summary,
    data.frame(
      alpha = a,
      lambda_min = cvfit$lambda.min,
      lambda_1se = cvfit$lambda.1se,
      cvm_min = min(cvfit$cvm)
    )
  )
}

best_alpha <- cv_summary$alpha[which.min(cv_summary$cvm_min)]
best_cvfit <- cv_list[[as.character(best_alpha)]]

# =========================================================
# 7. WT LOAO prediction
# =========================================================
loao_pred <- rep(NA_real_, length(y_wt))

for (id in unique(id_wt)) {
  tr <- id_wt != id
  te <- id_wt == id
  
  X_tr <- X_wt[tr, , drop = FALSE]
  y_tr <- y_wt[tr]
  id_tr <- id_wt[tr]
  foldid_tr <- as.integer(factor(id_tr))
  
  cvfit_tr <- cv.glmnet(
    x = X_tr,
    y = y_tr,
    family = "gaussian",
    alpha = best_alpha,
    foldid = foldid_tr,
    standardize = TRUE,
    type.measure = "mse"
  )
  
  fit_tr <- glmnet(
    x = X_tr,
    y = y_tr,
    family = "gaussian",
    alpha = best_alpha,
    lambda = cvfit_tr$lambda.1se,
    standardize = TRUE
  )
  
  loao_pred[te] <- as.numeric(predict(fit_tr, newx = X_wt[te, , drop = FALSE]))
}

wt_perf <- data.frame(
  sample_id = meta$sample_id[meta$group == "WT"],
  animal_id = id_wt,
  actual_age = y_wt,
  predicted_age = loao_pred
)

wt_r <- cor(wt_perf$actual_age, wt_perf$predicted_age, method = "pearson")
wt_rmse <- sqrt(mean((wt_perf$predicted_age - wt_perf$actual_age)^2))
wt_mae <- mean(abs(wt_perf$predicted_age - wt_perf$actual_age))

# =========================================================
# 8. fit final model on all WT and predict all samples
# =========================================================
final_fit <- glmnet(
  x = X_wt,
  y = y_wt,
  family = "gaussian",
  alpha = best_alpha,
  lambda = best_cvfit$lambda.1se,
  standardize = TRUE
)

all_pred <- as.numeric(predict(final_fit, newx = X_clr))

res <- meta %>%
  mutate(predicted_microbiome_age = all_pred)

# =========================================================
# 9. age-matched WT centroid distance
# =========================================================
X_clr_df <- as.data.frame(X_clr)
X_clr_df$sample_id <- rownames(X_clr)

feature_cols <- setdiff(colnames(X_clr_df), "sample_id")

wt_centroids <- X_clr_df %>%
  left_join(meta %>% dplyr::select(sample_id, group, age_month), by = "sample_id") %>%
  filter(group == "WT") %>%
  group_by(age_month) %>%
  summarise(across(all_of(feature_cols), mean), .groups = "drop")

dist_to_age_centroid <- sapply(seq_len(nrow(X_clr)), function(i) {
  this_age <- meta$age_month[i]
  centroid <- wt_centroids %>%
    filter(age_month == this_age) %>%
    dplyr::select(all_of(feature_cols)) %>%
    as.numeric()
  
  sqrt(sum((X_clr[i, ] - centroid)^2))
})

res$dist_to_age_matched_WT <- dist_to_age_centroid

# =========================================================
# 10. helper functions
# =========================================================
group_colors <- c(WT = "#2EC7C9", RTT = "#E91E63")

n_label_df <- res %>%
  dplyr::count(group, age_month) %>%
  mutate(label = paste0("n=", n))

wilcox_by_age <- function(data, value_col) {
  ages <- sort(unique(data$age_month))
  out <- lapply(ages, function(a) {
    dat_a <- data %>% filter(age_month == a)
    p <- tryCatch(
      wilcox.test(dat_a[[value_col]] ~ dat_a$group)$p.value,
      error = function(e) NA_real_
    )
    data.frame(age_month = a, p_value = p)
  }) %>% bind_rows()
  
  out %>%
    mutate(
      p_label = case_when(
        is.na(p_value) ~ "ns",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

sig_pred <- wilcox_by_age(res, "predicted_microbiome_age")
sig_dist <- wilcox_by_age(res, "dist_to_age_matched_WT")

# =========================================================
# 11. WT model performance
# =========================================================
anno_text <- paste0(
  "r = ", round(wt_r, 2), "\n",
  "RMSE = ", round(wt_rmse, 2), "\n"
)

p <- ggplot(wt_perf, aes(x = actual_age, y = predicted_age, color = factor(actual_age))) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey45") +
  annotate("text", x = 1.25, y = max(wt_perf$predicted_age), label = anno_text,
           hjust = 0, vjust = 1, size = 4.2) +
  scale_x_continuous(breaks = c(1, 3, 6, 9)) +
  scale_y_continuous(breaks = c(1, 3, 6, 9)) +
  labs(
    title = "WT microbiome age model",
    subtitle = "Leave-one-animal-out validation",
    x = "Chronological age (months)",
    y = "Predicted microbiome age (months)",
    color = "Age"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
ggsave("fig1_K.pdf", p, width = 5, height = 4)

# =========================================================
# 12. predicted microbiome age
# =========================================================
pred_y <- res %>%
  group_by(age_month) %>%
  summarise(y = max(predicted_microbiome_age) + 0.45, .groups = "drop") %>%
  left_join(sig_pred, by = "age_month")

p <- ggplot(res, aes(x = factor(age_month), y = predicted_microbiome_age, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.3,
               position = position_dodge(width = 0.72)) +
  geom_jitter(aes(color = group),
              size = 2.5, alpha = 0.45, stroke = 0.2,
              position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.72)) +
  geom_text(data = n_label_df,
            aes(x = factor(age_month), y = -0.65, label = label, color = group),
            position = position_dodge(width = 0.72),
            size = 3.6, show.legend = FALSE) +
  geom_text(data = pred_y,
            aes(x = factor(age_month), y = y, label = p_label),
            inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = c(-1, max(pred_y$y) + 0.3)) +
  labs(
    title = "Predicted microbiome age",
    x = "Chronological age (months)",
    y = "Predicted microbiome age (months)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
ggsave("fig1_L.pdf", p, width = 5, height = 4)

# =========================================================
# 13. age-matched compositional deviation
# =========================================================
dist_y <- res %>%
  group_by(age_month) %>%
  summarise(y = max(dist_to_age_matched_WT) + 5, .groups = "drop") %>%
  left_join(sig_dist, by = "age_month")

p <- ggplot(res, aes(x = factor(age_month), y = dist_to_age_matched_WT, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.3,
               position = position_dodge(width = 0.72)) +
  geom_jitter(aes(color = group),
              size = 2.5, alpha = 0.45, stroke = 0.2,
              position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.72)) +
  geom_text(data = dist_y,
            aes(x = factor(age_month), y = y, label = p_label),
            inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Age-matched compositional deviation",
    x = "Chronological age (months)",
    y = "Distance to age-matched WT centroid"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )
ggsave("fig1_M.pdf", p, width = 5, height = 4)

