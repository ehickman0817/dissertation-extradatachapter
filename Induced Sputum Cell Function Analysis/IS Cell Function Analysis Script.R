rm(list = ls(all.names = TRUE))

# If needed, install and load packages
library(readxl) # for data import 
library(table1) # to make tables
library(tidyverse) # for data cleaning and organization
library(ggplot2) # for bar plot
library(ggsignif) # for bar plot
library(ggpubr) # for bar plot

##############################################################
####### 1. Import data
##############################################################

all_data <- data.frame(read_excel("Induced Sputum Cell Function Analysis/Input Data/2021_12_28 Sputum Cell Functional Data.xlsx"))

##############################################################
####### 2. Demographics table
##############################################################

# Change labels
all_data$Group <- factor(all_data$Group, 
                         levels = c("Nonvaper", "4th Gen"),
                         labels = c("Non-vapers", "4th Gen E-Cig Users"))

all_data$Sex <- factor(all_data$Sex, 
                       levels=c("F", "M"),
                       labels=c("Female", "Male"))

all_data$Race <- factor(all_data$Race, 
                        levels=c("W", "B", "AZN"),
                        labels=c("White", "Black", "Asian"))

all_data$Hispanic <- factor(all_data$Hispanic, 
                            levels=c("NO", "YES"),
                            labels=c("No", "Yes"))

all_data$Device <- factor(all_data$Device, 
                          levels=c("JUUL", "JUUL/Disposable", "Disposable", "None"))

# Determine whether continuous variables are normally distributed
shapiro.test(all_data$Age) 
shapiro.test(all_data$BMI)

# Function for adding p-values to table
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a non-parametric t-test (Wilcox Test)
    p <- wilcox.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a Fisher's exact test (or Chi squared if there are expected frequencies >5)
    # Simulate p value = TRUE added because groups are too small
    p <- fisher.test(table(y, g), simulate.p.value = TRUE)$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

# Make demographics table
table1(~ Sex + Race + Hispanic + Age + BMI + Device | Group,
       data = all_data, 
       render.missing = NULL,
       overall = NULL,
       extra.col = list(`P-value` = pvalue))

##############################################################
####### 3. Sputum cell characteristics table
##############################################################

# Statistical test of normality
diffs_normality <- apply(all_data[11:16], 2, shapiro.test)

# Create data frame to summarize results
diffs_normality <- do.call(rbind.data.frame, diffs_normality )
diffs_normality <- format(diffs_normality, scientific = FALSE)

# Add column for normality conclusion
diffs_normality <- diffs_normality %>% mutate(normal = ifelse(p.value < 0.05, F, T)) # majority are non-normal

# Function for custom table so that Mean (SEM) is shown for continuous variables
my.render.cont.mean.sem <- function(x) {
  s <- stats.default(x)
  s$SEM <- with(s, SD/sqrt(N))
  with(stats.apply.rounding(s), c("",
                                  "Mean (SEM)"=sprintf("%s (%s)", MEAN, SEM)))
}

# Relabel rows
label(all_data$SelectSampleWeight) <- "Select Sample Weight"
label(all_data$TotalCellCount) <- "Total Cell Count"
label(all_data$CellsPerMG) <- "Cells/mg"
label(all_data$PercPMN) <- "% PMN"
label(all_data$PercMac) <- "% Macrophage"
label(all_data$PercViable) <- "% Viable"

# Make table
table1(~ SelectSampleWeight + TotalCellCount + CellsPerMG + PercPMN + PercMac + PercViable | Group,
       data = all_data, 
       render.missing = NULL,
       render.continuous = my.render.cont.mean.sem,
       overall = NULL,
       extra.col = list(`P-value` = pvalue))

##############################################################
####### 4. Overall demographics assessment 
##############################################################

## JC-1
shapiro.test(all_data$JC1)

JC1_sex_wilcox <- wilcox.test(JC1 ~ Sex, all_data)
JC1_sex_wilcox$p.value

JC1_sex_ttest <- t.test(JC1 ~ Sex, all_data)
JC1_sex_ttest$p.value

JC1_spearman_age <- cor.test(all_data$Age, all_data$JC1, method = "spearman")
JC1_spearman_age

JC1_spearman_BMI <- cor.test(all_data$BMI, all_data$JC1, method = "spearman")
JC1_spearman_BMI

## Phagocytosis
shapiro.test(all_data$Phagocytosis)

phagocytosis_sex <- wilcox.test(Phagocytosis ~ Sex, all_data)
phagocytosis_sex$p.value

phagocytosis_spearman_age <- cor.test(all_data$Age, all_data$Phagocytosis, method = "spearman")
phagocytosis_spearman_age

phagocytosis_spearman_BMI <- cor.test(all_data$BMI, all_data$Phagocytosis, method = "spearman")
phagocytosis_spearman_BMI

##############################################################
####### 5. Group Comparisons
##############################################################

## JC-1

# Statistical tests
JC1_group_ttest <- t.test(JC1 ~ Group, all_data)
JC1_group_ttest$p.value

JC1_group_wilcox <- wilcox.test(JC1 ~ Group, all_data)
JC1_group_wilcox$p.value

# Set theme
theme_set(theme_bw())

# Open graphical device
pdf("Induced Sputum Cell Function Analysis/Output Figures/SputumCellJC1Groups.pdf",
    colormodel = "cmyk",
    width = 7,
    height = 5.5)

set.seed(0817)

# Make graph
ggplot(all_data, aes(x = Group, y = JC1, fill = Group)) + 
  scale_fill_manual(values = c("gray74", "gray34")) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(0.15),
              shape = 19,
              size = 3) +
  ggsignif::geom_signif(comparisons = list(c("Non-vapers", "4th Gen E-Cig Users")),
                        map_signif_level=TRUE,
                        test = "t.test",
                        size = 1,
                        textsize = 16,
                        margin_top = 0.2) +
  ylim(0, 5) +
  ylab("JC-1 MFI (Red/Green)\n") +
  theme(legend.position = "None",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, face = "bold", color = "black"),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 0.5),
        axis.text.y = element_text(size = 18, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# Close graphical device
dev.off()

## PHAGOCYTOSIS
phagocytosis_group_wilcox <- wilcox.test(Phagocytosis ~ Group, all_data)
phagocytosis_group_wilcox$p.value

# Open graphical device
pdf("Induced Sputum Cell Function Analysis/Output Figures/SputumCellPhagocytosisGroups.pdf",
    colormodel = "cmyk",
    width = 7,
    height = 5.5)

set.seed(0817)

ggplot(all_data, aes(x = Group, y = Phagocytosis, fill = Group)) + 
  scale_fill_manual(values = c("gray74", "gray34")) +
  geom_boxplot(color = "black",
               outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.15),
              shape = 19,
              size = 3) +
  ggsignif::geom_signif(comparisons = list(c("Non-vapers", "4th Gen E-Cig Users")),
                        map_signif_level=TRUE,
                        test = "t.test",
                        size = 1,
                        textsize = 16,
                        margin_top = 0.2) +
  ylim(0, 60000) +
  ylab("Phagocytosis (MFI)\n") +
  theme(legend.position = "None",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, face = "bold", color = "black"),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 0.5),
        axis.text.y = element_text(size = 18, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# Close graphical device
dev.off()

##############################################################
####### 6. JC-1 with CCCP Control Analysis
##############################################################

# Test normality
shapiro.test(all_data$JC1withCCCP)

# Compare groups with paired t test
wilcox.test(all_data$JC1, all_data$JC1withCCCP, paired = TRUE, alternative = "two.sided")

# Filter data to include only what we want to graph
all_data_forgraph <- all_data[c("SubjectID", "JC1", "JC1withCCCP")]
all_data_forgraph <- na.omit(all_data_forgraph)
all_data_forgraph <- all_data_forgraph %>% pivot_longer(!SubjectID, names_to = "Group", values_to = "MFI")

# Make labels prettier
all_data_forgraph$Group <- factor(all_data_forgraph$Group,
                                  levels = c("JC1", "JC1withCCCP"),
                                  labels = c("JC-1", "JC-1 + CCCP"))

# Open graphical device
pdf("Induced Sputum Cell Function Analysis/Output Figures/SputumCellJC1vCCCP.pdf",
    colormodel = "cmyk",
    width = 5,
    height = 5.5)

# Graph
ggplot(all_data_forgraph, aes(x = Group, y = MFI)) +
  scale_fill_manual(values = c("cadetblue", "darkseagreen")) +
  geom_boxplot(aes(fill = Group), color = "black") +
  stat_compare_means(comparisons = list(c("JC-1", "JC-1 + CCCP")),
                     map_signif_level=TRUE,
                     label = "p.signif",
                     paired = TRUE,
                     size = 16,
                     bracket.size = 1,
                     label.y = 3.9) +
  ylim(0, 5) +
  geom_line(aes(group = SubjectID)) +
  geom_point(size = 2) +
  ylab("MFI (Red/Green)\n") +
  theme(legend.position = "None",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, face = "bold", color = "black"),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 0.5),
        axis.text.y = element_text(size = 18, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# Close graphical device
dev.off()
