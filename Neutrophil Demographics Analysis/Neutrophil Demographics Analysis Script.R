rm(list = ls(all.names = TRUE))

# If needed, install and load packages
library(readxl) # for data import 
library(table1) # to make tables
library(tidyverse) # for data cleaning and organization
library(ggplot2) # for bar plot
library(ggsignif) # for bar plot
library(corrplot) # for correlation plot
library (RColorBrewer) # for correlation plot

##############################################################
####### 1. Import data
##############################################################

## SEAHORSE DATA

# Import data
seahorse_data <- data.frame(read_excel("Neutrophil Demographics Analysis/Input Data/2021_12_14 Seahorse Demographics Analysis.xlsx"))

# Make sure variables are numeric that should be numeric. 
# They sometimes read in as characters because of how the study coordinators extract them from the database.
seahorse_data$Age <- as.numeric(seahorse_data$Age)
seahorse_data$BMI <- as.numeric(seahorse_data$BMI)

## PHAGOCYTOSIS DATA

# Import data
phagocytosis_data <- data.frame(read_excel("Neutrophil Demographics Analysis/Input Data/2021_12_14 Phagocytosis Demographics Analysis.xlsx"))

# Make sure variables are numeric that should be numeric. 
# They sometimes read in as characters because of how the study coordinators extract them from the database.
phagocytosis_data$Age <- as.numeric(phagocytosis_data$Age)
phagocytosis_data$BMI <- as.numeric(phagocytosis_data$BMI)


##############################################################
####### 2. Make demographics table
##############################################################

## SEAHORSE DATA

# Change labels
seahorse_data$Race <- factor(seahorse_data$Race, 
                             levels=c("W", "B", "A", "MO"), 
                             labels=c("White", "Black", "Asian", "Mixed/Other"))

seahorse_data$Sex <- factor(seahorse_data$Sex, 
                            levels=c("F", "M"),
                            labels=c("Female", "Male"))

# Make table.
table1(~ Sex + Race + Hispanic + Age + BMI,
       data = seahorse_data, 
       render.missing = NULL)

## PHAGOCYTOSIS DATA

# Change labels
phagocytosis_data$Race <- factor(phagocytosis_data$Race, 
                                 levels=c("W", "B", "A", "MO"), 
                                 labels=c("White", "Black", "Asian", "Mixed/Other"))

phagocytosis_data$Sex <- factor(phagocytosis_data$Sex, 
                                levels=c("F", "M"),
                                labels=c("Female", "Male"))

# Make table
table1(~ Sex + Race + Hispanic + Age + BMI,
       data = phagocytosis_data, 
       render.missing = NULL)


##############################################################
####### 3. Assess normality of phagocytosis and seahorse data
##############################################################

## SEAHORSE DATA

# Statistical test of normality
seahorse_normality <- apply(seahorse_data[9:16], 2, shapiro.test)

# Create data frame to summarize results
seahorse_normality <- do.call(rbind.data.frame, seahorse_normality )
seahorse_normality <- format(seahorse_normality, scientific = FALSE)

# Add column to adjust for multiple hypothesis testing
seahorse_normality$p.value.adj <- p.adjust(seahorse_normality$p.value, "BH")

# Add column for normality conclusion
seahorse_normality <- seahorse_normality %>% mutate(normal = ifelse(p.value.adj < 0.05, F, T)) 

# Make vector of names for Seahorse metrics that are normally distributed.
normal_seahorse_metrics <- filter(seahorse_normality, normal == 'TRUE')
normal_seahorse_metrics = rownames(normal_seahorse_metrics)

## PHAGOCYTOSIS DATA
phagocytosis_normality <- shapiro.test(phagocytosis_data$Phagocytosis)
phagocytosis_normality$p.value # p = 0.1772016 means normally distributed

##############################################################
####### 4. Statistical Analysis
##############################################################

## SEAHORSE 

# T-test for normally distributed data and Wilcox test for non-normally distributed data

# Create results data frame
seahorse_sex = data.frame()

# For loop that runs a t-test if the data are normally distributed and a wilcox test if the data are non-normally distributed.
for (i in 9:ncol(seahorse_data)) {
  if (names(seahorse_data[i]) %in% normal_seahorse_metrics == TRUE) {
    test <- t.test(as.formula(paste0(names(seahorse_data)[i], "~", "Sex", sep = "")), seahorse_data)
  } else {test <- wilcox.test(as.formula(paste0(names(seahorse_data)[i], "~", "Sex", sep = "")), seahorse_data) 
  }
  
  pval <- data.frame(test$p.value, row.names = noquote(paste0(names(seahorse_data[i]))))
  pval$method <- test$method
  
  seahorse_sex <- bind_rows(seahorse_sex, pval)
  
}

# View results
seahorse_sex

# Save results
seahorse_sex$metric <- rownames(seahorse_sex)
write.table(seahorse_sex, file = "Neutrophil Demographics Analysis/Output Tables/SeahorseSexTtest.csv", sep = "\t", row.names = FALSE)

# Correlation testing for age and BMI. Spearman's correlation chosen because it can detect both linear and monotonic correlations, 
# whereas Pearson's correlation only detects linear correlations. 

seahorse_spearman_age <- apply(seahorse_data[, 9:16], 2, cor.test, seahorse_data$Age, method = "spearman")
seahorse_spearman_age <- as.data.frame(do.call(rbind, seahorse_spearman_age))
seahorse_spearman_age

seahorse_spearman_age$metric <- rownames(seahorse_spearman_age)
seahorse_spearman_age <- apply(seahorse_spearman_age,2,as.character)
write.table(seahorse_spearman_age, file = "Neutrophil Demographics Analysis/Output Tables/SeahorseAgeSpearman.csv", sep = "\t", row.names = FALSE)

seahorse_spearman_bmi <- apply(seahorse_data[, 9:16], 2, cor.test, seahorse_data$BMI, method = "spearman")
seahorse_spearman_bmi <- as.data.frame(do.call(rbind, seahorse_spearman_bmi))
seahorse_spearman_bmi

seahorse_spearman_bmi$metric <- rownames(seahorse_spearman_bmi)
seahorse_spearman_bmi <- apply(seahorse_spearman_bmi,2,as.character)
write.table(seahorse_spearman_bmi, file = "Neutrophil Demographics Analysis/Output Tables/SeahorseBMISpearman.csv", sep = "\t", row.names = FALSE)

## PHAGOCYTOSIS 

phagocytosis_sex <- t.test(Phagocytosis ~ Sex, phagocytosis_data)
phagocytosis_sex$p.value

phagocytosis_spearman_age <- cor.test(phagocytosis_data$Age, phagocytosis_data$Phagocytosis, method = "spearman")
phagocytosis_spearman_age

phagocytosis_spearman_BMI <- cor.test(phagocytosis_data$BMI, phagocytosis_data$Phagocytosis, method = "spearman")
phagocytosis_spearman_BMI

## PLOT SEX DIFFERENCES IN PHAGOCYTOSIS

# Set theme
theme_set(theme_bw())

# Open graphical device
pdf("Neutrophil Demographics Analysis/Output Figures/PhagocytosisSexGraph.pdf",
    colormodel = "cmyk",
    width = 5,
    height = 5.5)

# Set seed so that jitters/points look the same each time
set.seed(8016)

# Make graph
ggplot(phagocytosis_data, aes(x = Sex, y = Phagocytosis, fill = Sex)) + 
  scale_fill_manual(values = c("tomato", "steelblue1")) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_jitter(0.15),
              shape = 19,
              size = 3) +
  geom_signif(comparisons = list(c("Female", "Male")),
              map_signif_level=TRUE,
              size = 1,
              textsize = 16,
              margin_top = 0.2) +
  ylim(0, 80000) +
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
####### 5. Correlation between phagocytosis and seahorse metrics - all subjects together
##############################################################

## PREPARE DATA

# Make phagocytosis data frame for merging
phagocytosis_data_trimmed <- subset(phagocytosis_data, TRUE, c("SubjectNumber", "Phagocytosis"))

# Combined Data
combined_data <- transform(merge(phagocytosis_data_trimmed, seahorse_data, by = "SubjectNumber"))

# Make data frame for just correlations
combined_data_trimmed <- combined_data[-c(3:9)]
combined_data_trimmed <- data.frame(combined_data_trimmed[, -1], row.names = combined_data_trimmed$SubjectNumber)

## MAKE DEMOGRAPHICS TABLE

table1(~ Sex + Race + Hispanic + Age + BMI,
       data = combined_data, 
       render.missing = NULL)

## MAKE CORRELATION PLOT

# Convert to data matrix
combined_data_trimmed_dm <- as.matrix(combined_data_trimmed)

# Create correlation data
correlations <- cor(combined_data_trimmed_dm, method = "spearman")

# To add significance to correlation plot, run this function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat_correlations <- cor.mtest(combined_data_trimmed_dm)

# Open graphical device
pdf("Neutrophil Demographics Analysis/Output Figures/NeutrophilFunctionCorrelationsAllSubjects.pdf",
    colormodel = "cmyk")

# Make plot
corrplot(correlations, 
         method = "circle", 
         type = "upper",
         order = "original",
         tl.col = "black",
         p.mat = p.mat_correlations,
         sig.level = 0.05,
         tl.srt = 45,
         insig = "label_sig",
         pch.col = 'white',
         diag = FALSE,
         col=brewer.pal(n=10, name="RdYlBu"))

# Close Graphical Device
dev.off()

##############################################################
####### 6. Correlation between phagocytosis and seahorse metrics - sexes separated
##############################################################

### MALES

# Filter to only males
combined_data_male <- filter(combined_data, Sex == "Male")

# Make data frame with only variables of interest for correlating
combined_data_male_trimmed <- combined_data_male[-c(3:9)]
combined_data_male_trimmed <- data.frame(combined_data_male_trimmed[, -1], row.names = combined_data_male_trimmed$SubjectNumber)

# Convert to data matrix
combined_data_male_trimmed_dm <- as.matrix(combined_data_male_trimmed)

# Create correlation data
correlations_male <- cor(combined_data_male_trimmed_dm, method = "spearman")

# Significance of correlations
p.mat_correlations_male <- cor.mtest(combined_data_male_trimmed_dm)

# PLOT CORRELATION MATRIX

# Open graphical device
pdf("Neutrophil Demographics Analysis/Output Figures/NeutrophilFunctionCorrelations_Male.pdf",
    colormodel = "cmyk")

corrplot(correlations_male, 
         method = "circle", 
         type = "upper",
         order = "original",
         tl.col = "black",
         p.mat = p.mat_correlations_male,
         sig.level = 0.05,
         tl.srt = 45,
         insig = "label_sig",
         pch.col = 'white',
         diag = FALSE,
         col=brewer.pal(n=10, name="RdYlBu"))

# Close Graphical Device
dev.off()

### FEMALES

# Repeat same code but for the female data

combined_data_female <- filter(combined_data, Sex == "Female")

# Make data frame with only variables of interest for correlating
combined_data_female_trimmed <- combined_data_female[-c(3:9)]
combined_data_female_trimmed <- data.frame(combined_data_female_trimmed[, -1], row.names = combined_data_female_trimmed$SubjectNumber)

# Convert to data matrix
combined_data_female_trimmed_dm <- as.matrix(combined_data_female_trimmed)

# Create correlation data
correlations_female <- cor(combined_data_female_trimmed_dm, method = "spearman")

# Significance of correlations
p.mat_correlations_female <- cor.mtest(combined_data_female_trimmed_dm)

# PLOT CORRELATION MATRIX

# Open graphical device
pdf("Neutrophil Demographics Analysis/Output Figures/NeutrophilFunctionCorrelations_Female.pdf",
    colormodel = "cmyk")

corrplot(correlations_female, 
         method = "circle", 
         type = "upper",
         order = "original",
         tl.col = "black",
         p.mat = p.mat_correlations_female,
         sig.level = 0.1,
         tl.srt = 45,
         insig = "label_sig",
         pch.col = 'white',
         diag = FALSE,
         col=brewer.pal(n=10, name="RdYlBu"))

# Close Graphical Device
dev.off()

## PLOT SPECIFIC CORRELATIONS

# Open graphical device
pdf("Neutrophil Demographics Analysis/Output Figures/Phagocytosis versus OCR Time to Max.pdf",
    colormodel = "cmyk")

# Scatterplot with regression line
theme_set(theme_bw())

ggplot(combined_data, aes(x= Phagocytosis, y=OCRTimeToMax, color=Sex, shape=Sex)) +
  geom_point(size = 5) + 
  scale_color_manual(values = c("tomato", "steelblue1")) +
  geom_smooth(method=lm, se=FALSE, linetype = "dashed") +
  ylab("OCR Time to Max (minutes)\n") +
  xlab("\nPhagocytosis (MFI)") +
  theme( axis.title.x = element_text(size = 20, face = "bold", color = "black"),
         axis.title.y = element_text(size = 20, face = "bold", vjust = 0.5),
         axis.text.y = element_text(size = 18, color = "black"),
         axis.text.x = element_text(size = 18, color = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(colour = "black", fill=NA, size=1),
         legend.text = element_text(size = 16),
         legend.position = c(0.25, 0.85),
         legend.title = element_blank(),
         legend.background = element_rect(fill = "white", 
                                          size = 0.3, 
                                          linetype = "solid",
                                          color = "black"))

# Close Graphical Device
dev.off()

# Open graphical device
pdf("Neutrophil Demographics Analysis/Output Figures/Phagocytosis versus ECAR Time to Max.pdf",
    colormodel = "cmyk")

# Scatterplot with regression line
theme_set(theme_bw())

ggplot(combined_data, aes(x= Phagocytosis, y=ECARTimeToMax, color=Sex, shape=Sex)) +
  geom_point(size = 5) + 
  scale_color_manual(values = c("tomato", "steelblue1")) +
  geom_smooth(method=lm, se=FALSE, linetype = "dashed") +
  ylab("ECAR Time to Max (minutes)\n") +
  xlab("\nPhagocytosis (MFI)") +
  theme( axis.title.x = element_text(size = 20, face = "bold", color = "black"),
         axis.title.y = element_text(size = 20, face = "bold", vjust = 0.5),
         axis.text.y = element_text(size = 18, color = "black"),
         axis.text.x = element_text(size = 18, color = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(colour = "black", fill=NA, size=1),
         legend.text = element_text(size = 16),
         legend.position = c(0.25, 0.85),
         legend.title = element_blank(),
         legend.background = element_rect(fill = "white", 
                                          size = 0.3, 
                                          linetype = "solid",
                                          color = "black"))

# Close Graphical Device
dev.off()


