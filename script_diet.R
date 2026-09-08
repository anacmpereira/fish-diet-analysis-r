##Diet Analysis, Comparison Between the Diets of Phalloceros harpagos (native species) and Poecilia reticulata (non-native species)

rm(list=ls()) #Clearing the R memory#

#1: Calculating niche width (Levins' Index) and Occurrence Frequencies#
#1.1: Load data and import libraries#
library(dplyr) #For data manipulation#
library(tidyverse)
dados <- read.csv("diet_database.csv") #Loading the database# 

#1.2: Defining variable types#
str(dados)
dados$sp = factor(dados$sp) #Defining the species variable 'sp' as a factor#
dados$item = factor(dados$item) #Defining the 'item' variable as a factor# 

#1.3: Dividing the matrix into two parts: one for each species
##The first 31 rows are for Poecilia reticulata (pr)
##and rows 32 to 77 are for Phalloceros harpagos (ph)
pr_abundance <- dados[1:31,-1]
ph_abundance <- dados[32:77,-1]

#1.4: Calculating occurrence frequency and niche width for Poecilia reticulata (pr)
##Occurrence frequency of each item in the diet
freq_pr <- colSums(pr_abundance > 0) / nrow(pr_abundance)
print(freq_pr) #Printing occurrence frequency data#
mean(rowSums(pr_abundance > 0))

##Sum of the squared frequencies
soma_pi2_pr <- sum(freq_pr^2)

##Niche width using Levins' index
largura_nicho_pr <- 1 / soma_pi2_pr

##Print results#
cat("Niche width for Poecilia reticulata:", largura_nicho_pr, "\n")

#1.5: Calculating occurrence frequency and niche width for Phalloceros harpagos (ph)
##Occurrence frequency of each item in the diet
freq_ph <- colSums(ph_abundance > 0) / nrow(ph_abundance)
print(freq_ph) #Printing occurrence frequency data#
mean(rowSums(ph_abundance > 0))

##Sum of the squared frequencies
soma_pi2_ph <- sum(freq_ph^2)
##Niche width using Levins' index
largura_nicho_ph <- 1 / soma_pi2_ph

##Printing results
cat("Niche width for Phalloceros harpagos:", largura_nicho_ph, "\n")

###########################################################################
#2: Using NMDS and PERMANOVA technique for graphical visualization of overlap#

#2.1: Loading required libraries
library(vegan)    # For NMDS analysis
library(ggplot2)  # For visualization
library(tidyr)    # For data manipulation
library(dplyr)    # For data manipulation

#2.2: Randomizing the data#
set.seed(123) #To ensure reproducibility#

#2.3: Separating the species column and abundance columns
traco <- dados %>% select(sp)
composi <- dados %>% select(Tricoptera:Vermes)

#2.4: Normalizing the abundance matrix and calculating the dissimilarity matrix
compos.rel <- decostand(composi, method = "hellinger")  # Hellinger transformation
compos.dist <- vegdist(compos.rel, method = "bray")  # Bray-Curtis distance

#2.5: Fitting the NMDS model
nmds_result <- metaMDS(compos.dist, k = 2, maxit = 999, trymax = 1000)
stress_value <- nmds_result$stress
print(stress_value)

#2.6: Verifying NMDS fit
## Calculating goodness-of-fit
gof <- goodness(nmds_result)
stressplot(nmds_result, main = "Shepard Plot")
plot(nmds_result, type = "t", main = "Goodness of Fit")
points(nmds_result, display = "sites", cex = gof * 200)

#2.7: Plotting the NMDS graph#
## Adding species information to NMDS results
data_plot <- as.data.frame(scores(nmds_result, display = "sites"))  # Site scores (individuals)
data_plot$sp <- traco$sp  # Adding species column to results

## Creating convex hulls
group_pr <- data_plot[data_plot$sp == "pr", ][chull(data_plot[data_plot$sp == "pr", c("NMDS1", "NMDS2")]), ]  # Convex hull for pr
group_ph <- data_plot[data_plot$sp == "ph", ][chull(data_plot[data_plot$sp == "ph", c("NMDS1", "NMDS2")]), ]  # Convex hull for ph

## Combining hull data
hull_data <- rbind(group_pr, group_ph)  # Combining both hulls data

## Creating the NMDS plot with ggplot2
ggplot(data_plot, aes(x = NMDS1, y = NMDS2, color = sp, shape = sp)) + 
  geom_point(size = 4, alpha = 0.7) +  # Adding points to the plot
  geom_polygon(data = hull_data, aes(fill = sp, group = sp), alpha = 0.5, color = NA) +  # Adding convex polygons
  scale_color_manual(values = c("ph" = "blue", "pr" = "red"), 
                     labels = c("ph" = expression(italic("Phalloceros harpagos")), 
                                "pr" = expression(italic("Poecilia reticulata")))) +  # Defining labels in italics
  scale_fill_manual(values = c("ph" = "blue", "pr" = "red"), guide = "none") +  # Defining fill colors
  scale_shape_manual(values = c("ph" = 16, "pr" = 17), 
                     labels = c("ph" = expression(italic("Phalloceros harpagos")), 
                                "pr" = expression(italic("Poecilia reticulata")))) +  # Defining shape labels in italics
  labs(x = "NMDS1", y = "NMDS2", 
       color = "Species", shape = "Species",
       title = NULL) +  # Plot title
  theme_minimal() +  # Minimal theme
  theme(plot.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.title = element_text(size = 14))  # Removing grid lines

## Saving the plot in high quality#
ggsave("grafico_nmds2.jpeg", plot = last_plot(), dpi = 300, width = 8, height = 6)

#2.8: Performing PERMANOVA for variation analysis
ru.permanova <- adonis2(composi ~ sp, data = traco, method = "bray")
print(ru.permanova)

###################################################################################
#3: Using GLMM analysis to check the most important food items

#3.1: Loading libraries
library(tidyverse)
library(lme4)
library(ggplot2)
library(performance)
library(DHARMa)
library(car)

#3.2: Transforming the table from wide to long format
dados_long <- dados %>%
  gather(key = "item", value = "abundance", Tricoptera:Vermes) %>% # Transforming columns into rows
  mutate(ind = seq_len(nrow(.))) # Adding a column for individuals

#3.3: Exploratory data analysis
hist(dados_long$abundance) #not normal, many zeros
boxplot(dados_long$abundance ~ dados_long$sp)
boxplot(dados_long$abundance ~ dados_long$item)

#3.4: Fitting GLMM models with Poisson distribution
## Defining reference categories for the analysis to better respond to research objectives
dados_long$item <- relevel(factor(dados_long$item), ref = "Ephemeroptera")
dados_long$sp <- relevel(factor(dados_long$sp), ref = "pr") 

m1 <- glmer(abundance ~ item + sp + (1|ind), data = dados_long, family = "poisson")
summary(m1) #Analyzing m1 coefficients

m2 <- glmer(abundance ~ item * sp + (1|ind), data = dados_long, family = "poisson")
summary(m2) #Analyzing m2 coefficients

## Creating the null model to compare if variation is explained by variables or randomness
m_nulo <- glmer(abundance ~ 1 + (1|ind), data = dados_long, family = "poisson")
summary(m_nulo) #Analyzing m_nulo coefficients

#Note: The Poisson family was chosen because the response variable consists of
#count data, and the data does not follow a normal distribution. The fixed 
#factors are "item" and "sp," which are the variables to be estimated directly.
#The variable 'ind,' representing the individual, was included as a random factor
#to account for uncontrolled variation. This allows the influence of individual 
#variation to be considered without estimating direct coefficients, assuming 
#that each fish may exhibit its own variation in diet.

#Note 2: Two different models were used (in addition to the null model - m_nulo).
#In m1, the model was created by adding the main effects, estimating "item" 
#and "sp" (species) separately on abundance. This means that the average prey
#abundance varies among the different types of items, regardless of fish species,
#and the average prey abundance varies among species without interaction with 
#the type of "item." In m2, the model includes the interaction between "item"
#and "sp." This model allows testing whether certain species consume more of
#specific items than others

#3.5: Checking for overdispersion using the function
check_overdispersion <- function(model) {
  res <- residuals(model, type = "pearson")
  phi <- sum(res^2) / df.residual(model)
  message("Overdispersion factor: ", round(phi, 2))
  if (phi > 1.5) {
    message("Overdispersion detected! Consider a quasi-Poisson or negative binomial distribution.")
  } else {
    message("No significant overdispersion detected.")
  }
}

check_overdispersion(m1) #no overdispersion
check_overdispersion(m2) #no overdispersion
check_overdispersion(m_nulo) #no overdispersion

#3.6: Checking homogeneity
## Important to ensure that the model is fitted correctly

plot(m1)
library(car)
res_m1 <- residuals(m1, type = "pearson")
res_m2 <- residuals(m2, type = "pearson")
res_m_nulo <- residuals(m_nulo, type = "pearson")

fitted_vals_m1 <- fitted(m1)
fitted_vals_m2 <- fitted(m2)
fitted_vals_m_nulo <- fitted(m_nulo)

leveneTest(res_m1 ~ sp, data = dados_long)
leveneTest(res_m2 ~ sp, data = dados_long)
leveneTest(res_m_nulo ~ sp, data = dados_long)

##Note: A high p-value indicates that there are no significant differences in the variances between groups,
#the residual variances are homogeneous for the models.

#3.7: Comparing models using AICc
## Installing and loading the necessary library
install.packages("AICcmodavg") 
library(AICcmodavg)

models <- list(m1, m2, m_nulo)
resultado_aicc <- aictab(models)
print(resultado_aicc)

#Conclusion: m1 is sufficient to explain the variability in abundance without the need to include interactions between the two variables

##3.8: Creating a plot with GLMM coefficients
dados_graf <- data.frame(
  Variable = c("(Intercept)", "Diptera", "Hymenoptera", "Megaloptera", "Odonata", "Tricoptera", "Vermes", "sp: Phalloceros harpagos"),
  Estimate = c(-4.1240, 4.2185, -0.7252, -0.7287, 1.0544, 0.6670, -0.7205, -0.4211),
  Std_Err = c(0.7244, 0.7152, 1.2296, 1.2311, 0.8197, 0.8688, 1.2288, 0.2350),
  p_value = c(1.25e-08, 3.68e-09, 0.5554, 0.5539, 0.1983, 0.4427, 0.5577, 0.0731)
)

##Create 'Color' column to color the points based on p_value values
dados_graf$Color <- ifelse(dados_graf$p_value < 0.01, "red", "black")

##Set the desired order for the 'Variable' factor from top to bottom, reversing the levels
dados_graf$Variable <- factor(dados_graf$Variable, 
                              levels = rev(c("(Intercept)", "Diptera", "Hymenoptera", "Megaloptera", 
                                             "Odonata", "Tricoptera", "Vermes", "sp: Phalloceros harpagos")))

## Create plot with ggplot2
graf_high_quality <- ggplot(dados_graf, aes(x = Estimate, y = Variable)) +
  geom_point(aes(color = Color), size = 8) +  # Points colored according to significance
  geom_errorbarh(aes(xmin = Estimate - Std_Err, xmax = Estimate + Std_Err), height = 0.2) +  # Error bars
  scale_color_manual(values = c("black", "red"), 
                     labels = c("p < 0.01", "p > 0.01")) +  # Customize legend labels
  theme_minimal() +
  labs(
    title = NULL,
    x = "Regression coefficients estimates",
    y = "Groups",
    color = "Significance"
  ) +
  theme(
    axis.text.x = element_text(size = 16),  # Increase the font size of the X-axis labels
    axis.text.y = element_text(size = 16),  # Increase the font size of the Y-axis labels
    axis.title.x = element_text(size = 16, face = "bold"),  # Bold and larger font for X-axis title
    axis.title.y = element_text(size = 16, face = "bold"),  # Bold and larger font for Y-axis title
    plot.title = element_text(hjust = 0.5),  # Center the title (it doesn't make a difference without a title)
    legend.position = "none"  # Completely remove the legend
  )

print(graf_high_quality)

## Saving the plot
ggsave("graf_high_quality5.png", plot = graf_high_quality, width = 12, height = 8, dpi = 300)

##3.9: Creating a plot with log-transformed mean abundance values

summary_data <- dados_long %>%
  group_by(sp, item) %>%
  summarise(mean_abundance = mean(abundance), 
            se_abundance = sd(abundance) / sqrt(n())) %>%
  ungroup()

graf_abundance <- ggplot(summary_data, aes(x = item, y = mean_abundance, fill = sp)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
                position = position_dodge(0.9), width = 0.2) +
  scale_y_log10() +
  labs(title = "Mean abundance of food items per species",
       x = "Food item", 
       y = "Mean abundance (log)",
       fill = "Species") +  
  scale_fill_manual(values = c("ph" = "skyblue", "pr" = "lightcoral"), 
                    labels = c("ph" = expression(italic("Phalloceros harpagos")),
                               "pr" = expression(italic("Poecilia reticulata")))) + 
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Arial", size = 12, face = "bold"),
    axis.title.y = element_text(family = "Arial", size = 11, face = "italic"),
    axis.title.x = element_text(family = "Arial", size = 11, face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid = element_blank())

print(graf_abundance) #Print the plot

##Saving the high-quality plot
ggsave("graf_abundance.jpeg", plot = last_plot(), dpi = 300, width = 8, height = 6)

############################################################

# Calculating Pianka Index and Schoener Index

# Pianka Index using relative abundance #

rm(list = ls())  # clearing workspace memory
# Loading the file: Import > Import text (readr) > diet_database

dados <- diet_database # Renaming the imported dataset

# 1: Loading required packages
library(dplyr)
library(tibble)

# 2: Aggregating abundance by species
diet_species <- dados %>% 
  group_by(sp) %>% 
  summarise(across(where(is.numeric), sum))  
# Summing all food items for all individuals of each species

# 3: Converting abundance values into proportions
diet_prop <- diet_species %>% 
  column_to_rownames("sp")  # Setting species column as row names

diet_prop <- sweep(diet_prop, 1, rowSums(diet_prop), "/")
# Dividing each row by its sum to obtain relative abundance

# 4: Creating Pianka Index function
pianka_index <- function(p1, p2) {
  sum(p1 * p2) / sqrt(sum(p1^2) * sum(p2^2))
  # Formula: Ojk = Σ(pij * pik) / sqrt(Σ(pij²) * Σ(pik²))
}

# 5: Calculating Pianka Index between the two species
pianka <- pianka_index(
  diet_prop["ph", ],  # species 1
  diet_prop["pr", ]   # species 2
)

pianka  # Pianka index result

# 6: Creating Schoener Index function
schoener_index <- function(p1, p2) { 
  1 - 0.5 * sum(abs(p1 - p2), na.rm = TRUE)
  # Formula: D = 1 - 0.5 * Σ|pij - pik|
}

# Calculating Schoener Index between the same two species
schoener <- schoener_index(
  diet_prop["ph", ],  # species 1
  diet_prop["pr", ]   # species 2
)

schoener  # Schoener index result