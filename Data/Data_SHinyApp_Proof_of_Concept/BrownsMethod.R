#Brown's Method

#p.values are really small using fisher 
#because the indicators are derived from the same data and correlated
#this is only useful as a relative ranking
#use Brown’s method (an extension of Fisher’s that adjusts the χ² degrees of freedom for correlation)

#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("EmpiricalBrownsMethod", type = "source")
library(EmpiricalBrownsMethod)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)



#populate matrix using Brown's method
merged_coefficients_df<-read.csv(here::here("Data/Data_SHinyApp_Proof_of_Concept/Scaled_Slopes_all_indicators.csv"))  

#make df wide

#separate by region, and make df wide 
USA.p<-subset(merged_coefficients_df, merged_coefficients_df$Region =="USA")
Canada.p<-subset(merged_coefficients_df, merged_coefficients_df$Region =="Canada")

USA.p_wide <- USA.p %>%
  select(Region, Indicator, p.value) %>%     
  pivot_wider(
    names_from = Indicator,                   # new columns = each indicator
    values_from = p.value                    # fill with corresponding estimates
  )
Canada.p_wide <- Canada.p %>%
  select(Region, Indicator, p.value) %>%    
  pivot_wider(
    names_from = Indicator,                   
    values_from = p.value                    
  )
#range edge values have different names because candaa is the leading and usa is the trailing
names(Canada.p_wide)
names(Canada.p_wide)[c(4,5)] <- c("Range_Edge_East", "Range_Edge_North")
names(USA.p_wide)
names(USA.p_wide)[c(4,5)] <- c("Range_Edge_East", "Range_Edge_North")

#Step 1: Brown’s Method function
#Standard version: small combined p-values mean a stronger joint relationship between indicator pairings.
combine_Browns <- function(p1, p2, indicator_data) {
  combined_p <- EmpiricalBrownsMethod::empiricalBrownsMethod(
    data_matrix = indicator_data,
    p_values = c(p1, p2)
  )
  return(combined_p)
}

#combine_fisher <- function(p1, p2) {
#  stat <- -2 * (log(p1) + log(p2))
#  df <- 4  # 2 * number of p-values combined
#  combined_p <- 1 - pchisq(stat, df)
#  return(combined_p)
#}

#Step 2: Function to combine indicators within a single-region dataframe
calculate_Browns_pairs <- function(df) {
  # Combine east/north indicators into single mean-based metrics
  df <- df %>%
    mutate(
      COG = rowMeans(select(., COG_East, COG_North), na.rm = TRUE),
      Range_Edge = rowMeans(select(., Range_Edge_East, Range_Edge_North), na.rm = TRUE)
    ) %>%
    select(-COG_East, -COG_North, -Range_Edge_East, -Range_Edge_North)
  # Create all unique indicator pairs
  pairs <- combn(names(df)[!names(df) %in% "Region"], 2, simplify = FALSE)
  # Compute combined p-values for each pair
  Browns_results <- map_dfr(pairs, function(x) {
    p1 <- df[[x[1]]]
    p2 <- df[[x[2]]]
    #combined <- combine_Browns(p1, p2)
    combined <- combine_Browns(p1, p2, df[, c(x[1], x[2])])#estiamate covariance using the original data
    data.frame(Indicator1 = x[1],
               Indicator2 = x[2],
               Combined_P = combined)
  }) %>%
    arrange(Combined_P)
  
  return(Browns_results)
}

canada_Browns <- calculate_Browns_pairs(Canada.p_wide)
usa_Browns <- calculate_Browns_pairs(USA.p_wide)
head(canada_Browns)

#3. Add some new fields: significance and rankings
#classify the combined p.values
usa_Browns <- usa_Browns %>%
  mutate(
    p.significant = case_when(
      Combined_P < 0.001 ~ "Very strong",
      Combined_P < 0.01  ~ "Strong",
      Combined_P < 0.05  ~ "Moderate",
      TRUE ~ "Not significant"
    )
  )
canada_Browns <- canada_Browns %>%
  mutate(
    p.significant = case_when(
      Combined_P < 0.001 ~ "Very strong",
      Combined_P < 0.01  ~ "Strong",
      Combined_P < 0.05  ~ "Moderate",
      TRUE ~ "Not significant"
    )
  )
# make sure p.significant is an ordered factor
usa_Browns$p.significant <- factor(
  usa_Browns$p.significant,
  levels = c("Not significant", "Marginal", "Strong", "Very strong")
)
canada_Browns$p.significant <- factor(
  canada_Browns$p.significant,
  levels = c("Not significant", "Marginal", "Strong", "Very strong")
)
#rank the combined p.values from smalles to largest
usa_Browns <- usa_Browns %>%
  mutate(Rank = rank(Combined_P, ties.method = "first"))
canada_Browns <- canada_Browns %>%
  mutate(Rank = rank(Combined_P, ties.method = "first"))
str(usa_Browns)

usa_Browns$Combined_P<-format(usa_Browns$Combined_P, scientific = FALSE)
canada_Browns$Combined_P<-format(canada_Browns$Combined_P, scientific = FALSE)

#Step 4: Plotting/ Interpretation
library(ggplot2)

# rename the indicators so it plots easier
unique(canada_Browns$Indicator1)
unique(canada_Browns$Indicator2)

canada_Browns <- canada_Browns %>%
  mutate(
    Indicator1 = recode(Indicator1,
                        "COG" = "COG",                     
                        "Distance to Border"  = "DtoB",    
                        "Area Occupied" = "AO" ,          
                        "Range_Edge"  =  "RE",          
                        "Abundance Weighted Depth" = "AWD",
                        "Abundance" = "Abun"),
    Indicator2 = recode(Indicator2,
                        "COG" = "COG",                     
                        "Distance to Border"  = "DtoB",    
                        "Area Occupied" = "AO" ,          
                        "Range_Edge"  =  "RE",          
                        "Abundance Weighted Depth" = "AWD",
                        "Abundance" = "Abun"))
usa_Browns <- usa_Browns %>%
  mutate(
    Indicator1 = recode(Indicator1,
                        "COG" = "COG",                     
                        "Distance to Border"  = "DtoB",    
                        "Area Occupied" = "AO" ,          
                        "Range_Edge"  =  "RE",          
                        "Abundance Weighted Depth" = "AWD",
                        "Abundance" = "Abun"),
    Indicator2 = recode(Indicator2,
                        "COG" = "COG",                     
                        "Distance to Border"  = "DtoB",    
                        "Area Occupied" = "AO" ,          
                        "Range_Edge"  =  "RE",          
                        "Abundance Weighted Depth" = "AWD",
                        "Abundance" = "Abun"))

#Plot Canada Matrix
#factor levels for even plotting
all_levels <- union(
  unique(na.omit(canada_Browns$Indicator1)),
  unique(na.omit(canada_Browns$Indicator2))
)
canada_Browns$Indicator1 <- factor(canada_Browns$Indicator1, levels = all_levels)
canada_Browns$Indicator2 <- factor(canada_Browns$Indicator2, levels = all_levels)
#some of the relations need to be swapped for plotting
rows_to_swap <- c(6, 8)  # for example

for (i in rows_to_swap) {
  tmp <- canada_Browns$Indicator1[i]
  canada_Browns$Indicator1[i] <- canada_Browns$Indicator2[i]
  canada_Browns$Indicator2[i] <- tmp
}

ggplot(canada_Browns, aes(x = Indicator1, y = Indicator2, fill = p.significant)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Rank), size = 4) +
  scale_x_discrete(drop = FALSE) +   # keep all levels
  scale_y_discrete(drop = FALSE) + 
  scale_fill_manual(
    values = c(
      "Not significant" = "#e0e0e0",  # light grey
      "Marginal" = "#b3cde3",         # light blue
      "Strong" = "#6497b1",           # medium blue
      "Very strong" = "#005b96"       # dark blue
    )
  ) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Significance",
    title = "Ranked Indicator Pair Relationships"
  )

#Plot USA Matrix
usa_Browns$Indicator1 <- factor(usa_Browns$Indicator1, levels = all_levels)
usa_Browns$Indicator2 <- factor(usa_Browns$Indicator2, levels = all_levels)
#some of the relations need to be swapped for plotting
rows_to_swap <- c(2, 8)  # for example

for (i in rows_to_swap) {
  tmp <- usa_Browns$Indicator1[i]
  usa_Browns$Indicator1[i] <- usa_Browns$Indicator2[i]
  usa_Browns$Indicator2[i] <- tmp
}

ggplot(usa_Browns, aes(x = Indicator1, y = Indicator2, fill = p.significant)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Rank), size = 4) +
  scale_x_discrete(drop = FALSE) +   # keep all levels
  scale_y_discrete(drop = FALSE) + 
  scale_fill_manual(
    values = c(
      "Not significant" = "#e0e0e0",  # light grey
      "Marginal" = "#b3cde3",         # light blue
      "Strong" = "#6497b1",           # medium blue
      "Very strong" = "#005b96"       # dark blue
    )
  ) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Significance",
    title = "Ranked Indicator Pair Relationships"
  )




