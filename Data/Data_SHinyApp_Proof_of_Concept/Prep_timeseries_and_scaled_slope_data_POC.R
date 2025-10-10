
#Data handover/prep for Jake Lawlor
library(dplyr)
library(ggplot2)
library(patchwork)
#Bring in all the data and pair it down until it contains only what is needed for handover----

#SECTION 1: timeseries data
#SECTION 2: scaled slopes

#SECTION 1: timeseries data
#1: ABUNDANCE----
  #Data from 3.1 Data_prep.R, step 1: Generated stratified abundance and standard error estimates
  abundance_ind_Region<-read.csv(here::here("2025-04-23/Output/IndexAbundance/abundance_ind_Region.csv"))  
  str(abundance_ind_Region)

  #subset to Spring and remove "all" (Canada and USA only)
abundance_ind_Region_Sp<-subset(abundance_ind_Region, abundance_ind_Region$Season=="Spring" & abundance_ind_Region$Index_Region=="Canada" | 
                                  abundance_ind_Region$Season=="Spring" & abundance_ind_Region$Index_Region=="USA")

ggplot(abundance_ind_Region_Sp, 
       aes(x = Year, y = Index_Estimate, colour = Index_Region, group = Index_Region)) +
  geom_line() +
  geom_ribbon(aes(ymin = Index_Estimate - Index_SD,
                  ymax = Index_Estimate + Index_SD,
                  fill = Index_Region), 
              alpha = 0.25, colour = NA) +
  theme_bw() +
  labs(x = "Year", y = "Index Estimate",
       colour = "Region", fill = "Region")

head(abundance_ind_Region_Sp)
#remove the columns that are not useful
  abundance_ind_Region_Sp$Time <- NULL
  abundance_ind_Region_Sp$Category  <- NULL
  abundance_ind_Region_Sp$YearGroup  <- NULL
  abundance_ind_Region_Sp$Date  <- NULL
  abundance_ind_Region_Sp$Season  <- NULL
#rename some columns for clarity
  names(abundance_ind_Region_Sp)[c(2,3,4)] <- c("Region", "Estimate", "SD")
#reset X column 
  abundance_ind_Region_Sp$X <- 1:102
  str(abundance_ind_Region_Sp)
#save
write.csv(abundance_ind_Region_Sp,(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_Abundance.csv")), row.names = FALSE)

#2 AREA OCCUPIED----
  #Data from 3.1 Data_prep.R, steps 2&3: 
    #get the abundance estimates per grid location, Add season, Year, depth, and the area (km2) of the Stratum
      #then 8.1 Deepening.R adds a depth field
      abdest<- read.csv(here::here("2025-04-23/Output/IndexAbundance/ForShiftAnalysis/AbundanceEstimates_GridCentriods_Reg_wDepth.csv"))
      dim(abdest) 
      str(abdest) 
      abdest$Region<-factor(abdest$Stratum)

#Select for Spring
abdest.spr<- abdest %>%
  filter(Season == "Spring")
str(abdest.spr)

# Calculate total annual spring abundance by region and year
abdest.spr <- abdest.spr %>%
  group_by(Region, Year) %>%
  mutate(TotalSpringAbundance = sum(Abundance))

## Calculating Area Thresholds
#We'll calculate the area containing different percentages of abundance (50%, 75%, 90%, 95%) for each region and year:
# Function to calculate area for different abundance thresholds 
#from EAOrmarkdown.Rmd
calculate_areas <- function(data, thresholds = c(50, 75, 90, 95)) {
  result_list <- list()
  
  for (threshold in thresholds) {
    threshold_result <- data %>%
      group_by(Region, Year) %>%
      mutate(Total_Abundance = sum(Abundance)) %>%
      mutate(Total_Area = sum(Area_km2)) %>%#
      arrange(Region, Year, desc(Abundance/Area_km2)) %>%
      mutate(
        Cumulative_Abundance = cumsum(Abundance),
        Percent_Abundance = Cumulative_Abundance / Total_Abundance * 100,
        Cumulative_Area = cumsum(Area_km2)
      ) %>%
      filter(Percent_Abundance <= threshold) %>%
      summarize(
        Threshold = threshold,
        Area_Threshold = sum(Area_km2),
        Total_Area = first(Total_Area),#
        Percent_Area_Used = Area_Threshold / Total_Area * 100,#
        Total_Abundance = first(Total_Abundance),
        #Total_Area = sum(Area_km2, na.rm = TRUE),
        #Percent_Area_Used = Area_Threshold / sum(Area_km2, na.rm = TRUE) * 100,
        n_cells = n()
      )
    
    result_list[[as.character(threshold)]] <- threshold_result
  }
  
  # Combine all results
  bind_rows(result_list)
}

# Apply the function to calculate areas for different thresholds
area_thresholds <- calculate_areas(abdest.spr)
summary(area_thresholds)
area_thresholds <- area_thresholds %>%
  mutate(Area_Efficiency = Area_Threshold / Total_Abundance)


#Plot Area Occupied: The total area containing 90%,of the abundance
#subset 90%
area_thresholds_sub<-subset(area_thresholds, area_thresholds$Threshold ==90)
AO1<-ggplot(area_thresholds_sub, 
       aes(x = Year, y = Area_Threshold, colour = Region, group = Region)) +
  geom_line() +
  theme_bw() +
  labs(x = "Year", y = "Total Area occupied (90% Abundance)",
       colour = "Region", fill = "Region")+
  theme(legend.position = "none")

AO2<-ggplot(area_thresholds_sub, 
       aes(x = Year, y = Percent_Area_Used  , colour = Region, group = Region)) +
  geom_line() +
  theme_bw() +
  labs(x = "Year", y = "Percent of Area Used (90% Abundance)",
       colour = "Region", fill = "Region")+
  theme(legend.position = "none")

AO3<-ggplot(area_thresholds_sub, 
            aes(x = Year, y =Area_Efficiency , colour = Region, group = Region)) +
  geom_line() +
  theme_bw() +
  labs(x = "Year", y = "Area Efficiency (90% Abundance)",
       colour = "Region", fill = "Region")
AO1 + AO2 + AO3 + plot_layout(widths = c(1,1,1.2))

#save
write.csv(area_thresholds,(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_AreaOccupied.csv")), row.names = FALSE)

#3 ABUNDANCE-WEIGHTED DEPTH----
  #Data from 8.1Deepening.R:
    #Depth field is added to the abundance estimates per grid location (AbundanceEstimates_GridCentriods_Reg_wDepth.csv),
    # Data are grouped by year/season/region
      #calculates the mean, median, Q5, and Q95 depth, weighted by estimated abundance values 
      D_data_Reg<- read.csv(here::here("2025-04-23/Output/Shift_Indicators/Seasonal_Deepening_Reg.csv"))
      str(D_data_Reg) 

#subset spring
  D_data_Reg<- subset(D_data_Reg, D_data_Reg$Season == "Spring")
  names(D_data_Reg)[names(D_data_Reg) == "Stratum"] <- "Region" #rename some columns for clarity
  D_data_Reg$Period  <- NULL #remove unneeded field
  D_data_Reg$Season  <- NULL #remove unneeded field
#these should all be negative values 
  D_data_Reg$Depth_Mean<- D_data_Reg$Depth_Mean*-1
  D_data_Reg$Depth_Median<- D_data_Reg$Depth_Median*-1
  D_data_Reg$Depth_Q5<- D_data_Reg$Depth_Q5*-1
  D_data_Reg$Depth_Q95<- D_data_Reg$Depth_Q95*-1

#plot
ggplot(D_data_Reg, 
       aes(x = Year, y = Depth_Median , colour = Region, group=Region)) +
  geom_line() +
  geom_ribbon(aes(ymin = Depth_Q5,
                  ymax = Depth_Q95,
                  fill = Region), 
              alpha = 0.25, colour = NA) +
  theme_bw() +
  labs(x = "Year", y = "Depth",
       colour = "Region", fill = "Region")
write.csv(D_data_Reg,(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_AWD.csv")), row.names = FALSE)

#4 RANGE EDGE----
  #from 7.1RangeEdge.R: Calculates 5th/50th/95th percentile of the spatial distribution (Weighted by abundance est,  quantile of the coordinate values)
  RE_DAT<- read.csv(here::here("R/DataforFinalFigs/Edge_df_NSreshp.csv"))
  str(RE_DAT)
  
#these data area already subset to spring and represent the whole study area, but we can remove some useless fields
RE_DAT$YearGroup   <- NULL
RE_DAT$Season   <- NULL
write.csv(RE_DAT,(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_RangeEdge.csv")), row.names = FALSE)

#5 CENTRE OF GRAVITY----
  #data from: 5.1Centre_of Gravity.R, 
    #which takes the abundance estimates per grid location data (AbundanceEstimates_GridCentriods_Reg_wDepth.csv), from 3.1 Data_prep.R
      #and calculates the weighted mean of lon/lat based on the abundance values,
      centroid_data_Reg<- read.csv(here::here("2025-04-23/Output/Shift_Indicators/seasonal_centroid_data_region.csv"))
      str(centroid_data_Reg)
      
#subset for spring,rename some columns for clarity, and remove unneeded
  centroid_data_Reg<-subset(centroid_data_Reg, centroid_data_Reg$Season=="Spring")
  names(centroid_data_Reg)[names(centroid_data_Reg) == "Stratum"] <- "Region" 
  centroid_data_Reg$Season  <- NULL 
  str(centroid_data_Reg)
  write.csv(centroid_data_Reg,(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_COG.csv")), row.names = FALSE)


#6 DISTANCE TO SHARED BORDER----
  #data from 6.1Distance_From_Hague.R
    #identifies the closest point on the border to each COG in the time series for each season/Region grouping,
      #before calculating the distance (km) between each "closest" point and the corresponding COG 
      dist_hague_Reg<- read.csv(here::here("2025-04-23/Output/Shift_Indicators/dist_hague_Reg_seasonal.csv"))
      str(dist_hague_Reg)

#subset for spring,rename some columns for clarity, and remove unneeded
      dist_hague_Reg<-subset(dist_hague_Reg, dist_hague_Reg$Season=="Spring")
      names(dist_hague_Reg)[names(dist_hague_Reg) == "Stratum"] <- "Region" 
      dist_hague_Reg$Season  <- NULL       
      str(dist_hague_Reg)
#multiply USA by -1 so they are on the left side of the border
      dist_hague_Reg2 <- dist_hague_Reg %>%
        mutate(
          Dist_Mean = if_else(Region == "USA", Dist_Mean * -1, Dist_Mean),
          Dist_Med = if_else(Region == "USA", Dist_Med * -1, Dist_Med),
          Dist_Q5 = if_else(Region == "USA", Dist_Q5 * -1, Dist_Q5),
          Dist_Q95 = if_else(Region == "USA", Dist_Q95 * -1, Dist_Q95)
        )
#plot      
ggplot(dist_hague_Reg2, 
  aes(x = Year, y = Dist_Mean, colour = Region, group = Region)) +
 geom_line() +
        theme_bw() +
        labs(x = "Year", y = "Distance to Border (km)",
             colour = "Region", fill = "Region")  +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8)

write.csv(dist_hague_Reg2,(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_DtoB.csv")), row.names = FALSE)
#End section1----


#SECTION 2: scaled slopes----
#read in the timeseries data
POC_Abundance<-read.csv(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_Abundance.csv"))
POC_AO<-read.csv(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_AreaOccupied.csv"))
POC_AWD<-read.csv(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_AWD.csv"))
POC_RE<- read.csv(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_RangeEdge.csv"))
POC_COG<-read.csv(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_COG.csv"))
POC_DtoB<-read.csv(here::here("Data/Data_SHinyApp_Proof_of_Concept/POC_DtoB.csv"))
# We are only calculating the slope on the period since 2006 when the water began warming 
#add the time period info 
add_period <- function(df) {
  df$Period <- NULL
  df$Period[df$Year<2006]<-"Before Warming"
  df$Period[df$Year>2005]<-"During Warming"
  return(df)
}
POC_Abundance<-add_period(POC_Abundance)
POC_AO<-add_period(POC_AO)
POC_AWD<-add_period(POC_AWD)
POC_RE<-add_period(POC_RE)
POC_COG<-add_period(POC_COG)
POC_DtoB<-add_period(POC_DtoB)

#calculate scaled slopes by Region and Period
library(dplyr)
library(broom)

#ABUNDANCE
#estimates scaled change in estimated abundance per year
str(POC_Abundance)
Abundance_coefficients_df <- POC_Abundance %>%
  group_by(Region,Period) %>%
  do({
    model <- lm(scale(Estimate) ~ scale(Year), data = .)
    data.frame(t(coef(model)))
    tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
  }) %>%
  ungroup()
Abundance_coefficients_df<-Abundance_coefficients_df%>%
  filter(term == "scale(Year)", Period == "During Warming") %>%
  mutate(Indicator = "Abundance")
Abundance_coefficients_df

#ARE OCCUPIED
#estimates scaled change in the total area  required to reach abundance Threshold==90 per year
str(POC_AO)
POC_AO_sub<-subset(POC_AO, POC_AO$Threshold ==90)
AO_coefficients_df <- POC_AO_sub %>%
  group_by(Region,Period) %>%
  do({
    model <- lm(scale(Area_Threshold) ~ scale(Year), data = .)
    data.frame(t(coef(model)))
    tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
  }) %>%
  ungroup()
AO_coefficients_df<-AO_coefficients_df%>%
  filter(term == "scale(Year)", Period == "During Warming") %>%
  mutate(Indicator = "Area Occupied")
AO_coefficients_df

#ABUNDANCE-WEIGHTED DEPTH
#estimated scaled change in average depth per year
str(POC_AWD)
AWD_coefficients_df <- POC_AWD %>%
  group_by(Region,Period) %>%
  do({
    model <- lm(scale(Depth_Mean) ~ scale(Year), data = .)
    data.frame(t(coef(model)))
    tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
  }) %>%
  ungroup()
AWD_coefficients_df<-AWD_coefficients_df%>%
  filter(term == "scale(Year)", Period == "During Warming") %>%
  mutate(Indicator = "Abundance Weighted Depth")
AWD_coefficients_df

#RANGE EDGE
str(POC_RE)
  #Trailing Edge(located in USA)
    #Trailing edge East
      RE_Trail_East_coefficients_df <- POC_RE %>%
        group_by(Period) %>%
        do({
          model <- lm(scale(Estimate_km_E_quantile_0.05) ~ scale(Year), data = .)
          data.frame(t(coef(model)))
          tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
        }) %>%
        ungroup()
      RE_Trail_East_coefficients_df<-RE_Trail_East_coefficients_df%>%
        filter(term == "scale(Year)", Period == "During Warming") %>%
        mutate(Indicator = "Trailing_Edge_East",
                Region = "USA")%>%
        select(Region, everything())
      RE_Trail_East_coefficients_df
      
      #Trailing edge North
      RE_Trail_North_coefficients_df <- POC_RE %>%
        group_by(Period) %>%
        do({
          model <- lm(scale(Estimate_km_N_quantile_0.05) ~ scale(Year), data = .)
          data.frame(t(coef(model)))
          tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
        }) %>%
        ungroup()
      RE_Trail_North_coefficients_df<-RE_Trail_North_coefficients_df%>%
        filter(term == "scale(Year)", Period == "During Warming") %>%
        mutate(Indicator = "Trailing_Edge_North",
               Region = "USA")%>%
        select(Region, everything())
      RE_Trail_North_coefficients_df
      
  #Leading Leading(located in Canada)
    #Leading edge East
      RE_Lead_East_coefficients_df <- POC_RE %>%
        group_by(Period) %>%
        do({
          model <- lm(scale(Estimate_km_E_quantile_0.95) ~ scale(Year), data = .)
          data.frame(t(coef(model)))
          tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
        }) %>%
        ungroup()
      RE_Lead_East_coefficients_df<-RE_Lead_East_coefficients_df%>%
        filter(term == "scale(Year)", Period == "During Warming") %>%
        mutate(Indicator = "Leading_Edge_East",
                Region = "Canada")%>%
        select(Region, everything())
      RE_Lead_East_coefficients_df
      
      #Leading edge North 
      RE_Lead_North_coefficients_df <- POC_RE %>%
        group_by(Period) %>%
        do({
          model <- lm(scale(Estimate_km_N_quantile_0.95) ~ scale(Year), data = .)
          data.frame(t(coef(model)))
          tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
        }) %>%
        ungroup()
      RE_Lead_North_coefficients_df<-RE_Lead_North_coefficients_df%>%
        filter(term == "scale(Year)", Period == "During Warming") %>%
        mutate(Indicator = "Leading_Edge_North",
               Region = "Canada")%>%
        select(Region, everything())
      RE_Lead_North_coefficients_df
    
#CENTRE OF GRAVITY
str(POC_COG)
  #east
  COG_E_coefficients_df <- POC_COG %>%
    group_by(Region,Period) %>%
    do({
      model <- lm(scale(centroid_longitude) ~ scale(Year), data = .)
      data.frame(t(coef(model)))
      tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
    }) %>%
    ungroup()
  COG_E_coefficients_df<-COG_E_coefficients_df%>%
    filter(term == "scale(Year)", Period == "During Warming") %>%
    mutate(Indicator = "COG_East")
  COG_E_coefficients_df
  #north
  COG_N_coefficients_df <- POC_COG %>%
    group_by(Region,Period) %>%
    do({
      model <- lm(scale(centroid_latitude) ~ scale(Year), data = .)
      data.frame(t(coef(model)))
      tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
    }) %>%
    ungroup()
  COG_N_coefficients_df<-COG_N_coefficients_df%>%
    filter(term == "scale(Year)", Period == "During Warming") %>%
    mutate(Indicator = "COG_North")
  COG_N_coefficients_df

#DISTANCE TO SHARED BORDER 
str(POC_DtoB)
DtoB_coefficients_df <- POC_DtoB %>%
  group_by(Region,Period) %>%
  do({
    model <- lm(scale(Dist_Mean) ~ scale(Year), data = .)
    data.frame(t(coef(model)))
    tidy(model, conf.int = TRUE) # Includes coefficients with 95% CI by default
  }) %>%
  ungroup()
DtoB_coefficients_df<-DtoB_coefficients_df%>%
  filter(term == "scale(Year)", Period == "During Warming") %>%
  mutate(Indicator = "Distance to Border")
DtoB_coefficients_df

#the dfs all have the same structure so merge them into one
merged_coefficients_df <- bind_rows(Abundance_coefficients_df, 
                       AO_coefficients_df, 
                       RE_Lead_East_coefficients_df, RE_Lead_North_coefficients_df, 
                       RE_Trail_East_coefficients_df, RE_Trail_North_coefficients_df,
                       COG_E_coefficients_df, COG_N_coefficients_df,
                       DtoB_coefficients_df, 
                       AWD_coefficients_df)

merged_coefficients_df$Period<-"2006-2023"
merged_coefficients_df<-merged_coefficients_df%>%
  select(Indicator, everything())

merged_coefficients_df <- merged_coefficients_df %>%
  mutate(
    p.significant = case_when(
      p.value < 0.001 ~ "Very strong",
      p.value < 0.01  ~ "Strong",
      p.value < 0.05  ~ "Moderate",
      TRUE ~ "Not significant"
    )
  )
 str(merged_coefficients_df)
 merged_coefficients_df$p.significant <- factor(merged_coefficients_df$p.significant, 
                          levels = c("Very strong", "Strong","Moderate",  "Not significant"))
 merged_coefficients_df<-read.csv(here::here("Data/Data_SHinyApp_Proof_of_Concept/Scaled_Slopes_all_indicators.csv"))
 
 #plot
 ggplot(merged_coefficients_df, 
        aes(x = estimate, 
            y = Indicator, 
            colour = Region,
            size = p.significant)) +
   geom_point(position = position_dodge(width = 0.6)) +
   scale_size_manual(
     values = c(
       "Very strong" = 4,   
       "Strong" = 3,
       "Moderate" = 2,
       "Not significant" = 1
     )
   ) +
   geom_errorbar(aes(xmin = conf.low, xmax = conf.high), 
                 width = 0.2, 
                 position = position_dodge(width = 0.6),
                 linewidth = 0.6) +  # fixed bar thickness
   geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
   theme_minimal(base_size = 14) +
   labs(
     x = "Estimated slope (±95% CI)",
     y = "",
     colour = "Region",
     size = "Significance",
     title = "2006–2023 Trends"
   ) +
   theme(
     panel.grid.minor = element_blank(),
     axis.title.y = element_text(margin = margin(r = 10))
   )
 
write.csv(merged_coefficients_df,(here::here("Data/Data_SHinyApp_Proof_of_Concept/Scaled_Slopes_all_indicators.csv")), row.names = FALSE)
