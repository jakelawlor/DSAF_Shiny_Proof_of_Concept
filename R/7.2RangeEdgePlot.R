######################
#Range Edge Mapped------
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggplot2)
library(dplyr)
library(geosphere)
library(ggpmisc)

RE_DAT<- read.csv(here::here("R/DataforFinalFigs/Edge_df_NSreshp.csv"))
names(RE_DAT)
# Get high-res coastline
coast <- ne_coastline(scale = "large", returnclass = "sf")
######################
#test on subset----
# Isolate coordinates
Est_0.5 <- RE_DAT[, c("Year", "Estimate_km_E_quantile_0.5", "Estimate_km_N_quantile_0.5")]

# Convert from km to m
Est_0.5$Estimate_km_E_quantile_0.5 <- Est_0.5$Estimate_km_E_quantile_0.5 * 1000
Est_0.5$Estimate_km_N_quantile_0.5 <- Est_0.5$Estimate_km_N_quantile_0.5 * 1000

# Make sf object in local UTM (meters)
Est_0.5_pts <- st_as_sf(Est_0.5,
                        coords = c("Estimate_km_E_quantile_0.5", "Estimate_km_N_quantile_0.5"),
                        crs = "EPSG:32621")

# Reproject to WGS84
Est_0.5_pts2 <- st_transform(Est_0.5_pts, crs = 4326)

# Plot
ggplot() +
  geom_sf(data = coast, color = "gray50") +
  geom_sf(data = Est_0.5_pts2, color = "blue", size = 1) +
  #coord_sf(xlim = c(-80, -40), ylim = c(35, 60), expand = FALSE) +
  theme_minimal() +
  labs(title = "Sample Projected Coordinates in Northwest Atlantic",
       subtitle = "Reprojected from km-based custom CRS to WGS84")
######################
#function to reproject for plotting----
#separate each set of coordinates into their own objects
Est_0.05<-RE_DAT[ , c("Year","Estimate_km_E_quantile_0.05", "Estimate_km_N_quantile_0.05" )]
Est_0.5<-RE_DAT[ , c("Year", "Estimate_km_E_quantile_0.5",  "Estimate_km_N_quantile_0.5")] 
Est_0.95<-RE_DAT[ , c("Year", "Estimate_km_E_quantile_0.95",  "Estimate_km_N_quantile_0.95")] 
Std_0.05<-RE_DAT[ , c("Year", "Std_Dev_km_E_quantile_0.05",  "Std_Dev_km_N_quantile_0.05")]
Std_0.5<-RE_DAT[ , c("Year", "Std_Dev_km_E_quantile_0.5", "Std_Dev_km_N_quantile_0.5")] 
Std_0.95<-RE_DAT[ , c("Year", "Std_Dev_km_E_quantile_0.95", "Std_Dev_km_N_quantile_0.95")] 

#function
convert_km_to_wgs84 <- function(df, east_col, north_col) {
  # Convert from km to m
  df[[east_col]] <- df[[east_col]] * 1000
  df[[north_col]] <- df[[north_col]] * 1000
  # Make sf object in local CRS
  pts_local <- st_as_sf(df, coords = c(east_col, north_col), crs = "EPSG:32621")
  # Reproject to WGS84
  pts_wgs84 <- st_transform(pts_local, crs = 4326)
  return(pts_wgs84)
}

#Trailing edge Estimates
Est_0.05_pts <- convert_km_to_wgs84(
  Est_0.05,
  east_col = "Estimate_km_E_quantile_0.05",
  north_col = "Estimate_km_N_quantile_0.05"
)
#Mean Estimates
Est_0.5_pts <- convert_km_to_wgs84(
  Est_0.5,
  east_col = "Estimate_km_E_quantile_0.5",
  north_col = "Estimate_km_N_quantile_0.5"
)
#Leadind edge Estimates
Est_0.95_pts <- convert_km_to_wgs84(
  Est_0.95,
  east_col = "Estimate_km_E_quantile_0.95",
  north_col = "Estimate_km_N_quantile_0.95"
)
#Trailing edge Standard Deviation
Std_0.05_pts <- convert_km_to_wgs84(
  Std_0.05,
  east_col = "Std_Dev_km_E_quantile_0.05",
  north_col = "Std_Dev_km_N_quantile_0.05"
)
#Mean Standard Deviation
Std_0.5_pts <- convert_km_to_wgs84(
  Std_0.5,
  east_col = "Std_Dev_km_E_quantile_0.5",
  north_col = "Std_Dev_km_N_quantile_0.5"
)
#Leading edge Standard Deviation
Std_0.95_pts <- convert_km_to_wgs84(
  Std_0.95,
  east_col = "Std_Dev_km_E_quantile_0.95",
  north_col = "Std_Dev_km_N_quantile_0.95"
)
######################
# RangeEdge_Mapped: Plot Leading and Trailing edge on a map-------
#first run Sup2DataPlots.R up to line 39 to get the shapefiles 
Land <- ne_countries(scale = "large", returnclass = "sf")

RangeEdge_Mapped <-  ggplot() +
#start with mapping shapefiles 
  geom_sf(data = contours, color="lightblue", linewidth = .4) +
  geom_sf(data = NAFO, color="dimgrey", fill = NA) +
  geom_sf(data = Hague, color="black", linewidth = 1) +
  geom_sf(data = EEZ, color="black", linetype = "dashed", linewidth  = 1) +
  geom_sf(data = Land , color = "black", fill = "cornsilk") +
  #coord_sf(xlim = c(-72, -47), ylim = c(41, 50), expand = FALSE) 
  
#Add Trailing Edge
  # Path connecting points by year
  geom_path(data = Est_0.05_pts %>% 
              st_coordinates() %>% 
              as.data.frame() %>% 
              cbind(Year = Est_0.05_pts$Year) %>% 
              arrange(Year),
            aes(x = X, y = Y, color = Year),
            linewidth = 0.8) +
  # Regular points for 0.05 quantile
  geom_sf(data = Est_0.05_pts, 
          aes(color = Year),
          na.rm = TRUE, 
          alpha = 1, 
          size = 2, 
          shape = 19)+
  # Larger terminal point for 0.05 quantile (most recent year)
  geom_sf(data = Est_0.05_pts[Est_0.05_pts$Year == max(Est_0.05_pts$Year, na.rm = TRUE), ],
             aes(color = Year),
             #position = pos_jitter, 
             na.rm = TRUE, 
             alpha = 1, 
             size = 5, 
             shape = 19) +
  # Annotation for 0.05 quantile terminal point
  annotate("text",
           x=-64.5, y=42.1,
           label = "Trailing Edge",
           color = "black",
           size = 5,
           family = "serif") +
#add Leading Edge 
  # Path connecting points by year
  geom_path(data = Est_0.95_pts %>% 
              st_coordinates() %>% 
              as.data.frame() %>% 
              cbind(Year = Est_0.95_pts$Year) %>% 
              arrange(Year),
            aes(x = X, y = Y, color = Year),
            linewidth = 0.8) +
  # Regular points for 0.95 quantile
  geom_sf(data = Est_0.95_pts, 
          aes(color = Year),
          na.rm = TRUE, 
          alpha = 1, 
          size = 2, 
          shape = 19)+
  # Larger terminal point for 0.95 quantile (most recent year)
  geom_sf(data = Est_0.95_pts[Est_0.95_pts$Year == max(Est_0.95_pts$Year, na.rm = TRUE), ],
          aes(color = Year),
          #position = pos_jitter, 
          na.rm = TRUE, 
          alpha = 1, 
          size = 5, 
          shape = 19) + 
  # Annotation for 0.95 quantile terminal point  
  annotate("text",
           x=-53, y=45,
           label = "Leading Edge",
           color = "black",
           size = 5,
           family = "serif") +
# Blue to orange gradient for Year
  scale_color_gradient2(low ="#005A9C",       ##"#0072B2", "#003366",  #"steelblue3"
                        mid = "#FFD700",      #yellow", 
                        high =  "orangered",  #"#D55E00",
                        midpoint = 2005, 
                        name = "Year") + 
# Labels and plot range
  xlab("Range Edge W to E") +
  ylab("Range Edge S to N") +
  coord_sf(xlim = c(-72, -47), ylim = c(41, 50), expand = FALSE) +
# Theme
  theme_bw() +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 10, 10, 10),
    axis.text.x = element_text(size = 14, 
                               family = "serif"),
    axis.text.y = element_text(size = 14, 
                               family = "serif"),
    axis.title.x = element_text(size = 14, 
                                hjust = 0.5, 
                                vjust = -2, 
                                family = "serif"),
    axis.title.y = element_text(size = 14, 
                                hjust = 0.5, 
                                vjust = 4, 
                                angle = 90, 
                                family = "serif"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, family = "serif"),
    strip.text = element_text(size = 14, 
                              family = "serif", 
                              angle = 0),
    strip.background = element_rect(colour = "black", 
                                    fill = "white"),
    panel.grid.minor = element_blank(), panel.grid.major =  element_blank(),
  )

RangeEdge_Mapped
######################



 
  

  

  

  
