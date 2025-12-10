# Global options for Transboundary Jusisdictional Framework Shiny



# libraries ---------------------------------------------------------------
# load the following packages
library(shiny)
library(bs4Dash)
library(dplyr)
library(readr)
library(plotly)
library(leaflet)
library(DT)
library(fresh)
library(shinyWidgets)
library(scales)
library(reactable)


# source UI tabs ----------------------------------------------------------
# The content for each tab is stored in a separate file. 
# Source all .R files in the current directory that start with "ui_":
sapply(list.files(
  pattern = "^ui_.*\\.R$",
  path = ".",
  full.names = TRUE
),
source)



# source text -------------------------------------------------------------
# script containing all app text in a list, which is called 
# for every text entry on various pages
source("app_text.R")


# set theme ---------------------------------------------------------------
# make some minimal theme changes - color of primary and secondary items
theme <- create_theme(
  bs4dash_status(
    primary = "#E1EDED",
    secondary = "#a7effa"
  )
)


# upload premade plots for spp X -----------------------------------------------------------
sppx_plots <- list(
  p_pie_Region1 = readRDS("www/premade_plots/pie_Region1.rds"),
  p_pie_Region2 = readRDS("www/premade_plots/pie_Region2.rds"),
  p_depth =  readRDS("www/premade_plots/depth_plot.rds"),
  p_abundance = readRDS("www/premade_plots/abundance_plot.rds"),
  p_area = readRDS("www/premade_plots/area_plot.rds"),
  p_map = readRDS("www/premade_plots/mapplot_NEW.rds"),
  p_dist = readRDS("www/premade_plots/distance_to_border_plot.rds"),
  p_slopes = readRDS("www/premade_plots/slope_plot.rds"),
  Region1_table = readRDS("www/premade_plots/Region1_table.rds"),
  Region2_table = readRDS("www/premade_plots/Region2_table.rds"),
  p_pie_legend = readRDS("www/premade_plots/pie_legend.rds")
)


# trend function for circle plots --------------------------------------------
trend_pal <- col_numeric(
  c("tomato2", "white", "cornflowerblue"),
  domain = c(-1, 1)
)


# upload nafo zone names for download button choices ---------------------
nafo_zone_names <- readRDS("www/nafo_zone_names.rds")



# upload study area plot --------------------------------------------------
study_area_imgs <- list(
  jurisdictional = "premade_plots/study_area.png"
)


