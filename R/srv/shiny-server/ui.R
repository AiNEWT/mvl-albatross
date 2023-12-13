library(shiny)
library(shinymanager)
library(shinythemes)
library(shinyWidgets)
library(shinydashboard)
library(shinycssloaders)
library(rintrojs)
library(shinyBS)
library(bsplus)

ui <- fluidPage(
  # auth_ui(
  #   id = "auth",
  #   # add image on top ?
  #   tags_top = 
  #     tags$div(
  #       tags$h4("Albatross Analytics", style = "align:center"),
  #       tags$img(
  #         src = "https://www.r-project.org/logo/Rlogo.png", width = 100
  #     )
  #   ),
  #   # add information on bottom ?
  #   tags_bottom = tags$div(
  #     tags$p(
  #       "For any question, please contact ",
  #       tags$a(
  #         href = "mailto:info@ainewt.ai?Subject=Albatross%20aAnalytics",
  #         target="_top", "info"
  #       )
  #     )
  #   ),
  #   background = NULL,
  #   lan = use_language("en")
  # ),
  # authentication module
  verbatimTextOutput("auth_output")
)

ui <- navbarPage(
    title = "Albatross Analytics",
    # theme = shinytheme("simplex"),
    theme = shinytheme("flatly"),

    source('ui/ui_dataImport.R', local = TRUE)[[1]],
    source("ui/ui_dataManagement.R", local = TRUE)[[1]],
    source("ui/ui_basicAnalysis.R", local = TRUE)[[1]],
    source("ui/ui_regression.R", local = TRUE)[[1]],
    source("ui/ui_randomEffects.R", local = TRUE)[[1]],
    source("ui/ui_survivalAnalysis.R", local = TRUE)[[1]],
    source("ui/ui_multiResponse.R", local = TRUE)[[1]],
    source("ui/ui_albatrossMaterials.R", local = TRUE)[[1]]
    # source("ui/ui_albatrossQA.R", local = TRUE)[[1]]
)

ui <- secure_app(ui)
