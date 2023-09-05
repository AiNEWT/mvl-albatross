library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinydashboard)
library(shinycssloaders)
library(rintrojs)
library(shinyBS)
library(bsplus)

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