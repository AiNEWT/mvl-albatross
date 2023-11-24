if (!require(shiny)) {install.packages("shiny"); library(shiny)}
if (!require(shiny)) {install.packages("shinymanager"); library(shinymanager)}

source("./ui.R", local = TRUE)  
source("./server.R", local = TRUE)  

shinyApp(ui, server)
