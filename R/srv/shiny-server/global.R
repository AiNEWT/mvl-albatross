if (!require(shiny)) {install.packages("shiny"); library(shiny)}

source("./ui.R", local = TRUE)  
source("./server.R", local = TRUE)  

shinyApp(ui, server)
