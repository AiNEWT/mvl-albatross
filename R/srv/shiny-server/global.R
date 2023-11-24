# define some basic credentials (on data.frame)
credentials <- data.frame(
  user = c("albatross", "admin"), # mandatory
  password = c("albatross", "admin"), # mandatory
  start = c("2023-11-01"), # optinal (all others)
  expire = c(NA, "2019-12-31"),
  admin = c(FALSE, TRUE),
#  comment = "Simple and secure authentification mechanism for single ‘Shiny’ applications.",
  stringsAsFactors = FALSE
)

if (!require(shiny)) {install.packages("shiny"); library(shiny)}
if (!require(shiny)) {install.packages("shinymanager"); library(shinymanager)}

source("./ui.R", local = TRUE)  
source("./server.R", local = TRUE)  

shinyApp(ui, server)
