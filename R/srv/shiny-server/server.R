library(editData)
library(tidyverse)
library(psych)
library(MASS)
library(lmtest)
library(car)        # leveneTest
library(ggfortify)

library(sandwich)   # vcovHC
library(margins)    # margins

library(GGally)     # Scatter Plot in Correlation Analysis
library(ggcorrplot) # Corrplot in Correlation Analysis

library(survival)   # coxph in cox model
library(survminer)  # ggsurvplot in K-M estimates
library(boot)       # glm.diag()

library(ERSA)
library(prediction)
library(Hmisc)
library(dplyr)    
library(dplyrAssist)
library(lattice)
library(zoo)
library(xtable)
library(GauPro)

library(ggpubr)     # ggarrange

library(lavaan)     # SEM
library(semTools)   # SEM
library(semPlot)    # SEM
library(diagram)    # SEM
#library(Rcmdr)

pdf(NULL)

# data2 is main data frame in albatross
data2 <<- NULL

g_sheetindex <<- 1
g_env <<- new.env()
g_openfile <<- FALSE
g_resetresult <<- FALSE
g_datatype <<- NULL

source("module/chooser.R")
source("module/frailtyHLAR1.R")
source("module/jmfit.R")
source("module/jointfit_20190515.R")
source("module/dhglmjoint.R")
source("module/Module1_2.R")
source("module/VersionKorea3.R")
source("module/dhglmfit_20180603.R")
source("module/dhglmfit_run3.R")
source("module/ggplotAlbatross.R")
source("module/ggsem.R")

options(shiny.maxRequestSize=1024^3)
options(java.parameters = "-Xmx8192m")

# Password Rules
validate_pwd_custom <- function(pwd) {
  all(vapply(
    X = c("[0-9]+", "[a-z]+", "[A-Z]+", "[[:punct:]]+", ".{8,}"),
    FUN = grepl, x = pwd, FUN.VALUE = logical(1)
  ))
}

server <- function(input, output, session) {

    # call the server part
    # check_credentials returns a function to authenticate users
    res_auth <- secure_server(
        check_credentials = check_credentials(credentials)
    )

    output$auth_output <- renderPrint({
        reactiveValuesToList(res_auth)
    })

    # Default Apps
    source("logic/logic_dataImport.R", local = T)
    source("logic/logic_dataManagement.R", local = T)
    source("logic/logic_basicAnalysis.R", local = T)
    source("logic/logic_regression.R", local = T)
    source("logic/logic_randomEffects.R", local = T)
    source("logic/logic_survivalAnalysis.R", local = T)
    source("logic/logic_multiResponse.R", local = T)
    source("logic/logic_albatrossMaterials.R", local = T)

}

