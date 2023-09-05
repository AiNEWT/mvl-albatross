navbarMenu(
    "Survival Analysis",
    source("ui/SurvivalAnalysis/ui_saKM.R", local = TRUE)[[1]],
    source("ui/SurvivalAnalysis/ui_saCox.R", local = TRUE)[[1]],
    source("ui/SurvivalAnalysis/ui_saFM.R", local = TRUE)[[1]],
    source("ui/SurvivalAnalysis/ui_saCR.R", local = TRUE)[[1]]
)