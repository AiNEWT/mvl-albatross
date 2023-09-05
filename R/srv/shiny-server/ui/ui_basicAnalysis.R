navbarMenu(
    "Basic Analysis",
    source("ui/BasicAnalysis/ui_baDescriptive.R", local = TRUE)[[1]],
    source("ui/BasicAnalysis/ui_baTtest.R", local = TRUE)[[1]],
    source("ui/BasicAnalysis/ui_baAnova.R", local = TRUE)[[1]],
    source("ui/BasicAnalysis/ui_baFrequency.R", local = TRUE)[[1]],
    source("ui/BasicAnalysis/ui_baCorrelation.R", local = TRUE)[[1]]
)
