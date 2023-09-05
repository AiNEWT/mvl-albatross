navbarMenu(
    "Random Effect Model",
    source("ui/RandomEffects/ui_reLmm.R", local = TRUE)[[1]],
    source("ui/RandomEffects/ui_reGlmm.R", local = TRUE)[[1]],
    source("ui/RandomEffects/ui_reHglm.R", local = TRUE)[[1]],
    source("ui/RandomEffects/ui_reDhglm.R", local = TRUE)[[1]],
    source("ui/RandomEffects/ui_reCubicSpline.R", local = TRUE)[[1]]
)