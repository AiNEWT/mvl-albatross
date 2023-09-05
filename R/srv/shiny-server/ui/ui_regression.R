navbarMenu(
    "Regression",
    source("ui/Regression/ui_rgLm.R", local = TRUE)[[1]],
    source("ui/Regression/ui_rgGlm.R", local = TRUE)[[1]],
    source("ui/Regression/ui_rgLogit.R", local = TRUE)[[1]],
    source("ui/Regression/ui_rgProbit.R", local = TRUE)[[1]],
    source("ui/Regression/ui_rgLoglinear.R", local = TRUE)[[1]],
    source("ui/Regression/ui_rgJglm.R", local = TRUE)[[1]]
)
