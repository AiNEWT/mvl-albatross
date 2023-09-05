navbarMenu(
    "Multiple Response Analysis",
    source('ui/MultiResponse/ui_mrMdhglm.R', local = TRUE)[[1]],
    source('ui/MultiResponse/ui_mrJointModel.R', local = TRUE)[[1]],
    source('ui/MultiResponse/ui_mrFactor.R', local = TRUE)[[1]],
    source('ui/MultiResponse/ui_mrSem.R', local = TRUE)[[1]],
#    source('ui/MultiResponse/ui_mrDSem.R', local = TRUE)[[1]],
    source('ui/MultiResponse/ui_mrDSem2.R', local = TRUE)[[1]]
)
