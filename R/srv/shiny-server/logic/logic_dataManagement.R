source("logic/DataManagement/logic_dmMake.R", local = T)
source("logic/DataManagement/logic_dmMakeInterval.R", local = T)
source("logic/DataManagement/logic_dmConvert.R", local = T)
source("logic/DataManagement/logic_dmSelect.R", local = T)
source("logic/DataManagement/logic_dmRename.R", local = T)
source("logic/DataManagement/logic_dmDelete.R", local = T)
source("logic/DataManagement/logic_dmMerge.R", local = T)

# Tab Change ####

observeEvent(input$dm_inputtabset, {
    dmSelectetab = NULL
    if (input$dm_inputtabset == "dm_mv_tab") {
        dmSelectetab = "dm_dftab"
    }  else if (input$dm_inputtabset == "dm_ct_tab") {
        dmSelectetab = "dm_datab"
    } else if (input$dm_inputtabset == "dm_sr_tab") {
        dmSelectetab = "dm_srptab"
    } 
    updateTabsetPanel(session, "dm_resulttabset", selected = dmSelectetab)
})

# Data Handling ####

output$dm_datatitle <- renderText({
    paste("Data Summary : ", strong(input$file$name))
})

output$dm_dataframe <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_mvi_run
    input$dm_ct_run
    input$dm_sr_run  
    input$dm_rn_run
    input$dm_del_run
    input$dm_md_run
    
    head(data2, 20)
}, bordered = TRUE)

dataAttribute <- reactive({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_mvi_run
    input$dm_ct_run
    input$dm_sr_run  
    input$dm_rn_run
    input$dm_del_run
    input$dm_md_run
    
    attributeTable<-NULL
    dataAttributes<-sapply(data2, class)
    attributeTable$Name<-names(data2)
    attributeTable$Type<-dataAttributes
    attributeTable
})

output$dm_dataattribute <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_mvi_run
    input$dm_ct_run
    input$dm_sr_run  
    input$dm_rn_run
    input$dm_del_run
    input$dm_md_run
    
    if(is.null(dataAttribute()))
        return()
    
    dataAttribute()
}, bordered = TRUE)

output$dm_selectrowspreview <- renderTable({
    input$resetData
    input$di_option_run
    input$dm_sr_preview
    
    if (is.null(g_sr_previewdata)) {
        return()
    }
    g_sr_previewdata
}, bordered = TRUE)

# Merge Dataset ####

output$dm_mergetitle <- renderText({
    paste("Merge", strong(input$md_file$name), "into", strong(input$file$name))
})

output$dm_firstdataframe <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_mvi_run
    input$dm_ct_run
    input$dm_sr_run  
    input$dm_rn_run
    input$dm_del_run
    input$dm_md_run
    
    
    head(data2, 10)
}, caption = "Main Dataset", bordered = TRUE, caption.placement = getOption("xtable.caption.placement", "top"))

output$dm_seconddataframe <- renderTable({
    input$md_file
    
    head(mergedata(), 10)
}, caption = "Second Dataset", bordered = TRUE, caption.placement = getOption("xtable.caption.placement", "top"))

output$dm_mergepreview <- renderTable({
    input$resetData
    input$di_option_run
    input$dm_md_preview
    
    if (is.null(g_md_previewdata)) {
        return()
    }
    g_md_previewdata
})