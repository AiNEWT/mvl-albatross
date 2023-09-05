mergedata <- reactive({
    if(is.null(input$md_file)) {return()}
    md_data = read.table(file=input$md_file$datapath, sep=input$md_sep, header=input$md_header, stringsAsFactors = TRUE,fileEncoding='EUC-Kr')
#    for (i in 1:length(md_data)) {
#        names(md_data)[i] <- paste0(names(md_data[i]), ".1")
#    }
    return(md_data)
})

# Excel Type File Sheet Setting ####

output$dm_md_sheetindex <- renderUI({
    fluidRow(
        column(
            4,
            selectInput(
                "dm_md_sheetindex",
                "Sheet Index", 
                choices = 1:5,
                #choices = (1:l_sheetnum), 
                multiple = FALSE
            )
        ),
        column(
            2,
            actionButton("dm_md_sheetindex_run", "Choose", icon = icon("play")),
            style = "text-align:left; padding:15px"
        )
    )
})

# Merge Type Setting ####

output$dm_md_mainall <- renderUI({
    input$md_file
    
    md_data <- mergedata()
    nameValue = names(md_data)
    
    selectInput(
        "dm_md_mainall",
        "All Observarions in the Main Datasets Included",
        choices = as.list(c(FALSE,TRUE)), 
        selected = FALSE, 
        multiple = FALSE
    )
})


output$dm_md_secondall <- renderUI({
    input$md_file
    
    md_data <- mergedata()
    nameValue = names(md_data)
    
    selectInput(
        "dm_md_secondall",
        "All Observarions in the Second Datasets Included",
        choices = as.list(c(FALSE,TRUE)), 
        selected = FALSE, 
        multiple = FALSE
    )
})


# Merge Dataset Preview ####

observeEvent(input$dm_md_preview, {
    if (input$dm_md_varname1 == "" || input$dm_md_varname2 == "") {
        return()
    }
    md_data <<- mergedata()
    key1 = input$dm_md_varname1
    key2 = input$dm_md_varname2
    g_md_countrows <<- 0
    lapply(data2[[key1]], function(i) {
        g_md_countrows <<- g_md_countrows + sum(i == md_data[[key2]])
    }) 
    if (g_md_countrows > nrow(data2)) {
        showNotification("Too many overlap data. Check key variable", type="warning")
        return()
    }

    arg1=FALSE
    arg2=FALSE
    if (input$dm_md_mainall=="TRUE") arg1=TRUE 
    if (input$dm_md_secondall=="TRUE") arg2=TRUE 

    g_md_previewdata <<- merge(data2, md_data, by.x = key1, by.y = key2,all.x=arg1,all.y=arg2)
})

# Merge Dataset Run ####

observeEvent(input$dm_md_run, {
    if (input$dm_md_varname1 == "" || input$dm_md_varname2 == "") {
        return()
    }
    md_data <<- mergedata()
#    md_data <- readr::read_delim(input$md_file$datapath, comment="#", col_names=input$md_header, delim=input$md_sep)
#    for (i in 1:length(md_data)) {
#        names(md_data)[i] <- paste0(names(md_data[i]), ".1")
#    }
    key1 = input$dm_md_varname1
    key2 = input$dm_md_varname2
    g_md_countrows <<- 0;
    lapply(data2[[key1]], function(i) {
        g_md_countrows <<- g_md_countrows + sum(i == md_data[[key2]])
    }) 
    if (g_md_countrows > nrow(data2)) {
        showNotification("Too many overlap data. Check key variable", type="warning")
        return()
    }
    arg1=FALSE
    arg2=FALSE
    if (input$dm_md_mainall=="TRUE") arg1=TRUE 
    if (input$dm_md_secondall=="TRUE") arg2=TRUE 

    data2 <<- merge(data2, md_data, by.x = key1, by.y = key2,all.x=arg1,all.y=arg2)

    md_attributes <- sapply(md_data, class)
    md_names <- names(md_data)
    md_whichnames <- which(md_names==key2)
    md_colattr <- attr(md_data, "spec")
    for (i in 1:length(md_names)) {
        if (i == md_whichnames) {
            next
        }
        g_attributes <<- c(g_attributes, md_attributes[[i]])
        g_colattr$cols[[md_names[i]]] <<- md_colattr$cols[[i]]
    }
})

# Merge Dataset Components ####

output$dm_md_selectvarname1 <- renderUI({
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
    
    nameValue = names(data2)
    
    selectInput(
        "dm_md_varname1",
        "Reference Variable for Main Dataset",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = FALSE
    )
})

output$dm_md_selectvarname2 <- renderUI({
    input$md_file
    
    md_data <- mergedata()
    nameValue = names(md_data)
    
    selectInput(
        "dm_md_varname2",
        "Reference Variable for Second Dataset",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = FALSE
    )
})