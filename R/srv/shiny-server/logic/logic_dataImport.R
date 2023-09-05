
# reactive maindata ####
# for attribute(data type) on csv type
maindata <- reactive({
    
    if(is.null(input$file) || is.null(input$file$datapath)) {return()}    
    
    l_datapath = strsplit(input$file$datapath, split = "[.]")[[1]]
    g_datatype <<- toupper(l_datapath[length(l_datapath)])
    
    
        
    if(g_datatype == "CSV" || g_datatype == "TXT") {
        if(is.null(input$sep) || is.null(input$header)) {return()}
        return(read.table(file=input$file$datapath, sep=input$sep, header=input$header, stringsAsFactors = TRUE, fileEncoding='EUC-Kr'))
        
    } else if(g_datatype == "XLSX" || g_datatype == "XLS") {
        if(is.null(g_sheetindex) || is.null(input$header)) {return()}
        return(readxl::read_excel(path=input$file$datapath, sheet = g_sheetindex, col_names = input$header))
        
    } else if(g_datatype == "RDATA" && !is.null(input$di_option_rdata_index)) {
        if(is.null(g_env) || is.null(input$di_option_rdata_index)) {return()}
        return(g_env[[input$di_option_rdata_index]])
    } else if(g_datatype == "RDS") {
        return(as.data.frame(readRDS(file = input$file$datapath)))
    } else if(g_datatype == "SAV") {
        return(haven::read_sav(file = input$file$datapath))
    }
    else {
        return()
    }
})

# File Option ####

output$di_option <- renderUI({
    input$file
    
    if(is.null(input$file) || is.null(input$file$datapath)) {
        return()
    }
    
    l_datapath = strsplit(input$file$datapath, split = "[.]")[[1]]
    g_datatype <<- toupper(l_datapath[length(l_datapath)])
    
    if(g_datatype == "CSV" || g_datatype == "TXT") {
        uiOutput("di_option_csv")
    } else if(g_datatype == "XLSX" || g_datatype == "XLS") {
        uiOutput("di_option_xlsx")
    } else if(g_datatype == "RDATA") {
        uiOutput("di_option_rdata")
    } else {
        return()
    }
})

output$di_option_header <-renderUI({
    checkboxInput(
        inputId = 'header',
        label = 'Header',
        value = TRUE
    )
})

output$di_option_run <- renderUI({
    actionButton("di_option_run", "File Open", icon = icon("play"))
})

output$di_option_csv <- renderUI({
    input$file
    input$resetData
    
    tagList(
        uiOutput("di_option_header"),
        fluidRow(
            column(
                width = 3, 
                selectInput(
                    "sep",
                    "Separator",
                    choices = c(
                        Comma = ',',
                        Semicolon = ';',
                        Tab = '\t',
                        Space = ' '
                    ),
                    multiple = FALSE
                )
                
            ),
            column(
                width = 3,
                uiOutput("di_option_run"),
                style = "text-align:left; padding:15px"
            )
        )
    )
})

output$di_option_xlsx <- renderUI({
    input$file
    input$resetData
    
    #l_datapath = strsplit(input$file$datapath[i], split = "[.]")[[1]]
    #g_datatype <<- l_datapath[length(l_datapath)]
    
    if(g_datatype == "XLSX" || g_datatype == "XLS") {
        l_sheetlist = readxl::excel_sheets(input$file$datapath)
        l_sheetnum = length(l_sheetlist)
        l_sheetchoice = 1:l_sheetnum
        names(l_sheetchoice) = l_sheetlist
        
        tagList(
            uiOutput("di_option_header"),
            fluidRow(
                column(
                    width = 3,
                    selectInput(
                        "di_option_xlsx_index",
                        "Sheet Index", 
                        choices = l_sheetchoice,
                        multiple = FALSE
                    )
                ),
                column(
                    width = 3,
                    uiOutput("di_option_run"),
                    style = "text-align:left; padding:15px"
                )
            )
        )
    } else {return()}
})

output$di_option_rdata <- renderUI({
    input$file
    input$resetData
    
    if(is.null(g_env)) {return()}
    
    if(g_datatype == "RDATA") {
        load(file = input$file$datapath, envir = g_env)
        l_envnames = names(g_env)
        l_envnames = l_envnames[unlist(lapply(l_envnames, function(x) {is.data.frame(g_env[[x]]) || is.matrix(g_env[[x]])}))]
        
        
        if(length(l_envnames) == 0) {return()}
        
        fluidRow(
            column(
                width = 3,
                selectInput(
                    "di_option_rdata_index",
                    "Choose Data names", 
                    choices = l_envnames,
                    multiple = FALSE
                )
            ),
            #column(
            #    width = 4,
            #    uiOutput("di_option_rdata_write")
            #),
            column(
                width = 3,
                uiOutput("di_option_run"),
                style = "text-align:left; padding:15px"
            )
        )
    } else {return()}
   
})

observeEvent(input$di_option_run, {
    g_openfile <<- TRUE
    g_resetresult <<- FALSE
})


#observeEvent(input$di_sheetindex_run, {
#    g_sheetindex_change <<- TRUE
#    g_resetresult <<- FALSE
#})



# editUI ####
# this function run at import data, change data, push reset button.
# Set all variables to initial values.
# g_attributes      - data attributes
# g_resetattributes - data attributes use to reset
# g_error           - Check the error when trycatch used.
# g_resetresult     - Reset the result when the user changes the file.

output$editUI=renderUI({
    input$file
    input$di_option_run
    input$resetData
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$file))
        return()
    
    data<-uiname<-result<-mylist<-textname<-list()
    count=length(input$file$datapath)
    
    
    
    #if(!g_openfile) {return()}
    if(is.null(maindata())) {return()}
    
    
    
    # reset global variables
    if (is.null(data2) || g_filepath != input$file$datapath || g_reset) {
        g_datatype <<- NULL
        g_sheetindex <<- 1
        #g_env <<- new.env()
        g_attributes <<- sapply(maindata(), class)
        g_resetattributes <<- sapply(maindata(), class)
        g_error <<- FALSE
        g_resetresult <<- FALSE
        g_sr_previewdata <<- NULL
        g_md_previewdata <<- NULL
        g_ba_an_interaction <<- NULL
        
        g_rg_lm_interaction <<- NULL
        g_rg_lm_comparisonmodel <<- NULL
        g_rg_glm_interaction <<- NULL
        g_rg_glm_comparisonmodel <<- NULL
        g_rg_logit_interaction <<- NULL
        g_rg_probit_interaction <<- NULL
        g_rg_loglinear_interaction <<- NULL
        
        g_rg_jglm_m_interaction <<- NULL
        g_rg_jglm_p_interaction <<- NULL
        g_rg_jglm_r_comparisonmodel_1 <<- NULL
        g_rg_jglm_r_comparisonmodel_2 <<- NULL
        
        g_re_lmm_m_interaction <<- NULL
        g_re_lmm_m_randinteraction <<- NULL
        g_re_lmm_r_comparisonmodel_1 <<- NULL
        g_re_lmm_r_comparisonmodel_2 <<- NULL
        
        g_re_glmm_m_interaction <<- NULL
        g_re_glmm_m_randinteraction <<- NULL
        g_re_glmm_r_comparisonmodel_1 <<- NULL
        g_re_glmm_r_comparisonmodel_2 <<- NULL
        
        g_re_hglm_m_interaction <<- NULL
        g_re_hglm_m_randinteraction <<- NULL
        g_re_hglm_p_interaction <<- NULL
        g_re_hglm_l_interaction <<- NULL
        g_re_hglm_r_comparisonmodel_1 <<- NULL
        g_re_hglm_r_comparisonmodel_2 <<- NULL
        g_re_hglm_r_comparisonmodel_3 <<- NULL
        
        g_re_dhglm_m_interaction <<- NULL
        g_re_dhglm_m_randinteraction <<- NULL
        g_re_dhglm_p_interaction <<- NULL
        g_re_dhglm_p_randinteraction <<- NULL
        g_re_dhglm_l_interaction <<- NULL
        g_re_dhglm_l_randinteraction <<- NULL
        g_re_dhglm_r_comparisonmodel_1 <<- NULL
        g_re_dhglm_r_comparisonmodel_2 <<- NULL
        g_re_dhglm_r_comparisonmodel_3 <<- NULL
        
        g_re_cs_m_interaction <<- NULL
        g_re_cs_p_interaction <<- NULL

        g_mr_mdhglm_m_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_mdhglm_m_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_mdhglm_p_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_mdhglm_p_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_mdhglm_l_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_mdhglm_l_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_mdhglm_r_comparisonmodel_1 <<- NULL
        g_mr_mdhglm_r_comparisonmodel_2 <<- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_mdhglm_r_comparisonmodel_3 <<- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_mdhglm_r_comparisonmodelcount <<- 0
        
        g_mr_factor_m_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_factor_m_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_factor_p_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_factor_p_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_factor_l_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_factor_l_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)      
        
        g_mr_shared_m_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_shared_m_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_shared_p_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_shared_p_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_shared_l_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_shared_l_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)      
        
        g_mr_joint_m_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_joint_m_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_joint_p_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_joint_p_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_joint_l_interaction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        g_mr_joint_l_randinteraction <<- list(NULL, NULL, NULL, NULL, NULL, NULL)
        
        g_mr_joint_s1_interaction <<- NULL
        g_mr_joint_s2_interaction <<- NULL
        
        g_sa_cox_interaction <<- NULL
        g_sa_fm_interaction <<- NULL
        g_sa_fm_randinteraction <<- NULL
        g_sa_fm_comparisonmodel <<- NULL
        g_sa_cr_s1_interaction <<- NULL
        g_sa_cr_s2_interaction <<- NULL
    }
    
    
    #if datapath is not empty, load data
    if(count>0) {
        for(i in  1:count){
            # save data2 & change dataframe attribute 
            if (is.null(data2) || g_filepath != input$file$datapath[i] || g_openfile || g_reset ) {
                #g_env <- new.env()
                g_openfile <<- FALSE
                g_reset <<- FALSE
                
                l_datapath = strsplit(input$file$datapath[i], split = "[.]")[[1]]
                g_datatype <<- toupper(l_datapath[length(l_datapath)])
                
                if(g_datatype == "CSV" || g_datatype == "TXT") {
                    tempdata <<- readr::read_delim(input$file$datapath[i], delim = input$sep, comment="#", col_names = input$header, locale=locale(encoding='EUC-Kr'))
                    
                    g_colattr <<- attr(tempdata, "spec")
                    g_resetcolattr <<- attr(tempdata, "spec")
                    g_filepath <<-input$file$datapath
                    g_attributes <<- sapply(maindata(), class) #temporarly
                    lapply(1:length(tempdata), function(k) {
                        if (class(tempdata[[k]]) == "character") {
                            tempdata[[k]] <<- factor(tempdata[[k]])
                        }
                        if (g_attributes[[k]] == "integer") {
                            tempdata[[k]] <<- as.integer(tempdata[[k]])
                        }
                    })
                    data[[i]]<-tempdata
                } else if(g_datatype == "XLSX" || g_datatype == "XLS") {
                    g_sheetindex <<- as.integer(input$di_option_xlsx_index)
                    tempdata <- readxl::read_excel(path=input$file$datapath, sheet = g_sheetindex, col_names = input$header)
                    
                    g_colattr <<- attr(tempdata, "spec")
                    g_resetcolattr <<- attr(tempdata, "spec")
                    g_filepath <<-input$file$datapath[i]
                    lapply(1:length(tempdata), function(k) {
                        if (class(tempdata[[k]]) == "character") {
                            tempdata[[k]] <<- factor(tempdata[[k]])
                        }
                    })
                    data[[i]]<-tempdata
                } else if(g_datatype == "RDATA") {
                    tempdata <- g_env[[input$di_option_rdata_index]]
                    
                    g_colattr <<- attr(tempdata, "spec")
                    g_resetcolattr <<- attr(tempdata, "spec")
                    g_filepath <<-input$file$datapath[i]
                    lapply(1:length(tempdata), function(k) {
                        if (class(tempdata[[k]]) == "character") {
                            tempdata[[k]] <<- factor(tempdata[[k]])
                        }
                    })
                    data[[i]]<-tempdata
                } else if(g_datatype == "RDS") {
                    tempdata <- as.data.frame(readRDS(file = input$file$datapath[i]))
                    
                    g_colattr <<- attr(tempdata, "spec")
                    g_resetcolattr <<- attr(tempdata, "spec")
                    g_filepath <<-input$file$datapath[i]
                    lapply(1:length(tempdata), function(k) {
                        if (class(tempdata[[k]]) == "character") {
                            tempdata[[k]] <<- factor(tempdata[[k]])
                        }
                    })
                    data[[i]]<-tempdata
                } else if(g_datatype == "SAV") {
                    tempdata <- haven::read_sav(file = input$file$datapath[i])
                    
                    g_colattr <<- attr(tempdata, "spec")
                    g_resetcolattr <<- attr(tempdata, "spec")
                    g_filepath <<-input$file$datapath[i]
                    lapply(1:length(tempdata), function(k) {
                        if (class(tempdata[[k]]) == "character") {
                            tempdata[[k]] <<- factor(tempdata[[k]])
                        }
                    })
                    data[[i]]<-tempdata
                }
            } else {
                attr(data2, "spec") <<- g_colattr
                data[[i]]<-data2
            }
            uiname[[i]]<-paste0("table",i)
            title=NULL
            mylist[[3*i-2]]<-h2(title)
            mylist[[3*i-1]]<-editableDTUI(uiname[[i]])
            textname[[i]]=paste0("text",i)
            mylist[[3*i]]<-uiOutput(textname[[i]])
            
            # editUI (Delete, Add New, EditData)
            local({
                j<-i
                result[[j]]=callModule(editableDT,uiname[[j]],data=reactive(data[[j]]))
                
                output[[textname[[j]]]]=renderUI({
                    data2<<-result[[j]]()
                    data2<<-type_convert(data2, attr(data2, "spec"))
                    lapply(1:length(data2), function(k) {
                        if (class(data2[[k]]) == "character") {
                            data2[[k]] <<- factor(data2[[k]])
                        }
                        if (g_attributes[[k]] == "integer") {
                            data2[[k]] <<- as.integer(data2[[k]])
                        }
                    })
                    return()
                })
            })
        }
        do.call(tagList, mylist)
    }
})

# Download & Reset ####

# Download Button Action
output$downloadData <- downloadHandler(
    filename = function() {
        paste("data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
    },
    content = function(downloadfile) {
        write.csv(data2, downloadfile, row.names = FALSE)
    }
)

output$download_rdata <- downloadHandler(
    filename = function() {
        paste("data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData", sep = "")
    },
    content = function(downloadfile) {
        save_env <- new.env()
        save_env$Albatross = data2
        save(Albatross, file = downloadfile, envir = save_env)
    }
)

# Reset Button Action
observeEvent(input$resetData, {
    g_reset <<- TRUE
})

# Download Plot
downloadPlot <- function(plot) {
    downloadHandler(
        filename = function() {
            paste("plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, plot = plot, width = 8, type = "cairo")
        }
    )
}

# Download K-M Plot
downloadSurvPlot <- function(plot) {
    downloadHandler(
        filename = function() {
            paste("plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, print(plot), width = 8, type = "cairo")
        }
    )
}
