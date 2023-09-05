# Rename Preview ####

#observeEvent(input$dm_rn_preview, {
#  if (is.null(input$dm_rn_expression)) {
#    return()
#  }
#  if (is.null(input$dm_rn_varname1)) {
#    return()
#  }  
#  names(data2)[names(data2)==input$dm_rn_varname1] <-input$dm_rn_expression 
#})

# Rename Run ####

observeEvent(input$dm_rn_run, {
  if (is.null(input$dm_rn_expression)) {
    return()
  }
  if (is.null(input$dm_rn_varname1)) {
    return()
  }  
  names(data2)[names(data2)==input$dm_rn_varname1]<<-input$dm_rn_expression 
})

# Select Rows Components ####

output$dm_rn_selectvarname1 <- renderUI({
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
    "dm_rn_varname1",
    "Variable",
    choices = as.list(c("", nameValue)), 
    selected = NULL, 
    multiple = FALSE
  )
})


output$dm_rn_expression<-renderUI({
  input$dm_rn_varname1
  textAreaInput("dm_rn_expression","New Variable name",value=NULL)
})


