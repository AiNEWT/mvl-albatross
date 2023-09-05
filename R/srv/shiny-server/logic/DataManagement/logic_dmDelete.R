# Delete variable ####

observeEvent(input$dm_del_run, {
  if (is.null(input$dm_del_varname1)) {
    return()
  }
  length1=length(input$dm_del_varname1)
  for (i in 1:length1) {
       loc= which(names(data2)==input$dm_del_varname1[i])
       data2 <<- data2[-loc]
  }
})

# Select Rows Components ####

output$dm_del_selectvarname1 <- renderUI({
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
    "dm_del_varname1",
    "Variable",
    choices = as.list(c("", nameValue)), 
    selected = NULL, 
    multiple = TRUE
  )
})


