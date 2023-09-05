# Make Interval Variable Run ####

observeEvent(input$dm_mvi_run, {
  if (input$dm_mvi_newvariablename == "") {
    return()
  }
  
  newVariableName <- input$dm_mvi_newvariablename
  nameAccepted = which(names(data2)==newVariableName)
  if (length(nameAccepted) != 0) {
    showNotification("Type another variable name.", type="warning")
    return()
  }
  
  interval_method="intervals"
  if (input$dm_mvi_method=="Equal-count") interval_method="proportions"
  if (input$dm_mvi_method=="Natual breaks (K-means clustering)") interval_method="natural"
  new_interval=bin.var(data2[input$dm_mvi_varname][[1]],bins=input$dm_mvi_interval_count,method=interval_method)
  
#  tryCatch(expr = data2<<-mutate(data2, newVariable=!!rlang::parse_expr(input$dm_mvi_expression)),
  tryCatch(expr = data2<<-mutate(data2, newVariable=new_interval),
             error = function(e) {
             showNotification("Wrong Expression.", type="warning")
             g_error <<- TRUE
           }
  )
  
  if (g_error == TRUE) {
    g_error <<- FALSE
    return()
  }
#  print(new_interval)
#  data2<<-mutate(data2,newVariable=new_interval)
  names(data2)[length(data2)]<<-newVariableName
  
  if (class(data2[[length(data2)]]) == "integer") {
    g_attributes <<- c(g_attributes, "integer")
    g_resetattributes <<- c(g_resetattributes, "integer")
  } else {
    g_attributes <<- c(g_attributes, "numeric")
    g_resetattributes <<- c(g_resetattributes, "numeric")
  }
  
  g_colattr$cols[[newVariableName]] <<- col_double()
  g_resetcolattr$cols[[newVariableName]] <<- col_double()
  
  # run mutate -> attribute removed
  # if attribute == NULL editData is not working
  # make attribute
  attr(data2, "spec") <<- g_colattr
})

# Make Interval Variable Components ####

output$dm_mvi_newvariablename <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  
  textAreaInput("dm_mvi_newvariablename","New Variable Name",value=NULL)
})

output$dm_mvi_expression<-renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$dm_mvi_varname
  input$dm_mvi_method
  input$dm_mvi_interval_count
  
  expressionValue=NULL
  interval_method="intervals"
  if (input$dm_mvi_method=="Equal-count") interval_method="proportions"
  if (input$dm_mvi_method=="Natural breaks (K-means clustering)") interval_method="natural"
  expressionValue=paste0("bin.var","(",input$dm_mvi_varname,",","bins=",input$dm_mvi_interval_count,",","method=","'",interval_method,"'",")")
  textAreaInput("dm_mvi_expression","R code",value=expressionValue)
})

output$dm_mvi_check_tools <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  
  checkboxInput("dm_mvi_check_tools", "Selection for Making Interval Variable", value = FALSE)
})

output$dm_mvi_numberinterval <- renderUI({
#  if(is.null(input$dm_mvi_check_tools) || input$dm_mvi_check_tools == FALSE)
#      return()
  sliderInput(
    inputId = sprintf("dm_mvi_interval_count"),
    label = "Number of Intervlas",
    value = 3,
    min = 1,
    max = 20
  )
})

output$dm_mvi_selectvarname <- renderUI({
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
  
#  if(is.null(input$dm_mvi_check_tools) || input$dm_mvi_check_tools == FALSE)
#    return()
  
  div(
    selectInput(
      "dm_mvi_varname",
      "Variable",
      choices = as.list(c("", names(select_if(data2, is.numeric)))), 
      selected = NULL, 
      multiple = TRUE
    ),
    bsTooltip(
      "dm_mvi_varname", 
      "Numeric Only",
      "right", 
      options = list(container = "body")
    )
  )
})

output$dm_mvi_selectmethod <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  
#  if(is.null(input$dm_mvi_check_tools) || input$dm_mvi_check_tools == FALSE)
#    return()
  
  Methods=c("Equal-width","Equal-count","Natual breaks (K-means clustering)")
  selectInput("dm_mvi_method","Method", choices=Methods, selected="Equal-width", multiple = FALSE)
})