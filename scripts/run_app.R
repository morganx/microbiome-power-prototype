args <- commandArgs(trailingOnly=TRUE)
name <- if(length(args)) args[1] else "case-control"
if(!name %in% c("case-control","continuous")) stop("Choose case-control or continuous")
if(!file.exists("DESCRIPTION")) stop("Run from the repository root")
e <- new.env(parent=globalenv())
app <- source(file.path("apps",name,"app.R"),local=e)$value
port <- if(length(args)>1) as.integer(args[2]) else if(name=="case-control") 3850L else 3851L
shiny::runApp(app,host="127.0.0.1",port=port,launch.browser=FALSE)
