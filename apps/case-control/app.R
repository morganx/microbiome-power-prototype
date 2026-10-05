# Launch from repository root: Rscript scripts/run_app.R case-control
library(shiny)
source("R/training_data.R",local=TRUE)
source("R/case_control_power.R",local=TRUE)
source("R/case_control_grant_text.R",local=TRUE)
training <- load_applet_training()
ui <- fluidPage(
  titlePanel("Case–control power • cross-sectional"),
  sidebarLayout(sidebarPanel(
    selectInput("dataset","Training feature distributions",unique(training$key),selected="Original: HMP1-2 taxonomy stool"),
    numericInput("cases","Number of independent cases",100,min=2,step=1),
    numericInput("controls","Number of independent controls",100,min=2,step=1),
    sliderInput("target","Target power",min=.5,max=.99,value=.8,step=.01),
    numericInput("family_alpha","Family-wise significance level",.05,min=.0001,max=.2),
    numericInput("features","Features tested",260,min=1,step=1),
    numericInput("tests_per_feature","Tests per feature in the correction family",1,min=1,step=1),
    helpText("One independent observation per person. Two real groups are required."),
    selectInput("grant_quartile","Feature SD for grant summary",c("75th percentile"="75%","Median"="50%","25th percentile"="25%"),selected="75%"),
    helpText("Grant examples use a fixed 1% relative abundance in controls."),
    downloadButton("download","Download estimates"),
    downloadButton("downloadSummary","Download grant text")
  ),mainPanel(
    uiOutput("design"),
    tableOutput("estimates"),
    plotOutput("curves"),
    h4("Grant-ready summary"),
    uiOutput("grantSummary"),
    p("The detectable effect is a difference in mean arcsine-square-root transformed feature abundance between cases and controls. Curves use the 25th, 50th and 75th percentiles of training SDs."),
    p("Assumes independent groups and a common within-group SD. Training SDs exclude zeros and are planning proxies, not measured model residual SDs. This is an unadjusted two-sample t-test calculation; it does not model covariate adjustment or repeated measurements."),
    p("Multiplicity correction uses the tested feature count and explicit tests per feature. Covariates are not automatically counted as hypotheses. Sequencing depth and read-count translations are not modeled.")
  ))
)
server <- function(input,output,session) {
  selected <- reactive({req(input$dataset);training[training$key==input$dataset,,drop=FALSE]})
  observeEvent(input$dataset,{updateNumericInput(session,"features",value=selected()$total_features)},ignoreInit=FALSE)
  settings <- reactive({
    whole <- function(x) length(x)==1 && is.finite(x) && x==floor(x)
    validate(need(whole(input$cases) && input$cases>=2,"Enter at least two real cases."),
      need(whole(input$controls) && input$controls>=2,"Enter at least two real controls. No control group is imputed."),
      need(whole(input$features) && input$features>=1,"Feature count must be a positive integer."),
      need(whole(input$tests_per_feature) && input$tests_per_feature>=1,"Tests per feature must be a positive integer."),
      need(is.finite(input$family_alpha) && input$family_alpha>0 && input$family_alpha<1,"Significance must be between zero and one."))
    alpha <- input$family_alpha/input$features/input$tests_per_feature
    validate(need(is.finite(input$target) && input$target>alpha && input$target<1,"Target power must exceed per-test alpha and be below one."))
    list(cases=input$cases,controls=input$controls,alpha=alpha)
  })
  estimates <- reactive({
    t<-selected();s<-settings();sd<-as.numeric(t[1,c("sd_25","sd_50","sd_75")])
    data.frame(SD_quartile=c("25%","50%","75%"),Training_SD=sd,
      Detectable_difference=do.call(case_control_detectable,c(list(feature_sd=sd,target=input$target),s)))
  })
  grant_text <- reactive({
    req(input$grant_quartile)
    validate(need(input$grant_quartile %in% c("25%","50%","75%"),"Choose a feature SD percentile."))
    make_case_control_grant_text(selected(),settings(),input$target,input$features,
      input$tests_per_feature,input$family_alpha,input$grant_quartile)
  })
  output$grantSummary <- renderUI(tags$textarea(readonly="readonly",`aria-label`="Grant summary",
    style="width:100%;min-height:300px;font-size:16px;line-height:1.5;padding:12px;",grant_text()))
  output$downloadSummary <- downloadHandler(filename=function()"case_control_grant_summary.txt",
    content=function(file)writeLines(grant_text(),file,useBytes=TRUE))
  output$design <- renderUI({t<-selected();s<-settings();tagList(
    h3(t$dataset),p(t$definition),p(t$class_features," features retained; ",t$excluded," excluded. Training samples: ",t$samples,"."),
    p(s$cases," cases and ",s$controls," controls. Per-test alpha: ",format(s$alpha,scientific=TRUE),"."))})
  output$estimates <- renderTable({
    d<-estimates();d$From_1_percent_to<-case_control_abundance_example(d$Detectable_difference)
    names(d)<-c("Feature SD percentile","Training SD","Detectable transformed difference","Illustration: 1% becomes (%)")
    d
  },digits=6,na="Outside range")
  output$curves <- renderPlot({
    e<-estimates();s<-settings();x<-seq(0,max(e$Detectable_difference)*1.3,length.out=160)
    y<-sapply(e$Training_SD,function(sd)do.call(case_control_power,c(list(delta=x,feature_sd=sd),s)))
    matplot(x,y,type="l",lty=c(2,1,3),lwd=2,col=c("#6689A6","#176B83","#A65326"),
      xlab="Difference in transformed means",ylab="Power",ylim=c(0,1))
    abline(h=input$target,lty=2,col="gray50")
    legend("bottomright",legend=c("25% SD","50% SD","75% SD"),lty=c(2,1,3),col=c("#6689A6","#176B83","#A65326"),bty="n")
  })
  output$download <- downloadHandler(filename=function()"case_control_power.csv",content=function(file){
    s<-settings();d<-estimates();d$reference_abundance_percent<-1;d$illustrative_case_abundance_percent<-case_control_abundance_example(d$Detectable_difference);d$dataset<-input$dataset;d$cases<-s$cases;d$controls<-s$controls
    d$family_alpha<-input$family_alpha;d$features_tested<-input$features;d$tests_per_feature<-input$tests_per_feature
    d$per_test_alpha<-s$alpha;d$target_power<-input$target;write.csv(d,file,row.names=FALSE)
  })
}
shinyApp(ui,server)
