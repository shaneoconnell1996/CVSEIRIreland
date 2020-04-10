# Seir / shiny app
# Last updated,
#
# Thu Apr  2 15:06:17 IST 2020 Brian (B.OSULLIVAN19@nuigalway.ie)
# Initial test version.
#
# Fri Apr  5 14:19:46 IST 2020 Brian (B.OSULLIVAN19@nuigalway.ie)
# Update from Nick adding get_I0 function (in associated support.R) to 
# finds the value of I0 such that an initial SEIR vector of
# (S=1, E=0, I=I0, R=0) at t=Start evolves to hit the right
# number of cases at t=hosts_guess_t.
# Also general updates to UI to tidy it up and make inputs resetable.
# Change legends and add deaths to plot.
#
# Mon Apr  6 21:58:11 IST 2020 Brian (B.OSULLIVAN19@nuigalway.ie)
# Remoce line "hse = hse[-14,]" (that was left in by mistake).

# Tues Apr 7 Nick (nick.tosh@nuigalway.ie)
# Tweaked plot axis labels, added plot titles, changed range of times
# plotted, changed colour of E plot, reduced weight of dashed lines,
# tweaked mins and step sizes for input parameters, and updated HSE data.

if (!requireNamespace("shiny", quietly = TRUE))
  BiocManager::install("shiny")

if (!requireNamespace("deSolve", quietly = TRUE))
  BiocManager::install("deSolve")

if (!requireNamespace("reshape2", quietly = TRUE))
  BiocManager::install("reshape2")

if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")

if (!requireNamespace("gridExtra", quietly = TRUE))
  BiocManager::install("gridExtra")

source('Seir_model.R')
source('support_modified.R')

library(deSolve)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(shiny)
# user interface.
ui <- fluidPage(
  
  # App title ----
  #  titlePanel("Seir Model (COVID-19 Ireland), last updated 5/4/2020"),
  titlePanel(title=div(img(src="seir_logo.png", width="100%"))),
  
  HTML(paste0('<div align="left"style="padding-top: 0%; padding-bottom: 0%; font-size: xx-small;">last updated 3/4/2020</div>')),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      uiOutput('resetable_inputs'),
      HTML('<div align="center" style="padding-top: 2%;">'),
      actionLink("reset_input", "reset inputs"),
      HTML('</div>'),
      width = 4
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(   
      tabsetPanel(type = "tabs",
                  tabPanel("Epidemiological Model", plotly::plotlyOutput(outputId = "seirPlot1")),
                  tabPanel("HSE Data", plotly::plotlyOutput(outputId = "seirPlot2"))
                  # tabPanel("Table", tableOutput("table"))
                  
      )
    )
  )
)

# Server
server <- function(input, output) {
  
  # We setup the UI here so we can recall it if the user wants to reset to initial values.
  data <- reactive({
    if(is.null(input$hosts_count_guess))
    {
      # print("got a null in reactivePlot...")
      return(NULL)
    }
    
    
    # Time measures: March 1st is t=1. March 2nd is t=2, etc.
    df <- create_df(-18,R0_pre=input$R0_pre,R0_post=input$R0_post,R0_post2=input$R0_post2,mu=input$mu,CV=input$CV,ICU_rate=input$ICU_rate,
                    hosts_count_guess=input$hosts_count_guess)[[1]]
    #p1_t_start = -10
  })
  data_hse <- reactive({
    df <- create_df(-18,R0_pre=input$R0_pre,R0_post=input$R0_post,R0_post2=input$R0_post2,mu=input$mu,CV=input$CV,ICU_rate=input$ICU_rate,
                    hosts_count_guess=input$hosts_count_guess)[[3]]
  })
  data2 <- reactive({
    if(is.null(input$hosts_count_guess))
    {
      # print("got a null in reactivePlot...")
      return(NULL)
    }
    df <- create_df(-18,R0_pre=input$R0_pre,R0_post=input$R0_post,R0_post2=input$R0_post2,mu=input$mu,CV=input$CV,ICU_rate=input$ICU_rate,
                    hosts_count_guess=input$hosts_count_guess)[[2]]
    
  })
  output$resetable_inputs <- renderUI({
    input$reset_input
    
    
    div(
      tags$head(
        tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle;padding: 0%;} 
                 .inline .form-group{display: table-row;}
                 .inline .form-control{width: 100%;}")
      ),
      
      HTML('<div align="left">'),
      
      HTML(paste0('<div style="align: center; padding-top: 0%; padding-bottom: 2%;">Epidemiological parameters:</div>')),
      
      HTML('<div style="display:inline-block; padding-right: 8px; padding-bottom: 5%; bottom: 15px; font-size: 11px">'),
      # Number of hosts (exposed + infectious) at hosts_guess_t
      tags$div(class = "inline", numericInput("hosts_count_guess", label=HTML("Estimated number of affected individuals (E+I) on 1/3/2020"), value = 650, min = 50, max = 10000, step = 50, width = "100px")),
      HTML('</div>'),
      
      
      HTML('<div style="display:inline-block; width: 32%; bottom: 15px; font-size: 11px">'),
      numericInput("R0_pre", HTML("R<sub>0</sub> pre 1<sup>st</sup> intervention"), value = 3, min = 0, max = 1000, step = 0.05, width = "80%"),
      HTML('</div>'),
      
      HTML('<div style="display:inline-block; width: 32%; bottom: 15px; font-size: 11px">'),
      numericInput("R0_post", HTML("R<sub>0</sub> post 1<sup>st</sup> intervention"), value = 1.05, min = 0, max = 1000, step = 0.05, width = "80%"),
      HTML('</div>'),
      
      
      HTML('<div style="display:inline-block; width: 32%; bottom: 15px; font-size: 11px">'),
      numericInput("R0_post2", HTML("R<sub>0</sub> post 2<sup>nd</sup> intervention"), value = 0.95, min = 0, max = 1000, step = 0.05, width = "80%"),
      HTML('</div>'),
      
      HTML('</div>'),
      
      
      
      
      
      
      # Advanced options.
      
      # Gamma-distributed delay to ICU
      
      HTML('<div style="border-radius: 5px; padding: 5px; background: #ebebeb; border: 1px solid black;">'),
      
      HTML(paste0('<div style="align: center; padding-top: 0; padding-bottom: 2%; font-size: 14px">Clinical parameters:</div>')),
      
      HTML('<div align="left">'),
      
      HTML('<div style="display:inline-block; padding-right: 16px; bottom: 15px; font-size: 11px">'),
      # Mean
      numericInput("mu", "mean:", value = 14, min = 0.5, max = 30, step = 0.5, width = "65px"),
      HTML('</div>'),
      
      HTML('<div style="display:inline-block; padding-right: 16px; bottom: 15px; font-size: 11px">'),
      #coeffcient of variation
      numericInput("CV", "coeff.var.:", value = 0.4, min = 0.1, max = 1, step = 0.1, width = "65px"),
      HTML('</div>'),
      
      HTML('<div style="display:inline-block; padding-right: 0px; bottom: 15px; font-size: 11px">'),
      # Estimated ICU admissions
      # highly uncertain; user may wish to experiment
      numericInput("ICU_rate", "ICU rate", value = 0.013, min = 0, max = 1, step = 0.001, width = "80px"),
      HTML('</div>'),
      
      
      
      HTML('</div>'),
      
      HTML('</div>')
      
    )
  })
  
  output$seirPlot1 <- plotly::renderPlotly({
    p1 <- plotly::ggplotly(
      
      ggplot(data(),aes(time,value)) + geom_line(aes(colour = series)) +
        xlab('time (days)') +
        ylab('number of people') +
        geom_vline(xintercept=13,color="red",linetype="dashed",size=0.4) +
        geom_vline(xintercept=27,color="red",linetype="dashed",size=0.4) +
        geom_vline(xintercept=as.numeric(difftime(Sys.Date(),as.Date("2020-3-1"))),color="black",linetype="dashed", size=0.4) +
        scale_colour_manual(name = "Breakdown", values=c("orange", "green","blue")) +
        ggtitle('Epidemiological model') +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=12),
              legend.justification=c(1,1),legend.position=c(1,1),legend.title=element_blank()
        ) )
    
    
    p1
    # OK, debounce or throttle don't work the way I want them to here.
    # Shelve it for now.
    # TODO!!!! find out how to take only the last event when a number a fired in quick succession.   
    #db_reactivePlot()
    # reactivePlot()
    
    
  })
  output$seirPlot2 <- plotly::renderPlotly({
    p2 <- plotly::ggplotly(
      ggplot(data = data2(),aes(time,value)) +
        geom_line(aes(colour="C.ICU")) +
        #      geom_point(data=hse,aes(day,deaths, colour = "Deaths")) +
        geom_point(data=data_hse(),aes(day,icu, colour = "ICU")) +
        #      scale_colour_manual(name = "Breakdown", values=c("blue", "black","red")) +
        scale_colour_manual(name = "Breakdown", values=c("blue","red")) + geom_line(color="black") + 
        geom_vline(xintercept=13,color="red",linetype="dashed", size=0.4) +
        geom_vline(xintercept=27,color="red",linetype="dashed", size=0.4) +
        geom_vline(xintercept=as.numeric(difftime(Sys.Date(),as.Date("2020-3-1"))),color="black",linetype="dashed", size=0.4) +
        xlab('time (days)') +
        ggtitle('HSE data') +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              axis.title.y = element_blank(),
              legend.text=element_text(size=12),
              legend.justification=c(1,1),legend.position=c(1,1),legend.title=element_blank()
        ) )
    p2
  })
}
# Run the application
shinyApp(ui = ui, server = server)
