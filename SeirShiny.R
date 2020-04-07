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
source('support.R')

library(deSolve)
library(reshape2)
library(ggplot2)
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
      # Seir model output.
      plotOutput(outputId = "seirPlot")
    )
  )
)

# Server
server <- function(input, output) {
  
  # We setup the UI here so we can recall it if the user wants to reset to initial values.
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
      tags$div(class = "inline", numericInput("hosts_count_guess", label=HTML("Estimated number of affected individuals on 1/3/2020"), value = 650, min = 0, max = 1000, step = 10, width = "100px")),
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
      numericInput("mu", "mean:", value = 14, min = 0, max = 1000, step = 0.1, width = "65px"),
      HTML('</div>'),
      
      HTML('<div style="display:inline-block; padding-right: 16px; bottom: 15px; font-size: 11px">'),
      #coeffcient of variation
      numericInput("CV", "coeff.var.:", value = 0.4, min = 0, max = 10, step = 0.1, width = "65px"),
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
  
  reactivePlot <- reactive({
    # Make sure the user interface is there before trying to plot.
    # Only need to check the first variable, if thats there the rest should be also.
    if(is.null(input$hosts_count_guess))
    {
      # print("got a null in reactivePlot...")
      return(NULL)
    }
    
    
    # Time measures: March 1st is t=1. March 2nd is t=2, etc.
    Start = -18
    hosts_guess_t = 1     # user will supply guess re number of hosts at this time
    
    intervention_t = 13 # 1st intervention
    intervention_t2 = 27 # 2nd intervention
    End = 60
    dt = 1/20
    today = as.numeric(difftime(Sys.Date(),as.Date("2020-3-1")))
    
    
    
    # Other model parameters (empirical)
    
    # User may also want to experiment with these values, but
    # app should make clear that the default settings are well grounded
    # in the scientific literature
    infectious_period = 2.9 # from Goh epidemic simulator
    latent_period = 5.2     # from epidemic simulator (will need better source)
    
    #Seir model parameters based on above assumptions
    alpha = 1 / latent_period
    gamma = 1 / infectious_period
    beta1 = gamma * input$R0_pre
    beta2 = gamma * input$R0_post
    beta3 = gamma * input$R0_post2
    parameter_list = c (alpha, beta1, beta2, gamma, beta3, intervention_t, intervention_t2)
    
    # Compute initial SEIR state
    N = 4921500    # population of Republic (Wikipedia)
    # Number of hosts (exposed + infectious) at hosts_guess_t / Population
    hosts_guess = input$hosts_count_guess / N
    
    # compute I0 such that SEIR state (S=1, E=0, I=I0, R=0) at t=Start
    # will evolve, under SEIR dynamics, to a state with E+I = hosts_guess
    # at t = hosts_guess_t
    I0 = get_I0(hosts_guess, hosts_guess_t, Start, dt, seir_model, parameter_list)
    initial_values = c(S = 1, E = 0, I = I0, R = 0)
    
    # Solve.
    timepoints = seq (Start, End, by=dt)
    output = lsoda (initial_values, timepoints, seir_model, parameter_list, rtol=1e-7, atol=1e-7)
    
    # Estimate ICU admissions based on user-inputted ICU_rate, k, mu and CV.
    # filter width: number of mus to scan back. User doesn't need to play with this.
    k = 4
    
    # Compute ICU admissions per timestep, and hence also cumulative admissions.
    # Admissions rate now will be proportional to a weighted average of E values
    # over the past k*mu days.
    filter_weights <- input$ICU_rate*N*dt*dgamma(seq(0, k*input$mu, dt), shape=1/input$CV^2, rate=1/(input$mu * input$CV^2))
    admissions <- zp_filter(alpha*dt*output[,'E'], filter_weights)
    cumulative_admissions <- cumsum(admissions)
    
    # HSE data from date of first intervention + 3.
    hse = read.csv('hse_data.csv')
    
    # Number datapoints in the csv file.
    icupoints = nrow(hse)
    
    # Get per day increases on ICU admissions.
    hse = cbind(day=((intervention_t+3):(intervention_t+3+icupoints-1)),hse)
    daily_icu = c()
    daily_icu[2:icupoints] = hse$icu[2:icupoints] - hse$icu[1:(icupoints-1)]
    daily_icu[1] = hse$icu[1]
    hse = cbind(hse,daily_icu)
    
    # Get per day deaths increases.
    daily_deaths = c()
    daily_deaths[2:icupoints] = hse$deaths[2:icupoints] - hse$deaths[1:(icupoints-1)]
    daily_deaths[1] = hse$deaths[1]
    hse = cbind(hse,daily_deaths)
    
    
    # Plots
    disp = c(3,4)
    df = data.frame(time=output[,1],output[,disp]*N, ICU=cumulative_admissions)
    df = melt(df,id.vars = 'time', variable.name = 'series')
    p1 = ggplot(df, aes(time,value)) + geom_line(aes(colour = series)) +
      geom_vline(xintercept=intervention_t,color="red",linetype="dashed") +
      geom_vline(xintercept=intervention_t2,color="red",linetype="dashed") +
      geom_vline(xintercept=today,color="black",linetype="dashed") +
      scale_colour_manual(name = "Breakdown", values=c("red", "green","blue")) +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            legend.text=element_text(size=12),
            legend.justification=c(1,1),legend.position=c(1,1),legend.title=element_blank()
      ) 
    
    
    
    df = data.frame(time=output[,1], ICU= cumulative_admissions)
    df = melt(df,id.vars = 'time', variable.name = 'series')
    
    p2 = ggplot()+
      geom_line(data=df,aes(y=value,x=time,colour="C.ICU")) +
      geom_point(data=hse,aes(day,deaths, colour = "Deaths")) +
      geom_point(data=hse,aes(day,icu, colour = "ICU")) +
      scale_colour_manual(name = "Breakdown", values=c("blue", "black","red")) +
      geom_line(color="black")+geom_vline(xintercept=intervention_t,color="red",linetype="dashed") +
      geom_vline(xintercept=intervention_t2,color="red",linetype="dashed") +
      geom_vline(xintercept=today,color="black",linetype="dashed") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            legend.text=element_text(size=12),
            legend.justification=c(1,1),legend.position=c(1,1),legend.title=element_blank()
      ) 
    
    grid.arrange(p1,p2,nrow=1) 
    
  })
  
  #    db_reactivePlot <- debounce(reactivePlot, 200)
  
  output$seirPlot <- renderPlot({
    # OK, debounce or throttle don't work the way I want them to here.
    # Shelve it for now.
    # TODO!!!! find out how to take only the last event when a number a fired in quick succession.   
    #db_reactivePlot()
    reactivePlot()
    
  })
  
}
# Run the application
shinyApp(ui = ui, server = server)
