# Script for Shiny app

library("shiny")

# Define UI

ui <- fluidPage(

  # App title ----
  titlePanel("Score and visualize a variant's EVE scores"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Select a file ----
      fileInput("eveData", "EVE Data",
                multiple = FALSE,
                accept = ".vcf"),

      tags$p("Upload a vcf file from https://evemodel.org/ corresponding to the
             gene you would like to analyze.  The vcf file for the gene of
             interest can be found at https://evemodel.org/download/protein and
             searching for the gene."),

      # Horizontal line ----
      tags$hr(),

      # Input: Select a file ----
      fileInput("variant1", "Variant 1 data",
                multiple = FALSE,
                accept = ".csv"),

      tags$p("Upload a csv file containing single nucleotide variants (SNVs) for
             your chosen gene."),

      # add checkbox about if in protein form

      # Horizontal line ----
      tags$hr(),

      # Input: Select a file ----
      fileInput("variant2", "Variant 2 data (optional)",
                multiple = FALSE,
                accept = ".csv"),

      tags$p("Optionally, upload a second csv file containing single nucleotide
      variants (SNVs) for your chosen gene."),

      # add checkbox about if in protein form
      # Horizontal line ----
      tags$hr(),

      # TODO: Change so its a slider to select gene range
      # Input: Slider for the number of observations to generate ----
      sliderInput("n",
                  "Number of observations:",
                  value = 500,
                  min = 1,
                  max = 1000)

      ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ variant 1 plot, variant 2 plot, and two variants
      # overlapped ----
      tabsetPanel(type = "tabs",
                  tabPanel("Variant 1", plotOutput("plot")),
                  tabPanel("Variant 2", verbatimTextOutput("summary")),
                  tabPanel("Overlapped", tableOutput("contents"))
      )
    )
  )
)


# TODO: Change the server logic
# Define server logic
server <- function(input, output) {

  # processedVar1Data <- reactive({
  #   req(input$variant1)
  #   variant1Data <- processVariantData(filePath = readRDS(input$variant1$datapath))
  # })
  # req(input$variant1)

  # variant1Data <- processVariantData(filePath = readRDS(input$variant1$datapath))
  # output$contents <- renderTable({
  #   req(input$variant1)
  #   variant1Processed <- processVariantData(filePath = input$variant1$datapath)
  # })

  output$plot <- renderPlot({
    req(input$eveData)
    req(input$variant1)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath)
    scoredVariant <- getEveScores(eveDataProcessed, variant1Processed, TRUE)
    visualizeVariant(scoredVariant, "NRXN1")
  })

  # output$contents <- renderPlot({
  #   visualizeVariant(eveData, "NRXN1", FALSE)
  # })

  # output$contents <- renderTable({
  #
  #   # input$file1 will be NULL initially. After the user selects
  #   # and uploads a file, head of that data file by default,
  #   # or all rows if selected, will be shown.
  #
  #   req(input$variant1)
  #
  #   # when reading semicolon separated files,
  #   # having a comma separator causes `read.csv` to error
  #   tryCatch(
  #     {
  #       df <- read.csv(input$variant1$datapath)
  #     },
  #     error = function(e) {
  #       # return a safeError if a parsing error occurs
  #       stop(safeError(e))
  #     }
  #   )
  #
  #   return(df)
  #
  # })
}


# Create Shiny app
shiny::shinyApp(ui, server)

# [END]
