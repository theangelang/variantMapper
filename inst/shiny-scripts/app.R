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

      textInput("geneName", "Gene name", placeholder = "Name of gene"),

      tags$p("Enter in the name of the gene of interest."),

      # Horizontal line ----
      tags$hr(),

      # Input: Select a file ----
      fileInput("variant1", "Variant 1 data",
                multiple = FALSE,
                accept = ".csv"),

      tags$p("Upload a csv file containing single nucleotide variants (SNVs) for
             your chosen gene."),

      # Input: Select protein or genomic form ----
      radioButtons("form1", "Form of variant 1' data",
                   choices = c("Protein" = TRUE,
                               "Genomic" = FALSE),
                   selected = character(0)),

      # TODO: Add blurb about what it means to be in protein and genomic form

      # Horizontal line ----
      tags$hr(),

      # Input: Select a file ----
      fileInput("variant2", "Variant 2 data (optional)",
                multiple = FALSE,
                accept = ".csv"),

      tags$p("Optionally, upload a second csv file containing single nucleotide
      variants (SNVs) for your chosen gene."),

      # Input: Select protein or genomic form ----
      radioButtons("form2", "Form of variant 2's data",
                   choices = c("Protein" = TRUE,
                               "Genomic" = FALSE),
                   selected = character(0)),

      # TODO: Add blurb about what it means to be in protein and genomic form

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
                  tabPanel("Variant 1", plotOutput("variantPlot1")),
                  tabPanel("Variant 2", plotOutput("variantPlot2")),
                  tabPanel("Overlapped", tableOutput("contents"))
      )
    )
  )
)


# TODO: Change the server logic
# Define server logic
server <- function(input, output) {

  print(input)
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

  output$variantPlot1 <- renderPlot({
    req(input$eveData)
    req(input$variant1)
    req(input$form1)
    req(input$geneName)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
    scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
    visualizeVariant(scoredVariant1, input$geneName)
  })

  output$variantPlot2 <- renderPlot({
    req(input$eveData)
    req(input$variant2)
    req(input$form2)
    req(input$geneName)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
    scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
    visualizeVariant(scoredVariant2, input$geneName)
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
