# Script for Shiny app

library("shiny")

# Define UI

ui <- fluidPage(

  # App title ----
  titlePanel("Score and visualize a variant's EVE scores"),

  # TODO: Add description about EVE and overall app

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
      radioButtons("form1", "Form of variant 1's data",
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

      actionButton(inputId = "button1",
                   label = "Run")
    ),


      # TODO: Add blurb about what it means to be in protein and genomic form

      # Horizontal line ----
      # tags$hr(),
      #
      # # TODO: Change so its a slider to select gene range
      # # Input: Slider for the number of observations to generate ----
      # sliderInput("n",
      #             "Number of observations:",
      #             value = 500,
      #             min = 1,
      #             max = 1000)
      #
      # ),


    # Main panel for displaying outputs ----
    mainPanel(

      # add filter to non-zero only
      # add average eve score

      # Output: Tabset w/ variant 1 plot, variant 2 plot, and two variants
      # overlapped ----
      tabsetPanel(type = "tabs",
                  tabPanel("About/Help",
                           tags$h3("Welcome to variantMapper's Shiny app!"),
                           tags$p("This Shiny app aims to make .... with EVE
                                  data easier"),
                           tags$h3("Instructions"),
                           tags$p("To use this Shiny app...")
                           ),
                  tabPanel("Variant 1",
                           plotOutput("variantPlot1"),
                           uiOutput("sliderValues1"),
                           h4("Average EVE score"),
                           textOutput("averageEveScore1"),
                           h4("Variant 1's information"),
                           tableOutput("variant1Data")
                           ),
                  tabPanel("Variant 2", plotOutput("variantPlot2"), uiOutput("sliderValues2")),
                  tabPanel("Overlapped", plotOutput("multiVariantPlot"), uiOutput("sliderValues3"))
      )
    )
  )
)


# TODO: Add reactive
# Define server logic
server <- function(input, output) {

  # change so that uses eventReactive and will only display if it is there

  print(input)

  output$sliderValues1 <- renderUI({
    req(input$eveData)
    req(input$variant1)
    req(input$form1)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
    scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
    firstRes <- min(scoredVariant1[, "resPos"])
    lastRes <- max(scoredVariant1[, "resPos"])
    sliderInput("numRes1", "Protein residue range", value = c(firstRes, lastRes), min = firstRes, max = lastRes)
  })
  output$sliderValues2 <- renderUI({
    req(input$eveData)
    req(input$variant2)
    req(input$form2)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
    scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
    firstRes <- min(scoredVariant2[, "resPos"])
    lastRes <- max(scoredVariant2[, "resPos"])
    sliderInput("numRes2", "Protein residue range", value = c(firstRes, lastRes), min = firstRes, max = lastRes)
  })
  output$sliderValues3 <- renderUI({
    req(input$eveData)
    req(input$variant2)
    req(input$form2)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
    scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
    firstRes <- min(scoredVariant2[, "resPos"])
    lastRes <- max(scoredVariant2[, "resPos"])
    sliderInput("numRes3", "Protein residue range", value = c(firstRes, lastRes), min = firstRes, max = lastRes)
  })
  output$variantPlot1 <- renderPlot({
    req(input$eveData)
    req(input$variant1)
    req(input$form1)
    req(input$geneName)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
    scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
    filteredScoredVariant1 <- dplyr::filter(scoredVariant1, resPos >= input$numRes1[1] & resPos <= input$numRes1[2])
    visualizeVariant(filteredScoredVariant1, input$geneName)
  })

  output$variantPlot2 <- renderPlot({
    req(input$eveData)
    req(input$variant2)
    req(input$form2)
    req(input$geneName)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
    scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
    filteredScoredVariant2 <- dplyr::filter(scoredVariant2, resPos >= input$numRes2[1] & resPos <= input$numRes2[2])
    visualizeVariant(filteredScoredVariant2, input$geneName)
  })

  output$multiVariantPlot <- renderPlot({
    req(input$eveData)
    req(input$variant1)
    req(input$form1)
    req(input$variant2)
    req(input$form2)
    req(input$geneName)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
    scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
    variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
    scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
    filteredScoredVariant1 <- dplyr::filter(scoredVariant1, resPos >= input$numRes3[1] & resPos <= input$numRes3[2])
    filteredScoredVariant2 <- dplyr::filter(scoredVariant2, resPos >= input$numRes3[1] & resPos <= input$numRes3[2])
    visualizeVariant2(filteredScoredVariant1, filteredScoredVariant2, input$geneName)
  })

  output$averageEveScore1 <- renderText({
    req(input$eveData)
    req(input$variant1)
    req(input$form1)
    req(input$geneName)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
    scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
    filteredScoredVariant1 <- dplyr::filter(scoredVariant1, resPos >= input$numRes1[1] & resPos <= input$numRes1[2])
    # header <- paste("Average EVE score from residue", input$numRes1[1], "to", input$numRes1[1], sep = " ")
    # h4(header)
    avgEveScore1 <- scoreVariant(filteredScoredVariant1$eveScores)
  })

  output$variant1Data <- renderTable({
    req(input$eveData)
    req(input$variant1)
    req(input$form1)
    req(input$geneName)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
    scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
    filteredScoredVariant1 <- dplyr::filter(scoredVariant1, resPos >= input$numRes1[1] & resPos <= input$numRes1[2])
  })
}

# Create Shiny app
shiny::shinyApp(ui, server)

# [END]
