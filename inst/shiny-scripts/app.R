# Script for Shiny app

library("shiny")

# Define UI

ui <- fluidPage(

  # App title ----
  titlePanel("Score and visualize a variant's EVE scores"),
  helpText("For information about this Shiny app and how to use it, please see
           the About/Help tab."),

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

      tags$p("For example data please see the About/Help tab for more
             information."),

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

      tags$p("For example data and details about file format please see the
             About/Help tab."),

      # Input: Select protein or genomic form ----
      radioButtons("form1", "Form of variant 1's data",
                   choices = c("Protein" = TRUE,
                               "Genomic" = FALSE),
                   selected = character(0)),

      tags$p("Protein form indicates the SNVs are listed as missense mutations
      at particular resiude positions.  Genomic form has SNVs listed with their
      genomic position and nucleotide changes.  The file types for 'Protein' and
      'Genomic' form are different, please see the About/Help tab for more
      information about file format and how to determine which type of file you
      have."),

      # Horizontal line ----
      tags$hr(),

      # Input: Select a file ----
      fileInput("variant2", "Variant 2 data (optional)",
                multiple = FALSE,
                accept = ".csv"),

      tags$p("Optionally, upload a second csv file containing single nucleotide
      variants (SNVs) for your chosen gene."),

      tags$p("For example data and details about file format please see the
             About/Help tab."),

      # Input: Select protein or genomic form ----
      radioButtons("form2", "Form of variant 2's data",
                   choices = c("Protein" = TRUE,
                               "Genomic" = FALSE),
                   selected = character(0)),

      tags$p("Protein form indicates the SNVs are listed as missense mutations
      at particular resiude positions.  Genomic form has SNVs listed with their
      genomic position and nucleotide changes.  The file types for 'Protein' and
      'Genomic' form are different, please see the About/Help tab for more
      information about file format and how to determine which type of file you
      have."),

      actionButton(inputId = "runCalculations",
                   label = "Run")
    ),


    # Main panel for displaying outputs ----
    mainPanel(

      # add filter to non-zero only
      # add average eve score

      # Output: Tabset w/ variant 1 plot, variant 2 plot, and two variants
      # overlapped ----
      tabsetPanel(type = "tabs",
                  tabPanel("About/Help",
                           tags$h3("Welcome to variantMapper's Shiny app!"),
                           tags$p("This Shiny app aims to make variant
                                  pathogenicity predictions from the
                                  Evolutionary model of Variant Effect (EVE)
                                  more accessible allowing for one source of
                                  variant pathogenicity predictions to be easily
                                  accessed. EVE is an unsupervised machine
                                  learning model shown to be accurate in its
                                  predictions and doesn’t rely on knowledge of
                                  protein function. It uses multiple sequence
                                  alignments to predict pathogenicity of
                                  missense variants. You can find out more about
                                  EVE here (https://evemodel.org/).
                                  With this Shiny app you will be able to
                                  asign an EVE score to your variants of
                                  interest and visualize their EVE scores by
                                  residue position.  You can also compare two
                                  variants of the same gene simultaneously on
                                  one plot."),
                           tags$h4("Variant 1 tab"),
                           tags$p("The Variant 1 tab corresponds to analysis and
                                  visualization from the data entered in the
                                  Variant 1 section in the input panel on the
                                  left."),
                           tags$h4("Variant 2 tab"),
                           tags$p("The Variant 2 tab corresponds to analysis and
                                  visualization from the data entered in the
                                  Variant 2 section in the input panel on the
                                  left."),
                           tags$h4("Overlapped tab"),
                           tags$p("The Overlapped panel shows analysis and
                                  visualizations of both variants
                                  ssimultaneously."),
                           tags$h3("Instructions"),

                           tags$blockquote("Note, you do not have to upload two
                                           sets of variant data in order to use
                                           this Shiny app.  The minimum required
                                           fields to use this Shiny app are
                                           'EVE Data', 'Gene name', 'Variant 1
                                           data', and 'Form of variant 1's
                                           data'.  The second variant option
                                           enables you to look at the EVE scores
                                           for two variants separately and
                                           overlapped onto the gene at one time.
                                           "),
                           tags$h4("To use this Shiny app to analyze one
                                   variant:"),
                           tags$ol(
                             tags$li("Upload EVE Data."),
                             tags$li("Enter Gene name."),
                             tags$li("Upload csv file with SNVs."),
                             tags$li("Select the form of variant 1's data i.e.
                                     Protein or Genomic.")
                           ),
                           tags$h4("To use this Shiny app to analyze two
                                  variants:"),
                           tags$ol(
                             tags$li("Do steps 1 to 4 in the above section for
                                     analyzing one variant."),
                             tags$li("Upload csv file with SNVs."),
                             tags$li("Select the form of variant 2's data i.e.
                                     Protein or Genomic.")
                           ),
                           tags$blockquote("Once all the fields are filled out
                           the calculations will begin and you can navigate
                           to the appropriate page.  If any of the fields
                           in the side panel are changed the tabs will be
                                           updated accordingly as well."),
                           tags$h3("Example data"),
                           tags$h4("EVE Data"),
                           tags$p("EVE Data is a vcf file downloaded from
                           (https://evemodel.org/).  A sample of what this file
                           should look like is 'NRX1B_HUMAN_SUBSET.vcf'.  It can
                           be found at
                           (https://github.com/theangelang/variantMapper/tree/
                           master/inst/extdata)."),
                           tags$p("Please note this is a subset of NRXN1 and
                           does not contain all the rows."),
                           tags$h4("Protein form"),
                           tags$p("Protein form is a csv file with three
                           columns.The first being 'wt_aa' with the single amino
                           acid code for the wild type amino acid.  The second
                           is 'position' containing integers representing
                           residue position.  The third column is 'mt_aa'which
                           is the amino acid variant at this particular residue
                           position.  It also contains the one letter amino acid
                           code.  For an example please see
                           'variant_data_protein.csv' at https://github.com
                           /theangelang/variantMapper/tree/master/inst/extdata"
                                  ),
                           tags$h4("Genomic form"),
                           tags$p("Genomic form is a csv file with five columns.
                           The first being 'chrom' which are either integers or
                           X or Y representing chromosome location.  The second
                           being 'start' which has integers indicating start
                           position of the SNV.  The third is 'end', containing
                           integers representing the end position of the SNV.
                           The fourth column is 'ref_allele' which has
                           characters representing the wildtype nucleotide.
                           The fifth column is 'alt_allele' which has characters
                           representing the variant nucleotide.  For an example
                           please see 'variant_data_genomic.csv' at
                                  https://github.com/theangelang/variantMapper/
                                  tree/master/inst/extdata")
                           ),
                  tabPanel("Variant 1",
                           uiOutput("variant1PlotTitle"),
                           plotOutput("variantPlot1"),
                           uiOutput("sliderValues1"),
                           uiOutput("averageScore1Title"),
                           textOutput("averageEveScore1"),
                           uiOutput("variant1TableTitle"),
                           uiOutput("variant1TableDesc"),
                           tableOutput("variant1Data")
                           ),
                  tabPanel("Variant 2",
                           uiOutput("variant2PlotTitle"),
                           plotOutput("variantPlot2"),
                           uiOutput("sliderValues2"),
                           uiOutput("averageScore2Title"),
                           textOutput("averageEveScore2"),
                           uiOutput("variant2TableTitle"),
                           uiOutput("variant2TableDesc"),
                           tableOutput("variant2Data")
                          ),
                  tabPanel("Overlapped",
                           uiOutput("variantMultiPlotTitle"),
                           plotOutput("variantPlotMulti"),
                           uiOutput("sliderValuesMulti"),
                           uiOutput("averageScoreMultiTitle"),
                           # have 2 columns here, add more labels
                           fluidRow(column(6,
                                           textOutput("averageEveScoreMulti1")),
                                    column(6,
                                           textOutput("averageEveScoreMulti2"))),
                           uiOutput("variantMultiTableTitle"),
                           uiOutput("variantMultiTableDesc"),
                           # have 2 columns here, add more labels
                           fluidRow(column(6,
                                           tableOutput("variant1MultiData")),
                                    column(6,
                                           tableOutput("variant2MultiData")))
                  )
      )
    )
  )
)


# TODO: Add reactive and logic for using button
# Define server logic
server <- function(input, output, session) {

  #----------- For the Variant 1 tab

  variant1Calculation <- reactive({
    req(input$eveData)
    req(input$variant1)
    req(input$form1)
    req(input$geneName)
    # req(input$runCalculations)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath,
                                            as.logical(input$form1))
    scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed,
                                   as.logical(input$form1))
    scoredVariant1
    # firstRes <- min(scoredVariant1[, "resPos"])
    # lastRes <- max(scoredVariant1[, "resPos"])
  })

  slider1Calculation <- reactive({
    # req(input$runCalculations)
    variant1Data <- variant1Calculation()
    firstRes <- min(variant1Data[, "resPos"])
    lastRes <- max(variant1Data[, "resPos"])
    c(firstRes, lastRes)
  })

  output$sliderValues1 <- renderUI({
    # variant1Data <- variant1Calculation()
    # firstRes <- min(variant1Data[, "resPos"])
    # lastRes <- max(variant1Data[, "resPos"])
    # req(input$runCalculations)
    slider1Values <- slider1Calculation()
    print(slider1Values)
    sliderInput("numRes1",
                "Protein residue range",
                value = c(slider1Values[[1]], slider1Values[[2]]),
                min = slider1Values[[1]], max = slider1Values[[2]])
    })

  variant1Filtering <- reactive({
    # req(input$runCalculations)
    req(input$numRes1)
    scoredVariant1 <- variant1Calculation()
    filteredScoredVariant1 <- dplyr::filter(scoredVariant1,
                                            resPos >= input$numRes1[1] &
                                              resPos <= input$numRes1[2])
  })

  output$variant1PlotTitle <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      header <- paste("Plot of EVE score vs residue position for",
                      input$geneName,
                      sep = " ")
      h4(header)
    }
  })

  output$variantPlot1 <- renderPlot({
    # req(input$runCalculations)
    variant1Data <- variant1Calculation()
    filteredVariant1Data <- variant1Filtering()
    visualizeVariant(filteredVariant1Data, input$geneName)
  })

  output$averageScore1Title <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      h4("Average EVE score")
      header <- paste("Average EVE score from residue",
                      input$numRes1[1],
                      "to",
                      input$numRes1[2],
                      "for",
                      input$geneName,
                      sep = " ")
      h4(header)
    }
  })

  output$averageEveScore1 <- renderText({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      filteredVariant1Data <- variant1Filtering()
      scoreVariant(filteredVariant1Data$eveScores)
    }
  })

  output$variant1TableTitle <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      header <- paste("Variant 1 of ",
                      input$geneName,
                      "'s EVE score and amino acid information",
                      sep = "")
      h4(header)
    }
  })

  output$variant1TableDesc <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      p("The columns...")
    }
  })

  output$variant1Data <- renderTable({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      filteredVariant1Data <- variant1Filtering()
      filteredVariant1Data
    }
  })

  #----------- For the Variant 2 tab

  variant2Calculation <- reactive({
    req(input$eveData)
    req(input$variant2)
    req(input$form2)
    req(input$geneName)
    # req(input$runCalculations)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant2Processed <- processVariantData(filePath = input$variant2$datapath,
                                            as.logical(input$form2))
    scoredVariant2 <- getEveScores(eveDataProcessed,
                                   variant2Processed,
                                   as.logical(input$form2))
    scoredVariant2
    # firstRes <- min(scoredVariant1[, "resPos"])
    # lastRes <- max(scoredVariant1[, "resPos"])
  })

  slider2Calculation <- reactive({
    # req(input$runCalculations)
    variant2Data <- variant2Calculation()
    firstRes <- min(variant2Data[, "resPos"])
    lastRes <- max(variant2Data[, "resPos"])
    return(c(firstRes, lastRes))
  })

  output$sliderValues2 <- renderUI({
    # variant1Data <- variant1Calculation()
    # firstRes <- min(variant1Data[, "resPos"])
    # lastRes <- max(variant1Data[, "resPos"])
    # req(input$runCalculations)
    slider2Values <- slider2Calculation()
    print(slider2Values)
    sliderInput("numRes2",
                "Protein residue range",
                value = c(slider2Values[[1]], slider2Values[[2]]),
                min = slider2Values[[1]], max = slider2Values[[2]])
  })

  variant2Filtering <- reactive({
    # req(input$runCalculations)
    req(input$numRes2)
    scoredVariant2 <- variant2Calculation()
    filteredScoredVariant2 <- dplyr::filter(scoredVariant2,
                                            resPos >= input$numRes2[1] &
                                              resPos <= input$numRes2[2])
  })

  # maybe don't need
  output$variant2PlotTitle <- renderUI({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      header <- paste("Plot of EVE score vs residue position for",
                      input$geneName,
                      sep = " ")
      h4(header)
    }
  })

  output$variantPlot2 <- renderPlot({
    # req(input$runCalculations)
    variant2Data <- variant2Calculation()
    filteredVariant2Data <- variant2Filtering()
    visualizeVariant(filteredVariant2Data, input$geneName)
  })

  output$averageScore2Title <- renderUI({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      h4("Average EVE score")
      header <- paste("Average EVE score from residue",
                      input$numRes2[1],
                      "to",
                      input$numRes2[2],
                      "for",
                      input$geneName,
                      sep = " ")
      h4(header)
    }
  })

  output$averageEveScore2 <- renderText({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      filteredVariant2Data <- variant2Filtering()
      scoreVariant(filteredVariant2Data$eveScores)
    }
  })

  output$variant2TableTitle <- renderUI({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      header <- paste("Variant 2 of ",
                      input$geneName,
                      "'s EVE score and amino acid information",
                      sep = "")
      h4(header)
    }
  })

  output$variant2TableDesc <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      p("The columns...")
    }
  })

  output$variant2Data <- renderTable({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      filteredVariant2Data <- variant2Filtering()
      filteredVariant2Data
    }
  })

  #----------- For the overlap tab

  # sliderMultiCalculation <- reactive({
  #   # req(input$runCalculations)
  #   variant2Data <- variant2Calculation()
  #   firstRes <- min(variant2Data[, "resPos"])
  #   lastRes <- max(variant2Data[, "resPos"])
  #   return(c(firstRes, lastRes))
  # })

  output$sliderValuesMulti <- renderUI({
    # variant1Data <- variant1Calculation()
    # firstRes <- min(variant1Data[, "resPos"])
    # lastRes <- max(variant1Data[, "resPos"])
    # req(input$runCalculations)
    # use values from variant 2 since need both variants for this page

    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      slider2Values <- slider2Calculation()
      print(slider2Values)
      sliderInput("numResMulti",
                "Protein residue range",
                value = c(slider2Values[[1]], slider2Values[[2]]),
                min = slider2Values[[1]], max = slider2Values[[2]])
    }
  })

  variant1FilteringMulti <- reactive({
    # req(input$runCalculations)
    req(input$numResMulti)
    scoredVariant1 <- variant1Calculation()
    filteredScoredVariant1Multi <- dplyr::filter(scoredVariant1,
                                                 resPos >= input$numResMulti[1] &
                                                   resPos <= input$numResMulti[2])
  })

  variant2FilteringMulti <- reactive({
    # req(input$runCalculations)
    req(input$numResMulti)
    scoredVariant2 <- variant2Calculation()
    filteredScoredVariant2Multi <- dplyr::filter(scoredVariant2,
                                            resPos >= input$numResMulti[1] &
                                              resPos <= input$numResMulti[2])
  })

  # maybe don't need
  output$variantMultiPlotTitle <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      header <- paste("Plot of EVE score vs residue position for Variant 1 and Variant 2 of",
                      input$geneName,
                      sep = " ")
      h4(header)
    }
  })

  output$variantPlotMulti <- renderPlot({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      filteredVariant1Data <- variant1FilteringMulti()
      filteredVariant2Data <- variant2FilteringMulti()
      visualizeVariant2(filteredVariant1Data,
                        filteredVariant2Data,
                        input$geneName)
    }
  })

  output$averageScoreMultiTitle <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      h4("Average EVE score")
      header <- paste("Average EVE score from residue",
                      input$numRes2[1],
                      "to",
                      input$numRes2[2],
                      "for",
                      input$geneName,
                      sep = " ")
      h4(header)
    }
  })

  output$averageEveScoreMulti2 <- renderText({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      filteredVariant2Data <- variant2FilteringMulti()
      scoreVariant(filteredVariant2Data$eveScores)
    }
  })

  output$averageEveScoreMulti1 <- renderText({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      filteredVariant1Data <- variant1FilteringMulti()
      scoreVariant(filteredVariant1Data$eveScores)
    }
  })

  output$variantMultiTableTitle <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      header <- paste("Variants 1 and 2 of ",
                      input$geneName,
                      "'s EVE score and amino acid information",
                      sep = "")
      h4(header)
    }
  })

  output$variantMultiTableDesc <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      p("The columns...")
    }
  })

  output$variant1MultiData <- renderTable({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      filteredVariant1Data <- variant1FilteringMulti()
      filteredVariant1Data
    }
  })

  output$variant2MultiData <- renderTable({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      filteredVariant2Data <- variant2FilteringMulti()
      filteredVariant2Data
    }
  })



  print(input)

  # output$sliderValues1 <- renderUI({
  #   req(input$eveData)
  #   req(input$variant1)
  #   req(input$form1)
  #   eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
  #   variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
  #   scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
  #   firstRes <- min(scoredVariant1[, "resPos"])
  #   lastRes <- max(scoredVariant1[, "resPos"])
  #   sliderInput("numRes1", "Protein residue range", value = c(firstRes, lastRes), min = firstRes, max = lastRes)
  # })
  # output$sliderValues2 <- renderUI({
  #   req(input$eveData)
  #   req(input$variant2)
  #   req(input$form2)
  #   eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
  #   variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
  #   scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
  #   firstRes <- min(scoredVariant2[, "resPos"])
  #   lastRes <- max(scoredVariant2[, "resPos"])
  #   sliderInput("numRes2", "Protein residue range", value = c(firstRes, lastRes), min = firstRes, max = lastRes)
  # })
  # output$sliderValues3 <- renderUI({
  #   req(input$eveData)
  #   req(input$variant2)
  #   req(input$form2)
  #   eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
  #   variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
  #   scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
  #   firstRes <- min(scoredVariant2[, "resPos"])
  #   lastRes <- max(scoredVariant2[, "resPos"])
  #   sliderInput("numRes3", "Protein residue range", value = c(firstRes, lastRes), min = firstRes, max = lastRes)
  # })
  # output$variantPlot1 <- renderPlot({
  #   req(input$eveData)
  #   req(input$variant1)
  #   req(input$form1)
  #   req(input$geneName)
  #   eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
  #   variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
  #   scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
  #   filteredScoredVariant1 <- dplyr::filter(scoredVariant1, resPos >= input$numRes1[1] & resPos <= input$numRes1[2])
  #   visualizeVariant(filteredScoredVariant1, input$geneName)
  # })
  #
  # output$variantPlot2 <- renderPlot({
  #   req(input$eveData)
  #   req(input$variant2)
  #   req(input$form2)
  #   req(input$geneName)
  #   eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
  #   variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
  #   scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
  #   filteredScoredVariant2 <- dplyr::filter(scoredVariant2, resPos >= input$numRes2[1] & resPos <= input$numRes2[2])
  #   visualizeVariant(filteredScoredVariant2, input$geneName)
  # })
  #
  # output$multiVariantPlot <- renderPlot({
  #   req(input$eveData)
  #   req(input$variant1)
  #   req(input$form1)
  #   req(input$variant2)
  #   req(input$form2)
  #   req(input$geneName)
  #   eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
  #   variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
  #   scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
  #   variant2Processed <- processVariantData(filePath = input$variant2$datapath, as.logical(input$form2))
  #   scoredVariant2 <- getEveScores(eveDataProcessed, variant2Processed, as.logical(input$form2))
  #   filteredScoredVariant1 <- dplyr::filter(scoredVariant1, resPos >= input$numRes3[1] & resPos <= input$numRes3[2])
  #   filteredScoredVariant2 <- dplyr::filter(scoredVariant2, resPos >= input$numRes3[1] & resPos <= input$numRes3[2])
  #   visualizeVariant2(filteredScoredVariant1, filteredScoredVariant2, input$geneName)
  # })
  #
  # output$averageEveScore1 <- renderText({
  #   req(input$eveData)
  #   req(input$variant1)
  #   req(input$form1)
  #   req(input$geneName)
  #   eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
  #   variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
  #   scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
  #   filteredScoredVariant1 <- dplyr::filter(scoredVariant1, resPos >= input$numRes1[1] & resPos <= input$numRes1[2])
  #   # header <- paste("Average EVE score from residue", input$numRes1[1], "to", input$numRes1[1], sep = " ")
  #   # h4(header)
  #   avgEveScore1 <- scoreVariant(filteredScoredVariant1$eveScores)
  # })
  #
  # output$variant1Data <- renderTable({
  #   req(input$eveData)
  #   req(input$variant1)
  #   req(input$form1)
  #   req(input$geneName)
  #   eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
  #   variant1Processed <- processVariantData(filePath = input$variant1$datapath, as.logical(input$form1))
  #   scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed, as.logical(input$form1))
  #   filteredScoredVariant1 <- dplyr::filter(scoredVariant1, resPos >= input$numRes1[1] & resPos <= input$numRes1[2])
  # })
}

# Create Shiny app
shiny::shinyApp(ui, server)

# [END]
