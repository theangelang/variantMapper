# Script for Shiny app

library("shiny")

# Define UI

ui <- fluidPage(

  # App title ----
  titlePanel("Score and visualize a variant's EVE scores"),
  helpText("For information about this Shiny app and how to use it, please see
           the About/Help tab."),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Select a vcf file for EVE data ----
      fileInput("eveData", "EVE Data",
                multiple = FALSE,
                accept = ".vcf"),

      tags$p("Upload a vcf file from https://evemodel.org/ corresponding to the
             gene you would like to analyze.  The vcf file for the gene of
             interest can be found at https://evemodel.org/download/protein and
             searching for the gene."),

      tags$p("For example data please see the About/Help tab for more
             information."),

      # Input: Enter gene name ----
      textInput("geneName", "Gene name", placeholder = "Name of gene"),

      tags$p("Enter in the name of the gene of interest."),

      # Horizontal line ----
      tags$hr(),

      # Input: Select a csv file for variant 1----
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

      # Input: Select a csv file for variant 2----
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

    ),


    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ about/help section, variant 1 information, variant 2
      # information, and information for two variants simultaneously ----
      tabsetPanel(type = "tabs",
                  # About/Help tab ----
                  tabPanel("About/Help",
                           tags$h3("Welcome to variantMapper's Shiny app!"),
                           tags$p("This Shiny app aims to make variant
                                  pathogenicity predictions from the
                                  Evolutionary model of Variant Effect (EVE)
                                  more accessible allowing for one source of
                                  variant pathogenicity predictions to be easily
                                  accessed.  Moreover, this Shiny app allows for
                                  visualization of variants and their EVE
                                  scores."),
                          tags$h4("About EVE"),
                          tags$p("EVE is an unsupervised machine
                                  learning model shown to be accurate in its
                                  predictions and doesn’t rely on knowledge of
                                  protein function. It uses multiple sequence
                                  alignments to predict pathogenicity of
                                  missense variants. An EVE score is assigned to
                                  missense variants from a continuous interval
                                  of 0 to 1 with 0 being benign and 1 being most
                                  pathogenic.  You can find out more about
                                  EVE here (https://evemodel.org/)."),
                          tags$p("With this Shiny app you will be able to
                                  assign an EVE score to your variants of
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
                                           this Shiny app.  The minimum is
                                           filling out the 'EVE Data', 'Gene
                                           name' and inputting the information
                                           for one variant.  The variant
                                           information will show up on the
                                           tab corresponding to the input
                                           location.  More specifically,
                                           the minimum required fields to use
                                           this Shiny app are either:"),
                          tags$ul(
                            tags$li("'EVE Data', 'Gene name', 'Variant 1
                                    data', and 'Form of variant 1's data'"),
                            tags$li("'EVE Data', 'Gene name', 'Variant 2
                                             data', and 'Form of variant 2's
                                             data'")
                          ),
                          tags$blockquote("Inputing variant information for two
                                          variants enables you to look at the
                                          EVE scores for two variants separately
                                          and overlapped onto the gene at one
                                          time."),
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
                  # Variant 1 tab ----
                  tabPanel("Variant 1",
                           uiOutput("variant1PlotTitle"),
                           plotOutput("variantPlot1"),
                           uiOutput("sliderValues1"),
                           uiOutput("averageEveScore1"),
                           uiOutput("variant1TableTitleDesc"),
                           tableOutput("variant1Data")
                           ),
                  # Variant 2 tab ----
                  tabPanel("Variant 2",
                           uiOutput("variant2PlotTitle"),
                           plotOutput("variantPlot2"),
                           uiOutput("sliderValues2"),
                           uiOutput("averageEveScore2"),
                           uiOutput("variant2TableTitleDesc"),
                           tableOutput("variant2Data")
                          ),
                  # Variant 1 and 2 displayed simultaneously tab ----
                  tabPanel("Overlapped",
                           uiOutput("variantMultiPlotTitle"),
                           plotOutput("variantPlotMulti"),
                           uiOutput("sliderValuesMulti"),
                           uiOutput("averageScoreMultiTitle"),
                           fluidRow(column(6,
                                           uiOutput("averageEveScoreMulti1")),
                                    column(6,
                                           uiOutput("averageEveScoreMulti2"))),
                           uiOutput("variantMultiTableTitleDesc"),
                           fluidRow(column(6,
                                           uiOutput("variant1MultiTableTitle")),
                                    column(6,
                                           uiOutput("variant2MultiTableTitle"))),
                           fluidRow(column(6,
                                           tableOutput("variant1MultiData")),
                                    column(6,
                                           tableOutput("variant2MultiData")))
                  )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {

  #----------- For the Variant 1 tab

  variant1Calculation <- reactive({
    req(input$eveData)
    req(input$variant1)
    req(input$form1)
    req(input$geneName)
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant1Processed <- processVariantData(filePath = input$variant1$datapath,
                                            as.logical(input$form1))
    scoredVariant1 <- getEveScores(eveDataProcessed, variant1Processed,
                                   as.logical(input$form1))
    scoredVariant1
  })

  slider1Calculation <- reactive({
    if (! is.null(variant1Calculation)) {
      variant1Data <- variant1Calculation()
      firstRes <- min(variant1Data[, "resPos"])
      lastRes <- max(variant1Data[, "resPos"])
      c(firstRes, lastRes)
    }
  })

  output$sliderValues1 <- renderUI({
    if (! is.null(slider1Calculation)) {
      slider1Values <- slider1Calculation()
      sliderInput("numRes1",
                  "Protein residue range",
                  value = c(slider1Values[[1]], slider1Values[[2]]),
                  min = slider1Values[[1]], max = slider1Values[[2]])
      }
    })

  variant1Filtering <- reactive({
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
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))) {
      variant1Data <- variant1Calculation()
      filteredVariant1Data <- variant1Filtering()
      visualizeVariant(filteredVariant1Data, input$geneName)
    }
  })

  output$averageEveScore1 <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      filteredVariant1Data <- variant1Filtering()
      header <- paste("Average EVE score from residue",
                      input$numRes1[1],
                      "to",
                      input$numRes1[2],
                      "for",
                      input$geneName,
                      sep = " ")
      tagList(
        tags$h4(header),
        scoreVariant(filteredVariant1Data$eveScores)
      )
    }
  })

  output$variant1TableTitleDesc <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1Filtering))){
      header <- paste("Variant 1 of ",
                      input$geneName,
                      "'s EVE score and amino acid information",
                      sep = "")
      tagList(
        tags$h4(header),
        tags$p("Below is a description of column titles:"),
        tags$ul(
          tags$li("eveScores: The EVE score of the amino acid at this residue
                  position in the variant."),
          tags$li("resPos: The residue position."),
          tags$li("wtAa: The wildtype amino acid at this residue position."),
          tags$li("varPos: The variant amino acid at this residue position in
                  the protein variant.")
        )
      )
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
    eveDataProcessed <- processEveData(filePath = input$eveData$datapath)
    variant2Processed <- processVariantData(filePath = input$variant2$datapath,
                                            as.logical(input$form2))
    scoredVariant2 <- getEveScores(eveDataProcessed,
                                   variant2Processed,
                                   as.logical(input$form2))
    scoredVariant2
  })

  slider2Calculation <- reactive({
    variant2Data <- variant2Calculation()
    firstRes <- min(variant2Data[, "resPos"])
    lastRes <- max(variant2Data[, "resPos"])
    return(c(firstRes, lastRes))
  })

  output$sliderValues2 <- renderUI({
    slider2Values <- slider2Calculation()
    sliderInput("numRes2",
                "Protein residue range",
                value = c(slider2Values[[1]], slider2Values[[2]]),
                min = slider2Values[[1]], max = slider2Values[[2]])
  })

  variant2Filtering <- reactive({
    req(input$numRes2)
    scoredVariant2 <- variant2Calculation()
    filteredScoredVariant2 <- dplyr::filter(scoredVariant2,
                                            resPos >= input$numRes2[1] &
                                              resPos <= input$numRes2[2])
  })

  output$variant2PlotTitle <- renderUI({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      header <- paste("Plot of EVE score vs residue position for",
                      input$geneName,
                      sep = " ")
      h4(header)
    }
  })

  output$variantPlot2 <- renderPlot({
    variant2Data <- variant2Calculation()
    filteredVariant2Data <- variant2Filtering()
    visualizeVariant(filteredVariant2Data, input$geneName)
  })

  output$averageEveScore2 <- renderUI({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      filteredVariant2Data <- variant2Filtering()
      header <- paste("Average EVE score from residue",
                      input$numRes2[1],
                      "to",
                      input$numRes2[2],
                      "for",
                      input$geneName,
                      sep = " ")
      tagList(
        tags$h4(header),
        scoreVariant(filteredVariant2Data$eveScores)
      )
    }
  })

  output$variant2TableTitleDesc <- renderUI({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      header <- paste("Variant 2 of ",
                      input$geneName,
                      "'s EVE score and amino acid information",
                      sep = "")
      tagList(
        tags$h4(header),
        tags$p("Below is a description of column titles:"),
        tags$ul(
          tags$li("eveScores: The EVE score of the amino acid at this residue
                  position in the variant."),
          tags$li("resPos: The residue position."),
          tags$li("wtAa: The wildtype amino acid at this residue position."),
          tags$li("varPos: The variant amino acid at this residue position in
                  the protein variant.")
        )
      )
    }
  })

  output$variant2Data <- renderTable({
    if ((! is.null(variant2Calculation())) & (! is.null(variant2Filtering))){
      filteredVariant2Data <- variant2Filtering()
      filteredVariant2Data
    }
  })

  #----------- For the overlap tab

  output$sliderValuesMulti <- renderUI({
    # use values from variant 2 since need both variants for this page

    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      slider2Values <- slider2Calculation()
      sliderInput("numResMulti",
                "Protein residue range",
                value = c(slider2Values[[1]], slider2Values[[2]]),
                min = slider2Values[[1]], max = slider2Values[[2]])
    }
  })

  variant1FilteringMulti <- reactive({
    req(input$numResMulti)
    scoredVariant1 <- variant1Calculation()
    filteredScoredVariant1Multi <- dplyr::filter(scoredVariant1,
                                                 resPos >= input$numResMulti[1] &
                                                   resPos <= input$numResMulti[2])
  })

  variant2FilteringMulti <- reactive({
    req(input$numResMulti)
    scoredVariant2 <- variant2Calculation()
    filteredScoredVariant2Multi <- dplyr::filter(scoredVariant2,
                                            resPos >= input$numResMulti[1] &
                                              resPos <= input$numResMulti[2])
  })

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

  output$averageEveScoreMulti2 <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      filteredVariant2Data <- variant2FilteringMulti()
      tagList(
        tags$h5("Variant 2"),
        scoreVariant(filteredVariant2Data$eveScores)
      )
    }
  })

  output$averageEveScoreMulti1 <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      filteredVariant1Data <- variant1FilteringMulti()
      tagList(
        tags$h5("Variant 1"),
        scoreVariant(filteredVariant1Data$eveScores)
      )
    }
  })

  output$variantMultiTableTitleDesc <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      header <- paste("Variants 1 and 2 of ",
                      input$geneName,
                      "'s EVE score and amino acid information",
                      sep = "")
      tagList(
        h4(header),
        tags$p("Below is a description of column titles:"),
        tags$ul(
          tags$li("eveScores: The EVE score of the amino acid at this residue
                  position in the variant."),
          tags$li("resPos: The residue position."),
          tags$li("wtAa: The wildtype amino acid at this residue position."),
          tags$li("varPos: The variant amino acid at this residue position in
                  the protein variant.")
        )
      )
    }
  })

  output$variant1MultiTableTitle <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      tags$h5("Variant 1")
    }
  })

  output$variant2MultiTableTitle <- renderUI({
    if ((! is.null(variant1Calculation())) & (! is.null(variant1FilteringMulti)) &
        (! is.null(variant2Calculation())) & (! is.null(variant2FilteringMulti))){
      tags$h5("Variant 2")
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

}

# Create Shiny app
shiny::shinyApp(ui, server)

# [END]
