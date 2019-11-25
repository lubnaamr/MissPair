#' A graphical user interface for the package MissPair
#'
#' This function provides a graphical user interface (GUI) for
#' calculating statistical tests in matched pairs with missing
#' values in a single arm or both arms.
#'
#' The function produces a GUI for the calculation of the test statistics.
#'  Data can be loaded via the "Browse" button. The formula,
#'  the missing pattern, testing hypothesis, testing method,
#'  number of permutations (default: 10,000), number of bootstraps
#'  (default: 1000) and the significance level alpha
#'  (default: 0.05) need to be specified.
#'
#' @aliases viewGUI
#'
#' @export
viewGUI <- function() {
  requireNamespace("shiny", quietly = TRUE)
  if (!("package:shiny" %in% search())) {
    attachNamespace("shiny")
  }


  ui <- fluidPage(theme = shinythemes::shinytheme("cerulean"),
                  shinyjs::useShinyjs(),
                  titlePanel("Tests for Matched Pairs with Missing Values"),
                  sidebarLayout(
                    sidebarPanel(
                      splitLayout(
                        fileInput("infile", "Choose CSV File",
                                  accept = c(
                                    "text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
                        checkboxInput("header", "Header", TRUE)

                      ),

                      splitLayout(

                      textInput("formula", "Formula of the form  ~x+y", "~ + "),

                      radioButtons("Missing", "Missing Pattern:",
                                   c("Both variables" = "botharms",
                                     "One variable" = "onearm"), inline = T)
                    ),
                      selectInput("Hypothesis", "Select Testing Hypothesis:",
                                  c("Testing Equality of Means (H0m)" = "h0m",
                                    "Testing H0p hypothesis" = "h0p")),

                      selectInput("Method", "Select Testing Method:",
                                  c("Weighted Permutation Test" = "wptest",
                                    "Multiplication combination test" = "mctest",
                                    "Bootsraping Paired data Test" = "ppd")),


                      splitLayout(

                        numericInput("nperm", "nperm/nbsp", value = 10000),

                        numericInput("alpha", "Alpha", value = 0.05, min = 0, max = 1),
                        radioButtons("alternative", "Alternative:",
                                     c("two.sided" = "two.sided",
                                       "less" = "less", "greater" = "greater"))
                      ),

                      actionButton("process", "Calculate", class = "btn-primary")
                    ),
                    mainPanel(

                      verbatimTextOutput("result")


                    )
                  )
  )









  server <- function(input, output, session) {

    datasetInput <- reactive({

      req(input$infile)

      if (is.null(input$infile))
        return(NULL)
      read.csv(input$infile$datapath, header = input$header)
    })


    observeEvent(input$Missing, {

      if (input$Missing == "onearm") {
        shinyjs::hide("alternative")
        updateNumericInput(session, "nperm", "nperm/nbsp", value = 1000)

        updateSelectInput(session, "Hypothesis",
                          label = paste("Testing Hypothesis:"),
                          choices = c("Testing Equality of Means (H0m)" = "h0m"))

        observeEvent(input$Hypothesis, {
          updateSelectInput(session, "Method",
                            label = paste("Testing Method:"),
                            choices = c("Bootsraping Paired data Test" = "ppd"))
        })## observeevent


      }


      if (input$Missing == "botharms") {
        shinyjs::show("alternative")
        updateNumericInput(session, "nperm", "nperm/nbsp", value = 10000)

        updateSelectInput(session, "Hypothesis",
                          label = paste("Select Testing Hypothesis:"),
                          choices = c("Testing Equality of Means (H0m)" = "h0m",
                                      "Testing H0p hypothesis" = "h0p"))
        observeEvent(input$Hypothesis, {
          if (input$Hypothesis == "h0m")
            updateSelectInput(session, "Method",
                              label = paste("Select Testing Method:"),
                              choices =  c("Weighted Permutation Test" = "wptest",
                                           "Multiplication combination test" = "mctest"))


          if (input$Hypothesis == "h0p")

            updateSelectInput(session, "Method",
                              label = paste("Testing Method:"),
                              choices =  c("Multiplication combination test" = "mctest"))

        })## observeevent
      }

    }) ## observeevent



    observeEvent(input$process, {
      if (input$formula != "~ + ") {
      if (all(all.vars(as.formula(input$formula)) %in% colnames(datasetInput()))) {
        data <- model.frame(as.formula(input$formula),
                            as.data.frame(datasetInput()), na.action = NULL)

      }
      }
if (input$formula == "~ + " || length(as.formula(input$formula)) != 2L || length(data) != 2L) {

  output$result <- renderPrint({
  "'formula' missing or invalid"
  })

} else {





      isolate(if (input$Missing == "onearm") {

        output$result <- renderPrint({
          testPBT <- withProgress(value = 0.8, expr = {
          PBT(isolate(c(data[, 1])),
              isolate(c(data[, 2])),
              isolate(input$nperm),
              isolate(input$alpha))
        }
        , message = "calculation is in progress...");
        list(PBT = testPBT$PBT, Descriptive = testPBT$Descriptive,
                                          Atrr = class(testPBT),
                                          Input = testPBT$Input) })


      } else {   #####Missing in both arms



        if (input$Method == "wptest") {

          output$result <- renderPrint({
            testWPT <- WPT(isolate(c(data[, 1])),
                           isolate(c(data[, 2])),
                           isolate(input$alternative),
                           isolate(input$nperm),
                           isolate(input$alpha));
          list(WPT = testWPT$WPT, Descriptive = testWPT$Descriptive,
               TestDecision = testWPT$TestDecision,
               Atrr = class(testWPT),
               Input = testWPT$Input) })



        } else { ## Method multiplication

          output$result <- renderPrint({
            testMCT <- MCT(isolate(c(data[, 1])),
                         isolate(c(data[, 2])),
                         isolate(input$Hypothesis),
                         isolate(input$alternative),
                         isolate(input$nperm),
                         isolate(input$alpha));
          list(MCT = testMCT$MCT, Descriptive = testMCT$Descriptive,
               TestDecision = testMCT$TestDecision,
               Atrr = class(testMCT),
               Input = testMCT$Input) })






        } ## end of method multiplication

      } ##end of missing in both arms
      ) ### end of isolate for input$method


} ## end of else checking formula

    }) ### end of observe event


  }
  shinyApp(ui, server)
} ### end of calculating GUI function
