# Load R packages
library(shiny)
library(shinythemes)
library(shinyFiles)
library(readxl)
library(ggplot2)
library(writexl)
library(here)

#Functions
HWE <- function(SNPtable, alpha = 0.05){
  # Observed genotype frequencies in the dataset
  HWTable <- table(SNPtable[, 2])
  
  # Observed allele frequencies in the dataset
  HWsum <- 2*(HWTable[[1]] + HWTable[[2]] + HWTable[[3]])
  p <- (2*HWTable[[1]] + HWTable[[2]])/HWsum
  q <- (2*HWTable[[3]] + HWTable[[2]])/HWsum
  
  # Theoretical probabilities of genotype frequencies
  HWexp <- c(p^2, 2*p*q, q^2)
  
  # Adding the expected genotype frequencies to the table
  HWTable <- rbind(HWTable, HWexp*sum(HWTable))
  rownames(HWTable) <- c("Observed", "Expected")
  
  # Printing the table with expected genotype counts
  print(HWTable) 
  
  # Finally, we can perform the chi squared test
  # We use just observed values for the testing. To determine the theoretical allele frequencies,
  # we use the 'p' argument of the chisq.test() function.
  pvalue <- chisq.test(HWTable[1,], p = HWexp)$p.value
  
  # Output of the function
  if(pvalue < alpha){
    pasting = paste0("The population is not at Hardy-Weinberg equilibrium (p-value = ", round(pvalue, 4), ")")
  }else{
    pasting = paste0("The population is at Hardy-Weinberg equilibrium (p-value = ", round(pvalue, 4), ")")
  }
  result <-  list(pasting,HWTable)
  
  return(result)
}

# Tests
SA_tests <- function(data){
  test_data <- SA_counts(data)
  return(list(
    fisher.test(test_data),
    chisq.test(test_data)
  ))
}

# Define UI
ui <- fluidPage(
  theme = shinytheme("darkly"),
  navbarPage(
    "Allelyzer",
    tabPanel("Input file and settings",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("select_ui_SNP_select"),
                 tags$h3("Tool selection:"),
                 checkboxGroupInput("checkboxes", "Select options:",
                                    choices = c("Hardy-Weinberg equilibrium", "CHi-Sq", "Heatmap"),
                                    selected = NULL),
                 textOutput("selected_text"),
                 tags$h4("Use selected tools"),
                 actionButton("tool_button", "Use Tools")
               ),
               mainPanel(
                 tags$h3("Input:"),
                 fileInput("file1", "Choose your file please", accept = ".xlsx"),
                 textOutput("HWE_out"),
                 div(
                   tableOutput("con_table"),
                   style = "margin-left: 100px"
                 )
               )
             )
    ), # tabPanel 1
    
    tabPanel("File overview", 
             sidebarLayout(
               sidebarPanel(
                 tags$h3("File viewer"),
                 uiOutput("select_ui2"),
                 uiOutput("select_ui"),
                 uiOutput("select_ui3")
               ),
               mainPanel(
                 tableOutput("data_table")
               )
             ) # sidebarLayout
    ) # tabPanel 2
  ) # navbarPage
) # fluidPage

# Define server function
server <- function(input, output) {
  # Read Excel file and create selectInput
  observeEvent(input$file1, {
    req(input$file1)
    # Read Excel file
    data <- read_xlsx(input$file1$datapath)
    
    # Extract unique values from "diagnosis" column
    unique_diagnosis <- unique(data$diagnosis)
    SNPs <- names(data)[-1]
    
    # Create selectInput with unique values
    output$select_ui <- renderUI({
      checkboxGroupInput("selected_diagnosis", "Select a diagnosis:", 
                         choices = c(unique_diagnosis), selected = NULL)
    })
    
    output$select_ui2 <- renderUI({
      selectInput("selected_SNP", "Select a SNP:", choices = c("ALL", SNPs))
    })
    
    output$select_ui3 <- renderUI({
      selected_SNP_for_search <- input$selected_SNP
      print(selected_SNP_for_search)
      unique_allels_of_SNP <- unique(data[, selected_SNP_for_search])
      unique_allels_of_SNP <- unique_allels_of_SNP[!is.na(unique_allels_of_SNP)]
      
      if (input$selected_SNP != "ALL") {
        selectInput("selected_allel", "Select a genotype:", 
                    choices = c("ALL", "A", "G", "T", "C",unique_allels_of_SNP))
      }
    })
    output$select_ui_SNP_select <- renderUI({
      selectInput("selected_SNP_for_analysis", "Select SNP to analyze:", choices = c(SNPs))
    })
  })    
  
  # Render data table based on selected diagnosis
  output$data_table <- renderTable({
    req(input$file1)
    req(input$selected_diagnosis)
    req(input$selected_SNP)
    req(input$selected_allel)
    
    # Read Excel file
    data <- read_xlsx(input$file1$datapath)
    
    # Filter data based on selected SNP
    if (input$selected_SNP != "ALL") {
      string_var <- as.character(input$selected_SNP)
      data <- data[, c("diagnosis", string_var)]
    }
    else{
      data <- data
    }
    
    # Filter data based on selected diagnosis
    if(input$selected_SNP != "ALL"){
      data <- data[data$diagnosis %in% input$selected_diagnosis & !is.na(data[,input$selected_SNP]), ]
    }
    else{
      data <- data[data$diagnosis %in% input$selected_diagnosis, ]
    }
    
    # Filter data based on selected allel
    if (input$selected_allel != "ALL") {
      data <- data[grepl(input$selected_allel, data[[input$selected_SNP]]) & !is.na(data[,input$selected_SNP]), ]
    }
    data
  })
  output$selected_text <- renderText({
    selected_options <- input$checkboxes
    if (length(selected_options) == 0) {
      "No options selected"
    } else {
      paste("Selected Tools:", paste(selected_options, collapse = ", "))
    }
  })
  #HWE function start
  output$HWE_out <- renderText({
    req(input$file1)
    data <- read_xlsx(input$file1$datapath)
    req(input$tool_button)
    if (input$checkboxes[1] == "Hardy-Weinberg equilibrium") {
      SNP <- cbind(data[, 1], data[, input$selected_SNP_for_analysis])
      
      # Cleaning the SNP table
      SNP <- SNP[SNP[, 2] %in% c("AA", "GG", "CC","TT", "AT", "AC", "AG", "CT", "CG", "GT", NA), ]
      
      # Factorizing the data
      SNP[, 1] <- as.factor(SNP[, 1])
      SNP[, 2] <- as.factor(SNP[, 2])
      list_from_HWE <- HWE(SNP)
      
      # Return the text from HWE function
      list_from_HWE[[1]]
    }
  })
  
  output$con_table <- renderTable({
    req(input$file1)
    data <- read_xlsx(input$file1$datapath)
    req(input$tool_button)
    if (input$checkboxes[1] == "Hardy-Weinberg equilibrium") {
      SNP <- cbind(data[, 1], data[, input$selected_SNP_for_analysis])
      
      # Cleaning the SNP table
      SNP <- SNP[SNP[, 2] %in% c("AA", "GG", "CC","TT", "AT", "AC", "AG", "CT", "CG", "GT", NA), ]
      
      # Factorizing the data
      SNP[, 1] <- as.factor(SNP[, 1])
      SNP[, 2] <- as.factor(SNP[, 2])
      list_from_HWE <- HWE(SNP)
      
      # Return the table from HWE function
      list_from_HWE[[2]]
    }
  }, rownames = TRUE)
  
}
# Run the application
shinyApp(ui = ui, server = server)