# Load R packages
library(shiny)
library(shinythemes)
library(shinyFiles)
library(readxl)
library(ggplot2)
library(writexl)
library(here)
library(epitools)
library(gridGraphics)

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

ChiSq <- function(data, significant_only = FALSE, alpha = 0.05){
  
  # Creating and printing contingency table
  conTable <- table(data)
  
  # For each diagnosis, we want to test if there is difference in the genotype frequencies
  # Empty data frame
  pData <- data.frame(matrix(nrow = 0, ncol = 2))
  
  # Chi-square for each diagnosis
  for(x in 1:length(rownames(conTable))){
    
    # Single diagnosis
    t <- conTable[x, ]
    
    # Observed allele frequencies in the dataset
    HWsum <- 2*(t[[1]] + t[[2]] + t[[3]])
    p <- (2*t[[1]] + t[[2]])/HWsum
    q <- (2*t[[3]] + t[[2]])/HWsum
    
    # Theoretical probabilities of genotype frequencies
    HWexp <- c(p^2, 2*p*q, q^2)
    
    # Chi-square test
    if(sum(t) > 0){
      pData <- rbind(pData, c(rownames(conTable)[x], round(chisq.test(t, p = HWexp)$p.value, 6)))
    }else{
      pData <- rbind(pData, c(rownames(conTable)[x], NA))
    }
  }
  
  # Renaming colums
  colnames(pData) <- c("Diagnosis", "p-value")
  
  # Changing p-values to numeric
  pData[,2] <- as.numeric(pData[,2])
  
  # Return the output table
  if(significant_only == TRUE){
    return(na.omit(pData[pData[,2] < alpha,]))
  }else{
    return(pData)
  }
}

# SNPHeatmap <- function(data, scaled = TRUE){
#   
#   # Create data for visualization
#   name <- colnames(data)[2]
#   print(name)
#   t <- table(data)
#   t <- t[rowSums(t) > 0,]
#   print(t)
#   # Plot the data
#   if(scaled == TRUE){
#     return(heatmap(t, scale = "row", main = name))
#   }else if(scaled == FALSE){
#     return(heatmap(t, scale = "none"))
#   }
# }

SNPHeatmap <- function(data){
  name <- colnames(data)[2]
  # Create a table of counts
  count_table <- table(data)
  
  # Convert table to a data frame
  df <- as.data.frame(as.table(count_table))
  names(df) <- c("diagnosis", "SNP", "count")
  
  # Create the heatmap using ggplot2
  heatmap_plot <- ggplot(df, aes(x = SNP, y = diagnosis, fill = count, label = count)) +
    geom_tile(color = "white") +
    geom_text(size = 5, color = "black") +  # Add text labels
    scale_fill_gradient(low = "cornflowerblue", high = "darkred") +  # Color gradient
    labs(title = "SNP Heatmap", x = name, y = "diagnosis") +  # Axis labels and title
    theme_minimal() +  # Minimal theme
    
    # Rotate x-axis labels for better readability
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot
  print(heatmap_plot)
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
                 ),
                 uiOutput("HW_download_button"),
                 div(
                   tableOutput("Chi_sq_table"),
                   style = "margin-left: 100px"
                   
                 ),
                 uiOutput("ChiSq_download_button"),
                 
                 plotOutput("SNP_heat"),
                 uiOutput("heatMap_download_button")
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
                         choices = c(unique_diagnosis), selected = c(unique_diagnosis))
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
  
  
  
  
  ########################################################
  observeEvent(input$tool_button, {
    req(input$file1)
    req(input$tool_button)
    data <- read_xlsx(input$file1$datapath)
    SNP <- cbind(data[, 1], data[, input$selected_SNP_for_analysis])
    SNP <- SNP[SNP[, 2] %in% c("AA", "GG", "CC","TT", "AT", "AC", "AG", "CT", "CG", "GT", NA), ]
    
    if ("Hardy-Weinberg equilibrium" %in% input$checkboxes) {
      output$HWE_out <- renderText({
        # Factorizing the data
        SNP[, 1] <- as.factor(SNP[, 1])
        SNP[, 2] <- as.factor(SNP[, 2])
        list_from_HWE <- HWE(SNP)
        
        # Return the text from HWE function
        list_from_HWE[[1]]
      })
      
      output$con_table <- renderTable({
        output$HW_download_button <- renderUI({
          downloadButton('downloadHW', 'Download HW table')
        })
        # Factorizing the data
        SNP[, 1] <- as.factor(SNP[, 1])
        SNP[, 2] <- as.factor(SNP[, 2])
        list_from_HWE <- HWE(SNP)
        
        # Return the table from HWE function
        list_from_HWE[[2]]
      }, rownames = TRUE)
      
      output$downloadHW <- downloadHandler(
        filename = function() { "HW_Table.xlsx" },
        content = function(file) {
          SNP[, 1] <- as.factor(SNP[, 1])
          SNP[, 2] <- as.factor(SNP[, 2])
          list_from_HWE <- HWE(SNP)
          write_xlsx(as.data.frame(list_from_HWE[[2]]),path = file,)
        })
    }else {
      output$HWE_out <- renderText({
        return(NULL)
      })
      output$con_table <- renderTable({
        return(NULL)
      })
      output$HW_download_button <- renderUI({
        return(NULL)
      })
    }
    ########################################################
    
    if ("CHi-Sq" %in% input$checkboxes) {
      output$Chi_sq_table <- renderTable({
        output$ChiSq_download_button <- renderUI({
          downloadButton('downloadChiSq', 'Download CHi-Sq table')
        })
        CHi_sq <- ChiSq(SNP)
      }, rownames = TRUE)
      
      output$downloadChiSq <- downloadHandler(
        filename = function() "CHi-Sq_Table.xlsx" ,
        content = function(file) {
          xlsx_out <- write_xlsx(as.data.frame(ChiSq(SNP)),path = file)
        }
      )
    }else {
      output$Chi_sq_table <- renderTable({
        return(NULL)
      })
      output$ChiSq_download_button <- renderUI({
        return(NULL)
      })
    }
    
    ########################################################
    if ("Heatmap" %in% input$checkboxes) {
      output$SNP_heat <- renderPlot({
        heat <- SNPHeatmap(SNP)
        output$heatMap_download_button <- renderUI({
          downloadButton('downloadPlot', 'Download Plot')
        })
      })
      
      output$downloadPlot <- downloadHandler(
        filename = function() { "test.png" },
        content = function(file) {
          if ("Heatmap" %in% input$checkboxes) {
            png(file)
            print(SNPHeatmap(SNP))
            dev.off()
          }
        })
    }else {
      output$SNP_heat <- renderTable({
        return(NULL)
        })
      output$heatMap_download_button <- renderUI({
        return(NULL)
      })
    }
  })
}
# Run the application
shinyApp(ui = ui, server = server)