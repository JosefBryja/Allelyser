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
  result <-  list(pasting, HWTable)
  
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

#####################################
###### Test for single alleles ######
#####################################

# Get the allele count for each diagnosis
SA_counts <- function(data){
  data <- data[!is.na(data[,2]),]
  diagnoses <- levels(data[,1])
  names <- levels(as.factor(unlist(strsplit(as.character(data[,2]), ""))))
  
  # Initialize empty dataframe
  out <- data.frame(matrix(nrow = 0, ncol = length(names)))
  
  for(diag in diagnoses){
    subdata <- data[data[,1] == diag,]
    if(nrow(subdata) == 0) next  # Skip if no observations for this diagnosis
    single <- unlist(strsplit(as.character(subdata[,2]), ""))
    single <- single[!is.na(single)]
    Allele_table <- as.data.frame(table(single))
    
    # Append to output dataframe
    out <- rbind(out, Allele_table[,2])
  }
  
  # Set row names for diagnoses with observations
  rownames(out) <- diagnoses[sapply(diagnoses, function(x) any(data[,1] == x))]
  colnames(out) <- names
  return(out)
}

# Tests
SA_tests <- function(data){
  test_data <- SA_counts(data)
  return(list(
    fisher.test(test_data, workspace = 2e7),
    chisq.test(test_data)
  ))
}

odds_ratio <-function(data){
  counts <- as.matrix(SA_counts(data))
  rows <- rownames(counts)
  out <- as.data.frame(matrix(nrow = 0, ncol = 4))
  for(row in rows){
    test_table <- counts[row,]
    test_table <- rbind(test_table, colSums(counts[!(rownames(counts) == row),]))
    rownames(test_table) <- c(row, "Others")
    ftest <- fisher.test(test_table)
    #print(test_table)
    out <- rbind(out, c(ftest$estimate[[1]], ftest$conf.int[1], ftest$conf.int[2], ftest$p.value))
  }
  colnames(out) <- c("Odds Ratio", "Lower", "Upper", "p-value")
  rownames(out) <- rows
  return(out)
}

mosaicplots <-function(data, diagnosis){
  counts <- as.matrix(SA_counts(data))
  row <- diagnosis
  print(row)
  test_table <- counts[row,]
  test_table <- rbind(test_table, colSums(counts[!(rownames(counts) == row),]))
  rownames(test_table) <- c(row, "Others")
  mosaicplot(test_table, main = row, color = c("cornflowerblue", "darkred"), cex.axis = 1.5)
}


# Define UI
ui <- fluidPage(
  theme = shinytheme("darkly"),
  navbarPage(
    "Allelyzer",
    tabPanel("Input file and settings",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("main_sidebar1"),
                 uiOutput("main_sidebar2"),
                 uiOutput("main_sidebar3"),
                 uiOutput("main_sidebar4"),
                 uiOutput("main_sidebar5"),
                 uiOutput("main_sidebar6"),
                 uiOutput("main_sidebar7"),
               ),
               mainPanel(
                 tags$h3("Input:"),
                 fileInput("file1", "Choose your file please", accept = ".xlsx"),
                 div(
                   textOutput("selected_SNP_text"),
                   style = "font-size: 20px; color: darkred;"
                 ),
                 div(
                 textOutput("filter_out_text"),
                 style = "font-size: 15px;"
                 ),
                 verbatimTextOutput("filtered_out"),
                 ##########---HWE---###########
                 div(
                   textOutput("HWE_head_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 textOutput("HWE_out"),
                 div(
                   tableOutput("con_table"),
                   style = "margin-left: 10px"
                 ),
                 div(
                   textOutput("HWE_tail_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 uiOutput("HW_download_button"),
                 
                 ##########---Chi-sq---###########
                 div(
                   textOutput("Chi_sq_head_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 tableOutput("Chi_sq_table"),
                 div(
                   textOutput("Chi_sq_tail_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 uiOutput("ChiSq_download_button"),
                 
                 ##########---Heatmap---###########
                 div(
                   textOutput("Heatmap_head_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 uiOutput("SNP_heat_UI"),
                 div(
                   textOutput("Heatmap_tail_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 uiOutput("heatMap_download_button"),
                 
                 ##########---SA-tests---###########
                 div(
                   textOutput("SA_tests_head_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 verbatimTextOutput("fisher_out"),
                 verbatimTextOutput("chi_out"),
                 div(
                   textOutput("SA_tests_tail_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 
                 ##########---Odds_ratio---###########
                 div(
                   textOutput("Odds_ratio_head_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 tableOutput("odds_ratio_table"),
                 div(
                   textOutput("Odds_ratio_tail_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 uiOutput("odds_ratio_download_button"),
                 
                 ##########---Mosaicplot---###########
                 div(
                   textOutput("Mosaicplot_head_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 plotOutput("mosaic_table"),
                 div(
                   textOutput("Mosaicplot_tail_text"),
                   style = "font-size: 15px;color: green"
                 ),
                 uiOutput("mosaic_download_button")
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
    
    #############_Main_sidebar###########
    output$main_sidebar1 <- renderUI({
      uiOutput("select_ui_SNP_select")
    })
    output$main_sidebar2 <- renderUI({
      tags$h3("Tool selection:")
    })
    output$main_sidebar3 <- renderUI({
      checkboxGroupInput("checkboxes", "Select options:",
                         choices = c("Hardy-Weinberg equilibrium", "CHi-Sq", "Heatmap", "SA_tests", 
                                     "odds_ratio", "mosaicplots"),
                         selected = NULL)
    })
    output$main_sidebar4 <- renderUI({
      textOutput("selected_text")
    })
    output$main_sidebar5 <- renderUI({
      tags$h4("Use selected tools")
    })
    output$main_sidebar6 <- renderUI({
      actionButton("tool_button", "Use Tools")
    })
    output$main_sidebar7 <- renderUI({
      uiOutput("mosaic_select")
    })
    #####################################
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
  
  observeEvent(input$checkboxes,{
    req(input$file1)
    data <- read_xlsx(input$file1$datapath)
    unique_diagnosis <- unique(data$diagnosis)
    if ("mosaicplots" %in% input$checkboxes){
      output$mosaic_select <- renderUI({
        selectInput("mosaic_select_for_analysis", "Select diagnosis for mosaicplots", choices = c(unique_diagnosis))
      })
    }
  })
  
  
  
  ########################################################
  observeEvent(input$tool_button, {
    req(input$file1)
    req(input$tool_button)
    data <- read_xlsx(input$file1$datapath)
    text_data <- input$selected_SNP_for_analysis
    unique_diagnosis <- unique(data$diagnosis)
    
    SNP <- cbind(data[, 1], data[, input$selected_SNP_for_analysis])
    SNP_not_in_filter <- SNP[!(SNP[, 2] %in% 
                                 c("AA", "GG", "CC","TT", "AT", "AC", "AG", "CT", "CG", "GT", NA)), ]
    
    SNP <- SNP[SNP[, 2] %in% c("AA", "GG", "CC","TT", "AT", "AC", "AG", "CT", "CG", "GT", NA), ]
    
    if (!is.null(input$checkboxes)){
      output$selected_SNP_text <- renderText({
        req(input$tool_button)
        paste("Selected SNP for analysis:", text_data)
       })
      
      output$filter_out_text <- renderText({
        req(input$tool_button)
        paste("Values that were filtered out of the dataset for  ", text_data,"  SNP:")
      })
      
      output$filtered_out <- renderPrint({
        # Convert the result to a character vector
        result_text <- capture.output(print(SNP_not_in_filter))
        
        # Print the result
        cat(result_text, sep = "\n")
      })
      
    }else{
      output$selected_SNP_text <- renderText({
        return(NULL)
      })
      output$filtered_out <- renderText({
        return(NULL)
      })
    }
    
    
    if ("Hardy-Weinberg equilibrium" %in% input$checkboxes) {
      output$HWE_out <- renderText({
        # Factorizing the data
        SNP[, 1] <- as.factor(SNP[, 1])
        SNP[, 2] <- as.factor(SNP[, 2])
        list_from_HWE <- HWE(SNP)
        
        # Return the text from HWE function
        list_from_HWE[[1]]
      })
      output$HWE_head_text <- renderText({
        req(input$tool_button)
        paste("Test if the population is in the Hardy-Weinberg equilibrium for SNP", text_data,":")
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
      
      output$HWE_tail_text <- renderText({
        req(input$tool_button)
        paste("If p-value is greater than 0.05, the population is in the Hardy-Weinberg equilibrium")
      })
      
      output$downloadHW <- downloadHandler(
        filename = function() { 
          paste("HWE_Table",text_data,Sys.Date(),".xlsx" , sep = "_")
        },
        content = function(file) {
          SNP[, 1] <- as.factor(SNP[, 1])
          SNP[, 2] <- as.factor(SNP[, 2])
          list_from_HWE <- HWE(SNP)
          write_xlsx(as.data.frame(list_from_HWE[[2]]),path = file,)
        })
    }else {
      output$HWE_head_text <- renderText({
        return(NULL)
      })
      output$HWE_out <- renderText({
        return(NULL)
      })
      output$con_table <- renderTable({
        return(NULL)
      })
      output$HW_download_button <- renderUI({
        return(NULL)
      })
      output$HWE_tail_text <- renderText({
        return(NULL)
      })
    }
    ########################################################
    
    if ("CHi-Sq" %in% input$checkboxes) {
      
      output$Chi_sq_head_text <- renderText({
        req(input$tool_button)
        paste("Test if the population is in the Hardy-Weinberg equilibrium for each diagnosis in SNP", text_data,":")
      })
      
      output$Chi_sq_table <- renderTable({
        output$ChiSq_download_button <- renderUI({
          downloadButton('downloadChiSq', 'Download CHi-Sq table')
        })
        CHi_sq <- ChiSq(SNP)
      }, rownames = TRUE)
      
      output$Chi_sq_tail_text <- renderText({
        req(input$tool_button)
        paste("The p-value from Chi-Square test is presented for each diagnosis in selected SNP. 
              If this value is greater than 0.05, the population is in the Hardy-Wwinberg equilibrium in selected SNP and diagnosis.")
      })
      
      output$downloadChiSq <- downloadHandler(
        filename = function() { 
          paste("Chi_sq_Table",text_data,Sys.Date(),".xlsx" , sep = "_")
        },
        content = function(file) {
          xlsx_out <- write_xlsx(as.data.frame(ChiSq(SNP)),path = file)
        }
      )
    }else {
      output$Chi_sq_head_text <- renderText({
        return(NULL)
      })
      output$Chi_sq_table <- renderTable({
        return(NULL)
      })
      output$ChiSq_download_button <- renderUI({
        return(NULL)
      })
      output$Chi_sq_tail_text <- renderText({
        return(NULL)
      })
    }
    
    ########################################################
    if ("Heatmap" %in% input$checkboxes) {
      
      output$Heatmap_head_text <- renderText({
        req(input$tool_button)
        paste("Check the distribution of the data with heatmap representation")
      })
      ### - No gap plot
      output$SNP_heat <- renderPlot({
        heat <- SNPHeatmap(SNP)
      })
      
      output$SNP_heat_UI <- renderUI({
        plotOutput("SNP_heat")
      })
      #################################
      
      output$heatMap_download_button <- renderUI({
        downloadButton('downloadPlot', 'Download Plot')
      })
      
      output$Heatmap_tail_text <- renderText({
        req(input$tool_button)
        paste("You can check how are the genotypes distributed in different diagnoses in selected SNP.")
      })
      
      output$downloadPlot <- downloadHandler(
        filename = function() { 
          paste("Heatmap",text_data,Sys.Date(),".png" , sep = "_")
        },
        content = function(file) {
          png(file)
          print(SNPHeatmap(SNP))
          dev.off()
        })
    }else {
      output$Heatmap_head_text <- renderText({
        return(NULL)
      })
      output$SNP_heat_UI <- renderTable({
        return(NULL)
      })
      output$heatMap_download_button <- renderUI({
        return(NULL)
      })
      output$Heatmap_tail_text <- renderText({
        return(NULL)
      })
    }
    ######################################################
    if("SA_tests" %in% input$checkboxes){
      
      output$SA_tests_head_text <- renderText({
        req(input$tool_button)
        paste("Chi-squared test and Fisher's exact test for SPP", text_data,":")
      })
      
      SNP[, 1] <- as.factor(SNP[, 1])
      SNP[, 2] <- as.factor(SNP[, 2])
      output$fisher_out <- renderPrint({
        # Call the SA_tests function to get the Fisher's test output
        fisher_result <- SA_tests(SNP)
        
        # Convert the result to a character vector
        result_text <- capture.output(print(fisher_result[[1]]))
        
        # Print the result
        cat(result_text, sep = "\n")
      })
      output$chi_out <- renderPrint({
        # Call the SA_tests function to get the Fisher's test output
        fisher_result <- SA_tests(SNP)
        
        # Convert the result to a character vector
        result_text <- capture.output(print(fisher_result[[2]]))
        
        # Print the result
        cat(result_text, sep = "\n")
      })
      
      output$SA_tests_tail_text <- renderText({
        req(input$tool_button)
        paste("These tests are showing you if there is any relation ship between the alleles and diagnoses. 
              If the p-value is lower than 0.05, you may conclude that there is some relationship. 
              If you have data with very low number of observations, it is better to use Fisher's exact test. 
              If you have High number of observations, you can use Chi-Squared test.")
      })
      
    }else{
      output$chi_out <- renderTable({NULL})
      output$fisher_out <- renderTable({NULL})
      output$SA_tests_head_text <- renderText({
        return(NULL)
      })
      output$SA_tests_tail_text <- renderText({
        return(NULL)
      })
    }
    ######################################################
    if ("odds_ratio" %in% input$checkboxes) {
      
      output$Odds_ratio_head_text <- renderText({
        req(input$tool_button)
        paste("Check what are the odds of beeing diagnosed with particular disease based on allele")
      })
      
      
      SNP[, 1] <- as.factor(SNP[, 1])
      SNP[, 2] <- as.factor(SNP[, 2])
      output$odds_ratio_table <- renderTable({
        output$odds_ratio_download_button <- renderUI({
          downloadButton('download_odds_ratio', 'Download odds_ratio table')
        })
        odds_ratio_out <- odds_ratio(SNP)
        odds_ratio_out
      }, rownames = TRUE)
      
      output$Odds_ratio_tail_text <- renderText({
        req(input$tool_button)
        paste("Check the p-value or confidence interval to see if this result is statistically significant.")
      })
      
      output$download_odds_ratio <- downloadHandler(
        filename = function() { 
          paste("Odds_ratio",text_data,Sys.Date(),".xlsx" , sep = "_")
        }, 
        content = function(file) {
          odds_ratio_df <- odds_ratio(SNP) # Get odds ratio result as data frame
          odds_ratio_df <- cbind(rownames(odds_ratio_df), odds_ratio_df)  # Add row names as first column
          colnames(odds_ratio_df)[1] <- ""  # Rename the first column to "RowNames"
          write_xlsx(odds_ratio_df, path = file)  # Write data frame to Excel
        }
      )
    }else {
      output$Odds_ratio_head_text <- renderText({
        return(NULL)
      })
      output$odds_ratio_table <- renderTable({
        return(NULL)
      })
      output$Odds_ratio_tail_text <- renderText({
        return(NULL)
      })
      output$odds_ratio_download_button <- renderUI({
        return(NULL)
      })
    }
    ####################################################
    if ("mosaicplots" %in% input$checkboxes) {
      output$Mosaicplot_head_text <- renderText({
        req(input$tool_button)
        paste("Check the distribution of alleles for each disease compared to the rest of the dataset in selected SNP.",
              input$mosaic_select_for_analysis,"__",":")
      })
      
      SNP[, 1] <- as.factor(SNP[, 1])
      SNP[, 2] <- as.factor(SNP[, 2])
      output$mosaic_table <- renderPlot({
        heat <- mosaicplots(SNP,input$mosaic_select_for_analysis)
        output$mosaic_download_button <- renderUI({
          downloadButton('downloadPlot_mosaic', 'Download Plot')
        })
      })
      
      output$Mosaicplot_tail_text <- renderText({
        req(input$tool_button)
        paste("You can compare this graphical representation with the odds ratio results.",
              input$mosaic_select_for_analysis,"__",":")
      })
      
      output$downloadPlot_mosaic <- downloadHandler(
        filename = function() { 
          paste("Mosiac_plot",text_data,input$mosaic_select_for_analysis
                ,Sys.Date(),".png" , sep = "_")
        },
        content = function(file) {
          if ("mosaicplots" %in% input$checkboxes) {
            png(file)
            print(mosaicplots(SNP,input$mosaic_select_for_analysis))
            dev.off()
          }
        })
    }else {
      output$Mosaicplot_head_text <- renderText({
        return(NULL)
      })
      output$mosaic_table <- renderTable({
        return(NULL)
      })
      output$Mosaicplot_tail_text <- renderText({
        return(NULL)
      })
      output$mosaic_download_button <- renderUI({
        return(NULL)
      })
    }
  })
}
# Run the application
shinyApp(ui = ui, server = server)