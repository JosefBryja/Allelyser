# List of required libraries
required_libraries <- c("shiny", "shinythemes", "shinyFiles", "readxl", 
                        "ggplot2", "writexl", "here", "epitools", "gridGraphics")

# Function to check and install missing libraries
check_and_install_libraries <- function(required_libraries) {
  missing_libraries <- required_libraries[!(required_libraries %in% installed.packages())]
  
  if (length(missing_libraries) > 0) {
    message("Installing missing libraries...")
    install.packages(missing_libraries, dependencies = TRUE)
    message("Installation complete.")
  } else {
    message("All required libraries are already installed.")
  }
}

# Call the function
check_and_install_libraries(required_libraries)
