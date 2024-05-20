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