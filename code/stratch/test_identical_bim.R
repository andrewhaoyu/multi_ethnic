
# Load required library
library(data.table)
# Define the list of chromosomes
chromosomes <- paste0("chr", 1:22)
# Function to compare files
compare_files <- function(chr) {
  ancestries <- c("EUR", "AFR", "AMR", "EAS", "SAS")
  
  # List to store file contents
  file_contents <- vector("list", length = length(ancestries))
  
  # Read data for each ancestry and chromosome
  for (i in seq_along(ancestries)) {
    # File path for the current ancestry and chromosome
    #file_path <- paste0("/data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh37/", ancestries[i], "/", chr, ".bim")
    file_path <- paste0("/data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh38/", ancestries[i], "/", chr, ".bim")
    # Read content of the file using fread
    file_contents[[i]] <- fread(file_path, header = FALSE)
  }
  
  # Check if all files are identical
  identical_status <- sapply(2:length(file_contents), function(i) identical(file_contents[[i]], file_contents[[1]]))
  
  identical_status
}


# Loop through chromosomes
for (chr in chromosomes) {
  identical_status <- compare_files(chr)
  
  # Output results
  cat("Chromosome", chr, ":\n")
  for (i in seq_along(identical_status)) {
    if (identical_status[[i]]) {
      cat("Ancestry", c("EUR", "AFR", "AMR", "EAS", "SAS")[i], "files are identical.\n")
    } else {
      cat("Ancestry", c("EUR", "AFR", "AMR", "EAS", "SAS")[i], "files are not identical.\n")
    }
  }
  cat("\n")
}



#Define the list of chromosomes
chromosomes <- paste0("chr", 1:22)

# Function to check if second column contains "."
check_dot <- function(chr) {
  # List of ancestries
  ancestries <- c("EUR", "AFR", "EAS", "AMR", "SAS")
  
  for (ancestry in ancestries) {
    # File path for the current ancestry and chromosome
    file_path <- file.path("/data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh37", ancestry, paste0(chr, ".bim"))
    
    # Read content of the file using fread
    file_content <- fread(file_path, header = FALSE)
    # Check if any element in the second column contains "."
    idx <- which(file_content$V2==".")
    if(length(idx)>0){
      cat("Chromosome", chr, "for", ancestry, "contains '.' in the second column.\n")
    }else {
      cat("Chromosome", chr, "for", ancestry, "does not contain '.' in the second column.\n")
    }
   
   
  }
}

# Loop through chromosomes
for (chr in chromosomes) {
  check_dot(chr)
}
