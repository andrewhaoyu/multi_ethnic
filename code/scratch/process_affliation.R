# Sample R code to process a list of authors and their affiliations

# Assuming 'authors_raw.txt' contains the raw list of authors and their affiliations in a structured format
# Load necessary library
library(dplyr)
library(stringr)

# Read in the data as a single string
authors_raw <- readLines("./data/author_raw.txt")  # Replace with your actual file path

# Split the raw string by each author entry based on the pattern for affiliations
authors_split <- str_split(authors_raw, ",\\s*(?=\\w+\\s+\\w+\\d+\\*)")[[1]]

# Initialize a dataframe to store parsed data
authors_data <- data.frame(FirstName = character(),
                           MiddleName = character(),
                           LastName = character(),
                           Affiliation = character(),
                           stringsAsFactors = FALSE)

# Loop over each entry in authors_split
for (entry in authors_split) {
  # Remove the asterisk (*) and extra spaces
  entry <- str_remove(entry, "\\*")
  
  # Extract name part and affiliation number
  name_affil <- str_match(entry, "^(.+?)(\\d+)$")
  
  # Split the name part into first, middle, and last names
  name_part <- str_trim(name_affil[2])
  affiliation_number <- as.integer(name_affil[3])
  
  # Split names into components (first, middle, last)
  name_parts <- str_split(name_part, "\\s+")[[1]]
  first_name <- name_parts[1]
  last_name <- name_parts[length(name_parts)]
  middle_name <- ifelse(length(name_parts) > 2, name_parts[2], "")
  
  # Select the first affiliation based on the number
  affiliation <- paste0("Affiliation", affiliation_number)  # Placeholder; use a lookup table if you have affiliations mapped
  
  # Add to the dataframe
  authors_data <- authors_data %>%
    add_row(FirstName = first_name,
            MiddleName = middle_name,
            LastName = last_name,
            Affiliation = affiliation)
}

# Write the results to a TSV file
write.table(authors_data, "authors_affiliations_parsed.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


