# FEBA.R based on Wetmore et al. 2025
# Estimate gene fitness values for methylamine conditions

# STEP 1: Load and Prepare Data
# --------------------------------------------------------------------------------------------

# Load gene count data: contains counts of sequencing reads mapped to genes.
genes <- read.csv("FEBA/MA/genes.csv", header = TRUE, sep = ",")

# Load gene count data for post-condition (countCond) and initial (countT0) experiments.
# These tables contain counts of sequencing reads for each gene across multiple experimental conditions.
countCond <- read.csv("FEBA/MA/countCond.csv", header = TRUE, row.names = 1, sep = ",") 
countT0 <- read.csv("FEBA/MA/countT0.csv", header = TRUE, row.names = 1, sep = ",")

# Load metadata describing experimental conditions and replicates.
expsUsed <- read.csv("FEBA/MA/expsUsed.csv", header = TRUE, row.names = 1, sep = ",")



# STEP 2: Define Functions
# ---------------------------------------------------------------------------------------
# - Develop functions to calculate fitness as a measure of gene importance under different conditions.
# - Normalize fitness values to account for variability across experiments or scaffolds.

# Function to calculate fitness: measures the relative growth advantage/disadvantage of mutants.
calculate_fitness <- function(countCond, countT0) {
  # Log2 ratio of post-condition counts to initial counts. 
  # A positive value indicates increased representation from JLW8 single to co-culture; negative indicates loss.
  log2(1 + countCond) - log2(1 + countT0)
}

# Function to normalize fitness values: centers fitness values for meaningful comparisons.
normalize_fitness <- function(fitness, genes, minToUse = 1) {
  # Normalize fitness globally if there are sufficient genes available.
  if (nrow(genes) < minToUse) {
    stop("Too few genes for normalization.")
  }
  fitness - median(fitness, na.rm = TRUE)  # Center fitness values around the median.
}

# Function to aggregate fitness results: summarizes data for each gene.
aggregate_fitness <- function(fitness, genes) {
  # Summarize fitness values across all experimental conditions.
  data.frame(
    locusId = genes$locusId,  # Gene identifiers.
    fit = rowMeans(fitness, na.rm = TRUE),  # Mean fitness per gene.
    se = apply(fitness, 1, function(x) sd(x, na.rm = TRUE)),  # Standard error for variability.
    n = rowSums(!is.na(fitness))  # Number of experiments with valid fitness values.
  )
}



# STEP 3: Execute Workflow
# ---------------------------------------------------------------------------------------
# - Identify genes of interest by filtering for sufficient data quality.
# - Calculate and normalize fitness values for each experimental condition.
# - Aggregate fitness results to determine gene-level significance.

# Define thresholds for filtering genes based on counts in the initial condition (JLW8 single cultures).
minTotalCounts <- 1   # Minimum total counts across all conditions.
minNonZeroConditions <- 1  # Minimum number of conditions with at least one count.

# Apply filtering criteria to retain only high-confidence genes.
genesUsed <- (rowSums(countT0) >= minTotalCounts) & 
  (rowSums(countT0 > 0) >= minNonZeroConditions)
filtered_genes <- genes[genes$locusId %in% rownames(countT0)[genesUsed], ]  
countT0 <- countT0[genesUsed, ]  
countCond <- countCond[genesUsed, ]  

# Calculate fitness for each experiment (column) using the fitness function.
fitness_results <- lapply(colnames(countT0), function(exp_name) {
  # Calculate fitness for the experiment.
  fitness <- calculate_fitness(countCond[, exp_name], countT0[, exp_name])
  # Normalize the fitness values.
  normalized_fitness <- normalize_fitness(fitness, filtered_genes)
  # Return results as a data frame.
  data.frame(locusId = filtered_genes$locusId, fit = normalized_fitness)
})
names(fitness_results) <- colnames(countT0)

# Combine and summarize fitness results across all experiments.
all_fitness_MA <- do.call(cbind, lapply(fitness_results, function(df) df$fit))  # Merge fitness results.
aggregated_results_MA <- aggregate_fitness(all_fitness_MA, filtered_genes)  # Summarize results.


# STEP 4: Add Gene Product Information
# ---------------------------------------------------------------------------------------
# - Annotate fitness results with gene product descriptions for biological interpretation.

# Load gene information data and rename columns for clarity.
gene_info <- read.csv("FEBA/MA/gene_information.csv", header = TRUE, sep = ",")

# Merge product information into the aggregated fitness results.
aggregated_results_MA <- merge(
  aggregated_results_MA,  
  gene_info[, c("locusId", "product")], 
  by = "locusId",  
  all.x = TRUE  
)



# STEP 5: Visualize Fitness Data
# ---------------------------------------------------------------------------------------
# - Identify and highlight genes with the most extreme fitness values (e.g., highly negative).

# Select the 20 genes with the most negative fitness values.
top_negative_fit_MA <- aggregated_results_MA %>%
  arrange(fit) %>%  
  slice_head(n = 20)  

# Create a bar plot to visualize these fitness values.
library(ggplot2)
ggplot(top_negative_fit_MA, aes(x = reorder(product, -fit), y = fit)) +  
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +  
  labs(
    title = "Most Negative Fitness Values (MA conditions)",  
    x = "JLW8 gene mutant product",  
    y = "Normalized gene fitness"  
  ) +
  theme_minimal() +
  coord_flip()  


# STEP 6: Check counts of individual experiments for a given gene
# ---------------------------------------------------------------------------------------
get_gene_counts <- function(gene_name, condition1, condition2, condition3, gene_count_cds, metadata) {
  if (!gene_name %in% rownames(gene_count_cds)) {
    stop("Gene name not found in gene_count_cds")
  }
  gene_counts <- gene_count_cds[gene_name, ]
  counts_df <- data.frame(Position = colnames(gene_count_cds), Count = as.numeric(gene_counts))
  metadata_with_counts <- merge(metadata, counts_df, by = "Position")
  filtered_data <- metadata_with_counts[
    metadata_with_counts$Condition %in% c(condition1, condition2, condition3), 
  ]
  result <- filtered_data[, c("Condition", "Replicate", "Count")]
  result <- result[order(result$Condition, result$Replicate), ]  
  return(result)
}

# Usage
condition1 <- "J_MeOH_La_P3"  # Replace with the first condition
condition2 <- "J_MeOH_La_P5" # Replace with the second condition
condition3 <- "J_MA_La_P3" # Replace with the third condition
gene_of_interest <- "MMOL_RS10940"  # Replace with the gene name you're interested in

# Call the function
result <- get_gene_counts(gene_of_interest, condition1, condition2, condition3, gene_count_cds, metadata)

# View the result
print(result)

