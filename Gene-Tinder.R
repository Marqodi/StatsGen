# ==============================================================================
# FILE: Gene_Tinder.R
# AUTHOR: Marqus Odisho
# ==============================================================================

# --- 1. LOAD LIBRARIES ---
# Check for and install missing packages to ensure portability.
# 'future' & 'future.apply': Essential for parallelizing the independent simulation runs.
# 'dplyr' & 'readr': Standard tools for fast data manipulation and file I/O.
required_packages <- c("future", "future.apply", "dplyr", "readr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(future)       
library(future.apply) 
library(dplyr)        
library(readr)        

# ==============================================================================
# PRIMARY FUNCTION: Gene_Tinder
# ==============================================================================

#' GENE-TINDER: Forward Genetic Simulation of Assortative Mating
#' GENEtic Traits INfluencing Divergence & Ecotype Reproduction
#'
#' This function runs a multi-generation assortative mating simulation with overlapping 
#' generations and fitness selection based on phenotype. Individuals can survive and reproduce for 2 generations.
#' 
#' Individuals are sorted into two distinct populations that can interbreed to produce hybrids.
#' Simulates the eco-evolutionary dynamics of two diverging populations (Species A & B) 
#' in a 3D environment. The model integrates:
#' 1. Habitat Isolation (Spatial segregation)
#' 2. Assortative Mating (Sexual selection based on phenotype/genotype)
#' 3. Natural Selection (Fitness penalties)
#' 4. Density-Dependent Growth (Discrete Logistic Regulation)
#'
#' @param experiment_name String. The name of the output directory where CSV files 
#'   will be saved (if running in batch mode). Default is "GeneTinder_Logistic".
#' @param num_runs Integer. The number of independent simulations to run. 
#'   If 1, returns the population object directly (Interactive Mode). 
#'   If >1, runs in parallel batch mode and saves CSVs to disk.
#' @param parallel Boolean. If TRUE and num_runs > 1, utilizes available CPU cores to run simulation in parallel.
#' 
#' @param weight_dist Numeric. The weight of physical Euclidean distance in mate choice.
#'   High values simulate low dispersal or preference for local mates.
#' @param weight_q Numeric. The weight of genetic similarity (q-score) in mate choice.
#' @param weight_p Numeric. The weight of phenotypic similarity in mate choice.
#' @param k_dist Numeric. Sensitivity coefficient for distance decay.
#' @param k_q Numeric. Sensitivity coefficient for genetic difference decay.
#' @param k_p Numeric. Sensitivity coefficient for phenotype difference decay.
#' 
#' @param min_fitness_scalar Numeric (0.0 - 1.0). Post-zygotic isolation strength.
#'   The relative survival probability of the least fit phenotype.
#' @param species_A_ratio Numeric (0.0 - 1.0). The initial proportion of Species A.
#'   Species A is homozygous for allele 0; Species B is homozygous for allele 2.
#' 
#' @param fine_scale_assortment Boolean. If TRUE, applies a non-linear boost to 
#'   mating probabilities for immediate neighbors (micro-allopatry).
#' @param fine_scale_multiplier Numeric. The magnitude of the neighbor boost.
#' 
#' @param pheno_loci_indices Vector of Integers. Which genome columns drive the phenotype?
#' @param pheno_heritability Numeric (0.0 - 1.0). The proportion of phenotypic 
#'   variance ($V_P$) attributable to additive genetic variance ($V_A$).
#' @param weight_random Numeric. Noise floor to prevent zero-probability mating matrices.
#' 
#' @param num_generations Integer. Duration of the simulation.
#' @param initial_pop Integer. Starting population size ($N_0$).
#' @param num_loci Integer. Genome length per individual.
#' @param max_population_size Integer. Environmental Carrying Capacity ($K$).
#' @param intrinsic_growth_rate Numeric. The intrinsic rate of increase ($r$) for 
#'   the discrete logistic equation. 
#' 
#' @param x_sd Numeric. Spatial variance for X-coordinate initialization.
#' @param y_sd Numeric. Spatial variance for Y-coordinate initialization.
#' @param z_sd Numeric. Spatial variance for Z-coordinate initialization.
#'
#' @return If num_runs = 1, returns a matrix of the final population. 
#'   If num_runs > 1, returns the file path to the aggregated results CSV.
#' @export

Gene_Tinder <- function(
    # --- Execution Parameters ---
  experiment_name = "GeneTinder_Logistic", 
  num_runs = 1,                        
  parallel = TRUE,                     
  
  # --- Mating Preference Parameters ---
  weight_dist,         
  weight_q,            
  weight_p,            
  k_dist,              
  k_q,                 
  k_p,                 
  
  # --- Selection Parameters ---
  min_fitness_scalar,  
  species_A_ratio,     
  
  fine_scale_assortment = FALSE, 
  fine_scale_multiplier = 10,    
  
  # --- Genetics & Environment ---
  pheno_loci_indices = 1:10,     
  pheno_heritability = 0.3,      
  weight_random = 0.1,           
  
  # --- Population Dynamics ---
  num_generations = 5,
  initial_pop = 2000,
  num_loci = 10000,
  max_population_size = 30000,   
  intrinsic_growth_rate = 0.5,   
  
  x_sd = 50, y_sd = 50, z_sd = 50 
) {
  
  # ============================================================================
  # INTERNAL WORKER: The Biological Engine
  # ============================================================================
  # This function runs a SINGLE simulation entirely within its own scope.
  # This encapsulation is crucial for 'future_lapply' to distribute the work
  # to different CPU cores without sharing memory states.
  simulation_worker <- function(run_id, save_dir = NULL) {
    
    # 1. SETUP OUTPUT
    if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    
    # 2. INITIALIZE POPULATION
    # Biological Context: We start with two divergent lineages (Species A & B).
    # This mimics secondary contact zones where previously isolated populations meet.
    
    n_species_A <- round(species_A_ratio * initial_pop)
    n_species_B <- initial_pop - n_species_A
    
    # --- Create Genotypes ---
    # We represent the genome as a matrix of alleles (0, 1, or 2).
    # 0 = Homozygous Species A allele
    # 2 = Homozygous Species B allele
    geno_mat <- matrix(0, nrow = initial_pop, ncol = num_loci)
    
    # Fill Species B alleles (if any exist in the initial ratio)
    if (n_species_B > 0) {
      geno_mat[(n_species_A + 1):initial_pop, ] <- 2
    }
    
    # --- Calculate Traits ---
    # Q-score (Hybrid Index): A continuous metric from 0.0 (Pure A) to 1.0 (Pure B).
    # Useful for tracking introgression and hybridization.
    q_scores <- rowSums(geno_mat) / (num_loci * 2)
    
    # Phenotype Calculation (Quantitative Genetics):
    # P = G + E
    # Phenotype = (Genotype * Heritability) + (Environment * (1 - Heritability))
    # We only sum specific loci ('pheno_loci_indices') to simulate specific trait architecture.
    p_gen <- rowSums(geno_mat[, pheno_loci_indices, drop = FALSE]) / (2 * length(pheno_loci_indices))
    pheno_scores <- (pheno_heritability * p_gen) + ((1 - pheno_heritability) * runif(initial_pop))
    
    starting_genotype <- cbind(geno_mat, q_score = q_scores, phenotype_score = pheno_scores)
    
    # --- Spatial Distribution (Genotype-Environment Correlation) ---
    # Biological Context: Divergent species often prefer different micro-habitats.
    # Species A prefers location 0; Species B prefers location 100.
    # Hybrids naturally fall into intermediate spatial niches.
    x_mean <- q_scores * 100; y_mean <- (1 - q_scores) * 100; z_mean <- (1 - q_scores) * 100
    
    # Place individuals stochastically around their preferred niche
    x_c <- round(pmax(0, rnorm(initial_pop, x_mean, x_sd)))
    y_c <- round(pmax(0, rnorm(initial_pop, y_mean, y_sd)))
    z_c <- round(pmax(0, rnorm(initial_pop, z_mean, z_sd)))
    
    # Primary Dataframe: Holds Genetics, Traits, Location, and Demographics
    pop_data <- cbind(starting_genotype, 
                      x_coordinate = x_c, y_coordinate = y_c, z_coordinate = z_c, 
                      age = rep(1, initial_pop),
                      id = 1:initial_pop,
                      maternal_id = rep(NA, initial_pop),
                      paternal_id = rep(NA, initial_pop))
    
    current_gen <- 1
    
    # ==========================================================================
    # MAIN GENERATION LOOP
    # ==========================================================================
    while (current_gen <= num_generations) {
      
      # --- A. FITNESS SELECTION (Viability Selection) ---
      # Biological Context: Individuals with maladaptive phenotypes are removed
      # from the gene pool before they can reproduce.
      n_before <- nrow(pop_data)
      if (n_before > 0) {
        # Fitness Function: Linear gradient.
        # If min_fitness_scalar = 0.7, Phenotype 0 has 70% survival probability vs Phenotype 1.
        fit_scores <- (1.0 - min_fitness_scalar) * pop_data[, "phenotype_score"] + min_fitness_scalar
        
        # Stochastic Survival Roll
        survivors <- which(runif(n_before) < fit_scores)
        pop_data <- pop_data[survivors, , drop = FALSE]
      }
      
      n_ind <- nrow(pop_data)
      if (n_ind < 2) break # Extinction threshold
      
      # ========================================================================
      # B. POPULATION GROWTH: DISCRETE LOGISTIC EQUATION
      # ========================================================================
      # Biological Context: Populations cannot grow indefinitely. They are regulated
      # by density-dependent factors (food, space, predation).
      # This model produces an S-shaped growth curve that stabilizes at K.
      
      # 1. Define Parameters
      N_t <- n_ind                    # Current Density (Post-selection)
      r   <- intrinsic_growth_rate    # Max growth rate (per capita)
      K   <- max_population_size      # Carrying Capacity
      
      # 2. Calculate Growth Increment (dN)
      # The discrete logistic difference equation:
      # dN = r * N_t * (1 - N_t / K)
      # If N < K, growth is positive. If N > K, growth is negative.
      growth_increment <- r * N_t * (1 - (N_t / K))
      
      # 3. Calculate Target Size (N_t+1)
      N_target_next_gen <- N_t + growth_increment
      
      # 4. Recruitment Calculation
      # We assume overlapping generations where parents (Age 1) survive to Age 2.
      # Offspring Needed = Target Total - Surviving Adults
      n_survivors <- N_t
      n_offspring_needed <- round(N_target_next_gen - n_survivors)
      
      # Safety: Cannot recruit negative offspring
      n_offspring_needed <- max(0, n_offspring_needed)
      
      # ========================================================================
      
      # --- C. MATING POOL FORMATION ---
      # We split the population into Females (Mothers) and Males (Fathers).
      # Currently assumes a 1:1 operational sex ratio.
      mat_size <- floor(n_ind / 2)
      if (mat_size < 1) break
      
      all_idx <- 1:n_ind
      mat_idx <- sample(all_idx, mat_size) # Mothers
      pat_idx <- all_idx[-mat_idx]         # Fathers
      
      moms <- pop_data[mat_idx, , drop=FALSE]
      dads <- pop_data[pat_idx, , drop=FALSE]
      n_moms <- nrow(moms)
      
      actual_offspring <- 0
      offspring_pop <- NULL
      
      # --- D. MATE CHOICE (Assortative Mating) ---
      # Biological Context: Females choose males based on "Gene Tinder" rules:
      # 1. Are they close by? (Spatial)
      # 2. Do they look like me? (Phenotypic Assortment)
      # 3. Are they genetically similar? (Genetic Assortment)
      # 
      
      if (n_moms > 0 && n_offspring_needed > 0) {
        
        # Distribute reproductive load: Each mother gets a share of the total offspring needed.
        base_spawn <- floor(n_offspring_needed / n_moms)
        remainder <- n_offspring_needed %% n_moms
        spawns <- rep(base_spawn, n_moms)
        if (remainder > 0) spawns[sample(n_moms, remainder)] <- spawns[sample(n_moms, remainder)] + 1
        
        # 1. Vectorized Trait Comparisons
        # 'outer' computes the difference matrix for every female-male pair instantly.
        x_d <- abs(outer(moms[,"x_coordinate"], dads[,"x_coordinate"], "-"))
        y_d <- abs(outer(moms[,"y_coordinate"], dads[,"y_coordinate"], "-"))
        z_d <- abs(outer(moms[,"z_coordinate"], dads[,"z_coordinate"], "-"))
        q_d <- abs(outer(moms[,"q_score"], dads[,"q_score"], "-"))
        p_d <- abs(outer(moms[,"phenotype_score"], dads[,"phenotype_score"], "-"))
        
        # 2. Calculate 3D Spatial Distance (Euclidian Distance)
        dist_mat <- sqrt(x_d^2 + y_d^2 + z_d^2)
        
        # 3. Preference Functions (Exponential Decay)
        # Simulates a preference kernel: Attraction drops exponentially as difference increases.
        s_dist <- exp(-k_dist * dist_mat)
        s_q <- exp(-k_q * q_d)
        s_p <- exp(-k_p * p_d)
        
        # 4. Weighted Suitability Score
        total_score <- (weight_dist * s_dist) + (weight_q * s_q) + (weight_p * s_p)
        
        # 5. Micro-Allopatry / Fine Scale Assortment
        if (fine_scale_assortment) {
          is_close <- (x_d <= 1) | (y_d <= 1) | (z_d <= 1)
          total_score[is_close] <- total_score[is_close] * fine_scale_multiplier
        }
        
        # 6. Convert Scores to Probabilities
        total_score <- total_score + weight_random # Noise floor prevents division by zero
        row_t <- rowSums(total_score)
        if(any(row_t == 0)) { total_score[row_t==0,] <- 1; row_t <- rowSums(total_score) }
        probs <- total_score / row_t
        
        # 7. Mate Selection
        # Each mother picks ONE father based on her specific probability row.
        dad_picks <- sapply(1:n_moms, function(i) sample(nrow(dads), 1, prob=probs[i,]))
        
        # --- E. REPRODUCTION (Mendelian Inheritance) ---
        mom_reps <- rep(1:n_moms, times=spawns)
        actual_offspring <- length(mom_reps)
        
        if (actual_offspring > 0) {
          m_genes <- moms[mom_reps, 1:num_loci]
          d_genes <- dads[dad_picks[mom_reps], 1:num_loci]
          
          # Fast Matrix Inheritance:
          # 1. Homozygous 2 -> Passes 1 allele (100%)
          # 2. Heterozygous 1 -> Passes 1 allele (50%) or 0 allele (50%)
          m_c <- matrix(0, actual_offspring, num_loci)
          m_c[m_genes==2] <- 1
          het_m <- m_genes==1; m_c[het_m] <- sample(0:1, sum(het_m), replace=TRUE)
          
          d_c <- matrix(0, actual_offspring, num_loci)
          d_c[d_genes==2] <- 1
          het_d <- d_genes==1; d_c[het_d] <- sample(0:1, sum(het_d), replace=TRUE)
          
          off_genes <- m_c + d_c
          
          # Calculate Offspring Traits
          off_q <- rowSums(off_genes)/(num_loci*2)
          off_p_gen <- rowSums(off_genes[,pheno_loci_indices,drop=F])/(2*length(pheno_loci_indices))
          off_pheno <- (pheno_heritability*off_p_gen) + ((1-pheno_heritability)*runif(actual_offspring))
          
          # Assign Offspring Coordinates (Eco-evolutionary feedback)
          off_x <- round(pmax(0, rnorm(actual_offspring, off_q*100, x_sd)))
          off_y <- round(pmax(0, rnorm(actual_offspring, (1-off_q)*100, y_sd)))
          off_z <- round(pmax(0, rnorm(actual_offspring, (1-off_q)*100, z_sd)))
          
          # Assign Offspring IDs and Parentage
          current_max_id <- max(pop_data[, "id"])
          off_ids <- (current_max_id + 1):(current_max_id + actual_offspring)
          mom_ids <- moms[mom_reps, "id"]
          dad_ids <- dads[dad_picks[mom_reps], "id"]
          
          offspring_pop <- cbind(off_genes, q_score=off_q, phenotype_score=off_pheno, 
                                 x_coordinate=off_x, y_coordinate=off_y, z_coordinate=off_z, 
                                 age=rep(1, actual_offspring),
                                 id=off_ids, maternal_id=mom_ids, paternal_id=dad_ids)
        }
      }
      
      # --- F. DEMOGRAPHY UPDATE ---
      # 1. Aging: Surviving parents become Age 2
      survivors <- pop_data[pop_data[,"age"] < 2, , drop=FALSE]
      if (nrow(survivors) > 0) survivors[,"age"] <- survivors[,"age"] + 1
      
      # 2. Merging: Combine survivors with new offspring
      if (nrow(survivors) > 0 && !is.null(offspring_pop)) {
        colnames(offspring_pop) <- colnames(survivors)
        pop_data <- rbind(survivors, offspring_pop)
      } else if (nrow(survivors) > 0) { pop_data <- survivors
      } else if (!is.null(offspring_pop)) { pop_data <- offspring_pop
      } else { break }
      
      # --- G. SAVE SNAPSHOT (Batch Only) ---
      if (!is.null(save_dir)) {
        het <- rowSums(pop_data[, 1:num_loci, drop=FALSE] == 1) / num_loci
        # We strip the heavy genomic columns to save disk space, keeping only metadata.
        df <- data.frame(generation=current_gen, q_score=pop_data[,"q_score"], 
                         phenotype_score=pop_data[,"phenotype_score"], heterozygosity=het,
                         x=pop_data[,"x_coordinate"], y=pop_data[,"y_coordinate"], z=pop_data[,"z_coordinate"], 
                         age=pop_data[,"age"],
                         id=pop_data[,"id"], maternal_id=pop_data[,"maternal_id"], paternal_id=pop_data[,"paternal_id"])
        write_csv(df, file.path(save_dir, paste0("generation_data_G", current_gen, ".csv")))
      }
      
      current_gen <- current_gen + 1
    }
    return(pop_data)
  }
  
  
  # ============================================================================
  # LOGIC: SINGLE RUN VS PARALLEL BATCH
  # ============================================================================
  
  if (num_runs == 1) {
    # --- SINGLE RUN MODE (Interactive) ---
    cat(paste("\n--- Running Single Gene_Tinder Simulation: ", experiment_name, " ---\n"))
    result <- simulation_worker(run_id = 1, save_dir = NULL)
    cat("Simulation complete. Returning final population matrix.\n")
    return(result)
    
  } else {
    # --- BATCH MODE (Parallelized) ---
    cat(paste("\n--- Starting Gene_Tinder Batch:", experiment_name, "| Runs:", num_runs, "---"))
    if (!dir.exists(experiment_name)) dir.create(experiment_name, recursive = TRUE)
    
    # Configure Parallel Backend
    if (parallel) {
      num_cores <- availableCores() - 2
      if (num_cores < 1) num_cores <- 1
      plan(multisession, workers = num_cores)
      cat(paste("\n--- Mode: PARALLEL (Cores:", num_cores, ") ---\n"))
    } else {
      plan(sequential)
      cat(paste("\n--- Mode: SEQUENTIAL (Debug) ---\n"))
    }
    
    start_time <- Sys.time()
    
    # Distribute simulations across cores
    future_lapply(1:num_runs, function(i) {
      run_dir <- file.path(experiment_name, paste0("run_", i))
      tryCatch({
        simulation_worker(run_id = i, save_dir = run_dir)
      }, error = function(e) { print(paste("Error in run", i, ":", e$message)) })
      return(NULL)
    }, future.seed = TRUE) # Ensures reproducible randomness
    
    # Clean up parallel workers to free RAM
    plan(sequential)
    duration <- Sys.time() - start_time
    cat(paste("\n--- Batch Complete. Duration:", round(duration, 2), units(duration), "---\n"))
    
    # --- DATA AGGREGATION ---
    # Combines thousands of small generation files into one primary dataset
    cat("Aggregating results... ")
    file_pattern <- "generation_data_G\\d+\\.csv"
    csv_files <- list.files(path = experiment_name, pattern = file_pattern, recursive = TRUE, full.names = TRUE)
    
    if (length(csv_files) > 0) {
      all_data <- lapply(csv_files, function(fp) {
        rid <- as.integer(gsub(".*?run_(\\d+).*", "\\1", fp))
        d <- read_csv(fp, show_col_types = FALSE)
        d$run_id <- rid
        d
      })
      final_df <- bind_rows(all_data)
      out_file <- file.path(experiment_name, "aggregated_results_primary.csv")
      write_csv(final_df, out_file)
      cat(paste("Saved to:", out_file, "\n"))
      return(out_file)
    } else {
      warning("No data generated.")
      return(NULL)
    }
  }
}

# Parallel Batch Run (Saves CSVs):
final_csv <- Gene_Tinder(
   experiment_name = "Gene-Tinder_IBD",
   num_runs = 10, parallel = TRUE,
   weight_dist = 1.0, weight_q = 0.0, weight_p = 0.0, weight_random = 1.0,
   k_dist = 1, k_q = 1, k_p = 1,
   min_fitness_scalar = 0.7, 
   species_A_ratio = 0.5
 )
