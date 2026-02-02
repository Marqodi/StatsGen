# ==============================================================================
# FILE: Visualize_Results_Filtered.R
# DESC: Visualizes simulation variance vs empirical target (Excluding Fixation)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# --- 1. CONFIGURATION ---
input_file <- "Gene-Tinder_Random/aggregated_results_primary.csv"
empirical_alpha <- 7.1
empirical_beta  <- 2.8

# --- 2. LOAD & PREPARE DATA ---
if (!file.exists(input_file)) stop("Aggregated CSV not found.")
data <- read_csv(input_file, show_col_types = FALSE)

# Filter for the final generation only
final_gen <- max(data$generation)
final_data <- data %>% filter(generation == final_gen)

cat(paste("Processing", n_distinct(final_data$run_id), "runs from generation", final_gen, "...\n"))

# --- 3. FIT BETA DISTRIBUTIONS (FILTERING FIXATION) ---

# Custom function: Removes 0s and 1s, then fits Beta via Method of Moments
estimate_beta_filtered <- function(x) {
  # FILTER STEP: Remove fixed individuals (0.0 and 1.0)
  # We use a small epsilon tolerance to catch floating point near-matches
  x_hybrid <- x[x > 0.00001 & x < 0.99999]
  
  # Safety check: Need enough data points to calculate variance
  if (length(x_hybrid) < 2) return(c(NA, NA, 1.0)) # 1.0 = 100% fixation
  
  # Calculate Fixation Rate (for reporting purposes)
  fixation_rate <- 1 - (length(x_hybrid) / length(x))
  
  mu <- mean(x_hybrid)
  var <- var(x_hybrid)
  
  if (var < 1e-10) return(c(NA, NA, fixation_rate))
  
  # Method of Moments
  alpha <- mu * ((mu * (1 - mu) / var) - 1)
  beta  <- (1 - mu) * ((mu * (1 - mu) / var) - 1)
  
  return(c(alpha = alpha, beta = beta, fixation = fixation_rate))
}

# Calculate parameters for every run
run_stats <- final_data %>%
  group_by(run_id) %>%
  summarise(
    stats = list(estimate_beta_filtered(q_score)), 
    .groups = "drop"
  ) %>%
  mutate(
    alpha = sapply(stats, `[`, 1),
    beta  = sapply(stats, `[`, 2),
    fix_rate = sapply(stats, `[`, 3)
  )

# Report on how much data was removed
avg_fixation <- mean(run_stats$fix_rate, na.rm = TRUE)
cat(paste0("Average Fixation Rate (Removed Data): ", round(avg_fixation * 100, 2), "%\n"))

# Remove runs that failed to fit (e.g., 100% fixation)
valid_runs <- run_stats %>% filter(!is.na(alpha))

# --- 4. GENERATE DENSITY CURVES ---
x_grid <- seq(0.001, 0.999, length.out = 500) # Avoid 0 and 1 for dbeta

sim_densities <- lapply(1:nrow(valid_runs), function(i) {
  a <- valid_runs$alpha[i]
  b <- valid_runs$beta[i]
  data.frame(
    x = x_grid,
    y = dbeta(x_grid, shape1 = a, shape2 = b),
    run_id = valid_runs$run_id[i]
  )
}) %>% bind_rows()

# --- 5. AGGREGATE SIMULATION STATS ---
sim_summary <- sim_densities %>%
  group_by(x) %>%
  summarise(
    mean_y = mean(y),
    sd_y   = sd(y),
    ymin   = pmax(0, mean_y - sd_y),
    ymax   = mean_y + sd_y
  )

# --- 6. GENERATE EMPIRICAL DATA ---
empirical_df <- data.frame(
  x = x_grid,
  y = dbeta(x_grid, shape1 = empirical_alpha, shape2 = empirical_beta)
)

# --- 7. PLOTTING ---
col_sim <- "#4E79A7"
col_emp <- "#E15759"

p <- ggplot() +
  # Simulation Variance (Shadow)
  geom_ribbon(data = sim_summary, 
              aes(x = x, ymin = ymin, ymax = ymax), 
              fill = col_sim, alpha = 0.2) +
  
  # Simulation Mean
  geom_line(data = sim_summary, 
            aes(x = x, y = mean_y, color = "Simulated Distribution"), 
            linewidth = 1.2) +
  
  # Empirical Target
  geom_line(data = empirical_df, 
            aes(x = x, y = y, color = "Empirical Distribution"), 
            linewidth = 1.2, linetype = "dashed") +
  
  scale_color_manual(name = NULL, 
                     values = c("Simulated Distribution" = col_sim, 
                                "Empirical Distribution" = col_emp)) +
  labs(
    title = "Hybrid q-Score Distribution",
    #subtitle = paste0("Comparing Hybrids Only | Avg Fixation Removed: ", round(avg_fixation*100, 1), "%"),
    x = "Hybrid Index (q-Score)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.3, 0.9), # Place legend inside plot top-right
    legend.background = element_rect(fill="white", color="gray90"),
    plot.title = element_text(face = "bold", size = 14)
  )

print(p)
ggsave("Beta_Comparison_Filtered.png", plot = p, width = 8, height = 6)