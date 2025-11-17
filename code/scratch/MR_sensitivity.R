# --- Setup ---
# Load required packages
library(MendelianRandomization)

# --- 1. Simulate Example Data (Replace with your real data in practice) ---
set.seed(123)
n_snps <- 50
b_exp   <- rnorm(n_snps, mean = 0.1, sd = 0.05)         # SNP effects on exposure
se_exp  <- runif(n_snps, min = 0.01, max = 0.03)        # SE of exposure effects
causal_effect <- 0.5
b_out   <- rnorm(n_snps, mean = b_exp * causal_effect, sd = 0.08)  # SNP effects on outcome
se_out  <- runif(n_snps, min = 0.02, max = 0.05)        # SE of outcome effects
snp_names <- paste0("rs", 1:n_snps)

# Create MR input object
mr_data <- mr_input(bx = b_exp, bxse = se_exp,
                    by = b_out, byse = se_out,
                    snps = snp_names)

# --- 2. Perform MR Analysis ---
mr_ivw_res   <- mr_ivw(mr_data)
mr_egger_res <- mr_egger(mr_data)

# --- 3. Visualization ---

# 3.1 Forest Plot of SNP-level causal estimates
png("forest_plot_ivw.png", width = 800, height = 600, res = 100)
mr_forest(mr_data)
dev.off()

# 3.2 Funnel Plot to assess directional pleiotropy
png("funnel_plot.png", width = 800, height = 600, res = 100)
mr_funnel(mr_data)
dev.off()

# 3.3 MR Scatter Plot
png("mr_scatter_plot.png", width = 800, height = 600, res = 100)
mr_plot(mr_data)
dev.off()

# 3.4 Leave-One-Out Sensitivity Plot (Optional but recommended)


png("leave_one_out_plot.png", width = 800, height = 600, res = 100)
mr_loo(mr_data)
dev.off()
