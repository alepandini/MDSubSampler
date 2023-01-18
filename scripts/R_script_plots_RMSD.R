###############################################################################
# Read replica 5 and different samples 
###############################################################################
R5 <- read.table("5_0.5%_RMSD_Bhattacharyya.dat")
R5_0.05 <- read.table("5_0.05%_RMSD_sample_Bhattacharyya.dat")
R5_0.1 <- read.table("5_0.1%_RMSD_sample_Bhattacharyya.dat")
R5_0.25 <- read.table("5_0.25%_RMSD_sample_Bhattacharyya.dat")
R5_0.5 <- read.table("5_0.5%_RMSD_sample_Bhattacharyya.dat")
R5_1 <- read.table("5_1%_RMSD_sample_Bhattacharyya.dat")
R5_2.5 <- read.table("5_2.5%_RMSD_sample_Bhattacharyya.dat")
R5_5 <- read.table("5_5%_RMSD_sample_Bhattacharyya.dat")
R5_10 <- read.table("5_10%_RMSD_sample_Bhattacharyya.dat")
R5_20 <- read.table("5_20%_RMSD_sample_Bhattacharyya.dat")
R5_25 <- read.table("5_25%_RMSD_sample_Bhattacharyya.dat")
R5_50 <- read.table("5_50%_RMSD_sample_Bhattacharyya.dat")

###############################################################################
# Read replica 1 to 5 and different samples 
###############################################################################
R1_R5 <- read.table("1_5_0.5%_RMSD_Bhattacharyya.dat")
R1_R5_0.05 <- read.table("1_5_0.05%_RMSD_Bhattacharyya.dat")
R1_R5_0.1 <- read.table("1_5_0.1%_RMSD_Bhattacharyya.dat")
R1_R5_0.25 <- read.table("1_5_0.25%_RMSD_Bhattacharyya.dat")
R1_R5_0.5 <- read.table("1_5_0.5%_RMSD_sample_Bhattacharyya.dat")
R1_R5_1 <- read.table("1_5_1%_RMSD_Bhattacharyya.dat")
R1_R5_2.5 <- read.table("1_5_2.5%%_RMSD_Bhattacharyya.dat")
R1_R5_5 <- read.table("1_5_5%_RMSD_sample_Bhattacharyya.dat")
R1_R5_10 <- read.table("1_5_10%_RMSD_sample_Bhattacharyya.dat")
R1_R5_20 <- read.table("1_5_20%_RMSD_sample_Bhattacharyya.dat")
R1_R5_25 <- read.table("1_5_25%_RMSD_sample_Bhattacharyya.dat")
R1_R5_50 <- read.table("1_5_50%_RMSD_sample_Bhattacharyya.dat")

###############################################################################
# Read replica 6 to 10 and different samples 
###############################################################################
R6_R10 <- read.table("6_10_0.5%_RMSD_Bhattacharyya.dat")
R6_R5_0.05 <- read.table("6_10_0.05%_RMSD_sample_Bhattacharyya.dat")
R6_R5_0.1 <- read.table("6_10_0.1%_RMSD_sample_Bhattacharyya.dat")
R6_R5_0.25 <- read.table("6_10_0.25%_RMSD_sample_Bhattacharyya.dat")
R6_R5_0.5 <- read.table("6_10_0.5%_RMSD_sample_Bhattacharyya.dat")
R6_R5_1 <- read.table("6_10_1%_RMSD_sample_Bhattacharyya.dat")
R6_R5_2.5 <- read.table("6_10_2.5%_RMSD_sample_Bhattacharyya.dat")
R6_R10_5 <- read.table("6_10_5%_RMSD_sample_Bhattacharyya.dat")
R6_R10_10 <- read.table("6_10_10%_RMSD_sample_Bhattacharyya.dat")
R6_R10_20 <- read.table("6_10_20%_RMSD_sample_Bhattacharyya.dat")
R6_R5_25 <- read.table("6_10_25%_RMSD_sample_Bhattacharyya.dat")
R6_R5_50 <- read.table("6_10_50%_RMSD_sample_Bhattacharyya.dat")

###############################################################################
# Read replica 1 to 10 and different samples 
###############################################################################
R1_R10 <- read.table("1_10_0.5%_RMSD_Bhattacharyya.dat")
R1_R5_0.05 <- read.table("1_10_0.05%_RMSD_sample_Bhattacharyya.dat")
R1_R5_0.1 <- read.table("1_10_0.1%_RMSD_sample_Bhattacharyya.dat")
R1_R5_0.25 <- read.table("1_10_0.25%_RMSD_sample_Bhattacharyya.dat")
R1_R5_0.5 <- read.table("1_10_0.5%_RMSD_sample_Bhattacharyya.dat")
R1_R5_1 <- read.table("1_10_1%_RMSD_sample_Bhattacharyya.dat")
R1_R5_2.5 <- read.table("1_10_2.5%_RMSD_sample_Bhattacharyya.dat")
R1_R10_5 <- read.table("1_10_5%_RMSD_sample_Bhattacharyya.dat")
R1_R10_10 <- read.table("1_10_10%_RMSD_sample_Bhattacharyya.dat")
R1_R10_20 <- read.table("1_10_20%_RMSD_sample_Bhattacharyya.dat")
R1_R5_25 <- read.table("1_10_25%_RMSD_sample_Bhattacharyya.dat")
R1_R5_50 <- read.table("1_10_50%_RMSD_sample_Bhattacharyya.dat")

###############################################################################
# Function that generates plots with distributions
###############################################################################
ref_hist <- function(ref_prop_values, prop_values, plot_title, n_breaks = 50){
  min_value = min(c(min(ref_prop_values),min(prop_values)))
  max_value = max(c(max(ref_prop_values),max(prop_values)))
  breaks_seq = seq(min_value, max_value, (max_value - min_value)/n_breaks)
  hist(
    ref_prop_values,
    freq=FALSE, 
    breaks = breaks_seq,
    ylim = c(0,0.8),
    border = NA,
    col = rgb(1,0,0,0.3),
    xlab = '',
    main = plot_title,
  )
  hist(
    prop_values, 
    freq=FALSE, 
    breaks=breaks_seq, 
    ylim = c(0,0.8), 
    border = NA,
    col = rgb(0,0,1,0.3),
    xlab = '',
    main = '',
    add=TRUE
  )  
  legend(
    'topright',
    c("reference", "sample"),
    pch = 15,
    col = c(
      rgb(1,0,0,0.3),
      rgb(0,0,1,0.3)
    ),
    bty = 'n'
  )
}

###############################################################################
# Create plots in pdf file for replica 5
###############################################################################
pdf(file= "R5.pdf" )
par(mfrow = c(2, 2))
ref_hist(R5$V2, R5_0.5$V2, "0.5%")
ref_hist(R5$V2, R5_5$V2, "5%")
ref_hist(R5$V2, R5_10$V2, "10%")
ref_hist(R5$V2, R5_20$V2, "20%")
dev.off()

###############################################################################
# Create plots in pdf file for replica 1 to 5
###############################################################################
pdf(file= "R1_R5.pdf" )
par(mfrow = c(2, 2))
ref_hist(R1_R5$V2, R1_R5_0.5$V2, "0.5%")
ref_hist(R1_R5$V2, R1_R5_5$V2, "5%")
ref_hist(R1_R5$V2, R1_R5_10$V2, "10%")
ref_hist(R1_R5$V2, R1_R5_20$V2, "20%")
dev.off()

###############################################################################
# Create plots in pdf file for replica 6 to 10
###############################################################################
pdf(file= "R6_R10.pdf" )
par(mfrow = c(2, 2))
ref_hist(R6_R10$V2, R6_R10_0.5$V2, "0.5%")
ref_hist(R6_R10$V2, R6_R10_5$V2, "5%")
ref_hist(R6_R10$V2, R6_R10_10$V2, "10%")
ref_hist(R6_R10$V2, R6_R10_20$V2, "20%")
dev.off()

###############################################################################
# Create plots in pdf file for replica 1 to 10
###############################################################################
pdf(file= "R1_R10.pdf" )
par(mfrow = c(2, 2))
ref_hist(R1_R10$V2, R1_R10_0.5$V2, "0.5%")
ref_hist(R1_R10$V2, R1_R10_5$V2, "5%")
ref_hist(R1_R10$V2, R1_R10_10$V2, "10%")
ref_hist(R1_R10$V2, R1_R10_20$V2, "20%")
dev.off()