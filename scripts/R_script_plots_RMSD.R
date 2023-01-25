###############################################################################
# Read replica 5 and different samples 
###############################################################################
R5 <- read.table("5_RMSD_Bhattacharyya.dat")
R5_0.05 <- read.table("5_0.05_RMSD_Bhattacharyya.dat")
R5_0.1 <- read.table("5_0.1_RMSD_Bhattacharyya.dat")
R5_0.25 <- read.table("5_0.25_RMSD_Bhattacharyya.dat")
R5_0.5 <- read.table("5_0.5_RMSD_Bhattacharyya.dat")
R5_1 <- read.table("5_1.0_RMSD_Bhattacharyya.dat")
R5_2.5 <- read.table("5_2.5_RMSD_Bhattacharyya.dat")
R5_5 <- read.table("5_5.0_RMSD_Bhattacharyya.dat")
R5_10 <- read.table("5_10.0_RMSD_Bhattacharyya.dat")
R5_20 <- read.table("5_20.0_RMSD_Bhattacharyya.dat")
R5_25 <- read.table("5_25.0_RMSD_Bhattacharyya.dat")
R5_50 <- read.table("5_50.0_RMSD_Bhattacharyya.dat")

###############################################################################
# Read replica 1 to 5 and different samples 
###############################################################################
R1_R5 <- read.table("1_5_RMSD_Bhattacharyya.dat")
R1_R5_0.05 <- read.table("1_5_0.05_RMSD_Bhattacharyya.dat")
R1_R5_0.1 <- read.table("1_5_0.1_RMSD_Bhattacharyya.dat")
R1_R5_0.25 <- read.table("1_5_0.25_RMSD_Bhattacharyya.dat")
R1_R5_0.5 <- read.table("1_5_0.5_RMSD_Bhattacharyya.dat")
R1_R5_1 <- read.table("1_5_1.0_RMSD_Bhattacharyya.dat")
R1_R5_2.5 <- read.table("1_5_2.5_RMSD_Bhattacharyya.dat")
R1_R5_5 <- read.table("1_5_5.0_RMSD_Bhattacharyya.dat")
R1_R5_10 <- read.table("1_5_10.0_RMSD_Bhattacharyya.dat")
R1_R5_20 <- read.table("1_5_20.0_RMSD_Bhattacharyya.dat")
R1_R5_25 <- read.table("1_5_25.0_RMSD_Bhattacharyya.dat")
R1_R5_50 <- read.table("1_5_50.0_RMSD_Bhattacharyya.dat")

###############################################################################
# Read replica 6 to 10 and different samples 
###############################################################################
R6_R10 <- read.table("6_10_RMSD_Bhattacharyya.dat")
R6_R10_0.05 <- read.table("6_10_0.05_RMSD_Bhattacharyya.dat")
R6_R10_0.1 <- read.table("6_10_0.1_RMSD_Bhattacharyya.dat")
R6_R10_0.25 <- read.table("6_10_0.25_RMSD_Bhattacharyya.dat")
R6_R10_0.5 <- read.table("6_10_0.5_RMSD_Bhattacharyya.dat")
R6_R10_1 <- read.table("6_10_1.0_RMSD_Bhattacharyya.dat")
R6_R10_2.5 <- read.table("6_10_2.5_RMSD_Bhattacharyya.dat")
R6_R10_5 <- read.table("6_10_5.0_RMSD_Bhattacharyya.dat")
R6_R10_10 <- read.table("6_10_10.0_RMSD_Bhattacharyya.dat")
R6_R10_20 <- read.table("6_10_20.0_RMSD_Bhattacharyya.dat")
R6_R10_25 <- read.table("6_10_25.0_RMSD_Bhattacharyya.dat")
R6_R10_50 <- read.table("6_10_50.0_RMSD_Bhattacharyya.dat")

###############################################################################
# Read replica 1 to 10 and different samples 
###############################################################################
R1_R10 <- read.table("1_10_RMSD_Bhattacharyya.dat")
R1_R10_0.05 <- read.table("1_10_0.05_RMSD_Bhattacharyya.dat")
R1_R10_0.1 <- read.table("1_10_0.1_RMSD_Bhattacharyya.dat")
R1_R10_0.25 <- read.table("1_10_0.25_RMSD_Bhattacharyya.dat")
R1_R10_0.5 <- read.table("1_10_0.5_RMSD_Bhattacharyya.dat")
R1_R10_1 <- read.table("1_10_1.0_RMSD_Bhattacharyya.dat")
R1_R10_2.5 <- read.table("1_10_2.5_RMSD_Bhattacharyya.dat")
R1_R10_5 <- read.table("1_10_5.0_RMSD_Bhattacharyya.dat")
R1_R10_10 <- read.table("1_10_10.0_RMSD_Bhattacharyya.dat")
R1_R10_20 <- read.table("1_10_20.0_RMSD_Bhattacharyya.dat")
R1_R10_25 <- read.table("1_10_25.0_RMSD_Bhattacharyya.dat")
R1_R10_50 <- read.table("1_10_50.0_RMSD_Bhattacharyya.dat")

###############################################################################
# Function that generates plots with distributions
###############################################################################
ref_hist <- function(ref_prop_values, prop_values, plot_title, n_breaks = 30){
  min_value = min(c(min(ref_prop_values),min(prop_values)))
  max_value = max(c(max(ref_prop_values),max(prop_values)))
  breaks_seq = seq(min_value, max_value, (max_value - min_value)/n_breaks)
  hist(
    ref_prop_values,
    freq=FALSE, 
    breaks = breaks_seq,
    ylim = c(0,0.05),
    border = NA,
    col = rgb(1,0,0,0.3),
    xlab = '',
    main = plot_title,
  )
  hist(
    prop_values, 
    freq=FALSE, 
    breaks=breaks_seq, 
    ylim = c(0,0.05), 
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
par(mfrow = c(2,2))
ref_hist(R5$V2, R5_0.05$V2, "0.05%    1.22327")
ref_hist(R5$V2, R5_0.1$V2, "0.1%    0.86817")
ref_hist(R5$V2, R5_0.25$V2, "0.25%    0.62658")
ref_hist(R5$V2, R5_0.5$V2, "0.5%    0.33794")
ref_hist(R5$V2, R5_1$V2, "1%    0.18858")
ref_hist(R5$V2, R5_2.5$V2, "2.5%    0.05975")
ref_hist(R5$V2, R5_5$V2, "5%    0.02770")
ref_hist(R5$V2, R5_10$V2, "10%    0.01491")
ref_hist(R5$V2, R5_20$V2, "20%    0.00799")
ref_hist(R5$V2, R5_25$V2, "25%    0.00561")
ref_hist(R5$V2, R5_50$V2, "50%    0.00278")
dev.off()

###############################################################################
# Create plots in pdf file for replica 1 to 5
###############################################################################
pdf(file= "R1_R5.pdf" )
par(mfrow = c(2,2))
ref_hist(R1_R5$V2, R1_R5_0.05$V2, "0.05%    0.52670")
ref_hist(R1_R5$V2, R1_R5_0.1$V2, "0.1%    0.36049")
ref_hist(R1_R5$V2, R1_R5_0.25$V2, "0.25%    0.18348")
ref_hist(R1_R5$V2, R1_R5_0.5$V2, "0.5%    0.07094")
ref_hist(R1_R5$V2, R1_R5_1$V2, "1%    0.02361")
ref_hist(R1_R5$V2, R1_R5_2.5$V2, "2.5%    0.01027")
ref_hist(R1_R5$V2, R1_R5_5$V2, "5%    0.00574")
ref_hist(R1_R5$V2, R1_R5_10$V2, "10%    0.00227")
ref_hist(R1_R5$V2, R1_R5_20$V2, "20%    0.00129")
ref_hist(R1_R5$V2, R1_R5_25$V2, "25%    0.00118")
ref_hist(R1_R5$V2, R1_R5_50$V2, "50%    0.00070")
dev.off()

###############################################################################
# Create plots in pdf file for replica 6 to 10
###############################################################################
pdf(file= "R6_R10.pdf")
par(mfrow = c(2, 2))
ref_hist(R6_R10$V2, R6_R10_0.05$V2, "0.05%    0.48561")
ref_hist(R6_R10$V2, R6_R10_0.1$V2, "0.1%    0.35384")
ref_hist(R6_R10$V2, R6_R10_0.25$V2, "0.25%    0.18276")
ref_hist(R6_R10$V2, R6_R10_0.5$V2, "0.5%    0.07635")
ref_hist(R6_R10$V2, R6_R10_1$V2, "1%    0.02893")
ref_hist(R6_R10$V2, R6_R10_2.5$V2, "2.5%    0.00953")
ref_hist(R6_R10$V2, R6_R10_5$V2, "5%    0.00562")
ref_hist(R6_R10$V2, R6_R10_10$V2, "10%    0.00263")
ref_hist(R6_R10$V2, R6_R10_20$V2, "20%    0.00145")
ref_hist(R6_R10$V2, R6_R10_25$V2, "25%    0.00104")
ref_hist(R6_R10$V2, R6_R10_50$V2, "50%    0.00048")
dev.off()

###############################################################################
# Create plots in pdf file for replica 1 to 10
###############################################################################
pdf(file= "R1_R10.pdf" )
par(mfrow = c(2,2))
ref_hist(R1_R10$V2, R1_R10_0.05$V2, "0.05%    0.34288")
ref_hist(R1_R10$V2, R1_R10_0.1$V2, "0.1%    0.20030")
ref_hist(R1_R10$V2, R1_R10_0.25$V2, "0.25%    0.06485")
ref_hist(R1_R10$V2, R1_R10_0.5$V2, "0.5%    0.04043")
ref_hist(R1_R10$V2, R1_R10_1$V2, "1%    0.01741")
ref_hist(R1_R10$V2, R1_R10_2.5$V2, "2.5%    0.00558")
ref_hist(R1_R10$V2, R1_R10_5$V2, "5%    0.00365")
ref_hist(R1_R10$V2, R1_R10_10$V2, "10%    0.00143")
ref_hist(R1_R10$V2, R1_R10_20$V2, "20%    0.00078")
ref_hist(R1_R10$V2, R1_R10_25$V2, "25%    0.00061")
ref_hist(R1_R10$V2, R1_R10_50$V2, "50%    0.00035")
dev.off()
