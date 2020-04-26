setwd("c:/00_Berman_lab_project/side_projects/tet2-mutant/tet2_human_sgRNAs/New_data/analysis")

data <- read.csv('translations.csv')


sgRNA1_data <- data[data$sgRNA == 'sgR1',]
sgRNA2_data <- data[data$sgRNA == 'sgR2',]
sgRNA3_data <- data[data$sgRNA == 'sgR3',]

eff_sgR1 <- mean(sgRNA1_data$Protein_size < 200)
eff_sgR2 <- mean(sgRNA2_data$Protein_size < 200)
eff_sgR3 <- mean(sgRNA3_data$Protein_size < 200)

efficiencies <- data.frame(sgRNA = c("sgRNA-1", "sgRNA-2", "sgRNA-3"), efficiency = c(eff_sgR1, eff_sgR2, eff_sgR3))

efficiencies

# > efficiencies
# sgRNA efficiency
# 1 sgRNA-1  0.4280841
# 2 sgRNA-2  0.5823075
# 3 sgRNA-3  0.3858569