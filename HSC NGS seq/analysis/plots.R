library(ggplot2)

setwd("c:/00_Berman_lab_project/side_projects/tet2-mutant/tet2_human_sgRNAs/New_data/analysis")

data <- read.csv("first_match_coordinates.csv")

ggplot(data, aes(x= FirstCoordinate, fill=Result)) +
  geom_histogram(alpha=0.5, bins = 60, position="identity", color="black") +
  geom_density(alpha=0.5) +
  scale_y_sqrt() +
  labs(x= "FirstCoordinate", y = "Number of events") +
  guides(fill=guide_legend(title="Plot of mutation event coordinates"))+ 
  facet_wrap(~sgRNA)