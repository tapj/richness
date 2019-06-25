rat_data = read.csv2("rat_richness_start.csv")



source("rat_mixer.r")

rat_groups = rat_mixer(weight=rat_data$Weight,microbiota=rat_data$Bacteroides, filiation  = rat_data$Filiation, cage_provider=as.factor(as.character(rat_data$Cage_provider)), groups_nb=7,nb_filiation_max = 2, nb_cage_max=3, plot=TRUE)

#statistical check
pairwise.wilcox.test(rat_data$Weight,rat_groups, p.adjust="none")
pairwise.wilcox.test(rat_data$Bacteroides,rat_groups, p.adjust="none")


write.csv2(data.frame(rat_data, rat_groups), file="rat_data_groups.csv")