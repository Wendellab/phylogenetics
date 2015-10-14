#Chloroplast phylogenetics/diversity for Wake Island cotton
#This is a cleaned up version of wilkes.phylogenetics.R
#Basically, this code takes in a NEXUS file of concatenated sequences, converts
#it to a distance matrix, and makes a NJ tree (to show phylogeny) and a PCOA
#plot to show diversity of accessions

library(ape)#sequence handling
library(ggplot2)#visualization

#read in sequence file
cp.nexus <- read.nexus.data("AD7_cp_concatenated.no_ps24.nex")
dist.matrix <- dist.dna(as.DNAbin(cp.nexus),model="K80",pairwise.deletion=TRUE, as.matrix=TRUE)

#determine which sequences can actually be compared and keep just those
keepers <- complete.cases(dist.matrix)
reduced.dist.matrix <- dist.matrix[keepers,keepers]

#generate Principal Coordinates Analysis object and plot
pcoa <- pcoa(reduced.dist.matrix)
accessions <- gsub("'","",dimnames(pcoa$vector)[[1]])
temp.list <- strsplit(accessions, split = "_")
species <- sapply(temp.list,`[`,1)
genome <- c(rep("A",3),rep("ADX",16),rep("AD7",14))
ggplot(as.data.frame(pcoa$vectors),aes(x=Axis.3,y=Axis.4,colour=species,shape = genome)) + 
  geom_jitter(size = 4, position = position_jitter(height = 0.00001,width = 0.00001)) +
  xlab(label = "Axis 3") + ylab(label = "Axis 4") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24)) +
  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16),legend.title = element_text(size = 16)) +
  scale_shape_discrete(guide = FALSE)
