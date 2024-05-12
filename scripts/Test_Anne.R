setwd("C:/Users/annew/Documents/One Health/Network diagrams/")#change to working drive

library(igraph)
library(harrietr)
library(dplyr) 


# read in SNP table row.names and col.names indicate that these are not part of the table data  
data <- read.delim(file = "~/Documents/UTS/APG/Retro/Manuscript_1/APG-OHEC-Retro-M1/analysis/SKA/ST69_dists/ST69_dist_100_SNP.distances.tsv", check.names = FALSE, sep = "\t")

#### MC
data <- data %>% select(`Sample 1`, `Sample 2`, SNPs) %>% rename("iso1" = `Sample 1`, "iso2" = `Sample 2`, "dist" = SNPs)



# convert SNP distances to binary ) for no link 1 for link 
edge_all <- melt_dist(data) #Melt SNP matrix removes mirrored and self pairs 
edge <- edge_all # keeping a copy bc I'm a hoarder and afraid of breaking things 
# DO THIS BEFORE CONVERTING ANY NON MATCHES TO 0. 
edge$dist [edge$dist <90] <- 1 #convert any SNP less than what ever your threshold for linked is to a 1. This needs to be < LESS THAN the SNP threshold to ensure complete matches are kept 
edge$dist [edge$dist >1] <- 0 # convert anything above the SNP threshold to 0 = no link. if you do this before converting the linked cases you will lose complete matches 
edge1 <- subset(edge, dist == 1) # keep only those with a link in the network map 
edge1$dist <- NULL # remove distance column for igraph edge list 
edge2 <- as.matrix(edge1) # convert to matrix
edge3 <- graph.edgelist(data, directed = FALSE) # convert to igraph edgelist 

# V(edge3) #view nodes
# E(edge3) #view edges 

nodes <- igraph::as_data_frame(data, what = "vertices") #convert nodes in edgelist to df, this is important bc you need the samples in the order they appear in the node list to apply colour 

#read in large metadata set 
whole_metadata <- read.delim(file = "APG-OHEC-Retro-M1/delims/genometa_n5631.txt", sep = "\t")
ST_meta <- subset(whole_metadata, select = c(name, Revised_Source_Niche)) #subset down to the columns that are needed 
nodes <- inner_join(nodes, ST_meta) #dont use merge as it will alter the order of the samples 

nodes$colour <- "" #create new column to add colour. igraph needs colours to be added as a vector individually for each node in the order of the node list 

# assign colours for each possible node type in the dataset
nodes$colour<- ifelse(nodes$Revised_Source_Niche =="Captive Animal", "#b54673", as.character(nodes$colour))
nodes$colour<- ifelse(nodes$Revised_Source_Niche =="Companion Animal", "#59398d", as.character(nodes$colour))
nodes$colour<- ifelse(nodes$Revised_Source_Niche =="Environmental", "#709b46", as.character(nodes$colour))
nodes$colour<- ifelse(nodes$Revised_Source_Niche =="Food", "#c09f3d", as.character(nodes$colour))
nodes$colour<- ifelse(nodes$Revised_Source_Niche =="Human", "#48c595", as.character(nodes$colour))
nodes$colour<- ifelse(nodes$Revised_Source_Niche =="Livestock", "#c26bbc", as.character(nodes$colour))
nodes$colour<- ifelse(nodes$Revised_Source_Niche =="Waste", "#6d83da", as.character(nodes$colour))
nodes$colour<- ifelse(nodes$Revised_Source_Niche =="Wild Animal", "#b9553d", as.character(nodes$colour))

#turn metadata dataframe into vectors to colour graph by 
##These MUST be in the order of the node list in the edge graph 
source = nodes[["Revised_Source_Niche"]] #not 100% sure this one is needed 
colour = nodes[["colour"]] 


#### adding weighted edges by SNP threshold

# convert SNP matrix into long format 
g <- melt_dist(data)
#convert SNPs into weighted distances 
data$dist [data$dist >= 0 & data$dist <= 20] <- 3 #thickest link 
data$dist [data$dist >90] <- 0 # convert any SNP greater than X to a 0 (no lnik) MAKE SURE THIS MATCHES THE THRESHOLD YOU USED TO MAKE THE EDGE LIST 
data$dist [data$dist >50] <- 1 #thinnest link 
data$dist [data$dist>20] <- 2 # medium link .. you get the idea 

#Keep only those with a weight 
# THIS SHOULD BE IN THE SAME ORDER AS THE NODE LIST! 
g1 <- subset(data, dist > 0)

g1$iso1 <- g1$iso2 <- NULL # remove everything but the weight. left with weights in order of the nodes to weight the links 

#plot network 

pdf(file = "ST69 network with ID.pdf", width=9, height=9, useDingbats=FALSE) 

plot (edge3, vertex.size=3,  vertex.color=colour, vertex.label = NA, edge.width=g1$dist, rescale = TRUE)
title(main = "E. coli ST69 SNP network")
legend("bottomleft",
       legend = c(#"Captive Animal",
         #"Companion Animal", 
         #"Environmental", 
         #"Food", 
         "Human", 
         "Livestock", 
         #"Waste", 
         "Wild Animal"),
       fill = c(#"#b54673",
         #"#59398d",
         #"#709b46",
         #"#c09f3d",
         "#48c595",
         "#c26bbc",
         #"#6d83da",
         "#b9553d"),
       pch=19, bty = "n", pt.cex = 1, cex = 1, text.col="black" , horiz = FALSE, inset(0, 0))

dev.off()

