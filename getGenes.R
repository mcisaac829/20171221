setwd('/Users/scott/Google Drive/Calico Experiments/20171214_AgingSheet')
library(tidyverse)
library(seqinr)

#load data
df <- read.delim(file = 'ESR_aging_masterSheet.txt', header = TRUE)
df <- as.tibble(df)

# load full promoter list
promoters <- read.fasta('YeastPromoters.txt')

#remove 0-length promoters
promoters <- promoters[as.vector(unlist((lapply(promoters, function(x){x[[1]][1] != ">"}))))] 


#get categories
category_matrix <- df %>% group_by(category) %>% summarise()
category_vec <- category_matrix$category
category_vec <- as.vector(category_vec[which(category_vec != "")])

#make fold folder for each category
for(j in 1:length(category_vec)){
  dir.create(category_vec[j])
}

#write data for each category to its own folder
for(j in 1:length(category_vec)){
  mydir <- paste0("/Users/scott/Google Drive/Calico Experiments/20171214_AgingSheet/",category_vec[j])
  setwd(mydir)
  
  data <- df %>% 
    filter(category == category_vec[j])
  
  genes_only <- df %>%
    filter(category == category_vec[j]) %>%
    select(gene_id)
  genes_only <- as.character(genes_only$gene_id)
  
  write.table(data,file=paste0(category_vec[j],'.txt'),sep="\t",quote=FALSE) #write data for category
  
  subset_promoters <- promoters[attributes(promoters)$names %in% genes_only]
  write.fasta(subset_promoters, names = names(subset_promoters),file.out = paste0(category_vec[j],'_promoters.txt'))
  
  system(paste0("export PATH=$HOME/meme/bin:$PATH && meme ",category_vec[j],"_promoters.txt -dna -revcomp -nmotifs 5 -maxw 12 -oc ", category_vec[j], "_meme -maxsize 300000"),intern=FALSE)
  
}

