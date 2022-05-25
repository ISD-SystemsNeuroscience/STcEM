library(FNN)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
#Test on real data
merfish_b2s20 = readRDS("data/merfish_b3s20.rds")

#Subset tolesioned area
lesion_b2s20 = subset(merfish_b2s20, region.lesion %in% c("lesion.core", "lesion.edge.inner", "lesion.edge.outer"))

DimPlot(lesion_b2s20, reduction = "spatial")

"%notin%" = Negate("%in%")

#--------------------------------------Calculate enrichment of k neighbors vs all cells per one section--------------------------------------------
"%notin%" = Negate("%in%")
query.cell.types = c(unique(lesion_b3s20$cell.cluster_manual)) #will search for neighbors of these cells

AllCells = lesion_b3s20@meta.data %>% select(center_x, center_y, sample.name, cell.type_manual, cell.cluster_manual)
head(AllCells)
results.list = list()
k = 6 #specify number of neighbors you want to search + 1 Specifying k=11, will find 10 neighbors (one self-neighbor is deleted in the process)
for(query.cell.type in query.cell.types){
  if(all(query.cell.types %in% lesion_b3s20$cell.cluster_manual)){
    
  }else {
    stop("One or more query cell types cannot be found in the data")
  }
  
  query.cells = lesion_b3s20@meta.data %>% filter(cell.cluster_manual == query.cell.type) %>% 
    select(center_x, center_y, sample.name, cell.type_manual, cell.cluster_manual)
  head(query.cells)
  
  result = get.knnx(data = AllCells[,1:2], query = query.cells[,1:2], k = k) #find k nearest neighbors of query.cells among AllCells
  
  all.result = list()
  for(i in 1:nrow(query.cells)){
    filter = result$nn.index[i,]
    neighbors = AllCells[filter[-1], ] # get dataframe of neighbors (exclude 1st row cause that's the same cell)
    all.result[[i]] = neighbors
  }
  names(all.result) = rownames(query.cells)
  abc = bind_rows(all.result) #merge neighborhing cells into one dataframe
  temp1 <- as.data.frame(table(abc$cell.cluster_manual)) #count occurences of clusters in neighborhoods
  temp1$Frequency = temp1$Freq/sum(temp1$Freq) #add column of frequencies of occurences of clusters in neighborhoods
  head(temp1)
  
  #Calculate statistics
  allcells = lesion_b3s20@meta.data %>% 
    select(center_x, center_y, sample.name, cell.type_manual, cell.cluster_manual)
  head(allcells)
  intersect(allcells$sample.name, abc$sample.name)  #must be empty
  temp2 <- as.data.frame(table(allcells$cell.cluster_manual)) #count occurences of clusters in neighborhoods
  temp2$Frequency = temp2$Freq/sum(temp2$Freq) #add column of frequencies of occurences of clusters in neighborhoods
  head(temp2)
  
  temp3 = left_join(temp1, temp2, by = "Var1")
  temp3$enrichment = temp3$Frequency.x / temp3$Frequency.y
  temp3$log2FC.enrichment = log2(temp3$enrichment) #log2fold change in frequency of cell type in neighborhood vs non-neighborhood
  colnames(temp3) = c("cell", "Count.neighbors", "Freq.neighbors", "Count.allcells", "Freq.allcells", "Enrichment", "log2FC.enrichment")
  head(temp3)
  
  abc$neighbor = rep("neighbor", nrow(abc))
  allcells$neighbor = rep("non.neighbor", nrow(allcells))
  
  for(cell.type in unique(temp3$cell)){
    abc[[cell.type]] = rep(ifelse(abc$cell.cluster_manual == cell.type, cell.type, "Other_cell"))#Ref
    allcells[[cell.type]] = rep(ifelse(allcells$cell.cluster_manual  == cell.type, cell.type, "Other_cell"))#Ref
  }
  
  fish.dat = rbind(abc, allcells)
  
  fisher.results = c()
  for(cell.type in unique(temp3$cell)){
    test = fisher.test(fish.dat$neighbor, fish.dat[[cell.type]], alternative = "greater")
    res.df = data.frame(cell = cell.type, p.val = test$p.value, odds.ratio = test$estimate)
    fisher.results = rbind(fisher.results, res.df)
  }
  results.df = left_join(temp3, fisher.results, by = "cell")
  results.list[[query.cell.type]] = results.df
}


results.measured = bind_rows(results.list, .id = "query_cell")

results.measured$p_log10 = -log10(results.measured$p.val)
write.table(results.measured, file = "results_fisher.txt", sep = "\t", quote = F, row.names = F)


#--------------------------------------Bootstrapped enrichment of k neighbors---------------------------------------------
library(parallel)
query.cell.types = c("Oligo_Interferon", "T-cells", "Micro_Homeostatic", "Micro_Activated_1", "Micro_Activated_2", "Micro_Interferon") #will search for neighbors of these cells
AllCells = lesion_b3s20@meta.data %>% select(center_x, center_y, sample.name, cell.type_manual, cell.cluster_manual)
head(AllCells)
results.list = list()
k = 6 #specify number of neighbors you want to search + 1 Specifying k=11, will find 10 neighbors (one self-neighbor is deleted in the process)
nperm = 10000
#library(parallel)
bootstrapped.res.list = list()
bootstrapped.res.list = mclapply(c(1:nperm),mc.cores = 14, function(x){

results.list = list()
AllCells$cell.cluster_manual = AllCells$cell.cluster_manual[sample(nrow(AllCells))]
print(paste0("working on permutation", x))

for(query.cell.type in query.cell.types){
  query.cells = AllCells %>% filter(cell.cluster_manual == query.cell.type) %>% 
    select(center_x, center_y, sample.name, cell.type_manual, cell.cluster_manual)
  #head(query.cells)
  #plot(query.cells$center_x, query.cells$center_y)
  result = get.knnx(data = AllCells[,1:2], query = query.cells[,1:2], k = k) #find k neares neighbors of query.cells among AllCells
  
  all.result = list()
  for(i in 1:nrow(query.cells)){
    filter = result$nn.index[i,]
    neighbors = AllCells[filter[-1], ] # get dataframe of neighbors (exclude 1st row cause that's the same cell)
    all.result[[i]] = neighbors
  }
  names(all.result) = rownames(query.cells)
  abc = bind_rows(all.result) #merge neighborhing cells into one dataframe
  temp1 <- as.data.frame(table(abc$cell.cluster_manual)) #count occurences of clusters in neighborhoods
  temp1$Frequency = temp1$Freq/sum(temp1$Freq) #add column of frequencies of occurences of clusters in neighborhoods
  #head(temp1)
  
  allcells = lesion_b3s20@meta.data %>% 
    select(center_x, center_y, sample.name, cell.type_manual, cell.cluster_manual)
  #print(head(allcells))
  temp2 <- as.data.frame(table(allcells$cell.cluster_manual)) #count occurences of clusters in neighborhoods
  temp2$Frequency = temp2$Freq/sum(temp2$Freq) #add column of frequencies of occurences of clusters in neighborhoods
  #head(temp2)
  
  temp3 = left_join(temp1, temp2, by = "Var1")
  temp3$enrichment = temp3$Frequency.x / temp3$Frequency.y
  temp3$log2FC.enrichment = log2(temp3$enrichment) #log2fold change in frequency of cell type in neighborhood vs non-neighborhood
  colnames(temp3) = c("cell", "Count.neighbors", "Freq.neighbors", "Count.allcells", "Freq.allcells", "Enrichment", "log2FC.enrichment")
  #head(temp3)
  
  abc$neighbor = rep("neighbor", nrow(abc))
  allcells$neighbor = rep("non.neighbor", nrow(allcells))
  
  for(cell.type in unique(temp3$cell)){
    abc[[cell.type]] = rep(ifelse(abc$cell.cluster_manual == cell.type, cell.type, "Other_cell"))#Ref
    allcells[[cell.type]] = rep(ifelse(allcells$cell.cluster_manual  == cell.type, cell.type, "Other_cell"))#Ref
  }
  
  fish.dat = rbind(abc, allcells)
  
  fisher.results = c()
  for(cell.type in unique(temp3$cell)){
    test = fisher.test(fish.dat$neighbor, fish.dat[[cell.type]], alternative = "greater")
    res.df = data.frame(cell = cell.type, p.val = test$p.value, odds.ratio = test$estimate)
    fisher.results = rbind(fisher.results, res.df)
  }
  results.df = left_join(temp3, fisher.results, by = "cell")
  results.list[[query.cell.type]] = results.df
}
  bootstrapped.res.list[[x]] = results.list
  #rm(temp1, temp2, temp3, allcells, abc, results.df, results.list, fisher.results, fish.dat, res.df, test)
})

#Process bootstrapping results-------------------------------
names(bootstrapped.res.list) = paste0("perm", seq(1:length(bootstrapped.res.list)))
bootstrapped.results = list()
bootstrapped.results = lapply(names(bootstrapped.res.list), function(x){
  bootstrapped.results[[x]] = bind_rows(bootstrapped.res.list[[x]], .id = "query.cell")
})

names(bootstrapped.results) = names(bootstrapped.res.list)
bootstrapped.results = bind_rows(bootstrapped.results, .id = "permutation")
temp4 = list()
temp5 = list()
for(query.cell.type in query.cell.types){
  temp4[[query.cell.type]] = bootstrapped.results %>% filter(query.cell == query.cell.type)
  query.result = results.measured %>% filter(query_cell == query.cell.type) # Get actual measured results of neighborhood composition for given query cell type. "results" here is from code before bootstrapping
  temp5[[query.cell.type]] = temp4[[query.cell.type]] %>% group_by(cell) %>%
    summarize(Count.neighbors.average.permuted = mean(Count.neighbors, na.rm=TRUE),
              Freq.neighbors.average.permuted = mean(Freq.neighbors, na.rm=TRUE)) %>% as.data.frame() #get average count and frequency of cells in permuted neighborhoods
  temp5[[query.cell.type]] = left_join(query.result, temp5[[query.cell.type]], by = "cell") #join actual and permuted data together
  temp5[[query.cell.type]]$avg.log2FC.enrichment.relative.to.permuted = log2(temp5[[query.cell.type]]$Freq.neighbors / temp5[[query.cell.type]]$Freq.neighbors.average.permuted) #calculate log2FC enrichment of real neighborhood relative to average permuted neighborhood
  
  #count number of times the actual measured frequency in neighborhood is larger than permuted neighborhood for each neighboring cell type
  def.res = c()
  for(neighbor.cell.type in temp5[[query.cell.type]]$cell){
    def = temp4[[query.cell.type]] %>% filter(cell == neighbor.cell.type, Freq.neighbors >= query.result[query.result$cell==neighbor.cell.type, "Freq.neighbors"]) %>% 
      nrow()
    def = data.frame(query_cell = query.cell.type,
                         cell = neighbor.cell.type,
                         times.larger.frequency.in.permuted = def)
    def.res = rbind(def.res, def)
  }
  temp5[[query.cell.type]] = left_join(temp5[[query.cell.type]], def.res, by = c("query_cell","cell")) #paste this count to previous result
  #calculate p-value, by dividing number of times permuted frequencies are larger by number of permutations
  #temp5 is the final resulting list of dataframes
  temp5[[query.cell.type]]$p.val.permuted = temp5[[query.cell.type]]$times.larger.frequency.in.permuted / nperm
  
  #Optionally filter out low-occuring cell types where results will be unreliable
  temp5[[query.cell.type]] = temp5[[query.cell.type]] %>% filter(Count.allcells >= 5 & Count.neighbors >= 3)
}

results.all = bind_rows(temp5, .id = "query_cell")


results.all$p_log10.permuted = -log10(results.all$p.val.permuted)
write.table(results.all, file = "neighbors.lesion.clusterlevel.txt", sep = "\t", quote = F, row.names = F)


