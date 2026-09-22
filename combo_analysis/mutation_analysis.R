####
#
# Code to Analyze Mutational Data and Relationship to Synergy
#
###
library(ggpubr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(glmnet)

rename.ccls <- read.csv("./data/CL_data/ccl_map.csv",sep=";")

colnames(rename.ccls) <- c("idSample","model_name")


mutations <- read.csv("./data/CL_data/mutations_all_20250318.csv")
mutations <- mutations %>% select(gene_symbol, model_id)
model.info <- read.csv("./data/CL_data/model_list_20250423.csv")
model.info <- model.info %>% select(model_id, model_name) %>% 
	filter(model_name %in% rename.ccls$model_name)


scores = read.csv("./results/median_synergy_scores.csv",row.names=1)



model.info <- model.info %>% merge(rename.ccls, by = "model_name") %>% 
			merge(mutations, by = "model_id")


binarized_mutation = c()

all.genes <- unique(model.info$gene_symbol)
n.genes <- length(all.genes)
ccls <- unique(model.info$idSample)

for(ccl in ccls){
	 
	ccl.mutations <-  model.info %>% filter(idSample==ccl) %>% pull(gene_symbol)
	# print(ccl.mutations)
	row <- rep(0,n.genes)
	# print(row)
	row[match(ccl.mutations,all.genes)]=1
	# res <- c(ccl,row)

	binarized_mutation <- rbind(binarized_mutation,row)
}
print(binarized_mutation)
colnames(binarized_mutation) <- all.genes
rownames(binarized_mutation) <- ccls
binarized_mutation <- as.data.frame(binarized_mutation)
print(binarized_mutation)
keep.genes <- names(which(colSums(binarized_mutation)>=20))
df <- as.data.frame(colSums(binarized_mutation)) %>%tibble::rownames_to_column("Gene Symbol")
colnames(df)[2] <- "Num. CLs w/ mut."
write.csv(df,"./results/all_mutation_count.csv")
binarized_mutation
score.mut <- binarized_mutation %>% tibble::rownames_to_column("idSample")%>%
	select(c("idSample",keep.genes))%>%
	merge(scores, on= "idSample")

syn.cor <- as.data.frame(cor(score.mut[,c(keep.genes,"Bliss","ZIP")]))%>%
	tibble::rownames_to_column("Gene Symbol")


labeled.cor <- merge(syn.cor[c("Gene Symbol","Bliss","ZIP")],df,by="Gene Symbol")

z <- max(labeled.cor[,c("Bliss","ZIP")])

write.csv(labeled.cor,"./results/all_mutation_correlation.csv")


#####
#
#
# Boxplots
#
########
mutation_status = c()
for(ccl in unique(model.info$idSample)){
	temp <- model.info %>% filter(idSample==ccl)
	pik_status = ifelse("PIK3CA" %in% temp$gene_symbol,"MUT","WT")
	p53_status = ifelse("TP53" %in% temp$gene_symbol,"MUT","WT")
	both_status = ifelse("TP53" %in% temp$gene_symbol & "PIK3CA" %in% temp$gene_symbol,
		"MUT","WT")
	res = c(ccl, pik_status,p53_status,both_status)
	mutation_status <- rbind(mutation_status,res)
}
colnames(mutation_status)<- c("idSample","PIK3CA_mut", "TP53_mut", "BOTH_mut")

mutation_status = as.data.frame(mutation_status)
rownames(mutation_status) <-NULL
scores = read.csv("./results/median_synergy_scores.csv",row.names=1)
print(scores)
ann_synergy <- mutation_status %>% merge(scores, by = "idSample")


# mutation_status[c("PIK3CA_mut", "TP53_mut", "BOTH_mut")] = sapply(mutation_status["PIK3CA_mut", "TP53_mut", "BOTH_mut"],numeric)
my_comparisons =  list( c("MUT","WT"))




mut = "PIK3CA_mut"

for(ref.model in c("ZIP", "Bliss")){
	for(mut in c("PIK3CA_mut", "TP53_mut","BOTH_mut")){
		mut.gene <- strsplit(mut,"_")[[1]][1]
		title.string <- paste0(ref.model," Synergy vs. ", mut.gene," Mutation Status")
		
		# stop()

		fname <- file.path("./results",paste0(mut.gene,"_",ref.model,".pdf"))
		pdf(fname)
		print(
			ggplot(ann_synergy,mapping = aes(x=.data[[mut]],y=.data[[ref.model]],fill=.data[[mut]])) +
 			 	geom_boxplot(outlier.shape = NA) + 
 				geom_hline(yintercept = 0,lty=2,col="red") +
  				geom_jitter(width = 0.2) +
  				ggtitle(title.string) +
 				stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 7.4) +
  # scale_fill_manual(values=PIK3CA_mut) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x =element_text(size=14,face="bold"))

)
dev.off()

	}
}

