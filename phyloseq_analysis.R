library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")

setwd("C:/Users/apardo/Documents/metagenomica/scripts-github/DADA2_out/")

## Importing OTUs and taxa tables

# Importing an exisitng OTUs table from txt 
########################
otu_facu <- read.csv("EX_otutable.txt", row.names=1)
OTU = otu_table(otu_facu , taxa_are_rows = TRUE)
head(OTU)
physeq = phyloseq(OTU)
########################

# Importing an exisitng otu table from RDS 
# Open files
taxa <- readRDS("tax_example.rds")
seqtab.nochim <- readRDS("seqtab_nochim_example.rds")

# Creation of the phyloseq objetc and adding atributes

samples.out <- rownames(seqtab.nochim)

onlySeasson <- function(stringe){
  tmp=strsplit(stringe,"_")
  seasson=tmp[[1]][2]
  return(seasson)
}
onlyLocation <- function(stringe){
  tmp=strsplit(stringe,"_")
  name=tmp[[1]][1]
  return(name)
}

location <- sapply(samples.out, onlyLocation)
seasson <- sapply(samples.out, onlySeasson)

samdf <- data.frame(Location=location,Seasson=seasson)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
sample_data(ps)

# Selection samples by atributes
# ps_interest= subset_samples(ps, sample_data(ps)$Location %in% Seasson)
# ps_interest= subset_samples(ps, sample_data(ps)$Seasson %in% Location)
# sample_data(ps_interest)

# Filter samples by reads count (Example 5000 reads)
sort(sample_sums(ps))
ps_reads_filter <- prune_samples(sample_sums(ps)>=5000, ps)
sample_data(ps_reads_filter)

# Show available ranks in the dataset
rank_names(ps)
# Create table, number of features for each phyla
# The taxonomic rank could be changed to "Kingdom","Phylum","Class","Order","Family","Genus","Species"
(table(tax_table(ps)[, "Kingdom"], exclude = NULL))

# I select only the Bacterias Kingdom. 
# You could use subset_taxa to filter by any Taxonomic rank
ps0 = subset_taxa(ps_reads_filter, Kingdom == "Bacteria")
# mapping=sample_data(ps0)
# write.csv(mapping, "/home/agustin/laura_Raiger_amplicones/cecarData/cecarOut/mapping.txt")
# check that the selection was satisfactory
(table(tax_table(ps0)[, "Kingdom"], exclude = NULL))


## Export OTUs and taxa table to files

# Export the OTUs table and taxa with the sequences as anames in csv format
otuTable_original=t(otu_table(ps0))
write.csv(otuTable_original, "otu_table.csv")
taxaTable_original=tax_table(ps0)
write.csv(taxaTable_original, "taxa_table.csv")

# # Export the OTU table and taxa changing the name by SV# in csv format
ps0_SV=ps0
taxa_names(ps0_SV) <- paste0("SV", seq(ntaxa(ps0_SV)))
# Exporto otu y taxa table
otuTable=t(otu_table(ps0_SV))
write.csv(otuTable, "otu_table_SV.csv")
taxaTable=tax_table(ps0_SV)
write.csv(taxaTable, "taxa_table_SV.csv")

# Export a file with the refercen between the sequence and the SV name
sequences=taxa_names(ps0)
id=taxa_names(ps0_SV)
reference_id=cbind(id ,sequences)
write.csv(reference_id, "reference_taxa_otu.txt",row.names=FALSE)

## Singletons removing
ps0_rare_filter=prune_taxa(taxa_sums(ps0) > 1, ps0)

## Rarefaction curves
source("C:/Users/apardo/Documents/metagenomica/ggrare.R", local = TRUE)
p <- ggrare(ps0_rare_filter, step = 5000, color = "Location", label = "Location", se = FALSE)
p + theme(legend.position="none") +  ggtitle("")
#p + facet_wrap(~Seasson)+ theme(legend.position="none")

## Subset selection
# Select by Seasson atributte
subset_seasson=subset_samples(ps0_rare_filter,sample_data(ps0_rare_filter)$Seasson == "Fall2016")
sample_names(subset_seasson)
# Select by Location atributte
subset_location=subset_samples(ps0_rare_filter,sample_data(ps0_rare_filter)$Location == "Y1A")
sample_names(subset_location)

# Example case
# I will continue to work with Fall2016 seasson subset
subset_ps0=subset_seasson

# Remove NA from each taxonomic rank
subset_ps0_noNAPhylum <- subset_taxa(subset_ps0, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
subset_ps0_noNAClass <- subset_taxa(subset_ps0_noNAPhylum, !is.na(Class) & !Class %in% c("", "uncharacterized"))
subset_ps0_noNAOrder <- subset_taxa(subset_ps0_noNAClass, !is.na(Order) & !Order %in% c("", "uncharacterized"))
subset_ps0_noNAFamily <- subset_taxa(subset_ps0_noNAOrder, !is.na(Family) & !Family %in% c("", "uncharacterized"))
subset_ps0_noNAGenus <- subset_taxa(subset_ps0_noNAFamily , !is.na(Genus) & !Genus %in% c("", "uncharacterized"))
# Check the NAs remotion
table(tax_table(subset_ps0_noNAPhylum)[, "Phylum"], exclude = NULL)

#Export otu and taxa table of the subset
GPr  = transform_sample_counts(subset_ps0, function(x) x / sum(x) )
otuTable=otu_table(GPr)
otuTable_tranpose=t(otuTable)
write.csv(otuTable_tranpose, "Fall2016OT.csv") 
taxaTable=tax_table(GPfr)
write.csv(taxaTable, "Fall2016TT.csv")

####################
#     BAR PLOTs    #
####################

# Taxonomic filtering. Filter by abundance of OTU in the samples.
# Hint: When you create a subset, the remaining OTUs are there with no read counts.
# That is why if you apply very low filters you still have a reduction, because 
# the abundance of that OTUs is 0.

GPr_Phylum  = transform_sample_counts(subset_ps0_noNAPhylum, function(x) x / sum(x) )
# Filter by abundancy: 10e-100 (play with the number)
GPfr_Phylum = filter_taxa(GPr_Phylum, function(x) mean(x) > 10e-100, TRUE)
# Black lines in the milde of the bars. Solved in the lines below.

GPr_Class  = transform_sample_counts(subset_ps0_noNAClass, function(x) x / sum(x) )
GPfr_Class = filter_taxa(GPr_Class, function(x) mean(x) > 10e-200, TRUE)
GPr_Order  = transform_sample_counts(subset_ps0_noNAOrder, function(x) x / sum(x) )
GPfr_Order = filter_taxa(GPr_Order, function(x) mean(x) > 10e-50, TRUE)
GPr_Family = transform_sample_counts(subset_ps0_noNAFamily, function(x) x / sum(x) )
GPfr_Family= filter_taxa(GPr_Family, function(x) mean(x) > 0.0008, TRUE)
GPr_Genus = transform_sample_counts(subset_ps0_noNAGenus, function(x) x / sum(x) )
GPfr_Genus= filter_taxa(GPr_Genus, function(x) mean(x) > 0.002, TRUE)

# Hint: Change large names that dont fit in the frame of the plots
# Is up to you to add more if sentences.
temp = tax_table(GPfr_Genus)
len=dim(a[,6])
total=len[1]
vec=c(1:total)
for (i in vec) {
  if(a[,6][i]=="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"){
    a[,6][i]="ANPR"
    print(i)
  }
  if(a[,6][i]=="Burkholderia-Caballeronia-Paraburkholderia"){
    a[,6][i]="BCP"
    print(i)
  }
}

tax_table(GPfr_Genus) = temp
# Reecalculate the abundancy (avoiding that the bars don't reach the 100% percent)
GPfr_Genus= transform_sample_counts(GPfr_Genus, function(x){x / sum(x)})

## Plots
# Plot with out black lines (abundancy of the taxonomic lower leves inside the bars)
phylumGlommed_phylum = tax_glom(GPfr_Phylum, "Phylum")
plot_bar(phylumGlommed_phylum, fill ="Phylum") + theme(legend.position="none")

# Export the data
table(tax_table(phylumGlommed_phylum )[, "Phylum"])
otuTable=otu_table(phylumGlommed_phylum )
write.csv(otuTable, "Phylum_otuTable.csv") 
taxaTable=tax_table(phylumGlommed_phylum)
write.csv(taxaTable, "Phylum_taxa_table.csv") 

phylumGlommed_class = tax_glom(GPfr_Class, "Class")
plot_bar(phylumGlommed_class, fill ="Class") #+ theme(legend.position="none")
phylumGlommed_order = tax_glom(GPfr_Order , "Order")
plot_bar(phylumGlommed_order, fill ="Order") #+ theme(legend.position="none")
phylumGlommed_family = tax_glom(GPfr_Family, "Family")
plot_bar(phylumGlommed_family, fill ="Family") #+ theme(legend.position="none")

# Tener en cuenta de que si quedan samples vacias( por ul filtrado muy restrictivo,
# se rompe el siguiente comando)
phylumGlommed_genus = tax_glom(GPfr_Genus, "Genus")
plot_bar(phylumGlommed_genus, fill ="Genus")
# Si pasa lo comentado arriba hacer este grafico..
plot_bar(GPfr_Genus, fill ="Genus")

# If I want to select expecfic taxa(Example with Phylum):
########################################
taxa_target <- c("Acidobacteria","Proteobacteria","Bacteroidetes","Chloroflexi",
                 "Cyanobacteria","Firmicutes","Nitrospirae")
select_phylum_ps=subset_taxa(subset_ps0 , Phylum  %in% taxa_target)
GPr = transform_sample_counts(select_phylum_ps, function(x) x / sum(x) )
phylumGlommed_phylum = tax_glom(GPr, "Phylum")
plot_bar(phylumGlommed_phylum, fill ="Phylum") 

# Also, I can merge samples (using merge_samples). 
# In this case by Seasson, but it can be with Location too.
subset_ps0_merge=merge_samples(subset_ps0,"Seasson")
taxa_target <- c("Acidobacteria","Proteobacteria","Bacteroidetes","Chloroflexi",
                 "Cyanobacteria","Firmicutes","Nitrospirae")
select_phylum_ps=subset_taxa(subset_ps0_merge , Phylum  %in% taxa_target)
GPr = transform_sample_counts(select_phylum_ps, function(x) x / sum(x) )
phylumGlommed_phylum = tax_glom(GPr, "Phylum")
plot_bar(phylumGlommed_phylum, fill ="Phylum") 


####################
# Ordination plots #
####################

#   ORDINATION PLOT, TODOS   #
########################

sample_data(ps0_rare_filter)
GP.ord_NMDS <- ordinate(ps0_rare_filter , "PCoA", "bray")
p1 = plot_ordination(ps0_rare_filter , GP.ord_NMDS , type="samples", color="Location",shape = "Seasson")
p1 = p1  + geom_point(size=5) + ggtitle("")
p1


#temporadas=list("Spring2014","Fall2015","Fall2016","Spring2016", "Spring2017", "Winter2017")
#locaciones=list("Planta1","Planta4","Yucra1", "Yucra4", "Yucra6" , "Severino","Referencia1", "Referencia2")
temporadas=list(S14=spring2014,S16=spring2016,S17=spring2017,F15=fall2015,F16=fall2016,W17=winter2017)
locaciones=list(P1=planta_1,P4=planta_4,Y1=yucra_1,Y4=yucra_4,Y6=yucra_6,Ref1=referencia1,S=severino)
# Saque Ref2 porque solo esta en una temporada

# ** Todos ** #
# subset_ps0 = subset_samples(ps0, sample_names(ps0) %in% completeSubsets)
#Ordinate function NMDS
sample_data(ps0_RF)
GP.ord_NMDS <- ordinate(ps0_RF , "PCoA", "bray")
p1 = plot_ordination(ps0_RF , GP.ord_NMDS , type="samples", color="Location", shape = "Seasson")
p1 = p1  + geom_point(size=5) #+ ggtitle(subset_name)
p1

# ** Para temporadas ** #
for (temporada in temporadas) {
  subset_ps0 = subset_samples(ps0, sample_names(ps0) %in% spring2017)
  owner_name="bia"

  temporada_name=strsplit(temporada[1],"_")
  temporada_name=temporada_name[[1]][2]
  
  path="/home/agustin/laura_Raiger_amplicones/Analisis/web/ordinationPlots/temporada/"
  a <- paste(path,owner_name,sep = "")
  b <- paste(a,"", sep = "")
  
  #Ordinate function NMDS
  GP.ord_NMDS <- ordinate(subset_ps0 , "NMDS", "bray")
  #Ordination NMDS plot
  c <- paste(b,temporada_name,sep = "_")
  d <- paste(c,"OP_NMDS_area",sep = "_")
  e <- paste(d,"svg",sep = ".")
  svg(e)
  p1 = plot_ordination(subset_ps0 , GP.ord_NMDS , type="samples", color="Location")
  p1 = p1 + geom_polygon(aes(fill=Location)) + geom_point(size=5) #+ ggtitle(subset_name)

  print(p1)
  dev.off()
  
  #Sin rellenar areas
  c <- paste(b,temporada_name,sep = "_")
  d <- paste(c,"OP_NMDS",sep = "_")
  e <- paste(d,"svg",sep = ".")
  svg(e)
  p1 = plot_ordination(subset_ps0 , GP.ord_NMDS , type="samples", color="Location")
  p1 = p1  + geom_point(size=5) #+ ggtitle(subset_name)
  print(p1)
  dev.off()
  
  #Ordinate function PCoA
  GP.ord_PCoA <- ordinate(subset_ps0 , "PCoA", "bray")
  #Ordination PCoA plot
  c <- paste(b,temporada_name,sep = "_")
  d <- paste(c,"OP_PCoA_area",sep = "_")
  e <- paste(d,"svg",sep = ".")
  svg(e)
  p2 = plot_ordination(subset_ps0 , GP.ord_PCoA, type="samples", color="Location")
  p2 = p2 + geom_polygon(aes(fill=Location)) + geom_point(size=5) #+ ggtitle(subset_name)
  print(p2)
  dev.off()
  
  #Sin rellenar areas
  c <- paste(b,temporada_name,sep = "_")
  d <- paste(c,"OP_PCoA",sep = "_")
  e <- paste(d,"svg",sep = ".")
  svg(e)
  p2 = plot_ordination(subset_ps0 , GP.ord_PCoA  , type="samples", color="Location")
  p2 = p2  + geom_point(size=5) #+ ggtitle(subset_name)
  print(p2)
  dev.off()
  
}

# ** Para Locacion ** #
for (locacion in locaciones) {
  subset_ps0 = subset_samples(ps0, sample_names(ps0) %in% locacion)
  owner_name="bia"
  
  location_name=strsplit(locacion[1],"_")
  location_name=location_name[[1]][1]
  location_name=substr(location_name,0,nchar(location_name)-1)

  path="/home/agustin/laura_Raiger_amplicones/Analisis/web/ordinationPlots/locacion/"
  a <- paste(path,owner_name,sep = "")
  b <- paste(a,"", sep = "")
  
  #Ordinate function NMDS
  GP.ord_NMDS <- ordinate(subset_ps0 , "NMDS", "bray")

  #Ordination NMDS plot
  c <- paste(b,location_name,sep = "_")
  d <- paste(c,"OP_NMDS_area",sep = "_")
  e <- paste(d,"svg",sep = ".")
  svg(e)
  p1 = plot_ordination(subset_ps0 , GP.ord_NMDS , type="samples", color="Seasson")
  p1 = p1 + geom_polygon(aes(fill=Seasson)) + geom_point(size=5) #+ ggtitle(subset_name)
  p1
  print(p1)
  dev.off()
  
  #Sin rellenar areas
  c <- paste(b,location_name,sep = "_")
  d <- paste(c,"OP_NMDS",sep = "_")
  e <- paste(d,"svg",sep = ".")
  svg(e)
  p1 = plot_ordination(subset_ps0 , GP.ord_NMDS , type="samples", color="Seasson")
  p1 = p1  + geom_point(size=5) #+ ggtitle(subset_name)
  print(p1)
  dev.off()
  
  #Ordinate function PCoA
  GP.ord_PCoA <- ordinate(subset_ps0 , "PCoA", "bray")
  #Ordination PCoA plot
  c <- paste(b,location_name,sep = "_")
  d <- paste(c,"OP_PCoA_area",sep = "_")
  e <- paste(d,"svg",sep = ".")
  svg(e)
  p2 = plot_ordination(subset_ps0 , GP.ord_PCoA, type="samples", color="Seasson")
  p2 = p2 + geom_polygon(aes(fill=Seasson)) + geom_point(size=5) #+ ggtitle(subset_name)
  print(p2)
  dev.off()
  
  #Sin rellenar areas
  c <- paste(b,location_name,sep = "_")
  d <- paste(c,"OP_PCoA",sep = "_")
  e <- paste(d,"svg",sep = ".")
  svg(e)
  p2 = plot_ordination(subset_ps0 , GP.ord_PCoA  , type="samples", color="Seasson")
  p2 = p2  + geom_point(size=5) #+ ggtitle(subset_name)
  print(p2)
  dev.off()
}

### Todos los Ordination Methods Juntos ###
###Supported Ordination Methods###
library("plyr"); packageVersion("plyr")

dist = "bray"

# All methods
# "DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA"
# DPCoa necesita Arbol filogenetico

ord_meths = c( "DCA", "CCA", "RDA",  "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="SampleType")
},subset_ps0, dist)

names(plist) <- ord_meths

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=SampleType, fill=SampleType))
p = p + geom_point(size=4) + geom_polygon()
p = p + facet_wrap(~method, scales="free")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + scale_colour_brewer(type="qual", palette="Set1")+ ggtitle(subset_name)
p

############################
##Alpha diversity graphics##
############################
subset_ps0 = subset_samples(ps0, sample_data(ps0)$Seasson %in% "Spring2017")
subset_name="Spring2017"
#Prepare data
subset_ps0_preprare <- prune_species(speciesSums(subset_ps0) > 0, subset_ps0)
plot_richness(subset_ps0_preprare,x="Location", color = "Location")  + ggtitle(subset_name)
plot_richness(subset_ps0_preprare,x="SampleType",color = "SampleType") + ggtitle(subset_name)

plot_richness(subset_ps0_preprare,x="Location", color = "Location",measures=c("Chao1")) + geom_point(size=2)+theme(legend.position="none",axis.title.x=element_blank())
plot_richness(subset_ps0_preprare,x="Location",  color = "Location",measures=c("Shannon")) + geom_point(size=2)+theme(legend.position="none",axis.title.x=element_blank())

# Mecanizado #
temporadas=list(S14=spring2014,S16=spring2016,S17=spring2017,F15=fall2015,F16=fall2016,W17=winter2017)
locaciones=list(P1=planta_1,P4=planta_4,Y1=yucra_1,Y4=yucra_4,Y6=yucra_6,Ref1=referencia1,Ref2=referencia2,S=severino)
# Temporada #
for (temporada in temporadas) {
  subset_ps0 = subset_samples(ps0, sample_names(ps0) %in% spring2017)
  owner_name="bia"
  temporada_name=strsplit(temporada[1],"_")
  temporada_name=temporada_name[[1]][2]
  
  path="/home/agustin/laura_Raiger_amplicones/Analisis/web/alphaPlots/temporada/"
  a <- paste(path,owner_name,sep = "")
  b <- paste(a,temporada_name,sep = "_")
  c <- paste(b,"alpha_Chao1",sep = "_")
  d <- paste(c,"svg",sep = ".")
  subset_ps0_preprare <- prune_species(speciesSums(subset_ps0) > 0, subset_ps0)
  svg(d)
  p1=plot_richness(subset_ps0_preprare,x="Location", color = "Location",measures=c("Chao1")) + geom_point(size=2)+theme(legend.position="none",axis.title.x=element_blank())
  p1
  print(p1)
  dev.off()
  
  c <- paste(b,"alpha_Shannon",sep = "_")
  d <- paste(c,"svg",sep = ".")
  svg(d)
  p2=plot_richness(subset_ps0_preprare,x="Location",  color = "Location",measures=c("Shannon")) + geom_point(size=2)+theme(legend.position="none",axis.title.x=element_blank())
  print(p2)
  dev.off()
}
# Locacion #
# Resolver problemas de orden de nombres....
miOrden_nombres=c("Spring2014", "Fall2015",   "Spring2016", "Fall2016",   "Winter2017", "Spring2017")
for (locacion in locaciones) {
  subset_ps0 = subset_samples(ps0, sample_names(ps0) %in% planta_4)
  owner_name="bia"
  
  sample_data(subset_ps0)
  
  location_name=strsplit(locacion[1],"_")
  location_name=location_name[[1]][1]
  location_name=substr(location_name,0,nchar(location_name)-1)
  
  path="/home/agustin/laura_Raiger_amplicones/Analisis/web/alphaPlots/locacion/"
  a <- paste(path,owner_name,sep = "")
  b <- paste(a,location_name,sep = "_")
  c <- paste(b,"alpha_Chao1",sep = "_")
  d <- paste(c,"svg",sep = ".")
  subset_ps0_preprare <- prune_species(speciesSums(subset_ps0) > 0, subset_ps0)
  
  sample_data(subset_ps0_preprare)
  
  svg(d)
  p1=plot_richness(subset_ps0_preprare,x="Order", color = "Seasson",measures=c("Chao1")) + geom_point(size=2)+theme(legend.position="none",axis.title.x=element_blank())
  p1 + scale_x_discrete(labels = miOrden_nombres) 
  print(p1)
  dev.off()
  
  c <- paste(b,"alpha_Shannon",sep = "_")
  d <- paste(c,"svg",sep = ".")
  svg(d)
  p2=plot_richness(subset_ps0_preprare,x="Seasson",  color = "Seasson",measures=c("Shannon")) + geom_point(size=2)+theme(legend.position="none",axis.title.x=element_blank())
  print(p2)
  dev.off()
}

#Datos de estimacion Tabla
subset_ps0_seasson = subset_samples(ps0, sample_names(ps0) %in% spring2017) 
subset_ps0=merge_samples(subset_ps0_seasson ,"SampleType")
sample_names(subset_ps0)
alphaDiversity_table=estimate_richness(subset_ps0, split = TRUE, measures = NULL)
salida="/home/agustin/laura_Raiger_amplicones/Analisis/web/Spring17/alphaDiversity_tableSpring17.csv"
write.csv(alphaDiversity_table,salida )