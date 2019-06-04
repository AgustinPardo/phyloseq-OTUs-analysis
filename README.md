# phyloseq-OTUs-app

### Data requirements

This workflow is develop to use OTUs and taxa tables in RDS format. 
You could also use OTUs table from txt/csv and provide your own sample_data and tax_table, but you should fits a tittle the first steps of the pipeline.

### Starting point

An example data set that cames straightforwar from the pipeline [**BacSixTnS**](https://github.com/AgustinPardo/BacSixTnS) is used. This data set is an OTU and taxa table from a 16S sequencing samples of diferents areas and seassons of a farming area.

### Things that you could do
* Creation phyloseq object from OTUs, taxa table in a RDS and csv format.
* Add to phyloseq object metadata from the samples to filter and play with it
* Counts the reads per sample
* Filter the samples per reads
* Analyse available taxonomic ranks in the dataset
* OTUs singletons removing
* Export new parsed OTUs and taxa table
* Plot rarefaction curves
![rare](https://github.com/AgustinPardo/phyloseq-OTUs-app/blob/master/rarefaction_curve.png)
* Subset selecion by the metada (e.g: Location, Season)

### Bar plots
* Create Bar plots of the diferetns taxonomic ranks, "Kingdom","Phylum","Class","Order","Family","Genus","Species".
![BP_all](https://github.com/AgustinPardo/phyloseq-OTUs-app/blob/master/barplot_phylum_all.png)
* Select especific taxa of a rank. For example: in the case of Phylum, you could select "Acidobacteria","Proteobacteria","Bacteroidetes","Chloroflexi", "Cyanobacteria","Firmicutes","Nitrospirae".
![BP_selected](https://github.com/AgustinPardo/phyloseq-OTUs-app/blob/master/barplot_phylum_selected.png)
### Ordination plots
* Coming soon...
![OP_all](https://github.com/AgustinPardo/phyloseq-OTUs-app/blob/master/ordination_plot_all.png)

### Data input examples
https://drive.google.com/drive/folders/1Zd_ErYQro95F-oArd_CoXHLOP9eQTZiG?usp=sharing
