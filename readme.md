# Detecting the effect of biological categories on genomes composition

This method was developed and optimized by Fantin Mesny and Nathan Vannier (Max Planck Institute For Plant Breeding Research - Cologne, Germany).

## Used in

Hage et al. "Gene family expansions and transcriptome signatures uncover fungal adaptations to wood decay." *Environmental Microbiology* (2021).

Miyauchi et al. "Large-scale genome sequencing of mycorrhizal fungi provides insights into the early evolution of symbiotic traits." *Nature communications* 11, no. 1 (2020): 1-17.


## About

A common difficulty in comparative genomics is to correct for phylogenetic signals that are often confounding factors when trying to detect genomic markers of a biological process.

When doing comparative genomics, approaches like effectors annotation and orthology prediction often result in count tables. For instance, orthology prediction with tools like OrthoFinder results in an *orthogroups* table where each row corresponds to one genome, and each column to the number of genes contained in an orthogroup.

To assess the link between genomic contents (i.e. orthogroup counts, effector counts, etc) and biological categories of interest, we have developed a method taking into account the phylogenetic relationships between organisms. This method allows first to test if there is a detectable effect of the biological categories of interest on the genomes content in genes. If there is one, its second step consists in comparing the genomic contents of the different categories. 

While this approach can be used for a wide variety of biological categories (virtually any category that is susceptible to affect genomic composition in genes), we used it to detect the effect of eukaryotic organisms’ lifestyle on genomic contents. 

## Methodology: main steps

- **Convert the phylogeny into continuous variables**
Pairwise phylogenetic distances between organisms are extracted from the input phylogenetic tree (using the function tree.distance() from Biopython's Phylo module), then formatted as a distance matrix.
This matrix is used to build a Principal Component Analysis (PCA). The coordinates of organisms on the PC1 and PC2 axes can then be used as continuous variables in statistical models 
(NB: Our script only uses the two first PCs. This can be modified to include more principal components. According to our experience, two PCs are generally sufficient to cover >80% of variance attributed to phylogeny.)


- **Build a distance matrix on the genomes composition in gene families**
Jaccard distances are calculated for each pair of organisms in the dataset, using the function vegdist(method="jaccard") from *Vegan* R package.
The input of this vegdist function is a simple dataframe with organisms as rows, and gene families as columns.

- **Testing the effects of phylogeny and lifestyles on the compositions in gene families**
PERMANOVA can be used to assess the effects of any variable on a distance matrix.
To assess both the effects of phylogeny and lifestyle categories on the content in gene families, the function *adonis2* (from *Vegan* R package) is run using the following model: distMatrix~phyloPC1+phyloPC2+Lifestyle. The PERMANOVA model explains the variance of the distance matrix in a sequential manner. The PERMANOVA tries to explain the variance with the first factor, then it tries to explain the variance that is left unexplained with the second factor, etc.


   The output of the test can be seen in the file *permanova.txt*. To interpret the results, first look at the *p*-values for each of the factors. If these are significant, have a look at the percentage of variance explained by each factor (R²). If the *lifestyle* factor significantly explains differences in genome contents, the next step can be interpreted, otherwise stop here.

- **Comparison of lifestyles**
If lifestyle categories explain part of the variance in gene families in the genomes, they can be compared to each other. Pairwise comparisons of lifestyles are carried out using the function *pairwise.perm.manova* (from R package *RVaidememoire*). This test returns a matrix of pairwise p-values. When a p-value is significant, it means the two categories have significantly different genome contents. The results can therefore be used to tell which categories are different.

## Software

The Python script in this repository allows the rapid execution of our method. Graphical outputs are returned, allowing a quick analysis of the results.

```bash
python run.py -i example/data.csv -t example/tree.nwk -o /output/directory -colors lifestyle1:blue,lifestyle2:green
```

### Output files

##### Graphical output files
- **pca.pdf**: Principal Component Analysis showing the genome compositions in gene families
- **pvalMatrix.pdf**: Graphical representation of lifestyle pairwise differences (red if significantly different)
- **network.pdf**: Simple network where lifestyles are connected if their composition in gene families is similar.

##### Other output files
- **distMatrix.csv**: Jaccard distances calculated on the genome compositions in gene families
- **permanova.txt**: output of the PERMANOVA test, containing the *Rsquared* and *p* values
- **pairwiseComparisons.csv**: output of the post-hoc testing, showing the *p*-values of pairwise lifestyle comparisons


### Installation


The following dependencies must be installed to run the ```run.py``` script.

- Python 3.x
- R (executable in command line using ```Rscript```)
- Python libraries: Pandas, Seaborn, Biopython, Sklearn, Networkx
- R packages: Vegan, RVAideMemoire

The ```run.py``` script can be downloaded from this repository.
To check it is able to run correctly, the example files can be used.

```bash
python run.py -i example/data.csv -t example/tree.nwk -o /output/directory 
```
NB: Example files were randomly generated and are not biologically relevant. Thus, results generated upon the execution of the command line above give negative results.


### Execution

##### Arguments

- ```--i```: CSV dataframe containing gene counts. One genome per row (names in column *genome*), one gene family per column. One column must be called *lifestyle*, and should contain the categories of interest. See ```example/data.csv```
- ```--t```: Phylogenetic tree (Newick format), with leaf name matching the names in column *genome* of the input dataframe (```--i```)
- ```--o```: Output/working directory
- ```--colors``` (facultative): To specify colors to associate to lifestyles. 
    For instance: ```--colors lifestyle1:red,lifestyle2:blue,lifestyle3:#00FF00```
