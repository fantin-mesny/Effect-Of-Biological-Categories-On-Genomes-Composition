# Effect of lifestyles on the genomic composition in gene families

This method was developed and optimized by Fantin Mesny and Nathan Vannier (Max Planck Institute For Plant Breeding Research - Cologne, Germany).

## About

When doing comparative genomics, approches like effectors annotation and orthology prediction often result in gene-count tables. For instance, orthology prediction with a tool like OrthoFinder results in an *orthogroups* table where each row correspond to one genome, and each column to the number of genes contained in an orthogroup.

To assess the link between the composition in genes and biological categories of interest - that we will call here *lifestyles* - we have developed a method taking into account the phylogenetic relationships between organisms. This method allows both to test if there is an effect of lifestyles on the composition in genes, and if there is one, to compare the content in genes of the different lifestyles.

## Methodology: main steps

- **Convert the phylogeny into continuous variables**
Pairwise phylogenetic distances between samples are extracted from the input phylogenetic tree (using the function tree.distance() from Biopython's Phylo module), then formatted as a distance matrix.
This matrix is used to build a Principal Component Analysis (PCA). PC1 and PC2 coordinates can then be used as continuous variables in statistical models 
(NB: Our script only use the two first PCs. It can be modified to include more. According to our experience, two PCs are generally sufficient to cover >80% variance of the phylogeny.)


- **Build a distance matrix on the compositions in gene families**
Jaccard distances are calculated for each pair of organisms in the dataset, using the function vegdist(method="jaccard") from *Vegan* R package.
The input of this vegdist function is a simple dataframe with each row being one organism, and each column being one gene family.

- **Testing the effects of phylogeny and lifestyles on the compositions in gene families**
PERMANOVA can be used to assess the effects of any variable on a distance matrix.
To assess both the effects of phylogeny and lifestyle categories on the content in gene families, the function *adonis2* (from *Vegan* R package) is run using the following model: distMatrix~phyloPC1+phyloPC2+Lifestyle.
It is important to know that PERMANOVA is sequential. The test explains the maximum variance by the first factor, then the variance that is left by the second one, etc.

    The output to this test can be seen in file 'permanova.txt'
    Important results of this test are the p-values for each of the factors, as well as their Rsquared.
    If the *lifestyle* factor significantly explains differences in genome contents, next step can be interpreted.


- **Comparison of lifestyles**
If lifestyle categories explain part of the differences in gene family compositions, they can be compared for these compositions.
Pairwise comparisons of lifestyles is carried out using the function *pairwise.perm.manova* (from R package *RVaidememoire*. This test returns a matrix of pairwise p-values. When a p-value is significant, it means the two categories are significantly different. It can therefore be used to tell which categories are similar.


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
- ```--colors``` (facultative): To specify colors to associate to lifestyles. ```--colors lifestyle1:red,lifestyle2:blue,lifestyle3:#00FF00```
