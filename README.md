# GSSG : Combine SNP-gene linking strategy with Gene Sets 

A pipeline to perform the following steps. 

- Generate gene programs for analysis
- Link SNPs to genes in a gene program to generate a genome-wide annotation
- Perform disease heritability enrichment analysis using S-LDSC


## Citation

If you use the data or code from this repository, please cite these following papers

<ins> **Blood enhancer-regulated and master-regulator gene programs** </ins>
Dey, K.K., Gazal, S., van de Geijn, B., Kim, S.S., Nasser, J., Engreitz, J.M. and Price, A.L., 2022. SNP-to-gene linking strategies reveal contributions of enhancer-related and candidate master-regulator genes to autoimmune disease. Cell Genomics, 2(7), p.100145.

<ins> **Cell-type gene programs (sc-linker)** </ins>
Jagadeesh, K.J.\*, Dey, K.K.\* et al bioRxiv. 2021. Identifying disease-critical cell types and cellular processes across the human body by integration of single-cell profiles and human genetics.[Link](https://www.biorxiv.org/content/10.1101/2021.03.19.436212v2)

If you use the ABC S2G strategies please cite Nasser, J., Engreitz, J. et al. unpublished data. 2020. and [Fulco et al, 2019](https://www.nature.com/articles/s41588-019-0538-0). If you use PC-HiC S2G strategies, please cite [Javierre et al 2016 Cell](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123897/). If you use ATAC (Yoshida), cite [Yoshida et al 2019 Cell](https://www.cell.com/cell/pdf/S0092-8674(18)31650-7.pdf). If you use Roadmap S2G links, cite the references [here](https://ernstlab.biolchem.ucla.edu/roadmaplinking/). 


## Code Tutorial

### STEP 1: Generate a Gene Program

An example of how a gene program looks  ( see demo ~/code/GeneSet_toS2G/Test_GeneSets/geneset_example.txt )

```head -n 5 geneset_example.txt ```

```
RITA1	0
FAM225B	0.00320672144802546
TBC1D23	0.0716823315695179
SIRPG	0
L3MBTL4-AS1	0
```

It is a 2 column file with the first column being the names of genes (HGNC names). If the gene names are in Emntrez or
Ensembl, please convert them first to HGNC. 


How to generate gene scores related to our papers?

`code/calc_gene_scores` : Directory with codes to generate gene programs 

  - `calc_gene_scores_enhancers.R`: Generate enhancer-regulated gene scores in blood (Dey et al 2022 Cell Genomics)
  
  - `calc_gene_scores_masterreg.R`: Generate candidate master-regulator gene scores in blood (Dey et al 2022 Cell Genomics)
  
  - `process_sclinker_output.R`: Generate cell type programs from Alzheimers disease progression (sclinker pipeline)
  
`code/calc_PPI_scores`: Directory with codes to run PPI-informed gene programs 

  - `ppi_RWR.R`: Gene Programs using Random Walk with Restart on a general edgelist graph, see `example.R`
  
  - `ppi_string_RWR.R`: Gene Programs using Random Walk with Restart on the STRING PPI network. 

See `genesets` for all gene programs corresponding to Dey et al 2022 Cell Genomics and `sclinker_genescores` for all gene programs 
corresponding to Jagadeesh\*, Dey\* et al (sc-linker). 

*Illustration code for sclinker gene program generation*

```{r}
cd GSSG
Rscript process_sclinker_output.R --inputcell sclinker --prefix Alzheimers --outcell sclinker_genescores
```



### STEP 2: Gene set/program to SNP annotation

For this, please download all contents from the directories below

`https://alkesgroup.broadinstitute.org/LDSCORE/Dey_Enhancer_MasterReg/processed_data/`  and 
`https://alkesgroup.broadinstitute.org/LDSCORE/Jagadeesh_Dey_sclinker/extras/`


and save them into the `processed_data` directory. 

For batch download, please use the `gsutil` command (https://cloud.google.com/storage/docs/gsutil), 
and use the commands here (thanks to Krzysztof Polanski for pointing it out)

` gsutil cp -r gs://broad-alkesgroup-public/LDSCORE/Dey_Enhancer_MasterReg/processed_data .  `

` gsutil cp -r gs://broad-alkesgroup-public/LDSCORE/Jagadeesh_Dey_sclinker/extras .  `



`code/GeneSet_toS2G` :  Directory with codes to combine gene set with S2G strategy to create annotation

We first run the script that converts the gene set file to probabilistic .bed format files (bed-graph) for different S2G strategies.

First go to the folder

```
cd ~/GSSG/code/GeneSet_toS2G

```

### STEP 2A: Gene sets to bedgraph format files using region-gene linking


For **blood specific data (Dey et al 2020)** with 10 different S2G strategies, run the following script

```
geneset_dir=Test_GeneSets
geneset=geneset_example
bed_dir=Test_Beds
bash geneset_to_bed.sh $geneset_dir $bed_dir $geneset
```

For **sc-linker (Jagadeesh, Dey et al 2021)** with Roadmap-U-ABC enhancer-gene strategy for a tissue-specific enhancer,
please run the following script

```
enhancer_tissue="BLD" (for blood, can also be "BRN", "GI", "LNG", "LIV", "KID", "SKIN",
                      "FAT", "HRT" for brain, intestine, lung, liver, kidney, skin, adipose and heart)
bash geneset_to_bed_sclinker.sh $geneset_dir $bed_dir $geneset $enhancer_tissue
```

*Illustration code for sclinker gene program to bedgraph*

```{r}
Rscript geneset_to_bed_sclinker.R --genescore_dir sclinker_genescores/Alzheimers/ --bed_dir sclinker_beds/Alzheimers/ --anot_name Disease_Endothelial_L2 --enhancer BLD
```

### STEP 2B: Cleaning the bedgraph format files

We postprocess the created bedgraph files by removing overlapping intervals. The user needs BEDTOOLS and BEDOPS for the
following script. Check `clean_bed.sh` for details.

```
bash clean_bed.sh $bed_dir/$geneset 
bash clean_bed.sh $sclinker_beds/$celltype
```

*Illustration code for sclinker cleaning of bedgraph files*

```{r}
bash clean_bed.sh sclinker_beds/Disease_Endothelial_L2
```

### STEP 2C: From bedgraph to annotation files

Next we use a .bim file for a list of variants to annotate. We here use the .bim files widely used for S-LDSC analysis.

```
bimfile_path="~/GSSG/processed_data/BIMS"
annot_path=Test_Annots
bash bed_to_annot.sh $bed_dir $bimfile_path $annot_path $geneset
```

If you processed the bed files with **all blood S2G strategies Dey et al 2020)** using `geneset_to_bed.sh`, then for the `geneset_example.txt` file, the above codes will create a folder structure of the form 

```
- Test_Annots
  - geneset_example 
    - 100kb
    - 5kb
    - ABC
    - FinemapBloodeQTL (eQTL S2G in paper)
    - PCHiC
    - Roadmap_Enhancer
    - Yoshida (ATAC S2G in paper)
```

This includes 7 of 10 S2G strategies discussed in the manuscript. For the other strtegies, we look at TSS, Promoter and Coding regions for each gene set using the UCSC nomenclature and the already set-up baseline-LD model. For this, download the baseline-LD annotations from here `https://data.broadinstitute.org/alkesgroup/LDSCORE/`. Say these baseline-LD annotations are in a directory called `baseline_cell`.

Then run the following script to generate the S2G annotations for the remaining 3 strategies.

```
bash geneset_coding_tss_promoter.sh $annot_path $geneset $baseline_cell
```


If you processed the bed files with **tissue-specific enhancer-gene S2G strategies (Jagadeesh, Dey et al 2021)** using `geneset_to_bed_sclinker.sh`, then for the `geneset_example.txt` file, the above codes will create a folder structure 
of the form 

```
- Test_Annots
  - geneset_example 
    - ABC_Road_{enhancer_tissue}
    - 100kb
```

Example code

*Illustration code for sclinker generating annotation files*

```{r}
bash bed_to_annot.sh sclinker_beds processed_data/BIMS sclinker_annots Disease_Endothelial_L2
```



### Disease heritability signal of annotations using S-LDSC

   
1. For running the S-LDSC analysis, go to `GSSG/code/ldsc`. 

    - `create_ldscores.sh` : Create LD-scores from a reference panel
    
    - `run_ldsc_reg.sh`: Run regression model in S-LDSC
    
	  
For running these codes, you will need some files that you can save in two new directories -  LDSC_PATH and SUMSTATS_PATH.

1) Download the LDSC from git (https://github.com/bulik/ldsc/wiki/Partitioned-Heritability)
2) Download the baselineLD_v2.1 annotations from Broad webpage (https://alkesgroup.broadinstitute.org/LDSCORE/baselineLD_v2.1_annots/)
3) Download 1000G plink files (https://alkesgroup.broadinstitute.org/Variant_effects/1000G_EUR_Phase3_plink/)
4) HAPMAP3 SNPs (https://alkesgroup.broadinstitute.org/LDSCORE/hapmap3_snps.tgz)
5) Download the weights file for 1000G (https://alkesgroup.broadinstitute.org/LDSCORE/weights_hm3_no_hla.tgz)
5) Download the baseline frq file for 1000G available (https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_frq.tgz)
6) Download all sumstats you want to analyze and put them in SUMSTATS_PATH (https://alkesgroup.broadinstitute.org/LDSCORE/all_sumstats/)


You can postprocess the LDSC output across traits using following codes

	  - `ldsc_postprocess_taustar.R`: calculate tau-star metric for S2G annotations of a gene set
	  
	  - `ldsc_postprocess_enrichment.R` : calculate the enrichment metric for S2G annotations of a gene set
	  
	  - `ldsc_postprocess_joint_taustar.R` : Joint S-LDSC analysis tau-star
	  
	  - `ldsc_postprocess_joint_enrichment.R` : Joint S-LDSC analysis enrichment 
	  
	  - `ldsc_postprocess_combined_taustar.R` : Compute combined tau-star metric for the annotations 
    
## Calculating E-scores

In sc-linker, we are proposing E-score as a metric to assess disease association of a gene program against all protein-coding genes.
The E-score is defined as

**E-score (program, trait) := S-LDSC Enrichment ( SNP annotation from Gene Program x Enhancer-gene strategy ) - 
                                S-LDSC Enrichment (SNP annotation from All protein-coding genes x Enhancer-gene )**

We compute the E-score (program, trait) for each program across a large collection of traits (60).
We define an estimate of standard error of E-score as

**sE-score (program, trait) := sqrt(S1^2 + S2^2 + 2cor(SNP annotation from Gene Program x Enhancer-gene strategy, 
                                                   SNP annotation from All protein-coding genes x Enhancer-gene)*S1*S2)**

where

*S1 = S-LDSC s.e. Enrichment (SNP annotation from Gene program x Enhancer-gene )* 

*S2 = S-LDSC s.e. Enrichment (SNP annotation from All protein-coding genes x Enhancer-gene )* 

For small gene programs, you can alternatively use

**sE-score (program, trait) := sqrt(S1^2 + S2^2)**

We then define

**zE-score :=. E-score/sE-score**

The distribution of zE-score is not always N(0, 1) across traits as traits have inherent correlation and pleiotropy, but you can use z-score distribution across many traits to get a p-value. 
		

## Gene Sets

Probabilistic gene sets (.txt files with 2 columns, first column gene symbol HGNC, the second symbol
the probabilitic membership of the gene in the gene set. Equals to 1 for binary gene sets. Genes with 0 weight
can be ignored)

- Enhancer-driven gene scores: ABC9_G_top10_0_015.txt (ABC-distal), EDS_Binary_top10.txt (EDS-binary),
                              eQTL_CTS_prob.txt (eQTL-CTS), Expecto_MVP.txt (ExPecto-MVP), HOMOD_prob.txt (ATAC-distal)
                              PCHiC_binary.txt (PC-HiC-binary), SEG_GTEx_top10.txt (SEG-blood)  [`genesets/`]
                              
- Master-regulator gene scores:  master_regulator.txt (Trans-master), TF_genes_curated.txt (TF-genes)  [`genesets/`]

- PPI-informed gene scores:  PPI_Enhancer.txt (PPI-enhancer), PPI_Master.txt (PPI-master),
                                    PPI_All.txt (PPI-all), PPI_control.txt (PPI-control)  [`genesets/`]
                                    
- Additional gene scores: pLI_genes2.txt (pLI), RegNet_Enhancer.txt (RegNet-Enhancer), Trans_Reg_genes.txt (Trans-regulated)  [`genesets/`]

- Example sc-linker gene scores: Alzheimeers disease progression example - all cell types (`sclinker_genescores/Alzheimers`)

## Annotations

All gene set-S2G annotations studied in the companion papers can be found at `https://alkesgroup.broadinstitute.org/LDSCORE/Dey_Enhancer_MasterReg/` (Dey et al Cell Genomics 2022) and `https://alkesgroup.broadinstitute.org/LDSCORE/Jagadeesh_Dey_sclinker/` (Jagadeesh\*, Dey\* et al, to appear, Nat Genet). 


## How to use these annotations?

```
ANNOT FILE header (*.annot):

CHR -- chromosome
BP -- physical position (base pairs)
SNP -- SNP identifier (rs number)
CM -- genetic position (centimorgans)
all additional columns -- Annotations
```

NOTE: Although one would expect the genetic position to be non-negative for all 1000G SNPs, we have checked that
in fact the genetic position is negative for a handful of 1000G SNPs that have a physical position that is smaller
than the smallest physical position in the genetic map. The genetic positions were obtained by running PLINK on
the Oxford genetic map (http://www.shapeit.fr/files/genetic_map_b37.tar.gz).

MORE DETAIL CAN BE OBTAINED FROM https://github.com/bulik/ldsc/wiki/LD-File-Formats


```
LD SCORE FILE header (*.l2.ldscore):

CHR -- chromosome
BP -- physical position (base pairs)
SNP -- SNP identifier (rs number)
all additional columns -- LD Scores

```

## Contact 

In case of any questions, please open an issue or send an email to me at `kdey@hsph.harvard.edu`.









