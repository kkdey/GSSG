=====================================================
Steps for converting gene set (GS) to S2G annotations
=====================================================

- First step:
  Start with a gene score file. This is a standard .txt file (no column name, no row name), with gene name in the first column, and
  the score in the second column (it does not automatically work with Ensembl and Entrez ids, so please convert to gene names for the 
  Ensembl ids)
  
 - See an example gene score file 
- INPUT_FOLDER=/n/groups/price/kushal/TEST_GS2S2G; geneset_dir=$INPUT_FOLDER/Test_GeneSets;  head -n 5 $INPUT_FOLDER/Test_GeneSets/geneset_example.txt`

RITA1						   0
FAM225B						   0.00320672144802546
TBC1D23						   0.0716823315695179
SIRPG						   0
L3MBTL4-AS1					   0

- Next create a folder where you want to save the bed (graph) format files created from this gene score file
- bed_dir=$INPUT_FOLDER/Test_Beds`

- Run the first script to save the bed(graph) files for ABC, Roadmap, and ABC U Roadmap S2G strategies for the given gene set. If you have multiple
  gene set files in the same geneset_dir, then also the following code should work
- bash geneset_to_bed.sh $geneset_dir $bed_dir 'BLD'

This script will create a number of folders as number of gene sets in $bed_dir with bed files corresponding to ABC, Roadmap and ABC Roadmap union
- In this case, it will look like the following
- Test_Beds/
  - geneset_example/
    - ABCI_BLD.bed
    - ABC_Road_I_BLD.bed
    - Road_I_BLD.bed

- Now take for the gene set you are interested (in case of multiple gene sets), run the following cleaning operation that uses bedops to clean the 
  bedgraph files from previous step.
- geneset=geneset_example
- bash clean_bed.sh $bed_dir/$geneset

As a result of running this script, the bed files in the previous step will get modified, no extra files will be created.

- Now, you are all set to generate annotations from the cleaned bed (graph) files. First, create  folder where you want to save your annotations 
  data.
- annot_path=$INPUT_FOLDER/Test_Annots

- Then run the following code 
- bimfile_path=/n/groups/price/kushal/DATA/BIMS
- bash bed_to_annot.sh $bed_dir $bimfile_path $annot_path $geneset

- This script will create annotation files of the following form
- Test_Annots/
  - geneset_example
     - ABCI_BLD
     - ABC_Road_I_BLD
     - Road_I_BLD

Each of the sub-directories contains the S2G annotations for the gene set linked by the S2G strategy (ABC, ABC+Roadmap and Roadmap) respectively.

The full pipeline:


INPUT_FOLDER=/n/groups/price/kushal/TEST_GS2S2G
geneset_dir=$INPUT_FOLDER/Test_GeneSets
bed_dir=$INPUT_FOLDER/Test_Beds

bash geneset_to_bed.sh $geneset_dir $bed_dir 'BLD'

geneset=geneset_example
bash clean_bed.sh $bed_dir/$geneset

bimfile_path=/n/groups/price/kushal/DATA/BIMS
annot_path=$INPUT_FOLDER/Test_Annots
bash bed_to_annot.sh $bed_dir $bimfile_path $annot_path $geneset


 
