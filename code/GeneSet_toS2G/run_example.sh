geneset_dir=Test_GeneSets
bed_dir=Test_Beds
geneset=geneset_example.txt

#bash geneset_to_bed.sh $geneset_dir $bed_dir $geneset

bash geneset_to_bed.sh $geneset_dir $bed_dir 

#geneset=geneset_example.txt
geneset=geneset_example
bash clean_bed.sh $bed_dir/$geneset

bimfile_path=../../processed_data/BIMS
annot_path=Test_Annots
bash bed_to_annot.sh $bed_dir $bimfile_path $annot_path $geneset

#baseline_cell=/n/groups/price/kushal/Enhancer_MasterReg/data/ANNOTATIONS/Baselines

baseline_cell=[BaselineLD annotations Directory; add from alkesgrp to ../../processed_data]
bash geneset_coding_tss_promoter.sh $annot_path $geneset $baseline_cell

 
