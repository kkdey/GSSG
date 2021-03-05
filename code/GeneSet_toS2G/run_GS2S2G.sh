INPUT_FOLDER=/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/code/GeneSet_toS2G
geneset_dir=$INPUT_FOLDER/Test_GeneSets
bed_dir=$INPUT_FOLDER/Test_Beds

bash geneset_to_bed.sh $geneset_dir $bed_dir 'BLD'

geneset=geneset_example
bash clean_bed.sh $bed_dir/$geneset

bimfile_path=/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/data/BIMS
annot_path=$INPUT_FOLDER/Test_Annots
bash bed_to_annot.sh $bed_dir $bimfile_path $annot_path $geneset


 
