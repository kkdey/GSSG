genescore_cell=/n/groups/price/kushal/singlecellLDSC/data/Gene_Scores/Modules/disease/celltypeenriched
bed_cell=/n/groups/price/kushal/singlecellLDSC/data/BEDFILES/Modules/disease/celltypeenriched
enhancer_tissue=BLD

module load gcc/6.2.0
module load R

IFS="
"

TASKFILE=/n/groups/price/kushal/singlecellLDSC/data/Tnames.txt

for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
   temp=`echo $line | awk '{print $1}'`
   bed_dir=$bed_cell/$temp
   if [ ! -d $bed_dir ]
   then
       mkdir $bed_dir
   fi
   genescore_dir=$genescore_cell/$temp
   if [ ! -d $genescore_dir ]
   then
       mkdir $genescore_dir
   fi
   for ll in `ls -1 $genescore_dir | sed 's/\.txt//g' | awk '{print $1}' | sort | uniq`;
   do
      annot_name=`echo $ll | awk '{print $1}'`
      echo $temp $annot_name
#     if [ ! -f $bed_dir/$annot_name/100kb.bed ]
#      then
	  cmd="Rscript build_module_annotations.R  $genescore_dir $bed_dir $annot_name $enhancer_tissue"
	  sbatch --time=90:00 --mem=20000 --output=geneS2G.out --error=geneS2G.err -p short -c 1 --wrap="$cmd"
#      fi
   done
done

