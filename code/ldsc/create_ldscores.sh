annot_cell=/n/groups/price/kushal/Enhancer_MasterReg/data/ANNOTATIONS/GENE_SCORES2
ldsc_path=/n/groups/price/kushal/LDSC/ldsc
bfile_path=/n/groups/price/kushal/LDSC/1000G_EUR_Phase3_plink
hapmap_path=/n/groups/price/kushal/LDSC/hapmap3_snps

IFS="
"

TASKFILE=/n/groups/price/kushal/Enhancer_MasterReg/data/temp.txt

module load conda2
source activate ldsc

for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
   annot_module=`echo $line | awk '{print $1}'`
   echo $annot_cell $annot_module
   for ll in `ls $annot_cell/$annot_module | awk '{print $1}' | sort | uniq`;
   do
       annot_dir=`echo $ll | awk '{print $1}'`
       echo $annot_dir
       if [ ! -d $annot_cell/$annot_module/$annot_dir ]
       then
	   mkdir $annot_cell/$annot_module/$annot_dir
       fi
       for chrom in {1..22}
       do
       if [ ! -f $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom.l2.ldscore.gz ]
       then
           cmd="~/.conda/envs/ldsc/bin/python $ldsc_path/ldsc.py --bfile $bfile_path/1000G.EUR.QC.$chrom --l2 --ld-wind-cm 1 --yes-really --annot $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom.annot.gz --print-snps $hapmap_path/hm.$chrom.snp --out $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom"
           sbatch --time=90:00 --mem=20000 --output=mega.out --error=mega.err -p short -c 1 --wrap="$cmd"
       fi
    done
  done
done

