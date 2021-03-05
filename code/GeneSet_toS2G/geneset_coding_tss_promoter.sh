module load gcc/6.2.0
module load R

annot_cell=$1
annot_name=$2
baseline_cell=$3

IFS="
"

Rscript geneset_coding_tss_promoter.R $annot_cell $annot_name $baseline_cell
#cmd="Rscript geneset_coding_tss_promoter.R $annot_cell $annot_name $baseline_cell"
#sbatch --time=60:00 --mem=20000 --output=geneset.out --error=geneset.err -p short -c 1 --wrap="$cmd"

