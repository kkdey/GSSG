module load gcc/6.2.0
module load conda2
source activate ldsc

bedfile_path=$1
bimfile_path=$2
annot_path=$3
geneset=$4
FOLDER=~/GSSG/code/GeneSet_toS2G

IFS="
"

if [ ! -d $annot_path/$geneset ]
then
    mkdir $annot_path/$geneset
fi


for bedline in `ls $bedfile_path/$geneset/ | cat | sort | uniq | cut -f 1 -d '.'`;
do
	bedname=`echo $bedline | awk '{print $1}'`
	if [ ! -d $annot_path/$geneset/$bedname ]
	then
	    mkdir $annot_path/$geneset/$bedname
	fi
	if [ ! -f $annot_path/$geneset/$bedname/$bedname.22.annot.gz ]
	then
	    python  $FOLDER/bedgraph_to_annot.py --bedname $bedname --bedfile_path $bedfile_path/$geneset --bimfile_path $bimfile_path --annot_path $annot_path/$geneset/$bedname
	    #cmd="python  $FOLDER/bedgraph_to_annot.py --bedname $bedname --bedfile_path $bedfile_path/$geneset --bimfile_path $bimfile_path --annot_path $annot_path/$geneset/$bedname"
	    #sbatch --time=150:00 --mem=20000 --output=annot.out --error=annot.err -p short -c 1 --wrap="$cmd"
	fi
done
