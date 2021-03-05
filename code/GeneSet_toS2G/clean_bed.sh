bedops_cell=[BEDOPS_DIRECTORY]/BEDOPS/bin
bedtools_cell=[BEDTOOLS_DIRECTORY]/BEDTOOLS/bedtools2/bin
input_cell=$1
names=`ls $input_cell | cut -f 1 -d '.'`

for name in $names
do
$bedtools_cell/bedtools sort -i $input_cell/$name.bed > $input_cell/$name.2.bed
$bedtools_cell/bedtools merge -i $input_cell/$name.2.bed -c 4 -o max > $input_cell/$name.3.bed
mv $input_cell/$name.3.bed $input_cell/$name.bed
rm $input_cell/$name.2.bed
done
