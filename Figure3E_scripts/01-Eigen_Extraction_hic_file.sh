#Eigen extraction from .hic data
#Used in Miura et al., Nat.Genet. 2019 Fig3e

#CPUs for parallel run
CPU=2

#juicer norm
norm="NONE"
juicer_jar_file="juicer_tools_linux_0.8.jar"

size=$1
name=$2
genome_file=$3
genome=$4

#chr list w/o chrY & chrM#
chr_list=(`cut -f1 ${genome_file} | grep -e chrY -e chrM -v`)

#reference/mm9.200kb.bins
bin_file=$5

#reference/mm9.200kb.bins.gene.density.bigWig
GD_track=$6

#eigenvector <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]
#files=(`find juicer_hic -name '*.hic'`)

#out dir
#out_dir="juicer_hic"
out_dir=$7

####################################only 200kb#########################################################
#for ((i=0;i<${#files[*]};i++))
#do
#res=${names[0]}
#res2=${sizes[0]}
#res3="200kb"
#file=${files[i]}
res=${name}
res2=${size}
res3=${name}

file=$8
out="${out_dir}/${name_hic}"

#name_hic=`basename ${files[i]}`
name_hic=`basename ${file}`
out="${out_dir}/${name_hic}"
#1
parallel -j ${CPU} "java -jar ${juicer_jar_file} -p eigenvector ${norm} \
${file} \
{} BP ${res2} \
${out}_${res}_eigen_{}.txt" ::: ${chr_list[@]}
for chr in ${chr_list[@]}
do
#2
grep $chr$'\t' ${bin_file} | paste - ${out}_${res}_eigen_${chr}.txt | awk '{OFS="\t"}{if ($4 != "NaN"){print}}'> ${out}_${res}_eigen_${chr}.bg
grep $chr$'\t' ${bin_file}  | paste - ${out}_${res}_eigen_${chr}.txt | awk '{OFS="\t"}{print}'> ${out}_${res}_eigen_${chr}_allbins.bg

bedGraphToBigWig ${out}_${res}_eigen_${chr}.bg ${genome_file} ${out}_${res}_eigen_${chr}.bw
#3
wigCorrelate ${GD_track} ${out}_${res}_eigen_${chr}.bw
cor=`wigCorrelate ${GD_track} ${out}_${res}_eigen_${chr}.bw | cut -f3`
awk -v cor=$cor '{OFS="\t"}{if (cor < 0){print $1,$2,$3,$4*-1}else{print}}' ${out}_${res}_eigen_${chr}.bg > ${out}_${res}_eigen_${chr}_correct.bg
awk -v cor=$cor '{OFS="\t"}{if (cor < 0){if ($4=="NaN"){print}else{print $1,$2,$3,$4*-1}}else{print}}' ${out}_${res}_eigen_${chr}_allbins.bg > ${out}_${res}_eigen_${chr}_correct_allbins.bg
done
cat ${out}_${res}_eigen_*_correct.bg | bedtools sort -i - > ${out}_${res}_eigen_${genome}_correct.bedGraph
cat ${out}_${res}_eigen_*_correct_allbins.bg | bedtools sort -i - > ${out}_${res}_eigen_${genome}_correct_allbins.bedGraph
rm ${out}_${res}_eigen_*.bg ${out}_${res}_eigen_*.bw ${out}_${res}_eigen_*.txt
#done


