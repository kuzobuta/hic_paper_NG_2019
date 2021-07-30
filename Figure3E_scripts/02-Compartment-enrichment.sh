#Compartment enrichment analysis
#Used in Miura et al., Nat.Genet. 2019 Fig3e
#Python ver. 2.7

#juicer norm
norm="NONE"
juicer_jar_file="juicer_tools_linux_0.8.jar"

genome_file=$3
genome=$4

#chr list#
chr_list=`cut -f1 ${genome_file} | grep -e chrM -e chrY -v`

#eigenvector <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]
#files=(`find juicer_hic -name '*.hic'`)

#out dir
#out_dir="juicer_hic"
out_dir=$5

#res2="200000"
res2=$1
#res="200kb"
res=$2

file=$6
#for ((i=0;i<${#files[*]};i++))
#do
#name_hic=`basename ${files[i]}`
name_hic=`basename ${file}`
out="${out_dir}/${name_hic}"
#file=${files[i]}
for chr in ${chr_list}
do
if [ "$chr" != "chrM" ]
then
java -jar ${juicer_jar_file} dump oe ${norm} \
${file} \
${chr} ${chr} BP ${res2} \
${out}_${res}_oe_IC_${chr}.txt
awk -v "chr"=$chr '{OFS="\t"}{print chr,$1,chr,$2,$3}' ${out}_${res}_oe_IC_${chr}.txt > ${out}_${res}_oe_IC_${chr}_name.txt
rm ${out}_${res}_oe_IC_${chr}.txt
fi
done
cat ${out}_${res}_oe_IC_*_name.txt > ${out}_${res}_oe_IC_ALL.txt
rm ${out}_${res}_oe_IC_*_name.txt

##python scripts##
#Compartment_file=${out}_${res}_eigen_${genome}_correct.bedGraph
Compartment_file=$7

python Compartment-Enrichment-juicer.py ${out}_${res}_eigen_${genome}_correct.bedGraph ${out}_${res}_oe_IC_ALL.txt ${out}_${res}
python Compartment-Enrichment-juicer-plot-npy.py ${out}_${res}_mean_med_std_oe_percentile_AtoB.npy ${out}_${res}

#Test Run 210730
#python Compartment-Enrichment-juicer-plot-npy.py juicer_hic/CBMS1_d0_ALL_200k_newfilter_IC_cis.hic_200kb

#python 10-Compartment-Enrichment-juicer-plot-percentile.py ${out}_${res}_eigen_mm9_correct.bedGraph ${out}_${res}

#done

