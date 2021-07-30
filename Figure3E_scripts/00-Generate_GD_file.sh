#Generate GD tracks
#Used in Miura et al., Nat.Genet. 2019 Fig3e

#Prep GD track 100kb ~ 1Mb
#sizes=(200000)
#names=("200kb")
size=$1
name=$2
#genome_file="references/mm9.chrom.sizes"
#ref_gene="references/mm9.refseq.sort.txt"
genome_file=$3
ref_gene=$4
ref_dir=$5
genome=$6

#for ((i=0;i<${#sizes[*]};i++))
#do
#size=${sizes[i]}
#name=${names[i]}
bedtools makewindows -g ${genome_file} -w ${size} > ${ref_dir}/${genome}.${name}.bins
bedtools intersect -c -a ${ref_dir}/mm9.${name}.bins -b ${ref_gene} | awk '{OFS="\t"}{print}' \
| LC_COLLATE=C sort -k1,1 -k2,2n - > ${ref_dir}/${genome}.${name}.bins.gene.density.bedGraph
bedGraphToBigWig ${ref_dir}/${genome}.${name}.bins.gene.density.bedGraph ${genome_file} ${ref_dir}/${genome}.${name}.bins.gene.density.bigWig
#done


