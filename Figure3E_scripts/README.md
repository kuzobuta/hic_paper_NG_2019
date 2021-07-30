# README

[210730]

**Custom scripts for compartment enrichment analysis**

Figure3E in Miura et al., Nature genetics (2019)


## Softwares

- **[parallel](https://www.gnu.org/software/parallel/) (GNU)**
- **[bedtools](https://bedtools.readthedocs.io/en/latest/) (v2.27.1)**
- **[UCSC tools](http://hgdownload.cse.ucsc.edu/admin/exe/)**
	- **bedGraphToBigWig**
	- **wigCorrelate**
- **[juicer tools](https://github.com/aidenlab/juicer/wiki/Download)** (Tested in `juicer_tools_linux_0.8.jar`)  
- **python libraries** (Tested in **python 2.7.16**)
	- *library names (tested ver.)* 
	- **matplotlib (2.0.2)**
	- **numpy (1.16.2)**
	- **joblib (0.14.1)**
	- **scipy (0.19.0)**

## Example scripts

```
#setup
git clone https://github.com/kuzobuta/hic_paper_NG_2019.git
cd hic_paper_NG_2019/
```

**0: Generation of gene density track**

```bash
./00-Generate_GD_file.sh 200000 200kb \
../data/reference/mm9.chrom.sizes \
../data/reference/mm9.refseq.sort.txt \
../data/reference \
mm9
```

**1: Generate A/B compartments (Eig)** 

```bash
./01-Eigen_Extraction_hic_file.sh 200000 200kb \
../data/reference/mm9.chrom.sizes \
mm9 \
../data/reference/mm9.200kb.bins \
../data/reference/mm9.200kb.bins.gene.density.bigWig \
../data/compartment_enrichment \
../data/hic_files/CBMS1_d0_ALL_200k_newfilter_IC_cis.hic
```

**2: Analysis of compartment enrichment (a.k.a. *saddle plot*)**

This run took a long time to finish..
`~30 min` by 2CPUs in MacBook Air

```bash
./02-Compartment-enrichment.sh 200000 200kb \
../data/reference/mm9.chrom.sizes \
mm9 \
../data/compartment_enrichment \
../data/hic_files/CBMS1_d0_ALL_200k_newfilter_IC_cis.hic \
../data/compartment_enrichment/CBMS1_d0_ALL_200k_newfilter_IC_cis.hic_200kb_eigen_mm9_correct.bedGraph
```

**3: Plot the heatmap in different directory (Optional)**

```
python Compartment-Enrichment-juicer-plot-npy.py \
../data/compartment_enrichment/CBMS1_d0_ALL_200k_newfilter_IC_cis.hic_200kb_mean_med_std_oe_percentile_AtoB.npy \
../data/Figures/CBMS1_d0_ALL_200k_newfilter_IC_cis.hic_200kb_contract_enrichment_in_compartments_scale
```



