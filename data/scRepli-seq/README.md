# README

scRepli-seq profiles at 80-kb bins  
[d0~d7 & EpiSCs, total 884 S-phase or G1 cells]

## Description of the files
 
### scRepliseq_884cells_Sphase_G1_data_set.Rdata:

Rdata format of binarized scRepli-seq profile (884 cells) at 80-kb bins.

### scRepliseq_884cells_Sphase_G1_data_set_all.txt.zip:

Txt format of binarized scRepli-seq profile at 80-kb bins.  
Same information as "scRepliseq_884cells_Sphase_G1_data_set.Rdata".  
-1: unreplicated, 1: replicated, 0: no data

### scRepliseq_884cells_Sphase_G1_data_set_filtered.txt.zip:

Filtered binarized scRepli-seq profile at 80-kb bins.  
See "Figure4\_01\_Filtering\_scRepliseq\_data.R"

### scRepliseq_884cells_Sphase_G1_data_set_numbers.txt:

The number of each sample set.

### scRepliseq_SPRING_pos.txt:

Output file by "Figure4\_02\_Run\_SPRING.py"  
Position of each single cell profile after SPRING analysis

### All_Sample_Info_filtering.txt:

QC of each single cell data.  
See "Figure4\_04\_Filtering\_by\_Manhattan\_Distance.R"

column information:  
repscores_woX => Repliscores without chrX [replicated bins/ [unreplicated + replicated bins]]
G1\_cells => G1 cell or not (1 means G1 cell)  

IQR\_filter => Outlier cells based on IQR filtering step using SPRING result  
(See "Figure4\_03\_SPRING\_data\_plot.R" or "Figure4\_04\_Filtering\_by\_Manhattan\_Distance.R")  
(1 means outlier cells)  
out\_5\_95\_rep => S phase samples with 0.05 < repliscores < 0.95 as 1  

Man\_filter => S phase samples filtered out by manhattan distance  
(See "Figure4\_04\_Filtering\_by\_Manhattan\_Distance.R")  
(This filtering was only applied for d0 & d7 samples)
