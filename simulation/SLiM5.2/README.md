## Simulation using SLiM 5.2
Below are example commands to simulate Chromosome 21 data for an eight-way admixed population using SLiM 5.2 (https://messerlab.org/slim/). The founder population samples are from the 1000 Genomes Project and the Human Genome Diversity Project (https://www.internationalgenome.org/data-portal/data-collection).

Download the 1000 Genomes Project and the Human Genome Diversity Project Chromosome 21 VCF and population label files:
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_chr21.recalibrated_variants.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_chr21.recalibrated_variants.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
wget https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz
wget https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz.tbi
wget https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/metadata/hgdp_wgs.20190516.metadata.txt
```

Remove multiallelics and indels for both VCF files using BCFTools (https://samtools.github.io/bcftools/):
```
bcftools view -m2 -M2 -v snps ./CCDG_13607_B01_GRM_WGS_2019-02-19_chr21.recalibrated_variants.vcf.gz -Oz -o ./tgp_chr21_no_multiallelics_no_indels.vcf.gz
bcftools index -t ./tgp_chr21_no_multiallelics_no_indels.vcf.gz
bcftools view -m2 -M2 -v snps ./hgdp_wgs.20190516.full.chr21.vcf.gz -Oz -o ./hgdp_chr21_no_multiallelics_no_indels.vcf.gz
bcftools index -t ./hgdp_chr21_no_multiallelics_no_indels.vcf.gz
```

Merge two VCF files:
```
bcftools merge --force-samples --output-type z --output ./tgp_hgdp_chr21_no_multiallelics_no_indels.vcf.gz ./tgp_chr21_no_multiallelics_no_indels.vcf.gz ./hgdp_chr21_no_multiallelics_no_indels.vcf.gz
```

Download genetic map:
```
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
```

Apply phasing and imputation to the merged VCF file using Beagle (https://faculty.washington.edu/browning/beagle/beagle.html):
```
java -Xmx1000g -jar ./beagle.22Jul22.46e.jar gt=./tgp_hgdp_chr21_no_multiallelics_no_indels.vcf.gz map=./plink.chr21.GRCh38.map out=./tgp_hgdp_chr21_no_multiallelics_no_indels_phased
```

Split the merged VCF file by founder populations:
```
awk '{print $1 > ($2 ".samples")}' ./tgp_hgdp_3433ind_8way_population_label.txt
bcftools view -S ./WAS.samples -Oz -o ./tgp_hgdp_WAS.vcf.gz ./tgp_hgdp_chr21_no_multiallelics_no_indels_phased.vcf.gz
bcftools view -S ./AFR.samples -Oz -o ./tgp_hgdp_AFR.vcf.gz ./tgp_hgdp_chr21_no_multiallelics_no_indels_phased.vcf.gz
bcftools view -S ./AMR.samples -Oz -o ./tgp_hgdp_AMR.vcf.gz ./tgp_hgdp_chr21_no_multiallelics_no_indels_phased.vcf.gz
bcftools view -S ./CAS.samples -Oz -o ./tgp_hgdp_CAS.vcf.gz ./tgp_hgdp_chr21_no_multiallelics_no_indels_phased.vcf.gz
bcftools view -S ./EAS.samples -Oz -o ./tgp_hgdp_EAS.vcf.gz ./tgp_hgdp_chr21_no_multiallelics_no_indels_phased.vcf.gz
bcftools view -S ./EUR.samples -Oz -o ./tgp_hgdp_EUR.vcf.gz ./tgp_hgdp_chr21_no_multiallelics_no_indels_phased.vcf.gz
bcftools view -S ./OCE.samples -Oz -o ./tgp_hgdp_OCE.vcf.gz ./tgp_hgdp_chr21_no_multiallelics_no_indels_phased.vcf.gz
bcftools view -S ./SAS.samples -Oz -o ./tgp_hgdp_SAS.vcf.gz ./tgp_hgdp_chr21_no_multiallelics_no_indels_phased.vcf.gz
```

Run simulation using SLiM 5.2:
```
slim -d seed=21 -d mu=1.72443e-08 -d re=1.29e-8 -d sp=0 -d ep=46709983 -d on1=765 -d "op1='./tgp_hgdp_AFR.vcf'" -d on2=408 -d "op2='./tgp_hgdp_AMR.vcf'" -d on3=29 -d "op3='./tgp_hgdp_CAS.vcf'" -d on4=727 -d "op4='./tgp_hgdp_EAS.vcf'" -d on5=658 -d "op5='./tgp_hgdp_EUR.vcf'" -d on6=28 -d "op6='./tgp_hgdp_OCE.vcf'" -d on7=657 -d "op7='./tgp_hgdp_SAS.vcf'" -d on8=161 -d "op8='./tgp_hgdp_WAS.vcf'" -d gen=15 -d an=200 -d rn=400 -d "ts='./slim_admix_8way_a200_r400.ts'" -d "vcf='./slim_admix_8way_a200_r400.vcf'" ./recipe_admix_8way.slim
```

Split the simulated VCF file by reference and admix panels:
```
bcftools query -l ./slim_admix_8way_a200_r400.vcf > ./slim_admix_8way_a200_r400_names.txt
python generate_two_columns.py ./slim_admix_8way_a200_r400_names.txt ./slim_admix_8way_a200_r400_names_mapped.txt AFR AMR CAS EAS EUR OCE SAS WAS ADMIX
bcftools reheader -s ./slim_admix_8way_a200_r400_names_mapped.txt -o ./slim_admix_8way_a200_r400_renamed.vcf ./slim_admix_8way_a200_r400.vcf
bcftools query -l ./slim_admix_8way_a200_r400_renamed.vcf > ./slim_admix_8way_a200_r400_names_renamed.txt
awk '/^ADMIX/ {print > "./slim_admix_8way_a200_r400_names_admix.txt"; next} {print > "./slim_admix_8way_a200_r400_names_reference.txt"}' ./slim_admix_8way_a200_r400_names_renamed.txt
bcftools view ./slim_admix_8way_a200_r400_renamed.vcf -Oz -o ./slim_admix_8way_a200_r400_renamed.vcf.gz
bcftools index --csi ./slim_admix_8way_a200_r400_renamed.vcf.gz
bcftools view -S ./slim_admix_8way_a200_r400_names_reference.txt -Ov -o ./slim_admix_8way_a200_r400_renamed_reference.vcf ./slim_admix_8way_a200_r400_renamed.vcf.gz
bcftools view -S ./slim_admix_8way_a200_r400_names_admix.txt -Ov -o ./slim_admix_8way_a200_r400_renamed_admix.vcf ./slim_admix_8way_a200_r400_renamed.vcf.gz
```

Create population label file of reference panel:
```
awk -F'_' 'BEGIN {print "#Sample_ID\tAncestry_Label"} {print $0 "\t" $1}' ./slim_admix_8way_a200_r400_names_reference.txt > ./slim_admix_8way_a200_r400_population_labels.txt
```

Extract ancestry ground truth using tskit (https://tskit.dev/software/tskit.html):
```
python ./get_al_from_ts.py 8 ./slim_admix_8way_a200_r400.ts ./slim_admix_8way_a200_r400_ground_truth.txt 8 0 AFR AMR CAS EAS EUR OCE SAS WAS
```
