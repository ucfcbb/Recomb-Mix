## Run Simulation Program
Belows are the example commands simulating the admixed population using SLiM 4.0 (https://messerlab.org/slim/). The input population samples are from the 1000 Genomes Project and the Human Genome Diversity Project (https://www.internationalgenome.org/data-portal/data-collection).

An example command to generate the three-way inter-continental datasets:
```
slim -d seed=18 -d mu=1.46455e-8 -d re=1.29e-8 -d sp=10105 -d ep=80262967 -d on1=108 -d "op1='./AFR_chr18.vcf'" -d on2=103 -d "op2='./EAS_chr18.vcf'" -d on3=99 -d "op3='./EUR_chr18.vcf'" -d gen=15 -d an=200 -d rn=500 -d "ts='./admix_3way_tree_sequence.trees'" -d "vcf='./admix_3way_individuals.vcf'" ./recipe_admix_3way.slim
```

An example command to generate the three-way intra-continental datasets:
```
slim -d seed=18 -d mu=1.46455e-8 -d re=1.29e-8 -d sp=10105 -d ep=80262967 -d on1=99 -d "op1='./FIN_chr18.vcf'" -d on2=91 -d "op2='./GBR_chr18.vcf'" -d on3=107 -d "op3='./TSI_chr18.vcf'" -d gen=15 -d an=200 -d rn=500 -d "ts='./admix_3way_tree_sequence.trees'" -d "vcf='./admix_3way_individuals.vcf'" ./recipe_admix_3way.slim
```

An example command to generate the seven-way inter-continental datasets:
```
slim -d seed=18 -d mu=1.46455e-8 -d re=1.29e-8 -d sp=10089 -d ep=80263251 -d on1=79 -d "op1='./AFR_chr18.vcf'" -d on2=185 -d "op2='./EAS_chr18.vcf'" -d on3=135 -d "op3='./EUR_chr18.vcf'" -d on4=49 -d "op4='./NAT_chr18.vcf'" -d on5=23 -d "op5='./OCE_chr18.vcf'" -d on6=179 -d "op6='./SAS_chr18.vcf'" -d on7=152 -d "op7='./WAS_chr18.vcf'" -d gen=15 -d an=200 -d rn=500 -d "ts='./admix_7way_tree_sequence.trees'" -d "vcf='./admix_7way_individuals.vcf'" ./recipe_admix_7way.slim
```
