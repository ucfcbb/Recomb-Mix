# Recomb-Mix
Author: Yuan Wei (yuan.wei@knights.ucf.edu)

## Run Recomb-Mix Program

Recomb-Mix program is compiled using GCC 9.1.0 with -Os optimization flag under a 64-bit Unix based operating system:
```
g++ -std=c++17 RecombMix.cpp -o RecombMix_v0.3 -Os
```

Recomb-Mix program has below parameters:
- p, or panel `<INPUT PANEL FILE>`, where `<INPUT PANEL FILE>` is the input reference panel path and file name.
- q, or query `<INPUT QUERY FILE>`, where `<INPUT QUERY FILE>` is the input admixture panel path and file name.
- g, or genetic `<INPUT GENETIC MAPPING FILE>`, where `<INPUT GENETIC MAPPING FILE>` is the input genetic mapping path and file name.
- a, or ancestry `<INPUT POPULATION ANCESTRY FILE>`, where `<INPUT POPULATION ANCESTRY FILE>` is the input population labels of reference panel path and file name.
- o, or output `<OUTPUT DIRECTORY PATH>`, where `<OUTPUT DIRECTORY PATH>` is the output directory path for all files.
- e, or weight `<WEIGHT>`, where `<WEIGHT>` is the weight of cross population penalty in cost function.
- u, or outputcompactpanel `<IDENTIFIER>`, where `<IDENTIFIER>` (0 or 1) specifies whether the program outputs compact reference panel (default is 0: no output).

An example command of running Recomb-Mix program:
```
./RecombMix_v0.3 -p ./reference_panel.vcf -q ./admixture_panel.vcf -a ./reference_panel_population_labels.txt -g ./recombination_map_GRCh37_chr18.txt -o ./result/ -e 1.5
```

The command to get the help of the program:
```
./RecombMix_v0.3 -h
```

## Generate Compact Reference Panel
Recomb-Mix program can generate compact panel from given reference panel and population labels of reference panel. Below is an example command:
```
./RecombMix_v0.3 -p ./reference_panel.vcf -a ./reference_panel_population_labels.txt -o ./result/ -u 1
```

The generated compact reference panel file and its population labels file are saved in given output folder. They can be reused for future ancestry inference queries. Below is an example command:
```
./RecombMix_v0.3 -p ./result/compact_reference_panel.vcf -q ./admixture_panel.vcf -a ./result/compact_reference_panel_population_labels.txt -g ./recombination_map.txt -o ./result/ -e 1.5
```

One can generate compact reference panel while making local ancestry inference calls on given queries against given reference panel. Above commands are equivalent to below one:
```
./RecombMix_v0.3 -p ./reference_panel.vcf -q ./admixture_panel.vcf -a ./reference_panel_population_labels.txt -g ./recombination_map.txt -o ./result/ -e 1.5 -u 1
```

## Input and Output Files
Four input files are required to run Recomb-Mix program: the reference panel file in VCF format, the admixture panel file in VCF format, the genetic map file in HapMap format, and the population labels of reference panel in text format. More than one individuals can be included in the admixture panel. To perform the local ancestry inference against the reference panel, the values of *CHROM* field and *POS* field in both reference panel and admixture panel VCF files should match. The genetic map file uses the HapMap format, whose description should be found in the first line of the file. The format of each line starting with the second line contains four tab-delimited fields: *Chromosome*, *Position(bp)*, *Rate(cM/Mb)*, and *Map(cM)*. Note that the value of the *Rate(cM/Mb)* field is not used. Instead, it is calculated based on *Position(bp)* and *Map(cM)* fields. If the genetic mapping of the physical position in VCF file is not found, interpolation is used to estimate the genetic distance of such position. The population labels of reference panel contains individual's population label per line in a comma-delimited fashion: *Individual id*, *Population label*.

The output file contains the inferred ancestry labels of each individual haplotypes in admixture panel in a tab-delimited text format. Each line represents result of one individual haplotype, starting with *Admixture individual haplotype id*, following by a list of inferred segments, having three fields: *Physical start position*, *Physical end position*, *Inferred ancestry label id*. The *Inferred ancestry label id* is a zero-based indexing of population labels in the order of appearance in input population labels of reference panel file.
