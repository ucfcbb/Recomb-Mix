# Recomb-Mix
Contact Author: Yuan Wei (yuan.wei@ucf.edu)

## Run Recomb-Mix Program

Recomb-Mix program uses C++ Boost Libraries (https://www.boost.org/), and is compiled using GCC 9.1.0 with -Os optimization flag under a 64-bit Unix-based operating system:
```
g++ -std=c++17 RecombMix.cpp -l boost_iostreams -o RecombMix_v0.6 -Os
```

Recomb-Mix program has below parameters:
- p, or panel `<INPUT PANEL FILE>`, where `<INPUT PANEL FILE>` is the input reference panel path and file name (required).
- q, or query `<INPUT QUERY FILE>`, where `<INPUT QUERY FILE>` is the input admixture panel path and file name (required).
- g, or genetic `<INPUT GENETIC MAPPING FILE>`, where `<INPUT GENETIC MAPPING FILE>` is the input genetic mapping path and file name (required).
- a, or ancestry `<INPUT POPULATION ANCESTRY FILE>`, where `<INPUT POPULATION ANCESTRY FILE>` is the input population labels of reference panel path and file name (required).
- o, or output `<OUTPUT DIRECTORY PATH>`, where `<OUTPUT DIRECTORY PATH>` is the output directory path for all files (optional; default is the current directory).
- i, or inferred `<OUTPUT INFERRED FILE NAME>`, where `<OUTPUT INFERRED FILE NAME>` is the output inferred local ancestry file name (optional; default is admix_inferred_ancestral_values_local.txt).
- e, or weight `<WEIGHT>`, where `<WEIGHT>` is the weight of cross population penalty in cost function (optional; default is 1.5).
- f, or frequency `<ALLELE FREQUENCY>`, where `<ALLELE FREQUENCY>` is the minor allele frequency threshold to exclude the allele value whose minor allele frequency is below the threshold (optional; default is 0).
- u, or outputcompactpanel `<IDENTIFIER>`, where `<IDENTIFIER>` (0 or 1) specifies whether the program outputs a compact reference panel (optional; default is 0: no output).

An example command of running the Recomb-Mix program:
```
./RecombMix_v0.6 -p ./test/reference_panel.vcf -q ./test/admixture_panel.vcf -a ./test/reference_panel_population_labels.txt -g ./maps/recombination_map_GRCh37_chr18.txt
```

The command to get the help of the program:
```
./RecombMix_v0.6 -h
```

## Generate Compact Reference Panel
Recomb-Mix program can generate a compact panel from a given reference panel and population labels of the reference panel. Below is an example command:
```
./RecombMix_v0.6 -p ./test/reference_panel.vcf -a ./test/reference_panel_population_labels.txt -o ./result/ -u 1
```

The generated compact reference panel file and its population labels file are saved in the given output folder. They can be reused for future ancestry inference queries. Below is an example command:
```
./RecombMix_v0.6 -p ./result/compact_reference_panel.vcf -q ./test/admixture_panel.vcf -a ./result/compact_reference_panel_population_labels.txt -g ./maps/recombination_map.txt -o ./result/
```

One can generate a compact reference panel while making local ancestry inference calls on given queries against a given reference panel. The above commands are equivalent to the below one:
```
./RecombMix_v0.6 -p ./test/reference_panel.vcf -q ./test/admixture_panel.vcf -a ./test/reference_panel_population_labels.txt -g ./maps/recombination_map.txt -o ./result/ -u 1
```

## Input and Output Files
Four input files are required to run Recomb-Mix program: the reference panel file in VCF or compressed VCF (\*.vcf or \*.vcf.gz) format, the admixture panel file in VCF or compressed VCF (\*.vcf or \*.vcf.gz) format, the genetic map file in HapMap text format, and the population labels of reference panel in text format. More than one individual can be included in the admixture panel. The genetic map file uses the HapMap format, whose description should be found in the first line of the file. The format of each line starting with the second line contains four tab-delimited fields: *Chromosome*, *Position(bp)*, *Rate(cM/Mb)*, and *Map(cM)*. Note that the value of the *Rate(cM/Mb)* field is not used. Instead, it is calculated based on the *Position(bp)* and *Map(cM)* fields. If the genetic mapping of the physical position in the VCF file is not found, interpolation is used to estimate the genetic distance of such position. The population labels of the reference panel contain individual's population label per line in a tab-delimited fashion: *Sample id*, *Population label*.

The output file contains the inferred ancestry labels of each individual haplotype in the admixture panel in a tab-delimited text format. Each line represents the result of one individual haplotype, starting with *Admixture individual haplotype id*, followed by a list of inferred segments, having three fields: *Physical start position*, *Physical end position*, *Inferred ancestry label id*. The *Inferred ancestry label id* is a zero-based indexing of population labels in the order of appearance in input population labels of the reference panel file, which can be found in the first line of the file.
