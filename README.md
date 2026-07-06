# Recomb-Mix
Contact Author: Yuan Wei (yuan.wei@ucf.edu)

## Run the Recomb-Mix Program

The Recomb-Mix program uses C++ Boost and OpenMP Libraries (https://www.boost.org/, https://www.openmp.org/), and is compiled using GCC 9.1.0 with -Os optimization flag under a 64-bit Unix-based operating system:
```
g++ -std=c++17 -fopenmp RecombMix.cpp -l boost_iostreams -o RecombMix -Os
```
The Recomb-Mix program has the following parameters:
- p, or panel `<INPUT PANEL FILE>`, where `<INPUT PANEL FILE>` is the input reference panel path and file name (required).
- q, or query `<INPUT QUERY FILE>`, where `<INPUT QUERY FILE>` is the input admixture panel path and file name (required).
- g, or genetic `<INPUT GENETIC MAPPING FILE>`, where `<INPUT GENETIC MAPPING FILE>` is the input genetic mapping path and file name (required).
- a, or ancestry `<INPUT POPULATION ANCESTRY FILE>`, where `<INPUT POPULATION ANCESTRY FILE>` is the input population labels of reference panel path and file name (required).
- o, or output `<OUTPUT DIRECTORY PATH>`, where `<OUTPUT DIRECTORY PATH>` is the output directory path for all files (optional; default is the current directory).
- i, or inferred `<OUTPUT INFERRED FILE NAME>`, where `<OUTPUT INFERRED FILE NAME>` is the output inferred local ancestry file name (optional; default is inferred_local_ancestral_values.txt).
- e, or weight `<WEIGHT>`, where `<WEIGHT>` is the weight of cross population penalty in cost function (optional; default is 1.5).
- f, or frequency `<ALLELE FREQUENCY>`, where `<ALLELE FREQUENCY>` is the minor allele frequency threshold to exclude the allele values for the markers whose minor allele frequencies are below the threshold (optional; default is 0). By default, it is assumed that the reference panel contains the sequence data. If the reference panel contains SNP-array-like data, it is recommended to use this parameter to filter out minor alleles for the markers based on the given threshold.
- s, or estimate `<MAXIMUM GAP PHYSICAL DISTANCE>`, where `<MAXIMUM GAP PHYSICAL DISTANCE>` is the maximum gap physical distance (number of markers) for local ancestry estimation (optional; default is 0: no estimation). A gap refers to a query region with at least one but no more than the maximum gap physical distance markers that are not present in the reference panel. An estimate is made only if the inferred ancestral labels of both the left and right adjacent regions are identical, and the shared ancestral label is used to smooth out the gap.
- t, or threads `<NUMBER OF THREADS>`, where `<NUMBER OF THREADS>` is the number of CPU cores to use (optional; default is the number of available CPU cores).

An example command of running the Recomb-Mix program:
```
./RecombMix -p ./test/reference_panel.vcf -q ./test/admixture_panel.vcf -a ./test/reference_panel_population_labels.txt -g ./maps/example/recombination_map_GRCh37_chr18.txt
```

The command to get the help of the program:
```
./RecombMix -h
```

## Compact Reference Panels
The Recomb-Mix program utilizes compact reference panels for local ancestry inference. A compact reference panel is space-efficient, as it includes only sample templates containing population-level information. The available compact reference panels (located in *./compact_panels/*) were generated in both GRCh37 and GRCh38 builds. The original panels were comprised of the 1000 Genomes Project (1000GP) and the Human Genome Diversity Project (HGDP) (https://www.internationalgenome.org/), and were phased and imputed using [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html). Below are the names and reference populations in the available compact reference panels.

| Name | Reference Populations |
| :--- | :--- |
| twoway_AFR_AMR | AFR (African), AMR (Admixed American) |
| twoway_AFR_EUR | AFR (African), EUR (European) |
| threeway_intra | FIN (Finnish in Finland), GBR (British in England and Scotland), TSI (Toscani in Italia) |
| threeway_inter | AFR (African), EAS (Eastern Asian), EUR (European) |
| sevenway_inter | AFR (African), AMR (Admixed American), EAS (Eastern Asian), EUR (European), OCE (Oceanian), SAS (South Asian), WAS (Western Asian) |
| eightway_inter | AFR (African), AMR (Admixed American), CAS (Central Asian), EAS (Eastern Asian), EUR (European), OCE (Oceanian), SAS (South Asian), WAS (Western Asian) |

## Input and Output Files
Four input files are required to run the Recomb-Mix program: the reference panel file in VCF or compressed VCF (\*.vcf or \*.vcf.gz) format, the admixture panel file in VCF or compressed VCF (\*.vcf or \*.vcf.gz) format, the genetic map file in HapMap text format, and the population labels of reference panel in text format. More than one individual can be included in the admixture panel. The genetic map file uses the HapMap format, whose description should be found in the first line of the file. The format of each line starting with the second line contains four tab-delimited fields: *Chromosome*, *Position(bp)*, *Rate(cM/Mb)*, and *Map(cM)*. Note that the value of the *Rate(cM/Mb)* field is not used. Instead, it is calculated based on the *Position(bp)* and *Map(cM)* fields. If the genetic mapping of the physical position in the VCF file is not found, interpolation is used to estimate the genetic distance of such position. The population labels of the reference panel contain individual's population label per line in a tab-delimited fashion: *Sample id*, *Population label*.

The output file contains three sections: *Ancestry*, *Position*, and *Result*. The *Ancestry* section contains the population labels and their IDs, with one population per line, in tab-delimited format. The *Position* section lists physical positions in ascending order, with one position per line. The *Result* section contains the inferred ancestry labels for each individual haplotype in the admixture panel, in tab-delimited format. Each line represents the result of one individual haplotype, starting with *Admixture individual haplotype id*, followed by a list of inferred segments, having three fields: *Physical start position*, *Physical end position*, *Inferred ancestry label id*. The *Inferred ancestry label id* is a zero-based index of population labels, and the corresponding population labels can be found in the *Ancestry* section.

## Example
Below is a tutorial on how to run Recomb-Mix using the example provided in the test folder in this repository. This example infers local ancestry labels of each site for three admixed individuals (data is in *./test/admixture_panel.vcf* file), using 30 reference individuals (data of 10 Africans, 10 Europeans, and 10 Asians is in *./test/reference_panel.vcf* file, and their population labels data is in *./test/reference_panel_population_labels.txt* file). The recombination rates used for the inference are loaded from a recombination map (data is in *./maps/example/recombination_map_GRCh37_chr18.txt* file).
1. Clone the repository in a local directory and compile
```
git clone https://github.com/ucfcbb/Recomb-Mix.git
cd Recomb-Mix
make
```
2. Run the Recomb-Mix example in the local directory
```
./RecombMix -p ./test/reference_panel.vcf -q ./test/admixture_panel.vcf -a ./test/reference_panel_population_labels.txt -g ./maps/example/recombination_map_GRCh37_chr18.txt
```
The inferred ancestry labels of each admixed individual haplotype per site are output to a file in the current directory (result data is in *./inferred_local_ancestral_values.txt*). The output file format can be found in [Input and Output Files](https://github.com/ucfcbb/Recomb-Mix/tree/main?tab=readme-ov-file#input-and-output-files) section.
