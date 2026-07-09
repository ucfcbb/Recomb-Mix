# Recomb-Mix utility

## Create Your Own Compact Reference Panel

The Compact Reference Panel Generator can convert a given original reference panel in VCF format to its compact version. The compact reference panel is space-efficient and can be reused for future ancestry inference queries. 

Below is an example command to compile:
```
g++ -std=c++17 ./utility/generate_compact_panel.cpp -o CRPGenerator -l boost_iostreams -Os
```

The Compact Reference Panel Generator program has the following parameters:
- p, or panel `<INPUT PANEL FILE>`, where `<INPUT PANEL FILE>` is the input panel path and file name (required).
- a, or ancestry `<INPUT POPULATION LABEL FILE>`, where `<INPUT POPULATION LABEL FILE>` is the input panel population label path and file name (required).
- o, or output `<OUTPUT COMPACT PANEL FILE>`, where `<OUTPUT COMPACT PANEL FILE>` is the output compact panel path and file name (required).
- c, or compact `<OUTPUT COMPACT POPULATION LABEL FILE>`, where `<OUTPUT COMPACT POPULATION LABEL FILE>` is the output compact panel population label path and file name (required).

An example command of running the Compact Reference Panel Generator program:
```
./CRPGenerator -p ./test/reference_panel.vcf -a ./test/reference_panel_population_labels.txt -o ./test/compact_reference_panel.vcf -c ./test/compact_reference_panel_population_labels.txt
```

## Convert Your Local Ancestry Calls to Other Formats

The Ancestry Call Converter can convert the local ancestry inference result from the run-length encoding format (the default output format of the Recomb-Mix program) to ancestry dosage and ancestry state formats for downstream analysis.

```
Run-length Encoding
├── Ancestry
│   ├── population_label
│   └── population_id
├── Position
│   └── physical_position
└── Result
    ├── individual_haplotype_id
    └── interval
        ├── start_physical_position
        ├── end_physical_position
        └── ancestry_id
```

```
Ancestry Dosage
└── Result
    ├── physical_position
    └── individual_diploid_ancestry_dosage
```

```
Ancestry State
└── Result
    ├── physical_position
    └── individual_haplotype_ancestry_identifier
```

Below is an example command to compile:
```
g++ -std=c++17 ./utility/convert_ancestry_call.cpp -l boost_iostreams -o ACConverter -Os
```

The Ancestry Call Converter program has the following parameters:
- i, or input `<INPUT RESULT FILE>`, where `<INPUT RESULT FILE>` is the input local ancestry inference result (run-length encoding format) path and file name (required).
- o, or output `<OUTPUT RESULT FILE>`, where `<OUTPUT RESULT FILE>` is the output local ancestry inference result path and file name (required).
- f, or format `<IDENTIFIER>`, where `<IDENTIFIER>` is the output local ancestry inference result file format (0: ancestry state; 1: ancestry dosage).
   
An example command of running the Ancestry Call Converter program, converting the run-length encoding format into the ancestry state format:
```
./ACConverter -i ./inferred_local_ancestral_values.txt -o ./inferred_local_ancestral_values_as_ancestry_state.txt -f 0
```

Another example command of running the Ancestry Call Converter program, converting the run-length encoding format into the ancestry dosage format:
```
./ACConverter -i ./inferred_local_ancestral_values.txt -o ./inferred_local_ancestral_values_as_ancestry_dosage.txt -f 1
```
