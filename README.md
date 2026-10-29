# PLINK Workflow GUI for Integrated Regional Analysis

This repository contains a Python graphical user interface (GUI) for exploratory regional analysis of genomic intervals identified through genome-wide association studies (GWAS). The workflow integrates regional marker extraction, pairwise linkage disequilibrium (LD) calculation, LD-decay visualisation, exploratory SNP–trait analysis, candidate LD region detection, graph-based LD clustering and pairwise genotype-combination analysis.

The current release was developed and tested on Windows 11.

## Main functions

The GUI provides the following modules:

- PLINK Extract
- LD Decay
- LD Decay Graph
- Exploratory SNP–Trait Analysis
- Candidate LD Regions
- Graph-Based LD Clusters
- Pairwise Genotype Combinations

PLINK 2.0 performs regional marker extraction and pairwise LD calculation. Python functions process the resulting files and generate tabular and graphical outputs.

## Repository structure

```text
Scripts
├── initialise_funcs.py
├── plink_funcs.py
├── plinkUtilsGUI.py
├── requirements.txt
├── LICENSE
├── README.md
└── plink
    ├── exe_dir
    │   └── plink2.exe
    ├── Images
    │   └── plink_workflow_logo.png
    ├──├── map_ped
    │   ├── example_chr7_region.map
    │   ├── example_chr7_region.ped
    │   └── Fe_Phene.csv
    ├
    └── setup_config_files
        └── plink_utils_config.txt

```

## Installation

1. Download or clone the repository and retain the supplied folder structure. Extract a downloaded ZIP archive completely.
2. Install Python 3.12 or a compatible environment.
3. Open PowerShell in the directory containing `requirements.txt` and run:

```powershell
python -m pip install -r requirements.txt
```

4. Confirm that `plink\exe_dir\plink2.exe` exists. PLINK documentation and releases are available from <https://www.cog-genomics.org/plink/2.0/>.

## Genotype input

The GUI directly accepts one matching PLINK PED/MAP pair. Place the `.ped` and `.map` files in `plink\map_ped`. Both files must have the same root filename, for example `example_dataset.ped` and `example_dataset.map`. If multiple matching pairs are present, the first detected pair is used.

VCF and PLINK BED/BIM/FAM files require prior conversion to PED/MAP with PLINK.

```powershell
plink2.exe --bfile input_dataset --export ped --out converted_dataset
plink2.exe --vcf input_dataset.vcf --export ped --out converted_dataset
```

## Phenotype input

The phenotype-based modules require `Fe_Phene.csv` in the selected output directory or its parent directory. It must contain `ID`, matching the first PED identifier, and `AWD_Fe57_grain_Y1`, the demonstration phenotype column. For another trait, change the `trait_column` value in `plink_funcs.py` and provide a column with the same name.
The example phenotype file is supplied as `plink\map_ped\Fe_Phene.csv`. Before running **Exploratory SNP–Trait Analysis** or **Pairwise Genotype Combinations** with the example dataset, copy this file into the output directory selected in the GUI.

## Example dataset

The `plink\map_ped` directory contains a matching regional PED/MAP pair for chromosome 7 (`example_chr7_region.ped` and `example_chr7_region.map`) and the corresponding phenotype file (`Fe_Phene.csv`). The genotype files contain 4,894 markers for 260 accessions. The phenotype file contains 282 records, including 239 accessions with matching genotype records and non-missing grain iron measurements.

To test the workflow, retain the example PED/MAP files in `plink\map_ped`, create a writable output directory, copy `Fe_Phene.csv` into that directory, launch the GUI and enter chromosome `7`, lower position `28780000` and upper position `29780000`.

## Running the GUI

Open and run `plinkUtilsGUI.py` in Spyder. The helper scripts do not need to be run separately.

```python
%runfile C:/path/to/Scripts/plinkUtilsGUI.py
```

Alternatively, use PowerShell:

```powershell
cd "C:\path\to\Scripts"
python plinkUtilsGUI.py
```

## Starting an analysis

1. Launch `plinkUtilsGUI.py`.
2. Select the chromosome and enter the lower and upper base-pair coordinates.
3. Click **Output Directory** and select a writable folder.
4. Enter a root filename such as `Plink_Output`.
5. Run the required modules.

Recommended order: PLINK Extract, LD Decay, LD Decay Graph, Exploratory SNP-Trait Analysis, Candidate LD Regions, Graph-Based LD Clusters and Pairwise Genotype Combinations.

The **Candidate LD Regions** module writes
`<root>_candidate_ld_regions.xlsx`. The **Graph-Based LD Clusters** module
reads that workbook and writes `<root>_graph_based_ld_clusters.xlsx`.
Corresponding backend functions are named `candidate_ld_regions()` and
`graph_based_ld_clusters()` to reflect that these procedures summarise LD and
do not infer haplotype blocks.

## Analytical notes

Regional extraction applies chromosome and position criteria but does not independently apply minor-allele-frequency, missingness, heterozygosity or Hardy-Weinberg-equilibrium filters. Appropriate quality control should be completed before regional analysis.

Candidate LD regions and graph-based LD clusters use same-chromosome pairs separated by no more than 1,000 bp with pairwise `r2 >= 0.6`. They are threshold-dependent summaries of regional LD rather than established haplotype blocks.

In phenotype-based modules, heterozygous and missing calls are treated as missing. Markers require homozygous calls for at least 90% of phenotype-eligible samples, and groups with fewer than three samples are excluded. The GUI one-way ANOVA does not adjust for population structure, genomic relatedness or experimental covariates and should be interpreted as exploratory screening.

Pairwise genotype combinations are constructed without estimating chromosomal phase and must not be interpreted as confirmed haplotypes.

## Troubleshooting

- **ModuleNotFoundError:** Keep the three Python files together and launch only `plinkUtilsGUI.py`.
- **PLINK executable not found:** Confirm that `plink\exe_dir\plink2.exe` exists.
- **PED/MAP files not found:** Place a matching pair in `plink\map_ped`.
- **Phenotype file not found:** Place `Fe_Phene.csv` in the selected output directory or its parent.
- **Blank output directory:** Click **Output Directory** and select a writable folder.

## Source code

<https://github.com/parthaomics/GWAS_Downstream_Analysis_Windows11>

## Licence

See `LICENSE` for the repository licence. PLINK is distributed under its own licence and remains the work of its respective developers.
