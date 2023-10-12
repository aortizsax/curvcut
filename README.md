# curvcut: README: Post-Qiime2 auto remove features/speices
### Adrian Ortiz-Velez


This project was created to proved a mathematical backing to eliminating sparse features in a count tables. 
Input: OTU table tsv or SFS data
Output: Zero filtered table ready for downstream analysis

In this directory, pertaining to usage and the data used for the paper
```
├── data
│   ├── building_data
│   │   └── 16S_OTU_Table_building_materials.tsv
│   ├── minor_allele_data
│   │   └── sfs-output.csv
│   ├── fernanda_data
│   │   └── ASV_table.tsv
│   └── periodontal_data
│       ├── 2020-09-18_PT_16S_OTU_Tabl.csv
│       └── 2020-09-18_PT_metagen_table.csv
├── environment.yml
├── README.md
└── scripts
    └── curvcut.py

```

Dependincies:
 * Pandas 1.3.3
 * Numpy 1.19.2
 * Matplotlib 3.3.4
 * SciPy 1.7.1
 * Jenkspy 0.3.3
 * enviroment.yml

To download 

```bash
git clone https://github.com/aortizsax/curvcut.git
```

```bash
cd curvcut
```

```bash
conda env create -f environment.yml
conda activate curvcut
```

To run in commandline, add script to path (instructions depend on the system). Usage below.

```bash
python3 curvcut.py -h
```

```bash
usage: curvcut.py [-h] [-ct FILE] [-af FILE] [-o OUTPUT_PREFIX]
                  [-u USER_CUTOFF]

optional arguments:
  -h, --help            show this help message and exit
  -ct FILE, --counttable FILE
                        Path to source file.
  -af FILE, --minorallelefreq FILE
                        Path to source file.
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix for output files [default=output].
  -u USER_CUTOFF, --user-cutoff USER_CUTOFF
                        Prefix for output files [default=-1].

```

Example:

	cd ./data/building_data/
	python3 curvcut.py -ct 16S_OTU_Table_building_materials.tsv
	
Output:
* Makes a new directory in the same as the input table
* Saves graphs of finding the max curvature
* Saves saves visual for user diagnosis
* Saves zero filtered table in the new directory
* Outputs code below for example
```
Cubic Spline Finshed
Execution time in seconds: 1.1464474201202393
Recommended  Cutoff: Trim Features that are present in less than 2.35 samples

Zero Filtered OTU table saved: ./output/16S_OTU_Table_building_materialstable.zerofiltered.csv
```
If there is a permission error trying to access your root, forgive me. Please to the following
* Try including './' in front of the data (ex:./16S_OTU_Table_building_materials.tsv) 

<![plot](./scripts/output/2020-09-18_PT_16S_OTU_Tablprocessingcutofffinal.svg)>


