# CASIO (Combining Annotation Software to Identify Orthologous genes)  

------------------------------------------------------------------------

## Description

**CASIO (Combining Annotation Software to Identify Orthologous genes)** is an automated and containerized bioinformatics workflow designed to identify, filter, and curate sets of orthologous genes from genome assemblies.

**CASIO** integrates multiple complementary gene annotation strategies, including *ab initio* gene prediction with **BRAKER2** (Brůna et al. 2021), homology-based annotation, and protein-to-genome alignment approaches implemented in **Miniprot** (Li 2023) and **Scipio** (Keller et al. 2008), while **BUSCO** (Manni et al. 2021) is used to retrieve conserved orthologous genes, in order to maximize ortholog recovery. 

The workflow harmonizes outputs from different annotation tools by standardizing formats, removing short or redundant sequences, filtering alternative isoforms and paralogs, and retaining only reliable orthologous groups. CASIO relies on **OrthoFinder** (Emms & Kelly 2015) for orthogroup inference and provides downstream steps for sequence alignment with **MACSE** (Ranwez et al. 2018), alignment cleaning with **HMMCleaner** (Di Franco et al. 2019) and **PhylteR** (Comte et al. 2023), and quality assessment. 

Implemented in **Snakemake** and distributed as a fully self-contained container, CASIO ensures reproducibility, scalability, and ease of deployment across computational environments. Its modular design allows users to customize annotation methods, filtering thresholds, and ortholog selection criteria through a simple configuration file, making CASIO adaptable to diverse taxonomic groups and research objectives in comparative genomics and phylogenomics.

------------------------------------------------------------------------

## Installation

CASIO can be installed using the provided installation script, which automatically sets up all required components.

1. Clone the repository (Needs [git](https://git-scm.com/install/) installed)
```bash
cd /path/to/your/desired/directory
git clone https://github.com/Gwenaelle-Vigo/CASIO
cd CASIO
```
2. Run INSTALL.sh to download apptainer image and configure repository. (Needs wget and md5sum installed)
```bash
bash INSTALL.sh
```
3. Install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (minimum version 8.20)
4. Install [apptainer](https://apptainer.org/docs/admin/main/installation.html)


Once completed, the pipeline can be launched by editing the configuration file, setting input files and running Snakemake as described in the **Usage** section.


------------------------------------------------------------------------
## Usage

CASIO is executed through a Snakemake workflow and is fully controlled via a configuration file. This design allows users to customize each step of the pipeline while ensuring full reproducibility.

### 1. Set input files (Genomes, Queries, Protein Database if Braker used)


CASIO uses 2 kind of input datas (+1 optional depending of annotation methods used).
  - **Genomes** : Define species genome assembly in FASTA format that you want to use to develop a lineage-specific ortholog dataset (often species from the same family). This input dataset takes the form of a FASTA file for each species, stored in a directory with the name Genus_species.fasta.  
The genomes must be in uppercase for the first letter of the Genus and lower case for the rest. The format must respect Genus and species separate by "_". And the extension must be .fasta.  
The INSTALL.sh script automatically creates a directory Data/Genomes/ directory in which you can store your species assemblies. Otherwise, store all your genomes in a directory and set it as assemblies_dir in configparam.yaml.


  - **Queries** : Define reference annotation that are used as homology based annotation and prtein-to-genome alignment. These species should have a good annotation and be from the same family as genomes. This input dataset takes the form of a FASTA file and a GFF file (with CDS feature type) for each reference species, stored in a directory with the name Genus_species.fasta and Genus_species.gff.  
  The genomes must be in uppercase for the first letter of the Genus and lower case for the rest. The format must respect Genus and species separate by "_". And the extension must be .fasta and .gff.  
  It is recommend to not set too many queries to not increase the running time of the workflow. We recommend to use 3 queries.  
The INSTALL.sh script automatically creates a directory Data/Query/ directory in which you can store your reference species assemblies and annotations. Otherwise, store all your FASTA and gff in a directory and set it as queries_dir in configparam.yaml.


  - (Required if Braker method used) **Protein_Database** : If you set Braker in annotation methods, Braker will need a protein database. Braker recommend to use [OrthoDB](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/). Download the most appropriate database for your dataset and store it in a directory.  
The INSTALL.sh script automatically creates a directory Data/Protein_Database/ directory in which you can store the protein database in FASTA format. Otherwise, store the database in a directory and set it as protein_dir in configparam.yaml.

### 2. Configure the pipeline

Before running CASIO, you need to edit the configuration file (`configparam.yaml`) to define your analysis parameters.  This file is fully documented and will guide you through the available options, including:

- Input and Output directory as explained above.
- Selection of the annotation programs to run (BRAKER2, BUSCO, Miniprot, Scipio).
- A list of genomes in a file that have already been masked (to avoid the masking step in the Braker method).
- The minimum length of contigs in genomes to be annotated (short contigs are not useful for annotation, and a genome that is too fragmented can cause bugs in Braker).
- Minimum protein length thresholds applied to CDS predictions.
- Choice of the BUSCO [lineage dataset](https://busco-data.ezlab.org/v5/data/lineages/).
- OrthoFinder parameters and tree inference method.
- The minimum number of species required to retain a gene or alignment at different filtering steps.
- Method to filter orthologous genes without paralogs.

### 2. Launch the workflow

Once the configuration file is set, run the pipeline from the CASIO directory using Snakemake:

```bash
snakemake -s CASIO.smk --configfile configparam.yaml -p --cores <N> --use-singularity --singularity-args /path/to/Data
```
Replace the number of CPU cores to allocate with `<N>`. If possible, use a machine with many cores, as the runtime of the CASIO workflow is long.

If the directory containing the input data is located outside the CASIO working directory, add the following option: `-B /path/to/data_directory/`, replacing the path with the appropriate location of your data directory.

### 3. Outputs

Upon completion, CASIO automatically generates:

- Annotated gene sets for each selected annotation method.
- Orthogroup inference results produced by OrthoFinder.
- Filtered and cleaned multiple sequence alignments of orthologous genes.
- Summary statistics and quality metrics for the recovered orthologs.
- Comparative figures of orthologous gene recovery across annotation methods.

All results are organized in a structured output directory to facilitate downstream analyses, such as phylogenomic inference or molecular evolution studies.

### 4. Reproducibility and portability

CASIO is distributed as a containerized Snakemake workflow. All software dependencies and versions are managed by apptainer, ensuring consistent results across different computing environments and enabling straightforward reuse on new datasets or taxonomic groups.

------------------------------------------------------------------------

## References of software used

- Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics, Volume 26, Issue 6, March 2010, Pages 841–842, https://doi.org/10.1093/bioinformatics/btq033

- Brůna T, Hoff KJ, Lomsadze A, Stanke M, Borodovsky M. 2021. BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database. NAR Genomics Bioinforma. 3:lqaa108. doi: 10.1093/nargab/lqaa108. 

- Comte A et al. 2023. PhylteR: Efficient Identification of Outlier Sequences in Phylogenomic Datasets. Mol. Biol. Evol. 40:msad234. doi: 10.1093/molbev/msad234. 

- Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  
(Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717

- Eddy SR. 2011. Accelerated Profile HMM Searches. PLoS Comput. Biol. 7:e1002195. doi: 10.1371/journal.pcbi.1002195.

- Emms DM, Kelly S. 2019. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol. 20:238. doi: 10.1186/s13059-019-1832-y.

- Flynn JM et al. 2020. RepeatModeler2 for automated genomic discovery of transposable element families. Proc. Natl. Acad. Sci. U. S. A. 117:9451–9457. doi: 10.1073/pnas.1921046117. 

- Keller O, Odronitz F, Stanke M, Kollmar M, Waack S. 2008. Scipio: Using protein sequences to determine the precise exon/intron structures of genes and their orthologs in closely related species. BMC Bioinformatics. 9:278. doi: 10.1186/1471-2105-9-278.

- Li H. 2023. Protein-to-genome alignment with miniprot. Bioinforma. Oxf. Engl. 39:btad014. doi: 10.1093/bioinformatics/btad014.

- Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. 2021. BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Mol. Biol. Evol. 38:4647–4654. doi: 10.1093/molbev/msab199.

- Minh BQ et al. 2020. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. Evol. 37:1530–1534. doi: 10.1093/molbev/msaa015.

- Mölder F et al. 2021. Sustainable data analysis with Snakemake. F1000Research. 10:33. doi: 10.12688/f1000research.29032.3.

- Pertea G, Pertea M. 2020. GFF Utilities: GffRead and GffCompare. doi: 10.12688/f1000research.23297.1.  

- Ranwez V, Chantret N, Delsuc F. 2021. Aligning Protein-Coding Nucleotide Sequences with MACSE. In: Multiple Sequence Alignment: Methods and Protocols. Katoh, K, editor. Springer US: New York, NY pp. 51–70. doi: 10.1007/978-1-0716-1036-7_4.

- Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F. 2018. MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. Mol. Biol. Evol. 35:2582–2584. doi: 10.1093/molbev/msy159.

- Rice P, Longden I, Bleasby A. 2000. EMBOSS: the European Molecular Biology Open Software Suite. Trends Genet. TIG. 16:276–277. doi: 10.1016/s0168-9525(00)02024-2. 

- Shen, W., Le, S., Li, Y. and Hu, F., 2016. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PloS one, 11(10), p.e0163962.

- Tange O. 2023. GNU Parallel 20230522 ( Charles’) [stable]. Zenodo doi: 10.5281/zenodo.7958356. 

- Tegenfeldt F et al. 2025. OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes. Nucleic Acids Res. 53:D516–D522. doi: 10.1093/nar/gkae987. 

- http://repeatmasker.org/

