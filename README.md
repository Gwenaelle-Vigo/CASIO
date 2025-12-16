# CASIO (Combining Annotation Software to Identify Orthologous genes)  

------------------------------------------------------------------------

## Description

**CASIO (Combining Annotation Software to Identify Orthologous genes)** is an automated and containerized bioinformatics workflow designed to identify, filter, and curate sets of orthologous genes from genome assemblies. CASIO integrates multiple complementary gene annotation strategies, including *ab initio* gene prediction with BRAKER2 (Brůna et al. 2021), homology-based annotation, and protein-to-genome alignment approaches implemented in Miniprot (Li 2023) and Scipio (Keller et al. 2008), while BUSCO (Manni et al. 2021) is used to retrieve conserved orthologous genes, in order to maximize ortholog recovery. The workflow harmonizes outputs from different annotation tools by standardizing formats, removing short or redundant sequences, filtering alternative isoforms and paralogs, and retaining only reliable orthologous groups. CASIO relies on OrthoFinder (Emms & Kelly 2015) for orthogroup inference and provides downstream steps for sequence alignment with MACSE (Ranwez et al. 2018), alignment cleaning with HMMCleaner (Di Franco et al. 2019) and PhylteR (Comte et al. 2023), and quality assessment. Implemented in Snakemake and distributed as a fully self-contained container, CASIO ensures reproducibility, scalability, and ease of deployment across computational environments. Its modular design allows users to customize annotation methods, filtering thresholds, and ortholog selection criteria through a simple configuration file, making CASIO adaptable to diverse taxonomic groups and research objectives in comparative genomics and phylogenomics.

------------------------------------------------------------------------

## Installation

CASIO can be installed using the provided installation script, which automatically sets up all required components.

1. Download the `INSTALL.sh` file from the repository.
2. Move to the directory where you want CASIO to be installed:
```bash
cd /path/to/your/desired/directory
```
3. Make the installation script executable:
```bash
chmod +x INSTALL.sh
```
4. Run the installation script:
```bash
./INSTALL.sh
```
The installation process will download and configure all necessary dependencies and prepare the CASIO workflow for use. Once completed, the pipeline can be launched by editing the configuration file and running Snakemake (installable [here](https://gist.github.com/RomainFeron/da9df092656dd799885b612fedc9eccd)) as described in the **Usage** section.

------------------------------------------------------------------------

## Usage

CASIO is executed through a Snakemake workflow and is fully controlled via a configuration file. This design allows users to customize each step of the pipeline while ensuring full reproducibility.

### 1. Configure the pipeline

Before running CASIO, edit the configuration file (`configparam.yaml`) to define your analysis parameters. The configuration file is fully documented and guides the user through the available options, including:

- Selection of the annotation programs to run (BRAKER2, BUSCO, Miniprot, Scipio)
- Choice of the BUSCO lineage dataset
- Minimum protein length thresholds applied to Miniprot and Scipio predictions
- Minimum number of species required to retain a gene or alignment at different filtering steps
- OrthoFinder parameters and tree inference method
- Computational resources (threads, memory, temporary directories)

### 2. Launch the workflow

Once the configuration file is set, run the pipeline from the CASIO directory using Snakemake:

```bash
snakemake -s Database_OG.smk --configfile configparam.yaml --cores <N> --use-singularity --singularity-args
```
Replace `<N>` with the number of CPU cores to allocate.

If the directory containing the input data is located outside the CASIO working directory, add the following option: `-B /path/to/data_directory/`, replacing the path with the appropriate location of your data directory.

### 3. Outputs

Upon completion, CASIO automatically generates:

- Annotated gene sets for each selected annotation method
- Orthogroup inference results produced by OrthoFinder
- Filtered and cleaned multiple sequence alignments of orthologous genes
- Summary statistics and quality metrics for the recovered orthologs
- Comparative figures of orthologous gene recovery across annotation methods

All results are organized in a structured output directory to facilitate downstream analyses, such as phylogenomic inference or molecular evolution studies.

### 4. Reproducibility and portability

CASIO is distributed as a containerized Snakemake workflow. All software dependencies and versions are managed internally, ensuring consistent results across different computing environments and enabling straightforward reuse on new datasets or taxonomic groups.

------------------------------------------------------------------------

## References of software used

Eddy SR. 2011. Accelerated Profile HMM Searches. PLoS Comput. Biol. 7:e1002195. doi: 10.1371/journal.pcbi.1002195.

Emms DM, Kelly S. 2019. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol. 20:238. doi: 10.1186/s13059-019-1832-y.

Keller O, Odronitz F, Stanke M, Kollmar M, Waack S. 2008. Scipio: Using protein sequences to determine the precise exon/intron structures of genes and their orthologs in closely related species. BMC Bioinformatics. 9:278. doi: 10.1186/1471-2105-9-278.

Li H. 2023. Protein-to-genome alignment with miniprot. Bioinforma. Oxf. Engl. 39:btad014. doi: 10.1093/bioinformatics/btad014.

Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. 2021. BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Mol. Biol. Evol. 38:4647–4654. doi: 10.1093/molbev/msab199.

Minh BQ et al. 2020. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. Evol. 37:1530–1534. doi: 10.1093/molbev/msaa015.

Ranwez V, Chantret N, Delsuc F. 2021. Aligning Protein-Coding Nucleotide Sequences with MACSE. In: Multiple Sequence Alignment: Methods and Protocols. Katoh, K, editor. Springer US: New York, NY pp. 51–70. doi: 10.1007/978-1-0716-1036-7_4.

Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F. 2018. MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. Mol. Biol. Evol. 35:2582–2584. doi: 10.1093/molbev/msy159.
