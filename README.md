# CASIO (Combining Annotation Software to Identify Orthologous genes)  
(An example of an application to the Papilionidae family)

## 1) Genome Annotation

Three different programs were used to annotate the genomes:

### BUSCO (Benchmarking Universal Single-Copy Orthologs) (version 5.7.0) (Manni et al. 2021)

Used to assess the integrity of genome assemblies and for gene prediction using **AUGUSTUS**. Detailed operation:

-   **BUSCO database**: Uses a library of single-copy orthologs, conserved and present in a unique form in most of the genomes of a taxonomic group. The library adapted to your group must be defined (`GROUPE_odb10`).
-   **Alignment and ortholog search**: BUSCO compares the sequences provided with the ortholog database using tools such as **BLAST+** or **HMMER**.
-   **Classification of orthologs** :
    -   *Complete*: orthologs found in their entirety.
    -   *Complete and duplicated*: orthologs present in several copies.
    -   *Fragmented*: orthologs partially found.
    -   *Missing*: orthologs not found.
-   **Integrity assessment**: Generation of reports showing the percentage of complete, duplicated, fragmented and missing genes.

BUSCO generates one FASTA file per gene. All the FASTA files for the complete classified genes (`single_copy_busco_sequences`) are compiled into a single FASTA file for each species.

### Miniprot (version 0.11_r234) (Li 2023)

Used to annotate genomes based on protein sequences. Detailed operation:

-   **Protein sequence pre-processing**: Miniprot uses well-annotated protein sequences for mapping.
-   **Protein mapping**: The sequences are aligned with the target genome to identify the corresponding regions.
-   **Exon and intron identification**: Miniprot identifies exons (coding) and introns (non-coding).
-   **Annotation transfer**: Miniprot transfers protein annotations to the target genome.
-   **Report generation** : Results are produced in GFF format.

### Scipio (scipio.1.4.1.pl) (Keller et al. 2008)

Scipio is used to identify coding genes and transfer functional annotations. Detailed operation:

-   **Protein alignment**: Scipio aligns protein sequences to the target genome.
-   **Exon and intron identification**: Identification of exons and introns in aligned regions.
-   **Gene annotation**: Scipio transfers functional annotations of proteins to the target genome.
-   **Results**: Scipio produces files in YAML format with gene annotations.

The genomes are annotated in three phases, with each phase using protein sequences from a different species. Each genome is therefore annotated three times.

The YAML files from Scipio are converted into GFF, and then several steps are carried out to process and clean up the results of the Miniprot and Scipio annotations:

### Processing GFF files

1.  **Copy** of the three GFF files.
2.  **Extraction of CDS** from GFF files (`listAlternativeTranscript2.py`).
3.  **Cleaning of GFF files** to eliminate alternative transcripts (`ExcludeCDSFromGFF.py`).
4.  **Reconciliation of cleaned CDS GFF** files.
5.  **Conversion of GFF files to FASTA** to obtain exonic sequences (`putIDinNameGFF.py`).
6.  **Concatenation** of exons into full transcripts (`concateExonFromBedtoolsGetFasta2_forCat.py`).

### Gene Selection and Translation

1.  **Gene sorting** to eliminate overlapping exons (`SelectOverlappingExonGFFsort_forCat.py`).
2.  **Recovery** of single transcripts (`copySeqList.py`).
3.  **Translation of transcripts** into protein sequences.
4.  **List of redundant** transcripts (`RenameTranscriptRedundants_forList.py`).
5.  **Renaming of redundant** transcripts (`RenameTranscriptRedundants_forFasta.py`).
6.  **Selection of final** transcripts (`SelectSeq.py`).
7.  **Transcript renaming** with species-specific information.

Bash scripts have been created to automate these steps for the outputs of **Miniprot** and **Scipio** : - `Annotation_Miniprot_concatenate.bash` - `Annotation_Scipio_concatenate.bash`

After annotation, we processed the results independently for each set of genes produced by **BUSCO**, **Miniprot** and **Scipio**.

------------------------------------------------------------------------

## 2)  Search for Orthologous Genes

Once the genomes were annotated, we used **OrthoFinder** (version 2.5.4) (Emms & Kelly 2015) to identify 1:1 orthologous genes between species. OrthoFinder searches for orthogroups and orthologs, constructs gene trees for orthogroups, and identifies gene duplication events. It generates a rooted species tree for the studied species and maps gene duplication events from the gene trees to the species tree. To run OrthoFinder, a set of protein sequence files for each species in FASTA format is required. Thus, we retrieved the BUSCO sequences, which were grouped into a single FASTA file per species, and the Scipio and Miniprot outputs, which were processed and converted to FASTA format.

### Selection of orthogroups

After running OrthoFinder, we selected the orthogroups for each annotation method by following several steps:

-   **Search for 1:1 orthologous genes shared by at least 5 different species**: The OG content description file generated by OrthoFinder (`Orthogroups.GeneCount.tsv`) was used with a Python program specifically created to retrieve the list of corresponding OGs (`SelectOrthogroupNamesNoParalog.py`). The resulting list is in the format of a text file (`List_orthogroup_names_no_paralog.txt`).
-   **Search for 1:1 orthologous genes shared by at least 5 different species with a paralog**: The same file (`Orthogroups.GeneCount.tsv`) was used to obtain the list of corresponding OGs using another program (`SelectOrthogroupNamesOneParalog.py`). The resulting list is in the format of a text file (`List_orthogroup_names_one_paralog.txt`).
-   **Identification of gene duplication events**: Among the orthogroups with one paralogous gene, we checked whether it was the result of a recent duplication. To do this, we used the gene trees generated by OrthoFinder, available in the `Resolved_Gene_Trees directory`, and from the list of orthogroups with a paralog, we obtained a list of those where the paralogous gene was in a sister or cousin branch of the gene from the same species (`SelectOrthogroupNamesOneParalogSister.py`). The resulting list is in the format of a text file (`List_orthogroup_names_one_paralog_sister.txt`). The paralogous gene, as well as the gene from which it originated, will both be removed from the alignments.

This search for orthologous genes provided the list of orthogroups to be retrieved for gene alignment (`List_orthogroup_for_alignment.txt`).

------------------------------------------------------------------------

## 3)  Alignment of orthologous gene sequences

To align the orthologous gene sequences, we used **SINGULARITY** (3.10.4) to run **MACSE** (Multiple Alignment of Coding Sequences) (omm_macse_v11.05b.sif) (Ranwez et al. 2021), a tool used to align gene sequences while taking codon structure into account, which is crucial for maintaining the consistency of reading frames in coding sequences. The program `omm_macse_v11.05b.sif` includes **HMMCLeaner** as an integrated cleaning tool. Again, these steps were performed independently for each genome annotation method (BUSCO, Miniprot, and Scipio).

Before performing the alignments, a data preparation phase from OrthoFinder outputs was necessary:

-   **Transfer of orthologous folders**: The orthogroup folders (`Orthogroup_Sequences`) from OrthoFinder, containing a file per gene sequence, were copied into a new folder (`Orthogroup_Sequences_faa`) in the alignment directory (`Alignment`), with an independent subdirectory for each set from the three annotation methods (`Align_Busco`; `Align_Miniprot`; `Align_Scipio`).
-   **Retrieval of nucleotide sequences**: The alignment requires nucleotide-level sequences and not amino acid sequences, so it was necessary to retrieve the corresponding sequences by going back to the gene annotation stage.
-   **Concatenation**: All retrieved nucleotide sequences were grouped into a single FASTA file for each different annotation set (`fasta_file.fna`).
-   **Creation of alignments for each orthologous gene**: The alignment is performed on a FASTA file containing all nucleotide sequences of the genes (here orthologous) from the different species. In our case, this corresponds to a FASTA file per orthogroup:
    -   **Copying .faa files while keeping only the gene names**: To retain the structure of the ".faa" files, which contain all genes of a given orthogroup, for each file, we copied the names of the sequences into another ".txt" file (`CopySeqFileByFile_Part1.sh`). This step allowed us to keep the names of the genes grouped in an orthogroup without their associated amino acid sequences for each previously selected orthogroup (`Orthogroup_Sequences_txt`).
    -   **Cleaning the .txt files**: In the ".txt" files, containing only the gene names of each orthogroup, we removed the gene names of the species possessing two paralogous genes. Then, we deleted the ".txt" files that had fewer than 5 genes.
    -   **Copying the .txt files and adding the nucleotide sequences**: To retrieve the nucleotide sequences associated with each gene name in each ".txt" file, we used a program that searches in a FASTA file for gene names (complete or incomplete) contained as a list in a TEXT file and copies those found with their corresponding sequences into an output FASTA file (`SelectSeq.py`; `CopySeqFileByFile_Part2.py`). This step allowed us to obtain a ".fna" file per orthogroup from the ".txt" files as gene lists and the FASTA file containing all nucleotide sequences (`Orthogroup_Sequences_fna`).
    - **Important check**: It is important to ensure that the gene names are not too long and to modify them if necessary, as the alignment program has a maximum number of characters for sequence names.

All the data formatting is necessary before executing the alignment, looped over all obtained ".fna" files. Below is an example of the command and its defined options for the alignment: 
`singularity run -H $PWD ~/bin/omm_macse_v11.05b.sif --out_dir ./Alignment/Align_Busco/Alignment_macse_Busco/${Orthogroup}_macse --out_file_prefix ${Orthogroup}_macse --in_seq_file ./Alignment/Align_Busco/Orthogroup_Sequences_fna/${Orthogroup}.fna --genetic_code_number 1 --alignAA_soft MAFFT --java_mem 10000m`

------------------------------------------------------------------------

## 4)  Alignment Cleaning

The alignments of the orthologous gene sequences were cleaned using the PhylteR program for each annotation method by following similar steps:

### Preparation of alignments for obtaining gene trees 

PhylteR detects aberrant gene family taxa in phylogenomic datasets, allowing for the removal of outlier genes. Thus, we sought to obtain gene trees from our orthologous gene alignments:

-   **Transfer of alignment files**: The nucleotide-level masked alignment files (`{OG_CODE}_macse_final_mask_align_NT.aln`) generated by MACSE were transferred into a directory containing a folder for each alignment (`Alignment_Phylter_Busco`;`Alignment_Phylter_Miniprot`;`Alignment_Phylter_Scipio`) within the directory specific to each annotation method (`PhylteR_Busco`; `PhylteR_Miniprot`; `PhylteR_Scipio`).
-   **Construction of gene trees**: The gene trees for each renamed alignment were constructed using IQTree and its iqtree2 program under the GTR+G4 model.

### Data preparation for PhylteR (Comte et al. 2023)

PhylteR takes as input a collection of phylogenetic trees (which it then converts into distance matrices), as well as a list of gene names ordered to correspond with the gene trees. In our case, these names correspond to those of our orthogroups:

-   **Collection of gene trees**: We concatenated all **IQTree** (2.2.2.6 COVID-edition) (Minh et al. 2020) output files (`{OG_CODE}_macse_final_mask_align_NT.aln.treefile`) in Newick format into a single file (`treefiles.tre`).
-   **Renaming of gene names**: The gene names in the file containing all gene trees (`treefiles.tre`) were standardized to retain only the species name.
-   **List of orthogroup names**: We also retrieved the list of orthogroup names associated with each ".treefile" (`list_treefiles.txt`).
-   **Execution of PhylteR**: To run PhylteR, we provided the corresponding script (`Script_PhylteR.r`) as input to the Rscript command.
-   **Removal of identified outliers**: PhylteR outputs a ".out" file (`phylter_Busco.out`; `phylter_Miniprot.out`; `phylter_Scipio.out`), which contains a list of all taxa detected as outliers and categorized by the name of the gene tree they belong to (the orthogroup names in our case). To process this output file and remove the outlier genes from our alignments, we created a Python program that takes the ".out" file and the directory containing copies of the alignments to filter as input (`PruneAlnFiles.sh`). This step allowed us to obtain alignments cleaned of taxa with aberrant values, renamed "pruned" or "filtered" according to the defined suffix.
-   **Removal of alignments with low species representation**: The "filtered" alignments containing fewer than five species after the removal of outlier species were discarded.

------------------------------------------------------------------------

## 5)  Comparison of alignments from the three annotation methods

After completing the first step of genome annotation, we independently processed the three sets of annotated genes, each produced by a different method: BUSCO, Miniprot, and Scipio. Although these three programs identified different numbers of orthologous genes, it is likely that some of these genes are common to two or all three gene sets.

To determine which alignment to prioritize among those of the orthologous genes found in common, we compared the alignments of the orthologous genes shared by the three annotation methods:

### Search for Orthologous Genes Common to BUSCO, Miniprot, and Scipio

We searched for alignments containing orthologous genes common between the BUSCO, Miniprot, and Scipio datasets based on reference genes from the species *Ornithoptera alexandrae*, *Papilio xuthus*, and *Parnassius apollo*. For these three species, the gene names were identical between the gene sets, as their annotation had not been performed by any of the three programs used. To achieve this, our approach was to create a dictionary for each gene set (B, M, and S) using Python. Each dictionary contained the OGA (Orthologous Gene Alignment) names as keys, and the reference gene names of the three species *O. alexandrae*, *P. xuthus*, and *P. apollo* as values, when at least one of them was present (`GenesBuscoMiniprotScipio.py`). Using this method, we then compared the values of the keys, one by one, between the dictionaries to match two OGAs and categorize them as common whenever a value from one key was identical to a value from another key. The OGA from a specific annotation method's gene set is called by the first letter of the corresponding program (OGA-B, OGA-M, and OGA-S), and when an OGA from one gene set is categorized as common with another, the two are labeled with the initials of the annotation programs of the gene sets they belong to (e.g., OGA-BM if OGA-B is common with OGA-M). To ensure this work was done properly, several sorting and search steps were necessary:

-   **Renaming Sequences**: First, we restored the full names of sequences in the filtered nucleotide alignments so we could use the gene names from *O. alexandrae*, *P. xuthus*, and *P. apollo* as comparison tools between the OGAs from different gene sets.
-   **Sorting Alignments**: We removed OGAs where the reference genes from *O. alexandrae*, *P. xuthus*, and *P. apollo* (values) were grouped differently between two gene sets. For example, if the *O. alexandrae* gene from an OGA-B was found in an OGA-M and the *P. xuthus* gene from the same OGA-B was found in a different OGA-M, they were removed from both the Busco and Miniprot gene sets. We built a Python program (`GenesBuscoMiniprotScipioToRemove.py`) that takes as input the name of the gene set to analyze and the three directories containing the annotation method alignments (OGA) and returns a list of alignments to remove from the specified gene set. We then removed the OGAs from the list of alignments to obtain filtered gene sets.
-   **Initial Categorization of Alignments**: Using the filtered OGA gene sets, we searched for OGAs common to several gene sets and those that were not. To do this:
    -   We retrieved the list of OGA names for each gene set using a Python program (`GenesBuscoMiniprotScipio.py`).
    -   Then, we searched for unique OGAs in each gene set using another program (`GenesBuscoMiniprotScipioUnique.py`), which takes as input the gene set name and the three directories containing the annotation method alignments and returns a list of alignments unique to the specified gene set.
    -   To search for OGAs common to all three annotation methods (referred to as OGA-BMS), we created a CSV file structured as a table listing OGA-BMS names in rows. The first column of this CSV was dedicated to renaming OGA-BMS in order to assign a unique name to OGAs from different gene sets but grouping the same genes. To do this, we used a Python program (`GenesBuscoMiniprotScipioCommonFor3.py`). From this file, we retrieved the list of OGA-BMS per gene set (B, M, and S) as a text file.
    -   We then searched for OGAs common between two gene sets (OGA-BM, OGA-BS, and OGA-MS). For this, a CSV file was created for each gene set pair using an adapted Python program (`GenesBuscoMiniprotScipioCommonFor2.py`). These files listed the OGAs categorized as common between the two specified gene sets. In these files, the first column is dedicated to the categorization "CommonToTwo" or "CommonToThree" for the two OGA names in the next two columns. Then, for each CSV file, we copied the CommonToTwo lines into another CSV to retrieve the OGA-BM/-BS/-MS lists separately from the OGA-BMS.
-   **Handling Issues**: Subsequent to the sorting phase, a number of inconsistencies were detected and addressed:
    -   **OGAs with Different Reference Genes**: We detected OGAs that were categorized as common, even though their reference genes were not consistent. The issue arose because they shared at least one reference gene. For instance, an OGA-S contains the reference genes for *O. alexandrae*, *P. xuthus*, and *P. apollo*. It shares *O. alexandrae* and *P. apollo* reference genes with an OGA-M, but only *P. xuthus* with an OGA-B. In this case, the OGAs were previously classified as OGA-BMS, since OGA-S shared genes with both OGA-M and OGA-B.

        |  | *O. alexandrae* | *P. xuthus* | *P. apollo* |
        |-----------|-----------|-----------|-----------|
        | **OGA-S** | Y | Y | Y |
        | **OGA-B** | abs | Y | **Z** |
        | **OGA-M** | Y | abs | Y |

        Above, recovered example of a categorisation problem (OGAs categorised as BMS when they should not have been.
        The resolution of these specific cases involved the implementation of a dedicated program written in the Bash language (`GenesBuscoMiniprotScipioCommonFor2Dif.sh`), designed for the purpose of searching for OGAs defined as -BMSs that possess a distinct reference gene by means of a pairwise comparison. A similar approach was adopted through the utilisation of a three-way comparison with a highly analogous program (`GenesBuscoMiniprotScipioCommonFor3Dif.sh`). The identification of OGAs was followed by their subsequent removal from the gene sets.

    -   **Recategorization**: After filtering, it was necessary to recategorize the OGAs using the various programs on the newly filtered gene sets (`GenesBuscoMiniprotScipio.py`; `GenesBuscoMiniprotScipioUnique.py`; `GenesBuscoMiniprotScipioCommonFor3.py`; `GenesBuscoMiniprotScipioCommonFor2_2.py`).
    -   **Missing OGA-BMS**: We detected a second problem in the categorization of OGA. Our method for categorizing OGA-BMS involves grouping three OGA from different gene sets whenever one OGA shares reference genes with the other two. As a result, the grouping is made even if two OGA among the three do not directly share reference genes with each other. However, we found that some OGA were present in two distinct categories of OGA common to two gene sets and were not classified under OGA-BMS. This is due to our method of detecting OGA-BMS, which is performed starting only from the OGA-B. Through this approach, we detect among the OGA-B those that share common reference genes with both OGA-M and OGA-S. However, we omit the OGA-M that share common reference genes with OGA-B and OGA-S, without OGA-B and OGA-S sharing reference genes with each other. Similarly, for OGA-S that share common reference genes with OGA-B and OGA-M, without these two sharing reference genes. To fix this problem, we added OGA-BS/-MS and -BM/-MS to the OGA-BMS. Example found of a categorization problem (OGA not categorized as -BMS, although they should have been).

        |  | *O. alexandrae* | *P. xuthus* | *P. apollo* |
        |-----------|-----------|-----------|-----------|
        | **OGA-S** | Y | Y | Y |
        | **OGA-M** | Y | Y | abs |
        | **OGA-B** | abs | abs | Y |

        Above, recovered example of a categorization problem (OGA not categorized as -BMS, although they should have been.

-   **Renaming Common OGAs**: The search for orthologous genes (Part 2) was carried out separately for each gene set from the different annotation programs. As a result, the names of the orthogroups (OGs) grouping the same orthologous genes vary from one gene set to another. At this stage, the names of the annotated orthogroups (OGAs) still did not allow us to recognize those that were common between the different gene sets. To resolve this issue, we renamed the OGAs following a specific scheme. The new name starts with the letters 'OG' for orthogroup, followed by the first letters of the relevant annotation programs: BMS for OGAs common to all three gene sets (OGAs-BMS); BMX, BXS, or XMS for OGAs common to two gene sets (OGAs-BM/-BS/-MS); and BXX, XMX, or XXS for OGAs specific to a single gene set (OGAs-B/-M/-S). Finally, the name ends with a sequence of five digits. For example, the first group of three OGAs-BMS will be named OGBMS00001. Similarly, for two OGAs-MS that appear in the 1234th position in the renaming order, the name will be OGXMS01234.

### Statistics on Common OGAs 

Now that we have identified and renamed the OGAs-BMS, the next step is to determine which of the three annotation methods to prioritize. To do this, we aimed to compare the quality of the OGA-BMS sequences by conducting two analyses:

-   **Sequence Statistics**: In order to compare the quality of sequences between the OGAs-BMS, we analyzed the sequence size, as well as the GAP rate and GC content, using at least one Python program (`seqStats.py`). This program collects these details and returns a CSV file per OGA. These files were then stripped of their headers and concatenated.
-   **dN/dS Calculation**: Next, we calculated the dN/dS values based on the idea that we are comparing the same genes, and therefore, a lower dN/dS ratio should indicate a better alignment. For the dN/dS calculation, we used the bppML and mapNH software:
    - **File Transfer**: The alignment files (already in nucleotide format) were copied into all subdirectories of the folder dedicated to this analysis.
    - **Sequence Renaming**: The alignment files were copied with the suffix “rename,” and the sequence names in these files were shortened to retain only the species name.
    - **Copying Configuration Files**: The necessary “.pp” configuration files for bppML and mapNH analyses were copied into all subdirectories of the alignments.
    - **Phylogenetic Tree Generation**: bppML and mapNH analyses require a species-specific tree for the alignment. To obtain a phylogenetic tree for each alignment, we created an R script (`TreeAlnFiles.R`) that uses functions from the ape library to prune the phylogenetic tree, which includes the 30 species of Papilionidae, based on the species present in each alignment.
    - **Running Calculations**: The bppML and mapNH analyses were run in parallel across all subdirectories containing a "rename" alignment file and its associated species tree, using a script to automate the process.
    - **Deleting Unsuccessful Mapping Directories**: Some subdirectories were excluded because the mapping process failed.
    - **Analysis of dN/dS**: The dN/dS analysis was performed using an R script (`ScriptR_Calculs_dN_dS.R`). This script reads the output files from the mapNH analysis (`counts_dN.dnd`; `counts_dS.dnd`). For each subdirectory, it retrieves the OGA name, species names, and their corresponding dN and dS values to construct a dataframe at the gene set level (BUSCO, Miniprot, or Scipio). From these dataframes, the average dN/dS was calculated for each gene set and for each species. Then, the dataframes were merged based on the common OGA-BMS names for case-by-case comparisons.

This section, which first focused on identifying the common orthologous genes between alignments from different annotation programs, allowed us to identify several of these genes common to two or all three gene sets, or unique to their respective annotation systems. Secondly, we were able to determine which alignment to prioritize when two or three alignments group the same genes but come from different annotation protocols. This step was crucial in establishing a database of orthologous genes for the purpose of searching for orthologous genes across all species of Papilionidae.

------------------------------------------------------------------------

## 6) Addition of new species to the alignments obtained by orthoPap

### Creation of HMM Profiles for OrthoPap Alignments

After aligning the orthologous sequences, we built HMM (Hidden Markov Models) profiles using **HMMER** (version 3.3.1) (Eddy 2011) for each orthologous alignment:

-   **Creating Profiles from Alignments**: The process of creating profiles was carried out on a collection of unmasked amino acid alignments from the MACSE outputs (`Alignment_macse_Busco`; `Alignment_macse_Miniprot`; `Alignment_macse_Scipio`).
-   **Collection of Alignments**: The alignments transferred to the directory dedicated to profile creation were selected from the lists of names of the OGAs to be prioritised (first Busco, then Miniprot and Scipio) and renamed according to the names of the common orthogroups (OGs).
-   **Collection of Profiles**: Once the HMM profiles were created using the `hmmbuild` command, we compiled them into a single file and indexed them using the `hmmpress` command for more efficient searching against target sequences.

This step allowed the creation of a robust orthologous gene database for the Papilionidae, providing a valuable tool for future phylogenetic studies and research on molecular evolution and adaptation within this taxonomic group.

### Annotation of new species

The same annotation procedure as for OrthoPap has been used in this example application.

### Search for orthologous genes of new species

To find the orthologous genes of the new species corresponding to those obtained using the OrthoPap method, we used the HMMER and the profiles produced previously.

-   **Search for profiles against sequences**: The `hmmsearch` command has been used between the annotated genomes of the new species and the profiles derived from the OrthoPap.

-   **Search for sequences against profiles**: The `hmmscan` command has been used between the annotated genes of the new species and the profiles derived from the OrthoPap.

-   **Verification of the consistency of the previous two searches**: Potential inconsistencies in the reciprocity search between the profiles corresponding to the new species sequences and the new species sequences corresponding to the profiles were resolved using the R script (`best_reciprocal_hmm.R`).

### Alignment of new species by enrichment

The **MACSE** (version 2.07) (Ranwez et al. 2018) `enrichAlignment` command was used to align new species with those used for OrthoPap. The principle is to add the sequences of the new species to the corresponding alignments according to the previous correspondence search performed on the HMM profiles.

------------------------------------------------------------------------

## Reference of software used

Eddy SR. 2011. Accelerated Profile HMM Searches. PLoS Comput. Biol. 7:e1002195. doi: 10.1371/journal.pcbi.1002195.

Emms DM, Kelly S. 2019. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol. 20:238. doi: 10.1186/s13059-019-1832-y.

Keller O, Odronitz F, Stanke M, Kollmar M, Waack S. 2008. Scipio: Using protein sequences to determine the precise exon/intron structures of genes and their orthologs in closely related species. BMC Bioinformatics. 9:278. doi: 10.1186/1471-2105-9-278.

Li H. 2023. Protein-to-genome alignment with miniprot. Bioinforma. Oxf. Engl. 39:btad014. doi: 10.1093/bioinformatics/btad014.

Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. 2021. BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Mol. Biol. Evol. 38:4647–4654. doi: 10.1093/molbev/msab199.

Minh BQ et al. 2020. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. Evol. 37:1530–1534. doi: 10.1093/molbev/msaa015.

Ranwez V, Chantret N, Delsuc F. 2021. Aligning Protein-Coding Nucleotide Sequences with MACSE. In: Multiple Sequence Alignment: Methods and Protocols. Katoh, K, editor. Springer US: New York, NY pp. 51–70. doi: 10.1007/978-1-0716-1036-7_4.

Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F. 2018. MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. Mol. Biol. Evol. 35:2582–2584. doi: 10.1093/molbev/msy159.
