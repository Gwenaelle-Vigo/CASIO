import re
import glob
import subprocess
import sys
import os
import gzip
import shutil
from pathlib import Path

# Pipeline

### Function to gzip file
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def compress_gz(filepath):
    with open(filepath, 'rb') as f_in:
            with gzip.open(filepath+".gz", 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
    out=filepath+".gz"
    return out

# gzip file in protein directory
#for fa in glob.glob(config["protein_dir"]+"/*.fa*"):
#    if not is_gz_file(fa):
#        newfile=compress_gz(fa)


ANNOTATIONS=list()

for a in config["annotation_tools"].replace(" ","").split(","):
    if(re.search('Busco', a, re.IGNORECASE)):
        ANNOTATIONS.append("Busco")
    if(re.search('Scipio', a, re.IGNORECASE)):
        ANNOTATIONS.append("Scipio")
    if(re.search('Miniprot', a, re.IGNORECASE)):
        ANNOTATIONS.append("Miniprot")
    if(re.search('Braker', a, re.IGNORECASE)):
        ANNOTATIONS.append("Braker")



def reference_to_query():
    # get all fasta files recursively under DIR_A
    fasta_indir = Path(config["queries_dir"]).glob('**/*fasta*')
    # pair query to infile path using a dictionary
    dquery2referencefa = {}
    dquery2referencegff = {}
    for f in fasta_indir:
        reference = re.sub('.fasta.*$', '', f.name)
        gff = re.sub('.fasta.*$', '.gff', str(f))
        if Path(gff).is_file():
            genus, species = reference.split("_")
            q = genus[0] + species[:3]
            dquery2referencefa[q.capitalize()] = str(f)
            dquery2referencegff[q.capitalize()] = str(gff)
        else:
            raise TypeError("\n\n#################\n\nError: Missing input file, "+ str(gff) +" not found.\nThe gff filename should be identic to fasta name : "+ str(f) +" found so search "+ str(gff) + ", but not found !\nDoes "+str(gff)+" exist ?\n\n#################\n")
    return dquery2referencefa, dquery2referencegff

fasta_query, gff_query = reference_to_query()

GENOMES, = glob_wildcards(config["assemblies_dir"] + "/{genome}.fasta")

QUERIES = list(fasta_query.keys())

print("Queries identify (based on Genus_species.fasta and Genus_species.gff) in the " +config["queries_dir"] + " directory are :\n")
print(QUERIES)
print("\n")

print("Genomes identify in the " +config["assemblies_dir"] + " directory are :\n")
print(GENOMES)
print("\n")

print("The pipeline will do the annotation with : ")
print(ANNOTATIONS)
print("\n")

# To avoid ambiguity between query and genome wildcards in rule group_genomefaa_orthofinder and group_queriesfaa_orthofinder
wildcard_constraints:
    query= '|'.join([re.escape(x) for x in QUERIES]),
    genome= '|'.join([re.escape(x) for x in GENOMES]),

config["results_dir"]=re.sub("/*$","",config["results_dir"])
config["assemblies_dir"]=re.sub("/*$","",config["assemblies_dir"])
config["queries_dir"]=re.sub("/*$","",config["queries_dir"])

### rule specifying the files to be generated
rule all:
    input:
        expand(config["results_dir"] + "/Queries/{query}.faa", query=QUERIES),
        expand(config["results_dir"] + "/Queries/{query}.fna", query=QUERIES),
        expand(config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_faa/{genome}.faa", genome=GENOMES, annotation=ANNOTATIONS),
        expand(config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_fna/{genome}.fna", genome=GENOMES, annotation=ANNOTATIONS),
        expand(config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_faa/{query}.faa", query=QUERIES, annotation=ANNOTATIONS),
        expand(config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_fna/{query}.fna", query=QUERIES, annotation=ANNOTATIONS),
        expand(config["results_dir"] + "/alignment/align_{annotation}/List_all_orthogroups_kept_{annotation}.txt", annotation=ANNOTATIONS),
        expand(config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_post_alignment.txt", annotation=ANNOTATIONS),
        expand(config["results_dir"] + "/cleaning/phylter_{annotation}/phylter_seqstats_{annotation}.out", annotation=ANNOTATIONS),
        expand(config["results_dir"] +"/comparison/comparison_{annotation}/", annotation=ANNOTATIONS),
        config["results_dir"] +"/comparison/Full_Table_OG_By_Query_and_Tool.txt",
        config["results_dir"] + "/profiling/orthologue_alignments.hmm",


### Prepare queries
rule generate_queries:
    input:
        reference = lambda wildcards: fasta_query[wildcards.query],
        annotation = lambda wildcards: gff_query[wildcards.query],
    output:
        query_fna = config["results_dir"] + "/Queries/{query}.fna",
        query_faa = config["results_dir"] + "/Queries/{query}.faa",
    params:
        output_dir = config["results_dir"] + "/Queries/",
        minlen = config["minLenCDS"],
    singularity:
        "Container/Braker/repeatM.sif"
    message:
        "Create CDS reference queries in NT and AA format by keeping longest isoform"
    log:
        config["results_dir"] +"/logs/Queries/{query}.log",
    shell:
        """
        agat_sp_keep_longest_isoform.pl --gff {input.annotation} -o {params.output_dir}/{wildcards.query}.gff &>> {log}
        
        gffread -l {params.minlen} -x {params.output_dir}/{wildcards.query}_tmp.fna -g {input.reference} {params.output_dir}/{wildcards.query}.gff &>> {log}
        sed -i "s/^>/>{wildcards.query}_/" {params.output_dir}/{wildcards.query}_tmp.fna
        
        transeq -sequence {params.output_dir}/{wildcards.query}_tmp.fna -outseq {params.output_dir}/{wildcards.query}_tmp.faa &>> {log}
        
        # Remove transcript with codon stop inside the CDS
        awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' {params.output_dir}/{wildcards.query}_tmp.faa |tail -n +2 | paste - - | awk -F "\\t" '{{if($2 !~ /\\*[A-Z*]/){{sub(/\\*$/, "", $2); print $1"\\t"$2}}}}' | tr "\\t" "\\n" > {output.query_faa}
        sed -i 's/_1$//' {output.query_faa}
        
        grep "^>" {output.query_faa} | sed 's/^>//' | sort > {params.output_dir}/list_transcript.txt
        seqkit grep -f {params.output_dir}/list_transcript.txt {params.output_dir}/{wildcards.query}_tmp.fna -w 0 -o {output.query_fna} &>> {log}
        
        rm {params.output_dir}/list_transcript.txt {params.output_dir}/{wildcards.query}_tmp.fna {params.output_dir}/{wildcards.query}_tmp.faa
        """

rule prepare_genome:
    input:
        assembly = config["assemblies_dir"] + "/{genome}.fasta",
    output:
        assembly_corrected = config["results_dir"] + "/Genomes/{genome}.fasta"
    params:
        minLenContigs = config["minLenContigs"],
        output_dir = config["results_dir"] + "/Genomes/"
    message:
        "To avoid too fragmented genomes, apply a filter on chromosomes/contigs by removing chromosomes/contigs with a length below {params.minLenContigs}. This rule rename chromosome/contigs name."
    singularity:
        "Container/Braker/repeatM.sif"
    log:
        config["results_dir"] +"/logs/Genomes/{genome}.log",
    shell:
        """
        seqkit seq -m {params.minLenContigs} -w 0 -o {params.output_dir}/{wildcards.genome}_tmp.fa {input.assembly} &>> {log}
        seqkit replace -p '^.*' -r '{wildcards.genome}_{{nr}}' -w 0 -o {output.assembly_corrected} {params.output_dir}/{wildcards.genome}_tmp.fa &>> {log}
        grep "^>" {output.assembly_corrected} | sed 's/^>//' > {params.output_dir}/{wildcards.genome}_Link_contigs_name_output.txt
        grep "^>" {params.output_dir}/{wildcards.genome}_tmp.fa | sed 's/^>//' > {params.output_dir}/{wildcards.genome}_Link_contigs_name_input.txt
        paste -d '\\t' {params.output_dir}/{wildcards.genome}_Link_contigs_name_input.txt {params.output_dir}/{wildcards.genome}_Link_contigs_name_output.txt > {params.output_dir}/{wildcards.genome}_Link_contigs_name_from_InputGenome_and_FilteredOutputGenome.txt
        rm {params.output_dir}/{wildcards.genome}_tmp.fa {params.output_dir}/{wildcards.genome}_Link_contigs_name_input.txt {params.output_dir}/{wildcards.genome}_Link_contigs_name_output.txt
        """



### Annotation
if "Busco" in ANNOTATIONS:
    ## BUSCO 
    # Download database odb
    rule download_busco:
        output:
            directory("busco_downloads/lineages/"+config["lineage_dataset_busco"]+"/")
        params:
            lineage=config["lineage_dataset_busco"],
        singularity:
            "Container/Busco/Busco.sif",
        message:
            "Download lineage database for busco",
        log:
            config["results_dir"] +"/logs/annotation/Busco/download_busco.log",
        shell:
            """
            busco --download {params.lineage} &> {log}
            """
    # rule for generating single copy sequences using the busco annotation program
    rule busco:
        input:
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
            database=rules.download_busco.output,
        output:
            summary=config["results_dir"] + "/annotation/Busco/Busco_{genome}/{genome}_summary.txt"
        params:
            output_dir = config["results_dir"] + "/annotation/Busco/",
            lineage_dataset = config["lineage_dataset_busco"]
        threads:
            max(int(workflow.cores / 5),1)
        singularity:
            "Container/Busco/Busco.sif"
        message:
            "Generate single copy genes for {params.lineage_dataset} lineage with busco program"
        log:
            config["results_dir"] +"/logs/annotation/Busco/Busco_{genome}/{genome}.log",
        shell:
            """
            busco -f -i {input.assembly} --out_path {params.output_dir} -o Busco_{wildcards.genome} --offline --download_path busco_downloads/ -m genome -c {threads} -l {params.lineage_dataset} --metaeuk &> {log}
            touch {output.summary}
            """
    # rule for generating one file per species
    rule busco_concate:
        input:
            rules.busco.output.summary
        output:
            faa=config["results_dir"] + "/annotation/Busco/Busco_{genome}/{genome}_final.faa",
            fna=config["results_dir"] + "/annotation/Busco/Busco_{genome}/{genome}_final.fna",
        params:
            output_dir = config["results_dir"] + "/annotation/Busco/Busco_{genome}",
            lineage_dataset = config["lineage_dataset_busco"]
        threads: 1
        singularity:
            "Container/Busco/Busco.sif"
        message:
            "By genome, generate a protein fasta file of all single copy genes found"
        shell:
            """
            for file in {params.output_dir}/run_{params.lineage_dataset}/busco_sequences/single_copy_busco_sequences/*.fna
            do
            sample=$(echo $file | awk -F "/" '{{print $NF}}' | sed 's/.fna//')
            sed "1s/.*/>{wildcards.genome}_Busco_${{sample}}/" "{params.output_dir}/run_{params.lineage_dataset}/busco_sequences/single_copy_busco_sequences/${{sample}}.fna" >> {output.fna}
            sed "1s/.*/>{wildcards.genome}_Busco_${{sample}}/" "{params.output_dir}/run_{params.lineage_dataset}/busco_sequences/single_copy_busco_sequences/${{sample}}.faa" >> {output.faa}
            done
            """
if "Braker" in ANNOTATIONS:
    rule repeat_modeler:
        input:
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
        output:
            TE_family=config["results_dir"] + "/annotation/Braker/TE_identification/{genome}/{genome}-families.fa",
        threads:
            max(int(workflow.cores / 5),1)
        params:
            output_dir=config["results_dir"] + "/annotation/Braker/TE_identification/{genome}/",
            skip_masking=config["skip_masking_step"] if config["skip_masking_step"] != "" else "NoGenomeMask",
        log:
            config["results_dir"] + "/logs/annotation/Braker/TE_identification/{genome}/{genome}_RM.log"
        singularity:
            "Container/Braker/repeatM.sif"
        message:
            "Build a database with Repeat Modeler to identify TE in genome {wildcards.genome}"
        shell:
            """
            if [ -f {params.skip_masking} ]; then
                if grep -qi "{wildcards.genome}" {params.skip_masking}; then
                    echo "Skip MASKING step for {wildcards.genome} because write in {params.skip_masking} file" 
                    ln -sr {input.assembly} {output.TE_family}
                else
                    BuildDatabase -name {params.output_dir}{wildcards.genome} {input.assembly} &>> {log}
                    RepeatModeler -database {params.output_dir}{wildcards.genome} -threads {threads} -LTRStruct &>> {log}
                fi
            else
                BuildDatabase -name {params.output_dir}{wildcards.genome} {input.assembly} &>> {log}
                RepeatModeler -database {params.output_dir}{wildcards.genome} -threads {threads} -LTRStruct &>> {log}
            fi
            """
    rule repeat_masker:
        input:
            TE_family=config["results_dir"] + "/annotation/Braker/TE_identification/{genome}/{genome}-families.fa",
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
        output:
            mask=config["results_dir"] + "/annotation/Braker/MASKING_genome/{genome}/{genome}.fasta.masked",
        threads:
            max(int(workflow.cores / 5),1)
        params:
            output_dir=config["results_dir"] + "/annotation/Braker/MASKING_genome/{genome}",
            skip_masking=config["skip_masking_step"] if config["skip_masking_step"] != "" else "NoGenomeMask",
        log:
            config["results_dir"] + "/logs/annotation/Braker/MASKING_genome/{genome}/{genome}_RM.log",
        singularity:
            "Container/Braker/repeatM.sif"
        message:
            "Identify and soft mask TE in genome {wildcards.genome} with Repeat Masker"
        shell:
            """
            if [ -f {params.skip_masking} ]; then
                if grep -qi "{wildcards.genome}" {params.skip_masking}; then
                    echo "Skip MASKING step for {wildcards.genome} because write in {params.skip_masking} file" 
                    ln -sr {input.assembly} {output.mask}
                else
                    t=$(awk -v threads="{threads}" 'BEGIN{{max=1; t=int(threads/4); t=(max > t ? max : t); print t}}')
                    RepeatMasker -xsmall -pa $t -e ncbi -lib {input.TE_family} -dir {params.output_dir} {input.assembly} &>> {log}
                fi
            else
                t=$(awk -v threads="{threads}" 'BEGIN{{max=1; t=int(threads/4); t=(max > t ? max : t); print t}}')
                RepeatMasker -xsmall -pa $t -e ncbi -lib {input.TE_family} -dir {params.output_dir} {input.assembly} &>> {log}
            fi
            """
    rule prot_braker:
        input:
            db_prot=config["protein_dir"],
            queries=expand(config["results_dir"] + "/Queries/{query}.faa", query=QUERIES),
        output:
            proteins=config["results_dir"] + "/annotation/Braker/proteins/proteins.faa",
        message:
            "Create a protein database for Braker2 by concatenating orthogroup proteins sequences in {input.db_prot} and proteins in queries {input.queries}"
        shell:
            """
            for file in {input.db_prot}/*
            do
                if file "$file" | grep -q gzip; then
                    zcat $file >> {output.proteins}
                else
                    cat $file >> {output.proteins}
                fi
            done
            cat {input.queries} >> {output.proteins}
            """
    rule braker:
        input:
            proteins=config["results_dir"] + "/annotation/Braker/proteins/proteins.faa",
            genome_mask=config["results_dir"] + "/annotation/Braker/MASKING_genome/{genome}/{genome}.fasta.masked",
        output:
            gtf=config["results_dir"] + "/annotation/Braker/ANNOTATE_genome/{genome}/braker.gtf",
        params:
            output_dir=config["results_dir"] + "/annotation/Braker/ANNOTATE_genome/{genome}/",
            log_dir=config["results_dir"] + "/logs/annotation/Braker/ANNOTATE_genome/{genome}/",
        singularity:
            "docker://teambraker/braker3:latest"
        threads:
            max(int(workflow.cores / 5),1)
        message:
            "Annotate {wildcards.genome} with Braker2 and {input.proteins} proteins database"
        log:
            config["results_dir"] + "/logs/annotation/Braker/ANNOTATE_genome/{genome}/{genome}_braker.log"
        shell:
            """
            if [ -d ".augustus/species/{wildcards.genome}/" ];then
                rm -r .augustus/species/{wildcards.genome}/
            fi
            if [ -d "{params.output_dir}" ];then
                rm -rf {params.output_dir}/*
            fi
            braker.pl --genome={input.genome_mask} --prot_seq={input.proteins} --threads={threads} --workingdir={params.output_dir} --species={wildcards.genome} &>> {log}
            cp {params.output_dir}/braker.log {params.log_dir}
            """
    rule braker_filter:
        input:
            gtf = config["results_dir"] + "/annotation/Braker/ANNOTATE_genome/{genome}/braker.gtf",
            assembly = config["results_dir"] + "/Genomes/{genome}.fasta",
        output:
            gtf_agat = config["results_dir"] + "/annotation/Braker/Braker_{genome}/{genome}_longest_isoform.gtf",
            gtf_fin = config["results_dir"] + "/annotation/Braker/Braker_{genome}/{genome}_CDS_clean_final.gtf",
            fna = config["results_dir"] + "/annotation/Braker/Braker_{genome}/{genome}_final.fna",
            faa = config["results_dir"] + "/annotation/Braker/Braker_{genome}/{genome}_final.faa",
        params:
            output_dir = config["results_dir"] + "/annotation/Braker/Braker_{genome}/",
            minlen = config["minLenCDS"],
        singularity:
            "Container/Braker/repeatM.sif"
        log:
            config["results_dir"] + "/logs/annotation/Braker/Braker_{genome}/{genome}.log",
        message:
            "Remove isoform in {wildcards.genome} Braker2 annotation by keeping longest isoform. Write fasta NT and AA for kept isoforms"
        shell:
            """
            agat_sp_keep_longest_isoform.pl -gff {input.gtf} -o {output.gtf_agat} &>> {log}
            
            grep -P "\\tCDS\\t" {output.gtf_agat} | bedtools sort > {output.gtf_fin}
            
            # transcript retrieval
            gffread -l {params.minlen} -x {params.output_dir}/{wildcards.genome}_all_transcript.fna -g {input.assembly} {output.gtf_fin} &>> {log}
            sed -i "s/>/>{wildcards.genome}_Braker_/" {params.output_dir}/{wildcards.genome}_all_transcript.fna
            
            transeq -sequence {params.output_dir}/{wildcards.genome}_all_transcript.fna -outseq {params.output_dir}/{wildcards.genome}_all_transcript.faa &>> {log}
            
            # deletion of sequences containing stop codons (*)
            awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' {params.output_dir}/{wildcards.genome}_all_transcript.faa |tail -n +2 | paste - - | awk -F "\\t" '{{if($2 !~ /\\*[A-Z*]/){{sub(/\\*$/, "", $2); print $1"\\t"$2}}}}' | tr "\\t" "\\n" > {output.faa}
            sed -i 's/_1$//' {output.faa}
            
            # generate fna from faa
            grep "^>" {output.faa} | sed 's/^>//' | sort > {params.output_dir}/list_transcript.txt
            seqkit grep -f {params.output_dir}/list_transcript.txt {params.output_dir}/{wildcards.genome}_all_transcript.fna -w 0 -o {output.fna} &>> {log}
            """
if "Miniprot" in ANNOTATIONS:
    # MINIPROT
    # rule for generating mpi files using the miniprot annotation program
    rule miniprot_mpi:
        input:
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
        output:
            mpi=config["results_dir"] + "/annotation/Miniprot/Miniprot_mpi/{genome}.mpi", 
        threads:
            max(int(workflow.cores / 5),1)
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        message:
            "Create an index with miniprot for the {wildcards.genome} genome"
        log:
            config["results_dir"] + "/logs/annotation/Miniprot/Miniprot_mpi/{genome}.log"
        shell:
            """
            miniprot -t {threads} -d {output.mpi} {input.assembly} 2> {log}
            """
    # rule for generating gff files using the miniprot annotation program
    rule miniprot_gff:
        input:
            query=config["results_dir"] + "/Queries/{query}.faa",
            mpi = rules.miniprot_mpi.output.mpi,
        output:
            gff=config["results_dir"] + "/annotation/Miniprot/Miniprot_{query}/{genome}.gff",
        threads:
            max(int(workflow.cores / 5),1)
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        log:
            config["results_dir"] + "/logs/annotation/Miniprot/Miniprot_{query}/{genome}_gff.log"
        message:
            "Protein alignment of {wildcards.query} query with miniprot to annotate {wildcards.genome}"
        shell:
            """
            miniprot -Iut{threads} --gff {input.mpi} {input.query} > {output.gff} 2> {log}
            """
    # rule for cleaning gff files from the miniprot annotation program
    rule miniprot_cds:
        input:
            gff=rules.miniprot_gff.output.gff,
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
        output:
            CDS_gff=config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}/{genome}_{query}_CDS.gff",
            CDS_fna=config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}/{genome}_{query}_CDS.fna",
        params:
            output_dir = config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}",
            minlen = config["minLenCDS"],
        threads: 1
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        message:
            "Extract {wildcards.genome} CDS from gff annotation file of {wildcards.query} query"
        log:
            config["results_dir"] + "/logs/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}/{genome}_{query}_summary.log"
        shell: 
            """
            # extract CDS and rename gene in gff
            grep -P "\\tCDS\\t" {input.gff} | sed "s/;Rank/_{wildcards.query};Rank/" | awk -F "\\t" '{{if($4 <= $5){{print}}}}' | bedtools sort > {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp1.gff
            
            # extract CDS and rename gene in fna
            gffread -l {params.minlen} -o {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp2.gff -x {output.CDS_fna} -g {input.assembly} {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp1.gff &>> {log}
            grep -P "\\tCDS\\t" {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp2.gff > {output.CDS_gff}
            rm {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp2.gff {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp1.gff
            """
    rule miniprot_concate_cds_fna:
        input:
            CDS_fna=expand(config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{{genome}}/{{genome}}_{query}_CDS.fna", query=QUERIES),
        output:
            CDS_ALL_fna=config["results_dir"] + "/annotation/Miniprot/Miniprot_{genome}/{genome}_all_CDS_clean.fna",
        threads: 1
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        message:
            "Concatenate {wildcards.genome} fasta Miniprot annotations from all queries" 
        shell:
            """
            cat {input.CDS_fna} | cut -d " " -f 1 > {output.CDS_ALL_fna}
            """
    rule miniprot_clean_cds:
        input:
            CDS_gff=config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}/{genome}_{query}_CDS.gff",
            CDS_ALL_fna=config["results_dir"] + "/annotation/Miniprot/Miniprot_{genome}/{genome}_all_CDS_clean.fna",
        output:
            CDS_clean_gff=config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}/{genome}_{query}_CDS_clean.gff",
        params:
            output_dir = config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}",
        threads: 1
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        message:
            "Remove isoform in {wildcards.genome} Miniprot annotation from {wildcards.query} proteins by keeping longest isoform."
        log:
            config["results_dir"] + "/logs/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}/{genome}_{query}_summary.log"
        shell:
            """
            bedtools window -a {input.CDS_gff} -b {input.CDS_gff} -w 0 > {params.output_dir}/tmp_{wildcards.query}.txt
            python scripts/listLongestAlternativeTranscript.py {params.output_dir}/tmp_{wildcards.query}.txt {input.CDS_ALL_fna} > {params.output_dir}/tmp2_{wildcards.query}.txt
            python scripts/ExcludeCDSFromGFF.py {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS.gff {params.output_dir}/tmp2_{wildcards.query}.txt > {output.CDS_clean_gff}
            rm {params.output_dir}/tmp2_{wildcards.query}.txt {params.output_dir}/tmp_{wildcards.query}.txt
            """
    rule miniprot_concate:
        input:
            CDS_clean_gff=expand(config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{{genome}}/{{genome}}_{query}_CDS_clean.gff", query=QUERIES),
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
            CDS_ALL_fna=config["results_dir"] + "/annotation/Miniprot/Miniprot_{genome}/{genome}_all_CDS_clean.fna",
        output:
            TRANSCRIPT_final_rename_cat_fna=config["results_dir"] + "/annotation/Miniprot/Miniprot_{genome}/{genome}_final.fna",
            TRANSCRIPT_final_rename_cat_faa=config["results_dir"] + "/annotation/Miniprot/Miniprot_{genome}/{genome}_final.faa"
        params:
            output_dir = config["results_dir"] + "/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}",
        threads: 1
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        message:
            "Remove isoform in {wildcards.genome} Miniprot annotation from concatenation of all queries proteins by keeping longest isoform. Write fasta NT and AA for kept isoforms"
        log:
            config["results_dir"] + "/logs/annotation/Miniprot/Miniprot_clean/gfftofasta_{genome}/{genome}_summary.log"
        shell:
            """
            # concatenation of the gff files for each query
            mkdir -p {params.output_dir}
            cat {input.CDS_clean_gff} | bedtools sort > {params.output_dir}/{wildcards.genome}_all_CDS.gff

            # remove CDS include in an other CDS
            bedtools window -a {params.output_dir}/{wildcards.genome}_all_CDS.gff -b {params.output_dir}/{wildcards.genome}_all_CDS.gff -w 0 > {params.output_dir}/tmp_cat.txt
            python scripts/listLongestAlternativeTranscript.py {params.output_dir}/tmp_cat.txt {input.CDS_ALL_fna} > {params.output_dir}/tmp2_cat.txt
            python scripts/ExcludeCDSFromGFF.py {params.output_dir}/{wildcards.genome}_all_CDS.gff {params.output_dir}/tmp2_cat.txt > {params.output_dir}/{wildcards.genome}_all_CDS_clean.gff
            rm {params.output_dir}/tmp2_cat.txt {params.output_dir}/tmp_cat.txt
            
            # information retrieval
            echo "INFO: The GFF file contains $(cut -f9 {params.output_dir}/{wildcards.genome}_all_CDS.gff | cut -d\\; -f1 | sort -u | wc -l) cds" >> {log}
            echo "INFO: The cleaned GFF file contains $(cut -f9 {params.output_dir}/{wildcards.genome}_all_CDS_clean.gff | cut -d\\; -f1 | sort -u | wc -l) cds" >> {log}
            
            # transcript retrieval
            gffread -x {params.output_dir}/{wildcards.genome}_all_transcript.fna -g {input.assembly} {params.output_dir}/{wildcards.genome}_all_CDS_clean.gff &>> {log}
            sed -i "s/>/>{wildcards.genome}_Miniprot_/" {params.output_dir}/{wildcards.genome}_all_transcript.fna
            
            transeq -sequence {params.output_dir}/{wildcards.genome}_all_transcript.fna -outseq {params.output_dir}/{wildcards.genome}_all_transcript.faa &>> {log}
            
            # deletion of sequences containing stop codons (*)
            awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' {params.output_dir}/{wildcards.genome}_all_transcript.faa |tail -n +2 | paste - - | awk -F "\\t" '{{if($2 !~ /\\*[A-Z*]/){{sub(/\\*$/, "", $2); print $1"\\t"$2}}}}' | tr "\\t" "\\n" > {output.TRANSCRIPT_final_rename_cat_faa}
            
            sed -i 's/_1$//' {output.TRANSCRIPT_final_rename_cat_faa}
            
            # generate fna from faa
            grep "^>" {output.TRANSCRIPT_final_rename_cat_faa} | sed 's/^>//' | sort > {params.output_dir}/list_transcript.txt
            seqkit grep -f {params.output_dir}/list_transcript.txt {params.output_dir}/{wildcards.genome}_all_transcript.fna -w 0 -o {output.TRANSCRIPT_final_rename_cat_fna} &>> {log}
            """
if "Scipio" in ANNOTATIONS:
    # SCIPIO
    # rule for generating yaml files using the scipio annotation program
    rule scipio_yaml:
        input:
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
            query=config["results_dir"] + "/Queries/{query}.faa"
        output:
            yaml=config["results_dir"] + "/annotation/Scipio/Scipio_{query}/{genome}.yaml"
        params:
            blat_psl=config["results_dir"] + "/annotation/Scipio/Scipio_{query}/{genome}.psl"
        log:
            config["results_dir"] + "/logs/annotation/Scipio/Scipio_{query}/{genome}.log"
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        threads:
            max(int(workflow.cores / 5),1)
        message:
            "Protein alignment of {wildcards.query} query with Scipio to annotate {wildcards.genome}"
        shell:
            """
            set +e
            scipio.1.4.1.pl --verbose --blat_output {params.blat_psl} {input.assembly} {input.query} > {output.yaml} 2> {log}
            exitcode=$? # captures the value of the exit code for the scipio command

            # Scipio can exitcode with a value of 1 whithout reason. Then, the exit code is modify based on the output yaml.
            # if yaml output is empty exit with error else continue correctly
            if [[ -s {output.yaml} ]]; then 
                exit 0
            else
                echo "The yaml is empty you must had an error with scipio go check the log file."
                exit 1
            fi
            """
    # rule for generating gff files using the Perl script yaml2gff.1.4.pl
    rule scipio_gff:
        input:
            yaml = rules.scipio_yaml.output.yaml,
        output:
            gff=config["results_dir"] + "/annotation/Scipio/Scipio_{query}/{genome}.gff"
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        message:
            "Convert to gff Scipio alignment of {wildcards.genome} from {wildcards.query} proteins"
        log:
            config["results_dir"] + "/logs/annotation/Scipio/Scipio_{query}/{genome}_gff.log"
        threads:
            max(int(workflow.cores / 5),1)
        shell:
            """
            yaml2gff.1.4.pl {input.yaml} > {output.gff} 2> {log}
            """
    # rule for cleaning gff files from the scipio annotation program
    rule scipio_cds:
        input:
            gff=rules.scipio_gff.output.gff,
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
        output:
            CDS_gff=config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{genome}/{genome}_{query}_CDS.gff",
            CDS_fna=config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{genome}/{genome}_{query}_CDS.fna",
        params:
            output_dir = config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{genome}",
            minlen = config["minLenCDS"],
        threads: 1
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        log:
            config["results_dir"] + "/logs/annotation/Scipio/Scipio_clean/gfftofasta_{genome}/{genome}_{query}_summary.log"
        message:
            "Extract {wildcards.genome} CDS from gff annotation file of {wildcards.query} query"
        shell: 
            """
            # extract CDS and rename gene
            grep -P "\\tprotein_match\\t" {input.gff} | sed -E "s/ID=([0-9]+);Query/Parent=SC\\1_{wildcards.query};Query/" | sed "s/protein_match/CDS/g" |awk -F "\\t" '{{if($4 <= $5){{print}}}}' | bedtools sort > {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp1.gff
            
            # extract CDS and rename gene in fna
            gffread -l {params.minlen} -o {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp2.gff -x {output.CDS_fna} -g {input.assembly} {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp1.gff &>> {log}
            grep -P "\\tCDS\\t" {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp2.gff > {output.CDS_gff}
            rm {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp2.gff {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS_tmp1.gff
            """
    rule scipio_concate_cds_fna:
        input:
            CDS_fna=expand(config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{{genome}}/{{genome}}_{query}_CDS.fna", query=QUERIES),
        output:
            CDS_ALL_fna=config["results_dir"] + "/annotation/Scipio/Scipio_{genome}/{genome}_all_CDS_clean.fna",
        threads: 1
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        message:
            "Concatenate {wildcards.genome} fasta Scipio annotations from all queries"
        shell:
            """
            cat {input.CDS_fna} | cut -d " " -f 1 > {output.CDS_ALL_fna}
            """
    rule scipio_clean_cds:
        input:
            CDS_gff=config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{genome}/{genome}_{query}_CDS.gff",
            CDS_ALL_fna=config["results_dir"] + "/annotation/Scipio/Scipio_{genome}/{genome}_all_CDS_clean.fna",
        output:
            CDS_clean_gff=config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{genome}/{genome}_{query}_CDS_clean.gff",
        params:
            output_dir = config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{genome}",
        threads: 1
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        log:
            config["results_dir"] + "/logs/annotation/Scipio/Scipio_clean/gfftofasta_{genome}/{genome}_{query}_summary.log"
        message:
            "Remove isoform in {wildcards.genome} Scipio annotation from {wildcards.query} proteins by keeping longest isoform. Write fasta NT and AA for kept isoforms"
        shell:
            """
            bedtools window -a {input.CDS_gff} -b {input.CDS_gff} -w 0 > {params.output_dir}/tmp_{wildcards.query}.txt
            python scripts/listLongestAlternativeTranscript.py {params.output_dir}/tmp_{wildcards.query}.txt {input.CDS_ALL_fna} > {params.output_dir}/tmp2_{wildcards.query}.txt
            python scripts/ExcludeCDSFromGFF.py {params.output_dir}/{wildcards.genome}_{wildcards.query}_CDS.gff {params.output_dir}/tmp2_{wildcards.query}.txt > {output.CDS_clean_gff}
            rm {params.output_dir}/tmp2_{wildcards.query}.txt {params.output_dir}/tmp_{wildcards.query}.txt
            """
    # rule for concatenating gff files, cleans up the concatenated gff file and converts it to fna and faa format
    rule scipio_concate:
        input:
            CDS_clean_gff=expand(config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{{genome}}/{{genome}}_{query}_CDS_clean.gff", query=QUERIES),
            assembly=config["results_dir"] + "/Genomes/{genome}.fasta",
            CDS_ALL_fna=config["results_dir"] + "/annotation/Scipio/Scipio_{genome}/{genome}_all_CDS_clean.fna",
        output:
            TRANSCRIPT_final_rename_cat_fna=config["results_dir"] + "/annotation/Scipio/Scipio_{genome}/{genome}_final.fna",
            TRANSCRIPT_final_rename_cat_faa=config["results_dir"] + "/annotation/Scipio/Scipio_{genome}/{genome}_final.faa",
        params:
            output_dir = config["results_dir"] + "/annotation/Scipio/Scipio_clean/gfftofasta_{genome}",
        threads: 1
        singularity:
            "Container/Miniprot_Scipio/miniprot_scipio.sif"
        log:
            config["results_dir"] + "/logs/annotation/Scipio/Scipio_clean/gfftofasta_{genome}/{genome}_summary.log"
        message:
            "Remove isoform in {wildcards.genome} Scipio annotation from concatenation of all queries proteins by keeping longest isoform. Write fasta NT and AA for kept isoforms"
        shell:
            """
            # concatenation of the gff files for each query
            mkdir -p {params.output_dir}
            cat {input.CDS_clean_gff} | bedtools sort > {params.output_dir}/{wildcards.genome}_all_CDS.gff
            
            # remove CDS include in an other CDS
            bedtools window -a {params.output_dir}/{wildcards.genome}_all_CDS.gff -b {params.output_dir}/{wildcards.genome}_all_CDS.gff -w 0 > {params.output_dir}/tmp_cat.txt
            python scripts/listLongestAlternativeTranscript.py {params.output_dir}/tmp_cat.txt {input.CDS_ALL_fna} > {params.output_dir}/tmp2_cat.txt
            python scripts/ExcludeCDSFromGFF.py {params.output_dir}/{wildcards.genome}_all_CDS.gff {params.output_dir}/tmp2_cat.txt > {params.output_dir}/{wildcards.genome}_all_CDS_clean.gff
            rm {params.output_dir}/tmp2_cat.txt {params.output_dir}/tmp_cat.txt
            
            # information retrieval
            echo "INFO: The GFF file contains $(cut -f9 {params.output_dir}/{wildcards.genome}_all_CDS.gff | cut -d\\; -f1 | sort -u | wc -l) cds" >> {log}
            echo "INFO: The cleaned GFF file contains $(cut -f9 {params.output_dir}/{wildcards.genome}_all_CDS_clean.gff | cut -d\\; -f1 | sort -u | wc -l) cds" >> {log}
            
            # transcript retrieval
            gffread -x {params.output_dir}/{wildcards.genome}_all_transcript.fna -g {input.assembly} {params.output_dir}/{wildcards.genome}_all_CDS_clean.gff &>>{log}
            sed -i "s/>/>{wildcards.genome}_Scipio_/" {params.output_dir}/{wildcards.genome}_all_transcript.fna
            
            transeq -sequence {params.output_dir}/{wildcards.genome}_all_transcript.fna -outseq {params.output_dir}/{wildcards.genome}_all_transcript.faa &>>{log}
            
            # deletion of sequences containing stop codons (*)
            awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' {params.output_dir}/{wildcards.genome}_all_transcript.faa |tail -n +2 | paste - - | awk -F "\\t" '{{if($2 !~ /\\*[A-Z*]/){{sub(/\\*$/, "", $2); print $1"\\t"$2}}}}' | tr "\\t" "\\n" > {output.TRANSCRIPT_final_rename_cat_faa}
            sed -i 's/_1$//' {output.TRANSCRIPT_final_rename_cat_faa}
            
            # generate fna from faa
            grep "^>" {output.TRANSCRIPT_final_rename_cat_faa} | sed 's/^>//' | sort > {params.output_dir}/list_transcript.txt
            seqkit grep -f {params.output_dir}/list_transcript.txt {params.output_dir}/{wildcards.genome}_all_transcript.fna -w 0 -o {output.TRANSCRIPT_final_rename_cat_fna} &>> {log}
            """

### End of annotation part
### Begin of orthology part

rule group_genomefaa_orthofinder:
    input:
        faa=config["results_dir"] + "/annotation/{annotation}/{annotation}_{genome}/{genome}_final.faa",
    output:
        output_faa = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_faa/{genome}.faa",
    singularity:
        "Container/Cleaning/cleaning_orthology.sif"
    message:
        "Prepare orthofinder by grouping AA fasta of all genome and queries in a same directory"
    shell:
        """
        ln -sr {input.faa} {output.output_faa}
        """


rule group_queriesfaa_orthofinder:
    input:
        query=config["results_dir"] + "/Queries/{query}.faa"
    output:
        output_annot = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_faa/{query}.faa",
    message:
        "Prepare orthofinder by grouping AA fasta of all genome and queries in a same directory"
    shell:
        """
        ln -sr {input.query} {output.output_annot}
        """

rule group_genomefna_orthofinder:
    input:
        fna=config["results_dir"] + "/annotation/{annotation}/{annotation}_{genome}/{genome}_final.fna",
    output:
        output_fna = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_fna/{genome}.fna",
    singularity:
        "Container/Cleaning/cleaning_orthology.sif"
    message:
        "Prepare orthofinder by grouping NT fasta of all genome and queries in a same directory"
    shell:
        """
        ln -sr {input.fna} {output.output_fna}
        """


rule group_queriesfna_orthofinder:
    input:
        query=config["results_dir"] + "/Queries/{query}.fna"
    output:
        output_annot = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_fna/{query}.fna",
    message:
        "Prepare orthofinder by grouping NT fasta of all genome and queries in a same directory"
    shell:
        """
        ln -sr {input.query} {output.output_annot}
        """


rule orthofinder:
    input:
        genome_faa = expand(config["results_dir"] + "/orthology/orthofinder/orthof_{{annotation}}/orthof_{{annotation}}_faa/{genome}.faa", genome=GENOMES),
        query_faa = expand(config["results_dir"] + "/orthology/orthofinder/orthof_{{annotation}}/orthof_{{annotation}}_faa/{query}.faa", query=QUERIES),
    output:
        genecount_ortholog = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_results/Orthogroups/Orthogroups.GeneCount.tsv",
        genelist_ortholog = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_results/Orthogroups/Orthogroups.tsv",
    params:
        outdir = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/Orthofinder/",
        resdir= config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_results/Orthogroups/",
        tempo = "orthof_{annotation}",
        dir_faa = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_faa/",
        orthofinder_method = config["orthofinder_method"],
    threads:
        max(int(workflow.cores / len(ANNOTATIONS)),1)
    message:
        "Identify orthologs and paralogs between all genomes and queries in {wildcards.annotation} with orthofinder"
    singularity:
        "Container/Cleaning/cleaning_orthology.sif"
    log:
        config["results_dir"] + "/logs/orthology/orthofinder/orthof_{{annotation}}/{{annotation}}.log"
    shell:
        """
        rm -f /tmp/test.ckp.gz
        rm -fr {params.outdir}/
        if [ "{wildcards.annotation}" = "Braker" ];then
            rm -rf RM_*
        fi
        orthofinder -f {params.dir_faa} -M msa -T {params.orthofinder_method} -o {params.outdir} -n {params.tempo} -t {threads} &>> {log}
        mkdir -p {params.resdir}
        mv {params.outdir}/Results_{params.tempo}/Orthogroups/* {params.resdir}
        """

### End to identify orthogroup
### Begin cleaning orthogroups


if config["orthogroup_gene"]=="ALL":        ## Version search orthogroup 1:1 in Query + Genome
    checkpoint select_ortholog:
        input:
            genecount_ortholog = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_results/Orthogroups/Orthogroups.GeneCount.tsv",
            genelist_ortholog = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_results/Orthogroups/Orthogroups.tsv",
        output:
            dir_og_fna = directory(config["results_dir"] + "/alignment/align_{annotation}/OG/Orthogroup_Sequences_fna/"),
            dir_og_faa = directory(config["results_dir"] + "/alignment/align_{annotation}/OG/Orthogroup_Sequences_faa/"),
            noparalog = config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_no_paralog.txt",
            oneparalog = config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_one_paralog.txt",
            oneparalog_sister = config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_one_paralog_sister.txt",
            og_list = config["results_dir"] + "/alignment/align_{annotation}/List_all_orthogroups_kept_{annotation}.txt",
        params:
            min_orthogroup_threshold = config["min_sp_ortholog"],
            res_dir=config["results_dir"] + "/alignment/align_{annotation}/",
            queries=QUERIES,
            genome=GENOMES,
            orthof_dir = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/Orthofinder/",
            tempo = "orthof_{annotation}",
            dir_genome_fna = directory(config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_fna/"),
            dir_genome_faa = directory(config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_faa/"),
        threads:
            max(int(workflow.cores / len(ANNOTATIONS)),1)
        singularity:
            "Container/Cleaning/cleaning_orthology.sif"
        message:
            "Filter orthogroups in {wildcards.annotation} by keeping orthogroups where each genome and queries are represent by 0 or 1 genes (orthologous gene 1:1 or with only one paralog). Only orthogroups that contain gene from {params.min_orthogroup_threshold} genomes and at least 1 query are kept. For kept Orthogroups the NT fasta is build"
        log:
            config["results_dir"] + "/logs/alignment/align_{annotation}/{annotation}.log"
        shell:
            """
            # Report orthologous gene 1:1 or with only one paralog
            
            awk -F "\\t" -v genome="{params.genome}" -v queries="{params.queries}" -v min_orthogroup="{params.min_orthogroup_threshold}" 'BEGIN{{
                max=1
                ngenome=split(genome,g," ")
                nqueries=split(queries,q," ")
                if(min_orthogroup<1){{
                    min_n_orthogroup=int((ngenome+nqueries)*min_orthogroup)
                    min_n_orthogroup = (max > min_n_orthogroup ? max : min_n_orthogroup)
                }}else{{
                    min_n_orthogroup=min_orthogroup
                }}
            }}
            {{
                if(NR==1){{
                    for(i=2;i<NF;i++){{
                        for (s in q){{
                            if($i == q[s]){{
                                toskip[i]==1
                            }}
                        }}
                    }}
                }}else{{
                    cnt=0
                    cnt_query=0
                    paralog=0
                    for(i=2;i<NF;i++){{
                        if (i in toskip == 0){{
                            if($i>2){{
                                print "At least a genome or a query have more than 2 genes ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                                next;
                            }}else if($i==2){{
                                paralog+=1
                                cnt+=1
                            }}else if($i==1){{
                                cnt+=1
                            }}
                        }}else{{
                            if($i>2){{
                                print "At least a genome or a query have more than 2 genes ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                                next;
                            }}else if($i==2){{
                                paralog+=1
                                cnt_query+=1
                            }}else if($i==1){{
                                cnt_query+=1
                            }}
                        }}
                    }}
                    if(cnt_query>0){{
                        if(paralog==0 && cnt+cnt_query>=min_n_orthogroup){{
                            print $1 >"{output.noparalog}"
                            print
                        }}else if(paralog==1 && cnt+cnt_query>=min_n_orthogroup){{
                            print $1 >"{output.oneparalog}"
                            print
                        }}else if(cnt+cnt_query<min_n_orthogroup){{
                            print "Not enough genome in orthogroup "cnt" < "min_n_orthogroup" ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                        }}else if(paralog>1){{
                            print "Too many paralogs : "paralog" genomes with a count equal to 2 ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                        }}
                    }}else{{
                        print "No gene in queries ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                    }}
                }}
            }}' {input.genecount_ortholog} > {params.res_dir}/GeneCount_all_orthogroups_kept.txt
            
            # Write all genome species in a file
            echo "{params.genome}" > {params.res_dir}/List_Genome.txt
            python scripts/SelectOrthogroupNamesOneParalogSister.py {output.oneparalog} {params.orthof_dir}/Results_{params.tempo}/Resolved_Gene_Trees {params.res_dir}/List_Genome.txt > {output.oneparalog_sister}
            
            
            # Write nucleotidic and proteic sequence of kept orthogroup.
            mkdir -p {params.res_dir}/OG/Orthogroup_Sequences_faa/
            
            # Keep all OG sequence in one directory
            for OG in $(cat {output.noparalog})
            do 
            ln -sr {params.orthof_dir}/Results_{params.tempo}/Orthogroup_Sequences/${{OG}}.fa {params.res_dir}/OG/Orthogroup_Sequences_faa/${{OG}}.faa; 
            done
            
            for OG in $(cat {output.oneparalog_sister})
            do
            ### Find sequence ID to remove (from paralogous species)
            grep $OG {input.genelist_ortholog} | awk -F "\\t" '{{for(i=2;i<=NF;i++){{n=split($i,a," ");if(n>1){{for(j=1;j<=n;j++){{print a[j]}}}}}}}}' | sed 's/,$//' > {params.orthof_dir}/Results_{params.tempo}/pattern_${{OG}}.txt
            ### Remove seqID
            awk '{{ if ((NR>1)&&($0~/^>/)) {{ printf("\\n%s", $0); }} else if (NR==1) {{ printf("%s", $0); }} else {{ printf("\\t%s", $0); }} }}' {params.orthof_dir}/Results_{params.tempo}/Orthogroup_Sequences/${{OG}}.fa | grep -v -Ff {params.orthof_dir}/Results_{params.tempo}/pattern_${{OG}}.txt - | tr "\\t" "\\n" > {params.res_dir}/OG/Orthogroup_Sequences_faa/${{OG}}.faa
            done
            
            mkdir -p {params.res_dir}/OG/Orthogroup_Sequences_fna/
            parallel -j {threads} --group 'echo {{}} >> {log}; grep "^>" {params.res_dir}/OG/Orthogroup_Sequences_faa/{{}}.faa | sed "s/^>//" > {params.res_dir}/OG/list_id_to_catch_{{}}.txt; seqkit grep -f {params.res_dir}/OG/list_id_to_catch_{{}}.txt {params.dir_genome_fna}/*.fna -o {params.res_dir}/OG/Orthogroup_Sequences_fna/{{}}.fna &>> {log}; rm {params.res_dir}/OG/list_id_to_catch_{{}}.txt' ::: $(find -L {params.res_dir}/OG/Orthogroup_Sequences_faa/ -type f | awk -F "/" '{{print $NF}}' |sed 's/.faa$//')
            
            # Build a file with list of OG kept for the analysis
            cat {output.noparalog} {output.oneparalog_sister} > {output.og_list}
            """
elif config["orthogroup_gene"]=="GENOME":        ## Version search orthogroup 1:1 in Genome only
    checkpoint select_ortholog:
        input:
            genecount_ortholog = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_results/Orthogroups/Orthogroups.GeneCount.tsv",
            genelist_ortholog = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_results/Orthogroups/Orthogroups.tsv",
        output:
            dir_og_fna = directory(config["results_dir"] + "/alignment/align_{annotation}/OG/Orthogroup_Sequences_fna/"),
            dir_og_faa = directory(config["results_dir"] + "/alignment/align_{annotation}/OG/Orthogroup_Sequences_faa/"),
            noparalog = config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_no_paralog.txt",
            oneparalog = config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_one_paralog.txt",
            oneparalog_sister = config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_one_paralog_sister.txt",
            og_list = config["results_dir"] + "/alignment/align_{annotation}/List_all_orthogroups_kept_{annotation}.txt",
        params:
            min_orthogroup_threshold = config["min_sp_ortholog"],
            res_dir=config["results_dir"] + "/alignment/align_{annotation}/",
            queries=QUERIES,
            genome=GENOMES,
            orthof_dir = config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/Orthofinder/",
            tempo = "orthof_{annotation}",
            dir_genome_fna = directory(config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_fna/"),
            dir_genome_faa = directory(config["results_dir"] + "/orthology/orthofinder/orthof_{annotation}/orthof_{annotation}_faa/"),
        threads:
            max(int(workflow.cores / len(ANNOTATIONS)),1)
        singularity:
            "Container/Cleaning/cleaning_orthology.sif"
        message:
            "Filter orthogroups in {wildcards.annotation} by keeping orthogroups where each genome (not queries) are represent by 0 or 1 genes (orthologous gene 1:1 or with only one paralog). Only orthogroups that contain gene from {params.min_orthogroup_threshold} genomes and at least 1 query are kept. For kept Orthogroups the NT fasta is build"
        log:
            config["results_dir"] + "/logs/alignment/align_{annotation}/{annotation}.log"
        shell:
            """
            # Report orthologous gene 1:1 or with only one paralog
            
            awk -F "\\t" -v genome="{params.genome}" -v queries="{params.queries}" -v min_orthogroup="{params.min_orthogroup_threshold}" 'BEGIN{{
                max=1
                ngenome=split(genome,g," ")
                nqueries=split(queries,q," ")
                if(min_orthogroup<1){{
                    min_n_orthogroup=int(ngenome*min_orthogroup)
                    min_n_orthogroup = (max > min_n_orthogroup ? max : min_n_orthogroup)
                }}else{{
                    min_n_orthogroup=min_orthogroup
                }}
            }}
            {{
                if(NR==1){{
                    for(i=2;i<NF;i++){{
                        for (s in q){{
                            if($i == q[s]){{
                                toskip[i]==1
                            }}
                        }}
                    }}
                }}else{{
                    cnt=0
                    cnt_query=0
                    paralog=0
                    for(i=2;i<NF;i++){{
                        if (i in toskip == 0){{
                            if($i>2){{
                                print "At least a genome have more than 2 genes ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                                next;
                            }}else if($i==2){{
                                paralog+=1
                                cnt+=1
                            }}else if($i==1){{
                                cnt+=1
                            }}
                        }}else{{
                            if($i>0){{
                                cnt_query+=1
                            }}
                        }}
                    }}
                    if(cnt_query>0){{
                        if(paralog==0 && cnt>=min_n_orthogroup){{
                            print $1 >"{output.noparalog}"
                            print
                        }}else if(paralog==1 && cnt>=min_n_orthogroup){{
                            print $1 >"{output.oneparalog}"
                            print
                        }}else if(cnt<min_n_orthogroup){{
                            print "Not enough genome in orthogroup "cnt" < "min_n_orthogroup" ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                        }}else if(paralog>1){{
                            print "Too many paralogs : "paralog" genomes with a count equal to 2 ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                        }}
                    }}else{{
                        print "No gene in queries ->\\t"$0 >"{params.res_dir}/GeneCount_NotKept_orthogroups.txt"
                    }}
                }}
            }}' {input.genecount_ortholog} > {params.res_dir}/GeneCount_all_orthogroups_kept.txt
            
            # Write all genome species in a file
            echo "{params.genome}" > {params.res_dir}/List_Genome.txt
            python scripts/SelectOrthogroupNamesOneParalogSister.py {output.oneparalog} {params.orthof_dir}/Results_{params.tempo}/Resolved_Gene_Trees {params.res_dir}/List_Genome.txt > {output.oneparalog_sister}
            
            
            # Write nucleotidic and proteic sequence of kept orthogroup.
            mkdir -p {params.res_dir}/OG/Orthogroup_Sequences_faa/
            
            # Keep all OG sequence in one directory
            for OG in $(cat {output.noparalog})
            do 
            ln -sr {params.orthof_dir}/Results_{params.tempo}/Orthogroup_Sequences/${{OG}}.fa {params.res_dir}/OG/Orthogroup_Sequences_faa/${{OG}}.faa; 
            done
            
            for OG in $(cat {output.oneparalog_sister})
            do
            ### Find sequence ID to remove (from paralogous species)
            grep $OG {input.genelist_ortholog} | awk -F "\\t" '{{for(i=2;i<=NF;i++){{n=split($i,a," ");if(n>1){{for(j=1;j<=n;j++){{print a[j]}}}}}}}}' | sed 's/,$//' > {params.orthof_dir}/Results_{params.tempo}/pattern_${{OG}}.txt
            ### Remove seqID
            awk '{{ if ((NR>1)&&($0~/^>/)) {{ printf("\\n%s", $0); }} else if (NR==1) {{ printf("%s", $0); }} else {{ printf("\\t%s", $0); }} }}' {params.orthof_dir}/Results_{params.tempo}/Orthogroup_Sequences/${{OG}}.fa | grep -v -Ff {params.orthof_dir}/Results_{params.tempo}/pattern_${{OG}}.txt - | tr "\\t" "\\n" > {params.res_dir}/OG/Orthogroup_Sequences_faa/${{OG}}.faa
            done
            
            mkdir -p {params.res_dir}/OG/Orthogroup_Sequences_fna/
            parallel -j {threads} --group 'echo {{}} >> {log} ; grep "^>" {params.res_dir}/OG/Orthogroup_Sequences_faa/{{}}.faa | sed "s/^>//" > {params.res_dir}/OG/list_id_to_catch_{{}}.txt; seqkit grep -f {params.res_dir}/OG/list_id_to_catch_{{}}.txt {params.dir_genome_fna}/*.fna -o {params.res_dir}/OG/Orthogroup_Sequences_fna/{{}}.fna &>> {log}; rm {params.res_dir}/OG/list_id_to_catch_{{}}.txt' ::: $(find -L {params.res_dir}/OG/Orthogroup_Sequences_faa/ -type f | awk -F "/" '{{print $NF}}' |sed 's/.faa$//')
            
            # Build a file with list of OG kept for the analysis
            cat {output.noparalog} {output.oneparalog_sister} > {output.og_list}
            """



rule omm_macse:
    input:
        og_list = config["results_dir"] + "/alignment/align_{annotation}/List_all_orthogroups_kept_{annotation}.txt",
        OG_fna = config["results_dir"] + "/alignment/align_{annotation}/OG/Orthogroup_Sequences_fna/{og}.fna",
    output:
        OG_align = directory(config["results_dir"] + "/alignment/align_{annotation}/MACSE_clean/{og}_macse"),
    params:
        prefix = "{og}_macse",
        macse_dir = config["results_dir"] + "/alignment/align_{annotation}/MACSE_clean/",
    threads:
        1
    container:
        "library://vranwez/default/omm_macse:v12.01"
    message:
        "Align with omm_macse all orthogroups from {wildcards.annotation}"
    log:
        config["results_dir"] + "/logs/alignment/align_{annotation}/MACSE_clean/{og}_macse.log"
    shell:
        """
        mkdir -p {params.macse_dir}
        bash /OMM_MACSE/S_OMM_MACSE_V12.01.sh --out_dir {output.OG_align} --out_file_prefix {params.prefix} --in_seq_file {input.OG_fna} --genetic_code_number 1 --alignAA_soft MAFFT --java_mem 10000m &>> {log} || true
        """

def get_orthogroup_names(wildcards):
    checkpoint_output = checkpoints.select_ortholog.get(**wildcards).output[0]
    dir_fna = expand(config["results_dir"] + "/alignment/align_{annotation}/MACSE_clean/{og}_macse", annotation=wildcards.annotation, og=glob_wildcards(os.path.join(checkpoint_output, "{og}.fna")).og)
    return dir_fna



rule omm_macse_list:
    input:
        get_orthogroup_names
    output:
        list_OG_align = config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_post_alignment.txt",
    threads:
        1
    message:
        "List orthogroup that pass alignment with omm_macse from {wildcards.annotation}"
    shell:
        """
        # Verify if MACSE work on the OG, write it in a list if yes.
        for d in {input}
        do
            nb_file=$(ls $d | wc -l)
            if [ "$nb_file" -ge 10 ]; then
                orthogroup=$(echo $d | awk -F "/" '{{print $NF}}' | sed 's/_macse//')
                if [ $(grep -c "^>" ${{d}}/${{orthogroup}}_macse_final_mask_align_NT.aln) -ge 3 ]; then
                    echo "$orthogroup" >> {output.list_OG_align};
                fi
            fi
        done
        """

checkpoint prepare_iqtree_phylter:
    input:
        list_OG_align = config["results_dir"] + "/alignment/align_{annotation}/List_orthogroup_names_post_alignment.txt",
    output:
        aln = directory(config["results_dir"] + "/cleaning/phylter_{annotation}/alignment_phylter_{annotation}/"),
    params:
        dir_results = config["results_dir"],
        queries=QUERIES,
        genome=GENOMES,
        list_Genome_Query = config["results_dir"] + "/List_Genome_Query.txt",
    message:
        "Prepare IQ-tree by renaming ID in align orthogroup fasta"
    shell:
        """
        echo "{params.queries} {params.genome}" | sed -e 's/ /\\n/g' > {params.list_Genome_Query}
        for OG_in_list in $(cat {input})
        do
        mkdir -p {output.aln}/${{OG_in_list}}_phylter/
        sed -i 's/!/-/g' {params.dir_results}/alignment/align_{wildcards.annotation}/MACSE_clean/${{OG_in_list}}_macse/${{OG_in_list}}_macse_final_mask_align_NT.aln 
        awk -F "\\t" 'fname != FILENAME {{fname = FILENAME; idx++}}
            idx == 1 {{list_genome[$1]=1}}
            idx == 2 {{if($0~/^>/){{for(g in list_genome){{if($0 ~ "^>"g){{print ">"g}}}}}}else{{print}}}}' {params.list_Genome_Query} {params.dir_results}/alignment/align_{wildcards.annotation}/MACSE_clean/${{OG_in_list}}_macse/${{OG_in_list}}_macse_final_mask_align_NT.aln > {output.aln}/${{OG_in_list}}_phylter/${{OG_in_list}}_final_mask_align_NT_rename.aln
        done
        """

rule iqtree:
    input:
        og_filter=config["results_dir"] + "/cleaning/phylter_{annotation}/alignment_phylter_{annotation}/{og_filter}_phylter/{og_filter}_final_mask_align_NT_rename.aln"
    output:
        tree=config["results_dir"] + "/cleaning/phylter_{annotation}/alignment_phylter_{annotation}/{og_filter}_phylter/{og_filter}_final_mask_align_NT_rename.aln.treefile"
    threads:
        1
    singularity:
        "Container/Cleaning/cleaning_orthology.sif"
    message:
        "Generate phylogenetic trees for othogroups from {wildcards.annotation}"
    log:
        config["results_dir"] + "/logs/cleaning/phylter_{annotation}/alignment_phylter_{annotation}/{og_filter}_phylter/{og_filter}_iqtree.log"
    shell:
        """
        iqtree2 -m GTR+G4 -s {input.og_filter} -T {threads} &>> {log}
        """


def get_og_filter(wildcards):
    checkpoint_output = checkpoints.prepare_iqtree_phylter.get(**wildcards).output[0]
    iqtree_out = expand(config["results_dir"] + "/cleaning/phylter_{annotation}/alignment_phylter_{annotation}/{og_filter}_phylter/{og_filter}_final_mask_align_NT_rename.aln.treefile", annotation=wildcards.annotation, og_filter=glob_wildcards(os.path.join(checkpoint_output, "{og_filter}_phylter/{og_filter}_final_mask_align_NT_rename.aln")).og_filter)
    return iqtree_out

rule phylter:
    input:
        get_og_filter,
    output:
        tre=config["results_dir"] + "/cleaning/phylter_{annotation}/treefiles_{annotation}.tre",
        list_tre=config["results_dir"] + "/cleaning/phylter_{annotation}/list_treefiles_{annotation}.txt",
        phylter=config["results_dir"] + "/cleaning/phylter_{annotation}/phylter_seqstats_{annotation}.out",
        phylter_outlier=config["results_dir"] + "/cleaning/phylter_{annotation}/phylter_outliers_{annotation}.txt",
    params:
        output_dir=config["results_dir"] + "/cleaning/phylter_{annotation}/",
    singularity:
        "Container/Cleaning/cleaning_orthology.sif"
    log:
        config["results_dir"] + "/logs/cleaning/phylter_{annotation}/phylter_outliers_{annotation}.log"
    message:
        "Define outliers in orthogroup from {wildcards.annotation}"
    shell:
        """
        for t in {input}
        do
            if [ -f "$t" ]; then
                cat "$t" >> {output.tre}
                echo "$t" | awk -F "/" '{{print $(NF-1)}}' | sed 's/_phylter//' >> {output.list_tre}
            fi
        done
        Rscript scripts/Script_PhylteR.r {output.tre} {output.list_tre} {output.phylter} &>> {log}
        
        # Then create for each OG the list of species outliers
        # END write OG with all the identified outlier.
        awk -F "\\t" '{{
            if($0 !~ /^#/){{
                if($1 in og2outlier){{
                    og2outlier[$1]=og2outlier[$1]"|"$2
                }}else{{
                    og2outlier[$1]=$2
                }}
            }}
        }}END{{
            for(o in og2outlier){{
                print o"\\t"og2outlier[o]
            }}
        }}' {output.phylter} > {output.phylter_outlier}
        
        """


if config["orthogroup_gene"]=="ALL":
    checkpoint remove_outgroup_filter_og:
        input:
            phylter_outlier=config["results_dir"] + "/cleaning/phylter_{annotation}/phylter_outliers_{annotation}.txt",
        output:
            pruned_aln=directory(config["results_dir"] +"/comparison/comparison_{annotation}/"),
            og_to_rerun=directory(config["results_dir"] + "/cleaning/phylter_{annotation}/OG_ToReRun/NT/"),
        params:
            dir_init_name = config["results_dir"] + "/alignment/align_{annotation}/MACSE_clean/",
            dir_og_phylter= config["results_dir"] + "/alignment/align_{annotation}/OG/Orthogroup_Sequences_fna/",
            dir_res_phylter= config["results_dir"] + "/cleaning/phylter_{annotation}/alignment_phylter_{annotation}/",
            min_orthogroup_threshold = config["min_sp_ortholog"],
            queries=QUERIES,
            genome=GENOMES,
        container:
            "library://vranwez/default/omm_macse:v12.01"
        message:
            "Keep orthogroups without outliers and prepare rerun of omm_macse on orthologs with outliers by filtering ID defined as outliers"
        shell:
            """
            mkdir -p {output.pruned_aln}/NT
            mkdir -p {output.pruned_aln}/AA
            mkdir -p {output.og_to_rerun}
            # BEGIN : Define minimum orthogroup threshold
            threshold=$(awk -F "\\t" -v genome="{params.genome}" -v queries="{params.queries}" -v min_orthogroup="{params.min_orthogroup_threshold}" 'BEGIN{{
                max=1
                ngenome=split(genome,g," ")
                nqueries=split(queries,q," ")
                if(min_orthogroup<1){{
                    min_n_orthogroup=int((ngenome+nqueries)*min_orthogroup)
                    min_n_orthogroup = (max > min_n_orthogroup ? max : min_n_orthogroup)
                }}else{{
                    min_n_orthogroup=min_orthogroup
                }}
                print min_n_orthogroup
            }}')
            echo "The threshold value is $threshold"
            
            for og_file in $(find {params.dir_init_name} -type f -name "*_final_align_NT.aln")
            do
            og_file_faa=$(echo $og_file | sed 's/_NT.aln$/_AA.aln/')
            og=$(echo $og_file | awk -F "/" '{{print $(NF-1)}}' | sed 's/_macse//')
            if grep -w -q $og {input.phylter_outlier}; then
                pattern=$(grep -w $og {input.phylter_outlier} | cut -f 2 | sed 's/|/|>/g' | sed 's/^/>/')
                queries=$(echo "^>{params.queries}" |sed 's/ /|^>/g')
                n=$(grep -v -E "$pattern" $og_file | grep -c "^>" ||true)
                if [ "$n" -lt "$threshold" ]; then
                    continue
                else
                    perl -pe '$. > 1 and /^>/ ? print "\\n" : chomp' {params.dir_og_phylter}/${{og}}.fna | paste - - | grep -v -E "$pattern" | tr "\\t" "\\n" > {output.og_to_rerun}/${{og}}_tokeep_phylter.fna
                fi
            else
                cp $og_file {output.pruned_aln}/NT/${{og}}_macse_final_mask_align_NT_pruned_complet.aln
                cp {params.dir_init_name}/${{og}}_macse/${{og}}_macse_final_unmask_align_AA.aln {output.pruned_aln}/AA/${{og}}_macse_final_unmask_align_AA_pruned_complet.aln
            fi
            done
            """
elif config["orthogroup_gene"]=="GENOME":
    checkpoint remove_outgroup_filter_og:
        input:
            phylter_outlier=config["results_dir"] + "/cleaning/phylter_{annotation}/phylter_outliers_{annotation}.txt",
        output:
            pruned_aln=directory(config["results_dir"] +"/comparison/comparison_{annotation}/"),
            og_to_rerun=directory(config["results_dir"] + "/cleaning/phylter_{annotation}/OG_ToReRun/NT/"),
        params:
            dir_init_name = config["results_dir"] + "/alignment/align_{annotation}/MACSE_clean/",
            dir_og_phylter= config["results_dir"] + "/alignment/align_{annotation}/OG/Orthogroup_Sequences_fna/",
            dir_res_phylter= config["results_dir"] + "/cleaning/phylter_{annotation}/alignment_phylter_{annotation}/",
            min_orthogroup_threshold = config["min_sp_ortholog"],
            queries=QUERIES,
            genome=GENOMES,
        container:
            "library://vranwez/default/omm_macse:v12.01"
        message:
            "Keep orthogroups without outliers and prepare rerun of omm_macse on orthologs with outliers by filtering ID defined as outliers"
        shell:
            """
            mkdir -p {output.pruned_aln}/NT
            mkdir -p {output.pruned_aln}/AA
            mkdir -p {output.og_to_rerun}
            # BEGIN : Define minimum orthogroup threshold
            threshold=$(awk -F "\\t" -v genome="{params.genome}" -v min_orthogroup="{params.min_orthogroup_threshold}" 'BEGIN{{
                max=1
                ngenome=split(genome,g," ")
                if(min_orthogroup<1){{
                    min_n_orthogroup=int(ngenome*min_orthogroup)
                    min_n_orthogroup = (max > min_n_orthogroup ? max : min_n_orthogroup)
                }}else{{
                    min_n_orthogroup=min_orthogroup
                }}
                print min_n_orthogroup
            }}')
            echo "The threshold value is $threshold"
            
            for og_file in $(find {params.dir_init_name} -type f -name "*_final_align_NT.aln")
            do
            og_file_faa=$(echo $og_file | sed 's/_NT.aln$/_AA.aln/')
            og=$(echo $og_file | awk -F "/" '{{print $(NF-1)}}' | sed 's/_macse//')
            if grep -w -q $og {input.phylter_outlier}; then
                pattern=$(grep -w $og {input.phylter_outlier} | cut -f 2 | sed 's/|/|>/g' | sed 's/^/>/')
                queries=$(echo "^>{params.queries}" |sed 's/ /|^>/g')
                n=$(grep -v -E "$pattern" $og_file | grep -v -E "$queries" | grep -c "^>" ||true)
                if [ "$n" -lt "$threshold" ]; then
                    continue
                else
                    perl -pe '$. > 1 and /^>/ ? print "\\n" : chomp' {params.dir_og_phylter}/${{og}}.fna | paste - - | grep -v -E "$pattern" | tr "\\t" "\\n" > {output.og_to_rerun}/${{og}}_tokeep_phylter.fna
                fi
            else
                cp $og_file {output.pruned_aln}/NT/${{og}}_macse_final_mask_align_NT_pruned_complet.aln
                cp {params.dir_init_name}/${{og}}_macse/${{og}}_macse_final_unmask_align_AA.aln {output.pruned_aln}/AA/${{og}}_macse_final_unmask_align_AA_pruned_complet.aln
            fi
            done
            """


rule rerun_omm_macse:
    input:
        og_rerun = config["results_dir"] + "/cleaning/phylter_{annotation}/OG_ToReRun/NT/{og_rerun}_tokeep_phylter.fna"
    output:
        OG_align = directory(config["results_dir"] + "/cleaning/phylter_{annotation}/OG_ToReRun/MACSE_clean/{og_rerun}_macse"),
    params:
        prefix = "{og_rerun}_macse",
        macse_dir = config["results_dir"] + "/cleaning/phylter_{annotation}/OG_ToReRun/MACSE_clean",
    threads:
        1
    container:
        "library://vranwez/default/omm_macse:v12.01"
    message:
        "Re-align orthogroups that contained outliers without outliers found by phylter"
    log:
        config["results_dir"] + "/logs/cleaning/phylter_{annotation}/OG_ToReRun/MACSE_clean/{og_rerun}.log"
    shell:
        """
        mkdir -p {params.macse_dir}
        bash /OMM_MACSE/S_OMM_MACSE_V12.01.sh --out_dir {output.OG_align} --out_file_prefix {params.prefix} --in_seq_file {input.og_rerun} --genetic_code_number 1 --alignAA_soft MAFFT --java_mem 10000m &>> {log} || true
        """

def get_rerun_phylter_names(wildcards):
    checkpoint_output = checkpoints.remove_outgroup_filter_og.get(**wildcards).output[1]
    rerun_aln = expand(config["results_dir"] + "/cleaning/phylter_{annotation}/OG_ToReRun/MACSE_clean/{og_rerun}_macse", annotation=wildcards.annotation, og_rerun=glob_wildcards(os.path.join(checkpoint_output, "{og_rerun}_tokeep_phylter.fna")).og_rerun)
    return rerun_aln

rule rerun_filt_omm_macse:
    input:
        get_rerun_phylter_names,
    output:
        list_OG_rerun = config["results_dir"] + "/cleaning/phylter_{annotation}/OG_ToReRun/List_orthogroup_pass_rerun.txt",
    params:
        comparison_dir = config["results_dir"] + "/comparison/comparison_{annotation}/",
    threads:
        1
    message:
        "Write in directory for comparison only AA and NT alignments that pass filter after omm_macse"
    shell:
        """
        for d in {input}
        do
            nb_file=$(ls $d | wc -l)
            if [ "$nb_file" -ge 10 ]; then
                orthogroup=$(echo $d | awk -F "/" '{{print $NF}}' | sed 's/_macse//')
                if [ $(grep -c "^>" ${{d}}/${{orthogroup}}_macse_final_mask_align_NT.aln) -ge 3 ]; then
                    cp ${{d}}/${{orthogroup}}_macse_final_mask_align_NT.aln {params.comparison_dir}/NT/${{orthogroup}}_macse_final_mask_align_NT_pruned_complet.aln
                    cp ${{d}}/${{orthogroup}}_macse_final_unmask_align_AA.aln {params.comparison_dir}/AA/${{orthogroup}}_macse_final_unmask_align_AA_pruned_complet.aln
                    echo "$orthogroup" >> {output.list_OG_rerun};
                fi
            fi
        done
        # If checkpoint as 0 results what happen -> to test
        if [ ! -f {output.list_OG_rerun} ]; then
            touch {output.list_OG_rerun}
        fi;
        """



checkpoint comparison_remove_gene:
    input:
        pruned_aln=expand(config["results_dir"] +"/comparison/comparison_{annotation}/", annotation=ANNOTATIONS),
        list_OG_rerun = expand(config["results_dir"] + "/cleaning/phylter_{annotation}/OG_ToReRun/List_orthogroup_pass_rerun.txt", annotation=ANNOTATIONS),
    output:
        full_table=config["results_dir"] +"/comparison/Full_Table_OG_By_Query_and_Tool.txt",
        newlist=config["results_dir"] +"/comparison/NewListGene_OG_GroupByAnnotationTools.tsv",
        upset=config["results_dir"] +"/comparison/UpsetPlot_OG_AnnotationTools.png",
        profile_dir=directory(config["results_dir"] +"/profiling/AAsequences/"),
    params:
        queries=QUERIES,
        annotation=ANNOTATIONS,
        output_dir=config["results_dir"] +"/comparison/",
    singularity:
        "Container/Cleaning/cleaning_orthology.sif"
    message:
        "Comparison between annotation for orthogroups that pass all filters. The comparison is done by regarding common queries between annotation. Generate new Orthogroup list of genes based on this comparisons. The new orthogroups names are based on the order of annotation tools set in parameter"
    shell:
        """        
        pattern=$(echo "^>{params.queries}" |sed 's/ /|^>/g')
        header="#Gene_Query"
        cnt=0;
        for a in {params.annotation}
        do
        cnt=$((cnt + 1));
        header+="\\t$a"
        if [[ $cnt -eq 1 ]]; then
            grep -E "$pattern" {params.output_dir}/comparison_${{a}}/NT/*.aln | awk -F ":" '{{sub(/^>/,"", $2);sub(/_macse_final_mask_align_NT_pruned_complet.aln$/,"", $1); n=split($1,a,"/"); print $2"\\t"a[n]}}' | sort -t $'\\t' -k1,1 > {params.output_dir}/join.tmp
        elif [[ $cnt -eq 2 ]]; then
            join -t $'\\t' -a 1 -a 2 -e X -o auto <(cat {params.output_dir}/join.tmp | sort -t $'\\t' -k1,1) <(grep -E "$pattern" {params.output_dir}/comparison_${{a}}/NT/*.aln | awk -F ":" '{{sub(/^>/,"", $2);sub(/_macse_final_mask_align_NT_pruned_complet.aln$/,"", $1); n=split($1,a,"/"); print $2"\\t"a[n]}}' | sort -t $'\\t' -k1,1) > {output.full_table}
            rm {params.output_dir}/join.tmp
        else
            join -t $'\\t' -a 1 -a 2 -e X -o auto <(cat {output.full_table} | sort -t $'\\t' -k1,1) <(grep -E "$pattern" {params.output_dir}/comparison_${{a}}/NT/*.aln | awk -F ":" '{{sub(/^>/,"", $2);sub(/_macse_final_mask_align_NT_pruned_complet.aln$/,"", $1); n=split($1,a,"/"); print $2"\\t"a[n]}}' | sort -t $'\\t' -k1,1) > {params.output_dir}/join.tmp
            rm {output.full_table}
            mv {params.output_dir}/join.tmp {output.full_table}
        fi
        done
        
        sed -i -E "1s/^/$header\\n/" {output.full_table}
        python scripts/Comparison_annotation.py {output.full_table} {params.output_dir}
        
        mkdir -p {output.profile_dir}
        while read newname tools og queries nqueries status
        do
            tokeep=$(echo $og |cut -d "|" -f 1);
            tool_tokeep=$(echo $tokeep | cut -d "_" -f 1)
            og_tokeep=$(echo $tokeep | cut -d "_" -f 2)
            sed 's/[*!]/-/g' {params.output_dir}/comparison_${{tool_tokeep}}/AA/${{og_tokeep}}_macse_final_unmask_align_AA_pruned_complet.aln > {output.profile_dir}/${{newname}}_AA_final.aln
        done <{output.newlist}
        """

rule hmmbuild:
    input:
        aln = config["results_dir"] + "/profiling/AAsequences/{og_results}_AA_final.aln"
    output:
        hmm = config["results_dir"] + "/profiling/orthologue_alignments_profiles/{og_results}.hmm",
    params:
        prefix = "{og_results}_macse",
        hmm_dir = config["results_dir"] + "/profiling/orthologue_alignments_profiles/",
    threads:
        2
    log:
        config["results_dir"] + "/log/profiling/orthologue_alignments_profiles/{og_results}.out"
    container:
        "Container/Cleaning/cleaning_orthology.sif"
    message:
        "Create hmm profile by orthogroup for new orthogroups genes list"
    log:
        config["results_dir"] + "/logs/profiling/orthologue_alignments_profiles/{og_results}_hmm.log"
    shell:
        """
        hmmbuild -n HMM_{wildcards.og_results} -o {log} --cpu {threads} --amino {output.hmm} {input.aln} &>> {log}
        """

def get_hmm_names(wildcards):
    checkpoint_output = checkpoints.comparison_remove_gene.get(**wildcards).output[3]
    profile_hmm = expand(config["results_dir"] + "/profiling/orthologue_alignments_profiles/{og_results}.hmm", og_results=glob_wildcards(os.path.join(checkpoint_output, "{og_results}_AA_final.aln")).og_results)
    return profile_hmm

rule hmmpress:
    input:
        get_hmm_names,
    output:
        cat_hmm = config["results_dir"] + "/profiling/orthologue_alignments.hmm",
        hmm_press = multiext(config["results_dir"] + "/profiling/orthologue_alignments.hmm", ".h3m", ".h3i", ".h3f", ".h3p"),
    threads:
        1
    container:
        "Container/Cleaning/cleaning_orthology.sif"
    message:
        "Merge and Build hmm profile for all orthogroups together"
    log:
        config["results_dir"] + "/logs/profiling/orthologue_alignments.log"
    shell:
        """
        cat {input} > {output.cat_hmm}
        hmmpress {output.cat_hmm} &> {log}
        """

