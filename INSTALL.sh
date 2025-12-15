#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

cd $SCRIPT_DIR

echo -e "\n"
echo "#################### Container preparation :"
echo "Create Container directory in $SCRIPT_DIR"
mkdir -p Container/Braker/
mkdir -p Container/Busco/
mkdir -p Container/Miniprot_Scipio/
mkdir -p Container/Cleaning/
echo -e "\n"

echo "Download apptainer image in appropriate Container directory :"
echo -e "\n"
echo -e "Downloading Busco image :\n"
wget https://data.indores.fr:443/api/access/datafile/31134 -O Container/Busco/Busco.sif
echo -e "\n"
echo -e "Downloading masking Tool image :\n"
wget https://data.indores.fr:443/api/access/datafile/31133 -O Container/Braker/repeatM.sif
echo -e "\n"
echo -e "Downloading Miniprot and Scipio image :\n"
wget https://data.indores.fr:443/api/access/datafile/31132 -O Container/Miniprot_Scipio/miniprot_scipio.sif
echo -e "\n"
echo -e "Downloading Orthology and cleaing tools image :\n"
wget https://data.indores.fr:443/api/access/datafile/31131 -O Container/Cleaning/cleaning_orthology.sif

echo -e "\n"
echo "Ensure that the apptainer image download was completed successfully."
md5sum -c checkmd5sum.md5
echo -e "\n\n\n"

echo "#################### Prepare working directory for CASIO pipeline :"
echo "Create directory to set genome at format Genus_species.fasta : Data/Genome/"
mkdir -p Data/Genome/
echo -e "\n"
echo "Now, you can download and put your Queries (Genus_species.fasta and Genus_species.gff) in Data/Query/"

echo -e "\n\n\n"
echo "Create directory to set query/queries at format Genus_species.fasta and Genus_species.gff : Data/Query/"
echo -e "\n"
mkdir -p Data/Query/

echo "Now you can download and put your Queries (Genus_species.fasta and Genus_species.gff) in Data/Query/"

echo -e "\n\n\n"
echo "Create directory to set OrthoDB Protein database for Braker (https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/): Data/Protein_Database/"
mkdir -p Data/Protein_Database/
echo -e "\n"
echo "Now, if you will use Braker you have to download a orthoDB protein (https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/) and put it in  Data/Protein_Database/"

echo -e "\n\n\n"
echo "#################### Last things to do :"
echo "- Download and put Genomes, queries and protein database like explain above."
echo "- Install snakemake and apptainer."
echo "- Configure the ${SCRIPT_DIR}/configparam.yaml file."
echo "- Execute the workflow with optionaly --singularity-args to bind data and results if there are not in CASIO directory:"
echo "snakemake -s CASIO.smk --configfile configparam.yaml -c {threads} --use-singularity (--singularity-args '-B /path/where/are/data -B /path/where/are/queries/ -B /path/where/are/protein_database -B /path/where/write/results/')"



