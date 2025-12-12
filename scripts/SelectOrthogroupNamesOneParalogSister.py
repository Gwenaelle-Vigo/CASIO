#! /usr/bin/env python3

# Ce script utilise le package ETE pour comparer les noms des feuilles dans les arbres de gène des orthogroupes contenant un gène paraloguà partir d'une liste de noms d'orthogroupes. Il identifie les orthogroupes où le nom d'une espèce apparait deux fois dans le même arbre de gène et déterminer si ces deux gènes paralogues appartiennent à des feuilles sœurs (provenant d'un même ancêtre) ou non.

import sys
from ete3 import Tree

# Vérifier si le nombre d'arguments est correct
if len(sys.argv) != 4:
    print("Utilisation : script.py List_orthogroup_names_one_paralog.txt pathway_to_directory_Resolved_Gene_Trees > Output_file.txt")
    sys.exit(1)

# Récupérer le chemin vers le fichier contenant la liste des orthogroups possédant un gène paralogue à partir des arguments en ligne de commande
list_Orthogroups_one_paralogs = sys.argv[1]

# Récupérer le chemin vers le répertoire contenant les arbres phylogénétiques résolus
pathway_to_directory_Resolved_Gene_Trees = sys.argv[2]

# Liste des noms d'espèces à parcourir
with open(sys.argv[3], "r") as genome_file:
    for line in genome_file:
        line = line.strip()
        species=line.split(" ")

# Initialiser une liste pour stocker les noms des feuilles correspondantes
leaf_names = []

# Initialiser une liste pour stocker les orthogroupes avec des feuilles sœurs
list_OG_paralog_sisters = []

# Parcourir les orthogroups
with open(list_Orthogroups_one_paralogs, 'r') as og_file:
    for OG in og_file:
        OG = OG.strip()  # Supprimer les espaces ou retours à la ligne
        # Récupérer le chemin vers le fichier contenant l'arbre Newick au format texte
        file_tree = f"{pathway_to_directory_Resolved_Gene_Trees}/{OG}_tree.txt"

        # Lire le contenu du fichier
        try:
            with open(file_tree, 'r') as file_t:
                tree_Newick = file_t.read()
        except FileNotFoundError:
            continue  # Passer au fichier suivant

        # Charger l'arbre phylogénétique depuis le contenu du fichier
        tree = Tree(tree_Newick, format=1)  # format 1 : flexible with internal node names

        # Parcourir les noms d'especes
        for sp in species:
            # Parcourir les feuilles de l'arbre
            for leaf1 in tree:
                if leaf1.is_leaf():
                    for leaf2 in tree:
                        if leaf2.is_leaf() and leaf1 != leaf2:
                            # Vérifier si le mot spécifié est présent dans le nom des deux feuilles
                            if sp in leaf1.name and sp in leaf2.name:
                                # Ajouter les noms des feuilles correspondantes pour cet OG à la liste principale
                                leaf_names.append((OG, leaf1.name, leaf2.name))
                                # Vérifier si les feuilles sont sœurs
                                common_ancestor = tree.get_common_ancestor(leaf1.name, leaf2.name)
                                if len(common_ancestor.get_leaves()) == 2 or len(common_ancestor.get_leaves()) == 3:
                                    # Les feuilles sont sœurs
                                    if OG not in list_OG_paralog_sisters:
                                        list_OG_paralog_sisters.append(OG)
                                        print(OG)
                                break # Passer à l'OG suivant dès qu'une paire est trouvée


