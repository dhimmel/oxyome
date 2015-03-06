
# coding: utf-8

# In[1]:

import os
import csv
import collections
import gzip
import itertools
import warnings
import re
import sys
import multiprocessing


# In[2]:

project_dir = '/home/dhimmels/Dropbox/lung/followup/interactome/'
ii_directory = os.path.join(project_dir, 'incomplete-interactome')
ii_source_dir = os.path.join(ii_directory, 'source')

sys.path.append(ii_source_dir)
import separation


# In[3]:

def read_graph(path):
    """ """
    graph = separation.read_network(path)
    separation.remove_self_links(graph)
    return graph


# In[4]:

path_interactome = os.path.join(ii_directory, 'data-supplement', 'DataS1_interactome.tsv')
interactome = read_graph(path_interactome)
interactome_genes = set(interactome.nodes())


# In[5]:

def read_ii_tsv(path):
    read_file = open(path)
    for line in read_file:
        if line.startswith('#'):
            fieldnames = line.lstrip('#').split('\t')
            fieldnames = map(str.strip, fieldnames)
        else:
            row = line.rstrip().split('\t')
            yield collections.OrderedDict(zip(fieldnames, row))
    read_file.close()


# In[6]:

disease_genes = dict()
path = os.path.join(ii_directory, 'data-supplement', 'DataS2_disease_genes.tsv')
for disease in read_ii_tsv(path):
    omim = set(disease['OMIM_genes'].split(';') if disease.get('OMIM_genes') else [])
    gwas = set(disease['GWAS_genes'].split(';') if disease.get('GWAS_genes') else [])
    
    omim &= interactome_genes
    gwas &= interactome_genes

    name = disease['disease']
    disease_genes[name] = {'omim': omim, 'gwas': gwas, 'all': omim | gwas}


# In[7]:

# Read GO annotations
def read_go_annotations(path):
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter='\t')
    for row in reader:
        row['genes'] = set(row['genes'].split(';'))
        yield row
    read_file.close()


# In[8]:

go_path = os.path.join(project_dir, 'gene-ontology', 'human-annotations-prop.tsv')
term_to_annotation = dict()
for row in read_go_annotations(go_path):
    row['genes'] &= interactome_genes
    size = len(row['genes'])
    if size < 5 or size > 200:
        continue
    term_to_annotation[row['go_id']] = row

len(term_to_annotation)


# In[9]:

def calculate_quantities(graph, genes_A, genes_B):
    """ """
    all_genes_in_network = set(graph.nodes())
    genes_A &= all_genes_in_network
    genes_B &= all_genes_in_network

    results = dict()
    results['size_A'] = len(genes_A)
    results['size_B'] = len(genes_B)

    distance_A = separation.calc_single_set_distance(graph, genes_A)
    distance_B = separation.calc_single_set_distance(graph, genes_B)
    results['distance_A'] = distance_A
    results['distance_B'] = distance_B

    distance_AB = separation.calc_set_pair_distances(graph, genes_A, genes_B)
    results['distance_AB'] = distance_AB

    results['separation_AB'] = distance_AB - (distance_A + distance_B) / 2.0 if all(map(lambda x: x is not None, [distance_A, distance_B, distance_AB])) else None

    return results

def compute_results(args):
    """ """
    disease, go_id, association_type = args
    genes_A = disease_genes[disease][association_type]
    genes_B = term_to_annotation[go_id]['genes']
    
    results = calculate_quantities(interactome, genes_A, genes_B)

    results['disease'] = disease
    results['association_type'] = association_type
    
    results['go_id'] = go_id
    results['go_term'] = term_to_annotation[go_id]['go_term']
    results['go_domain'] = term_to_annotation[go_id]['go_domain']

    return results


# In[10]:

diseases = set(disease_genes)
go_terms = set(term_to_annotation)
association_types = 'omim', 'gwas', 'all'

def generate_parameters(min_disease_size = 5, max_disease_size = 250, min_go_size = 5, max_go_size = 250):
    for disease, go_id, association_type in itertools.product(diseases, go_terms, association_types):
        disease_size = len(disease_genes[disease][association_type])
        go_size = len(term_to_annotation[go_id]['genes'])
        if min_disease_size > disease_size:
            continue
        if max_disease_size < disease_size:
            continue
        if min_go_size > go_size:
            continue
        if max_go_size < go_size:
            continue
        yield disease, go_id, association_type

parameters = generate_parameters()


# In[11]:

#warnings.filterwarnings('ignore')

results_fields = ['disease', 'association_type', 'go_id', 'go_term', 'go_domain',
                  'size_A', 'size_B', 'distance_A', 'distance_B', 'distance_AB', 'separation_AB']
results_path = os.path.join(project_dir, 'data', 'disease-go-pairs.tsv.gz')
results_file = gzip.open(results_path, 'w')
result_writer = csv.DictWriter(results_file, delimiter = '\t', fieldnames = results_fields)
result_writer.writeheader()

pool = multiprocessing.Pool(processes = 7)
result_generator = pool.imap_unordered(compute_results, parameters, chunksize = 500)
for result in result_generator:
    result_writer.writerow(result)
results_file.close()

