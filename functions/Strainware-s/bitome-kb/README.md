# bitome
A module for accumulating genomic features, storing them in as a knowledge base, and representing them in binary matrix
format based on genomic position.

This repository contains all the raw data and code required to construct a feature knowledgebase and binary feature
matrix (bitome) for E. coli strain K-12 substrain MG1655 (reference sequence ID NC_000913.3)

This module is written in Python 3 and requires the following third-party packages:
```
biopython>=1.73
matplotlib>=3.0.2
numpy>=1.16.1
pandas>=0.24.1
PyWavelets>=1.0.1
scipy>=1.2.0
tqdm>=4.31.1
```

## Module Usage

### Overview

The fundamental construct of this module is the Bitome class. It has a `load_data()` method that parses all of the raw
data contained in the repository into a set of objects and links them together, creating a knowledgebase. The
`load_matrix()` method then parses this knowledgebase into a binary feature matrix with rows representing features and
columns representing genome positions.

Here is a quick example for loading a bitome object (should be run from the top-level project directory):

```python
# import the Bitome object from the core sub-module
from bitome.core import Bitome

test_bitome = Bitome()
test_bitome.load_data() # this will take a few minutes
test_bitome.load_matrix() # this will take a few minutes
```

The following sections will highlight the attributes of the `test_bitome` object that are populated after these loading
methods are run.

### Knowledgebase (load_data)

`Bitome.load_data()` parses raw data (see Data Sources section for detail on data sources and parsing) and converts into
`BitomeFeature` objects or sub-objects. 

The `BitomeFeature` class is similar to Biopython's `SeqFeature` class, storing basic information about a genomic
feature such as its name, feature type, genomic location, and sequence.

This module extends this class with a number of feature specific subclasses, such as `Gene`, `TranscriptionUnit`, 
`Protein`, etc. These classes have the capability to link to other objects. For example, the `Gene` class has a
`.protein` attribute which can be populated with a `Protein` object (for coding genes)

The `load_data()` method not only creates these feature objects but links them together to create a graph/network/
knowledgebase of genomic features. For example, consider the following method of finding the regulons that a particular
protein is part of:

```python
from bitome.core import Bitome

test_bitome = Bitome()
test_bitome.load_data()

random_protein = test_bitome.proteins[100]
regulons_with_protein = []
for tu in random_protein.gene.transcription_units:
	regulons_with_protein += tu.promoter.regulons
```

### Bitome Matrix (load_matrix)

`Bitome.load_matrix()` uses the features in the knowledgebase created within a `Bitome` object by the `load_data()` 
method and parses it into a binary matrix with different genomic features in the rows and genomic positions in columns.

The matrix is a `scipy.sparse` matrix, and is actually stored in two different formats at two different attributes.
Because it is a sparse matrix, not a dataframe, the row labels (features) are stored separately.

```python
from bitome.core import Bitome

test_bitome = Bitome()
test_bitome.load_data()
test_bitome.load_matrix()

csr_matrix = test_bitome.matrix_csr # Compressed Sparse Row matrix, efficient for row slicing
csc_matrix = test_bitome.matrix_csc # Compressed Sparse Column matrix, efficient for column slicing
row_labels = test_bitome.matrix_labels
```

Some features are reading frame or strand sensitive; in these cases, multiple rows in `row_labels` correspond to the
same feature, just in different reading frames. For example, for coding genes, the following 6 row labels correspond
to genes of type CDS (coding):

- gene_CDS_(+1)
- gene_CDS_(+2)
- gene_CDS_(+3)
- gene_CDS_(-1)
- gene_CDS_(-2)
- gene_CDS_(-3)

So the reading frame for each "sub-row" is denoted with a suffix in the row label. Similarly, for features that are 
strand-specific but not reading frame-specific (e.g. promoters):

- promoter_(+1)
- promoter_(-1)

#### Knowledgebase Data Pruning for Matrix Inclusion

The following rules are observed when parsing knowledgebase features into the bitome matrix:

- genes of type `'CDS'` that do NOT have any associated transcription units are not included
- codons, proteins, and protein properties are only annotated for gene loci that code for proteins (and have at least
one known TU, as noted above)
- of non-gene features derived from the NCBI reference sequence (see Data Sources section for raw data parsing details),
only mobile element (insertion sequence), repeat regions, and ORI are included in the bitome
- only transcription units that contain a gene are used
- operons and promoters are drawn only from the above TUs (TUs with a gene)
- additional TU- or gene-linked regulatory features (e.g. attenuators, terminators) are only included if they are linked
to a TU that is included by the above criterion
- sigmulons, i-modulons and regulons are only included if they regulate TUs/promoters included by the above criterion
- TF binding sites specifically are ALL included even if they are NOT tied to a particular promoter 

## Data Sources

### Overview

The following data sources are required for construction of the E. coli K-12 MG1655 bitome. These data are stored in
the data/ directory of this repository and were downloaded directly from the relevant online database. See specific 
sections further down for a detailed overview of how genomic features are parsed from these sources.

_NCBI_
- Provides reference genome sequence and basic genome features
- source link: https://www.ncbi.nlm.nih.gov/nuccore/NC_000913
- repository link: https://github.com/SBRG/bitome/blob/master/data/NC_000913.3.gb

_RegulonDB_
- Provides transcriptional regulatory network data
- database link: http://regulondb.ccg.unam.mx/index.jsp
- repository link: https://github.com/SBRG/bitome/tree/master/data/regulon_db10.0

_ssbio package_
- Provides calculated protein properties
- ssbio package link: https://ssbio.readthedocs.io/en/latest/index.html
- repository link: https://github.com/SBRG/bitome/tree/master/data/local-gempro and 
https://github.com/SBRG/bitome/blob/master/data/global-gempro.csv

_SBRG_
- internal data from SBRG is included in this bitome
- yTF binding sites: https://github.com/SBRG/bitome/blob/master/data/ytf_binding_sites.xlsx
- i-modulon information: https://github.com/SBRG/bitome/blob/master/data/i_modulon.csv
