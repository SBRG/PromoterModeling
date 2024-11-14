# built-in modules
from collections import defaultdict
import os
from pathlib import Path
import re
from typing import Dict, List, Tuple, Union

# third-party modules
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqUtils import GC
from Bio.SeqUtils.MeltingTemp import Tm_NN
from Bio.motifs import matrix
from dna_features_viewer import CircularGraphicRecord, GraphicFeature, GraphicRecord
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, mannwhitneyu
import seaborn as sns
from sklearn.metrics import median_absolute_error

# local packages
from bitome.util import *

# -- INPUT PATTERNS ---
PATH_OR_STRING = Union[Path, str]
FILE_OR_DF = Union[PATH_OR_STRING, pd.DataFrame]

# --- REQUIRED TABLE COLUMNS ---
ALWAYS_REQUIRED = ['locus_tag']
LOCATION = ['left', 'right', 'strand']
LOCUS_LOCATION = ALWAYS_REQUIRED + LOCATION
REQUIRED_COLUMNS = {
    'gene': LOCUS_LOCATION,
    'protein': LOCUS_LOCATION,
    'tu': LOCUS_LOCATION + ['tss'],
    'operon': LOCUS_LOCATION,
    'tss': ALWAYS_REQUIRED + ['tss', 'strand'],
    'tts': ALWAYS_REQUIRED + ['tts', 'strand'],
    'tfbs': LOCUS_LOCATION,
    'terminator': LOCUS_LOCATION,
    'attenuator': LOCUS_LOCATION,
    'rbs': LOCUS_LOCATION,
    'riboswitch': LOCUS_LOCATION,
    'custom': LOCUS_LOCATION
}

# --- PLOTTING COLORS ---
FEATURE_TYPE_TO_COLOR = {
    'gene': 'tab:blue',
    'protein': 'tab:cyan',
    'tu': 'tab:purple',
    'attenuator': 'tab:orange',
    'terminator': 'tab:red',
    'tss': 'tab:pink',
    'tfbs': 'tab:green'
}

# --- EXPRESSION COLUMN NAMES ---
EXPRESSION_COLUMNS = ['fpkm', 'tpm']

# --- DATA FILE PATH ---
DATA_DIR_PATH = Path(Path(os.path.abspath(os.path.dirname(__file__))).parent, 'data')


class Bitome:

    def __init__(self, reference_sequence: PATH_OR_STRING, name: str = None,
                 use_genbank: bool = True,
                 origin: Tuple[int, int] = None, terminus: Tuple[int, int] = None,
                 cid_boundaries: List[int] = None,
                 gene_table: FILE_OR_DF = None, protein_table: FILE_OR_DF = None,
                 tu_table: FILE_OR_DF = None, operon_table: FILE_OR_DF = None,
                 tss_table: FILE_OR_DF = None, tts_table: FILE_OR_DF = None,
                 tfbs_table: FILE_OR_DF = None, terminator_table: FILE_OR_DF = None,
                 attenuator_table: FILE_OR_DF = None, rbs_table: FILE_OR_DF = None,
                 riboswitch_table: FILE_OR_DF = None,
                 additional_tables: Dict[str, FILE_OR_DF] = None):
        """
        A Bitome object parses, calculates, represents, and visualizes genome-scale
        information about a bacterial genome that can be associated with specific base pairs.
        This construct enables rapid and repeatable sequence-based calculations and analyses to be
        performed for any set of genome-scale experimental or computational data tied to a known
        genome sequence.

        The minimal input is a reference genome sequence in either FASTA or GenBank format.
        Additional, optional tables of genomic features may be included via keyword arguments (see
        the descriptions of these arguments for details about expected formats).

        All "Bitomic" data is inputted and internally represented as "feature tables" using pandas
        DataFrames. Certain columns will be required for feature tables. Other columns will be
        optional but yield consistent functionality if included. And any custom columns may be
        included as desired. See a specific input argument for the required/options columns. At a
        minimum, the columns 'left', 'right', and 'strand' are required to locate the feature. The
        location integers are 1-INDEXED AND INCLUSIVE. The locus_tag column is also required and
        will be used as the index.

        :param PATH_OR_STRING reference_sequence: the file path for a FASTA or GenBank file
            containing, at a minimum, the full reference genome sequence for the organism.
        :param str name: the common name to use for this Bitome; this name will be used on any
            auto-generated plots and figures
        :param bool use_genbank: indicates if any genomic features contained within the provided
            GenBank reference sequence should be parsed into feature table(s). Ignored if
            reference_sequence is a FASTA file.
        :param Tuple[int, int] origin: the location of the origin of replication for this genome,
            in the form of a (left, right) tuple (1-indexed, inclusive)
        :param Tuple[int, int] terminus: the location of the replication termination region for this
            genome, in the form of a (left, right) tuple (1-indexed, inclusive)
        :param List[int] cid_boundaries: the boundary points of known chromosome interacting
            domains (CIDs) in the genome; these data will typically come from a chromosome
            conformation capture (3C or Hi-C) experiment
        :param FILE_OR_DF gene_table: a pandas DataFrame (or file path to one) containing
            information about genes located on the provided reference genome. Genes, in this usage,
            are ANY sequences that are transcribed, not only protein-coding genes.
            Required columns: left, right, strand, locus_tag
            Optional columns: name, pseudo, type, essential, y-ome, cog
        :param FILE_OR_DF protein_table: a pandas DataFrame (or file path to one) containing
            information about proteins coded in the provided reference genome.
            Required columns: left, right, strand, locus_tag
            Optional columns: name
        :param FILE_OR_DF tu_table: a pandas DataFrame (or file path to one) containing
            information about transcription units coded in the provided reference genome.
            Required columns: left, right, strand, locus_tag, tss
            Optional columns: name, sigma_factors, regulons
        :param FILE_OR_DF operon_table: a pandas DataFrame (or file path to one) containing
            information about operons coded in the provided reference genome.
            Required columns: left, right, strand, locus_tag
        :param FILE_OR_DF tss_table: a pandas DataFrame (or file path to one) containing
            information about transcription start sites (TSS) coded in the provided reference genome
            Required columns: tss, strand, locus_tag
        :param FILE_OR_DF tts_table: a pandas DataFrame (or file path to one) containing
            information about transcription terminator sites (TTS) coded in the provided genome.
            Required columns: tts, strand, locus_tag
        :param FILE_OR_DF tfbs_table: a pandas DataFrame (or file path to one) containing
            information about transcription factor binding sites located on the provided reference
            genome.
            Required columns: left, right, strand, locus_tag, tf
        :param FILE_OR_DF terminator_table: a pandas DataFrame (or file path to one) containing
            information about transcriptional terminators coded in the provided reference genome.
            Required columns: left, right, strand, locus_tag
        :param FILE_OR_DF attenuator_table: a pandas DataFrame (or file path to one) containing
            information about transcriptional attenuators coded in the provided reference genome.
            Required columns: left, right, strand, locus_tag
        :param FILE_OR_DF rbs_table: a pandas DataFrame (or file path to one) containing
            information about Shine-Dalgarno sequences (ribosome binding sites) coded in the
            provided reference genome.
            Required columns: left, right, strand, locus_tag
        :param FILE_OR_DF riboswitch_table: a pandas DataFrame (or file path to one) containing
            information about riboswitches coded in the provided reference genome
            Required columns: left, right, strand, locus_tag
        :param FILE_OR_DF additional_tables: any custom data types not represented by a particular
            input option
        """

        reference_sequence_path = Path(reference_sequence)
        reference_sequence_file_type = reference_sequence_path.suffix.lower()[1:]
        if reference_sequence_file_type in ['fasta', 'gb']:
            # this Biopython function errors if the FASTA file doesn't have exactly one record
            reference_record = SeqIO.read(reference_sequence_path, reference_sequence_file_type)
        else:
            raise ValueError(f'{reference_sequence_file_type} file extension not recognized. Please'
                             f' provide either a FASTA file (.fasta) or a GenBank file (.gb)')

        self.reference_record = reference_record
        self.reference_id = reference_record.id
        self.description = reference_record.description
        self.sequence = reference_record.seq.upper()
        self.seq_length = int(len(self.sequence))
        # set this variable as a backup to ensure we always have a copy in case we do mutations
        self._reference_sequence = self.sequence

        if name is not None:
            self.name = name
        else:
            self.name = self.description

        self.origin = origin
        self.terminus = terminus
        if cid_boundaries is not None:
            self.cid_boundaries = sorted(cid_boundaries)
        else:
            self.cid_boundaries = None

        # gene and protein tables are prepared differently depending on whether we got a GenBank
        if use_genbank and reference_sequence_file_type == 'gb':
            self.gene_table, self.protein_table, self.misc_feature_table = \
                genbank_to_feature_tables(self.reference_record)
            # if the gene_table was also provided, treat it as a supplementary table, but use
            # GenBank as ground truth (i.e. a left join)
            gene_table_argument = _table_input_to_dataframe(gene_table, 'gene', fast_track=True)
            if not gene_table_argument.empty:
                self.gene_table = self.gene_table.merge(gene_table_argument, how='left',
                                                        left_index=True, right_index=True)
            protein_table_argument = _table_input_to_dataframe(protein_table, 'protein',
                                                               fast_track=True)
            if not protein_table_argument.empty:
                self.protein_table = self.protein_table.merge(protein_table_argument, how='left',
                                                              left_index=True, right_index=True)
        else:
            self.gene_table = _table_input_to_dataframe(gene_table, 'gene')
            self.protein_table = _table_input_to_dataframe(protein_table, 'protein')
            # dummy line that just makes an empty dataframe
            self.misc_feature_table = _table_input_to_dataframe(None, 'custom')

        self.tu_table = _table_input_to_dataframe(tu_table, 'tu')

        # add in a handy lookup of TU IDs to gene IDs if we have a gene_names column in TU table
        self._tu_to_genes = {}
        if 'gene_names' in self.tu_table:
            for tu_row in self.tu_table.itertuples():
                self._tu_to_genes[tu_row.Index] = list(
                    set(tu_row.gene_names.split(';')).intersection(
                        set(list(self.gene_table.index))
                    )
                )
        else:
            for tu_row in self.tu_table.itertuples():
                g_in_range = self.gene_table[
                    (self.gene_table['strand'] == tu_row.strand) &
                    (self.gene_table['left'] >= tu_row.left) &
                    (self.gene_table['right'] <= tu_row.right)
                    ]
                self._tu_to_genes[tu_row.Index] = list(g_in_range.index)

        # also prepare the reverse lookup, gene IDs to TU IDs
        self._gene_to_tus = defaultdict(list)
        for tu_id, tu_gene_ids in self._tu_to_genes.items():
            for tu_gene_id in tu_gene_ids:
                if tu_gene_id in self._gene_to_tus:
                    self._gene_to_tus[tu_gene_id].append(tu_id)
                else:
                    self._gene_to_tus[tu_gene_id] = [tu_id]

        self.operon_table = _table_input_to_dataframe(operon_table, 'operon')
        self.tss_table = _table_input_to_dataframe(tss_table, 'tss')
        self.tts_table = _table_input_to_dataframe(tts_table, 'tts')
        self.tfbs_table = _table_input_to_dataframe(tfbs_table, 'tfbs')
        self.terminator_table = _table_input_to_dataframe(terminator_table, 'terminator')
        self.attenuator_table = _table_input_to_dataframe(attenuator_table, 'attenuator')
        self.rbs_table = _table_input_to_dataframe(rbs_table, 'rbs')
        self.riboswitch_table = _table_input_to_dataframe(riboswitch_table, 'riboswitch')

        if additional_tables is not None:
            self.additional_tables = {k: _table_input_to_dataframe(v, 'custom')
                                      for k, v in additional_tables.items()}
        else:
            self.additional_tables = {}

        # internally, we need a nice way to access all of the tables that we have; use a dict
        self._all_tables = {
            'gene': self.gene_table,
            'protein': self.protein_table,
            'tu': self.tu_table,
            'operon': self.operon_table, 
            'tss': self.tss_table,
            'tts': self.tts_table,
            'tfbs': self.tfbs_table,
            'terminator': self.terminator_table,
            'attenuator': self.attenuator_table,
            'rbs': self.rbs_table,
            'riboswitch': self.riboswitch_table,
            'misc': self.misc_feature_table
        }
        self._all_tables = {**self._all_tables, **self.additional_tables}

        # determine if we have any expression data present in our gene table
        self._exp_col = None
        for exp_col_known in EXPRESSION_COLUMNS:
            if exp_col_known in self.gene_table.columns:
                self._exp_col = exp_col_known
                self._exp_col_log = f'log_{self._exp_col}'
                break

        # automatically add a log-transformed expression column if we have expression data
        if self._exp_col is not None:
            self.gene_table[self._exp_col_log] = np.log2(self.gene_table[self._exp_col] + 1)

        # add an empty attribute in case we compute sigma motifs
        self.sigma_motifs = None
        
        # initiate the 5mer shape lookup
        self.fivemer_shape_lookup = pd.read_csv(Path(DATA_DIR_PATH, '5mer_shape_lookup.csv'), index_col=0)

        # initiate empty hidden attributes for storing TSS shape information for TSS finding
        self._tss_search_shape_stats = None
        self._tss_search_shape_dict = None
        self._tss_search_n_upstream = None
        self._tss_reference_shape = None
        self._tss_reference_n_up = None
        self._tss_reference_n_down = None

    # --- BASIC SEQUENCE UTILITY METHODS ---

    def get_sequence(self, left: int, right: int, strand: Union[int, None]) -> Seq:
        """
        Given a left, right, and strand positions, return the sequence for the requested region
        
        NOTE: this handles a case where left < 1 or right > len(sequence); interprets as
        "wraparound" case and handles specially

        :param int left: the left end of the sequence to return
        :param int right: the right end of the sequence to return
        :param Union[int, None] strand: the strand end of the sequence to return; may be None
        """
        if left < 1:
            far_right_wraparound = FeatureLocation(
                self.seq_length - 1 + int(left), self.seq_length, int(strand)
            ).extract(self.sequence)
            far_left_continued = FeatureLocation(
                0, int(right), int(strand)
            ).extract(self.sequence)
            # have to stitch the sequences together differently depending on the strand
            if strand == 1:
                return far_right_wraparound + far_left_continued
            else:
                return far_left_continued + far_right_wraparound
        # BioPython's handy FeatureLocation object ONLY accepts Python integers; cast to these
        # make sure to handle the fact that Bitome uses 1-indexing
        else:
            return FeatureLocation(int(left)-1, int(right), int(strand)).extract(self.sequence)

    def get_feature_sequence(self, feature_row: pd.Series) -> Seq:
        """
        Given a feature row pulled from a feature table, return the sequence (in coding order,
        where possible) of the feature

        :param pd.Series feature_row: a row from a feature table for which to return the sequence
        :return Seq feature_sequence: the sequence of the requested feature, in coding order
        """
        left, right, strand = feature_row[['left', 'right', 'strand']]
        return self.get_sequence(left, right, strand)

    def gc_content(self, feature_row: pd.Series) -> float:
        """
        Given a feature row from a feature table, return the GC content of the feature's sequence

        :param pd.Series feature_row: the data table row for the feature
        :return float gc_content: the GC content for the feature as a ratio, NOT a percentage
        """
        feature_sequence = self.get_feature_sequence(feature_row)
        return GC(feature_sequence) / 100

    def mutate_sequence(self, mutations: Dict[int, str]):
        """
        Temporarily changes the raw sequence associated with this Bitome for the purposes of
        investigating the effect of mutation(s). Use restore_original_sequence() to undo any effects
        of calling this function. BE CAREFUL with this function, as it will affect all other
        functions that utilize the sequence UNTIL restore_original_sequence is called.

        If the original sequence is needed, it can be found at Bitome._reference_sequence

        :param Dict[int, str] mutations: a set of mutations to apply to the sequence, in the form
             of position (1-indexed) to new base key/value pairs. NOTE: these mutations are ALWAYS
             considered in the forward strand direction!
        """

        mutated_sequence = list(self.sequence)
        for pos, new_base in mutations.items():
            mutated_sequence[pos-1] = new_base
        self.sequence = Seq(''.join(mutated_sequence))

    def restore_original_sequence(self):
        """
        Undo the effects of calling mutate_sequence and reset the Bitome.sequence attribute to its
        original, reference state
        """
        self.sequence = self._reference_sequence

    # --- GENOME ORGANIZATION ---

    def genome_organization_table(self, n_genome_bins: int = 40, primary_tu: bool = False,
                                  tm_range: Tuple[int, int] = (-12, 5),
                                  tm_use_box: bool = True, n_box_10_ext: int = 3) -> pd.DataFrame:
        """
        A mega-function which creates a DataFrame of genome organization features for this Bitome's
        gene_table. Calls genome_organization_for_gene iteratively over the gene_table (see that
        function's documentation for details about features present).

        :param int n_genome_bins: the number of discrete bins in which to split the genome when
            approximating chromosome interacting domains
        :param bool primary_tu: indicates if only the primary tu (from the primary_tu column of the
            gene_table) should be considered for TU-related feature calculations
        :param Tuple[int, int] tm_range: indicates the range (relative to TSS at 0) that should be
            used for determining the promoter melting temperature
        :param bool tm_use_box: indicates if the upstream edge of the -10 box should be used instead
            of the default range for melting temperature range. Defaults to True.
        :param int n_box_10_ext: the number of bp upstream of the -10 box to consider as the
            extended -10 box (specific to the box_10_ext_gc feature)
        :return genome_organization_df pd.DataFrame: a DataFrame of the same index as the
            gene_table, containing a series of genome organization-related features in columns
        """

        genome_organization_rows = []
        for _, gene_row in self.gene_table.iterrows():
            genome_organization_rows.append(
                self.genome_organization_for_gene(gene_row, n_genome_bins=n_genome_bins,
                                                  primary_tu=primary_tu, tm_range=tm_range,
                                                  tm_use_box=tm_use_box, n_box_10_ext=n_box_10_ext)
            )
        genome_organization_df = pd.DataFrame(genome_organization_rows, index=self.gene_table.index)
        # drop any column that has absolutely no data
        genome_organization_df = genome_organization_df.dropna(how='all', axis=1)
        return genome_organization_df

    def genome_organization_for_gene(self, gene_row: pd.Series, n_genome_bins: int = 40,
                                     primary_tu: bool = False, tm_range: Tuple[int, int] = (-12, 5),
                                     tm_use_box: bool = True, n_box_10_ext: int = 3) -> dict:
        """
        This mega-function computes a suite of genome organization-related features for a gene row
        from this Bitome's gene table. The features are as follows:

        Basic Features (always computed)
        - Genome location (genome_loc): the location (as a discrete bin) of the gene within the
            genome; actually represented as genome_loc_sin and genome_loc_cos to capture genome's
            circularity. Control the number of bins with the n_genome_bins parameter

        ------

        Origin/Terminus-Dependent Features
        - replication region (rep_region): whether the gene falls on the leading or lagging strand
            w.r.t. genome replication. Only computed if origin and terminus are known.
        - distance to origin (origin_dist): the distance from the gene to the replication origin;
            again, only computed if origin is known

        ------

        TU-dependent Features
            All TU features come in 3 flavors:
            - None if the gene has no known TU
            - An average of the values for each of the gene's known TUs
            - The value for the gene's primary_tu ONLY (must set primary_tu to True and have a
                primary_tu column in the gene_table)
        - TU length (tu_len): the length of the TU(s) for the gene
        - TU GC content (tu_gc): the GC content of the gene's TU(s)
        - TU methylation sites (tu_gatc): indicates the presence of 1 or more GATC Dam methylation
            motifs in the core promoter region (-50 to +10 from TSS)
        - TU melting temp: the melting temperature (using nearest-neighbor from Biopython) of the
            region from the upstream end of the -10 box (if known; -12 otherwise) to +5
        - Order in TU (tu_order): the order of transcription of the gene within its TU(s)
        - TSS distance (tss_dist): the distance from the gene start to its TU's (or TUs') TSS
        - TSS base (tss_base): the literal DNA base of the TSS of the gene's TU(s); for genes with
            multiple TUs, the most common base is used
        - 5' UTR length (utr_len): the length of the 5'UTR for the gene's TU(s)
        - 5' UTR GC content (utr_gc) : the GC content of the 5'UTR for the gene's TU(s)
        - (if cid_boundaries is provided) also determine the actual CID that a gene belongs to, and
            also the distance from the gene's TSS to the nearest CID boundary

        ------

        Pribnow/-35 Box Features (illustrated for -10 box); requires prior annotation of -10 and -35
        box locations for each TU. Same primary_tu behavior as basic TU-dependent features
        - Box sequence (box_10_seq): the sequence of the -10 box (categorical)
        - Box distance to TSS (box_10_tss_dist): the distance from the -10 box to the TSS
        - Extended box GC content (box_10_ext_gc): the GC content of the nucleotides upstream of the
            -10 box (n of nucleotides controllable via n_box_10_ext parameter)

        :param pd.Series gene_row: the gene row (from the gene_table) for which to compute features
        :param int n_genome_bins: the number of discrete bins in which to split the genome when
            approximating chromosome interacting domains
        :param bool primary_tu: indicates if only the primary tu (from the primary_tu column of the
            gene_table) should be considered for TU-related feature calculations
        :param Tuple[int, int] tm_range: indicates the range (relative to TSS at 0) that should be
            used for determining the promoter melting temperature
        :param bool tm_use_box: indicates if the upstream edge of the -10 box should be used instead
            of the default range for melting temperature range. Defaults to True.
        :param int n_box_10_ext: the number of bp upstream of the -10 box to consider as the
            extended -10 box (specific to the box_10_ext_gc feature)
        :return genome_organization_dict dict: a dictionary of the genome organization features
            for this gene
        """

        go_row = {}

        # --- Genome location ---
        if gene_row.strand == 1:
            gene_start = gene_row.left
        else:
            gene_start = gene_row.right
        genome_loc = np.floor(gene_start / len(self.sequence) * n_genome_bins)
        go_row['genome_loc_sin'] = np.sin(2 * np.pi * genome_loc / n_genome_bins)
        go_row['genome_loc_cos'] = np.cos(2 * np.pi * genome_loc / n_genome_bins)

        # --- Replication region and distance to origin ---
        if self.origin is not None:
            go_row['rep_region'] = self.replication_region(gene_row)
            go_row['origin_dist'] = self.distance_to_origin(gene_row)

        # --- TU-dependent features ---

        # only proceed if we have at least 1 TU for the gene; also only use TUs with TSS
        if primary_tu:
            if not pd.isna(gene_row.primary_tu):
                tus_for_gene = [gene_row.primary_tu]
            else:
                tus_for_gene = []
        else:
            tus_for_gene = self._gene_to_tus[gene_row.name]
        tus_for_gene_df = self.tu_table.loc[tus_for_gene].dropna(subset=['tss'])

        if not tus_for_gene_df.empty:
            tu_features = ['tu_len', 'tu_gc', 'tu_order', 'tss_dist', 'tss_base', 'utr_len',
                           'utr_gc', 'box_10_seq', 'box_10_tss_dist', 'box_10_ext_gc',
                           'box_35_seq', 'box_35_tss_dist', 'spacer_len', 'spacer_gc',
                           'cid', 'cid_bound_dist', 'tu_tm', 'tu_gatc']
            tu_feat_dict_lst = {tu_f: [] for tu_f in tu_features}
            for g_tu_id, g_tu_row in tus_for_gene_df.iterrows():

                # we can get the TU length, GC content, TSS base, and TSS dist right from the TU
                # define some useful strand-specific values
                start_col = 'left' if g_tu_row.strand == 1 else 'right'
                tu_feat_dict_lst['tu_len'].append(g_tu_row.right - g_tu_row.left + 1)
                tu_feat_dict_lst['tu_gc'].append(self.gc_content(g_tu_row))
                g_tu_tss = g_tu_row.tss
                tss_base = str(self.get_sequence(g_tu_tss, g_tu_tss, g_tu_row.strand))
                tu_feat_dict_lst['tss_base'].append(tss_base)
                tu_feat_dict_lst['tss_dist'].append(abs(g_tu_tss - gene_row[start_col]))

                # get the melting temp of the region around -10 box and TSS
                if g_tu_row.strand == 1:
                    if tm_use_box and pd.notna(g_tu_row.box_10_left):
                        melt_left, melt_right = g_tu_row.box_10_left, g_tu_tss + tm_range[1]
                    else:
                        melt_left, melt_right = g_tu_tss + tm_range[0], g_tu_tss + tm_range[1]
                else:
                    if tm_use_box and pd.notna(g_tu_row.box_10_right):
                        melt_left, melt_right = g_tu_tss - tm_range[1], g_tu_row.box_10_right
                    else:
                        melt_left, melt_right = g_tu_tss - tm_range[1], g_tu_tss - tm_range[0]
                # handle crazy fim case where -10 is downstream of TSS
                if melt_left < melt_right:
                    tu_tm = Tm_NN(self.get_sequence(melt_left, melt_right, g_tu_row.strand))
                    tu_feat_dict_lst['tu_tm'].append(tu_tm)

                # check for GATC (or CTAG to handle opposite strand) methylation sites in promoter
                if g_tu_row.strand == 1:
                    core_prom_l, core_prom_r = g_tu_tss - 50, g_tu_tss + 10
                else:
                    core_prom_l, core_prom_r = g_tu_tss - 10, g_tu_tss + 50
                core_prom_seq = self.get_sequence(core_prom_l, core_prom_r, g_tu_row.strand)
                meth_sites = re.findall(r'(GATC|CTAG)', str(core_prom_seq))
                tu_feat_dict_lst['tu_gatc'].append(len(meth_sites) > 0)

                # for TU order, 5' UTR calculations, we need the gene order in the TU
                tu_genes = self._tu_to_genes[g_tu_id]
                tu_gene_df = self.gene_table.loc[tu_genes]
                sorted_genes, sorted_starts = zip(*sorted(
                    zip(tu_genes, tu_gene_df[start_col]),
                    key=lambda tup: tup[1], reverse=(g_tu_row.strand == -1)
                ))
                gene_order_dict = {g_id: i + 1 for i, g_id in enumerate(sorted_genes)}
                tu_feat_dict_lst['tu_order'].append(gene_order_dict[gene_row.name])
                if g_tu_row.strand == 1:
                    utr_l, utr_r = g_tu_tss, sorted_starts[0] - 1
                else:
                    utr_l, utr_r = sorted_starts[0] + 1, g_tu_tss
                # handle leaderless transcript case
                if utr_l > utr_r:
                    utr_seq = Seq('')
                else:
                    utr_seq = self.get_sequence(utr_l, utr_r, g_tu_row.strand)
                tu_feat_dict_lst['utr_len'].append(len(utr_seq))
                tu_feat_dict_lst['utr_gc'].append(GC(utr_seq) / 100)

                # for -10/-35box features, we need to have the annotated location of the box in TU
                if not pd.isna(g_tu_row.box_10_left):
                    b10_l, b10_r = g_tu_row.box_10_left, g_tu_row.box_10_right
                    box_10_seq = str(self.get_sequence(b10_l, b10_r, g_tu_row.strand))
                    tu_feat_dict_lst['box_10_seq'].append(box_10_seq)
                    b10_dist = np.min(np.abs(np.array([b10_l, b10_r]) - g_tu_tss))
                    tu_feat_dict_lst['box_10_tss_dist'].append(b10_dist)
                    if g_tu_row.strand == 1:
                        b10_ext_seq = self.get_sequence(b10_l - n_box_10_ext, b10_l - 1, 1)
                    else:
                        b10_ext_seq = self.get_sequence(b10_r + 1, b10_r + n_box_10_ext, -1)
                    tu_feat_dict_lst['box_10_ext_gc'].append(GC(b10_ext_seq) / 100)
                # handle case with no -35 box
                if 'box_35_left' in g_tu_row and not pd.isna(g_tu_row.box_35_left):
                    b35_l, b35_r = g_tu_row.box_35_left, g_tu_row.box_35_right
                    box_35_seq = str(self.get_sequence(b35_l, b35_r, g_tu_row.strand))
                    tu_feat_dict_lst['box_35_seq'].append(box_35_seq)
                    b35_dist = np.min(np.abs(np.array([b35_l, b35_r]) - g_tu_tss))
                    tu_feat_dict_lst['box_35_tss_dist'].append(b35_dist)

                    # we can do spacer calculations if we have both the -10/-35 box locations
                    if not pd.isna(g_tu_row.box_10_left):
                        if g_tu_row.strand == 1:
                            spacer_l, spacer_r = g_tu_row.box_35_right + 1, g_tu_row.box_10_left - 1
                        else:
                            spacer_l, spacer_r = g_tu_row.box_10_right + 1, g_tu_row.box_35_left - 1
                        # crazy case with fim TU where box10/35 are downstream of TSS, need invert?
                        if spacer_l < spacer_r:
                            spacer_seq = self.get_sequence(spacer_l, spacer_r, g_tu_row.strand)
                            tu_feat_dict_lst['spacer_len'].append(len(spacer_seq))
                            tu_feat_dict_lst['spacer_gc'].append(GC(spacer_seq) / 100)

                # if we have CID boundaries, with this bitome, we can get the CID of this gene
                # and also its distance to the nearest boundary (based on the TSS of TU)
                if self.cid_boundaries is not None:
                    cid_bounds = self.cid_boundaries
                    # need to be careful with wraparound when finding CID the TSS is in
                    cid_dists = []
                    for cid_bound in cid_bounds:
                        dist = genome_point_to_point(g_tu_tss, cid_bound, len(self.sequence))
                        cid_dists.append(dist)
                    min_idx = int(np.argmin(cid_dists))
                    # figure out which side of the closest CID boundary we're on to define CID
                    # we're in the edge case if we're higher than the largest index OR lower than
                    # the smallest index; define the wraparound CID to be cid_highest
                    high_idx = len(cid_dists) - 1
                    if min_idx == 0:
                        # the gene is to the right of the 0th boundary
                        if cid_bounds[0] <= g_tu_tss < cid_bounds[high_idx]:
                            tu_feat_dict_lst['cid'].append('cid_0')
                        # the gene is in the wraparound CID
                        else:
                            tu_feat_dict_lst['cid'].append(f'cid_{high_idx}')
                    elif min_idx == high_idx:
                        # the gene is not in the wraparound, in the one before
                        if cid_bounds[high_idx] <= g_tu_tss < cid_bounds[high_idx]:
                            tu_feat_dict_lst['cid'].append(f'cid_{high_idx - 1}')
                        # the gene is in the wraparound CID
                        else:
                            tu_feat_dict_lst['cid'].append(f'cid_{high_idx}')
                    # normal case don't have to worry about edges
                    else:
                        if g_tu_tss >= cid_bounds[min_idx]:
                            tu_feat_dict_lst['cid'].append(f'cid_{min_idx}')
                        else:
                            tu_feat_dict_lst['cid'].append(f'cid_{min_idx - 1}')
                    tu_feat_dict_lst['cid_bound_dist'].append(cid_dists[min_idx])

            tu_feat_dict_avg = {}
            for tu_f, f_values in tu_feat_dict_lst.items():
                if f_values:
                    # need the mode for collapsing categorical features
                    # TODO: handle ties, or handle this better generally
                    if tu_f in ['tss_base', 'box_10_seq', 'box_35_seq', 'cid', 'tu_gatc']:
                        modes, _ = mode(f_values)
                        tu_feat_dict_avg[tu_f] = modes[0]
                    else:
                        tu_feat_dict_avg[tu_f] = np.mean(f_values)
            go_row = {**go_row, **tu_feat_dict_avg}

        return go_row

    def inter_feature_distance(self, feature1_row: pd.Series, feature2_row: pd.Series,
                               feature1_by: str = 'start', feature2_by: str = 'start',
                               relative_to_origin: bool = False) -> float:
        """
        Given the feature table rows for two features, determine the distance between them on the
        genome. The reference point for each feature can be defined. The circularity of the genome
        will be taken into account. If the location of the replication origin is provided, this may
        be used to determine the relative distances of the features from the origin.
        NOTE: if no strand information is known for the feature, the midpoint will always be used.

        :param pd.Series feature1_row: 1st feature, pulled directly from a feature table; must have
            required columns left, right, strand to be located properly
        :param pd.Series feature2_row: 2nd feature, pulled directly from a feature table; must have
            required columns left, right, strand to be located properly
        :param str feature1_by: choose the reference point for feature 1; 'start' (default), 'end',
            or 'midpoint'
        :param str feature2_by: choose the reference point for feature 1; 'start' (default), 'end',
            or 'midpoint'
        :param bool relative_to_origin: indicates if distances will be computed relative to the
            origin (i.e., the difference between how far each feature is from the origin of
            replication); origin location must be known
        :return float distance: the distance between the 2 features
        """

        f1_left, f1_right, f1_strand = feature1_row[['left', 'right', 'strand']]
        f2_left, f2_right, f2_strand = feature2_row[['left', 'right', 'strand']]

        f1_point = location_to_point(f1_left, f1_right, f1_strand, feature1_by)
        f2_point = location_to_point(f2_left, f2_right, f2_strand, feature2_by)

        genome_length = len(self.sequence)

        if relative_to_origin:
            if self.origin is None:
                raise ValueError('Cannot compute distances relative to origin without origin.')
            else:
                origin_midpoint = (self.origin[1] + self.origin[0]) / 2
                f1_to_origin = genome_point_to_point(f1_point, origin_midpoint, genome_length)
                f2_to_origin = genome_point_to_point(f2_point, origin_midpoint, genome_length)
                distance = abs(f2_to_origin - f1_to_origin)
        else:
            distance = genome_point_to_point(f1_point, f2_point, genome_length)

        return distance

    def replication_region(self, feature_row: pd.Series) -> str:
        """
        Determine whether a given feature falls on the leading or lagging strand w.r.t. DNA
        replication; a third option, 'terminus', may be returned if the feature lies within the
        innermost termination sequences, making it impossible to guarantee in which direction the
        replication fork will pass the feature.
        Note: assumes that the reference sequence numbering is in typical, 5'->3' direction; in
        other words, coding direction is in the direction of increasing genome position when on +1
        strand

        :param pd.Series feature_row: a feature row pulled from a feature table
        :return str region: a string describing the location of the feature w.r.t. DNA replication;
            possible values are 'leading', 'lagging', and 'terminus'
        """

        if self.origin is None or self.terminus is None:
            raise ValueError('Origin and terminus must be known to determine replication region.')

        # make a temporary mapping such that the origin's left is located at genome position 1
        # subtracting 1 ensures that the origin's left position will end up at position 1 if this
        # offset is subtracted from the origin's left position (the goal)
        offset = self.origin[0] - 1

        def offset_position(position):
            if position > offset:
                new_position = position - offset
            else:
                new_position = len(self.sequence) - (offset - position)
            return new_position

        origin_right_offset = offset_position(self.origin[1])
        term_left_offset = offset_position(self.terminus[0])
        term_right_offset = offset_position(self.terminus[1])

        f_left, f_right, f_strand = feature_row[['left', 'right', 'strand']]
        if f_strand not in [1, -1]:
            raise ValueError(f'Strand unknown for feature from {f_left, f_right}; replication '
                             f'region cannot be determined.')
        f_start = location_to_point(f_left, f_right, f_strand, reference='start')
        f_start_offset = offset_position(f_start)

        # ensure 1-indexing is preserved when using Python's range function for convenience
        if f_start_offset in range(origin_right_offset, term_left_offset + 1):
            if f_strand == 1:
                region = 'leading'
            else:
                region = 'lagging'
        elif f_start_offset in range(term_left_offset, term_right_offset + 1):
            region = 'terminus'
        else:
            if f_strand == 1:
                region = 'lagging'
            else:
                region = 'leading'

        return region

    def distance_to_origin(self, feature_row: pd.Series) -> float:
        """
        Given a feature row, return the distance from that feature to the origin of replication
        (will return the closest possible distance, ignoring the feature's reading direction)

        :param pd.Series feature_row: the feature row for which to compute origin distance
        :return float origin_dist: the distance from the feature to the origin
        """
        l, r = feature_row.left, feature_row.right
        seq_max = len(self.sequence)
        return np.min([
            genome_point_to_point(self.origin[0], l, seq_max),
            genome_point_to_point(self.origin[0], r, seq_max),
            genome_point_to_point(self.origin[1], l, seq_max),
            genome_point_to_point(self.origin[1], r, seq_max)
        ])

    # --- FEATURE EXTRACTION AND VISUALIZATION ---

    def features_in_range(self, left: int, right: int) -> pd.DataFrame:
        """
        Given a range on the reference genome, return a single DataFrame containing all features
        found in the provided range

        :param int left: the left end of the range for which to extract features
        :param int right: the right end of the range for which to extract features
        """

        master_range_table = pd.DataFrame()

        for label, table in self._all_tables.items():
            label_range_table = table[
                ((left <= table['left']) & (table['left'] <= right)) |
                ((left <= table['right']) & (table['right'] <= right)) |
                ((left >= table['left']) & (right <= table['right']))
            ].copy()
            # TODO: handle GenBank misc features more cleanly
            if not label_range_table.empty:
                if label != 'misc':
                    label_range_table.loc[:, 'type'] = label
                master_range_table = master_range_table.append(label_range_table)

        return master_range_table

    def view_region(self, left: int, right: int) -> plt.Axes:
        """
        Given a genome range defined by a left and right position, display a genome browser-style
        figure highlighting all known features present in the region

        :param int left: the left end of the range for which to extract features
        :param int right: the right end of the range for which to extract features
        :return plt.Axes ax: the Axes containing the genome region plot
        """

        range_feature_table = self.features_in_range(left, right)

        graphic_features = []
        legend_types = set()

        for feature_row in range_feature_table.itertuples(index=False):

            color = FEATURE_TYPE_TO_COLOR.get(feature_row.type, 'gray')
            if isinstance(feature_row.name, str):
                if feature_row.type in ['protein']:
                    label = None
                else:
                    label = feature_row.name
            else:
                if feature_row.type in ['repeat_region', 'attenuator', 'terminator']:
                    label = None
                else:
                    label = feature_row.type

            graphic_feature = GraphicFeature(start=feature_row.left, end=feature_row.right,
                                             strand=feature_row.strand, label=label,
                                             color=color, box_color=color)
            graphic_features.append(graphic_feature)
            legend_types.add(feature_row.type)

        seq_length = (right-left+1)
        graphic_record = GraphicRecord(sequence=self.get_sequence(left, right, 1),
                                       sequence_length=seq_length, features=graphic_features,
                                       first_index=left, labels_spacing=2)
        ax, _ = graphic_record.plot(plot_sequence=(seq_length < 100))

        legend_handles = [Patch(color=FEATURE_TYPE_TO_COLOR.get(legend_type, 'gray'))
                          for legend_type in legend_types]
        ax.legend(legend_handles, list(legend_types), bbox_to_anchor=(1, 0.5))

        return ax

    def view_genome_scale_data(self, feature_types: List[str]) -> plt.Axes:
        """
        Plot all instances of a certain set of feature types in a circular genome plot

        :param List[str] feature_types: the feature type(s) to plot along the entire circular genome
        :return Axes ax: the Axes containing the circular genome plot
        """

        graphic_features = []
        for feature_type in feature_types:
            feature_table = self._all_tables[feature_type]
            for feature_row in feature_table.itertuples(index=False):
                graphic_feature = GraphicFeature(start=feature_row.left, end=feature_row.right,
                                                 strand=feature_row.strand,
                                                 color=FEATURE_TYPE_TO_COLOR[feature_type],
                                                 linecolor=FEATURE_TYPE_TO_COLOR[feature_type])
                graphic_features.append(graphic_feature)

        graphic_record = CircularGraphicRecord(sequence_length=len(self.sequence),
                                               features=graphic_features)
        ax, _ = graphic_record.plot()
        legend_handles = [Patch(color=FEATURE_TYPE_TO_COLOR.get(legend_type, 'gray'))
                          for legend_type in feature_types]
        ax.legend(legend_handles, list(feature_types), bbox_to_anchor=(1, 0.5))

        return ax

    # --- MOTIF SCORING AND PROMOTERS ---

    def featurize_promoter(self, tss: int, strand: int, minus_10_motif: Dict[str, List[float]],
                           minus_35_motif: Dict[str, List[float]],
                           minus_10_search: Tuple[int, int] = (-20, 0),
                           minus_35_search: Tuple[int, int] = (-45, -25),
                           usr: Tuple[int, int] = (-65, -45), dsr: Tuple[int, int] = (9, 20)
                           ) -> Dict[str, float]:
        """
        Given a transcription start site (TSS) (and its strand), looks in the core promoter region
        to find matches to provided -10 and -35 element motifs (sigma factor binding sites). The
        motifs should be provided as pre-made PSSMs (see argument for details of formatting). A
        linear search within the search space (customizable) is performed, and statistics for the
        best scoring sequence (if above a certain match threshold) is returned.

        :param int tss: the location of the TSS upstream of which to look for sigma factor motifs
        :param int strand: the strand of the TSS upstream of which to look for sigma factor motifs
        :param Dict[str, List[float]] minus_10_motif: the motif for the -10 element (Pribnow box)
            binding site of the sigma factor of interest. A PSSM in dictionary form is expected. The
            dictionary must have keys 'A', 'T', 'C', 'G', and each value must be a list of log-odds
            values for the corresponding nucleotide key at the appropriate index in the motif (in
            the 5' -> 3' direction). See Biopython's motifs.create function for a way to make one.
        :param Dict[str, List[float]] minus_35_motif: the motif for the -35 element binding site of
            the sigma factor of interest. A PSSM in dictionary form is expected. The dictionary must
            have keys 'A', 'T', 'C', 'G', and each value must be a list of log-odds values for the
            corresponding nucleotide key at the appropriate index in the motif (in the 5' -> 3'
            direction). See Biopython's motifs.create function for a way to make one.
        :param Tuple[int, int] minus_10_search: the sequence region in which to linearly search
            for the provided -10 motif upstream of the given TSS; defaults to -20 to 0 (expected to
            be w.r.t. TSS at 0)
        :param Tuple[int, int] minus_35_search: the sequence region in which to linearly search
            for the provided -35 motif upstream of the given TSSl defaults to -45 to -25 (expected
            to be w.r.t. TSS at 0)
        :param Tuple[int, int] usr: the sequence region to be considered the upstream sequence
            region (USR) for calculation of AT content, upstream of the given TSS; defaults to -65
            to -45 (expected to be w.r.t. TSS at 0)
        :param Tuple[int, int] dsr: the sequence region to be considered the downstream sequence
            region (DSR) for calculation of puring content; defaults to +9 to +20 (expected to be
            w.r.t. TSS at 0)
        :return Dict[str, float] promoter_feature_dict: a dictionary of features for the promoter
            region upstream of the TSS
        """

        m10_up, m10_down = minus_10_search
        m35_up, m35_down = minus_35_search
        usr_up, usr_down = usr
        dsr_up, dsr_down = dsr

        if strand == 1:
            m10_left, m10_right = tss+m10_up, tss+m10_down
            m35_left, m35_right = tss+m35_up, tss+m35_down
            usr_left, usr_right = tss+usr_up, tss+usr_down
            dsr_left, dsr_right = tss+dsr_up, tss+dsr_down
        else:
            m10_left, m10_right = tss-m10_down, tss-m10_up
            m35_left, m35_right = tss-m35_down, tss-m35_up
            usr_left, usr_right = tss-usr_down, tss-usr_up
            dsr_left, dsr_right = tss-dsr_down, tss-dsr_up

        m10_result = self.motif_search(m10_left, m10_right, strand, minus_10_motif).iloc[0, :]
        m35_result = self.motif_search(m35_left, m35_right, strand, minus_35_motif).iloc[0, :]

        # get the locations of the boxes in terms relative to the TSS
        m10_point = location_to_point(m10_result['left'], m10_result['right'], strand, 'midpoint')
        m10_location = -abs(tss - m10_point)
        m35_point = location_to_point(m35_result['left'], m35_result['right'], strand, 'midpoint')
        m35_location = -abs(tss - m35_point)

        # compute some properties of the spacer between the -10 and -35 elements
        if strand == 1:
            spacer_left, spacer_right = m35_result['right'], m10_result['left']
        else:
            spacer_left, spacer_right = m10_result['right'], m35_result['left']
        spacer_sequence = self.get_sequence(spacer_left, spacer_right, strand)
        spacer_at = 1 - GC(spacer_sequence) / 100

        usr_at = 1 - GC(self.get_sequence(usr_left, usr_right, strand)) / 100
        dsr_sequence = self.get_sequence(dsr_left, dsr_right, strand)
        dsr_ag = (dsr_sequence.count('A') + dsr_sequence.count('G')) / len(dsr_sequence)

        return {'m10_sequence': m10_result['match_sequence'], 'm10_score': m10_result['log_odds'],
                'm10_location': m10_location,
                'm35_sequence': m35_result['match_sequence'], 'm35_score': m35_result['log_odds'],
                'm35_location': m35_location,
                'spacer_sequence': spacer_sequence, 'spacer_length': len(spacer_sequence),
                'spacer_at': spacer_at, 'usr_at': usr_at, 'dsr_ag': dsr_ag}

    def promoter_motif_search(self, tss: int, strand: int, motif_pssm: Dict[str, List[float]],
                              n_up: int = 100, n_down: int = 50, n_best_matches: int = 1
                              ) -> pd.DataFrame:
        """
        Given a TSS and a sequence motif, search for that motif in a customizable range up- and/or
        downstream of the TSS. This function is a specialized wrapper for the more general function
        motif_search

        :param int tss: the location of the TSS around which to search for the motif
        :param int strand: the strand of the TSS; used to define "up-" and "down-" stream directions
        :param Dict[str, List[float]] motif_pssm: the position-specific scoring matrix (PSSM)
            describing the motif to look for in the range of the TSS
        :param int n_up: the number of base pairs upstream of the TSS to include in the search
            range; defaults to 100, not including TSS
        :param int n_down: the number of base pairs downstream of the TSS to include in the search
            range; defaults to 50, not including TSS
        :param int n_best_matches: the number of matched sequences to return; defaults to 1
        :return Dict[str, float] result_dict: a dictionary of result, with keys: match_sequence,
            log_odds, left, and right; left and right are in terms of the reference genome
        """

        if strand == 1:
            left, right = tss - n_up, tss + n_down
        else:
            left, right = tss - n_down, tss + n_up

        return self.motif_search(left, right, strand, motif_pssm, n_best_matches=n_best_matches)

    def motif_search(self, left: int, right: int, strand: int, motif_pssm: Dict[str, List[float]],
                        n_best_matches: int = 1) -> pd.DataFrame:
        
        """
        Given a sequence range to explore and a motif PSM to find, return statistics about the
        closest-matching subsequence within the given range, as a dict. None will be returned if no
        sufficiently-strong match is found.

        :param int left: the left end of the sequence range to search
        :param int right: the right end of the sequence range to search
        :param int strand: the strand of the sequence range to search
        :param Dict[str, List[float]] motif_pssm: the position-specific scoring matrix (PSSM)
            describing the motif to look for in the provided range
        :param int n_best_matches: the number of matched sequences to return; defaults to 1
        :return pd.DataFrame result_table: a DataFrame of results, with columns: match_sequence,
            log_odds, left, and right; left and right are in terms of the reference genome;
            number of rows is determined by n_best_matches
        """
            
        search_sequence = self.get_sequence(left, right, strand)

        motif_length = len(motif_pssm['A'])
        
        motif_pssm = matrix.PositionSpecificScoringMatrix(motif_pssm.keys(), motif_pssm)
    
        scores = motif_pssm.calculate(search_sequence)
            

        result_table = pd.DataFrame(columns=['match_sequence', 'log_odds', 'left', 'right'])

        best_matches_idx = (-scores).argsort()[:n_best_matches]
    
        for start_ind in best_matches_idx:
            # get the left/right ends of the sequence; note that the best_start_ind is w.r.t. the
            # SEQUENCE, which is in CODING order; so it's distance from the START of the sequence
            if strand == 1:
                best_left, best_right = left + start_ind, left + start_ind + motif_length - 1
            else:
                best_right, best_left = right - start_ind, right - start_ind - motif_length + 1
            
            sub_sequence_best = self.get_sequence(best_left, best_right, strand)

            result_table = pd.concat([result_table, pd.DataFrame({'match_sequence': str(sub_sequence_best),
                                                'log_odds': scores[start_ind],
                                                'left': best_left, 'right': best_right}, index=[0])],
                                                ignore_index=True)

        return result_table

    def create_sigma_motifs(self, cache: bool = True) -> dict:
        """
        This function will utilize -10/-35 box information from the tu_table (if present) to create
        motifs for the -10 box, -35 box, and spacer, for each sigma factor present in the tu_table

        :param bool cache: indicates if the resulting motif dict should be cached
        :return dict sigma_motif_dict: a dictionary with sigma factors as keys and further
            dictionaries as values. The value dictionaries contain the following keys:
            - "-10": the -10 box motif for the sigma factor
            - "-35": the -35 box motif
            - "spacer": the spacer motif
            - "spacer_lens": spacer lengths (i.e. the numbers of base pairs in the spacers)
            - "promoter": the full promoter motif (-10/spacer/-35 all as one block)
            - "b10_tss_lens": the distances between the downstream end of the -10 box and the tss;
                DIFFERENT from the spacer_lens, this is the mathematical difference (i.e.
                subtraction) between the positions of the TSS and the downstream -10 box end
        """

        if self.sigma_motifs is not None:
            return self.sigma_motifs

        # get unique sigma factors; assumes that all sigmas are represented with -10 box location
        sigmas_raw = self.tu_table['sigma_factor'].dropna().unique()
        sigmas = []
        for sigma_raw in sigmas_raw:
            sigmas_split = [sig.strip() for sig in sigma_raw.split(',')]
            sigmas += sigmas_split
        sigmas = sorted(list(set(sigmas)), reverse=True)

        # for each sigma factor, compute its motif; use a helper function to avoid code duplication
        def create_box_motif(box_name, sigma_name):
            box_prefix = f'box_{box_name}'
            tu_with_sigma_df = self.tu_table[
                (self.tu_table[f'{box_prefix}_left'].notna()) &
                (self.tu_table['sigma_factor'].str.contains(sigma_name))
            ]
            box_seqs = []
            for _, row in tu_with_sigma_df.iterrows():
                box_seq = self.get_sequence(
                    row[f'{box_prefix}_left'], row[f'{box_prefix}_right'], row.strand
                )
                box_seqs.append(box_seq)
            return create_motif(box_seqs)

        all_sigmas_motif_dict = {}
        for sigma in sigmas:

            # instantiate the sub-dicts for this individual sigma factor with the -10/-35 motifs
            sigma_motif_dict = {
                '-10': create_box_motif('10', sigma),
                '-35': create_box_motif('35', sigma)
            }

            # we want to record information about the motif of the spacer and the entire
            # -10/spacer/-35 region as well
            tu_both_box_df = self.tu_table[
                (self.tu_table['box_10_left'].notna()) &
                (self.tu_table['box_35_left'].notna()) &
                (self.tu_table['sigma_factor'].str.contains(sigma))
            ]
            spacer_seqs = []
            full_prom_seqs = []
            for tu_row in tu_both_box_df.itertuples():
                if tu_row.strand == 1:
                    spacer_l, spacer_r = tu_row.box_35_right + 1, tu_row.box_10_left - 1
                    full_prom_l, full_prom_r = tu_row.box_35_left, tu_row.box_10_right
                else:
                    spacer_l, spacer_r, = tu_row.box_10_right + 1, tu_row.box_35_left - 1
                    full_prom_l, full_prom_r = tu_row.box_10_left, tu_row.box_35_right
                # trxA weird case where there is no spacer at all (-10 and -35 overlap)
                if spacer_l < spacer_r:
                    spacer_seqs.append(self.get_sequence(spacer_l, spacer_r, tu_row.strand))
                full_prom_seqs.append(self.get_sequence(full_prom_l, full_prom_r, tu_row.strand))
            sigma_motif_dict['spacer'] = create_motif(spacer_seqs)
            spacer_lens = [len(spacer_seq) for spacer_seq in spacer_seqs]
            sigma_motif_dict['spacer_lens'] = spacer_lens
            sigma_motif_dict['promoter'] = create_motif(full_prom_seqs)

            # finally, get information about the -10 to TSS lens
            tu_b10_tss_df = self.tu_table[
                (self.tu_table['box_10_left'].notna()) &
                (self.tu_table['tss'].notna()) &
                (self.tu_table['sigma_factor'].str.contains(sigma))
            ]
            b10_tss_lens = []
            for tu_row in tu_b10_tss_df.itertuples():
                if tu_row.strand == 1:
                    b10_tss_len = tu_row.tss - tu_row.box_10_right
                else:
                    b10_tss_len = tu_row.box_10_left - tu_row.tss
                b10_tss_lens.append(b10_tss_len)
            sigma_motif_dict['b10_tss_lens'] = b10_tss_lens

            # add the completed motif dict for this sigma factor to the overall lookups
            all_sigmas_motif_dict[sigma] = sigma_motif_dict

        if cache:
            self.sigma_motifs = all_sigmas_motif_dict
        return all_sigmas_motif_dict

    def predict_tss(self, gene_row: pd.Series, n_upstream: int = 250,
                    box_weights: Tuple[float, float] = (0.6, 0.4),
                    spacer_len_sd_multiplier: float = 1.5,
                    b10_tss_sd_multiplier: float = 1.0,
                    n_up_shape_ref: int = 15, n_down_shape_ref: int = 5, n_shape_discern: int = 5,
                    true_tss_data: dict = None, plot: bool = False, plot_box_scores: bool = False,
                    plot_spacer_lens: bool = False):
        """
        This function attempts to locate a putative TSS upstream of a given gene row. It does so
        using a combination of motif matching of the -10/-35 box, taking into account spacer length,
        for all known sigma factors. Sigma factor motifs are calculated in the background. DNA
        shape is also utilized (shape reference is also computed in the background).

        :param pd.Series gene_row: a gene information row, taken from the gene_table, for which
            a TSS prediction is desired
        :param int n_upstream: # of bps upstream of the gene in which to search for a TSS
        :param Tuple[float, float] box_weights: the weights with which to consider the -10/-35
            scores (respectively) when assessing the score for a particular -10/spacer/-35 set
        :param float spacer_len_sd_multiplier: multiplier for standard deviation of spacer lengths
            to be used to determine wiggle around median spacer length
        :param float b10_tss_sd_multiplier: multiplier for standard deviation of -10 box to TSS
            lengths to be used to determine wiggle around median when testing TSS locations with
            shape downstream of -10 box
        :param int n_up_shape_ref: # of bps upstream of known TSSes in which to get reference shape
        :param int n_down_shape_ref: # of bps downstream of known TSSes in which to get ref shape
        :param int n_shape_discern: # of top promoter matches - from the motif search - for which
            to use the reference TSS shape to discern the best
        :param dict true_tss_data: a dictionary with keys 'tss', '-10', '-35' indicating actual
            positions of these features, if known, to be compared to the prediction. '-10' and '-35'
            entries are assumed to be the DOWNSTREAM edges
        :param bool plot: indicates if a plot of the prediction score should be outputted
        :param bool plot_box_scores: indicates if the scores for the individual box motifs should
            be plotted
        :param bool plot_spacer_lens: indicates if the individual score trajectories for different
            spacer lengths should be plotted; by default, the best score for each downstream -10 box
            end will be plotted as the score trajectory
        """

        # get the absolute genome positions in which we will search; orient up- to downstream
        strand = gene_row.strand
        if strand == 1:
            search_l, search_r = gene_row.left - n_upstream, gene_row.left
            search_positions = np.arange(search_l, search_r + 1)
            # record the name (left or right) that corresponds to upstream
            up_col, down_col = 'left', 'right'
        else:
            search_l, search_r = gene_row.right, gene_row.right + n_upstream
            search_positions = np.flip(np.arange(search_l, search_r + 1))
            up_col, down_col = 'right', 'left'
        # we want to have a strand-agnostic way to identify a position in the search range
        search_idx_dict = {pos: i for i, pos in enumerate(search_positions)}
        idx_to_pos_dict = {v: k for k, v in search_idx_dict.items()}

        # calculate the sigma factor motifs (returns a cached result if already called in session)
        all_sigma_motif_dict = self.create_sigma_motifs()

        # do the motif matching of all sigma factor motifs
        all_box_10_match_dict = {}
        all_box_35_match_dict = {}
        all_prom_match_dict = {}
        best_prom_match_dict = {}
        for sigma, sigma_motif_dict in all_sigma_motif_dict.items():

            # match the -10 box motif and -35 box motifs throughout the search range
            box_10_motif = sigma_motif_dict['-10']
            box_10_match_df = self.motif_search(
                search_l, search_r, strand, box_10_motif.pssm,
                n_best_matches=(n_upstream - len(box_10_motif))
            )
            box_35_motif = sigma_motif_dict['-35']
            box_35_match_df = self.motif_search(
                search_l, search_r, strand, box_35_motif.pssm,
                n_best_matches=(n_upstream - len(box_35_motif))
            )

            # add columns to the match DataFrames that locate the boxes in up/downstream manner
            for match_df in [box_10_match_df, box_35_match_df]:
                match_df['up'] = match_df[up_col].apply(lambda pos: search_idx_dict[pos])
                match_df['down'] = match_df[down_col].apply(lambda pos: search_idx_dict[pos])
                # duplicate the down column so we can merge on it but still retain distinct down
                # columns for each box after the merge
                match_df['search_idx'] = match_df['down'].copy()

            # record the individual box match tables in case we have to plot them
            all_box_10_match_dict[sigma] = box_10_match_df
            all_box_35_match_dict[sigma] = box_35_match_df

            # now we can merge the match DataFrames based on their search range index (downstream)
            merged_match_df = box_10_match_df.merge(box_35_match_df, how='inner', on='search_idx',
                                                    suffixes=("_b10", "_b35"))
            merged_match_df = merged_match_df.sort_values(by='search_idx')

            # extract the spacer lengths and set up the different spacer lengths to try out
            # NOTE: the spacer lengths include JUST the bps in the spacer
            spacer_lens = sigma_motif_dict['spacer_lens']
            spacer_med, spacer_sd = np.median(spacer_lens), np.std(spacer_lens)
            spacer_lens_to_consider = np.arange(
                int(spacer_med - spacer_sd * spacer_len_sd_multiplier),
                int(spacer_med + spacer_sd * spacer_len_sd_multiplier)
            )

            # iterate through each potential spacer length and get all promoter matching scores
            # for that spacer length throughout the search range
            prom_match_rows = []
            for spacer_len in spacer_lens_to_consider:
                # this is important; we want to offset the merged match so that we track through a
                # putative b10/b35 pair with a given spacer length; the search_idx is set up to be
                # the DOWNSTREAM edge of each box; thus, to ensure a given spacer, we need the INDEX
                # offset to be the spacer length PLUS the length of the -10 motif (because the
                # spacer is from the downstream edge of the -35 to the UPSTREAM edge of the -10)
                # also on top of that, keep in mind that the spacer length is  # of bps IN SPACER
                for b35_row, b10_row in zip(
                        # zip takes the shorter list, so the first list is auto cut off, as we want
                        merged_match_df.itertuples(),
                        merged_match_df.iloc[
                        (spacer_len + len(box_10_motif)):].itertuples()
                ):
                    combined_score = (
                        b10_row.log_odds_b10 * box_weights[0] +
                        b35_row.log_odds_b35 * box_weights[1]
                    )
                    prom_match_rows.append({
                        'left_b10': b10_row.left_b10, 'right_b10': b10_row.right_b10,
                        'left_b35': b35_row.left_b35, 'right_b35': b35_row.right_b35,
                        'up_b10': b10_row.up_b10, 'down_b10': b10_row.down_b10,
                        'up_b35': b35_row.up_b35, 'down_b35': b35_row.down_b35,
                        'b10_score': b10_row.log_odds_b10, 'b35_score': b35_row.log_odds_b35,
                        'score': combined_score, 'spacer_len': spacer_len
                    })

            prom_match_df = pd.DataFrame(prom_match_rows).sort_values(by='score', ascending=False)
            all_prom_match_dict[sigma] = prom_match_df

            # pull out the best prediction for the downstream edge of the -10 and record it
            best_b10_match = prom_match_df.iloc[0]['down_b10']
            best_prom_match_dict[sigma] = idx_to_pos_dict[best_b10_match]

        # also create the reference TSS shape; this will also populate the TSS search area
        # stats and raw shape in the search areas, which we will need later
        tss_reference_shape = self.tss_reference_shape(
            n_upstream=n_up_shape_ref, n_downstream=n_down_shape_ref, n_upstream_search=n_upstream
        )

        # get the shape score for the full range; since we called tss_reference_shape, that, as a
        # side effect (I know, I know), should have populated self._tss_search_shape_dict, where
        # we can get the current gene ID looked up; NOTE .name, for series, is the index (locus tag)
        search_range_shape_df = self._tss_search_shape_dict[gene_row.name]

        # score each of the potential TSS in the search range
        shape_ref_size = (n_up_shape_ref + n_down_shape_ref + 1)
        putative_tss_shape_scores = {}
        for up_idx in range(search_range_shape_df.shape[0] - (shape_ref_size - 1)):

            # compare the z-scored shapes of the ground truth and this putative range
            shape_range_df = search_range_shape_df.iloc[up_idx:(up_idx+shape_ref_size)]
            shape_mae_dict = {}
            for shape_col in tss_reference_shape.columns:
                true_col = tss_reference_shape[shape_col]
                putative_col = shape_range_df[shape_col]
                shape_mae_dict[shape_col] = median_absolute_error(true_col, putative_col)
            putative_tss_shape_scores[up_idx + n_up_shape_ref] = shape_mae_dict

        # now we can get the shape scores just for the ranges downstream from the -10 boxes of the
        # n_shape_discern motif matches
        shape_match_df_dict = {}
        for sigma, sigma_match_df in all_prom_match_dict.items():

            # get the actual distribution of -10 to TSS lengths for this sigma
            b10_tss_lens = self.sigma_motifs[sigma]['b10_tss_lens']
            b10_tss_med, b10_tss_sd = np.median(b10_tss_lens), np.std(b10_tss_lens)
            b10_tss_lens_to_consider = np.arange(
                int(b10_tss_med - b10_tss_sd * b10_tss_sd_multiplier),
                int(b10_tss_med + b10_tss_sd * b10_tss_sd_multiplier)
            )

            tss_match_rows = []
            for b10_tss_len in b10_tss_lens_to_consider:

                # for this particular -10 box to TSS length, look downstream of the best motif
                # matches to score the shape around those putative TSSes
                for _, match_row in sigma_match_df.iloc[:n_shape_discern].iterrows():

                    # the match row is oriented relative to the box 10 location; use the index of
                    # downstream end of the -10 to select from the shape range
                    down_b10_idx = match_row.down_b10

                    # can use no strand awareness here since the up/down orientation within the
                    # search range has already handled that for us; can get our shape scores too
                    # can't look to far ahead if we're near the edge of our search range
                    putative_tss_idx = down_b10_idx + b10_tss_len
                    # TODO: this kind of silently defines 5' UTR length, can this be better?
                    putative_tss_shape_score_dict = putative_tss_shape_scores.get(putative_tss_idx)
                    if putative_tss_shape_score_dict is None:
                        continue

                    # record all the information about this TSS match
                    tss_match_rows.append({
                        **dict(match_row),
                        'b10_tss_len': b10_tss_len, 'tss_idx': putative_tss_idx,
                        'tss': idx_to_pos_dict[putative_tss_idx],
                        **putative_tss_shape_score_dict,
                        'shape_med': np.median(list(putative_tss_shape_score_dict.values()))
                    })

            tss_match_df = pd.DataFrame(tss_match_rows)
            shape_match_df_dict[sigma] = tss_match_df

        # plot the scores of the matches (along with the actual TSS information if requested)
        if plot:

            # set up colors for actual features
            label_to_col = {'tss': 'gray', 'b10': 'blue', 'b35': 'green'}

            _, axs = plt.subplots(len(all_sigma_motif_dict), 1, figsize=(15, 15), sharex='col')

            for sigma, ax in zip(all_sigma_motif_dict.keys(), axs):

                # plot the individual -10/-35 motif match scores if asked
                if plot_box_scores:
                    sns.lineplot(x='search_idx', y='log_odds', data=all_box_10_match_dict[sigma],
                                 label='-10 Box', color='blue', ax=ax)
                    sns.lineplot(x='search_idx', y='log_odds', data=all_box_35_match_dict[sigma],
                                 label='-35 Box', color='green', ax=ax)

                # plot the combined score results, either separating the spacer lengths or not
                if plot_spacer_lens:
                    sns.lineplot(x='down_b10', y='score', data=all_prom_match_dict[sigma],
                                 hue='spacer_len', ax=ax)
                else:
                    all_spacer_prom_match_df = all_prom_match_dict[sigma]
                    best_by_down_rows = []
                    for down_b10, down_b10_df in all_spacer_prom_match_df.groupby('down_b10'):
                        best_by_down_rows.append(
                            down_b10_df.sort_values(by='score', ascending=False).iloc[0]
                        )
                    best_by_b10_df = pd.DataFrame(best_by_down_rows)
                    sns.lineplot(x='down_b10', y='score', data=best_by_b10_df, ax=ax,
                                 label='combined', color='orange')

                ax.set_ylabel(f'{sigma} Score', fontsize=14)

                # plot the shape results; the median of the median absolute errors; use a twin x
                ax_shape = ax.twinx()
                sigma_shape_match_df = shape_match_df_dict[sigma]
                sns.scatterplot(x='tss_idx', y='shape_med', data=sigma_shape_match_df, ax=ax_shape,
                                label='Shape Median MAE', color='black')
                ax_shape.set_ylabel('Shape Median MAE', fontsize=14)

                if true_tss_data is not None:
                    for label, locs in true_tss_data.items():
                        for loc in locs:
                            ax.axvline(search_idx_dict[loc], linestyle='--',
                                       label=f'Actual {label}', color=label_to_col[label],
                                       zorder=-10)

                # also get the best predicted promoter match from the combined scores and plot that
                ax.axvline(search_idx_dict[best_prom_match_dict[sigma]], color='red',
                           label='Predicted -10', zorder=-10)

                ax.legend(fontsize=14, bbox_to_anchor=(1, 0.5), loc='center left')

            axs[-1].set_xticks([i for i in range(len(search_positions)) if i % 10 == 0])
            axs[-1].set_xticklabels([pos for i, pos in enumerate(search_positions) if i % 10 == 0],
                                    rotation=45, ha='right')
            axs[-1].set_xlabel('Genome Position (bp)', fontsize=14)
            axs[0].set_title(f"{gene_row['name']}")

        return all_prom_match_dict, shape_match_df_dict

    # --- ONE-HOT ENCODING ---

    def encode_boxes(self, shape: bool = True, sequence: bool = True, length_is: int = 6
                     ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Encode the shape and sequence of the -10/-35 boxes in a position-dependent manner.

        :param bool shape: indicates if shape should be calculated for the -10/-35 boxes
        :param bool sequence: indicates if one-hot encoded sequence should be included
        :param int length_is: specifies the length filter to use for -10/-35 boxes
        :return List[pd.DataFrame, pd.DataFrame] box_pos_dfs: the DataFrames containing position-
            specific shape and/or sequence information for the -10/-35 boxes
        """

        # filter out only the TUs that have box information with the right length
        tu_table_to_use = self.tu_table[
            ((self.tu_table.box_10_right - self.tu_table.box_10_left) == length_is - 1) &
            ((self.tu_table.box_35_right - self.tu_table.box_35_left) == length_is - 1)
        ]

        # set up a separate table that just has the locations of the boxes
        box_df_rows = []
        box_df_index = []
        for tu_row in tu_table_to_use.itertuples():
            b10_l, b10_r = tu_row.box_10_left, tu_row.box_10_right
            b10_df_row = {'left': b10_l, 'right': b10_r, 'strand': tu_row.strand}
            box_df_rows.append(b10_df_row)
            box_df_index.append(f'{tu_row.Index}_b10')
            b35_l, b35_r = tu_row.box_35_left, tu_row.box_35_right
            b35_df_row = {'left': b35_l, 'right': b35_r, 'strand': tu_row.strand}
            box_df_rows.append(b35_df_row)
            box_df_index.append(f'{tu_row.Index}_b35')
        box_location_df = pd.DataFrame(box_df_rows, index=box_df_index)

        if shape:

            # add some wiggle room to the box DF to ensure we don't have edge effects
            box_loc_df_shape = box_location_df.copy()
            box_loc_df_shape['left'] = box_loc_df_shape['left'].copy() - 5
            box_loc_df_shape['right'] = box_loc_df_shape['right'].copy() + 5

            # use the get_dna_shape_for_features for all -10/35 box at once
            box_shape_dict = self.get_dna_shape_for_features(box_loc_df_shape)

            # make an internal helper to handle lookups for -10 and -35 boxes
            def make_shape_row(box_id, box_tu_row):
                box_shape_df_full = box_shape_dict[box_id]
                if 'b10' in box_id:
                    box_pos = np.arange(box_tu_row.box_10_left, box_tu_row.box_10_right + 1)
                else:
                    box_pos = np.arange(box_tu_row.box_35_left, box_tu_row.box_35_right + 1)
                # reverse the order for reverse strand TU so we have all boxes oriented
                if box_tu_row.strand == -1:
                    box_pos = np.flip(box_pos)
                box_shape_df = box_shape_df_full.loc[box_pos]
                box_shape_df.index = np.arange(1, length_is + 1)
                # this little nifty block of code from StackOverflow gives us the column names
                # we want (i.e. {shape_name}_{position index} in a single row for later
                box_shape_row = box_shape_df.stack()
                box_shape_row.index = box_shape_row.index.map('{0[1]}_{0[0]}'.format)
                return box_shape_row

            # now iterate again through the TU tables, and pick up the box shapes for the TU,
            # and turn them into flat vectors in a new DataFrame
            b10_shape_rows = []
            b10_shape_index = []
            b35_shape_rows = []
            b35_shape_index = []
            for tu_row in tu_table_to_use.itertuples():
                b10_id, b35_id = f'{tu_row.Index}_b10', f'{tu_row.Index}_b35'
                if b10_id in box_shape_dict:
                    b10_shape_row = make_shape_row(b10_id, tu_row)
                    b10_shape_rows.append(b10_shape_row)
                    b10_shape_index.append(tu_row.Index)
                if b35_id in box_shape_dict:
                    b35_shape_row = make_shape_row(b35_id, tu_row)
                    b35_shape_rows.append(b35_shape_row)
                    b35_shape_index.append(tu_row.Index)

            b10_shape_df = pd.DataFrame(b10_shape_rows, index=b10_shape_index)
            b35_shape_df = pd.DataFrame(b35_shape_rows, index=b35_shape_index)

        else:
            b10_shape_df = pd.DataFrame()
            b35_shape_df = pd.DataFrame()

        if sequence:

            box_one_hot_df = self.one_hot_encode_sequences(box_location_df)
            # re-separate the -10 and -35 box sequences
            b10_oh_df = box_one_hot_df[box_one_hot_df.index.str.contains('b10')]
            b10_oh_df.index = [idx.split('_')[0] for idx in b10_oh_df.index]
            b35_oh_df = box_one_hot_df[box_one_hot_df.index.str.contains('b35')]
            b35_oh_df.index = [idx.split('_')[0] for idx in b35_oh_df.index]

        else:
            b10_oh_df = pd.DataFrame()
            b35_oh_df = pd.DataFrame()

        # prepare the final DataFrames; by using an outer merge, we ensure that we keep the max
        # keys; this is a bit of a weird way of handling that one of these may be empty
        b10_final_df = b10_shape_df.merge(b10_oh_df, how='outer', left_index=True,
                                          right_index=True)
        b35_final_df = b35_shape_df.merge(b35_oh_df, how='outer', left_index=True,
                                          right_index=True)

        return b10_final_df, b35_final_df

    def one_hot_encode_sequences(self, feature_df: pd.DataFrame) -> pd.DataFrame:
        """
        Given any DataFrame containing left/right/strand information, return a DataFrame containing
        one-hot encoded sequences for the features in that DataFrame.

        :param pd.DataFrame feature_df: the feature table, as a DataFrame, for which one-hot
            encoding should be done en masse
        :return pd.DataFrame one_hot_df: the one-hot encoded feature table
        """

        one_hot_rows = []
        for row in feature_df.itertuples():
            one_hot_rows.append(self.one_hot_encode_sequence(row.left, row.right, row.strand))
        one_hot_df = pd.DataFrame(one_hot_rows, index=feature_df.index)
        return one_hot_df

    def one_hot_encode_sequence(self, left: int, right: int, strand: int) -> pd.Series:
        """
        Given a feature location range, return a one-hot encoded vector of that sequence. The
        alphabetical (ACGT) nucleotide order will be used.

        :param int left: the left end of the sequence to encode
        :param int right: the right end of the sequence to encode
        :param int strand: the strand end of the sequence to encode
        :return pd.Series one_hot_sequence: a pandas Series containing the one-hot-encoded sequence
        """

        sequence = self.get_sequence(left, right, strand)

        # not using pandas' get_dummies method ensures that we account for sequences that don't
        # have all 4 nucleotides, and also allows us to prepare a 1-dimensional one-hot vector
        one_hot_names = []
        one_hots = []
        for i, seq_base in enumerate(sequence):
            for base in 'ACGT':
                one_hot_names.append(f'{i}_{base}')
                if base == seq_base:
                    one_hots.append(1)
                else:
                    one_hots.append(0)

        one_hot_sequence = pd.Series(one_hots, index=one_hot_names)
        return one_hot_sequence

    def one_hot_encode_tu_promoters(self, local_align: bool = False,
                                    n_upstream: int = 50, n_downstream: int = 10,
                                    tss_radius: int = 3, box_10_radius: int = 6,
                                    box_35_radius: int = 6, shift_odd_len_upstream: bool = True
                                    ) -> pd.DataFrame:
        """
        High-level utility function to prepare a data matrix with all TUs (tu_table must be set)
        converted into one-hot sequences. Uses one_hot_encode_sequence under the hood.

        :param bool local_align: indicates if the one-hot encoding should include local regions
            centered around the specific annotated locations of key promoter features: TSS (CRE),
            -10 box, -35 box
        :param int n_upstream: the number of bases upstream from a TSS to include; only used if
            local_align is False; defaults to 50
        :param int n_downstream: the number of bases downstream from a TSS to include; only used if
            local_align is False, defaults to 10
        :param int tss_radius: the number of base pairs up and downstream of the TSS to one-hot
            encode (NOTE: this parameter is in ONE direction, so the final range will be double the
            size of this parameter, plus the TSS itself); only used if local_align is True
        :param int box_10_radius: the number of base pairs up and downstream of the center of the
            -10 box to one-hot encode; this box is typically a hexamer, so the center point will be
            between the 2 central nucleotides; if a box annotation has odd length, it will be
            shifted to ensure alignment with the even-length boxes. The direction of this shift can
            be controlled with shift_odd_len_upstream. Only used if local_align is True
        :param int box_35_radius: the number of base pairs up and downstream of the center of the
            -35 box to one-hot encode; this box is typically a hexamer, so the center point will be
            between the 2 central nucleotides; if a box annotation has odd length, it will be
            shifted to ensure alignment with the even-length boxes. The direction of this shift can
            be controlled with shift_odd_len_upstream. Only used if local_align is True
        :param bool shift_odd_len_upstream: in order to align odd length promoter boxes, indicates
            if they should be shifted 1 bp upstream (i.e. their "center" is 0.5 bp DOWNSTREAM of the
            central base pair) or 1 bp downstream (i.e. their "center" is 0.5 bp UPSTREAM of the
            central base pair); only used if local_align is True
        :return pd.DataFrame one_hot_tu_df: a DataFrame containing rows with one-hot-encoded TSS
            promoter sequences for all TUs with known TSS
        """

        if self.tu_table.empty:
            raise ValueError('No TU table available.')

        # create a utility that will help us rename one-hot row indices from one_hot_encode_sequence
        def reindex_one_hot(old_names, offset, prefix=''):
            new_names = []
            for old_name in old_names:
                ind, base = old_name.split('_')
                new_ind = int(ind) - offset
                new_name = f'{prefix}{new_ind}_{base}'
                new_names.append(new_name)
            return new_names

        if local_align:
            one_hot_index = []
            one_hot_rows = []
            tu_table_local_align = self.tu_table[
                (self.tu_table['tss'].notna()) &
                (self.tu_table['box_10_left'].notna()) &
                (self.tu_table['box_35_left'].notna())
            ]
            for tu_row in tu_table_local_align.itertuples():

                tss_l, tss_r = tu_row.tss - tss_radius, tu_row.tss + tss_radius

                # define a quick utility to appropriately align the (typically) hexamer boxes
                def get_box_mid(box_l, box_r, strand):
                    box_mid = location_to_point(box_l, box_r, strand, 'midpoint')
                    # adjust the midpoint if the box has an odd length
                    if (box_r - box_l) % 2 == 1:
                        if shift_odd_len_upstream:
                            box_mid = box_mid + 0.5 * strand
                        else:
                            box_mid = box_mid - 0.5 * strand
                    return box_mid

                # ceiling/floor handles the case of even length
                b10_mid = get_box_mid(tu_row.box_10_left, tu_row.box_10_right, tu_row.strand)
                b10_l, b10_r = np.ceil(b10_mid - box_10_radius), np.floor(b10_mid + box_10_radius)
                b35_mid = get_box_mid(tu_row.box_35_left, tu_row.box_35_right, tu_row.strand)
                b35_l, b35_r = np.ceil(b35_mid - box_35_radius), np.floor(b35_mid + box_35_radius)

                tss_oh = self.one_hot_encode_sequence(tss_l, tss_r, tu_row.strand)
                tss_oh.index = reindex_one_hot(tss_oh.index, tss_radius, prefix='tss_')
                box_10_oh = self.one_hot_encode_sequence(b10_l, b10_r, tu_row.strand)
                box_10_oh.index = reindex_one_hot(box_10_oh.index, box_10_radius, prefix='-10_')
                box_35_oh = self.one_hot_encode_sequence(b35_l, b35_r, tu_row.strand)
                box_35_oh.index = reindex_one_hot(box_35_oh.index, box_35_radius, prefix='-35_')
                full_oh_row = pd.concat([tss_oh, box_10_oh, box_35_oh])

                one_hot_index.append(tu_row.Index)
                one_hot_rows.append(full_oh_row)

            one_hot_tu_df = pd.DataFrame(one_hot_rows, index=one_hot_index)

        else:
            one_hot_index = []
            one_hot_rows = []
            tu_table_tss = self.tu_table[self.tu_table['tss'].notna()]
            for tu_row in tu_table_tss.itertuples():
                one_hot_index.append(tu_row.Index)
                if tu_row.strand == 1:
                    left, right = tu_row.tss - n_upstream, tu_row.tss + n_downstream
                else:
                    left, right = tu_row.tss - n_downstream, tu_row.tss + n_upstream
                one_hot_row = self.one_hot_encode_sequence(left, right, tu_row.strand)
                one_hot_rows.append(one_hot_row)

            one_hot_tu_df = pd.DataFrame(one_hot_rows, index=one_hot_index)

            # adjust the column names to reflect the locations relative to the TSS
            one_hot_tu_df.columns = reindex_one_hot(one_hot_tu_df.columns, n_upstream)

        return one_hot_tu_df

    # --- DNA shape-related functionalities ---

    def get_dna_shape_for_features(self,feature_table: pd.DataFrame
                                   ) -> Union[Dict[str, pd.DataFrame], None]:
        """
        Given a feature table (must have left, right, and strand columns at a minimum), compute the
            DNA shape parameters for the sequence of each feature in the table.
        Uses 5mer shape lookup computed with the DNAshapeR (https://rdrr.io/bioc/DNAshapeR/) package;
            see reference: https://academic.oup.com/bioinformatics/article/32/8/1211/1744520
        14 different shape parameters are returned: Opening, Rise, Stretch, Electrostatic
            Potential (EP), Tilt, Shear, Propeller Twist (ProT), Buckle, Shift, Helical Twist (HelT)
            Stagger, Slide, Minor Groove Width (MGW), and Roll. See the DNAshapeR package docs
            for detail about what each of these parameters means structurally

        :param pd.DataFrame feature_table: a Bitome-formatted feature table (i.e. with left, right
            and strand columns at a minimum) containing features for which the DNA shape of the
            features' sequences should be calculated.
        :return Dict[str, pd.DataFrame] shape_df_dict: a dictionary containing the DNA shape
            results; keys are the locus tags of the features from the feature table, and values
            are DataFrames of the form: feature sequence base pairs (rows) x shape features (cols)
        """
        return {idx: self.get_dna_shape(row.left, row.right, row.strand)\
                for idx, row in feature_table.iterrows()}

    def get_dna_shape(self, left: int, right: int, strand: int) -> pd.DataFrame:
        """
        Given a location in left, right, strand format, use the locally stored fivemer shape
        lookup to construct and return a DataFrame containing base-by-base shape in this region
        
        NOTE: because shape is 5mer-based and each shape value is centered within a 5mer, 
        this function will automagically and temporarily add the necessary overhang to 
        achieve the shape for the exact given sequence
        
        :param int left: the left end of the range for which to compute shape
        :param int right: the right end of the range for which to compute shape
        :param int strand: the strand from which to extract sequence for shape computation
        """
        
        # add the +/- 2 padding to ensure that we return shape as requested
        seq_for_shape = self.get_sequence(left - 2, right + 2, strand)
        seq_fivemers = [str(seq_for_shape[i:i+5]) for i in range(len(seq_for_shape) - 5 + 1)]
        shape_df = self.fivemer_shape_lookup.loc[seq_fivemers]
        
        # change the index of the raw shape DF to reflect the absolute positions given initially
        pos_index = list(range(left, right + 1))
        if strand == -1:
            pos_index = pos_index[::-1]
        shape_df.index = pos_index
        
        return shape_df
    
    def tss_search_shape_stats(self, n_upstream: int = 250, cache: bool = True):
        """
        A component of TSS finding machinery
        This utility function automatically computes the shape in the search ranges upstream of
        ALL genes

        :param int n_upstream: # of bps upstream of a gene start in which to compute shape
        :param bool cache: indicates if the resulting dict should be saved to the
            _tss_search_shape_stats attribute; NOTE: the cache will be busted if a different
            n_upstream is requested
        :return pd.DataFrame tss_search_shape_stat_df: a DataFrame containing mean/std columns for
            each shape feature, as computed within the search ranges
        """

        # return the cache if it is set
        if self._tss_search_shape_stats is not None:
            # bust the cache if n_upstream is different than the last request
            if n_upstream != self._tss_search_n_upstream:
                self._tss_search_shape_stats = None
            else:
                return self._tss_search_shape_stats

        # get the search ranges above all genes; skip if there's overlap as in some operons
        tss_search_shape_rows = []
        for gene_row in self.gene_table.itertuples():
            # find gene start and next upstream gene
            if gene_row.strand == 1:
                gene_start = gene_row.left
                upstream_gene_df = self.gene_table[self.gene_table['left'] < gene_start]
                upstream_gene_df = upstream_gene_df.sort_values(by='left', ascending=False)
                if not upstream_gene_df.empty and gene_row.left < upstream_gene_df.iloc[0, :].right:
                    continue
                # otherwise, define a search range in which to look for TSS
                else:
                    search_l, search_r = gene_start - n_upstream, gene_start
            else:
                gene_start = gene_row.right
                upstream_gene_df = self.gene_table[self.gene_table['right'] > gene_start]
                upstream_gene_df = upstream_gene_df.sort_values(by='right')
                if not upstream_gene_df.empty and gene_row.right > upstream_gene_df.iloc[0, :].left:
                    # ignore if next gene is on top of this one
                    continue
                # otherwise, define a search range in which to look for TSS
                else:
                    search_l, search_r = gene_start, gene_start + n_upstream
            tss_search_shape_rows.append({
                'left': search_l,
                'right': search_r,
                'strand': gene_row.strand,
                'locus_tag': gene_row.Index
            })

        # compute shape for all of these search areas
        tss_search_get_shape_df = pd.DataFrame(tss_search_shape_rows).set_index('locus_tag')
        tss_search_shape_df_dict = self.get_dna_shape_for_features(tss_search_get_shape_df)
        tss_search_shape_df = pd.concat(list(tss_search_shape_df_dict.values()))

        # compute the mean/SD and return/cache
        tss_search_shape_stat_df = pd.DataFrame(data={
            'mean': tss_search_shape_df.mean(axis=0),
            'std': tss_search_shape_df.std(axis=0)
        })
        if cache:
            self._tss_search_shape_stats = tss_search_shape_stat_df
            self._tss_search_shape_dict = tss_search_shape_df_dict
            self._tss_search_n_upstream = n_upstream

        return tss_search_shape_stat_df

    def tss_reference_shape(self, n_upstream: int = 15, n_downstream: int = 5,
                            n_upstream_search: int = 250, return_raw: bool = False,
                            cache: bool = True
                            ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, pd.DataFrame]]:
        """
        Prepare a reference shape fingerprint for the regions around known TSSes from the tu_table.
        Can be used to find TUs based (in part) on shape fingerprint in the TSS region.

        :param int n_upstream: # of bp upstream of the TSS to include in the shape fingerprint
        :param int n_downstream: # of bp downstream of the TSS to include in the shape fingerprint
        :param int n_upstream_search: # of bp upstream of gene WITHOUT TSS to use for creating the
            z-score stats necessary to scale the reference shape fingerprint
        :param bool return_raw: indicates if the raw (un-averaged) reference shape DF should be
            returned
        :param bool cache: indicates if the resulting reference shape should be saved to the
            private attribute _tss_reference_shape
        :return pd.DataFrame reference_shape_df: a DataFrame containing the reference shape for the
            known TSS regions, with shape features as columns and base positions relative to TSS as
            the row index
        """

        # use the cache if it is already set; bust if the n_upstream/downstream differ
        if self._tss_reference_shape is not None:
            if n_upstream != self._tss_reference_n_up or n_downstream != self._tss_reference_n_down:
                self._tss_reference_shape = None
            elif n_upstream_search != self._tss_search_n_upstream:
                self._tss_search_shape_dict = None
                self._tss_search_shape_stats = None
            else:
                return self._tss_reference_shape

        # set up the location information for which to compute shape around known TSSes
        rows_for_shape = []
        to_get_shape_index = []
        for tu_row in self.tu_table.itertuples():
            if pd.notna(tu_row.tss):
                if tu_row.strand == 1:
                    l, r = tu_row.tss - n_upstream, tu_row.tss + n_downstream
                else:
                    l, r = tu_row.tss - n_downstream, tu_row.tss + n_upstream
                rows_for_shape.append({
                    'left': int(l),
                    'right': int(r),
                    'strand': int(tu_row.strand)
                })
                to_get_shape_index.append(tu_row.Index)
        to_get_shape_df = pd.DataFrame(rows_for_shape, index=to_get_shape_index)

        # compute shape all at one for these TSS regions
        shape_df_dict = self.get_dna_shape_for_features(to_get_shape_df)
        # add in the bp information relative to the TSS
        for tu_id, tu_shape_df in shape_df_dict.items():
            tu_shape_df['bp'] = np.arange(-n_upstream, n_downstream + 1)
        # combine into a single giant table so we can compute averages by bp across all TSS
        combined_tss_shape_df = pd.concat(list(shape_df_dict.values()))

        # create the reference shape dataframe
        bp_avg_rows = []
        bp_index = []
        for bp, bp_df in combined_tss_shape_df.groupby('bp'):
            bp_avg_rows.append(bp_df.mean(axis=0))
            bp_index.append(bp)
        tss_shape_reference_df = pd.DataFrame(bp_avg_rows, index=bp_index)
        tss_shape_reference_df = tss_shape_reference_df.drop(columns='bp')

        # now we want to z-score these values (if asked); get the search area reference; this will
        # be cached also, unless the n_upstream_search is different
        z_score_stat_df = self.tss_search_shape_stats(n_upstream=n_upstream_search)
        for shape_col in tss_shape_reference_df.columns:
            mean, std = z_score_stat_df.loc[shape_col, ['mean', 'std']]
            tss_shape_reference_df[shape_col] = (
                (tss_shape_reference_df[shape_col].copy() - mean) / std
            )

        if cache:
            self._tss_reference_shape = tss_shape_reference_df
            self._tss_reference_n_up, self._tss_reference_n_down = n_upstream, n_downstream
        if return_raw:
            return tss_shape_reference_df, combined_tss_shape_df
        else:
            return tss_shape_reference_df

    # --- EXPRESSION WORKFLOW FUNCTIONS ---
    def expression_distribution(self) -> List[plt.Axes]:
        """
        Plot the gene expression distribution from expression data (if any) from this Bitome
        2 plots are generated:
            - the raw, un-transformed expression values
            - the log-transformed expression values

        Uses the expression_column value set upon creation of the Bitome

        :return List[plt.Axes] axes: returns the two Axes corresponding to the 2 subplots (raw and
            log-transformed expression, respectively)
        """

        _, (ax_raw, ax_log) = plt.subplots(1, 2, figsize=(15, 6))

        # plot the raw expression distribution
        sns.boxplot(x=self._exp_col, data=self.gene_table, fliersize=0, ax=ax_raw)
        sns.stripplot(x=self._exp_col, data=self.gene_table, color='gray', ax=ax_raw)
        ax_raw.set_xlabel(f'{self._exp_col.upper()}', fontsize=13)
        ax_raw.tick_params(axis='x', labelsize=12)
        ax_raw.set_title(f'{self.name} Gene Expression (n={self.gene_table.shape[0]})', fontsize=14)

        # plot the log-transformed expression distribution
        sns.boxplot(x=self._exp_col_log, data=self.gene_table, fliersize=0, ax=ax_log)
        sns.stripplot(x=self._exp_col_log, data=self.gene_table, color='gray', ax=ax_log)
        ax_log.set_xlabel(f'log({self._exp_col.upper()})', fontsize=13)
        ax_log.tick_params(axis='x', labelsize=12)
        ax_log.set_title(
            f'{self.name} Gene Expression (Log-Transformed) (n={self.gene_table.shape[0]})',
            fontsize=14
        )

        return [ax_raw, ax_log]

    def within_tu_expression(self) -> List[plt.Axes]:
        """
        Calculate and plot the distribution of gene expressions within TUs. This distribution will
        be compared to a random, bootstrapped distribution, drawn from an empirical distribution
        based on the observed TU sizes. Only TUs with multiple genes, of course, are considered

        :return List[plt.Axes] axes: the Axes objects of the generated plot; the TU sizes and
            within vs random TU expression variation, respectively
        """

        tu_expression_dict = {}
        gene_counts_multi = []
        gene_counts_all = []
        for tu_row in self.tu_table.itertuples():
            tu_gene_ids = self._tu_to_genes[tu_row.Index]
            n_tu_genes = len(tu_gene_ids)
            # only consider this TU if it has multiple genes with expression data
            if n_tu_genes > 1:
                gene_expressions = self.gene_table.loc[tu_gene_ids, self._exp_col_log]
                tu_expression_dict[tu_row.Index] = gene_expressions
                gene_counts_multi.append(n_tu_genes)
            if n_tu_genes > 0:
                gene_counts_all.append(n_tu_genes)

        sd_within_tu = [np.std(exps) for exps in tu_expression_dict.values()]

        # bootstrap by random sample groups of genes and getting their standard deviations;
        # have to random sample the number of genes too; sample from the same distribution
        tu_size_kde = gaussian_kde(gene_counts_multi)

        sds_random = []
        for i in range(len(sd_within_tu)):
            rand_n_genes = int(max(2, np.round(tu_size_kde.resample(1).flatten()[0])))
            gene_idxes = np.random.choice(len(sd_within_tu), rand_n_genes)
            exps = self.gene_table.iloc[gene_idxes, :][self._exp_col_log]
            sds_random.append(np.std(exps))

        u_stat, p_val = mannwhitneyu(sd_within_tu, sds_random, alternative='two-sided')

        tu_n_genes_df = pd.DataFrame(data={'# of Genes': gene_counts_all})

        tu_exp_sd_df = pd.concat([
            pd.DataFrame(data={f'SD of {self._exp_col_log}': sd_within_tu, 'Group': 'Within TU'}),
            pd.DataFrame(data={f'SD of {self._exp_col_log}': sds_random, 'Group': 'Random'})
        ])

        _, (ax_tu, ax_sd) = plt.subplots(1, 2, figsize=(14, 5))

        # plot the distribution of TU counts
        sns.histplot(x='# of Genes', data=tu_n_genes_df, discrete=True, ax=ax_tu)
        ax_tu.tick_params(axis='both', labelsize=13)
        ax_tu.set_xlabel('# of Genes', fontsize=13)
        ax_tu.set_ylabel('# of TUs', fontsize=13)
        ax_tu.set_title(f'{self.name} TU Gene Count Distribution (n={len(gene_counts_all)})',
                        fontsize=14)

        # plot the SD distributions for within TU/random
        sns.boxplot(x='Group', y=f'SD of {self._exp_col_log}', data=tu_exp_sd_df, fliersize=0,
                    ax=ax_sd)
        sns.swarmplot(x='Group', y=f'SD of {self._exp_col_log}', data=tu_exp_sd_df, color='gray',
                      ax=ax_sd)
        ax_sd.tick_params(axis='both', labelsize=13)
        ax_sd.set_xlabel('')
        ax_sd.set_ylabel(f'SD of log({self._exp_col.upper()})', fontsize=13)
        ax_sd.set_title(f'{self.name} Within-TU Expression Variation\n'
                        f'P={p_val:.2E}', fontsize=13)
        ax_sd.set_xticklabels([
            f"{lab.get_text()}\n"
            f"(n={tu_exp_sd_df[tu_exp_sd_df['Group'] == lab.get_text()].shape[0]})"
            for lab in ax_sd.get_xticklabels()
        ])

        return [ax_tu, ax_sd]


# --- LOCAL UTILITIES ---
def _table_input_to_dataframe(table_input: Union[FILE_OR_DF, None],
                              table_label: str, fast_track: bool = False) -> pd.DataFrame:
    """
    Given a table input that's either None, a pandas DataFrame, or a file path to a DataFrame,
    return a DataFrame (empty or otherwise)

    :param FILE_OR_DF table_input: the input to parse into a DataFrame (if necessary)
    :param str table_label: the specific type of table this is (e.g. 'gene', 'protein', etc.)
    :param bool fast_track: skip error-checking
    :return pd.DataFrame dataframe: the DataFrame output based on the provided input table argument
    """

    required_columns = REQUIRED_COLUMNS[table_label]

    if table_input is None:
        parsed_table = pd.DataFrame(columns=required_columns)
    elif isinstance(table_input, pd.DataFrame):
        parsed_table = table_input
    elif isinstance(table_input, str) or isinstance(table_input, Path):
        table_input_path = Path(table_input)
        parsed_table = pd.read_csv(table_input_path, index_col=0)
    else:
        raise ValueError(f'Provided input for {table_label} table is neither a file path nor a'
                         f' DataFrame.')

    if not fast_track and not set(required_columns).issubset(parsed_table.columns):
        raise ValueError(f'Provided input for {table_label} table does not contain all of the'
                         f' required columns: {required_columns}')

    parsed_table_reindex = parsed_table.set_index('locus_tag')

    # for specific tables that don't need left/right (i.e. single point tables), add these columns
    if table_label in ['tss', 'tts']:
        parsed_table_reindex['left'] = parsed_table_reindex[table_label]
        parsed_table_reindex['right'] = parsed_table_reindex[table_label]

    return parsed_table_reindex
