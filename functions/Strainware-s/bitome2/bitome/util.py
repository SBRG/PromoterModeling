# built-in modules
from typing import Dict, List, Tuple, Union

# third-party modules
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import pandas as pd
from scipy.stats import mode


def genbank_to_feature_tables(genbank_record: SeqRecord) -> \
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parse a GenBank record (pre-parsed into Biopython's SeqRecord format) into (up to) 3 feature
    tables (pandas DataFrames). The first table will list all genes, the second table will list
    proteins (CDS in GenBank nomenclature), and the third will list any miscellaneous additional
    features present in the record (e.g. mobile elements, repeat regions, origin of replication)

    :param SeqRecord genbank_record: a Biopython-parsed GenBank record containing genomic features
    :return Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame] genbank_feature_tables: 3 distinct
        pandas DataFrames containing tables parsed from the GenBank file, in the following order:
        gene_table, protein_table, misc_feature_table. Any/all of these may be completely empty,
        depending on the contents of the particular GenBank record
    """

    gene_rows = []
    protein_rows = []
    misc_feature_rows = []

    unique_id_counter = 1

    for feature in genbank_record.features:

        # feature types to skip: 'source' (ref sequence repeated), 'STS' (sequence-tagged site,
        # relevant only for biochemistry)
        if feature.type in ['source', 'STS']:
            continue

        feature_row_dict = {}

        left, right = feature.location.start.position, feature.location.end.position
        # Biopython's feature locations are 0-indexed; to convert to 1-indexed, add 1 to left
        feature_row_dict['left'], feature_row_dict['right'] = left + 1, right
        feature_row_dict['strand'] = feature.location.strand

        locus_tag = _get_seqfeature_qualifier(feature, 'locus_tag')
        if locus_tag is None:
            locus_tag = f'{genbank_record.id}_{unique_id_counter}'
            unique_id_counter += 1
        feature_row_dict['locus_tag'] = locus_tag

        if feature.type == 'gene':
            # TODO: there may be other ways that GenBank records specify a pseudogene
            feature_row_dict['pseudo'] = 'pseudo' in feature.qualifiers
            feature_row_dict['name'] = _get_seqfeature_qualifier(feature, 'gene')
            gene_rows.append(feature_row_dict)
        elif feature.type == 'CDS':
            feature_row_dict['name'] = _get_seqfeature_qualifier(feature, 'product')
            protein_rows.append(feature_row_dict)
        else:
            feature_row_dict['type'] = feature.type
            misc_feature_rows.append(feature_row_dict)

    # make sure to set the unique locus tags as the indices; ASSUME these are unique
    return (pd.DataFrame(gene_rows).set_index('locus_tag'),
            pd.DataFrame(protein_rows).set_index('locus_tag'),
            pd.DataFrame(misc_feature_rows).set_index('locus_tag'))


def _get_seqfeature_qualifier(seqfeature: SeqFeature, key: str) -> Union[str, None]:
    """
    Get a non-null attribute from a Biopython SeqFeature object

    :param SeqFeature seqfeature: the Biopython feature object from which to get a qualifier value
    :param str key:a key for the qualifiers attribute, which is a dictionary
    :return Union[str] value: the value stored in the provided feature's qualifiers dictionary for
        the given key; returns None if the key is missing or has an empty value
    """

    try:
        value = seqfeature.qualifiers[key]
    except KeyError:
        return None

    non_empty_str = value[0].strip()
    if non_empty_str == '':
        return None
    else:
        return non_empty_str


def genome_point_to_point(point1: float, point2: float, genome_length: int) -> float:
    """
    Compute the 1-D distance between 2 points on the genome, taking into account the genome's
    circularity

    :param float point1: the first genome position for which to compute an inter-point distance
    :param float point2: the second genome position for which to compute an inter-point distance
    :param int genome_length: the total length (linear) of the genome, a.k.a. the max position index
    :return float distance: the distance between the two points, taking into account genome
        circularity
    """

    within_distance = abs(point1 - point2)
    wraparound_distance = genome_length - max(point1, point2) + min(point1, point2)
    return min(within_distance, wraparound_distance)


def location_to_point(left: int, right: int, strand: Union[int, None], reference: str) -> float:
    """
    Return a single reference point for a feature based on its location and the desired reference

    :param int left: the left end of the feature
    :param int right: the right end of the feature
    :param int strand: the strand of the feature (may be None)
    :param str reference: the reference position to compute: 'start', 'end', or 'midpoint';
        'midpoint' will ALWAYS be used if strand is unknown
    :return float point: the point used to summarize the feature's location
    """

    if strand is None or reference == 'midpoint':
        feature_point = (right + left) / 2
    elif strand == 1:
        if reference == 'start':
            feature_point = left
        else:
            feature_point = right
    elif strand == -1:
        if reference == 'start':
            feature_point = right
        else:
            feature_point = left
    else:
        raise ValueError(f'Provided strand value {strand} is neither 1 nor -1')

    return feature_point


def score_motif_match(sequence: Seq, motif_pssm: Dict[str, List[float]]) -> float:
    """
    Given a candidate sequence and a position-specific scoring matrix (PSSM) in log-odds form,
    determine a log-odds score for the sequence's match to the consensus motif

    :param Seq sequence: the sequence to score against the given motif
    :param Dict[str, List[float]] motif_pssm: the PWM for the motif used to score the sequence
    :return float log_odds: the log-odds probability score for the sequence's match to the motif
    """

    if len(sequence) != len(motif_pssm['A']):
        raise ValueError(f"Sequence ({len(sequence)}) and PSSM ({len(motif_pssm['A'])}) do not have"
                         f" same lengths.")

    log_odds = 0
    for i, base in enumerate(sequence):
        # handles the case where a non-ACGT base pair is in the genome; just ignores it
        # TODO: score non-ACGT base pair codes
        if base in motif_pssm:
            log_odds += motif_pssm[base][i]

    return log_odds


def create_motif(sequences: List[Seq]) -> motifs.Motif:
    """
    Given a list of sequences, create a Biopython style motif object; automatically picks out the
    majority-length sequences if different lengths are given

    :param List[Seq] sequences: the sequences from which to generate a motif
    :return motifs.Motif motif: the Biopython Motif object generated from the provided sequences
    """
    # motifs have to be made from sequences of the same length for now
    sequence_lens = [len(seq) for seq in sequences]
    mode_len, _ = mode(sequence_lens)
    sequences = [seq for seq in sequences if len(seq) == mode_len]
    return motifs.create(sequences)


def one_hot_consensus(sequences: List[Seq], position_indices=None) -> pd.Series:
    """
    Given a list of aligned sequences, return a one-hot encoding

    :param List[Seq] sequences: the aligned sequences to one-hot encode
    :param position_indices: optional, alternative indices to use when creating the one hot row
    """

    base_counts = create_motif(sequences).counts

    if position_indices is None:
        position_indices = list(range(1, len(sequences[0] + 1)))

    base_rows = []
    for base in 'ATCG':
        base_row = pd.Series(base_counts[base], index=[f'{pos}_{base}' for pos in position_indices])
        base_rows.append(base_row)
    full_row = pd.concat(base_rows).divide(len(sequences))

    return full_row
