# built-in modules
from pathlib import Path
import re
from typing import List, Tuple, Union
from warnings import warn

# third-party modules
from Bio.SeqFeature import FeatureLocation
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import seaborn as sns

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']

DATA_DIR = Path(Path(__file__).parent.parent, 'data')


def is_int(potential_int):
    """
    Returns a boolean indicating if the given argument has type: int, np.int64

    :param potential_int: an input to check integer-ness for
    :return bool is_int: a boolean indicating if the input is a recognized integer type
    """
    return (
        isinstance(potential_int, int) or
        isinstance(potential_int, np.int64)
    )


def signed_int_str(input_int: int):
    """
    Converts an integer into a string that includes the integer's sign, + or -

    :param int input_int: the input integer to convert to a signed integer string
    :return str signed_int_str: the signed integer string
    """
    if input_int >= 0:
        return '+' + str(input_int)
    else:
        return str(input_int)


def get_reading_frame(location: FeatureLocation, num_positions: int) -> int:
    """
    Determine the reading frame of a genome location given the full length of the genome; position 1 in forward
    direction is considered reading frame +1, and the last position in the reverse direction is -1, etc.

    :param FeatureLocation location: a Biopython FeatureLocation object describing this feature's
    absolute start and end (0-indexed); strand may also be included. An alternative to manual providing of
    :param int num_positions: the total number of positions in the reference genome to which this location corresponds
    :return int reading_frame: the reading frame of the provided location
    """

    strand = location.strand
    if strand is None:
        raise ValueError('Strand must be known to determine reading frame')

    if strand == 1:
        start_for_reading_frame = location.start.position
    # assume strand is -1
    else:
        start_for_reading_frame = (num_positions - location.end.position)

    if start_for_reading_frame % 3 == 0:
        reading_frame = 3
    else:
        reading_frame = start_for_reading_frame % 3

    if strand == -1:
        reading_frame = reading_frame * -1

    return reading_frame


def to_strand_int(strand_str: str) -> int:
    """
    Converts a strand in string form into the integers +1 and -1 for forward and reverse strands, respectively
    :param str strand_str: a string description of a strand
    :return int strand_int: the integer version of the provided strand string
    """
    if strand_str == 'forward':
        strand_int = +1
    elif strand_str == 'reverse':
        strand_int = -1
    else:
        warn('Unrecognized strand string; defaulting to +1 (forward strand)')
        strand_int = 1

    return strand_int


def select_features_by_type(features: list, feature_type: str) -> list:
    """
    Selects features from a list of BitomeFeatures that have the provided type as their .type attribute
    :param List[BitomeFeature] features: the list of BitomeFeature objects (or sub-objects thereof) from which to
    extract features of a given type
    :param str feature_type: the type of feature to select
    :return list selected_features: the features selected from the list by type
    """
    return [feature for feature in features if feature.type == feature_type]


def parse_gempro(locus_tag: str) -> Union[pd.DataFrame, None]:
    """
    Given a locus tag, parses the appropriate gempro local amino acid information file, returning a DataFrame with
    amino acids of the protein as the index

    :param str locus_tag: the locus tag used to identify this protein on the genome
    :return pd.DataFrame local_gempro_df: the local GEM-PRO information about this protein, as a pandas DataFrame
    """
    # make sure the locus tag has a local gempro file
    full_file_path = Path(DATA_DIR, 'local-gempro', locus_tag + '.csv')
    if full_file_path.exists():
        return pd.read_csv(full_file_path, index_col=0)
    else:
        return None


def load_locus_tag_to_gene_names_lookup() -> dict:
    """
    Parses a raw text file with a lookup between E. coli gene b-numbers and all gene name synonyms into a dictionary

    :return dict tag_name_lookup: a dictionary serving as a mapping between locus tags and all synonym gene names
    """

    with open(Path(DATA_DIR, 'e_coli_gene_ids.txt'), 'r') as f:
        lines = f.readlines()

    header = True
    i = 0
    while header:
        # find first line that has a b-number as the first entry
        if re.search(r'^b\d{4}', lines[i]):
            header = False
        else:
            i += 1

    end_of_b_numbers = True
    j = len(lines) - 1
    while end_of_b_numbers:
        # find first line that has a b-number as the first entry
        if re.search(r'^b\d{4}', lines[j]):
            end_of_b_numbers = False
        else:
            j -= 1

    # order of the gene names matters; the first is considered the "true" name (from UniProt)
    tag_name_lookup = {}
    for line in lines[i:j+1]:
        b_number = re.findall(r'^b\d{4}', line)[0]
        gene_names = re.findall(r'[a-z]{3}[A-Z]?[A-z0-9]?', line)
        tag_name_lookup[b_number] = gene_names

    return tag_name_lookup


def load_locus_tag_cogs_lookup() -> dict:
    """
    Parses a raw file that contains associations between gene locus tags (b-numbers) and COG labels (single letters)

    :return dict locus_tag_cog_lookup: a dictionary with keys that are gene locus tags and values with the gene's COG
    """
    cogs_df = pd.read_csv(Path(DATA_DIR, 'cogs_ecoli_mg1655.csv'))
    locus_tag_cog_lookup = {locus_tag: cog for locus_tag, cog in zip(cogs_df['locus'], cogs_df['COG category primary'])}
    return locus_tag_cog_lookup


def load_essential_genes() -> list:
    """
    Parses a raw file containing annotation of gene essentiality from the Keio collection and 2 previous studies
    NOTE: only genes that are considered essential by all 3 studies (a score of 3 in the file) are noted as essential

    :return list essential_genes: a list of the locus tags of all essential genes
    """
    return list(pd.read_csv(Path(DATA_DIR, 'keio_essentiality_parsed.csv'))['b_number'])


def bits_per_bp_plot(
            bitome,
            ranges_by_feature: List[List[Tuple[int, int]]],
            names_by_feature: List[str],
            kde: bool = True,
            figsize: tuple = None,
            median: bool = False,
            split: bool = False,
            compare: bool = False,
            file_path: Union[str, Path] = None,
            return_axs: bool = False,
            show: bool = True
        ) -> Union[None, List[plt.Axes]]:
    """
    Create a plot showing the bits per bp (information density) of feature(s) for a bitome object

    :param bitome: a Bitome object from which to extract the feature ranges
    :param List[Tuple[int, int]] ranges_by_feature: list(s) of lists of tuples, where the lists of tuples correspond
    to genomic position ranges for the corresponding feature name in names_by_feature
    :param List[str] names_by_feature: the name(s) of the feature(s) whose ranges are passed in ranges_by_feature;
    must be index-matched to the former
    :param bool kde: indicates if a KDE should be included in the histogram
    :param tuple figsize: figsize to be passed to plt.subplots() argument of the same name
    :param bool median: indicate if the median of the distribution should be highlighted by a vertical line per subplot
    :param bool split: indicates if the features should be split into "small" and "large" groups around the median
    length, and each group's information density distribution plotted separately
    :param bool compare: indicates, if split is True, that the small/large groups' information densities should be
    statistically compared
    :param Union[Path, str] file_path: if set, a path to which the produced figure should be saved
    :param bool return_axs: indicates if any Axes generated should be returned
    :param bool show: indicates if the plt.show() command should be issued within this function
    :return Union[None, plt.Axes] none_or_axes: either None, or, if return_axs is True, the plt.Axes objects generated
    """

    bits_per_bp_all_features: List[Union[
        List[float],
        Tuple[List[float], List[float]]
    ]] = []
    for feature_ranges in ranges_by_feature:

        lengths_single_feature: List[float] = []
        bits_per_bp_single_feature: List[float] = []
        for feature_range in feature_ranges:
            feature_sub_matrix = bitome.extract(column_range=feature_range)
            feature_length = feature_sub_matrix.shape[1]
            lengths_single_feature.append(feature_length)
            bits_per_bp_single_feature.append(feature_sub_matrix.sum()/feature_length)

        if split:
            median_length = np.median(lengths_single_feature)
            small_half_indices = np.where(lengths_single_feature <= median_length)[0]
            large_half_indices = np.where(lengths_single_feature > median_length)[0]
            bits_per_bp_single_feature_small = np.array(bits_per_bp_single_feature)[small_half_indices]
            bits_per_bp_single_feature_large = np.array(bits_per_bp_single_feature)[large_half_indices]
            bits_per_bp_all_features.append((bits_per_bp_single_feature_small, bits_per_bp_single_feature_large))
        else:
            bits_per_bp_all_features.append(bits_per_bp_single_feature)

    if figsize is None:
        figsize = (4*len(names_by_feature), 8)
    fig, axs = plt.subplots(len(names_by_feature), 1, figsize=figsize, sharex='all')

    for bits_per_bp_list_or_tuple, name, ax in zip(bits_per_bp_all_features, names_by_feature, axs):
        title = f'{name} ({len(bits_per_bp_list_or_tuple)})'
        if split:
            bits_per_bp_small, bits_per_bp_large = bits_per_bp_list_or_tuple
            sns.distplot(
                bits_per_bp_small,
                ax=ax,
                kde=kde,
                color=(0, 0, 1, 0.4),
                kde_kws={'linewidth': 5},
                bins=np.arange(0, 21, 0.5)
            )
            sns.distplot(
                bits_per_bp_large,
                ax=ax,
                kde=kde,
                color=(1, 0, 0, 0.4),
                kde_kws={'linewidth': 5},
                bins=np.arange(0, 21, 0.5)
            )

            if compare:
                mann_whitney_u, p_value = mannwhitneyu(bits_per_bp_small, bits_per_bp_large, use_continuity=True)
                title = f'{title} (P: {p_value:.2f})'
        else:
            sns.distplot(
                bits_per_bp_list_or_tuple,
                ax=ax,
                kde=kde,
                kde_kws={'linewidth': 3},
                bins=np.arange(0, 21, 0.5)
            )
            if median:
                ax.axvline(np.median(bits_per_bp_list_or_tuple))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(axis='x', reset=True)
        ax.set_xticks([5, 10, 15, 20])
        ax.xaxis.set_ticks_position('bottom')
        ax.text(0.97, 0.15, title, fontsize=28, transform=ax.transAxes, horizontalalignment='right')
        ax.tick_params(axis='both', labelsize='25')

    if split:
        legend_elems = [
            Patch(facecolor=(0, 0, 1, 0.4), edgecolor=(0, 0, 1, 0.4), label='small'),
            Patch(facecolor=(1, 0, 0, 0.4), edgecolor=(1, 0, 0, 0.4), label='large')
        ]
        fig.legend(handles=legend_elems, prop={'size': 14}, loc='upper right')

    axs[-1].set_xlabel('Information density (bits/bp)', fontsize=34)
    plt.tight_layout()

    if file_path is not None:
        plt.savefig(f'{file_path}.svg')
    if show:
        plt.show()
    if return_axs:
        return axs
    else:
        return None
