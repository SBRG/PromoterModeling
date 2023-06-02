# built-in modules
from pathlib import Path
import re
from typing import List, Union

# third-party modules
import pandas as pd

# local modules
from bitome.features import Gene, TranscriptionFactor, TFBindingSite, IModulon

DATA_DIR = Path(Path(__file__).parent.parent, 'data')


def load_locus_tag_yome_lookup() -> dict:
    """
    Load y-ome information and return a dictionary with locus tags as keys and y-ome categorization as values

    :return dict y_ome_lookup:
    """
    y_ome_df = pd.read_csv(Path(DATA_DIR, 'y-ome', 'y-ome-genes.tsv'), sep='\t')
    y_ome_lookup = {}
    for y_ome_row in y_ome_df.itertuples(index=False):
        y_ome_lookup[y_ome_row.locus_id] = y_ome_row.category
    return y_ome_lookup


def load_ytfs_and_binding_sites(
            existing_tfs: List[TranscriptionFactor] = None
        ) -> List[Union[TranscriptionFactor, TFBindingSite]]:
    """
    Loads BitomeFeature objects for yTF binding sites observed experimentally in SBRG via ChIP-exo experiments

    :param List[TranscriptionFactor] existing_tfs: a list of pre-loaded TranscriptionFactor objects; if provided, and if
    a matching yTF is present, its binding sites will be added to this TF
    :return List[Union[TranscriptionFactor, TFBindingSite]] ytfs_and_binding_sites: a list of ytfs and binding sites
    """

    ytf_df = load_ytf_df()
    ytf_names = list(set(ytf_df['TF']))

    # ytfs_for_binding_sites includes all TFs, EVEN IF already existing in passed list; ytfs is the final return, NOT
    # including the pre-existing objects (all to avoid duplicates at higher level)
    ytfs_for_binding_sites = []
    ytfs = []
    for i, ytf_name in enumerate(ytf_names):

        if existing_tfs is not None:
            existing_tf_objs = [tf for tf in existing_tfs if tf.name == ytf_name]
            if existing_tf_objs:
                existing_tf = existing_tf_objs[0]
            else:
                existing_tf = None
        else:
            existing_tf = None

        if existing_tf is None:
            transcription_factor = TranscriptionFactor('yTF_' + str(i), ytf_name)
        else:
            transcription_factor = existing_tf

        ytf_site_df = ytf_df[ytf_df['TF'] == ytf_name]
        ytf_final_states = list(set(ytf_site_df['TF_FINAL_STATE']))

        for ytf_final_state in ytf_final_states:
            ytf_final_state_site_df = ytf_site_df[ytf_site_df['TF_FINAL_STATE'] == ytf_final_state]
            ytf_sites = []
            for ytf_site_row in ytf_final_state_site_df.itertuples(index=False):
                left_right = ytf_site_row.Start, ytf_site_row.End
                site_obj = TFBindingSite(
                    ytf_site_row.SITE_ID,
                    left_right
                )
                ytf_sites.append(site_obj)

            transcription_factor.add_binding_sites(ytf_sites, ytf_final_state)

        # don't append to the final ytf return list if we already had a TF object
        ytfs_for_binding_sites.append(transcription_factor)
        if existing_tf is None:
            ytfs.append(transcription_factor)

    ytf_binding_sites = []
    for ytf in ytfs_for_binding_sites:
        # TranscriptionFactor object stores binding sites in a dictionary, keys are conformation names, values are sites
        for conf_binding_sites in ytf.binding_sites.values():
            ytf_binding_sites += conf_binding_sites

    return ytfs + ytf_binding_sites


def load_ytf_df() -> pd.DataFrame:
    """
    Reads a file containing yTF binding sites and parses into a single dataframe

    :return pd.DataFrame ytf_df: a DataFrame listing binding sites for yTFs experimentally determined by ChIP-exo
    """

    # read_excel will return an OrderedDict with sheet_name: sheet_df items when sheet_name is None
    ytf_data_dict = pd.read_excel(Path(DATA_DIR, 'ytf_binding_sites.xlsx'), sheet_name=None)
    final_state_names = []
    tf_names = []
    dfs = []
    for ytf_name, ytf_df in ytf_data_dict.items():
        # capitalize the yTF name for cross-referencing later with RegulonDB; also swap underscore for dash to be
        # consistent with RegulonDB final conformation state nomenclature
        final_state_name = re.sub('_', '-', ytf_name[0].upper() + ytf_name[1:])
        tf_name = final_state_name[:4]
        final_state_names.append(final_state_name)
        tf_names.append(tf_name)
        dfs.append(ytf_df)
    for i, (final_name, base_tf_name) in enumerate(zip(final_state_names, tf_names)):
        dfs[i]['TF_FINAL_STATE'] = final_name
        dfs[i]['Transcription factors'] = base_tf_name

    full_ytf_df = pd.concat(dfs, ignore_index=True)
    full_ytf_df = full_ytf_df[['Transcription factors', 'Start', 'End', 'TF_FINAL_STATE']]
    full_ytf_df.rename(columns={'Transcription factors': 'TF'}, inplace=True)

    # add a pseudo ID column to identify each of these binding sites
    site_ids = ['ytf_site_' + str(i) for i in range(len(full_ytf_df))]
    full_ytf_df['SITE_ID'] = site_ids

    return full_ytf_df


def load_i_modulons(
            transcription_factors: List[TranscriptionFactor] = None,
            genes: List[Gene] = None
        ) -> List[IModulon]:
    """
    Loads i-Modulon objects, linking them to Gene objects they regulate (if gene objects provided)

    :param list[TranscriptionFactor] transcription_factors: TF objects already loaded to link to i-modulons
    :param list[Gene] genes: already loaded Gene objects to link to the i-modulons being loaded
    :return list[IModulon] i_modulons: i-modulon objects with information about the i-modulon and associated genes
    """

    i_modulon_metadata_df = pd.read_csv(Path(DATA_DIR, 'i_modulon.csv'), index_col=0)
    i_modulon_metadata_df = i_modulon_metadata_df[[
        'name',
        'Regulator'
    ]]
    # exclude rows that have nan in the Regulator column (not verified regulatory i-modulons)
    i_modulon_metadata_df = i_modulon_metadata_df.dropna(subset=['Regulator'])

    m_matrix_df = pd.read_csv(Path(DATA_DIR, 'i_modulon_m_matrix.csv'), index_col=0)

    i_modulons = []
    for i_modulon_row in i_modulon_metadata_df.itertuples():

        i_modulon_obj = IModulon(i_modulon_row.Index, i_modulon_row.name)

        if transcription_factors is not None:
            if isinstance(i_modulon_row.Regulator, str):
                i_modulon_tf_names = re.findall(r'[A-Z][a-z]{2}[A-Z]?', i_modulon_row.Regulator)
            else:
                i_modulon_tf_names = []
            i_modulon_tfs = []
            for transcription_factor in transcription_factors:
                if transcription_factor.name in i_modulon_tf_names:
                    i_modulon_tfs.append(transcription_factor)
            i_modulon_obj.link_transcription_factors(i_modulon_tfs)

        if genes is not None:
            i_modulon_gene_locus_tags = list(m_matrix_df.index[m_matrix_df[i_modulon_row.name]])
            i_modulon_genes = []
            for gene in genes:
                if gene.locus_tag in i_modulon_gene_locus_tags:
                    i_modulon_genes.append(gene)
            i_modulon_obj.link_genes(i_modulon_genes)

        i_modulons.append(i_modulon_obj)

    return i_modulons


def load_mrna_structure_ranges() -> List[tuple]:
    """
    Loads a list of tuples of mRNA folding energies as calculated by John Ooi/Troy Sandberg in SBRG (internal)
    Currently, only considers "tight" (top 10% of deltaG) mRNA folding across 100+ bp windows)

    :return List[tuple] mrna_structure_ranges: a list of start/end tuples for tight mRNA structures
    """

    raw_mrna_folding_df = pd.read_csv(Path(DATA_DIR, 'mrna_folding_energies.csv'))
    mrna_tight_df = raw_mrna_folding_df[raw_mrna_folding_df['Type'] == 'tight']

    # the locations are encoded weirdly; strings of two integers, separated by a dash, where the first number is the
    # start of the location range (100 bp window), and the last number is the last START of windows in this overall
    # mRNA region (so that last window + 100 is the last of the full range)
    location_strs = mrna_tight_df['Location']
    mrna_tight_structure_ranges = []
    for location_str in location_strs:
        start_end_list = [int(num.strip()) for num in location_str.split('-')]
        start = start_end_list[0]
        end = start_end_list[1] + 100
        mrna_tight_structure_ranges.append((start, end))

    return mrna_tight_structure_ranges


def load_term_seq() -> pd.DataFrame:
    """
    Loads a dataframe of 3' termination ends of transcription units

    :return pd.DataFrame term_seq_df: a DataFrame of termination sequence locations with strand information
    """
    term_seq_df = pd.read_csv(
        Path(DATA_DIR, 'Term_1629_single_tracked.gff'),
        sep='\t',
        header=None,
        names=['ref_genome', 'KNNC', 'TEP', 'TTS', 'right', 'strand', 'strand_alt', 'unknown', 'sample?']
    )
    # have validated that left/right are all the same
    term_seq_df = term_seq_df[['TTS', 'strand']].astype({'strand': int})
    return term_seq_df
