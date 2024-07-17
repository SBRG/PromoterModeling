# built-in modules
from math import isnan
from pathlib import Path
import re
from typing import List, Union

# third-party modules
from Bio.Seq import Seq
import pandas as pd

# local modules
from bitome.features import Gene, Operon, TranscriptionUnit, Promoter, ShineDalgarno, Terminator, Attenuator,\
                            Regulon, TranscriptionFactor, TFBindingSite, Riboswitch
from bitome.sbrg import load_term_seq
from bitome.utilities import load_locus_tag_to_gene_names_lookup, to_strand_int


REGULON_DB_FEATURE = Union[
    Gene, Operon, TranscriptionUnit, Promoter, ShineDalgarno, Terminator, Attenuator,
    Regulon, TranscriptionFactor, TFBindingSite, Riboswitch
]

DATA_DIR = Path(Path(__file__).parent.parent, 'data')


def load_regulon_db_features(reference_sequence: Seq, genbank_genes: List[Gene]) -> List[REGULON_DB_FEATURE]:
    """
    A master function that calls sub-functions to load genomic features stored in RegulonDB, and then links them
    together appropriately

    :param Bio.Seq.Seq reference_sequence: the reference GenBank sequence, as a Biopython Seq object
    :param List[Gene] genbank_genes: a list of pre-loaded GenBank gene objects that can be further linked to RegDB objs
    :return List[REGULON_DB_FEATURE] all_regulon_db_features: a list of all RegulonDB features loaded from the files
    """

    transcription_units = load_transcription_units(reference_sequence)
    operons = load_operons(reference_sequence)
    promoters = load_promoters(reference_sequence)

    # to the best of our ability, add RegulonDB IDs to the genbank genes
    genes = add_regulon_db_gene_ids(genbank_genes)

    # link transcription units to their promoters and operons by using the IDs stored in the TU objects
    tu_gene_link_df = load_regulon_db_df('tu_gene_link')
    for transcription_unit in transcription_units:

        # based on manual parsing of data before writing this code, every TU has either 0 or 1 promoter
        # NOTE: the reading_frame of the promoter, based on TSS, is copied to the TU
        tu_promoter_id = transcription_unit.promoter_id
        tu_promoter = select_features_by_id(promoters, tu_promoter_id, return_first=True)
        transcription_unit.link_promoter(tu_promoter)

        # based on manual parsing of data before writing this code, every TU has exactly 1 operon
        # reading frame of TU is copied to operon (important that it's taken from promoter TSS, above section)
        tu_operon_id = transcription_unit.operon_id
        tu_operon = select_features_by_id(operons, tu_operon_id, return_first=True)
        transcription_unit.link_operon(tu_operon)

        tu_id = transcription_unit.id
        tu_gene_ids = list(tu_gene_link_df[tu_gene_link_df['TRANSCRIPTION_UNIT_ID'] == tu_id]['GENE_ID'])
        tu_genes = select_features_by_id(genes, tu_gene_ids)
        transcription_unit.link_genes(tu_genes)

    # the following features are much sparser than their "parent" feature, so it's faster to loop over the small lists
    terminators = load_terminators()
    for terminator in terminators:
        terminator_tu = select_features_by_id(transcription_units, terminator.tu_id, return_first=True)
        terminator.link_transcription_unit(terminator_tu, reference_sequence=reference_sequence)

    # for attenuators, riboswitches, and terminators; since we are using the Genbank genes, we may not have loaded
    # the RegulonDB ID for genes that are only in RegulonDB; if these feature are linked to such genes, they won't
    # show up; thus we need to check that the linked gene is None and only call the link_gene method if not
    attenuators = load_attenuators(reference_sequence)
    for attenuator in attenuators:
        attenuator_gene = select_features_by_id(genes, attenuator.gene_id, return_first=True)
        if isinstance(attenuator_gene, Gene):
            attenuator.link_gene(attenuator_gene)

    riboswitches = load_riboswitches()
    for riboswitch in riboswitches:
        riboswitch_gene = select_features_by_id(genes, riboswitch.gene_id, return_first=True)
        if isinstance(riboswitch_gene, Gene):
            riboswitch.link_gene(riboswitch_gene, reference_sequence=reference_sequence)

    shine_dalgarnos = load_shine_dalgarnos()
    for shine_dalgarno in shine_dalgarnos:
        shine_dalgarno_gene = select_features_by_id(genes, shine_dalgarno.gene_id, return_first=True)
        if isinstance(shine_dalgarno_gene, Gene):
            shine_dalgarno.link_gene(shine_dalgarno_gene, reference_sequence=reference_sequence)

    # pull out all the binding sites from the transcription factors so we can link them to promoters
    transcription_factors = load_transcription_factors()
    tf_binding_sites = []
    for transcription_factor in transcription_factors:
        # TranscriptionFactor object stores binding sites in a dictionary, keys are conformation names, values are sites
        for conf_binding_sites in transcription_factor.binding_sites.values():
            tf_binding_sites += conf_binding_sites

    site_promoter_link_df = load_regulon_db_df('tf_site_promoter_link')
    for promoter in promoters:

        binding_site_ids = list(site_promoter_link_df[site_promoter_link_df['PROMOTER_ID'] == promoter.id]['SITE_ID'])
        binding_sites = select_features_by_id(tf_binding_sites, binding_site_ids)
        promoter.link_binding_sites(binding_sites)

    regulons = load_regulons()
    regulon_tf_link_df = load_regulon_db_df('regulon_tf_link')
    regulon_function_df = load_regulon_db_df('regulon_function')
    reg_func_prom_link_df = load_regulon_db_df('regulon_func_prom_link')
    for regulon in regulons:

        regulon_tf_ids = regulon_tf_link_df[regulon_tf_link_df['REGULON_ID'] == regulon.id]['TRANSCRIPTION_FACTOR_ID']
        regulon_tfs = select_features_by_id(transcription_factors, list(regulon_tf_ids))
        regulon.link_transcription_factors(regulon_tfs)

        # find the known "functions" of the regulon; mostly just 1, but some have two; these are regulatory modes
        this_regulon_function_df = regulon_function_df[regulon_function_df['REGULON_ID'] == regulon.id]
        for regulon_function_row in this_regulon_function_df.itertuples(index=False):

            reg_func_promoter_ids = reg_func_prom_link_df[
                reg_func_prom_link_df['REGULON_FUNCTION_ID'] == regulon_function_row.REGULON_FUNCTION_ID
            ]['PROMOTER_ID']
            reg_func_promoters = select_features_by_id(promoters, list(reg_func_promoter_ids))
            regulon_function_name = regulon_function_row.REGULON_FUNCTION_NAME
            regulon_function_name = re.sub(r'\s+', ' ', regulon_function_name.strip())

            regulon.link_promoters(reg_func_promoters, regulon_function_name)

    # assemble a massive list of all features to be returned; they can be selected at will by upper-level functions
    # based on the type attribute of each feature
    all_regulon_db_features = sum(
        [
            transcription_units, operons, promoters, terminators, attenuators,
            riboswitches, shine_dalgarnos, transcription_factors, tf_binding_sites, regulons
        ],
        []
    )

    return all_regulon_db_features


def select_features_by_id(
            features: List[REGULON_DB_FEATURE],
            feature_ids: Union[List[str], str],
            return_first: bool = False
        ) -> Union[List[REGULON_DB_FEATURE], REGULON_DB_FEATURE, None]:
    """
    Return the feature(s) from the provided feature list that have the given RegulonDB ID

    :param List[REGULON_DB_FEATURE] features: a list of BitomeFeature objects; must have an id attribute
    :param Union[List[str], str] feature_ids: the RegulonDB ID(s) to look for
    :param bool return_first: indicates if only the first instance of a feature with the given ID should be returned
    :return Union[List[REGULON_DB_FEATURE], REGULON_DB_FEATURE, None] found_features: the selected feature(s) with the\
    provided ID; returns None/empty list if none found (depending on value of return_first parameter)
    """

    if isinstance(feature_ids, list):
        feature_ids_list = feature_ids
    else:
        feature_ids_list = [feature_ids]

    found_features = []
    for feature in features:
        if feature.id in feature_ids_list:
            found_features.append(feature)
            if return_first:
                break

    if return_first:
        if found_features:
            return found_features[0]
        else:
            return None
    else:
        return found_features


def load_operons(reference_sequence: Seq) -> List[Operon]:
    """
    Load operon features from the operon CSV parsed from RegulonDB

    :param Seq reference_sequence: the reference GenBank sequence, as a Biopython Seq object
    :return List[Operon] operons: a list of Operon objects
    """

    operon_df = load_regulon_db_df('operon')

    operons = []
    for operon_row in operon_df.itertuples(index=False):

        left_right = operon_row.REGULATIONPOSLEFT, operon_row.REGULATIONPOSRIGHT
        strand = to_strand_int(operon_row.OPERON_STRAND)

        operon_obj = Operon(
            operon_row.OPERON_ID,
            left_right,
            strand,
            reference_sequence,
            name=operon_row.OPERON_NAME
        )
        operons.append(operon_obj)

    return operons


def load_transcription_units(reference_sequence: Seq) -> List[TranscriptionUnit]:
    """
    Load transcription unit features from the transcription unit CSV parsed from RegulonDB;
    this function does NOT take a reference sequence since we don't directly have strand information for the TUs
    (need to link to operons first; that is handled in load_regulon_db_features)

    :param Seq reference_sequence: the reference GenBank sequence, as a Biopython Seq object
    :return List[TranscriptionUnit] transcription_units: a list of TranscriptionUnit objects
    """

    tu_df = load_regulon_db_df('transcription_unit')
    term_seq_df = load_term_seq()

    tus = []
    for tu_row in tu_df.itertuples(index=False):

        left_right = tu_row.TU_POSLEFT, tu_row.TU_POSRIGHT
        tss_raw = tu_row.POS_1
        if isnan(tss_raw):
            tss = None
        else:
            tss = int(tss_raw)
        strand = to_strand_int(tu_row.TU_STRAND)

        # look for the nearest TTS; have parsed out the longest legit 3'UTR in RegulonDB is 726 bp long; thus, set
        # a cap at this length for Term-Seq data association to RegulonDB
        if strand == 1:
            term_seq_f_df = term_seq_df[term_seq_df['strand'] == 1]
            term_seq_right_df = term_seq_f_df[left_right[1] < term_seq_f_df['TTS']]
            tu_f_df = tu_df[tu_df['TU_STRAND'] == 'forward']
            tu_right_starts = list(tu_f_df[left_right[1] < tu_f_df['TU_POSLEFT']]['TU_POSLEFT'])
            closest_tts = min(list(term_seq_right_df['TTS']), default=left_right[1])
            next_tu_start = min(tu_right_starts, default=closest_tts+1)
            if closest_tts > next_tu_start or closest_tts > left_right[1] + 726:
                closest_tts = left_right[1]

        else:
            term_seq_r_df = term_seq_df[term_seq_df['strand'] == 1]
            term_seq_left_df = term_seq_r_df[term_seq_r_df['TTS'] < left_right[0]]
            tu_r_df = tu_df[tu_df['TU_STRAND'] == 'reverse']
            tu_left_starts = list(tu_r_df[tu_r_df['TU_POSRIGHT'] < left_right[0]]['TU_POSRIGHT'])
            closest_tts = max(list(term_seq_left_df['TTS']), default=left_right[0])
            next_tu_start = max(tu_left_starts, default=closest_tts-1)
            if closest_tts < next_tu_start or closest_tts < left_right[0]-726:
                closest_tts = left_right[0]

        # not guaranteed to have a promoter ID
        promoter_id_raw = tu_row.PROMOTER_ID
        if isinstance(promoter_id_raw, str):
            promoter_id = promoter_id_raw
        else:
            promoter_id = None

        tu_obj = TranscriptionUnit(
            tu_row.TRANSCRIPTION_UNIT_ID,
            tu_row.OPERON_ID,
            left_right,
            strand,
            reference_sequence,
            tss=tss,
            tts=closest_tts,
            promoter_id=promoter_id,
            name=tu_row.TRANSCRIPTION_UNIT_NAME
        )
        tus.append(tu_obj)

    return tus


def load_promoters(reference_sequence: Seq) -> List[Promoter]:
    """
    Load promoter features from the promoter CSV parsed from RegulonDB

    :param Seq reference_sequence: the reference GenBank sequence, as a Biopython Seq object
    :return List[Promoter] promoters: a list of Promoter objects
    """

    promoter_df = load_regulon_db_df('promoter')

    promoters = []
    for promoter_row in promoter_df.itertuples(index=False):

        tss = int(promoter_row.POS_1)
        strand = to_strand_int(promoter_row.PROMOTER_STRAND)

        # RegulonDB definitions for promoter range
        if strand == +1:
            left_right = tss - 60, tss + 20
        else:
            left_right = tss - 20, tss + 60

        sigma_factors_raw = promoter_row.SIGMA_FACTOR
        if isinstance(sigma_factors_raw, str):
            sigma_factors = [sig.strip() for sig in sigma_factors_raw.split(',')]
        else:
            sigma_factors = []

        # box regions are always together, so just check one to see if we have all
        if isnan(promoter_row.BOX_10_LEFT):
            box_10_left_right = None
            box_35_left_right = None
        else:
            # the DataFrame stores these columns as floats since there are some nan; int columns can't have nan
            box_10_left_right = int(promoter_row.BOX_10_LEFT), int(promoter_row.BOX_10_RIGHT)
            box_35_left_right = int(promoter_row.BOX_35_LEFT), int(promoter_row.BOX_35_RIGHT)

        promoter_obj = Promoter(
            promoter_row.PROMOTER_ID,
            tss,
            left_right,
            strand,
            reference_sequence,
            box_10_left_right=box_10_left_right,
            box_35_left_right=box_35_left_right,
            sigma_factors=sigma_factors,
            name=promoter_row.PROMOTER_NAME
        )

        promoters.append(promoter_obj)

    return promoters


def load_terminators() -> List[Terminator]:
    """
    Load terminator features from the terminator CSV parsed from RegulonDB
    this function does NOT take a reference sequence since we don't directly have strand information;
    (need to link to operons first; that is handled in load_regulon_db_features)

    :return List[Terminator] terminators: a list of Terminator objects
    """

    terminator_df = load_regulon_db_df('terminator')

    terminators = []
    for terminator_row in terminator_df.itertuples(index=False):

        left_right = terminator_row.TERMINATOR_POSLEFT, terminator_row.TERMINATOR_POSRIGHT
        terminator_obj = Terminator(
            terminator_row.TERMINATOR_ID,
            left_right,
            terminator_row.TERMINATOR_CLASS,
            terminator_row.TRANSCRIPTION_UNIT_ID
        )
        terminators.append(terminator_obj)

    return terminators


def load_attenuators(reference_sequence: Seq) -> List[Attenuator]:
    """
    Load attenuator features from the attenuator CSV parsed from RegulonDB

    :param Seq reference_sequence: the reference GenBank sequence, as a Biopython Seq object
    :return List[Attenuator] attenuators: a list of Attenuator objects
    """

    attenuator_df = load_regulon_db_df('attenuator')

    attenuators = []
    for attenuator_row in attenuator_df.itertuples(index=False):

        left_right = attenuator_row.A_TERMINATOR_POSLEFT, attenuator_row.A_TERMINATOR_POSRIGHT
        strand = to_strand_int(attenuator_row.ATTENUATOR_STRAND)

        attenuator_obj = Attenuator(
            attenuator_row.ATTENUATOR_ID,
            left_right,
            strand,
            reference_sequence,
            attenuator_row.ATTENUATOR_TYPE,
            attenuator_row.A_TERMINATOR_TYPE,
            attenuator_row.GENE_ID
        )
        attenuators.append(attenuator_obj)

    return attenuators


def load_shine_dalgarnos() -> List[ShineDalgarno]:
    """
    Load Shine-Dalgarno features from the Shine-Dalgarno CSV parsed from RegulonDB
    this function does NOT take a reference sequence since we don't directly have strand information;
    (need to link to operons first; that is handled in load_regulon_db_features)

    :return List[ShineDalgarno] shine_dalgarnos: a list of ShineDalgarno objects
    """

    shine_dalgarno_df = load_regulon_db_df('shine_dalgarno')

    shine_dalgarnos = []
    for shine_dalgarno_row in shine_dalgarno_df.itertuples(index=False):

        left_right = shine_dalgarno_row.SHINE_DALGARNO_POSLEFT, shine_dalgarno_row.SHINE_DALGARNO_POSLEFT
        shine_dalgarno_obj = ShineDalgarno(
            shine_dalgarno_row.SHINE_DALGARNO_ID,
            left_right,
            shine_dalgarno_row.GENE_ID
        )
        shine_dalgarnos.append(shine_dalgarno_obj)

    return shine_dalgarnos


def load_riboswitches() -> List[Riboswitch]:
    """
    Load RFAM features from the RFAM CSV parsed from RegulonDB; these are riboswitches
    this function does NOT take a reference sequence since we don't directly have strand information;
    (need to link to operons first; that is handled in load_regulon_db_features)

    :return List[Riboswitch] riboswitches: a list of Riboswitch objects
    """

    rfam_df = load_regulon_db_df('rfam')

    riboswitches = []
    for rfam_row in rfam_df.itertuples(index=False):

        left_right = rfam_row.RFAM_POSLEFT, rfam_row.RFAM_POSRIGHT
        riboswitch_obj = Riboswitch(rfam_row.RFAM_ID, left_right, rfam_row.RFAM_DESCRIPTION, rfam_row.GENE_ID)
        riboswitches.append(riboswitch_obj)

    return riboswitches


def load_transcription_factors() -> List[TranscriptionFactor]:
    """
    Loads transcription factors and their associated binding sites from RegulonDB files; returns just a list
    of transcription factor objects that have binding sites associated with them

    :return List[TranscriptionFactor] transcription_factors: a list of transcription factor objects that are themselves
    not associated with a particular genomic location but have a binding_sites attribute containing observed binding
    sites
    """

    transcription_factor_df = load_regulon_db_df('transcription_factor')
    conformation_df = load_regulon_db_df('tf_conformation')
    site_df = load_regulon_db_df('tf_binding_site')

    transcription_factors = []
    for tf_row in transcription_factor_df.itertuples(index=False):

        tf_id = tf_row.TRANSCRIPTION_FACTOR_ID

        if isnan(tf_row.SITE_LENGTH):
            site_length = None
        else:
            site_length = tf_row.SITE_LENGTH

        tf_obj = TranscriptionFactor(
            tf_row.TRANSCRIPTION_FACTOR_ID,
            tf_row.TRANSCRIPTION_FACTOR_NAME,
            site_length=site_length
        )

        # use the conformation-TF dataframe to get the conformation row that correspond to this TF
        tf_conformations_df = conformation_df[conformation_df['TRANSCRIPTION_FACTOR_ID'] == tf_id]

        for conformation_row in tf_conformations_df.itertuples(index=False):

            # get the binding sites for this conformation
            conformation_id = conformation_row.CONFORMATION_ID
            conf_sites_df = site_df[site_df['CONFORMATION_ID'] == conformation_id]

            sites = []
            for site_row in conf_sites_df.itertuples(index=False):

                left_right = site_row.SITE_POSLEFT, site_row.SITE_POSRIGHT
                site_obj = TFBindingSite(
                    site_row.SITE_ID,
                    left_right,
                    conformation_id=conformation_id
                )
                sites.append(site_obj)

            tf_obj.add_binding_sites(sites, conformation_row.FINAL_STATE)

        transcription_factors.append(tf_obj)

    return transcription_factors


def add_regulon_db_gene_ids(genbank_genes: List[Gene]) -> List[Gene]:
    """
    Add RegulonDB IDs to existing loaded GenBank gene objects, done by tying a RegulonDB gene name to its locus tag

    :param List[Gene] genbank_genes: a list of previously-loaded from Genbank Gene objects that can be used in lieu of
    making new ones from RegulonDB
    :return List[Gene] updated_genes: updated Genbank-loaded genes with RegulonDB IDs for linking to other objects
    """

    regulon_db_gene_df = load_regulon_db_df('gene')
    locus_tag_gene_name_lookup = load_locus_tag_to_gene_names_lookup()

    for genbank_gene in genbank_genes:

        # if none of the names associated with the locus tag get a hit in the RegulonDB, resort to trying just the
        # common name of the Genbank gene (but must match the exact primary RegulonDB name)
        possible_gene_names = locus_tag_gene_name_lookup.get(genbank_gene.locus_tag, [])
        possible_gene_names.append(genbank_gene.name)

        regulon_db_id = None

        if possible_gene_names is not None:
            for possible_gene_name in possible_gene_names:
                # RegulonDB gene names are unique, so this will be either 0 or 1 rows
                regulon_db_row = regulon_db_gene_df[regulon_db_gene_df['GENE_NAME'] == possible_gene_name]
                if not regulon_db_row.empty:
                    regulon_db_id = regulon_db_row['GENE_ID'].values[0]
                    break

        genbank_gene.id = regulon_db_id

    return genbank_genes


def load_regulons() -> List[Regulon]:
    """
    Loads Regulons from data files parsed from RegulonDB; a Regulon can contain one or more "functions", indicating
    regulatory modes, that are stored in the functions attribute (a dictionary with lists of regulated promoters
    for each function)

    :return List[Regulon] regulons: a list of Regulon objects
    """

    regulon_df = load_regulon_db_df('regulon')

    regulons = []
    for regulon_row in regulon_df.itertuples(index=False):
        regulon_obj = Regulon(regulon_row.REGULON_ID, regulon_row.REGULON_NAME)
        regulons.append(regulon_obj)

    return regulons


def load_regulon_db_df(file_name: str) -> pd.DataFrame:
    """
    Load a DataFrame from a CSV file in the appropriate directory, that was parsed from raw RegulonDB text files by
    assemble_regulon_db_dataframes()

    :param str file_name: the name of the file (NOT including .csv) for the desired DataFrame
    :return pd.DataFrame df: the DataFrame loaded from the CSV with the given file name
    """
    return pd.read_csv(Path(DATA_DIR, 'regulon_db_parsed', file_name + '.csv'))


def assemble_regulon_db_dataframes():
    """
    Assembles and saves merged DataFrames from raw RegulonDB files that have key information for each feature type;
    Only meant to be run to regenerate dataframes that are parsed by core functions
    """

    # set the Path base
    base_output_path = Path(DATA_DIR, 'regulon_db_parsed')

    # --- OPERONS ---

    operon_df = read_regulon_db_file('operon.txt')
    operon_df = operon_df[[
        'OPERON_ID',
        'OPERON_NAME',
        'REGULATIONPOSLEFT',
        'REGULATIONPOSRIGHT',
        'OPERON_STRAND'
    ]]
    operon_df = operon_df.astype({'REGULATIONPOSLEFT': int, 'REGULATIONPOSRIGHT': int})
    operon_df.to_csv(base_output_path.joinpath(Path('operon.csv')))

    # --- TRANSCRIPTION UNITS ---

    # bring in the transcription unit information; need to use a weird side dataframe to get TU left/right positions
    transcription_unit_df = read_regulon_db_file('transcription_unit.txt')
    transcription_unit_df = transcription_unit_df[[
        'TRANSCRIPTION_UNIT_ID',
        'PROMOTER_ID',
        'TRANSCRIPTION_UNIT_NAME',
        'OPERON_ID'
    ]]
    tu_objects_df = read_regulon_db_file('tu_objects_tmp.txt')
    tu_location_df = tu_objects_df[[
        'TRANSCRIPTION_UNIT_ID',
        'TU_POSLEFT',
        'TU_POSRIGHT',
        # space is a huge hack due to read_regulon_db_file not parsing properly with double-digit column ids
        ' TU_OBJECT_STRAND'
    ]].rename(columns={' TU_OBJECT_STRAND': 'TU_STRAND'})
    tu_location_df['TU_STRAND'] = tu_location_df['TU_STRAND'].apply(
        lambda s: 'forward' if s == 'F' else 'reverse'
    )
    tu_location_df = tu_location_df.drop_duplicates('TRANSCRIPTION_UNIT_ID')
    full_tu_df = transcription_unit_df.merge(tu_location_df, on='TRANSCRIPTION_UNIT_ID', how='outer')

    # --- PROMOTERS ---

    # prepare information about all promoters; don't want if they don't have a known TSS
    promoter_df = read_regulon_db_file('promoter.txt')
    promoter_df = promoter_df[[
        'PROMOTER_ID',
        'PROMOTER_NAME',
        'POS_1',
        'SIGMA_FACTOR',
        'PROMOTER_STRAND'
    ]].dropna(subset=['POS_1'])

    # squeeze the POS1 in with the TU df if we can
    tss_df = promoter_df[['PROMOTER_ID', 'POS_1']]
    full_tu_df_with_tss = full_tu_df.merge(tss_df, on='PROMOTER_ID', how='left')
    full_tu_df_with_tss = full_tu_df_with_tss.astype({'TU_POSLEFT': int, 'TU_POSRIGHT': int})
    full_tu_df_with_tss.to_csv(base_output_path.joinpath(Path('transcription_unit.csv')))

    promoter_feature_df = read_regulon_db_file('promoter_feature.txt')
    promoter_feature_df = promoter_feature_df[[
        'PROMOTER_ID',
        'BOX_10_LEFT',
        'BOX_10_RIGHT',
        'BOX_35_LEFT',
        'BOX_35_RIGHT'
    ]].dropna(subset=['BOX_10_LEFT', 'BOX_10_RIGHT', 'BOX_35_LEFT', 'BOX_35_RIGHT'])

    full_promoter_df = promoter_df.merge(promoter_feature_df, on='PROMOTER_ID', how='outer')
    full_promoter_df.to_csv(base_output_path.joinpath(Path('promoter.csv')))

    # --- TERMINATORS ---

    # move through the TU-terminator link to associate terminators with TUs
    terminator_df = read_regulon_db_file('terminator.txt')
    terminator_df = terminator_df[[
        'TERMINATOR_ID',
        'TERMINATOR_POSLEFT',
        'TERMINATOR_POSRIGHT',
        'TERMINATOR_CLASS'
    ]]
    tu_terminator_link_df = read_regulon_db_file('tu_terminator_link.txt')
    tu_terminator_link_df = tu_terminator_link_df[[
        'TRANSCRIPTION_UNIT_ID',
        'TERMINATOR_ID'
    ]]
    full_terminator_df = terminator_df.merge(tu_terminator_link_df, on='TERMINATOR_ID', how='inner')
    full_terminator_df = full_terminator_df.astype({'TERMINATOR_POSLEFT': int, 'TERMINATOR_POSRIGHT': int})
    full_terminator_df.to_csv(base_output_path.joinpath(Path('terminator.csv')))

    # --- GENES ---

    gene_df = read_regulon_db_file('gene.txt')
    gene_df = gene_df[[
        'GENE_ID',
        'GENE_NAME',
        'GENE_POSLEFT',
        'GENE_POSRIGHT',
        'GENE_STRAND'
    ]].dropna(subset=['GENE_NAME', 'GENE_POSLEFT', 'GENE_POSRIGHT'])
    gene_df = gene_df.astype({'GENE_POSLEFT': int, 'GENE_POSRIGHT': int})
    gene_df.to_csv(base_output_path.joinpath(Path('gene.csv')))

    tu_gene_link_df = read_regulon_db_file('tu_gene_link.txt')
    tu_gene_link_df = tu_gene_link_df[[
        'TRANSCRIPTION_UNIT_ID',
        'GENE_ID'
    ]]
    tu_gene_link_df.to_csv(base_output_path.joinpath(Path('tu_gene_link.csv')))

    # --- ATTENUATORS ---

    attenuator_link_df = read_regulon_db_file('attenuator.txt')
    attenuator_link_df = attenuator_link_df[[
        'ATTENUATOR_ID',
        'GENE_ID',
        'ATTENUATOR_TYPE',
        'ATTENUATOR_STRAND'
    ]]
    a_terminator_df = read_regulon_db_file('attenuator_terminator.txt')
    a_terminator_df = a_terminator_df[[
        'A_TERMINATOR_TYPE',
        'A_TERMINATOR_POSLEFT',
        'A_TERMINATOR_POSRIGHT',
        'A_TERMINATOR_ATTENUATOR_ID'
    ]]
    attenuator_df = attenuator_link_df.merge(
        a_terminator_df,
        left_on='ATTENUATOR_ID',
        right_on='A_TERMINATOR_ATTENUATOR_ID',
        how='inner'
    )
    attenuator_df = attenuator_df.astype({'A_TERMINATOR_POSLEFT': int, 'A_TERMINATOR_POSRIGHT': int})
    attenuator_df.to_csv(base_output_path.joinpath(Path('attenuator.csv')))

    # --- SHINE DALGARNO ---

    shine_dalgarno_df = read_regulon_db_file('shine_dalgarno.txt')
    shine_dalgarno_df = shine_dalgarno_df[[
        'SHINE_DALGARNO_ID',
        'GENE_ID',
        'SHINE_DALGARNO_POSLEFT',
        'SHINE_DALGARNO_POSRIGHT'
    ]]
    shine_dalgarno_df = shine_dalgarno_df.astype({'SHINE_DALGARNO_POSLEFT': int, 'SHINE_DALGARNO_POSRIGHT': int})
    shine_dalgarno_df.to_csv(base_output_path.joinpath(Path('shine_dalgarno.csv')))

    # --- RFAM ---

    rfam_df = read_regulon_db_file('rfam.txt')
    rfam_df = rfam_df[[
        'GENE_ID',
        'RFAM_ID',
        'RFAM_TYPE',
        'RFAM_DESCRIPTION',
        'RFAM_POSLEFT',
        'RFAM_POSRIGHT'
    ]]
    rfam_df = rfam_df.astype({'RFAM_POSLEFT': int, 'RFAM_POSRIGHT': int})
    rfam_df.to_csv(base_output_path.joinpath(Path('rfam.csv')))

    # --- TRANSCRIPTION FACTORS/REGULATION ---

    tf_df = read_regulon_db_file('transcription_factor.txt')
    tf_df = tf_df[[
        'TRANSCRIPTION_FACTOR_ID',
        'TRANSCRIPTION_FACTOR_NAME',
        'SITE_LENGTH'
    ]]
    tf_df.to_csv(base_output_path.joinpath(Path('transcription_factor.csv')))

    conformation_df = read_regulon_db_file('conformation.txt')
    conformation_df = conformation_df[[
        'CONFORMATION_ID',
        'TRANSCRIPTION_FACTOR_ID',
        'FINAL_STATE'
    ]]
    conformation_df.to_csv(base_output_path.joinpath(Path('tf_conformation.csv')))

    # now record the sites that actually have positions, and their associated conformation and/or promoter IDs
    tf_site_df = read_regulon_db_file('site.txt')
    tf_site_df = tf_site_df[[
        'SITE_ID',
        'SITE_POSLEFT',
        'SITE_POSRIGHT'
    ]]
    tf_site_df = tf_site_df[tf_site_df['SITE_POSLEFT'] != 0]
    tf_site_df = tf_site_df[tf_site_df['SITE_POSLEFT'] <= tf_site_df['SITE_POSRIGHT']]
    regulatory_interaction_df = read_regulon_db_file('regulatory_interaction.txt')
    regulatory_interaction_df = regulatory_interaction_df[[
        'CONFORMATION_ID',
        'SITE_ID'
    ]].dropna(subset=['CONFORMATION_ID', 'SITE_ID'])
    tf_site_interaction_df = tf_site_df.merge(regulatory_interaction_df, on='SITE_ID', how='left')
    tf_site_interaction_df.to_csv(base_output_path.joinpath(Path('tf_binding_site.csv')))

    tf_site_promoter_link_df = read_regulon_db_file('regulatory_interaction.txt')
    tf_site_promoter_link_df = tf_site_promoter_link_df[[
        'SITE_ID',
        'PROMOTER_ID'
    ]]
    tf_site_promoter_link_df.to_csv(base_output_path.joinpath(Path('tf_site_promoter_link.csv')))

    # --- REGULONS ---

    regulon_df = read_regulon_db_file('regulon_tmp.txt')
    regulon_df = regulon_df[[
        'REGULON_ID',
        'REGULON_NAME'
    ]]
    regulon_df.to_csv(base_output_path.joinpath(Path('regulon.csv')))

    regulon_func_df = read_regulon_db_file('regulon_function_tmp.txt')
    regulon_func_df = regulon_func_df[[
        'REGULON_FUNCTION_ID',
        'REGULON_ID',
        'REGULON_FUNCTION_NAME'
    ]]
    regulon_func_df.to_csv(base_output_path.joinpath(Path('regulon_function.csv')))

    regulon_tf_link_df = read_regulon_db_file('regulon_tf_link_tmp.txt')
    regulon_tf_link_df = regulon_tf_link_df[[
        'REGULON_ID',
        'TRANSCRIPTION_FACTOR_ID'
    ]]
    regulon_tf_link_df.to_csv(base_output_path.joinpath(Path('regulon_tf_link.csv')))

    regulon_prom_link_df = read_regulon_db_file('regulonfuncpromoter_link_tmp.txt')
    regulon_prom_link_df = regulon_prom_link_df[[
        'REGULON_FUNCTION_ID',
        'PROMOTER_ID'
    ]]
    regulon_prom_link_df.to_csv(base_output_path.joinpath(Path('regulon_func_prom_link.csv')))


def read_regulon_db_file(filename: Union[str, Path]) -> pd.DataFrame:
    """
    Given a filename for a raw TXT file from RegulonDB, parses into a DataFrame
    :param Union[str, Path] filename: the filename of the RegulonDB TXT file to read into a pandas DataFrame
    :return pd.DataFrame regulon_df: a pandas DataFrame parsed from a RegulonDB raw TXT file
    """

    full_filename = Path(DATA_DIR, 'regulon_db10.0', filename)

    with open(full_filename, 'r') as f:
        lines = f.readlines()

    comment = True
    i = 0
    while comment:
        if lines[i].startswith('#'):
            i += 1
        else:
            comment = False

    names = [line[5:-1] for line in lines if re.match(r'# \d', line)]
    df = pd.read_csv(full_filename, index_col=None, skiprows=i, sep='\t', header=None, names=names)

    return df.drop_duplicates()
