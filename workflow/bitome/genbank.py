# built-in modules
import re
from typing import List, Union
from warnings import warn

# third-party modules
from Bio import SeqFeature, SeqRecord
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from CAI import CAI, relative_adaptiveness

# local modules
from bitome.features import Gene, MobileElement, Origin, Protein, RepeatRegion, TRNA
from bitome.sbrg import load_locus_tag_yome_lookup
from bitome.utilities import load_essential_genes, load_locus_tag_cogs_lookup, parse_gempro


# define re-used input/output patterns
GENBANK_FEATURE = Union[Gene, MobileElement, Origin, Protein, RepeatRegion, TRNA]


def load_genbank_features(
            genbank_record: SeqRecord.SeqRecord,
            terminus: Union[SeqFeature.FeatureLocation] = None
        ) -> List[GENBANK_FEATURE]:
    """
    Parses a GenBank record and generates Bitome knowledgebase objects based on the features within the record
    Currently set up to create the following feature types:
        - Gene
        - Protein
        - TRNA
        - MobileElement
        - RepeatRegion
        - Origin

    :param SeqRecord.SeqRecord genbank_record: the Genbank record to parse
    :param Union[SeqFeature.FeatureLocation] terminus: the location of the terminus region for this genome; used to
    determine whether a GenBank feature is on the leading or lagging strand
    :return List[GENBANK_FEATURE] genbank_features: the GenBank-based knowledgebase objects for genomic features
    """

    # some functionality is limited to the E. coli K-12 MG1655 genome annotation currently; set a flag for that
    is_k12 = genbank_record.id == 'NC_000913.3'

    if is_k12:
        locus_tag_cog_lookup = load_locus_tag_cogs_lookup()
        locus_tag_yome_lookup = load_locus_tag_yome_lookup()
        essential_locus_tags = load_essential_genes()
    else:
        locus_tag_cog_lookup = {}
        locus_tag_yome_lookup = {}
        essential_locus_tags = []

    genome_seq = genbank_record.seq

    # separate the gene SeqFeatures and non-gene SeqFeatures from the GenBank record
    gene_seqfeatures = select_seqfeatures(genbank_record.features, 'gene')
    non_gene_seqfeatures = list(set(genbank_record.features).difference(set(gene_seqfeatures)))

    origin_seqfeatures = select_seqfeatures(non_gene_seqfeatures, 'rep_origin')
    origins: list = []
    for origin_seqfeature in origin_seqfeatures:
        origins.append(Origin(
            origin_seqfeature.location,
            genome_seq,
            name=get_seqfeature_qualifier(origin_seqfeature, 'note')
        ))

    genes: list = []
    proteins: list = []
    trnas: list = []
    for gene_seqfeature in gene_seqfeatures:

        locus_tag = get_seqfeature_qualifier(gene_seqfeature, 'locus_tag')
        gene_name = get_seqfeature_qualifier(gene_seqfeature, 'gene')
        gene_location = gene_seqfeature.location
        # pseudogenes have a 'pseudo' key (with empty value) in their qualifiers dictionary
        is_pseudo = 'pseudo' in gene_seqfeature.qualifiers

        # determine if feature is leading, lagging, or in terminus region (if we have that region provided)
        # assumes the origin is at a "higher" position in the linear numbering of the chromosome than terminus
        replication_strand = None
        origin_distance = None
        terminus_distance = None
        if len(origins) == 1 and terminus is not None:

            origin = origins[0]
            gene_start = gene_location.start.position
            gene_end = gene_location.end.position
            gene_strand = gene_location.strand

            # all below descriptions of conditional cases assume we're looking at the genome with origin at the top,
            # and the absolute position numbering goes clockwise; so the origin is 12:00, terminus is 5:30 - 6:30

            # the gene is in the 12:00 - 5:30 region; note, we're assuming that the wraparound point is here (i.e.
            # the spot where base 4.6M and base 1 are adjacent; also assuming that terminus region is 180 degrees
            # from origin, so that the clockwise direction will definitely be shorter
            if gene_start > origin.location.end.position or gene_start < terminus.start.position:

                if gene_strand == 1:
                    replication_strand = 'leading'
                else:
                    replication_strand = 'lagging'

                if gene_start > origin.location.end.position:
                    origin_distance = gene_start - origin.location.end.position
                    terminus_distance = terminus.start.position + (len(genome_seq) - gene_end)
                else:
                    origin_distance = (len(genome_seq) - origin.location.end.position) + gene_start
                    terminus_distance = terminus.start.position - gene_end

            # the gene is in the terminus region between 5:30 and 6:30; can't guarantee if it's leading or lagging
            # also don't assume which direction to origin is closer; distance to terminus is 0 since it's in there
            elif terminus.start.position < gene_start < terminus.end.position:
                replication_strand = 'terminus'
                origin_distance_1 = (len(genome_seq) - origin.location.end.position) + gene_start
                origin_distance_2 = origin.location.start.position - gene_end
                origin_distance = min(origin_distance_1, origin_distance_2)
                terminus_distance = 0

            # the gene is on the left of the clock (6:30 - 12:00)
            elif terminus.end.position < gene_start < origin.location.start.position:
                if gene_strand == 1:
                    replication_strand = 'lagging'
                else:
                    replication_strand = 'leading'
                origin_distance = origin.location.start.position - gene_end
                terminus_distance = gene_start - terminus.end.position

        # isolate the feature that this gene codes; GenBank record separates these; e.g. a gene and its 'CDS' (coding
        # sequence) will be distinct SeqFeature objects when parsed from the GenBank record
        # for coronavirus, we want to ignore mat_peptide for now
        if genbank_record.id == 'NC_045512.2':
            coded_seqfeature = find_locus_tag(locus_tag, non_gene_seqfeatures, ignore_types=['mat_peptide'])
        else:
            coded_seqfeature = find_locus_tag(locus_tag, non_gene_seqfeatures)

        if is_pseudo:
            gene_type = 'pseudo'
        # note; this ignores ONE gene in NC_000913.3; ralA, antisense toxin; TODO don't ignore this
        elif coded_seqfeature is None:
            warn(f'No coded feature found for {locus_tag}; no Gene object created')
            continue
        elif coded_seqfeature.type == 'ncRNA':
            ncrna_class = get_seqfeature_qualifier(coded_seqfeature, 'ncRNA_class')
            if ncrna_class == 'antisense_RNA':
                gene_type = 'antisense_RNA'
            else:
                gene_type = 'ncRNA'
        elif coded_seqfeature.type == 'mat_peptide':
            gene_type = 'CDS'
        # TODO don't ignore variation and mRNA features for lambda phage genome
        elif coded_seqfeature.type in ['variation', 'mRNA']:
            continue
        else:
            gene_type = coded_seqfeature.type

        # use the CDS location if the coded feature is a CDS; TODO this glosses over genes whose mRNA are altered to
        # make the CDS (see lambdap57 for an example)
        if gene_type == 'CDS':
            gene_location = coded_seqfeature.location

        gene = Gene(
            locus_tag,
            gene_type,
            gene_location,
            gene_name,
            genome_seq,
            # these lookups are non-empty only for GenBank record NC_000913.3 (E. coli K-12 MG1655)
            cog=locus_tag_cog_lookup.get(locus_tag, None),
            y_ome=locus_tag_yome_lookup.get(locus_tag, None),
            essential=(locus_tag in essential_locus_tags),
            replication_strand=replication_strand,
            origin_distance=origin_distance,
            terminus_distance=terminus_distance
        )
        genes.append(gene)

        if gene_type == 'CDS':
            protein_name = get_seqfeature_qualifier(coded_seqfeature, 'product')
            protein = protein_from_gene(gene, include_gempro=is_k12, name=protein_name)
            proteins.append(protein)
            gene.link_protein(protein)

        # if we have a gene coding for a tRNA, generate a TRNA object
        if gene_type == 'tRNA':
            trna_name = get_seqfeature_qualifier(coded_seqfeature, 'product')
            trna_note = get_seqfeature_qualifier(coded_seqfeature, 'note')
            if trna_note is None:
                trna_anticodon = None
            else:
                # assumes an anticodon will be somewhere in the note
                trna_anticodon = re.findall(r'[AUCGTaugct]{3}', trna_note)
            trna = TRNA(
                locus_tag,
                gene_location,
                gene.reading_frame,
                genome_seq,
                name=trna_name,
                anticodon=trna_anticodon
            )
            trnas.append(trna)
            gene.link_trna(trna)

    # add CAI for protein-coding genes
    cds_genes = [gene for gene in genes if gene.gene_type == 'CDS']
    calculate_and_add_cai(cds_genes)

    # load mobile element, repeat region, and origin of replication features
    mobile_element_seqfeatures = select_seqfeatures(non_gene_seqfeatures, 'mobile_element')
    mobile_elements: list = []
    for mobile_element_seqfeature in mobile_element_seqfeatures:
        mobile_elements.append(MobileElement(
            mobile_element_seqfeature.location,
            genome_seq,
            name=get_seqfeature_qualifier(mobile_element_seqfeature, 'mobile_element_type')
        ))

    repeat_region_seqfeatures = select_seqfeatures(non_gene_seqfeatures, 'repeat_region')
    repeat_regions: list = []
    for repeat_region_seqfeature in repeat_region_seqfeatures:
        repeat_regions.append(RepeatRegion(
            repeat_region_seqfeature.location,
            genome_seq,
            name=get_seqfeature_qualifier(repeat_region_seqfeature, 'note')
        ))

    all_genbank_features = genes + proteins + trnas + mobile_elements + repeat_regions + origins

    return all_genbank_features


def select_seqfeatures(seqfeatures: List[SeqFeature.SeqFeature], feature_type: str) -> List[SeqFeature.SeqFeature]:
    """
    Given a list of SeqFeatures and a specific type to isolate, returns a list of features of that type

    :param List[SeqFeature.SeqFeature] seqfeatures: a list of GenBank features in Biopython SeqFeature form
    :param str feature_type: a string for the type attribute of a GenBank SeqFeature
    :return List[SeqFeature.SeqFeature] features_of_type: a list of features of the requested type
    """
    return [seqfeature for seqfeature in seqfeatures if seqfeature.type == feature_type]


def find_locus_tag(
            locus_tag: str,
            seqfeatures: List[SeqFeature.SeqFeature],
            ignore_types: List[str] = None
        ) -> SeqFeature.SeqFeature:
    """
    Find and return the first instance of a Genbank SeqFeature with the given locus tag

    :param str locus_tag: the locus tag to hunt for
    :param List[SeqFeature.SeqFeature] seqfeatures: the list of SeqFeatures in which to look for the given locus tag
    :param List[str] ignore_types: use to restrict the types of features to look within
    :return SeqFeature.SeqFeature locus_feature: the SeqFeature object with the requested locus tag; note, will return
    None without error if no feature with the given locus tag is found in the provided SeqFeature list
    """
    matching_seqfeature = None
    for seqfeature in seqfeatures:
        if ignore_types is not None and seqfeature.type in ignore_types:
            continue
        if get_seqfeature_qualifier(seqfeature, 'locus_tag') == locus_tag:
            matching_seqfeature = seqfeature
            break

    return matching_seqfeature


def get_seqfeature_qualifier(seqfeature: SeqFeature.SeqFeature, key: str) -> Union[str, None]:
    """
    Get a non-null attribute from a Biopython SeqFeature object

    :param SeqFeature.SeqFeature seqfeature: the feature object from which to get a qualifier value
    :param str key: the name of the qualifier value to get; a key for the qualifiers attribute, a dictionary
    :return Union[str] value: the value stored in the provided feature's qualifiers dictionary for the given key
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


def protein_from_gene(gene: Gene, include_gempro: bool, name: str = None) -> Protein:
    """
    Create a Protein knowledgebase object based on a pre-created Gene object; computed properties of the protein will
    be included based on the GEM-PRO annotation pipeline by looking up the protein's locus tag (currently only possible
    for E. coli K-12 MG1655)

    :param Gene gene: a Bitome knowledgebase Gene object that codes for the protein for which a Protein object is
    desired
    :param bool include_gempro: indicates if protein annotation information from GEM-PRO should be accessed and included
    via looking up the protein's locus tag
    :param str name: the common name of the Protein
    :return Protein protein: a Protein object representing knowledge about the protein coded for by the provided gene
    """

    # note that the FeatureLocation object is 0-indexed
    locus_tag = gene.locus_tag
    gene_location = gene.location
    gene_left = gene_location.start.position
    gene_right = gene_location.end.position
    gene_strand = gene_location.strand
    if gene_strand == 1:
        protein_left = gene_left
        protein_right = gene_right - 3
    else:
        protein_left = gene_left + 3
        protein_right = gene_right
    protein_location = SeqFeature.FeatureLocation(protein_left, protein_right, gene_strand)

    # handle the error where we get a premature stop codon which is actually a selenocysteine
    coding_sequence = gene.sequence
    try:
        coded_aa_sequence = coding_sequence.translate(to_stop=True, table=11, cds=True)
    except TranslationError:
        coded_aa_sequence = ''

    if include_gempro:
        gempro_df = parse_gempro(locus_tag)
        has_gempro = gempro_df is not None
        if not has_gempro:
            warn(f'No GEM-PRO file found for {locus_tag}')

        if has_gempro and 'seq_residue' in gempro_df.columns:
            gempro_aa_sequence = Seq(''.join(list(gempro_df['seq_residue'])))
        else:
            gempro_aa_sequence = None

        # compare the GEM-PRO sequence to what the coding gene codes for when directly translated; we want to not use
        # any GEM-PRO annotations if the coded sequence is different from the GEM-PRO sequence; empirically have found
        # that the discrepancies are all on GEM-PRO's side (b0064, b0292, etc.)
        bad_genbank_sequence = False
        if gempro_aa_sequence is not None:
            if 'U' in gempro_aa_sequence:
                warn(f'Selenocysteine found in GEM-PRO sequence for {locus_tag}; Genbank translation audit skipped')
                final_aa_sequence = gempro_aa_sequence
            else:
                if gempro_aa_sequence != coded_aa_sequence:
                    warn(f'GEM-PRO sequence and coded sequence are not the same for {locus_tag}')
                    bad_genbank_sequence = True

                final_aa_sequence = coded_aa_sequence
        else:
            final_aa_sequence = coded_aa_sequence

        # pull out the prediction of an exposed or buried residue
        if has_gempro and 'seq_RSA-accpro' in gempro_df.columns and not bad_genbank_sequence:
            amino_acid_exposures = ''.join(list(gempro_df['seq_RSA-accpro']))
        else:
            amino_acid_exposures = None

        # pull out the secondary structure info
        if has_gempro and 'struct_SS-dssp' in gempro_df.columns and not bad_genbank_sequence:
            # ensure we don't have any nan; fill with the 'no structure' symbol, '-'
            secondary_structure = ''.join(list(gempro_df['struct_SS-dssp'].fillna('-')))
        else:
            secondary_structure = None
    else:
        final_aa_sequence = coded_aa_sequence
        amino_acid_exposures = None
        secondary_structure = None

    return Protein(
        locus_tag,
        protein_location,
        gene.reading_frame,
        final_aa_sequence,
        name=name,
        secondary_structure=secondary_structure,
        amino_acid_exposures=amino_acid_exposures
    )


def calculate_and_add_cai(all_genes: List[Gene]):
    """
    Calculates the codon adaptation index (CAI) of the DNA sequence for the provided genes. Uses Benjamin Lee's
    implementation that is based on Biopython's implementation in SeqUtils module but apparently corrects known issues
    https://cai.readthedocs.io/en/latest/api.html; no return value (modifies the gene objects in place)

    :param List[Gene] all_genes: the full reference list for which to calculate CAI and populate the CAI attribute
    """
    sequences = [str(gene.sequence) for gene in all_genes]
    weights = relative_adaptiveness(sequences=sequences, genetic_code=11)
    for gene in all_genes:
        # TODO don't calculate CAI if there's an 'N' in the sequence
        if 'N' not in gene.sequence:
            gene.cai = CAI(str(gene.sequence), weights=weights)
