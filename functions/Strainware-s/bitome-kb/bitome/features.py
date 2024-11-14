# built-in modules
from typing import List, Tuple, Union

# third-party modules
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from Bio.SeqUtils import CodonUsage, GC
import numpy as np

# local modules
from bitome.utilities import get_reading_frame, is_int


class BitomeFeature:

    def __init__(
        self,
        feature_type: str,
        location: Union[CompoundLocation, FeatureLocation] = None,
        left_right: Tuple[int, int] = None,
        strand: int = None,
        reference_sequence: Seq = None,
        name: str = None
    ):
        """
        A superclass used to hold common information for a feature included in the bitome; not meant to be directly
        instantiated

        :param str feature_type: the type of this feature
        :param Union[CompoundLocation, FeatureLocation] location: a Biopython FeatureLocation or CompoundLocation object
        describing this feature's absolute start and end (0-indexed); strand may also be included. An alternative to
        manual providing of left_right/strand, can only be passed if neither of those 2 parameters are set; defaults to
        None
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end);
        defaults to None
        NOTE: these positions are considered to be 1-indexed
        :param int strand: the strand on which the feature is located; +1 for forward, -1 for reverse; defaults to None
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq
        object; defaults to None
        :param str name: the name of this feature; defaults to None
        """

        # make sure that EITHER range OR location was used, NOT both
        if left_right is not None:
            if not isinstance(left_right, tuple):
                raise TypeError(f'Range parameter must be a tuple of 2 integers for {name}')
            elif len(left_right) != 2:
                raise ValueError(f'Range parameter must be a tuple of length 2 for {name}')
            elif not is_int(left_right[0]) or not is_int(left_right[1]):
                raise TypeError(f'Start and end values for range parameter must be integers for {name}')
            elif left_right[0] > left_right[1]:
                raise ValueError(f'First position of range tuple must be <= to second position for {name}')
            elif location is not None:
                raise ValueError(f'Use either the range + strand parameters, or the location parameter for {name}')
            else:
                start, end = left_right
                # FeatureLocation uses 0-indexing; also wants base ints; FeatureLocation will throw an error if strand
                # is not +1, -1, or None, let it do that
                start, end = int(start), int(end)
                self.location = FeatureLocation(start - 1, end, strand)
        elif location is not None:
            if not isinstance(location, FeatureLocation) and not isinstance(location, CompoundLocation):
                raise TypeError(f'Location parameter must be Biopython object for {name}')
            else:
                self.location = location
        else:
            self.location = None

        self.type = feature_type
        self.name = name

        self._set_sequence(reference_sequence=reference_sequence)

    def _set_sequence(self, reference_sequence: Seq = None):
        """
        Extract the sequence from the reference sequence, if present

        :param Seq reference_sequence: the reference sequence this feature corresponds to
        """
        # set the sequence of the feature if we have a location and reference; also set reading frame IFF we know strand
        if self.location is not None and reference_sequence is not None:
            self.sequence = self.location.extract(reference_sequence)
        else:
            self.sequence = None

    def add_strand(self, strand: int, reference_sequence: Seq = None):
        """
        A method to add strand specificity to a feature location if strand was not known at time of object creation;
        updates the location attribute and extracts a sequence if possible from the reference sequence

        :param int strand: the strand of this feature, as an integer; +1 for forward strand, -1 for reverse
        :param Seq reference_sequence: the reference sequence this feature corresponds to; after adding strand,
        reference sequence can be used to extract the exact sequence of this feature
        """
        self.location = FeatureLocation(self.location.start, self.location.end, strand)
        self._set_sequence(reference_sequence=reference_sequence)

    @property
    def gc_content(self):
        """
        Returns the GC-content of the sequence for this feature, rounded to 2 decimal places
        Returns None if the feature does not have a sequence associated
        """
        if self.sequence is not None:
            return np.around(GC(self.sequence), decimals=2)
        else:
            return None


class Gene(BitomeFeature):

    def __init__(
        self,
        locus_tag: str,
        gene_type: str,
        location: Union[FeatureLocation, CompoundLocation],
        name: str,
        reference_sequence: Seq,
        cog: str = None,
        y_ome: str = None,
        essential: bool = False,
        replication_strand: str = None,
        origin_distance: int = None,
        terminus_distance: int = None
    ):
        """
        An object representing a gene sequence along a genome. A gene may code for a protein, for a non-translated RNA
        element such as an rRNA or tRNA, or be a pseudogene. If it is a CDS, the gene object will have a protein
        attribute that contains a Protein object

        :param str locus_tag: the b-number for this gene
        :param str gene_type: the type of gene; could be coding, pseudo, rRNA, etc.
        :param Union[FeatureLocation, CompoundLocation] location: a Biopython FeatureLocation or CompoundLocation object
        describing this gene's absolute start and end (0-indexed) with strand information
        :param str name: the common name of the gene
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq object
        :param str cog: the cluster of orthologous genes (COG) group for the gene; describes the functional category
        :param str y_ome: the y-ome categorization for the gene; indicates if the gene is well-characterized (see y-ome
        publication by Zachary King)
        :param bool essential: indicates if this gene is essential for growth
        :param str replication_strand: indicates which strand this gene is on during chromosomal replication; valid
        values are 'leading', 'lagging', or 'terminus' (if the gene is in the 300 kb terminus region
        :param int origin_distance: the bp distance from this gene to the origin of replication
        :param int terminus_distance: the bp distance from this gene to the terminus region; 0 if the gene is in the
        terminus region
        """

        super().__init__('gene', location=location, reference_sequence=reference_sequence, name=name)
        self.reading_frame = get_reading_frame(self.location, len(reference_sequence))

        # if the gene is a coding sequence, it should have a multiple of 3 length; sequence is set by super init
        if gene_type == 'CDS' and len(self.location) % 3 != 0:
            raise ValueError(locus_tag + ': sequence should have multiple of 3 length if gene is coding')

        self.locus_tag = locus_tag
        self.gene_type = gene_type
        self.cog = cog
        self.y_ome = y_ome
        self.essential = essential
        self.replication_strand = replication_strand
        self.origin_distance = origin_distance
        self.terminus_distance = terminus_distance
        
        # only set by add_regulon_db_gene_ids
        self.id = None

        # only set after calculate_and_add_cai is run
        self.cai = None

        # only set after the appropriate linking functions are run
        self.protein = None
        self.trna = None
        self.transcription_units = []
        self.attenuators = []
        self.riboswitches = []
        self.shine_dalgarno = None
        self.i_modulons = []

    def link_protein(self, protein):
        """
        Link a Protein object to this gene; skips if the Gene already has a Protein
        :param Protein protein: the Protein object to add for this gene
        """
        if self.protein is None:
            self.protein = protein
            protein.link_gene(self)

    def link_trna(self, trna):
        """
        Link a TRNA object to this gene; skips if the Gene already has a Protein
        :param TRNA trna: the TRNA object to add for this gene
        """
        if self.trna is None:
            self.trna = trna
            trna.link_gene(self)

    def link_transcription_unit(self, transcription_unit):
        """
        Link a TranscriptionUnit object to this gene; skips if the TU is already in the list
        Note: meant to be called from TranscriptionUnit.link_genes, so does NOT modify the TranscriptionUnit object
        :param TranscriptionUnit transcription_unit: the TranscriptionUnit object to add for this gene
        """
        if transcription_unit.id not in [tu.id for tu in self.transcription_units]:
            self.transcription_units.append(transcription_unit)

    def link_attenuator(self, attenuator):
        """
        Link transcriptional/translational attenuator sequence present in this gene to the Gene object
        Note: called from Attenuator.link_gene, so does not update Attenuator object (just setting backlink)
        :param Attenuator attenuator: the Attenuator object tied to this gene
        """
        if attenuator.id not in [att.id for att in self.attenuators]:
            self.attenuators.append(attenuator)
            attenuator.link_gene(self)

    def link_riboswitch(self, riboswitch):
        """
        Links a Riboswitch that regulates this gene to the gene object
        Note: intended to be called by Riboswitch.link_gene, so assumes that the riboswitch is already linked to gene
        :param Riboswitch riboswitch: the Riboswitch object to link to this Gene object
        """
        if riboswitch.id not in [rs.id for rs in self.riboswitches]:
            self.riboswitches.append(riboswitch)

    def link_shine_dalgarno(self, shine_dalgarno):
        """
        Links a ShineDalgarno that is present in this gene to the gene object
        Note: intended to be called by ShineDalgarno.link_gene, so assumes that the S-D is already linked to gene
        :param ShineDalgarno shine_dalgarno: the S-D object to link to this Gene object
        """
        self.shine_dalgarno = shine_dalgarno

    def link_i_modulon(self, i_modulon):
        """
        Back-link an IModulon object to this Gene if the i-modulon regulates the Gene
        :param IModulon i_modulon: an IModulon object to link to this gene
        """
        if i_modulon not in self.i_modulons:
            self.i_modulons.append(i_modulon)

    @property
    def codon_usage(self):
        """
        Returns a Counter dictionary of codons and their usage counts within this gene
        """
        codons_dict = CodonUsage.CodonsDict.copy()
        codons = [str(self.sequence[i:i+3]) for i in range(0, len(self.sequence), 3)]
        for codon in codons:
            codons_dict[codon] += 1
        return codons_dict

    @property
    def position_in_operon(self):
        """
        Returns the position of this gene in its operon, w.r.t. transcription; if not a member of an operon, returns 1
        :return int operon_order: an integer indicating the position of the gene (1-indexed) in its operon
        """
        if self.transcription_units:
            tu_lengths = [len(tu.location) for tu in self.transcription_units]
            longest_tu = self.transcription_units[int(np.argmax(tu_lengths))]
            if longest_tu.location.strand == 1:
                gene_starts = sorted([gene.location.start.position for gene in longest_tu.genes])
                this_gene_start = self.location.start.position
            else:
                gene_starts = sorted([gene.location.end.position for gene in longest_tu.genes])
                gene_starts.reverse()
                this_gene_start = self.location.end.position
            position = np.where(np.array(gene_starts) == this_gene_start)[0][0] + 1
        else:
            position = 1
        return position


class Protein(BitomeFeature):

    def __init__(
        self,
        locus_tag: str,
        location: FeatureLocation,
        reading_frame: int,
        amino_acid_sequence: Seq,
        name: str = None,
        secondary_structure: str = None,
        amino_acid_exposures: str = None
    ):
        """
        A Bitome knowledgebase object containing information about a protein; information expected to come from
        GenBank and/or the GEM-PRO protein property annotation pipeline

        :param str locus_tag: the b-number for this protein
        :param FeatureLocation location: a Biopython FeatureLocation object describing this protein's absolute start and
         end (0-indexed) with strand information
        :param int reading_frame: the reading frame of this protein
        :param Seq amino_acid_sequence: the amino acid sequence of this protein, as a Biopython Seq object
        :param str name: the common name of the protein; defaults to None
        :param str secondary_structure: the predicted secondary structure - by amino acid - of this protein; the length
        of this string must match the length of amino_acid_sequence; defaults to None
        :param str amino_acid_exposures: the predicted exposure of the amino acids in this protein; the length of this
        string must match the length of amino_acid_sequence; defaults to None
        """

        super().__init__('protein', location=location, name=name)
        self.reading_frame = reading_frame
        self.locus_tag = locus_tag

        # make sure the sequence, secondary structure, and exposures (if any) are the same lengths
        amino_acid_length = len(amino_acid_sequence)
        aa_matched_lists = [secondary_structure, amino_acid_exposures]
        if all([aa_matched_list is not None for aa_matched_list in aa_matched_lists]):
            if not all(len(aa_matched_list) == amino_acid_length for aa_matched_list in aa_matched_lists):
                raise ValueError(f'Protein {locus_tag}: aa length ({amino_acid_length}) not match sec/exp length')

        self.amino_acid_sequence = amino_acid_sequence
        self.secondary_structure = secondary_structure
        self.amino_acid_exposures = amino_acid_exposures
        # gene is only set once the linking method from the Gene class is run
        self.gene = None

    def link_gene(self, gene: Gene):
        """
        Link the gene coding this protein; does not update if gene is already linked; meant to be called from
        Gene.link_protein(), so does NOT update Gene object

        :param Gene gene: the Gene object to link to this protein
        """
        if self.gene is None:
            self.gene = gene


class TRNA(BitomeFeature):

    def __init__(
        self,
        locus_tag: str,
        location: FeatureLocation,
        reading_frame: int,
        reference_sequence: Seq,
        name: str = None,
        anticodon: str = None
    ):
        """
        A Bitome knowledgebase object containing information about a protein

        :param str locus_tag: the b-number for this tRNA
        :param FeatureLocation location: a Biopython FeatureLocation object describing this tRNA's absolute start and
         end (0-indexed) with strand information
        :param int reading_frame: the reading frame of this tRNA
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq object
        :param str name: the common name of the tRNA; defaults to None
        :param str anticodon: the tRNA's anticodon (in RNA nucleotide form); defaults to None
        """

        super().__init__('tRNA', location=location, reference_sequence=reference_sequence, name=name)
        self.reading_frame = reading_frame
        self.locus_tag = locus_tag
        self.anticodon = anticodon
        self.gene = None

    def link_gene(self, gene: Gene):
        """
        Link the gene coding this tRNA; does not update if gene is already linked; meant to be called from
        Gene.link_trna(), so does NOT update Gene object

        :param Gene gene: the Gene object to link to this tRNA
        """
        if self.gene is None:
            self.gene = gene


class MobileElement(BitomeFeature):

    def __init__(
        self,
        location: FeatureLocation,
        reference_sequence: Seq,
        name: str = None
    ):
        """
        A Bitome knowledgebase object containing information about a mobile element (insertion sequence)

        :param FeatureLocation location: a Biopython FeatureLocation object describing this mobile element's absolute
        start and end (0-indexed), with strand information
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq object
        :param str name: the common name of the mobile element/insertion sequence; defaults to None
        """
        super().__init__('mobile_element', location=location, reference_sequence=reference_sequence, name=name)


class RepeatRegion(BitomeFeature):

    def __init__(
        self,
        location: FeatureLocation,
        reference_sequence: Seq,
        name: str = None
    ):
        """
        A Bitome knowledgebase object containing information about a repeat region

        :param FeatureLocation location: a Biopython FeatureLocation object describing this repeat region's absolute
        start and end (0-indexed), with strand information
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq object
        :param str name: the common name of the repeat region; defaults to None
        """
        super().__init__('repeat_region', location=location, reference_sequence=reference_sequence, name=name)


class Origin(BitomeFeature):

    def __init__(
        self,
        location: FeatureLocation,
        reference_sequence: Seq,
        name: str = None
    ):
        """
        A Bitome knowledgebase object containing information about an origin of replication

        :param FeatureLocation location: a Biopython FeatureLocation object describing this replication origin's
        absolute start and end (0-indexed), with strand information
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq object
        :param str name: the common name of the tRNA; defaults to None
        """
        super().__init__('origin', location=location, reference_sequence=reference_sequence, name=name)


class Operon(BitomeFeature):

    def __init__(
        self,
        operon_id: str,
        left_right: Tuple[int, int],
        strand: int, 
        reference_sequence: Seq, 
        name: str = None
    ):
        """
        Load an operon object that stores its name and position; link TranscriptionUnits to the object using the
        link_tu() method to set the appropriate attributes

        :param str operon_id: the RegulonDB unique ID for this operon; will be used to link TUs later
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute 
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end).
        NOTE: these positions are considered to be 1-indexed
        :param int strand: the strand on which the feature is located; +1 for forward, -1 for reverse; defaults to None
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq object
        :param str name: the name of this feature; defaults to None
        """

        super().__init__(
            'operon',
            left_right=left_right,
            strand=strand,
            reference_sequence=reference_sequence,
            name=name
        )
        self.id = operon_id
        self.reading_frame = get_reading_frame(self.location, len(reference_sequence))

        # this is only set after link_transcription_unit is run at least once
        self.transcription_units = []

    def link_transcription_unit(self, transcription_unit):
        """
        Add a transcription unit to this operon's list of TUs;
        NOTE: meant to be called from TranscriptionUnit.link_operon(), so this does NOT update the TU object

        :param TranscriptionUnit transcription_unit: a TranscriptionUnit object to link to this Operon
        """
        # silently don't do anything if this operon is already linked to a TU with the same ID
        if transcription_unit.id not in [tu.id for tu in self.transcription_units]:
            self.transcription_units.append(transcription_unit)


class TranscriptionUnit(BitomeFeature):

    def __init__(
        self,
        tu_id: str,
        operon_id: str,
        left_right: Tuple[int, int],
        strand: int,
        reference_sequence: Seq,
        tss: int = None,
        tts: int = None,
        promoter_id: str = None,
        name: str = None
    ):
        """
        Instantiate a TranscriptionUnit object storing the TU's position and name information; run the link_operon
        method to link to an operon object; NOTE: strand is not set until an operon is linked!

        :param str tu_id: the RegulonDB unique ID of this transcription unit; used to link to other objects
        :param str operon_id: the RegulonDB unique ID of the operon to which this TranscriptionUnit belongs, used for
        linkage establishment
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end).
        NOTE: these positions are considered to be 1-indexed
        :param int strand: the strand on which the feature is located; +1 for forward, -1 for reverse; defaults to None
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq object
        :param int tss: the TSS (transcription start site) for this TU, if annotated for this TU's promoter;
        will be compared to the left_right to determine the "true" left_right span of this TU before setting location
        :param int tts: the TTS (transcription termination site) for this TU, as determined experimentally by Term-Seq;
        will be compared to the left_right to determine the "true" left_right span of this TU before setting location
        :param str promoter_id: the RegulonDB unique ID of the promoter for this TranscriptionUnit, used for linking
        :param str name: the name of this feature; defaults to None
        """

        final_left = left_right[0]
        final_right = left_right[1]

        # always use an explicit TSS if we have it
        if tss is not None:
            if strand == 1:
                final_left = tss
            else:
                final_right = tss

        # compare the given TSS (if any) and TTS (if any) to left_right
        self._cho_tts = False
        if tts is not None:
            if strand == 1 and tts > left_right[1]:
                final_right = tts
                self._cho_tts = True
            elif strand == -1 and tts < left_right[0]:
                final_left = tts
                self._cho_tts = True

        if strand == 1:
            self.tts = final_right
        else:
            self.tts = final_left

        # set strand later when we link this to an operon
        super().__init__(
            'transcription_unit',
            left_right=(final_left, final_right),
            strand=strand,
            reference_sequence=reference_sequence,
            name=name
        )
        self.id = tu_id
        self.operon_id = operon_id
        self.promoter_id = promoter_id

        # only set after linking to operon, promoter, etc.
        self.operon = None
        self.promoter = None
        self.tss = None
        self.sigma_factors = []
        self.reading_frame = None
        self.genes = []
        self.terminators = []

    def link_operon(self, operon: Operon):
        """
        Link the object for this transcription unit's operon to this TU; also adds this TU to the given operon

        :param Operon operon: an Operon object
        """
        self.operon = operon
        operon.link_transcription_unit(self)

    def link_promoter(self, promoter):
        """
        Link this TU object to its promoter; can handle a value of None for promoter (does nothing); passes TSS from
        promoter object as shortcut for TranscriptionUnit
        :param Promoter promoter: a promoter object that is associated with this TU
        """
        if promoter is not None:
            self.promoter = promoter
            self.tss = promoter.tss
            self.sigma_factors = promoter.sigma_factors
            self.reading_frame = promoter.reading_frame
            promoter.link_transcription_unit(self)

    def link_genes(self, genes: List[Gene]):
        """
        Link this TU object to the genes it transcribes; can handle TUs that do not have a gene associated (empty list);
        sets a backlink in the genes for their transcription_units

        :param List[Gene] genes: a list of Gene objects to link to this TranscriptionUnit
        """

        # do a double-check to make sure we don't add duplicate genes
        for gene in genes:
            if gene.locus_tag is not None:
                if gene.locus_tag not in [gene.locus_tag for gene in self.genes]:
                    self.genes.append(gene)
                    gene.link_transcription_unit(self)
            elif gene.id is not None:
                if gene.id not in [gene.id for gene in self.genes]:
                    self.genes.append(gene)
                    gene.link_transcription_unit(self)

    def link_terminator(self, terminator):
        """
        Link this TU object to a Terminator within it
        Note; meant to be called from Terminator.link_transcription_unit (does not modify Terminator object)
        :param Terminator terminator: the Terminator object to link to this TU
        """
        if terminator.id not in [term.id for term in self.terminators]:
            self.terminators.append(terminator)


class Promoter(BitomeFeature):
    def __init__(
        self,
        promoter_id: str,
        tss: int,
        left_right: Tuple[int, int],
        strand: int,
        reference_sequence: Seq,
        box_10_left_right: Tuple[int, int] = None,
        box_35_left_right: Tuple[int, int] = None,
        sigma_factors: List[str] = None,
        name: str = None
    ):
        """
        A Promoter object that stores location information for a promoter that contains the transcription start site
        for a transcription unit. All promoters should eventually be linked to a transcription unit

        :param str promoter_id: the RegulonDB ID for this promoter (required for linking to other objects)
        :param int tss: the location of the transcription start site (TSS) for this promoter; should be 1-indexed
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute 
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end).
        NOTE: these positions are considered to be 1-indexed
        :param int strand: the strand on which the feature is located; +1 for forward, -1 for reverse; defaults to None
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq
        :param Tuple[int, int] box_10_left_right: the left and right ends of this promoter's -10 box, as a tuple of
        integers. Absolute range is used so these are the "left" and "right" ends of the feature, even if on reverse
        (i.e. start <= end). NOTE: these positions are considered to be 1-indexed
        :param Tuple[int, int] box_35_left_right: the left and right ends of this promoter's -35 box, as a tuple of
        integers. Absolute range is used so these are the "left" and "right" ends of the feature, even if on reverse
        (i.e. start <= end). NOTE: these positions are considered to be 1-indexed
        :param List[str] sigma_factors: a list of sigma factors known to bind to this promoter
        :param str name: the name of this feature; defaults to None
        """

        super().__init__(
            'promoter',
            left_right=left_right,
            strand=strand,
            reference_sequence=reference_sequence,
            name=name
        )
        self.id = promoter_id
        self.tss = tss
        self.tss_location = FeatureLocation(tss-1, tss, strand)
        if sigma_factors is None:
            self.sigma_factors = []
        else:
            self.sigma_factors = sigma_factors

        # set formal locations for the -10 and -35 elements only if present; ASSUMES that both are either set/unset
        if box_10_left_right is None:
            self.box_10_location = None
            self.box_35_location = None
        else:
            # FeatureLocation uses 0-indexing
            box_10_start, box_10_end = (box_10_left_right[0] - 1), box_10_left_right[1]
            self.box_10_location = FeatureLocation(box_10_start, box_10_end, strand)

            box_35_start, box_35_end = (box_35_left_right[0] - 1), box_35_left_right[1]
            self.box_35_location = FeatureLocation(box_35_start, box_35_end, strand)

        # set the promoter's reading frame based on the TSS
        self.reading_frame = get_reading_frame(self.tss_location, len(reference_sequence))

        # only set after linking functions are run
        self.transcription_unit = None
        self.tf_binding_sites = []
        self.regulons = []

    def link_transcription_unit(self, transcription_unit: TranscriptionUnit):
        """
        Link a TranscriptionUnit object to this promoter;
        NOTE: only mean to be called from TranscriptionUnit.link_promoter(); does NOT update the provided TU object
        :param TranscriptionUnit transcription_unit: a TU object that contains this Promoter
        """
        self.transcription_unit = transcription_unit

    def link_binding_sites(self, binding_sites):
        """
        Links TFBindingSite objects that reside within this promoter region
        :param list binding_sites: the binding sites to link to this promoter
        """
        # do a double-check to make sure we don't add duplicate binding sites
        for binding_site in binding_sites:
            if binding_site.id not in [bs.id for bs in self.tf_binding_sites]:
                self.tf_binding_sites.append(binding_site)
                binding_site.link_promoter(self)

    def link_regulon(self, regulon):
        """
        Link a Regulon object to this promoter;
        NOTE: only mean to be called from Regulon.link_promoters(); does NOT update the provided Regulon object

        :param Regulon regulon: a Regulon object that regulates this Promoter
        """
        if regulon.id not in [reg.id for reg in self.regulons]:
            self.regulons.append(regulon)


class Terminator(BitomeFeature):

    def __init__(
        self, 
        terminator_id: str, 
        left_right: Tuple[int, int], 
        terminator_class: str, 
        tu_id: str
    ):
        """
        Create an object for a terminator

        :param str terminator_id: the RegulonDB ID for this terminator (required for linking to other objects)
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute 
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end).
        NOTE: these positions are considered to be 1-indexed
        :param terminator_class:
        :param str tu_id: the RegulonDB ID for this terminator's TU (required for linking to other objects)
        """

        # we will set the strand later with set_strand; has to be taken from a TU that is linked to this object
        super().__init__('terminator', left_right=left_right)
        self.id = terminator_id
        self.terminator_class = terminator_class
        self.tu_id = tu_id

        # set after link_transcription_unit is run
        self.transcription_unit = None

    def link_transcription_unit(self, transcription_unit: TranscriptionUnit, reference_sequence: Seq = None):
        """
        Sets this terminator's transcription unit (assumed to have only one);

        :param TranscriptionUnit transcription_unit: the transcription unit object that contains this Terminator
        :param Seq reference_sequence: the reference sequence for this terminator; can be used with TU's
        strand to explicitly set the sequence of the terminator
        """
        self.transcription_unit = transcription_unit
        self.add_strand(transcription_unit.location.strand, reference_sequence=reference_sequence)
        transcription_unit.link_terminator(self)


class Attenuator(BitomeFeature):

    def __init__(
        self,
        attenuator_id: str,
        left_right: Tuple[int, int],
        strand: int,
        reference_sequence: Seq,
        attenuator_type: str,
        termination_type: str,
        gene_id: str
    ):
        """
        Creates an Attenuator-terminator object that stores the location of an attenuator-terminator for a given TU;
        overlapping terminator/anti-terminator/anti-anti-terminator elements should be considered SEPARATE objects

        :param str attenuator_id: the RegulonDB ID for this attenuator (required for linking to other objects)
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute 
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end).
        NOTE: these positions are considered to be 1-indexed
        :param int strand: the strand on which the feature is located; +1 for forward, -1 for reverse; defaults to None
        :param Seq reference_sequence: the reference sequence this feature corresponds to, as a Biopython Seq
        :param str attenuator_type: whether this is a transcriptional or translational attenuator
        :param str termination_type: the type of terminator loop this is; terminator, anti-, or anti-anti-terminator
        :param str gene_id: the RegulonDB ID for this attenuator's associated gene (required for linking)
        """

        super().__init__(
            'attenuator',
            left_right=left_right,
            strand=strand,
            reference_sequence=reference_sequence
        )

        self.id = attenuator_id
        self.attenuator_type = attenuator_type
        self.termination_type = termination_type
        self.gene_id = gene_id

        # populated by link_gene only
        self.gene = None

    def link_gene(self, gene: Gene):
        """
        Link the Gene this attenuator acts on to this Attenuator object; also adds this attenuator to the Gene

        :param Gene gene: a Gene object
        """
        self.gene = gene
        gene.link_attenuator(self)


class ShineDalgarno(BitomeFeature):

    def __init__(self, shine_dalgarno_id: str, left_right: Tuple[int, int], gene_id: str):
        """
        An object representing a Shine-Dalgarno sequence for a given gene

        :param str shine_dalgarno_id: the unique regulonDB ID for this Shine-Dalgarno sequence
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end).
        NOTE: these positions are considered to be 1-indexed
        :param str gene_id: the RegulonDB unique ID of the gene to which this Shine Dalgarno sequence belongs
        """
        super().__init__('shine_dalgarno', left_right=left_right)
        self.id = shine_dalgarno_id
        self.gene_id = gene_id

        # populated by link_gene only
        self.gene = None

    def link_gene(self, gene: Gene, reference_sequence: Seq = None):
        """
        Link the Gene this Shine-Dalgarno is part of to the S-D object; also adds S-D to the gene

        :param Gene gene: a Gene object
        :param Seq reference_sequence: the reference sequence for this S-D; can be used with gene's
        strand to explicitly set the sequence of the Shine-Dalgarno
        """
        self.gene = gene
        self.add_strand(gene.location.strand, reference_sequence=reference_sequence)
        gene.link_shine_dalgarno(self)


class Riboswitch(BitomeFeature):

    def __init__(self, riboswitch_id: str, left_right: Tuple[int, int], name: str, gene_id: str):
        """
        An object representing a regulatory RNA sequence for a given gene

        :param str riboswitch_id: the RegulonDB unique ID for this riboswitch
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end).
        NOTE: these positions are considered to be 1-indexed
        :param str name: the name of this regulatory RNA
        :param str gene_id: the RegulonDB unique ID of the gene this regulatory RNA regulates sequence belongs
        """
        super().__init__('riboswitch', left_right=left_right, name=name)
        self.id = riboswitch_id
        self.gene_id = gene_id
        # set after link_gene is run
        self.gene = None

    def link_gene(self, gene: Gene, reference_sequence: Seq = None):
        """
        Link the Gene that this riboswitch regulates to this Riboswitch object

        :param Gene gene: a Gene object to link this Riboswitch to
        :param Seq reference_sequence: the reference sequence for this terminator; can be used with TU's
        strand to explicitly set the sequence of the terminator
        """
        self.gene = gene
        self.add_strand(gene.location.strand, reference_sequence=reference_sequence)
        gene.link_riboswitch(self)


class TranscriptionFactor:

    def __init__(self, tf_id: str, name: str, site_length: int = None):
        """
        An object representing a transcription factor known to regulate a given genome; set up to contain lists of
        binding sites corresponding to different final conformations of this TF (i.e. with different bound metabolites)

        :param str tf_id: the unique RegulonDB ID for this transcription factor
        :param str name: the name of the transcription factor (excluding any conformation-specific suffixes)
        :param int site_length: the length of the binding site motif for this transcription factor
        """

        self.type = 'transcription_factor'
        self.id = tf_id
        self.name = name
        self.site_length = site_length

        # initialize the binding sites dictionary; we will always have the "base" case of binding sites for the
        # "unmodified" conformation (i.e. just the TF); so we can add that key now
        self.binding_sites = {name: []}
        self.regulons = []
        self.i_modulons = []

    def add_binding_sites(self, binding_sites, final_conformation: str = None):
        """
        Append binding sites to the binding_sites dictionary for this transcription factor

        :param list binding_sites: a list of already_created TFBindingSite objects
        :param str final_conformation: the name of the final conformation of the TF that the given binding sites
        correspond to; defaults to None, in which case the binding sites are assumed to be for the transcription factor
        without any effector molecules attached
        """

        if final_conformation is None:
            self.binding_sites[self.name] += binding_sites
        elif final_conformation in self.binding_sites:
            self.binding_sites[final_conformation] += binding_sites
        else:
            self.binding_sites[final_conformation] = binding_sites

        # also reverse_link the binding sites to this TF object
        for binding_site in binding_sites:
            binding_site.add_transcription_factor(self, final_conformation)

    def link_regulon(self, regulon):
        """
        Link this TranscriptionFactor to a Regulon object to which it belongs;
        NOTE: intended to be called from Regulon.link_transcription_factors, so assumes Regulon already points here

        :param Regulon regulon: the Regulon to link to this TF
        """
        if regulon.id not in [reg.id for reg in self.regulons]:
            self.regulons.append(regulon)

    def link_i_modulon(self, i_modulon):
        """
        Link this TranscriptionFactor to a IModulon object to which it belongs;
        NOTE: intended to be called from IModulon.link_transcription_factors, so assumes IModulon already points here

        :param IModulon i_modulon: the IModulon to link to this TF
        """
        if i_modulon.id not in [imod.id for imod in self.i_modulons]:
            self.i_modulons.append(i_modulon)


class Regulon:

    def __init__(self, regulon_id: str, name: str):
        """
        An object representing a Regulon known to regulate a given set of genes; set up to contain lists of
        regulated promoters sorted by the function of the regulon w.r.t. that promoter (activator, repressor)

        :param str regulon_id: the unique RegulonDB ID for this regulon
        :param str name: the name of this regulon
        """

        self.type = 'regulon'
        self.id = regulon_id
        self.name = name

        # running the link_promoters method will add entries here
        self.regulated_promoters = {}
        self.transcription_factors = []

    def link_transcription_factors(self, transcription_factors: List[TranscriptionFactor]):
        """
        Link transcription factors that enact this regulon's function to the Regulon

        :param List[TranscriptionFactor] transcription_factors: the TranscriptionFactor objects that form this Regulon
        """
        for transcription_factor in transcription_factors:
            if transcription_factor.id not in [tf.id for tf in self.transcription_factors]:
                self.transcription_factors.append(transcription_factor)
                transcription_factor.link_regulon(self)

    def link_promoters(self, promoters: List[Promoter], function_name: str):
        """
        Link a group of promoters to this Regulon that are regulated by a particular functional mode of the regulon

        :param List[Promoter] promoters: the list of Promoter objects to link to this Regulon
        :param str function_name: the name of the Regulon function that regulates the given promoters
        """
        if function_name not in self.regulated_promoters:
            self.regulated_promoters[function_name] = promoters
            for promoter in promoters:
                promoter.link_regulon(self)


class IModulon:

    def __init__(self, i_modulon_id: str, name: str):
        """
        An object representing an i-modulon calculated to regulate a given set of genes

        :param str i_modulon_id: the unique RegulonDB ID for this regulon
        :param str name: the name of this regulon
        """

        self.type = 'i_modulon'
        self.id = i_modulon_id
        self.name = name

        # running the link_genes method will add entries here
        self.transcription_factors = []
        self.genes = []

    def link_transcription_factors(self, transcription_factors: List[TranscriptionFactor]):
        """
        Link transcription factors that enact this i-modulon's function to the IModulon
        :param List[TranscriptionFactor] transcription_factors: the TranscriptionFactor objects that form this i-modulon
        """
        for transcription_factor in transcription_factors:
            if transcription_factor.id not in [tf.id for tf in self.transcription_factors]:
                self.transcription_factors.append(transcription_factor)
                transcription_factor.link_i_modulon(self)

    def link_genes(self, genes: List[Gene]):
        """
        Link genes that are regulated by this i-modulon

        :param List[Gene] genes: a list of gene objects to associate with this i-modulon
        """
        for gene in genes:
            if gene.locus_tag not in [gene.locus_tag for gene in self.genes]:
                self.genes.append(gene)
                gene.link_i_modulon(self)


class TFBindingSite(BitomeFeature):

    def __init__(self, site_id: str, left_right: Tuple[int, int], conformation_id: str = None):
        """
        A feature representing a genomic sequence that is bound by a transcription factor

        :param str site_id: the RegulonDB unique ID for this binding site
        :param Tuple[int, int] left_right: the left and right ends of this feature, as a tuple of integers. Absolute
        range is used so these are the "left" and "right" ends of the feature, even if on reverse (i.e. start <= end).
        NOTE: these positions are considered to be 1-indexed
        :param str conformation_id: the RegulonDB unique ID for the final state of a TF that binds to this binding site
        """

        super().__init__('tf_binding_site', left_right=left_right)

        self.id = site_id
        self.conformation_id = conformation_id

        # only set after add_transcription_factor is run
        self.transcription_factor = None
        self.tf_name = None
        self.tf_final_conformation = None

        # only set after link_promoter is run
        self.promoters = []

    def add_transcription_factor(self, tf_object: TranscriptionFactor, final_conformation: str = None):
        """
        Links a TranscriptionFactor object to this binding site, updating transcription_factor and tf_name

        :param TranscriptionFactor tf_object: a TranscriptionFactor object
        :param str final_conformation: the name of the final conformation of the transcription factor that binds site
        """

        self.transcription_factor = tf_object
        self.tf_name = tf_object.name
        if final_conformation is None:
            self.tf_final_conformation = tf_object.name
        else:
            self.tf_final_conformation = final_conformation

    def link_promoter(self, promoter: Promoter):
        """
        Link a promoter to this binding site if the binding site resides within the promoter;
        NOTE: meant to be called from Promoter.link_binding_sites, so assumes that promoter obj is already linked to BS

        :param Promoter promoter: the Promoter to back-link to this TFBindingSite
        """
        if promoter.id not in [prom.id for prom in self.promoters]:
            self.promoters.append(promoter)
