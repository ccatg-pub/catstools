import json
import logging
import os
from typing import Callable, Dict, Iterator, List, Optional, Tuple, Type

import pyliftover
import pysam

logger = logging.getLogger(__name__)
root_dir = os.path.dirname(os.path.realpath(__file__))

DEFAULT_AMPLIFICATION_LABELS = [
    "copyNumberAlterationType: amplification",
    "copyNumberAlterationType: gain",
    "copyNumberAlterationType: duplication"
]
DEFAULT_LOSS_LABELS = [
    "copyNumberAlterationType: loss",
    "copyNumberAlterationType: deletion",
    "copyNumberAlterationType: homozygous deletion"
]
DEFAULT_FUSION_LABELS = [
    "rearrangementType: gene fusion",
    "rearrangementType: gene fusion and frameshift variant",
    "rearrangementType: bidirectional gene fusion"
]
DEFAULT_INVERSION_LABELS = [
    "rearrangementType: inversion"
]
DEFAULT_DELETION_LABELS = [
    "rearrangementType: deletion"
]
DEFAULT_DUPLICATION_LABELS = [
    "rearrangementType: duplication",
    "rearrangementType: tandem duplication"
]
DEFAULT_TRUNCATION_LABELS = [
    "rearrangementType: truncation"
]
DEFAULT_EXON_SKIPPING_LABELS = [
    "rearrangementType: exon skipping"
]
DEFAULT_REARRANGEMENT_LABELS = [
    "rearrangementType: other"
]


class CatsMeta:
    """Meta Information Definition Class"""

    def __init__(self, json_cats: Dict, fasta_grch38: Optional[str] = None, fasta_grch37: Optional[str] = None) -> None:
        self.json_cats = json_cats
        self._load_chain_file()
        self._load_cytoband_file()
        self._load_fasta_file(fasta_grch38, fasta_grch37)
        self.configoptions = self.json_cats['metaData'].get('configOptions', {})
        self.cna_label = Label(self.configoptions, "copyNumberAlterationType")
        self.sv_label = Label(self.configoptions, "rearrangementType")

    def is_grch37_build(self) -> bool:
        build: str = self.json_cats['metaData']['referenceGenome']['grcRelease']
        return True if build.startswith('GRCh37') else False

    def _load_chain_file(self) -> None:
        """Reads a file containing information on the correspondence of gene coordinates between versions"""
        self.liftover = None
        if self.is_grch37_build():
            self.liftover = pyliftover.LiftOver(os.path.join(root_dir, "reference/hg19ToHg38.over.chain"))

    def _load_fasta_file(self, fasta_grch38: Optional[str], fasta_grch37: Optional[str]) -> None:
        """Read reference genome array file"""
        if fasta_grch38:
            self.fasta_grch38 = pysam.FastaFile(fasta_grch38)
        if fasta_grch37:
            self.fasta_grch37 = pysam.FastaFile(fasta_grch37)

    def _load_cytoband_file(self) -> None:
        """Read a file that defines Sideband information"""
        with open(os.path.join(root_dir, "reference/cytoBand.txt")) as f:
            items = [line.strip().split('\t') for line in f.readlines()]
            self.cyto_ref = [Bed(item[0], int(item[1]) + 1, int(item[2]), item[3]) for item in items]

    def get_hide_option(self, key) -> bool:
        """Extract hidden information flags"""
        return self.configoptions.get(key, False)

    def convert_position(self, chrom_org: str, position_org: int) -> Tuple[str, int]:
        """
        Convert genome coordinates

        If the genome version of the position information is GRCh37, the position is lifted over to GRCh38.
        If the genome version is GRCh38, return the original position information received.
        """
        if self.liftover is None:
            return chrom_org, position_org
        else:
            chrom, position = self.liftover.convert_coordinate(f"chr{chrom_org}", position_org)[0][:2]
            return chrom.replace('chr', ''), position

    def generate_sort_key(self, vcfline):
        """Generate Gene mutation information sort key"""
        return (self.fasta_grch38.references.index(vcfline[0]), vcfline[1])


class Label:
    def __init__(self, configoptions, label_key: str) -> None:
        self.configoptions = configoptions
        self.label_key = label_key
        self.amplification = self.get_label_values('typeLabelsInterpretedAsKbAmplification',
                                                   DEFAULT_AMPLIFICATION_LABELS)
        self.loss = self.get_label_values('typeLabelsInterpretedAsKbLoss', DEFAULT_LOSS_LABELS)
        self.deletion = self.get_label_values('typeLabelsInterpretedAsKbDeletion', DEFAULT_DELETION_LABELS)
        self.duplication = self.get_label_values('typeLabelsInterpretedAsKbDuplication', DEFAULT_DUPLICATION_LABELS)
        self.fusion = self.get_label_values('typeLabelsInterpretedAsKbGeneFusion', DEFAULT_FUSION_LABELS)
        self.inversion = self.get_label_values('typeLabelsInterpretedAsKbInversion', DEFAULT_INVERSION_LABELS)
        self.truncation = self.get_label_values('typeLabelsInterpretedAsKbTruncation', DEFAULT_TRUNCATION_LABELS)
        self.exon_skipping = self.get_label_values('typeLabelsInterpretedAsKbExonSkipping',
                                                   DEFAULT_EXON_SKIPPING_LABELS)
        self.rearrangement = self.get_label_values('typeLabelsInterpretedAsKbRearrangement',
                                                   DEFAULT_REARRANGEMENT_LABELS)
        self.sv_type_map = {
            'inversion': ['inversion', self.inversion],
            'deletion': ['long deletion', self.deletion],
            'duplication': ['duplication', self.duplication],
            'tandem duplication': ['duplication', self.duplication],
            'truncation': ['truncation', self.truncation],
            'exon skipping': ['exon skipping', self.exon_skipping],
            'other': ['rearrangements', self.rearrangement]
        }

    def get_label_values(self, configoption_key, default: List[str]) -> List[str]:
        """Extracting the genetic aberration label"""
        label_values = list()
        labels: List[str] = self.configoptions.get(configoption_key) or default
        for label in labels:
            key, value = map(lambda x: x.strip(" "), label.split(":"))
            if key == self.label_key:
                label_values.append(value)
        return label_values


class Builder:
    """Gene mutation information construction class"""

    def __init__(self, meta: CatsMeta,
                 variant_iterator: Callable[[CatsMeta, List[Callable]], Iterator[Type['Variant']]],
                 filters: List[Callable] = list()) -> None:
        self.meta = meta
        self.variant_iterator = variant_iterator
        self.filters = filters
        self.lines: List[List] = list()
        self.lines_json: Dict = dict()
        self.maps: Dict = dict()
        self._build()

    def _build(self):
        """Building Gene mutation Information"""
        for variant in self.variant_iterator(self.meta, self.filters):
            # Build output mutation information(VCF/BED/BEDPE)
            for std in variant.output_std:
                self.lines.append(std)

            # Construct output mutation information (JSON)
            for d in variant.output_std_json:
                for key, item in d.items():
                    self.lines_json.setdefault(key, list()).append(item)

            # Build mapping information
            for key, item in variant.output_map.items():
                self.maps[key] = item

    def write_std(self, file: str, header=None):
        """Output gene mutation information standardized conversion file"""
        with open(file, 'w') as f:
            if header:
                f.write(header)
            for line in self.lines:
                f.write("\t".join(map(str, line)) + "\n")


class Variant:
    """Common classes of Gene mutation information"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        self.meta = meta
        self.variant = variant

        # Sequencing Sample Information
        self.seq_samples: List[SequencingSample] = [sample for sample in iter_sequencing_samples(meta)]

    @property
    def itemid(self) -> str:
        return self.variant["itemId"]

    @property
    def reported(self) -> bool:
        return self.variant["reported"]

    @property
    def sampleitemid(self):
        return None

    @property
    def nucleic_acid(self):
        for sample in self.seq_samples:
            if self.sampleitemid == sample.itemid:
                return sample.nucleic_acid
        return None

    @property
    def output_std(self) -> List[List]:
        # Output information of mutation files (VCF/BED/BEDPE)
        return []

    @property
    def output_map(self) -> Dict:
        # Mapping file output information
        return {}

    @property
    def output_std_json(self) -> List[Dict]:
        # Output information of mutation files (JSON)
        return []


def write_json(map_file, maps: List[Dict]):
    """Output JSON file"""
    json_dict = dict()
    if os.path.exists(map_file):
        # If a mapping file already exists, read it in for appending
        with open(map_file, "r", encoding="utf-8") as f:
            json_dict = json.load(f)
    for m in maps:
        json_dict.update(m)
    with open(map_file, "w", encoding="utf-8") as f:
        json.dump(json_dict, f, indent=4, ensure_ascii=False)


def iter_variants(meta: CatsMeta, key: str):
    """Extracting genetic mutations"""
    variants_set: Dict = meta.json_cats.get('variants', {})
    for variant in variants_set.get(key, []):
        yield variant


def iter_top_node_elements(meta: CatsMeta, key: str):
    """Extract elements from the top-level node of the specified key"""
    for element in meta.json_cats.get(key, []):
        yield element


def iter_sequencing_samples(meta: CatsMeta):
    """Extract sequencingSamples information"""
    for sample in iter_top_node_elements(meta, key='sequencingSamples'):
        yield SequencingSample(sample)


def filter_reported(variant: Variant):
    """Provide a filter flag for report posting availability"""
    return True if variant.reported else False


def filter_dna(variant: Variant):
    """Provide DNA filter flags"""
    seq_itemids = [sample.itemid for sample in variant.seq_samples if sample.nucleic_acid == 'DNA']
    return True if variant.sampleitemid in seq_itemids else False


def filter_rna(variant: Variant):
    """Provide RNA filter flags"""
    seq_itemids = [sample.itemid for sample in variant.seq_samples if sample.nucleic_acid == 'RNA']
    return True if variant.sampleitemid in seq_itemids else False


def filter_variant(variant, filters: List[Callable]):
    """Filter for genetic mutations"""
    if filters:
        if all(f(variant) for f in filters):
            yield variant
    else:
        yield variant


class SequencingSample:
    def __init__(self, sample) -> None:
        self.sample = sample
        self.itemid = sample['itemId']
        self.tumor_or_normal = sample['tumorOrNormal']
        self.nucleic_acid = sample['nucleicAcid']


class GenomicRegion:
    """Genomic Region Information Store Class

    Attributes:
        chrom (str): chromosome which holds the region.
        nchrom (int): numerical representation for chromosome
        start (int): start position
        end (int): end position

    Todo:
        * define common interface function for subclasses (using ``NotImplemented``)
    """

    def __init__(self, chrom, start, end):
        if chrom[0:3] == 'chr':
            self.chrom = chrom[3:]
        else:
            self.chrom = chrom
        self.nchrom = chrom2num(chrom)
        if start < end:
            self.start = start
            self.end = end
        else:
            self.start = end
            self.end = start

    def __repr__(self):
        return '\t'.join([self.chrom, str(self.start), str(self.end)])


class Bed(GenomicRegion):
    """Bed Information Store Class

    Attributs:
        gene (str): name of gene.
        bed_type (int): 1 is mutation target region
                        2 is criterion region for cnv
                        3 is fusion target region
    """

    def __init__(self, chrom, start, end, gene):
        super(Bed, self).__init__(chrom, start, end)
        self.gene = gene
        self.bed_type = 1

    def __repr__(self):
        bed_repr = super(Bed, self).__repr__()
        return '\t'.join([bed_repr, self.gene, str(self.bed_type)])


def chrom2num(chrom):
    """convert chromosome to integral to export file in ascending order
       integral mean alignment sequence
    Args:
        chrom (str): name of chromosome you want to convert.
    Returns:
        int: integral which represents chromosome.
    """
    nchrom = 0
    if chrom.upper() in ['X', 'CHRX', 'CHR96']:
        nchrom = 96
    elif chrom.upper() in ['Y', 'CHRY', 'CHR97']:
        nchrom = 97
    elif chrom.upper() in ['M', 'MT', 'CHRM', 'CHRMT', 'CHR98']:
        nchrom = 98
    elif chrom != '*' and '_' not in chrom:
        if chrom[0:3].upper() == 'CHR':
            try:
                nchrom = int(chrom[3:])
            except Exception:
                nchrom = 99
        else:
            try:
                nchrom = int(chrom)
            except Exception:
                nchrom = 99
    else:
        nchrom = 99
    return nchrom
