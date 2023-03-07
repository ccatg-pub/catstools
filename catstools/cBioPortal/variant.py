from typing import Callable, Dict, List, Optional

from catstools.common import CatsMeta, Variant, filter_variant, iter_variants, iter_top_node_elements


def iter_snv(meta: CatsMeta, filters: List[Callable]):
    """Extracting SNV/insertion/deletion gene mutations"""
    for variant in iter_variants(meta, key='shortVariants'):
        snv = Snv(meta, variant)
        yield from filter_variant(snv, filters)


def iter_cna(meta: CatsMeta, filters: List[Callable]):
    """Extracting copyNumberAlterations gene mutations"""
    for variant in iter_variants(meta, key='copyNumberAlterations'):
        cna = Cna(meta, variant)
        yield from filter_variant(cna, filters)


def iter_sv(meta: CatsMeta, filters: List[Callable]):
    """Extracting rearrangements gene mutations"""
    for variant in iter_variants(meta, key='rearrangements'):
        sv = Sv(meta, variant)
        yield from filter_variant(sv, filters)


def iter_other_biomarkers(meta: CatsMeta, filters: List[Callable]):
    """Extracting otherBiomarkers"""
    for variant in iter_top_node_elements(meta, key='otherBiomarkers'):
        other_biomarker = OtherBiomarker(meta, variant)
        yield from filter_variant(other_biomarker, filters)


def iter_expressions(meta: CatsMeta, filters: List[Callable]):
    """Extracting expressions"""
    for sample in iter_top_node_elements(meta, key='expressions'):
        expression = Expression(meta, sample)
        yield from filter_variant(expression, filters)


def filter_fusion(sv: 'Sv'):
    """Provide filter flags for Gene rearrangement type (fusion)"""
    return True if sv.rearrangement_type_org in sv.meta.sv_label.fusion else False


class Snv(Variant):
    """Gene mutation definition class of SNV/insertion/deletion"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)

    def _set_position(self):
        self._chrom, self._pos = self.meta.convert_position(self.chrom_org, self.pos_org)

    @property
    def chrom(self):
        self._set_position()
        return self._chrom

    @property
    def pos(self):
        self._set_position()
        return self._pos

    @property
    def chrom_org(self):
        return self.variant['chromosome']

    @property
    def pos_org(self):
        return self.variant['position']

    @property
    def id(self):
        return "."

    @property
    def ref(self):
        return self.variant['referenceAllele']

    @property
    def alt(self):
        return self.variant['alternateAllele']

    @property
    def qual(self):
        return "."

    @property
    def filter(self):
        return "."

    @property
    def info(self):
        return "."

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']


class Cna(Variant):
    """Gene mutation definition class of copyNumberAlterations"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)
        self.cna_type_org = variant['copyNumberAlterationType']

    @property
    def transcripts(self):
        # Reconstruct transcript with the name of the first gene
        transcripts: Dict[str, Transcript] = dict()
        for ts_elem in self.variant['transcripts']:
            ts = Transcript(ts_elem)
            if ts.gene not in transcripts:
                transcripts[ts.gene] = ts
        return transcripts

    @property
    def cna_type(self):
        # Type CNA (amplification/deletion)
        if self.cna_type_org in self.meta.cna_label.amplification:
            return '2'
        elif self.cna_type_org in self.meta.cna_label.loss:
            return '-2'
        else:
            return ''

    @property
    def gene(self):
        return ','.join([gene for gene in self.transcripts.keys()])

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']


class Sv(Variant):
    """Gene mutation definition class of rearrangements"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)
        self.rearrangement_type_org = self.variant['rearrangementType']

    @property
    def breakends(self) -> List['Breakend']:
        return [Breakend(self.meta, be) for be in self.variant['breakends']]

    @property
    def rearrangement_type(self):
        if self.rearrangement_type_org == 'other':
            return 'rearrangements'
        else:
            return self.rearrangement_type_org

    @property
    def rearrangement_name(self):
        return self.variant.get('rearrangementNames', [''])

    @property
    def dna(self):
        seq_itemids = [sample.itemid for sample in self.seq_samples if sample.nucleic_acid == 'DNA']
        return 'yes' if self.sampleitemid in seq_itemids else 'no'

    @property
    def rna(self):
        seq_itemids = [sample.itemid for sample in self.seq_samples if sample.nucleic_acid == 'RNA']
        return 'yes' if self.sampleitemid in seq_itemids else 'no'

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']


class Breakend:
    """rearrangements breakend definition class"""

    def __init__(self, meta: CatsMeta, breakend: dict):
        self.meta = meta
        self.breakend = breakend
        self.set_position()

    @property
    def transcripts(self):
        # Reconstruct transcript with the name of the first gene
        transcripts: Dict[str, Optional[str]] = dict()
        for ts_elem in self.breakend['transcripts']:
            ts = Transcript(ts_elem)
            if ts.gene not in transcripts:
                transcripts[ts.gene] = ts.tsid
        return transcripts

    @property
    def gene(self):
        return ','.join([gene for gene in self.transcripts.keys()])

    @property
    def transcript_id(self):
        return ','.join([tsid for tsid in self.transcripts.values()])

    def set_position(self):
        self.chrom_org = self.breakend['chromosome']
        self.start_org = self.breakend['startPosition']
        self.chrom, self.start = self.meta.convert_position(self.chrom_org, self.start_org)


class OtherBiomarker(Variant):
    """Biomarker definition classes other than genetic mutation"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)

    @property
    def biomarker_type(self):
        return self.variant['biomarkerType']

    @property
    def metrics(self):
        return self.variant.get('biomarkerMetrics', [])

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']


class Expression(Variant):
    """Expression Information Definition Class"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)

    @property
    def transcripts(self):
        # Reconstruct transcript with the name of the first gene
        transcripts: Dict[str, Optional[str]] = dict()
        for ts_elem in self.variant['transcripts']:
            ts = Transcript(ts_elem)
            if ts.gene not in transcripts:
                transcripts[ts.gene] = ts.tsid
        return transcripts

    @property
    def gene(self):
        return ','.join([gene for gene in self.transcripts.keys()])

    @property
    def expression_level_metrics(self):
        return self.variant['expressionLevelMetrics']

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']


class Transcript:
    """Transcript Information Definition Class"""

    def __init__(self, transcript: dict):
        self.tsid: Optional[str] = transcript.get('transcriptId', '')
        self.gene: str = transcript.get('geneSymbol', '')
