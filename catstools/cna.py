import logging
from typing import Callable, Dict, List, Optional

from catstools.common import CatsMeta, Variant, filter_variant, iter_variants

logger = logging.getLogger(__name__)


def iter_cna(meta: CatsMeta, filters: List[Callable]):
    """Extracting copyNumberAlterations gene mutations"""
    for variant in iter_variants(meta, key='copyNumberAlterations'):
        cna = Cna(meta, variant)
        yield from filter_variant(cna, filters)


def filter_cna_type(cna: 'Cna'):
    """Provide CNA type filter flags"""
    return True if cna.cna_type else False


class Cna(Variant):
    """Gene mutation definition class of copyNumberAlterations"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)
        self.cna_type_org = variant['copyNumberAlterationType']

    def _set_position(self):
        # Chromosome number Start position
        self.chrom_org = self.variant.get('chromosome')
        self.start_org = self.variant.get('startPosition')
        if self.chrom_org is not None and self.start_org is not None:
            self._chrom, self._start = self.meta.convert_position(self.chrom_org, self.start_org)
        else:
            if self.chrom_org is not None:
                self._chrom = self.chrom_org
            else:
                self._chrom = "."
            if self.start_org is not None:
                self._start = self.start_org
            else:
                self._start = 0

        # End position
        self.end_org = self.variant.get('endPosition')
        if self.end_org is not None:
            if self.chrom_org is not None:
                self._end = self.meta.convert_position(self.chrom_org, self.end_org)[1]
            else:
                self._end = self.end_org
        else:
            self._end = 0

    @property
    def chrom(self):
        self._set_position()
        return self._chrom

    @property
    def start(self):
        self._set_position()
        return self._start

    @property
    def end(self):
        self._set_position()
        return self._end

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
            return 'Amplification'
        elif self.cna_type_org in self.meta.cna_label.loss:
            return 'Deletion'
        else:
            return None

    @property
    def genes(self):
        return [gene for gene in self.transcripts.keys()]

    @property
    def ratio(self):
        # Metrics
        metrics = self.variant.get('copyNumberMetrics')
        if self.meta.get_hide_option('hideCnaValue') or metrics is None:
            return "."
        else:
            return metrics[0]["value"]

    @property
    def strand(self):
        return '.'

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']

    @property
    def output_std(self) -> List[List]:
        gene = ",".join(self.genes)
        return [[f"chr{self.chrom}", self.start, self.end, f"{self.cna_type} {gene}", self.ratio, self.strand]]

    @property
    def output_map(self) -> Dict:
        map_dict = dict()
        pos_key = f"{str(self.start)}-{str(self.end)}"
        for gene in self.genes:
            map_dict["_".join([self.nucleic_acid, "somatic", self.chrom,
                              pos_key, "-", self.cna_type, gene])] = self.itemid
        return map_dict


class Transcript:
    """Transcript Information Definition Class"""
    def __init__(self, transcript: dict):
        self.gene: str = transcript['geneSymbol']
        self.cds: Optional[str] = transcript.get('cdsChange')
        self.amino: Optional[str] = transcript.get('aminoAcidsChange')
        self.effects: List[str] = transcript.get('calculatedEffects', [])
