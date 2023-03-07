import logging
import os
from typing import Callable, Dict, List, Optional

from catstools.common import CatsMeta, Variant, filter_variant, iter_variants

logger = logging.getLogger(__name__)


def iter_snv(meta: CatsMeta, filters: List[Callable]):
    """Extracting SNV/insertion/deletion gene mutations"""
    for variant in iter_variants(meta, key='shortVariants'):
        snv = Snv(meta, variant)
        yield from filter_variant(snv, filters)


def filter_somatic(snv: 'Snv'):
    """Provides a filter flag for somatic lineage mutations"""
    return True if snv.origin in ["somatic", "likely somatic"] else False


def filter_germline(snv: 'Snv'):
    """Provides a filter flag for germline mutations"""
    return True if snv.origin in ["germline", "likely germline"] else False


def make_header(meta: CatsMeta):
    """Generate header information"""
    sources = [f"{key}={value}" for key, value in meta.json_cats['testInfo'].items()]
    standard_meta_lines = '''\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=TRANSCRIPT,Number=.,Type=String,Description="Representative transcript id of this mutation">
''' + '\n'.join(f'##contig=<ID={ID},length={length}>' for ID, length in zip(meta.fasta_grch38.references, meta.fasta_grch38.lengths)) + '''
##reference=''' + os.path.basename(meta.fasta_grch38.filename.decode('UTF-8')) + '''
##source=<"''' + ",".join(sources) + '''">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
'''
    return standard_meta_lines


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
    def origin(self):
        return self.variant.get("variantOrigin", "somatic")

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
    def sampleitemid(self):
        return self.variant['sampleItemId']

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
    def af(self):
        # Allerle Frequency
        if not self.meta.get_hide_option('hideAlleleFrequency'):
            return self.variant['alternateAlleleFrequency']
        else:
            return None

    @property
    def dp(self):
        # Reading depth
        if not self.meta.get_hide_option('hideAlleleFrequency'):
            return self.variant.get('totalReadDepth', None)
        else:
            return None

    @property
    def output_std(self) -> List[List]:
        output = list()
        for ts in self.transcripts.values():
            info = []
            if self.af is not None:
                info.append("AF=" + str(self.af))
            if self.dp is not None:
                info.append("DP=" + str(self.dp))
            if self.transcripts is not None:
                info.append("TRANSCRIPT=" + ts.tsid)
            output.append([f"chr{self.chrom}", self.pos, self.id,
                           self.ref, self.alt, self.qual, self.filter, ";".join(info)])
        return output

    @property
    def output_map(self) -> Dict:
        origin = "somatic" if self.origin in ["somatic", "likely somatic"] else "germline"
        return {"_".join([self.nucleic_acid, origin, self.chrom, str(self.pos),
                          self.ref, self.alt]): self.itemid}


class Transcript:
    """Transcript Information Definition Class"""

    def __init__(self, transcript: dict):
        self.tsid: str = transcript['transcriptId']
        self.gene: str = transcript['geneSymbol']
        self.cds: Optional[str] = transcript['cdsChange']
        self.amino: Optional[str] = transcript['aminoAcidsChange']
        self.effects: List[str] = transcript.get('calculatedEffects', [])
