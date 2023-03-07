import itertools
import logging
from typing import Callable, Dict, List, Optional

from catstools.common import CatsMeta, Variant, chrom2num, filter_variant, iter_variants

logger = logging.getLogger(__name__)


def iter_sv(meta: CatsMeta, filters: List[Callable]):
    """Extracting rearrangements gene mutations"""
    for variant in iter_variants(meta, key='rearrangements'):
        if filter_fusion in filters:
            yield from filter_variant(SvFusion(meta, variant), filters)
        elif filter_except_fusion in filters:
            yield from filter_variant(SvExceptFusion(meta, variant), filters)
        else:
            # Build a common class for rearrangements that does not define output
            # specifications if there are neither fusion filters nor other than fusion filters
            yield from filter_variant(Sv(meta, variant), filters)


def filter_fusion(sv: 'Sv'):
    """Provide filter flags for Gene rearrangement type (fusion)"""
    return True if sv.sv_type_org in sv.meta.sv_label.fusion else False


def filter_except_fusion(sv: 'Sv'):
    """Provide filter flags for Gene rearrangement type (Other than fusion)"""
    labels = set(sv.meta.sv_label.amplification + sv.meta.sv_label.loss + sv.meta.sv_label.inversion +
                 sv.meta.sv_label.deletion + sv.meta.sv_label.duplication + sv.meta.sv_label.truncation +
                 sv.meta.sv_label.exon_skipping + sv.meta.sv_label.rearrangement) - set(sv.meta.sv_label.fusion)
    return True if sv.sv_type_org in labels else False


class Sv(Variant):
    """Gene mutation definition class of rearrangement"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)
        self.sv_type_org = self.variant['rearrangementType']

    @property
    def breakends(self):
        breakends: List[List[TranscriptSet, TranscriptSet]] = list()
        if not self.gene_pairs:
            # If orderedGenePairs tag is not present, geneSymbol combination is defined as orderedGenePairs
            _gene_pairs = itertools.product(self.breakends_org[0].transcripts.keys(),
                                            self.breakends_org[1].transcripts.keys())
        else:
            _gene_pairs = self.gene_pairs
        for gene_pair in _gene_pairs:
            if gene_pair[0] in self.breakends_org[0].transcripts:
                breakends.append([TranscriptSet(self.breakends_org[0], gene_pair[0]),
                                  TranscriptSet(self.breakends_org[1], gene_pair[1])])
            elif gene_pair[0] in self.breakends_org[1].transcripts:
                breakends.append([TranscriptSet(self.breakends_org[1], gene_pair[0]),
                                  TranscriptSet(self.breakends_org[0], gene_pair[1])])
            else:
                raise ValueError()
        return breakends

    @property
    def breakends_org(self) -> List['Breakend']:
        return [Breakend(self.meta, be) for be in self.variant['breakends']]

    @property
    def gene_pairs(self) -> List[List[str]]:
        return self.variant.get('orderedGenePairs', [])

    @property
    def sv_type(self):
        # Compare rearrangementType with the prioritized genetic aberration label and return the corresponding marker category
        for key, value in self.meta.sv_label.sv_type_map.items():
            if self.sv_type_org == key and self.sv_type_org in value[1]:
                return value[0]
        # Compares rearrangementType and genetic aberration labels in a round-robin fashion and returns the corresponding marker category
        for value in self.meta.sv_label.sv_type_map.values():
            if self.sv_type_org in value[1]:
                return value[0]
        # If none of the above matches, rearrangementType is returned as is
        return self.sv_type_org

    @property
    def strand(self):
        return '.' if len(self.gene_pairs) == 0 else '+'

    @property
    def id(self):
        return '.'

    @property
    def qual(self):
        return '.'

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']

    @property
    def output_map(self) -> Dict:
        map_dict = dict()
        rtype = self.sv_type
        if rtype in self.meta.sv_label.fusion:
            rtype = 'fusion'
        else:
            rtype = rtype.replace(' ', '-')
        for ts_5, ts_3 in self.breakends:
            chrom_key = f"{ts_5.breakend.chrom}:{ts_3.breakend.chrom}"
            pos_key = f"{str(ts_5.breakend.start)}-{str(ts_5.breakend.end)}:" \
                      f"{str(ts_3.breakend.start)}-{str(ts_3.breakend.end)}"
            if rtype != 'fusion':
                # In gene rearrangements other than fusion, the
                # If the gene name on the 5'side is null, the gene name on the 3'side is used as the key.
                genes = [f"{ts_5.gene}"] if ts_5.gene else [f"{ts_3.gene}"]
            else:
                genes = [f"{ts_5.gene}|{ts_3.gene}"]
                if self.strand == '.' and ts_5.gene != ts_3.gene:
                    # When the transcription direction is bidirectional
                    genes.append(f"{ts_3.gene}|{ts_5.gene}")
            for gene in genes:
                map_dict["_".join([self.nucleic_acid, "somatic", chrom_key,
                                   pos_key, "-", rtype, gene])] = self.itemid
        return map_dict


class SvFusion(Sv):
    """Gene mutation definition class of rearrangement(fusion)"""

    @property
    def output_std(self) -> List[List]:
        output = list()
        for ts_5, ts_3 in self.breakends:
            cytos = []
            for cyto_5, cyto_3 in zip(ts_5.breakend.cytobands, ts_3.breakend.cytobands):
                cytos.append("-".join([cyto_5, cyto_3]))
            cyto = ",".join(cytos)
            output.append([f"chr{ts_5.breakend.chrom}", ts_5.breakend.start, ts_5.breakend.end,
                           f"chr{ts_3.breakend.chrom}", ts_3.breakend.start, ts_3.breakend.end,
                           self.id, self.qual, self.strand, self.strand, "BND", "PASS",
                           "", "", "", "", "", "", ".", ".", "", "",
                           f"{ts_5.gene}::-{ts_3.gene}:: {cyto}"])
        return output


class SvExceptFusion(Sv):
    """Gene mutation definition class of rearrangement(Other than fusion)"""
    @property
    def skipped_exon(self):
        return [SkippedExon(ts) for ts in self.variant.get('skippedExonRanges', [])]

    @property
    def output_std_json(self) -> List[Dict]:
        output = list()
        for ts_5, ts_3 in self.breakends:
            output.append({
                'rearrangements': {
                    'breakends': [
                        {
                            'chromosome': ts.breakend.chrom,
                            'startPosition': ts.breakend.start,
                            'endPosition': ts.breakend.end,
                            'transcripts': [{
                                'geneSymbol': ts.gene
                            }]
                        } for ts in [ts_5, ts_3]
                    ],
                    'skippedExonRanges': [
                        {
                            'exonRange': se.range
                        }
                        for se in self.skipped_exon
                    ],
                    'rearrangementType': self.sv_type,
                    'cytoband': "chr" + str(ts_5.breakend.chrom) + ":" + str(ts_5.breakend.start) +
                                "-chr" + str(ts_3.breakend.chrom) + ":" + str(ts_3.breakend.end)
                }
            })
        return output


class TranscriptSet:
    def __init__(self, breakend: 'Breakend', gene: str) -> None:
        self.breakend = breakend
        self.gene = gene


class Breakend:
    """rearrangements breakend definition class"""

    def __init__(self, meta: CatsMeta, breakend: dict):
        self.meta = meta
        self.breakend = breakend
        self.set_position()
        self.stream: Optional[str] = breakend.get('matePieceLocation')
        self.set_transcripts()

    def set_position(self):
        self.chrom_org = self.breakend['chromosome']
        self.start_org = self.breakend['startPosition']
        self.end_org = self.breakend['endPosition']
        self.chrom, self.start = self.meta.convert_position(self.chrom_org, self.start_org)
        self.end = self.meta.convert_position(self.chrom_org, self.end_org)[1]

    def set_transcripts(self):
        # Reconstruct transcript with the name of the first gene
        self.transcripts = dict()
        for ts_elem in self.breakend['transcripts']:
            ts = Transcript(ts_elem)
            if ts.gene not in self.transcripts:
                self.transcripts[ts.gene] = ts

    @property
    def cytobands(self):
        cytobands = []
        for bed in self.meta.cyto_ref:
            if chrom2num(self.chrom) == bed.nchrom and self.chrom == bed.chrom:
                if self.start < bed.start:
                    break
                elif bed.start <= self.start and self.start <= bed.end:
                    cytobands.append(bed.gene)

        return cytobands


class TranscriptBase:
    """Transcript base Information Definition Class"""

    def __init__(self, transcript: dict):
        self.gene: str = transcript['geneSymbol']


class Transcript(TranscriptBase):
    """Transcript Information Definition Class"""

    def __init__(self, transcript: dict):
        super().__init__(transcript)
        self.cds: Optional[str] = transcript.get('cdsChange')
        self.amino: Optional[str] = transcript.get('aminoAcidsChange')
        self.effects: List[str] = transcript.get('calculatedEffects', [])


class SkippedExon(TranscriptBase):
    """Skipped ExonInformation Definition Class"""

    def __init__(self, skipped_exon: Dict):
        super().__init__(skipped_exon)
        self.range: List[int] = skipped_exon['exonRange']
