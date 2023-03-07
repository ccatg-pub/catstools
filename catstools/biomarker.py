import logging
from typing import Callable, Dict, List

from catstools.common import CatsMeta, Variant, filter_variant, iter_top_node_elements

logger = logging.getLogger(__name__)


def iter_other_biomarkers(meta: CatsMeta, filters: List[Callable]):
    """Extracting otherBiomarkers"""
    for variant in iter_top_node_elements(meta, key='otherBiomarkers'):
        other_biomarker = OtherBiomarker(meta, variant)
        yield from filter_variant(other_biomarker, filters)


def filter_biomarker_type(other_biomarker: 'OtherBiomarker'):
    """Provide filter flags for biomarker types"""
    return True if not other_biomarker.is_hide else False


def iter_composite_biomarkers(meta: CatsMeta, filters: List[Callable]):
    """Extracting compositeBiomarkers"""
    for variant in iter_top_node_elements(meta, key='compositeBiomarkers'):
        composite_biomarker = Variant(meta, variant)
        yield from filter_variant(composite_biomarker, filters)


def iter_expressions(meta: CatsMeta, filters: List[Callable]):
    """Extracting expressions"""
    for sample in iter_top_node_elements(meta, key='expressions'):
        expression = Expression(meta, sample)
        yield from filter_variant(expression, filters)


def iter_non_human_contents(meta: CatsMeta, filters: List[Callable]):
    """Extracting nonHumanContents"""
    for sample in iter_top_node_elements(meta, key='nonHumanContents'):
        expression = NonHumanContent(meta, sample)
        yield from filter_variant(expression, filters)


class OtherBiomarker(Variant):
    """Biomarker definition classes other than genetic mutation"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)

    @property
    def biomarker_type(self):
        return self.variant['biomarkerType']

    @property
    def is_hide(self) -> bool:
        if self.biomarker_type == 'MSI':
            return self.meta.get_hide_option('hideMsiValue')
        elif self.biomarker_type == 'TMB':
            return self.meta.get_hide_option('hideTmbValue')
        elif self.biomarker_type == 'LOH':
            return self.meta.get_hide_option('hideLohValue')
        return False

    @property
    def origin(self):
        return self.variant.get('biomarkerOrigin', '')

    @property
    def state(self):
        return self.variant.get('state', '')

    @property
    def metrics(self):
        metrics = self.variant.get('biomarkerMetrics', [])
        for metric in metrics:
            if metric['unit'] is None:
                metric['unit'] = ''
        return metrics

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']

    @property
    def output_map(self) -> Dict:
        return {self.biomarker_type: {"origin": self.origin,
                                      "state": self.state,
                                      "metrics": self.metrics}}


class Expression(Variant):
    """Expression Information Definition Class"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)

    @property
    def read_count(self):
        return self.variant.get('readCount', '')

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
    def expression_level_metrics(self):
        return self.variant['expressionLevelMetrics']

    @property
    def sampleitemid(self):
        return self.variant['sampleItemId']

    @property
    def output_std_json(self) -> List[Dict]:
        output = list()
        output.append({
            'expressions': {
                'readCount': self.read_count,
                'transcripts': [
                    {
                        'geneSymbol': gene
                    } for gene in self.transcripts.keys()
                ],
                'expressionLevelMetrics': self.expression_level_metrics
            }
        })
        return output


class NonHumanContent(Variant):
    """NonHumanContent Information Definition Class"""

    def __init__(self, meta: CatsMeta, variant: Dict) -> None:
        super().__init__(meta, variant)


class Transcript:
    """Transcript Information Definition Class"""

    def __init__(self, transcript: dict):
        self.gene: str = transcript['geneSymbol']
