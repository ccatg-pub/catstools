import itertools
from typing import Dict

from catstools.formatupdate.utils import check_version, iter_variants


def update_cna_metrics(json_cats: Dict):
    for variant in itertools.chain(iter_variants(json_cats, 'copyNumberAlterations')):
        for metric in variant.get('copyNumberMetrics', []):
            if metric['unit'] == 'absolute copy number':
                metric['unit'] = 'copy number'


def update_meta(json_cats: Dict):
    """Update metaData"""
    json_cats['metaData']['schemaVersion'] = "1.2.1"


def run(json_cats: Dict):
    check_version(json_cats, "1.2.0")
    update_cna_metrics(json_cats)
    update_meta(json_cats)
