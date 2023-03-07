from typing import Dict

from catstools.formatupdate.utils import check_version, iter_variants


def update_rearrangements(json_cats: Dict):
    """Update mutations in rearrangements"""
    for variant in iter_variants(json_cats, 'rearrangements'):
        if variant['rearrangementType'] == 'frameshift gene fusion':
            variant['rearrangementType'] = 'gene fusion and frameshift variant'


class Label:
    """Class for analyzing genetic aberration labels"""

    def __init__(self, label) -> None:
        self.key, self.value = [s.strip(" ") for s in label.split(":")]

    def __str__(self) -> str:
        return f"{self.key}: {self.value}"


def update_meta(json_cats: Dict):
    """Update metaData"""
    json_cats['metaData']['schemaVersion'] = "1.1.0"
    ref_genome: dict = json_cats['metaData']['referenceGenome']
    name = ref_genome['name']
    patch = ref_genome['patch']
    ref_genome['grcRelease'] = f"{name}.{patch}" if patch else name
    del ref_genome['patch']

    # Update genetic abnormality labels
    configoptions = json_cats['metaData'].get('configOptions')
    if configoptions:
        fusion_labels = configoptions.get('typeLabelsInterpretedAsKbGeneFusion', [])
        for index, label in enumerate(fusion_labels):
            label_obj = Label(label)
            if label_obj.value == 'frameshift gene fusion':
                label_obj.value = 'gene fusion and frameshift variant'
                fusion_labels[index] = str(label_obj)


def run(json_cats: Dict):
    check_version(json_cats, "1.0.1")
    update_meta(json_cats)
    update_rearrangements(json_cats)
