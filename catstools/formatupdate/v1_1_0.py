import itertools
from typing import Dict, List

from catstools.formatupdate.utils import check_version, iter_variants

MAP_ORIGIN = {
    "somatic": "tumor",
    "germline": "normal",
    "likely somatic": "tumor",
    "likely germline": "normal"
}

MAP_METHOD = {
    "DNA-seq": "DNA",
    "RNA-seq": "RNA"
}


class SequencingSample:
    def __init__(self, itemid: str, tumor_or_normal: str, nucleic_acid: str, adding: bool = False) -> None:
        self.itemid = itemid
        self.tumor_or_normal = tumor_or_normal
        self.nucleic_acid = nucleic_acid
        self.adding = adding


def associate_samples(json_cats: Dict):
    # Define dummies for sequencingSamples information
    samples: Dict[str, SequencingSample] = {
        "DNA_tumor": SequencingSample("dummy-sample-1", "tumor", "DNA"),
        "DNA_normal": SequencingSample("dummy-sample-2", "normal", "DNA"),
        "RNA_tumor": SequencingSample("dummy-sample-3", "tumor", "RNA"),
        "RNA_normal": SequencingSample("dummy-sample-4", "normal", "RNA"),
    }

    # Updated sequencingSamples information
    for sample_elem in json_cats.get("sequencingSamples", []):
        nucleic_acid = sample_elem["nucleicAcid"] = MAP_METHOD[sample_elem["testMethod"]]
        tumor_or_normal = sample_elem["tumorOrNormal"]
        del sample_elem["testMethod"]
        # Replace dummies with actual values
        samples[f'{nucleic_acid}_{tumor_or_normal}'] = SequencingSample(sample_elem["itemId"],
                                                                        tumor_or_normal,
                                                                        nucleic_acid,
                                                                        adding=True)

    # Associated with sampleItemId
    for variant in itertools.chain(iter_variants(json_cats, "shortVariants"),
                                   iter_variants(json_cats, "copyNumberAlterations"),
                                   iter_variants(json_cats, "rearrangements")):
        nucleic_acid = MAP_METHOD[variant["testMethod"]]
        tumor_or_normal = MAP_ORIGIN[variant.get("variantOrigin", "somatic")]
        del variant["testMethod"]
        sample = samples[f'{nucleic_acid}_{tumor_or_normal}']
        sample.adding = True
        variant["sampleItemId"] = sample.itemid

    # Associated with sampleItemId (otherBiomarker)
    for biomarker in json_cats.get("otherBiomarkers", []):
        nucleic_acid = MAP_METHOD['DNA-seq']
        tumor_or_normal = MAP_ORIGIN[biomarker.get("biomarkerOrigin", "somatic")]
        sample = samples[f'{nucleic_acid}_{tumor_or_normal}']
        sample.adding = True
        biomarker["sampleItemId"] = sample.itemid

    # If there is nothing to be added to sequencingSamples, a dummy definition of DNA/tumor is
    # targeted for addition
    if all([lambda x: x.adding is False for x in samples]):
        samples["DNA_tumor"].adding = True

    # Add missing information to sequencingSamples in CATS format as dummy
    sample_elems: List = json_cats.setdefault("sequencingSamples", list())
    for sample in samples.values():
        # If a dummy definition is associated, add
        if sample.itemid.startswith('dummy-') and sample.adding:
            sample_elems.append({
                "itemId": sample.itemid,
                "tumorOrNormal": sample.tumor_or_normal,
                "nucleicAcid": sample.nucleic_acid
            })


def edit_gene_pairs(json_cats: Dict):
    for variant in iter_variants(json_cats, "rearrangements"):
        genes = [breakend['transcripts'][0]['geneSymbol'] for breakend in variant['breakends']]
        gene_pairs = variant.get('genePairs', [])

        # If there is no genePairs tag, skip it.
        if not gene_pairs:
            continue

        # Change genePairs tag
        ordered_gene_pairs = list()
        if f"{genes[0]}-{genes[1]}" in gene_pairs:
            ordered_gene_pairs.append([genes[0], genes[1]])
        if f"{genes[1]}-{genes[0]}" in gene_pairs:
            ordered_gene_pairs.append([genes[1], genes[0]])
        del variant['genePairs']
        variant['orderedGenePairs'] = ordered_gene_pairs


def update_biomarker_metrics(json_cats: Dict):
    for biomarker in json_cats.get("otherBiomarkers", []):
        for metric in biomarker.get("biomarkerMetrics", []):
            metric['metricType'] = 'scalar'


def update_meta(json_cats: Dict):
    """Update metaData"""
    json_cats['metaData']['schemaVersion'] = "1.2.0"


def run(json_cats: Dict):
    check_version(json_cats, "1.1.0")
    associate_samples(json_cats)
    edit_gene_pairs(json_cats)
    update_biomarker_metrics(json_cats)
    update_meta(json_cats)
