import json
import os
from abc import ABCMeta, abstractmethod
from typing import Dict, List, Optional, Tuple


def convert_to_cats(dataset, grc_release, test_id, test_type, panel_name, panel_version) -> Dict:
    """Convert from cBioPortal input data to CATS format"""
    # Create objects for each mutation type from various data file elements
    variants = Variants()
    # Initialization of sequentially numbered elements of itemId
    stats = Stats(1)
    # Definition of metaData and testInfo tags
    json_dict = {
        "metaData": {
            "schemaVersion": "1.2.0",
            "referenceGenome": {
                "grcRelease": grc_release
            }
        },
        "testInfo": {
            "testId": test_id,
            "testType": test_type,
            "panelName": panel_name,
            "panelVersion": panel_version
        }
    }
    # Check for the existence of objects for each mutation type, and if so, create the corresponding tag
    sample_item_ids = {}

    if dataset.get("snv") is not None:
        variants.add(Snv(dataset["snv"], stats))
    if dataset.get("cna") is not None:
        variants.add(Cna(dataset["cna"], stats))
    if dataset.get("sv") is not None:
        variants.add(Rearrangement(dataset["sv"], stats))
    if variants.contents:
        json_dict["variants"] = variants.contents
        sample_item_ids.update(variants.sample_item_id)
    if dataset.get("other_biomarker") is not None:
        biomarker = Biomarkers(dataset["other_biomarker"], stats)
        if biomarker.content:
            json_dict["otherBiomarkers"] = biomarker.content
            sample_item_ids.update(biomarker.sample_item_id)
    if dataset.get("expression") is not None:
        expression = Expressions(dataset["expression"], stats)
        if expression.content:
            json_dict["expressions"] = expression.content
            sample_item_ids.update(expression.sample_item_id)

    # Creating sequencingSamples tag
    json_dict["sequencingSamples"] = create_sequencing_samples(sample_item_ids)

    return json_dict


def get_dataset(input_dir) -> Tuple[Optional[str], Dict]:
    """Obtain testId and various data file elements and return them in list type"""
    # Acquisition of various meta file elements
    meta_mutations_extended = get_metamap(input_dir, 'meta_mutations_extended.txt')
    meta_cna_discrete = get_metamap(input_dir, 'meta_cna_discrete.txt')
    meta_sv = get_metamap(input_dir, 'meta_sv.txt')
    meta_clinical_sample = get_metamap(input_dir, 'meta_clinical_sample.txt')
    meta_expression = get_metamap(input_dir, 'meta_expression.txt')

    # Acquisition of various data file elements
    datamap: Dict = {}
    datamap.update(get_datamap(input_dir, meta_mutations_extended.get("data_filename"), 'snv'))
    datamap.update(get_datamap(input_dir, meta_cna_discrete.get("data_filename"), 'cna'))
    datamap.update(get_datamap(input_dir, meta_sv.get("data_filename"), 'sv'))
    datamap.update(get_datamap(input_dir, meta_clinical_sample.get("data_filename"), 'other_biomarker'))
    datamap.update(get_datamap(input_dir, meta_expression.get("data_filename"), 'expression'))

    # Returns testId and various data file elements as a list
    return meta_clinical_sample.get("cancer_study_identifier"), datamap


def get_metamap(input_dir, filename: str) -> Dict:
    """Get the contents of the argument file and return it as a dictionary type"""
    meta: Dict = {}
    input_file = os.path.join(input_dir, filename)
    if os.path.isfile(input_file):
        with open(input_file, "r", encoding="utf-8") as f:
            for line in f.readlines():
                data = line.replace('\n', '').replace(' ', '').split(':')
                meta[data[0]] = data[1]
    return meta


def get_datamap(input_dir, filename: Optional[str], mutation_type) -> Dict:
    """Get the contents of the argument file and return it as a list type"""
    data: Dict = {}
    if filename is not None:
        input_file = os.path.join(input_dir, filename)
        if os.path.isfile(input_file):
            with open(input_file, "r", encoding="utf-8") as f:
                data[mutation_type] = f.readlines()
    return data


def create_sequencing_samples(sample_item_id: Dict):
    """Determines tumorOrNormal and nucleicAcid from arguments and returns elements of sequencingSamples"""
    sample: List = []
    for key, value in sample_item_id.items():
        if value[0] == 'germline':
            origin = 'normal'
        else:
            origin = 'tumor'
        sample.append({
            "itemId": key,
            "tumorOrNormal": origin,
            "nucleicAcid": value[1]
        })
    return sample


def write_cats(output, json_dict):
    """Output CATS format files"""
    with open(output, "w", encoding="utf-8") as f:
        json.dump(json_dict, f, indent=4)


class Variants(metaclass=ABCMeta):
    def __init__(self):
        self.contents = {}
        self.sample_item_id = {}

    def add(self, product: "Variant"):
        if not product.content:
            return
        if product.key in self.contents:
            self.contents[product.key].extend(product.content)
        else:
            self.contents[product.key] = product.content

        self.sample_item_id.update(product.sample_item_id)


class Stats:
    def __init__(self, count):
        self.count = count


class Variant:
    def __init__(self, stats, key):
        self.key = key
        self.content = []
        self.sample_item_id = {}
        self.stats: Stats = stats
        self.build()

    @abstractmethod
    def build(self):
        pass


class Snv(Variant):
    def __init__(self, item_list, stats):
        self.snv_list = item_list[2:]
        super().__init__(stats, "shortVariants")

    def build(self):
        for snv in self.snv_list:
            data = snv.replace('\n', '').split('\t')
            db_name = 'Ensembl' if str(data[37]).startswith('ENST') else 'RefSeq'
            snv_map = {
                "itemId": f"variant-{self.stats.count}",
                "chromosome": data[4],
                "position": int(data[5]),
                "referenceAllele": data[10],
                "alternateAllele": data[12],
                "alternateAlleleFrequency": float(data[76]) if data[76] else 0.0,
                "transcripts": [{
                    "transcriptId": data[37],
                    "transcriptDatabaseName": db_name,
                    "geneSymbol": data[0],
                    "strand": data[7],
                    "cdsChange": data[34],
                    "aminoAcidsChange": data[35],
                }],
                "sampleItemId": data[15]
            }
            if data[25]:
                snv_map["variantOrigin"] = data[25]
                self.sample_item_id.update({data[15]: [data[25], "DNA"]})
            else:
                self.sample_item_id.update({data[15]: ["somatic", "DNA"]})
            snv_map["reported"] = True
            self.content.append(snv_map)
            self.stats.count += 1


class Cna(Variant):
    def __init__(self, item_list: List[str], stats):
        self.test_ids = item_list[0].replace('\n', '').split('\t')[2:]
        self.cna_list: List[str] = item_list[1:]
        super().__init__(stats, "copyNumberAlterations")

    def build(self):
        for cna in self.cna_list:
            data = cna.replace('\n', '').split('\t')
            if data[2] == '-1' or data[2] == '-2':
                cna_type = 'loss'
            elif data[2] == '1' or data[2] == '2':
                cna_type = 'amplification'
            else:
                cna_type = ''
            cna_map = {
                "itemId": f"variant-{self.stats.count}",
                "copyNumberAlterationType": cna_type,
                "transcripts": [{
                    "geneSymbol": data[0],
                }],
                "sampleItemId": self.test_ids[0],
                "reported": True,
            }
            self.sample_item_id.update({self.test_ids[0]: ["somatic", "DNA"]})
            self.content.append(cna_map)
            self.stats.count += 1


class Rearrangement(Variant):
    def __init__(self, item_list: List[str], stats):
        self.sv_list: List[str] = item_list[1:]
        super().__init__(stats, "rearrangements")

    def build(self):
        for sv in self.sv_list:
            data = sv.replace('\n', '').split('\t')
            db_name_1 = 'Ensembl' if data[3].startswith('ENST') else 'RefSeq'
            db_name_2 = 'Ensembl' if data[10].startswith('ENST') else 'RefSeq'
            sv_map = {
                "itemId": f"variant-{self.stats.count}",
                "breakends": [{
                    "chromosome": data[5],
                    "startPosition": int(data[6]),
                    "endPosition": int(data[6]),
                    "transcripts": [{
                        "transcriptId": data[3],
                        "transcriptDatabaseName": db_name_1,
                        "geneSymbol": data[2],
                    }],
                }, {
                    "chromosome": data[12],
                    "startPosition": int(data[13]),
                    "endPosition": int(data[13]),
                    "transcripts": [{
                        "transcriptId": data[10],
                        "transcriptDatabaseName": db_name_2,
                        "geneSymbol": data[9],
                    }],
                }, ],
                "rearrangementType": data[32],
                "sampleItemId": data[0],
                "reported": True,
            }
            if data[17] == 'no' and data[18] == 'yes':
                self.sample_item_id.update({data[0]: ["somatic", "RNA"]})
            else:
                self.sample_item_id.update({data[0]: ["somatic", "DNA"]})
            self.content.append(sv_map)
            self.stats.count += 1


class Biomarkers(Variant):
    def __init__(self, item_list: List[str], stats):
        self.header = item_list[4].replace('\n', '').split('\t')
        self.biomarkers: List[str] = item_list[5:]
        super().__init__(stats, "otherBiomarkers")

    def build(self):
        biomarker_stats_count = 1
        for biomarker in self.biomarkers:
            data_map = dict(zip(self.header, biomarker.replace('\n', '').split('\t')))
            biomarker_map = dict()
            if 'TMB_NONSYNONYMOUS' in self.header and data_map["TMB_NONSYNONYMOUS"]:
                biomarker_map = {
                    "itemId": f"biomarker-{biomarker_stats_count}",
                    "biomarkerType": "TMB",
                    "biomarkerMetrics": [{
                        "value": float(data_map["TMB_NONSYNONYMOUS"]),
                        "unit": None,
                        "metricType": "TMB_NONSYNONYMOUS"
                    }],
                    "sampleItemId": data_map["SAMPLE_ID"],
                    "reported": True,
                }
            elif 'MSI_SENSOR_SCORE' in self.header and data_map["MSI_SENSOR_SCORE"]:
                biomarker_map = {
                    "itemId": f"biomarker-{biomarker_stats_count}",
                    "biomarkerType": "MSI",
                    "biomarkerMetrics": [{
                        "value": float(data_map["MSI_SENSOR_SCORE"]),
                        "unit": None,
                        "metricType": "MSI_SENSOR_SCORE"
                    }],
                    "sampleItemId": data_map["SAMPLE_ID"],
                    "reported": True,
                }
            else:
                break

            self.sample_item_id.update({data_map["SAMPLE_ID"]: ["somatic", "DNA"]})
            self.content.append(biomarker_map)
            biomarker_stats_count += 1


class Expressions(Variant):
    def __init__(self, item_list: List[str], stats):
        self.test_ids = item_list[0].replace('\n', '').split('\t')[2:]
        self.expressions: List[str] = item_list[1:]
        super().__init__(stats, "expressions")

    def build(self):
        expression_stats_count = 1
        for expression in self.expressions:
            data = expression.replace('\n', '').split('\t')
            self.content.append({
                "itemId": f"expression-{expression_stats_count}",
                "transcripts": [{
                    "geneSymbol": data[0]
                }],
                "expressionLevelMetrics": [{
                    "value": float(data[2]),
                    "unit": ""
                }],
                "sampleItemId": self.test_ids[0],
                "reported": True,
            })

            self.sample_item_id.update({self.test_ids[0]: ["somatic", "RNA"]})
            expression_stats_count += 1
