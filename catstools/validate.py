
import itertools
import json
import logging
import os
from typing import Dict

import jsonschema

from catstools.biomarker import (
    iter_composite_biomarkers, iter_expressions, iter_non_human_contents, iter_other_biomarkers)
from catstools.cna import iter_cna
from catstools.common import CatsMeta, filter_reported, root_dir
from catstools.snv import iter_snv
from catstools.sv import filter_except_fusion, filter_fusion, iter_sv

logger = logging.getLogger(__name__)


def check_format(json_cats: Dict):
    """Validate JSON schema in CATS format"""
    with open(os.path.join(root_dir, "reference/schema.json")) as f:
        schema = json.load(f)
        jsonschema.validate(json_cats, schema)


def check_itemids(meta: CatsMeta):
    """Check for duplicate itemId"""
    itemids = []
    for variant in itertools.chain(iter_snv(meta, filters=[filter_reported]),
                                   iter_cna(meta, filters=[filter_reported]),
                                   iter_sv(meta, filters=[filter_reported]),
                                   iter_other_biomarkers(meta, filters=[filter_reported]),
                                   iter_composite_biomarkers(meta, filters=[filter_reported]),
                                   iter_expressions(meta, filters=[filter_reported]),
                                   iter_non_human_contents(meta, filters=[filter_reported])):
        if variant.itemid in itemids:
            raise ValueError(f'ItemId({variant.itemid}) is duplicated.')
        itemids.append(variant.itemid)


def check_sampleitemids(meta: CatsMeta):
    """Check sampleItemId"""
    for variant in itertools.chain(iter_snv(meta, filters=[filter_reported]),
                                   iter_cna(meta, filters=[filter_reported]),
                                   iter_sv(meta, filters=[filter_reported]),
                                   iter_other_biomarkers(meta, filters=[filter_reported]),
                                   iter_expressions(meta, filters=[filter_reported])):
        seq_itemids = [sample.itemid for sample in variant.seq_samples]
        if variant.sampleitemid not in seq_itemids:
            raise ValueError(f'SampleItemId({variant.sampleitemid}) is not defined in sequencingSamples.')


def check_short_variants(meta: CatsMeta):
    """Check for SNV/insertion/deletion genetic mutations."""
    for snv in iter_snv(meta, filters=[filter_reported]):
        validate_position(meta, snv.chrom_org, snv.pos_org, snv.ref)


def validate_position(meta: CatsMeta, chrom: str, pos: int, ref: str):
    """Validate the position's reference base"""
    start = pos - 1
    end = start + len(ref)
    if meta.is_grch37_build():
        fasta_ref = meta.fasta_grch37.fetch(reference=f"{chrom}", start=start, end=end)
    else:
        fasta_ref = meta.fasta_grch38.fetch(reference=f"chr{chrom}", start=start, end=end)
    if ref != fasta_ref:
        raise ValueError(f"Reference allele does not match at chr{chrom}:{pos}. (Expect:{fasta_ref}, Found:{ref})")


def check_fusion(meta: CatsMeta):
    """Check for rearrangements(fusion) genetic mutations."""
    for fusion in iter_sv(meta, filters=[filter_reported, filter_fusion]):
        expected_genes = list(fusion.breakends_org[0].transcripts.keys()) \
            + list(fusion.breakends_org[1].transcripts.keys())
        ts_0 = fusion.breakends_org[0].transcripts
        ts_1 = fusion.breakends_org[1].transcripts
        for gene_pair in fusion.gene_pairs:
            # Confirmation that both genes comprising the fusion gene are defined by geneSymbol.
            # Check under the assumption that the fusion gene is not composed of the same gene.
            if gene_pair[0] == gene_pair[1]:
                raise ValueError(f"The genePairs consists of the same geneSymbol({gene_pair[0]}).")
            if not(gene_pair[0] in ts_0 and gene_pair[1] in ts_1 or gene_pair[0] in ts_1 and gene_pair[1] in ts_0):
                raise ValueError(f"The genePairs({gene_pair}) does not consist of the geneSymbol({expected_genes}).")


def check_except_fusion(meta: CatsMeta):
    """Check for rearrangements(Other than fusion) genetic mutations."""
    for except_fusion in iter_sv(meta, filters=[filter_reported, filter_except_fusion]):
        gene_0 = except_fusion.breakends[0][0].gene
        gene_1 = except_fusion.breakends[0][1].gene
        if not gene_0 and not gene_1:
            raise ValueError("Both breakend gene symbols consist of nulls.")


def check_map_keys(meta: CatsMeta):
    """Check mapping information"""
    maps: Dict[str, str] = dict()
    for variant in itertools.chain(iter_snv(meta, filters=[filter_reported]),
                                   iter_cna(meta, filters=[filter_reported]),
                                   iter_sv(meta, filters=[filter_reported]),
                                   iter_other_biomarkers(meta, filters=[filter_reported])):
        for key, itemid in variant.output_map.items():
            if key in maps:
                raise KeyError(f'Key({key}) already exists. (ItemId: {maps[key]}, {itemid})')
            maps[key] = itemid
