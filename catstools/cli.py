#!/usr/bin/env python3

import argparse
import json
import logging
import os
import sys
from copy import deepcopy
from typing import Dict

from catstools import __version__
from catstools.aggregation.main import run as run_aggregation
from catstools.biomarker import filter_biomarker_type, iter_expressions, iter_other_biomarkers
from catstools.cBioPortal.cats2portal import convert_to_portal, write_dataset
from catstools.cBioPortal.portal2cats import convert_to_cats, get_dataset, write_cats
from catstools.cna import filter_cna_type, iter_cna
from catstools.common import Builder, CatsMeta, filter_dna, filter_reported, filter_rna, write_json
from catstools.formatupdate.main import VERSION_LIST
from catstools.formatupdate.main import run as run_format_update
from catstools.snv import filter_germline, filter_somatic, iter_snv, make_header
from catstools.sv import filter_except_fusion, filter_fusion, iter_sv
from catstools.validate import (
    check_except_fusion, check_format, check_fusion, check_itemids, check_map_keys, check_short_variants)

logger = logging.getLogger(__name__)


def validate(args, json_cats: Dict):
    """Validate CATS format"""
    meta = CatsMeta(json_cats, args.fasta_grch38, args.fasta_grch37)
    check_format(json_cats)
    check_itemids(meta)
    check_short_variants(meta)
    check_fusion(meta)
    check_except_fusion(meta)
    check_map_keys(meta)


def format_update(args, json_cats: Dict):
    """Update CATS format"""
    output_cats = deepcopy(json_cats)
    is_updated = run_format_update(output_cats, args.to_version)
    if is_updated:
        write_json(args.output, maps=[output_cats])


def output_snv(args, json_cats: Dict):
    """Perform standardized conversions(SNV/insertion/deletion)"""
    filter_origin = filter_somatic if args.somatic else filter_germline
    meta = CatsMeta(json_cats, args.fasta_grch38)
    header = make_header(meta)
    builder = Builder(meta, iter_snv, filters=[filter_reported, filter_origin])
    builder.lines.sort(key=meta.generate_sort_key)
    builder.write_std(args.output_variant, header=header)
    write_json(args.output_map, maps=[builder.maps])


def output_cna(args, json_cats: Dict):
    """Perform standardized conversions(copyNumberAlterations)"""
    meta = CatsMeta(json_cats)
    builder = Builder(meta, iter_cna, filters=[filter_reported, filter_cna_type])
    builder.write_std(args.output_variant)
    write_json(args.output_map, maps=[builder.maps])


def output_rearrangement(args, json_cats: Dict):
    """Perform standardized conversions(rearrangements)"""
    # DNA/RNA filter
    filter_na = filter_dna if args.dna else filter_rna
    meta = CatsMeta(json_cats)
    if args.fusion:
        # rearrangements (fusion)
        builder = Builder(meta, iter_sv, filters=[filter_reported, filter_fusion, filter_na])
        builder.write_std(args.output_variant)
    else:
        # rearrangements (not fusion)
        builder = Builder(meta, iter_sv, filters=[filter_reported, filter_except_fusion, filter_na])
        write_json(args.output_variant, maps=[builder.lines_json])
    write_json(args.output_map, maps=[builder.maps])


def output_other_biomarker(args, json_cats: Dict):
    """Perform standardized conversions(otherBiomarkers)"""
    meta = CatsMeta(json_cats)
    builder = Builder(meta, iter_other_biomarkers, filters=[filter_reported, filter_biomarker_type])
    write_json(args.output_variant, maps=[builder.maps])


def output_expression(args, json_cats: Dict):
    """Perform standardized conversions(expressions)"""
    meta = CatsMeta(json_cats)
    builder = Builder(meta, iter_expressions, filters=[filter_reported])
    write_json(args.output_variant, maps=[builder.lines_json])


def output_portal(args, json_cats: Dict):
    """Perform conversion from CATS format to cBioPortal input data"""
    meta = CatsMeta(json_cats)
    output_data = convert_to_portal(meta)
    # Creation of cBioPortal input data
    write_dataset(args.output_dir, output_data)


def output_cats(args, json_cats: Dict):
    """Perform conversion from cBioPortal input data to CATS format"""
    test_id, dataset = get_dataset(args.input_dir)
    output_data = convert_to_cats(dataset, args.grc_release, test_id, args.test_type,
                                  args.panel_name, args.panel_version)
    # Creation of CATS format files
    write_cats(args.output, output_data)


def output_aggregation_results(args, json_cats: Dict):
    """Perform aggregate functions"""
    # Run aggregation
    run_aggregation(args.mode, args.input_pattern, args.output_dir)


def get_args():
    """Extract arguments."""
    parser = argparse.ArgumentParser(prog='catstools')
    parser.add_argument("--version", action="version", version='{prog} {version}'.format(
        prog=os.path.basename(sys.argv[0]), version=__version__))
    subparsers = parser.add_subparsers(dest='sub_command')
    subparsers.required = True
    input_help = 'input cats format json file'
    mode_help = 'graph aggregation mode\n' \
                '1: onco_plot\n' \
                '2: circos_plot\n' \
                '3: gene_analysis_plot_amino_acids\n' \
                '4: gene_analysis_plot_dna_except_fusion\n' \
                '5: gene_analysis_plot_cna\n' \
                '6: gene_analysis_plot_fusion'

    sp_update = subparsers.add_parser('format_update', formatter_class=argparse.RawTextHelpFormatter)
    sp_update.add_argument('--input', required=True, help=input_help)
    sp_update.add_argument('--output', required=True, help='output cats format json file')
    sp_update.add_argument('--to-version', help=f'updated cats format version\n[{VERSION_LIST}]')
    sp_update.set_defaults(func=format_update)

    sp_validate = subparsers.add_parser('validation')
    sp_validate.add_argument('--input', required=True, help=input_help)
    sp_validate.add_argument('--fasta-grch38', required=True, help='GRCh38 reference fasta file')
    sp_validate.add_argument('--fasta-grch37', required=True, help='GRCh37(hg19) reference fasta file')
    sp_validate.set_defaults(func=validate)

    sp_snv = subparsers.add_parser('short_variants')
    sp_snv.add_argument('--input', required=True, help=input_help)
    mutual_sp_origin = sp_snv.add_mutually_exclusive_group(required=True)
    mutual_sp_origin.add_argument('--somatic', action='store_true', help='extract somatic')
    mutual_sp_origin.add_argument('--germline', action='store_true', help='extract germline')
    sp_snv.add_argument('--fasta-grch38', required=True, help='GRCh38 reference fasta file')
    sp_snv.add_argument('--output-variant', required=True, help='output vcf file')
    sp_snv.add_argument('--output-map', required=True, help='output reference map file')
    sp_snv.set_defaults(func=output_snv)

    sp_cna = subparsers.add_parser('copy_number_alterations')
    sp_cna.add_argument('--input', required=True, help=input_help)
    sp_cna.add_argument('--output-variant', required=True, help='output cna bed file')
    sp_cna.add_argument('--output-map', required=True, help='output reference map file')
    sp_cna.set_defaults(func=output_cna)

    sp_sv = subparsers.add_parser('rearrangements')
    sp_sv.add_argument('--input', required=True, help=input_help)
    mutual_sv_fusion_or_not = sp_sv.add_mutually_exclusive_group(required=True)
    mutual_sv_fusion_or_not.add_argument('--fusion', action='store_true', help='extract fusion')
    mutual_sv_fusion_or_not.add_argument('--except-fusion', action='store_true', help='extract except fusion')
    mutual_sv_na = sp_sv.add_mutually_exclusive_group(required=True)
    mutual_sv_na.add_argument('--dna', action='store_true', help='extract DNA')
    mutual_sv_na.add_argument('--rna', action='store_true', help='extract RNA')
    sp_sv.add_argument('--output-variant', required=True, help='output rearrangement bedpe file')
    sp_sv.add_argument('--output-map', required=True, help='output reference map file')
    sp_sv.set_defaults(func=output_rearrangement)

    sp_other = subparsers.add_parser('other_biomarkers')
    sp_other.add_argument('--input', required=True, help=input_help)
    sp_other.add_argument('--output-variant', required=True, help='output other biomarker json file')
    sp_other.set_defaults(func=output_other_biomarker)

    sp_expression = subparsers.add_parser('expressions')
    sp_expression.add_argument('--input', required=True, help=input_help)
    sp_expression.add_argument('--output-variant', required=True, help='output expression json file')
    sp_expression.set_defaults(func=output_expression)

    sp_cats_cbioportal = subparsers.add_parser('cats2cBioPortal')
    sp_cats_cbioportal.add_argument('--input', required=True, help=input_help)
    sp_cats_cbioportal.add_argument('--output-dir', required=True, help='output directory path')
    sp_cats_cbioportal.set_defaults(func=output_portal)

    sp_cbioportal_cats = subparsers.add_parser('cBioPortal2cats')
    sp_cbioportal_cats.add_argument('--input-dir', required=True, help='input directory path')
    sp_cbioportal_cats.add_argument('--output', required=True, help='output cats format json file')
    sp_cbioportal_cats.add_argument('--grc-release', default='GRCh38', help='"GRCh37" or "GRCh38"')
    sp_cbioportal_cats.add_argument('--test-type', required=True, help='"tumor-only" or "tumor and matched-normal", '
                                    '"tumor-only (cell-free)", "tumor (cell-free) and matched-normal"')
    sp_cbioportal_cats.add_argument('--panel-name', default='defaultPanel',
                                    help='name of your genomic profiling (gene panel) test.')
    sp_cbioportal_cats.add_argument('--panel-version', default='0.1', help='version of your test')
    sp_cbioportal_cats.set_defaults(func=output_cats)

    sp_aggregation = subparsers.add_parser('aggregation', formatter_class=argparse.RawTextHelpFormatter)
    sp_aggregation.add_argument('--mode', required=True, help=mode_help)
    sp_aggregation.add_argument('--input-pattern', required=True, help='glob pattern to catch input files')
    sp_aggregation.add_argument('--output-dir', required=True, help='output directory path')
    sp_aggregation.set_defaults(func=output_aggregation_results)

    return parser.parse_args()


def main():
    args = get_args()
    json_cats: Dict = dict()
    if args.sub_command != 'cBioPortal2cats' and args.sub_command != 'aggregation':
        with open(args.input, 'r') as f:
            json_cats = json.load(f)
    if hasattr(args, 'func'):
        args.func(args, json_cats)


if __name__ == '__main__':
    main()
