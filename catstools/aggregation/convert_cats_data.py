#!/usr/bin/env python
# coding: utf-8
"""
Summary: Obtains various information from the CATS format,
        converts it to TSV format, and outputs it.
"""
import csv
import datetime
import glob
import json
import os

from .aggregation import Aggregation as Aggr

# Various consts
# Key to each dictionary
TEST_INFO = "testInfo"
PANEL_NAME = "panelName"
# Key for mutation type
VARIANTS = "variants"
OTHER_BIOMARKERS = "otherBiomarkers"
VARIANT_SNV = "shortVariants"
VARIANT_CNA = "copyNumberAlterations"
VARIANT_SV = "rearrangements"
# Key required to obtain
BIOMARKER_TYPE = "biomarkerType"
BIOMARKER_METRICS = "biomarkerMetrics"
MSI_VALUE = "value"
STATE = "state"
TMB_VALUE = "value"
VARIANT_ORIGIN = "variantOrigin"
CHROMOSOME = "chromosome"
SNV_POSITION = "position"
REFERENCE_ALLELE = "referenceAllele"
ALTERNATE_ALLELE = "alternateAllele"
START_POSITION = "startPosition"
END_POSITION = "endPosition"
ALLELE_FREQ = "alternateAlleleFrequency"
VARIANT_TYPE = "variantType"
TRANSCRIPTS = "transcripts"
GENE_SYMBOL = "geneSymbol"
AMINO_ACIDS_CHANGE = "aminoAcidsChange"
CDS_CHANGE = "cdsChange"
CNA_TYPE = "copyNumberAlterationType"
COPY_NUMBER_METRICS = "copyNumberMetrics"
COPY_VALUE = "value"
COPY_UNIT = "unit"
REARRANGEMENT_TYPE = "rearrangementType"
REARRANGEMENT_NAMES = "rearrangementNames"
BREAKENDS = "breakends"
GENE_PAIRS = "genePairs"
ORDERED_GENE_PAIRS = "orderedGenePairs"
# Target-determination key
REPORTED = "reported"
# Somatic variant
SOMATIC = "somatic"
# Mutation type
SNV = "SNV"
FUSION = "Fusion"
# CNA value type
FOLD_CHANGE = "fold-change"
ABSOLUTE_COPY_NUMBER = "absolute copy number"
# Value corresponding to various keys
MSI = "MSI"
TMB = "TMB"
CNA_AMPLIFICATION = "amplification"
CNA_DELETION = "loss"
GENE_FUSION = "gene fusion"
# CNA conversion dictionary
CNA_CONVERT_DICT = {
    CNA_AMPLIFICATION: "Amp",
    CNA_DELETION: "Loss",
}
# File extension
EXT_CATS = ".cats"
# Result file name
SUCCESS_CATS_TSV = "cats_file_list.tsv"
INVALID_CATS_TSV = "invalid_cats_file_list.tsv"
# Result file header
ITEM_NUMBER = "#"
HEADER_SAMPLE_ID = "sample_id"
HEADER_CATS_FILE = "cats_file_path"
# sample_id value
DIRECTLY_ARGUMENT = "DIRECTLY_ARGUMENT"


class CatsConverter:
    def __init__(self):
        # Initialize each variable with "-" to assume "-" for values
        # that do not exist for some mutation types.
        self.sample_id = "-"
        self.panel_name = "-"
        self.origin = "-"
        self.variant_type = "-"
        self.state = "-"
        self.marker_name = "-"
        self.amino_change = "-"
        self.cds_change = "-"
        self.ref_allele = "-"
        self.alt_allele = "-"
        self.chr1 = "-"
        self.start_pos1 = "-"
        self.end_pos1 = "-"
        self.value = "-"
        self.chr2 = "-"
        self.start_pos2 = "-"
        self.end_pos2 = "-"

    def create_variant_element_list(self):
        """
        List various elements of one mutation obtained from the CATS format.
        :return: List of various elements of one mutation retrieved
                from CATS format.
        """
        element_list = [self.sample_id, self.panel_name, self.origin,
                        self.variant_type, self.state, self.marker_name,
                        self.amino_change, self.cds_change, self.ref_allele,
                        self.alt_allele, self.chr1, self.start_pos1,
                        self.end_pos1, self.value, self.chr2,
                        self.start_pos2, self.end_pos2]

        return element_list

    def get_variant_origin(self, var_dict):
        """
        Obtain Somatic or Germline information from the CATS format.
        :param var_dict: Various information on one mutation.
        :return: Somatic or Germline information
        """
        if VARIANT_ORIGIN in var_dict.keys():
            self.origin = var_dict[VARIANT_ORIGIN]
        # If "variantOrigin" key does not exist,
        # the value of origin_type is set to somatic.
        else:
            self.origin = SOMATIC

        return self.origin

    def get_sv_marker_name(self, var_dict):
        """
        Get the marker name of SV mutation from CATS format.
        :param var_dict: Various information on one SV mutation.
        :return: Marker name for SV mutation
        """
        # If "genePairs" key exists
        if GENE_PAIRS in var_dict.keys():
            self.marker_name = var_dict[GENE_PAIRS][0]
        # If "orderedGenePairs" key exists
        elif ORDERED_GENE_PAIRS in var_dict.keys():
            a_marker_name = var_dict[ORDERED_GENE_PAIRS][0][0]
            b_marker_name = var_dict[ORDERED_GENE_PAIRS][0][1]
            self.marker_name = f"{a_marker_name}-{b_marker_name}"
        # other
        else:
            a_marker_name = var_dict[BREAKENDS][0][TRANSCRIPTS][0][GENE_SYMBOL]
            b_marker_name = var_dict[BREAKENDS][1][TRANSCRIPTS][0][GENE_SYMBOL]
            self.marker_name = f"{a_marker_name}-{b_marker_name}"

        return self.marker_name

    def get_sv_chromosome_position(self, var_dict, marker_name):
        """
        Obtain chromosome number and position information for the SV mutation
        from the CATS format.
        :param var_dict: Various information on one SV mutation.
        :param marker_name: Marker name for SV mutation
        :return: Chromosome number and position information for SV mutations
        """
        # Only gene names are extracted
        if "fusion" in marker_name:
            marker_name = marker_name.replace(" fusion", "")

        # Determination of genes A and B when "genePairs" key is present.
        a_marker_name = "-"
        b_marker_name = "-"
        if "-" in marker_name:
            split_marker_name = marker_name.split("-")
            a_marker_name = split_marker_name[0]
            b_marker_name = split_marker_name[1]
        elif "|" in marker_name:
            split_marker_name = marker_name.split("|")
            a_marker_name = split_marker_name[0]
            b_marker_name = split_marker_name[1]

        # Chromosome number, starting position, ending position
        for idx, fusion_element in enumerate(var_dict[BREAKENDS]):
            # Various information on gene A
            if a_marker_name == fusion_element[TRANSCRIPTS][0][GENE_SYMBOL]:
                self.chr1 = var_dict[BREAKENDS][idx][CHROMOSOME]
                self.start_pos1 = var_dict[BREAKENDS][idx][START_POSITION]
                self.end_pos1 = var_dict[BREAKENDS][idx][END_POSITION]

            # Various information on gene B
            if b_marker_name == fusion_element[TRANSCRIPTS][0][GENE_SYMBOL]:
                self.chr2 = var_dict[BREAKENDS][idx][CHROMOSOME]
                self.start_pos2 = var_dict[BREAKENDS][idx][START_POSITION]
                self.end_pos2 = var_dict[BREAKENDS][idx][END_POSITION]

        chr_pos_list = [self.chr1, self.start_pos1, self.end_pos1, self.chr2,
                        self.start_pos2, self.end_pos2]

        return chr_pos_list


def get_tsv_file_date():
    """
    Obtains the date to be appended to the name of the TSV file to be saved.
    :return: Date to be appended to the name of the TSV file to be saved
    """
    time_delta = datetime.timedelta(hours=9)
    jst = datetime.timezone(time_delta, 'JST')
    now_time = datetime.datetime.now(jst)
    conv_now_time = now_time.strftime('%Y%m%d%H%M%S')

    return conv_now_time


def create_cats_format_list(input_pattern):
    """
    Get the file path of CATS format from the input pattern and list it.
    :param input_pattern: input CATS format pattern
    :return: List containing the file path of CATS format
    """
    # If the input data pattern includes a ".cats" extension
    if check_extension_cats(input_pattern):
        cats_format_list = glob.glob(input_pattern, recursive=True)
    # If the input data condition does not include ".cats" extension
    else:
        cats_format_list = glob.glob(f"{input_pattern}/**/*{EXT_CATS}", recursive=True)

    return cats_format_list


def check_extension_cats(input_pattern):
    """
    Check if the extension "cats" is included in the pattern specification
    of the CATS format for graph aggregation.
    :param input_pattern: input CATS format pattern
    :return: The search pattern of the input data ends
            with the extension ".cats" or not.
    """
    return input_pattern.endswith(EXT_CATS)


def output_tsv_cats_distribution_result(result_tsv, target_cats_dict):
    """
    Outputs a list of CATS formats to be tabulated in TSV format.
    :param result_tsv: Name of TSV file to output the list of CATS formats
                        to be aggregated
    :param target_cats_dict: Dictionary of CATS formats to be included
                            in the listing output
                            (Key: CATS format path, Value: sample ID)
    """
    with open(result_tsv, mode="w", encoding="utf-8", newline="") as result_writer:
        writer = csv.writer(result_writer, delimiter="\t")
        result_header_list = [ITEM_NUMBER, HEADER_SAMPLE_ID, HEADER_CATS_FILE]
        writer.writerow(result_header_list)

        for idx, (cats_path, cats_id) in enumerate(target_cats_dict.items()):
            item_num = idx + 1
            result_element_list = [item_num, cats_id, cats_path]
            writer.writerow(result_element_list)


def get_panel_name(cats_dict):
    """
    Obtain the panel name from the CATS format.
    :param cats_dict: String read from the CATS format as a dictionary.
    :return: Panel name taken from CATS format.
    """
    return cats_dict[TEST_INFO][PANEL_NAME]


def get_cats_variant_data(cats_dict):
    """
    Obtain various information from CATS format.
    :param cats_dict: String read from the CATS format as a dictionary.
    :return: List containing various information obtained from CATS format.
    """
    # List to store mutation information
    var_msi_list = []
    var_tmb_list = []
    var_snv_list = []
    var_cna_list = []
    var_sv_list = []

    # Obtain various MSI and TMB values from "otherBiomarkers" key
    if OTHER_BIOMARKERS in cats_dict.keys():
        for biomarker in cats_dict[OTHER_BIOMARKERS]:
            # MSI
            if all([biomarker[BIOMARKER_TYPE] == MSI, biomarker[REPORTED]]):
                var_msi_list.append(biomarker)
            # TMB
            if all([biomarker[BIOMARKER_TYPE] == TMB, biomarker[REPORTED]]):
                var_tmb_list.append(biomarker)
    # Obtain various values from the "variants" key
    if VARIANTS in cats_dict.keys():
        for cats_var_type in cats_dict[VARIANTS]:
            # Generates a list containing various types of information
            # (one element contains various types of information
            # for one marker)
            # SNV
            if cats_var_type == VARIANT_SNV:
                var_snv_list = cats_dict[VARIANTS][cats_var_type]
            # CNA
            elif cats_var_type == VARIANT_CNA:
                var_cna_list = cats_dict[VARIANTS][cats_var_type]
            # SV(only fusion)
            elif cats_var_type == VARIANT_SV:
                var_sv_list = cats_dict[VARIANTS][cats_var_type]

    return var_msi_list, var_tmb_list, var_snv_list, var_cna_list, var_sv_list


def create_tsv_data(var_list, sample_data_id, panel_name, var_type):
    """
    Convert mutation information obtained from the CATS format
    into a string in TSV format.
    :param var_list: List of mutation information obtained from CATS format
    :param sample_data_id: sample ID from the CATS format file index
    :param panel_name: Panel name taken from CATS format.
    :param var_type: variant type
    :return: List of mutation information obtained from CATS format
            converted to TSV format strings.
    """
    # For storing return list
    cats_data_list = []

    # MSI
    if all([var_type == MSI, len(var_list) != 0]):
        for var_msi in var_list:
            # Generate instance
            cats_msi = CatsConverter()

            cats_msi.sample_id = sample_data_id
            cats_msi.panel_name = panel_name
            # Biomarker type(MSI)
            cats_msi.variant_type = MSI
            # MSI status
            cats_msi.state = var_msi[STATE]

            if BIOMARKER_METRICS in var_msi.keys():
                # MSI score
                cats_msi.value = str(var_msi[BIOMARKER_METRICS][0][MSI_VALUE])

            # Stores various elements of mutation in a return list
            msi_data = cats_msi.create_variant_element_list()
            cats_data_list.append(msi_data)

    # TMB
    if all([var_type == TMB, len(var_list) != 0]):
        for var_tmb in var_list:
            # Generate instance
            cats_tmb = CatsConverter()

            cats_tmb.sample_id = sample_data_id
            cats_tmb.panel_name = panel_name
            # Biomarker type(TMB)
            cats_tmb.variant_type = TMB
            # TMB score
            cats_tmb.value = str(var_tmb[BIOMARKER_METRICS][0][TMB_VALUE])

            # Stores various elements of mutation in a return list.
            tmb_data = cats_tmb.create_variant_element_list()
            cats_data_list.append(tmb_data)

    # SNV
    if all([var_type == VARIANT_SNV, len(var_list) != 0]):
        for var_snv in var_list:
            # reported = True only
            if not var_snv[REPORTED]:
                continue

            # Generate instance
            cats_snv = CatsConverter()

            cats_snv.sample_id = sample_data_id
            cats_snv.panel_name = panel_name
            cats_snv.variant_type = SNV
            cats_snv.marker_name = var_snv[TRANSCRIPTS][0][GENE_SYMBOL]
            # Somatic variant / Germline variant
            cats_snv.origin = cats_snv.get_variant_origin(var_snv)
            # Amino acid change (CDS change only in some cases)
            if var_snv[TRANSCRIPTS][0][AMINO_ACIDS_CHANGE] is not None:
                cats_snv.amino_change = var_snv[TRANSCRIPTS][0][AMINO_ACIDS_CHANGE]
            cats_snv.cds_change = var_snv[TRANSCRIPTS][0][CDS_CHANGE]
            # DNA sequence before change
            cats_snv.ref_allele = var_snv[REFERENCE_ALLELE]
            # DNA sequence after change
            cats_snv.alt_allele = var_snv[ALTERNATE_ALLELE]
            # Chromosome number
            cats_snv.chr1 = var_snv[CHROMOSOME]
            # Starting position
            cats_snv.start_pos1 = str(var_snv[SNV_POSITION])
            # Allele frequency
            cats_snv.value = str(var_snv[ALLELE_FREQ])

            if VARIANT_TYPE in var_snv.keys():
                # Mutation type(subdivision)
                cats_snv.state = var_snv[VARIANT_TYPE]

            # Stores various elements of mutation in a return list
            snv_data = cats_snv.create_variant_element_list()
            cats_data_list.append(snv_data)

    # CNA
    if all([var_type == VARIANT_CNA, len(var_list) != 0]):
        for var_cna in var_list:
            # reported = True only
            if not var_cna[REPORTED]:
                continue

            # Generate instance
            cats_cna = CatsConverter()

            cats_cna.sample_id = sample_data_id
            cats_cna.panel_name = panel_name
            base_variant_type = var_cna[CNA_TYPE]
            cats_cna.variant_type = CNA_CONVERT_DICT.get(base_variant_type, base_variant_type)
            cats_cna.marker_name = var_cna[TRANSCRIPTS][0][GENE_SYMBOL]
            # Somatic variant / Germline variant
            cats_cna.origin = cats_cna.get_variant_origin(var_cna)

            cna_value_dict = {}
            # Copy number(fold-change > absolute copy number)
            for cna_metrics in var_cna[COPY_NUMBER_METRICS]:
                copy_unit_value = cna_metrics[COPY_UNIT]
                cna_value_dict[copy_unit_value] = str(cna_metrics[COPY_VALUE])

                # Obtain values in the order fold-change ⇒ absolute copy number
                if FOLD_CHANGE in cna_value_dict:
                    cats_cna.value = cna_value_dict[FOLD_CHANGE]
                elif ABSOLUTE_COPY_NUMBER in cna_value_dict:
                    cats_cna.value = cna_value_dict[ABSOLUTE_COPY_NUMBER]

            # Chromosome number, starting position, ending position
            cats_cna.chr1 = var_cna.get(CHROMOSOME, "-")
            cats_cna.start_pos1 = var_cna.get(START_POSITION, "-")
            cats_cna.end_pos1 = var_cna.get(END_POSITION, "-")

            # Stores various elements of mutation in a return list.
            cna_data = cats_cna.create_variant_element_list()
            cats_data_list.append(cna_data)

    # SV(only fusion)
    if all([var_type == VARIANT_SV, len(var_list) != 0]):
        for var_sv in var_list:
            # reported = True and only fusion mutations are targeted
            if not all([var_sv[REPORTED], var_sv[REARRANGEMENT_TYPE] == GENE_FUSION]):
                continue

            # Generate instance
            cats_sv = CatsConverter()

            cats_sv.sample_id = sample_data_id
            cats_sv.panel_name = panel_name
            cats_sv.variant_type = FUSION
            cats_sv.marker_name = cats_sv.get_sv_marker_name(var_sv)
            # Somatic variant / Germline variant
            cats_sv.origin = cats_sv.get_variant_origin(var_sv)

            # get chromosome number and position information
            # for fusion mutations
            # (chr1, start_pos1, end_pos1, chr2, start_pos2, end_pos2)
            cats_sv.get_sv_chromosome_position(var_sv, cats_sv.marker_name)

            # Allele frequency
            if ALLELE_FREQ in var_sv.keys():
                cats_sv.value = str(var_sv[ALLELE_FREQ])

            # Stores various elements of mutation in a return list
            sv_data = cats_sv.create_variant_element_list()
            cats_data_list.append(sv_data)

    return cats_data_list


def create_tsv_file(cats_base_pattern, tsv_dir):
    """
    Reads CATS format and converts it to TSV file for output.
    :param cats_base_pattern: Search pattern for CATS format paths to be used
                            for graph aggregation.
    :param tsv_dir: Path to output TSV file
    :return: Name of output TSV file
    """
    # List of CATS format
    cats_file_list = create_cats_format_list(cats_base_pattern)
    # Define the name of the TSV file to be output
    tsv_date = get_tsv_file_date()
    result_file = os.path.join(tsv_dir, f"cats_input_data_{tsv_date}.tsv")
    # List of items that are TSV headers
    header_list = [Aggr.SAMPLE_ID, Aggr.PANEL_NAME, Aggr.ORIGIN_TYPE,
                   Aggr.TYPE, Aggr.VARIANT_TYPE_STATE, Aggr.MARKER_NAME,
                   Aggr.AMINO_ACIDS_CHANGE, Aggr.CDS_CHANGE,
                   Aggr.REFERENCE_ALLELE, Aggr.ALTERNATE_ALLELE,
                   Aggr.CHROMOSOME1, Aggr.START_POS1, Aggr.END_POS1,
                   Aggr.VALUE_COL, Aggr.CHROMOSOME2, Aggr.START_POS2,
                   Aggr.END_POS2]
    # For counting the number of Key Error
    err_cnt = 0
    # List to store elements of various mutation information
    variant_element_list = []

    # For storing a dictionary of CATS formats to be aggregated
    # (Key: CATS format path, Value: sample ID)
    success_cats_dict = {}
    # For storing a list of CATS formats that are not subject to aggregation
    # (Key: CATS format path, Value: sample ID)
    invalid_cats_dict = {}

    # Obtain various mutation information elements from the CATS format.
    for cats_idx, cats_format in enumerate(cats_file_list):
        # Sample ID
        sample_id = cats_idx + 1
        try:
            with open(cats_format, 'r', encoding='utf-8') as cats_reader:
                cats_file_dict = json.load(cats_reader)
                # Panel name
                panel = get_panel_name(cats_file_dict)
                # Variant list
                cats_variant_list = get_cats_variant_data(cats_file_dict)

            msi_list = cats_variant_list[0]
            tmb_list = cats_variant_list[1]
            snv_list = cats_variant_list[2]
            cna_list = cats_variant_list[3]
            sv_list = cats_variant_list[4]
            # List of data retrieved from CATS format
            msi_variant_list = create_tsv_data(msi_list, sample_id, panel, MSI)
            tmb_variant_list = create_tsv_data(tmb_list, sample_id, panel, TMB)
            snv_variant_list = create_tsv_data(snv_list, sample_id, panel,
                                               VARIANT_SNV)
            cna_variant_list = create_tsv_data(cna_list, sample_id, panel,
                                               VARIANT_CNA)
            sv_variant_list = create_tsv_data(sv_list, sample_id, panel,
                                              VARIANT_SV)
            # Stores elements of various mutation information in a list
            variant_element_list.append(msi_variant_list)
            variant_element_list.append(tmb_variant_list)
            variant_element_list.append(snv_variant_list)
            variant_element_list.append(cna_variant_list)
            variant_element_list.append(sv_variant_list)

            success_cats_dict[cats_format] = sample_id
        except KeyError as e:
            err_cnt += 1
            invalid_cats_dict[cats_format] = sample_id
            print(f"The key in the CATS format was not found correctly."
                  f"(file: {cats_format})")
            print(e)
            continue

    # Write data acquired from CATS format in TSV format
    with open(result_file, "w", encoding="utf-8", newline="") as tsv_writer:
        writer = csv.writer(tsv_writer, delimiter="\t")
        # Write header
        writer.writerow(header_list)
        # Writing of each element retrieved from CATS format
        for variant_element in variant_element_list:
            writer.writerows(variant_element)

    # Outputs a list of CATS-formatted text to be aggregated.
    success_cats_tsv_path = os.path.join(tsv_dir, SUCCESS_CATS_TSV)
    output_tsv_cats_distribution_result(success_cats_tsv_path, success_cats_dict)
    # Output a textual list of CATS formats that are no longer included
    # in the aggregation.
    invalid_cats_tsv_path = os.path.join(tsv_dir, INVALID_CATS_TSV)
    output_tsv_cats_distribution_result(invalid_cats_tsv_path, invalid_cats_dict)

    # Number of CATS formats to be aggregated
    cats_file_length = len(cats_file_list)

    print(f"\nFinished TSV file generation process."
          f"(Success: {cats_file_length - err_cnt}/{cats_file_length}, "
          f"Error: {err_cnt}/{cats_file_length})")

    return result_file
