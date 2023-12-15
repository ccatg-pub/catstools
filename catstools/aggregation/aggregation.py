#!/usr/bin/env python
# coding: utf-8

import os
import warnings

import matplotlib as mpl
import pandas as pd


class Aggregation():

    # data frame column name(for input data)
    SAMPLE_ID = 'Sample_id'
    PANEL_NAME = 'Panel_name'
    ORIGIN_TYPE = 'Origin_Type'
    AMINO_ACIDS_CHANGE = 'Amino_Acids_Change'
    CDS_CHANGE = 'Cds_Change'
    VALUE_COL = 'Value'
    CHROMOSOME1 = 'Chromosome1'
    CHROMOSOME2 = 'Chromosome2'
    STANDARD_GENE = 'standard_gene'
    ALIASES = 'aliases'

    # Capital letter
    REFERENCE_ALLELE_1STL = 'ReferenceAllele_1stL'
    REFERENCE_ALLELE = 'ReferenceAllele'
    REFERENCE_ALLELE_LEN = 'ReferenceAllele_len'
    ALTERNATE_ALLELE_LEN = 'AlternateAllele_len'
    ALTERNATE_ALLELE = 'AlternateAllele'
    ALTERNATE_ALLELE_1STL = 'AlternateAllele_1stL'
    RA_TF = 'R/A_TF'
    AMPLIFICATION = 'Amplification'
    LOSS = 'Loss'
    FUSION = 'Fusion'

    # various option value
    LEFT = 'left'
    RIGHT = 'right'
    CENTER = 'center'
    BOTTOM = 'bottom'
    RED = 'red'
    BLUE = 'blue'
    MEDIUM_BLUE = 'mediumblue'
    BLACK = 'black'
    ORANGE = 'orange'
    START_POS1 = 'Start1'
    START_POS2 = 'Start2'
    END_POS1 = 'End1'
    END_POS2 = 'End2'
    MUT_TYPE = 'Mut_type'
    TIGHT = 'tight'
    VARIANT_TYPE_STATE = 'Variant_Type/State'
    TYPE = 'Type'
    AMP = 'Amp'
    SNV = 'SNV'
    SNP = 'SNP'
    CNV = 'CNV'
    LINK = 'LINK'
    INSERTION = 'Insertion'
    DELETION = 'Deletion'
    DELINS = 'Delins'
    FREQUENCY = 'Frequency'
    GAINSBORO = 'gainsboro'
    MARKER_NAME = 'Marker_Name'
    CHROMOSOME = 'Chromosome'
    START_P = 'start_P'
    END_P = 'end_P'
    SEP_TAB = '\t'
    CONNECT_OR = ' or '
    CONNECT_AND = ' and '
    ORIGIN_SOMATIC = 'Origin_Type=="somatic"'
    NOT_TYPE_MSI = 'not Type=="MSI"'
    NOT_TYPE_TMB = 'not Type=="TMB"'
    TYPE_MSI = 'Type=="MSI"'
    TYPE_TMB = 'Type=="TMB"'
    TYPE_SNV = 'Type=="SNV"'
    TYPE_AMP = 'Type=="Amp"'
    TYPE_LOSS = 'Type=="Loss"'
    TYPE_FUSION = 'Type=="Fusion"'

    SEPALATOR_LINE = '-----------------------------------'
    NO_DATA = '-'
    GENE_ANALYSIS = 'gene_analysis_summary'
    INPUT_DATA_DIR = 'input_data'

    # determine backends of matplotlib
    mpl.use('Agg')
    warnings.filterwarnings('ignore')
    warnings.simplefilter('ignore')

    def __init__(self):
        # constructor
        pass

    def get_input_dir_path(self, file_name) -> str:
        return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            self.GENE_ANALYSIS, self.INPUT_DATA_DIR, file_name)

    def get_input_tid(self) -> str:
        return self.get_input_dir_path('transcriptionID.tsv')

    def get_input_strand(self) -> str:
        return self.get_input_dir_path('gene_strand_transcript.tsv')

    def get_input_aliases(self) -> str:
        return self.get_input_dir_path('gene_aliases.tsv')

    def get_input_chromosome(self) -> str:
        return self.get_input_dir_path('chromosome_length_GRCh37.tsv')

    def set_allele_info(self, in_df):
        df_type = in_df
        df_type[self.REFERENCE_ALLELE_LEN] = list(map(len, df_type[self.REFERENCE_ALLELE]))
        df_type[self.ALTERNATE_ALLELE_LEN] = list(map(len, df_type[self.ALTERNATE_ALLELE]))
        df_type[self.REFERENCE_ALLELE_1STL] = df_type[self.REFERENCE_ALLELE].str[:1]
        df_type[self.ALTERNATE_ALLELE_1STL] = df_type[self.ALTERNATE_ALLELE].str[:1]
        df_type[self.RA_TF] = (df_type[self.REFERENCE_ALLELE_1STL] == df_type[self.ALTERNATE_ALLELE_1STL]) * 1
        df_type[self.MUT_TYPE] = df_type.apply(lambda c: self.assign_mutation_type(c), axis=1)
        return df_type

    def assign_mutation_type(self, gene_df):
        """
        Assign mutation type based on the information
        of the reference allele and the alternate allele
        """
        variant_type = gene_df[self.TYPE]
        ref_len = gene_df[self.REFERENCE_ALLELE_LEN]
        alt_len = gene_df[self.ALTERNATE_ALLELE_LEN]
        first_ref_alt_match = gene_df[self.RA_TF]

        # Amplification
        if variant_type == self.AMP:
            mut_type = self.AMPLIFICATION
        # Loss
        elif variant_type == self.LOSS:
            mut_type = self.LOSS
        # Fusion
        elif variant_type == self.FUSION:
            mut_type = self.FUSION
        # short-variant
        else:
            # SNV
            if all([ref_len == 1, alt_len == 1, first_ref_alt_match == 0]):
                mut_type = self.SNV
            # Deletion
            elif all([ref_len >= 2, alt_len == 1, first_ref_alt_match == 1]):
                mut_type = self.DELETION
            # Insertion
            elif all([ref_len == 1, alt_len >= 2, first_ref_alt_match == 1]):
                mut_type = self.INSERTION
            # Delins
            else:
                mut_type = self.DELINS

        return mut_type

    def exclude_non_target_gene(self, gene_df, trans_df):
        """
        If an unexpected gene is included, it is excluded from the data frame
        and the unexpected gene is output to the console.
        """
        std_gene_df = self.create_standard_gene_df(gene_df)
        # For storing a list of unexpected genes
        no_marker_list = []

        # If no value exists, set the Chromosome, Start Position,
        # and End Position information that each gene has.
        for gene_idx in range(std_gene_df.shape[0]):
            if std_gene_df[self.TYPE].iat[gene_idx] == self.FUSION:
                continue

            target_marker_name = std_gene_df[self.MARKER_NAME].iat[gene_idx]
            # Obtain information on the gene of the transcript, if available.
            trans_info = trans_df[trans_df[self.MARKER_NAME] == target_marker_name]

            # If it is an unexpected gene, store it in the list
            # of unexpected genes
            if len(trans_info) == 0:
                no_marker_list.append(target_marker_name)

        # Output the names of markers that are not included in the aggregation
        # to the console.
        if len(no_marker_list) != 0:
            print('List of markers not included in the aggregation.\n'
                  '(because there is no information on the gene for the transcript.)')
            print(self.SEPALATOR_LINE)

            # Exclude gene duplications
            no_marker_list = sorted(list(dict.fromkeys(no_marker_list)))

            for no_marker in no_marker_list:
                print(f'{no_marker}')

            print(f'{self.SEPALATOR_LINE}\n')

        target_gene_df = std_gene_df.query(f"{self.MARKER_NAME} not in {no_marker_list}")

        return target_gene_df

    def set_gene_default_position(self, gene_df, trans_df):
        """
        For mutations that do not have position information,
        set the position information that is defined for each gene(CNA only).
        """
        # Maintains original data and generates data frames for return values
        conv_gene_df = gene_df.copy()

        # If no value exists, set the Chromosome, Start Position,
        # and End Position information that each gene has.
        for gene_idx in range(conv_gene_df.shape[0]):
            # Chromosome and position information is not completed
            # unless it is a CNA mutation.
            if not any([conv_gene_df[self.TYPE].iat[gene_idx] == self.AMP,
                        conv_gene_df[self.TYPE].iat[gene_idx] == self.LOSS]):
                continue

            target_marker_name = conv_gene_df[self.MARKER_NAME].iat[gene_idx]
            # Obtain chromosome and position as defined by the target gene
            gene_chr_pos_info = trans_df[trans_df[self.MARKER_NAME] == target_marker_name]

            # If the gene is not included in the aggregation,
            # proceed to the next gene.
            if any([conv_gene_df[self.CHROMOSOME1].iat[gene_idx] == self.NO_DATA,
                    conv_gene_df[self.START_POS1].iat[gene_idx] == self.NO_DATA,
                    conv_gene_df[self.END_POS1].iat[gene_idx] == self.NO_DATA]):
                # Chromosome1
                if conv_gene_df[self.CHROMOSOME1].iat[gene_idx] == self.NO_DATA:
                    conv_gene_df[self.CHROMOSOME1].iat[gene_idx] = gene_chr_pos_info[self.CHROMOSOME].iat[0]
                # Start1
                if conv_gene_df[self.START_POS1].iat[gene_idx] == self.NO_DATA:
                    conv_gene_df[self.START_POS1].iat[gene_idx] = gene_chr_pos_info[self.START_P].iat[0]
                # End1
                if conv_gene_df[self.END_POS1].iat[gene_idx] == self.NO_DATA:
                    conv_gene_df[self.END_POS1].iat[gene_idx] = gene_chr_pos_info[self.END_P].iat[0]

        return conv_gene_df

    def create_standard_gene_df(self, gene_df):
        """
        Generate a data frame that converts possible gene readings to standard gene names.
        """
        std_gene_df = gene_df.copy()
        gene_aliases = pd.read_csv(self.get_input_aliases(), sep=self.SEP_TAB)

        for gene_df_idx in range(gene_df.shape[0]):
            marker = gene_df[self.MARKER_NAME].iat[gene_df_idx]
            aliases_marker = gene_aliases[gene_aliases[self.ALIASES] == marker]
            aliases_marker_len = len(aliases_marker)

            # If multiple gene synonyms are present, output to log.
            if aliases_marker_len >= 2:
                print(f'{len(aliases_marker)} synonyms for {marker} exist and cannot be converted to '
                      f'standard gene names.')

            # If a gene synonym exists and can be converted to a standard gene name, it is converted.
            elif aliases_marker_len == 1:
                std_marker = aliases_marker[self.STANDARD_GENE].iat[0]
                std_gene_df[self.MARKER_NAME].iat[gene_df_idx] = std_marker

        return std_gene_df
