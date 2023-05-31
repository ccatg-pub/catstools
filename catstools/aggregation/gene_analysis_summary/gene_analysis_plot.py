#!/usr/bin/env python
# coding: utf-8

import os

import matplotlib.pyplot as plt
import numpy as np
from catstools.aggregation.aggregation import Aggregation


class GeneAnalysisPlot(Aggregation):
    # data frame column name(for input data)
    MARKER = 'Marker'
    CHR_NUM = 'chr_num'
    LENGTHS = 'Lengths'
    POSITION = 'position'
    START_POS = 'Start_P'
    START_POSITION = 'Start Position'
    END_POSITION = 'End Position'
    FUS_GENE1 = 'Gene1'
    FUS_GENE2 = 'Gene2'
    FUS_GENE1_GENE2 = 'Gene1-Gene2'
    FUS_GENE_NUM = 'Gene'
    SIZE = 'size'
    FREQ = 'freq'
    RANK = 'rank'
    STRAND_X = 'strand_x'
    STRAND_Y = 'strand_y'
    GENE1_STRAND = 'Gene1_strand'
    GENE2_STRAND = 'Gene2_strand'
    START_AND_END = 'SandE'
    START_SIDE = 'Start'
    END_SIDE = 'End'
    MARKER_NAME_NUMBER = 'Marker_Name_No'
    MUT_TYPE_NUMBER = 'Mut_type_No'
    MARKER_NUMBER = 'Marker_No'
    STRAND_NUMBER = 'strand_No'
    CHROMOSOME_NUMBER = 'Chromosome_No'
    POSITION_NUMBER = 'position_No'
    FREQ_NUMBER = 'freq_No'
    RANK_NUMBER = 'rank_No'
    GENE_NUMBER = 'Gene_No'
    START_AND_END_NUMBER = 'SandE_No'
    CHR_NUM_NUMBER = 'chr_num_No'

    # data frame column name(for strand)
    STRAND_GENE = 'gene'
    STRAND_DIRECTION = 'strand'
    ALIASES = 'aliases'

    OFFSET_POINTS = 'offset points'
    DENSE = 'dense'
    ROUND = 'round'
    ORANGE_RED = 'orangered'
    DARKBLUE = 'darkblue'
    BEST = 'best'
    MAJOR = 'major'
    FFILL = 'ffill'

    # Capital letter
    AMINO_ACIDS = 'AminoAcids'
    PROMOTER = 'promoter'
    SPLICE_SITE = 'splice site'
    AMINO_ACIDS_NUM = 'Amino_Acids_Num'
    PROTEIN_ID = 'ProteinID'
    PROTEIN_ID_No = 'ProteinID_No'
    AMINO_ACIDS_NUM_NO = 'Amino_Acids_Num_No'
    COLUMNS = 'columns'
    MINOR = 'minor'
    AMPLIFICATION_START_POSITION = 'Amplification_Start Position'
    AMPLIFICATION_END_POSITION = 'Amplification_End Position'
    LOSS_START_POSITION = 'Loss_Start Position'
    LOSS_END_POSITION = 'Loss_End Position'
    LENGTH = 'Length'
    TRANSCRIPT_ID_VERSION = 'TranscriptIDVersion'
    TRANSCRIPT_ID_VERSION_NUMBER = 'TranscriptIDVersion_No'
    START1_END1 = 'Start1-End1'
    START_P_PDIFF = 'start_P_pdiff'
    END_P_PDIFF = 'end_P_pdiff'
    END1_PDIFF = 'End1_pdiff'
    START1_PDIFF = 'Start1_pdiff'
    NONE = 'none'

    PROTEIN_ID_NUMBER = 'ProteinID_No'
    SIZE_NUMBER = 'size_No'
    AMINO_ACIDS_NUM_NUMBER = 'Amino_Acids_Num_No'

    CDS_CHANGE_NUM = 'Cds_Change_Num'
    P_DIFF = 'p_diff'
    P_DIFF_NO = 'p_diff_No'
    SIZE_NO = 'size_No'

    FONT_SIZE = 'font.size'
    LEGEND_FONTSIZE = 'legend.fontsize'
    LEGEND_HANDLELENGTH = 'legend.handlelength'
    OUT_EXTENSION = '.png'

    NOT_TYPE_FUSION = 'not Type=="Fusion"'
    MUT_TYPE_EQ = 'Mut_type=='
    REG_NUMBER_REP = r'(\d+)'
    NOT_CHROMOSOME1 = 'not Chromosome1=="-"'
    NOT_START1 = 'not Start1=="-"'
    NOT_END1 = 'not End1=="-"'

    # constructor
    def __init__(self):
        pass

    def get_input_tid(self) -> str:
        # get main.py current path
        main_py_path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(main_py_path, "input_data", 'transcriptionID.tsv')

    def get_input_strand(self) -> str:
        # get main.py current path
        main_py_path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(main_py_path, "input_data", 'gene_strand_transcript.tsv')

    def get_input_aliases(self) -> str:
        # get main.py current path
        main_py_path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(main_py_path, "input_data", 'gene_aliases.tsv')

    def get_input_chromosome(self) -> str:
        # get main.py current path
        main_py_path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(main_py_path, "input_data", 'chromosome_length_GRCh37.tsv')

    def draw_graphical_stem(self, mut_type_df, mut_type, mut_type_list):
        """
        Draw graphical stem(Plotting values)
        """
        if len(mut_type_df.freq) == 0:
            return None

        # plot values in graph
        if mut_type == self.SNV:
            plt.stem(mut_type_df.position, mut_type_df.freq)
        elif mut_type in mut_type_list:
            # Determine the color that will serve as the legend
            color_number = str(mut_type_list.index(mut_type) + 1)
            line_color = f"C{color_number}-"
            marker_color = f"C{color_number}o"

            plt.stem(mut_type_df.position, mut_type_df.freq, linefmt=line_color, markerfmt=marker_color)
        else:
            raise ValueError(f"Error: mut_type invalid. val={mut_type}")

    def draw_graphical_framework(self, max_freq):
        """
        Draw graphical framework(Set scale to Y-axis)
        """
        if max_freq == 0:
            return None
        elif 0 < max_freq < 6:
            plt.yticks(np.arange(0, max_freq, 1))
            plt.grid(which=self.MAJOR, axis="y", alpha=0.8, linestyle="--", linewidth=1)
        else:
            plt.yticks(np.arange(0, max_freq, 5))
            plt.minorticks_on()
            plt.grid(color=self.BLACK, axis="y", linestyle="--")
            plt.grid(which=self.MINOR, axis="y", alpha=0.8, linestyle="--", linewidth=1)

    def annotate_gene_graph(self, mut_type_df):
        """
        annotate text to gene graph plot(except fusion)
        """
        if len(mut_type_df.freq) == 0:
            return None

        for pos, freq in zip(mut_type_df[self.POSITION], mut_type_df[self.FREQ]):
            plt.annotate('{:.0f}'.format(pos), xy=(pos, freq), xytext=(0, 5),
                         textcoords=self.OFFSET_POINTS, ha=self.CENTER, fontsize=18)

    def annotate_gene_graphs(self, mut_type_df_list):
        for mut_type_df in mut_type_df_list:
            self.annotate_gene_graph(mut_type_df)

    def calculate_mode_frequency(self, *calc_val):
        """
        Calculate the mode frequency in graph
        """
        return max(calc_val)

    def draw_mutation_type_legend(self, mut_type_list, mut_type_max):
        """
        Determine the mutation type to be used in the legend
        """
        if 1 <= len(mut_type_list) <= mut_type_max:
            plt.legend(mut_type_list, fontsize=24)
        else:
            return None

    def calculate_frequency(self, mut_type_df):
        """
        Calculate the mode frequency for each mutation type
        """
        mut_type_freq = mut_type_df.freq
        if len(mut_type_df) == 0:
            return 0
        else:
            return int(max(mut_type_freq)) + 1

    def calculate_frequencies(self, mut_type_df_list):
        mcalc_list = []
        for mut_type_df in mut_type_df_list:
            mcalc_list.append(self.calculate_frequency(mut_type_df))
        return max(mcalc_list)

    def sort_mutation_type(self, mut_type, mut_type_list):
        """
        Sorting function of Mut_type
        """
        if mut_type in mut_type_list:
            return mut_type_list.index(mut_type) + 1
        else:
            return len(mut_type_list) + 1

    def create_data_frame(self, data_df, mut_type, rename_col):
        """
        Create data frames for each mutation type.
        """
        mut_type_df = data_df.query(f'{self.MUT_TYPE_EQ}"{mut_type}"')
        mut_type_df = mut_type_df.rename(columns={rename_col: self.POSITION, self.SIZE: self.FREQ})
        mut_type_df = mut_type_df.loc[:, [self.POSITION, self.FREQ]]

        return mut_type_df

    def get_last_loc(self, data_flame, key):
        return data_flame.loc[:, key].drop_duplicates().iloc[-1]

    def get_mut_type_pre(self, data_flame):
        mut_type_pre = data_flame.loc[:, self.MUT_TYPE].drop_duplicates().sort_index(ascending=False)
        mut_type_pre = mut_type_pre.reset_index().loc[:, self.MUT_TYPE]
        return mut_type_pre
