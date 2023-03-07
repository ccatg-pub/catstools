#!/usr/bin/env python
# coding: utf-8
import gc
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .gene_analysis_plot import GeneAnalysisPlot


class GeneAnalysisPlotDnaFusion(GeneAnalysisPlot):

    TYPE_FUSION = 'Type=="Fusion"'
    NOT_NO_VALUE_GENE = 'not Marker_Name.str.contains("None")'
    GENE1_START = 'Gene==1 & SandE=="Start"'
    GENE1_END = 'Gene==1 & SandE=="End"'
    GENE1_SANDE = 'Gene==1 & SandE=="SandE"'
    GENE2_START = 'Gene==2 & SandE=="Start"'
    GENE2_END = 'Gene==2 & SandE=="End"'
    GENE2_SANDE = 'Gene==2 & SandE=="SandE"'
    GENE3_START = 'Gene==3 & SandE=="Start"'
    GENE3_END = 'Gene==3 & SandE=="End"'
    GENE3_SANDE = 'Gene==3 & SandE=="SandE"'
    OUT_FILENAME = 'Fusion_geneplot_'

    def draw_gene_analysis_fusion(self, input_rawdata, output_for_gene_analysis):
        # PANEL
        base_panel_df = pd.read_csv(input_rawdata, sep=self.SEP_TAB)

        # Check if the number of data frame length is not zero
        if len(base_panel_df) == 0:
            print('No variant data in input data')
            return None

        panel_df = base_panel_df.query(self.CONNECT_AND.join([self.TYPE_FUSION, self.NOT_NO_VALUE_GENE]))

        # Check if the number of fusion data frame length is not zero
        if len(panel_df) == 0:
            print('No fusion variant data in input data')
            return None

        fus_plt = self.get_fus_plt(panel_df)

        # Specify the direction of gene
        gene_strand = pd.read_csv(self.get_input_strand(), sep=self.SEP_TAB)
        gene_fusion_plot = self.gene_fusion_plot(fus_plt, gene_strand)
        marker_name_num = gene_fusion_plot[self.MARKER_NAME_NUMBER]

        for index in range(1, (int(max(gene_fusion_plot[self.MARKER_NAME_NUMBER])) + 1)):
            self.create_graph(gene_fusion_plot[marker_name_num == index], output_for_gene_analysis)

    def create_graph(self, gene_fusion_graph, output_for_gene_analysis):
        # Drawing graphs
        # Draw the entire graph area
        fig = plt.figure(figsize=(20, 15), dpi=400)
        plt.rcParams[self.FONT_SIZE] = 14
        plt.rcParams[self.LEGEND_FONTSIZE] = 16
        plt.rcParams[self.LEGEND_HANDLELENGTH] = 3

        # Determine name to Y-axis
        plt.ylabel(self.FREQUENCY, fontsize=30)
        plt.tick_params(labelsize=30)

        # Gene: Start
        fus_gene_plt_s_list = [
            self.create_data_frame(gene_fusion_graph, self.GENE1_START),
            self.create_data_frame(gene_fusion_graph, self.GENE2_START)
        ]
        # Gene : End
        fus_gene_plt_e_list = [
            self.create_data_frame(gene_fusion_graph, self.GENE1_END),
            self.create_data_frame(gene_fusion_graph, self.GENE2_END)
        ]
        # Gene : Start-END
        fus_gene_plt_se_list = [
            self.create_data_frame(gene_fusion_graph, self.GENE1_SANDE),
            self.create_data_frame(gene_fusion_graph, self.GENE2_SANDE)
        ]

        if max(gene_fusion_graph[self.FUS_GENE_NUM]) == 3:
            ax_list = [plt.subplot2grid((3, 3), (0, 0)), plt.subplot2grid((3, 3), (0, 1)),
                       plt.subplot2grid((3, 3), (0, 2))]
            # Gene3: Start
            fus_gene_plt_s_list.append(self.create_data_frame(gene_fusion_graph, self.GENE3_START))
            # Gene3 : End
            fus_gene_plt_e_list.append(self.create_data_frame(gene_fusion_graph, self.GENE3_END))
            # Gene3 : Start-END
            fus_gene_plt_se_list.append(self.create_data_frame(gene_fusion_graph, self.GENE3_SANDE))
        else:
            ax_list = [plt.subplot2grid((2, 2), (0, 0)), plt.subplot2grid((2, 2), (0, 1))]

        # Display Gene
        self.display_gene_graph(ax_list, fus_gene_plt_s_list, fus_gene_plt_e_list,
                                fus_gene_plt_se_list, output_for_gene_analysis)

        # Open memory
        fig.clf()
        # Close a saved graph
        plt.close(fig)

        # Delete figure instance
        del fig
        # Open memory by garbage collection
        gc.collect()

    def create_data_frame_position(self, fus_df, start_end, gene_num):
        """
        Create data frames based on gene (start/end) position information
        """
        # Determining the column name to use,
        # depending on gene number and start_end_side
        if gene_num == 1:
            gene_col, strand_col, chr_col = self.FUS_GENE1, self.GENE1_STRAND, self.CHROMOSOME1
            # start side
            if start_end == self.START_SIDE:
                pos_col = self.START_POS1
            # end side
            elif start_end == self.END_SIDE:
                pos_col = self.END_POS1
            else:
                print("It is not possible to identify "
                      "whether it is on the start or end side.")
                return None
        elif gene_num == 2:
            gene_col, strand_col, chr_col = self.FUS_GENE2, self.GENE2_STRAND, self.CHROMOSOME2
            # start side
            if start_end == self.START_SIDE:
                pos_col = self.START_POS2
            # end side
            elif start_end == self.END_SIDE:
                pos_col = self.END_POS2
            else:
                print("It is not possible to identify "
                      "whether it is on the start or end side.")
                return None
        else:
            print("Can't assign a fusion gene number to one or two")
            return None

        conv_fus_df = fus_df.loc[:, [self.MARKER_NAME, self.MUT_TYPE,
                                     gene_col, strand_col, chr_col, pos_col, self.SIZE, self.RANK]]
        conv_fus_df = conv_fus_df.rename(
            columns={pos_col: self.POSITION, chr_col: self.CHROMOSOME, gene_col: self.MARKER,
                     strand_col: self.STRAND_DIRECTION})
        conv_fus_df[self.FUS_GENE_NUM] = gene_num
        conv_fus_df[self.START_AND_END] = start_end

        return conv_fus_df

    def create_data_frame_chromosome(self, fus_df, gene_chr_df, start_end, gene_num):
        """
        Create data frames based on chromosome length information
        """
        # Determining the column name to use, depending on gene number
        if gene_num == 1:
            gene_col, strand_col, chr_col = self.FUS_GENE1, self.GENE1_STRAND, self.CHROMOSOME1
        elif gene_num == 2:
            gene_col, strand_col, chr_col = self.FUS_GENE2, self.GENE2_STRAND, self.CHROMOSOME2
        else:
            print("Can't assign a fusion gene number to one or two")
            return None

        # Create the base data frame, depending on start_end_side
        if start_end == self.START_SIDE:
            conv_fus_df = fus_df.loc[:, [self.MARKER_NAME, self.MUT_TYPE,
                                         gene_col, strand_col, chr_col, self.START_POS]]
            conv_fus_df = conv_fus_df.rename(
                columns={self.START_POS: self.POSITION, chr_col: self.CHROMOSOME, gene_col: self.MARKER,
                         strand_col: self.STRAND_DIRECTION})
        elif start_end == self.END_SIDE:
            conv_fus_df = fus_df.loc[:, [self.MARKER_NAME, self.MUT_TYPE, gene_col, strand_col, chr_col]]
            conv_fus_df = conv_fus_df.rename(
                columns={chr_col: self.CHROMOSOME, gene_col: self.MARKER, strand_col: self.STRAND_DIRECTION})
        else:
            print("It is not possible to identify whether it is on the start or end side.")
            return None

        conv_fus_df = conv_fus_df.drop_duplicates()
        conv_fus_df[self.START_AND_END] = self.START_AND_END
        conv_fus_df[self.FUS_GENE_NUM] = gene_num
        conv_fus_df[self.SIZE] = 0

        # Obtain the chromosome length and create data frame.
        if start_end == self.END_SIDE:
            conv_fus_df[self.CHROMOSOME] = conv_fus_df.apply(lambda x: self.cast_chromosome(x.Chromosome), axis=1)
            conv_fus_df = conv_fus_df.astype({self.CHROMOSOME: str})
            conv_fus_df = pd.merge(conv_fus_df, gene_chr_df, how=self.LEFT, on=[self.CHROMOSOME])
            conv_fus_df = conv_fus_df.rename(columns={self.LENGTHS: self.POSITION})

        return conv_fus_df

    def cast_chromosome(self, x):
        """
        If chromosome is integer cast to stringer
        """
        # when the chromosome number is an integer
        if type(x) == int:
            return str(x)
        else:
            return x

    def combine_data_frame(self, fus_start_df, chr_start_df, chr_end_df, fus_end_df):
        """
        UNION-Combine start location and the number by base
        """
        fus_df = fus_start_df.append(chr_start_df, ignore_index=True)
        fus_df = fus_df.append(chr_end_df, ignore_index=True)
        fus_df = fus_df.append(fus_end_df, ignore_index=True)

        return fus_df

    def convert_gene_number(self, x, y):
        """
        If the chromosome of the target gene is not in the correct position
        and is in a different position than it should be,
        designate it as a third gene
        """
        # when the chromosome of the target gene is not in the correct position
        # and is in a different position than it should be
        if x == 2:
            return 3
        else:
            return y

    def add_number(self, fus_df):
        """
        Add the values corresponding to keys chr_num
        and fus_gene_num to the data frame.
        """
        fus_df[self.CHROMOSOME] = fus_df.apply(lambda x: self.cast_chromosome(x.Chromosome), axis=1)
        fus_df_chr_num = fus_df.groupby([self.MARKER_NAME, self.MARKER])[self.CHROMOSOME]
        fus_df[self.CHR_NUM] = fus_df_chr_num.rank(method=self.DENSE, ascending=True)

        fus_df[self.FUS_GENE_NUM] = fus_df.apply(lambda x: self.convert_gene_number(x.chr_num, x.Gene), axis=1)

        return fus_df

    def create_data_frame(self, fus_df, search_terms):
        """
        create input data frame for graph plotting
        """
        fus_query_df = fus_df.query(search_terms)

        # Check if the data frame matching the search terms is not zero
        if len(fus_query_df) == 0:
            print(f"No match data frame(terms: {search_terms})")
            return None

        graph_depiction_fus_df = fus_query_df.loc[:, [self.MARKER_NAME, self.MARKER, self.STRAND_DIRECTION, self.RANK,
                                                      self.CHROMOSOME, self.POSITION, self.FREQ]]

        return graph_depiction_fus_df

    def display_gene_graph(self, ax_list, plt_s_list, plt_e_list, plt_se_list, output_dir):
        """
        Display gene
        """
        for gene_num, (ax, plt_start, plt_end, plt_start_and_end) in enumerate(zip(ax_list, plt_s_list, plt_e_list,
                                                                                   plt_se_list), start=1):
            # Determine gene name
            marker_name = plt_start.loc[:, self.MARKER_NAME].drop_duplicates().iloc[-1]
            marker = plt_start.loc[:, self.MARKER].drop_duplicates().iloc[-1]
            chromosome = plt_start.loc[:, self.CHROMOSOME].drop_duplicates().iloc[-1]
            strand = plt_start.loc[:, self.STRAND_DIRECTION].drop_duplicates().iloc[-1]

            # Display gene name
            print(f"fusion{gene_num}_{marker}")

            # Draw Gene
            self.draw_gene(gene_num, ax, plt_start, plt_end, plt_start_and_end, marker_name, marker,
                           chromosome, strand)

            if gene_num == 2 or gene_num == 3:
                # Specify the save location
                savefig_pass = os.path.join(output_dir, f'{self.OUT_FILENAME}{marker_name}{self.OUT_EXTENSION}')
                plt.savefig(savefig_pass, bbox_inches=self.TIGHT, dpi=400)

    def draw_gene(self, gene_num, ax, plt_start, plt_end, plt_start_and_end, marker_name, marker, chromosome, strand):
        if gene_num == 1 or gene_num == 3:
            gene_start_values = [[-7, 12, 107, 107], [0, 13, 106, 89]]
            gene_end_values = [[5, 12, 107, 90], [6, 11, 108, 91]]
        elif gene_num == 2:
            gene_start_values = [[-5, 13, 106, 106], [-5, 12, 107, 107]]
            gene_end_values = [[7, 11, 108, 106], [5, 12, 107, 107]]

        # Draw Gene-START
        ax.stem(plt_start.position, plt_start.freq)
        for pos, freq, rank in zip(plt_start[self.POSITION], plt_start[self.FREQ], plt_start[self.RANK]):
            text_dict = None
            # Determine the coordinates of the text
            x_coordinate, y_coordinate = self.get_coordinate(strand, rank, gene_start_values)
            if gene_num == 3 and strand == '+':
                text_dict = dict(boxstyle=self.ROUND, fc=self.GAINSBORO, ec=self.MEDIUM_BLUE)
            # Determine transcription color
            ax.annotate('{:.0f}'.format(pos), xy=(pos, freq), xytext=(x_coordinate, y_coordinate), bbox=text_dict,
                        textcoords=self.OFFSET_POINTS, ha=self.RIGHT, va=self.BOTTOM, fontsize=10, color=self.BLUE)

        # Draw Gene-END
        markerline = ax.stem(plt_end.position, plt_end.freq, linefmt='C1-', markerfmt='C1o')[0]
        markerline.set_markerfacecolor(self.NONE)
        for pos, freq, rank in zip(plt_end[self.POSITION], plt_end[self.FREQ], plt_end[self.RANK]):
            text_dict = None
            # Determine the coordinates of the text
            x_coordinate, y_coordinate = self.get_coordinate(strand, rank, gene_end_values)
            if gene_num == 3:
                text_dict = dict(boxstyle=self.ROUND, fc=self.GAINSBORO, ec=self.ORANGE_RED)
            # Determine transcription color
            ax.annotate('{:.0f}'.format(pos), xy=(pos, freq), xytext=(x_coordinate, y_coordinate), bbox=text_dict,
                        textcoords=self.OFFSET_POINTS, ha=self.LEFT, va=self.BOTTOM, fontsize=10, color=self.ORANGE_RED)

        # Draw Gene-START&END
        ax.stem(plt_start_and_end.position, plt_start_and_end.freq, markerfmt='ko')
        for pos, freq in zip(plt_start_and_end[self.POSITION], plt_start_and_end[self.FREQ]):
            ax.annotate('{:.0f}'.format(pos), xy=(pos, freq), xytext=(0, 5),
                        textcoords=self.OFFSET_POINTS, ha=self.CENTER, fontsize=10)

        # Determine name to X-axis of axi
        ax.set_xlabel(f'Chromosome{chromosome} \n $\\it{marker}$ \n $\\it{marker_name}$ Fusion', fontsize=18)
        # Determine name to Y-axis of axi
        ax.set_ylabel(self.FREQUENCY, fontsize=18)
        y_max_g = int(max(plt_start.freq)) + 2

        ax.set_yticks(np.arange(0, y_max_g, 1))
        ax.grid(which=self.MAJOR, axis='y', alpha=0.8, linestyle='--', linewidth=1)
        ax.legend(fontsize=10, loc=self.BEST, labels=[self.START_POSITION, self.END_POSITION],
                  title=f'$\\it{marker}$ breakpoint', labelcolor=[self.BLUE, self.ORANGE_RED])

    def get_coordinate(self, strand, rank, gene_values):
        value_list = gene_values[0] if strand == '+' else gene_values[1]
        x_coordinate = value_list[0]
        # Determine y-coordinate
        if 1 <= rank <= 6:
            y_coordinate = (rank * 17) - value_list[1]
        elif 7 <= rank <= 10:
            y_coordinate = (rank * -17) + value_list[2]
        else:
            y_coordinate = value_list[3]
        return x_coordinate, y_coordinate

    def get_fus_plt(self, panel_df):
        target_fus_plt = panel_df.loc[:, [self.SAMPLE_ID, self.TYPE, self.VARIANT_TYPE_STATE,
                                          self.AMINO_ACIDS_CHANGE, self.CDS_CHANGE,
                                          self.MARKER_NAME, self.CHROMOSOME1, self.START_POS1,
                                          self.END_POS1, self.CHROMOSOME2, self.START_POS2,
                                          self.END_POS2]]
        fus_plt_pre = target_fus_plt.rename(columns={self.TYPE: self.MUT_TYPE})

        # Fusion
        # Case of baseDNA
        fus_plt_vb = fus_plt_pre[fus_plt_pre[self.MARKER_NAME].str.contains(r'\|')].copy()
        fus_plt_sl = fus_plt_pre[fus_plt_pre[self.MARKER_NAME].str.contains(r'\-')].copy()

        # Branching in case of fusion vertical bar
        if len(fus_plt_vb) != 0:
            fus_vb_marker_obj = fus_plt_vb[self.MARKER_NAME].str
            # fusion Gene1
            fus_plt_vb[self.FUS_GENE1] = fus_vb_marker_obj.extract(r'(.+)\|', expand=True)
            # fusion Gene2
            fus_vb_gene2_base = fus_vb_marker_obj.extract(r'\|(.+)', expand=True)
            fus_plt_vb[self.FUS_GENE2] = fus_vb_gene2_base.replace(['fusion', ' '], '', regex=True)
            # fusion Marker Name
            conv_fus_vb_marker_name = fus_vb_marker_obj.replace('|', '-')
            fus_plt_vb[self.MARKER_NAME] = conv_fus_vb_marker_name.str.replace(' fusion', '')

        if len(fus_plt_sl) != 0:
            fus_sl_marker_obj = fus_plt_sl.Marker_Name.str
            # fusion Gene1
            fus_plt_sl[self.FUS_GENE1] = fus_sl_marker_obj.extract(r'(.+)-', expand=True)
            # fusion Gene2
            fus_plt_sl[self.FUS_GENE2] = fus_sl_marker_obj.extract(r'-(.+)', expand=True)

        # UNION
        fus_plt = fus_plt_vb.append(fus_plt_sl, ignore_index=True)
        return fus_plt

    def gene_fusion_plot(self, fus_plt, gene_strand):
        target_gene_strand = gene_strand.loc[:, [self.STRAND_GENE, self.STRAND_DIRECTION]]
        gene_strand_pre = target_gene_strand.drop_duplicates()

        # Determine gene renaming and gene direction
        gene_aliases = pd.read_csv(self.get_input_aliases(), sep=self.SEP_TAB)
        gene_strand_ed = pd.merge(gene_strand_pre, gene_aliases, on=[self.STRAND_GENE], how=self.LEFT)
        gene_strand_ed = gene_strand_ed.loc[:, [self.STRAND_DIRECTION, self.STRAND_GENE, self.ALIASES]]
        gene_strand_ed = gene_strand_ed.fillna(method=self.FFILL, axis=1)
        gene_strand_ed = gene_strand_ed.loc[:, [self.ALIASES, self.STRAND_DIRECTION]]
        gene_strand_ed = gene_strand_ed.rename(columns={self.ALIASES: self.STRAND_GENE})

        # Gene chromosome
        gene_chromosome = pd.read_csv(self.get_input_chromosome(), sep=self.SEP_TAB)

        # Aggregate the number of cases by base
        base_fus_plt_df = fus_plt.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.FUS_GENE1,
                                          self.CHROMOSOME1, self.START_POS1, self.END_POS1,
                                          self.FUS_GENE2, self.CHROMOSOME2,
                                          self.START_POS2, self.END_POS2]]
        group_fus_plt_df = base_fus_plt_df.groupby([self.MARKER_NAME, self.MUT_TYPE, self.FUS_GENE1, self.FUS_GENE2,
                                                   self.CHROMOSOME1, self.START_POS1, self.END_POS1, self.CHROMOSOME2,
                                                   self.START_POS2, self.END_POS2], as_index=False)
        fus_plt_ed = group_fus_plt_df.size()
        fus_plt_ed[self.START_POS] = 0

        # Combine Gene1 and Gene2
        fus_plt_ed[self.FUS_GENE1_GENE2] = (fus_plt_ed[self.CHROMOSOME1].astype(str) + ':'
                                            + fus_plt_ed[self.START_POS1].astype(str) + '-'
                                            + fus_plt_ed[self.END_POS1].astype(str) + '&' +
                                            fus_plt_ed[self.CHROMOSOME2].astype(str) + ':'
                                            + fus_plt_ed[self.START_POS2].astype(str) + '-'
                                            + fus_plt_ed[self.END_POS2].astype(str))
        fus_plt_ed[self.FUS_GENE1_GENE2] = fus_plt_ed[self.FUS_GENE1_GENE2].str.replace(' ', '')
        group_fus_plt_ed = fus_plt_ed.groupby([self.MARKER_NAME])[self.FUS_GENE1_GENE2]
        fus_plt_ed[self.RANK] = group_fus_plt_ed.rank(ascending=True)

        # Combine
        gene_test_bp_p = pd.merge(fus_plt_ed, gene_strand_ed, left_on=self.FUS_GENE1,
                                  right_on=self.STRAND_GENE, how=self.LEFT)
        gene_test_bp_p2 = pd.merge(gene_test_bp_p, gene_strand_ed, left_on=self.FUS_GENE2,
                                   right_on=self.STRAND_GENE, how=self.LEFT)
        target_fus_plt_ed = gene_test_bp_p2.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.FUS_GENE1, self.FUS_GENE2,
                                                    self.CHROMOSOME1, self.START_POS1, self.END_POS1, self.CHROMOSOME2,
                                                    self.START_POS2, self.END_POS2, self.SIZE, self.START_POS,
                                                    self.FUS_GENE1_GENE2, self.RANK, self.STRAND_X, self.STRAND_Y]]
        fus_plt_ed = target_fus_plt_ed.rename(columns={self.STRAND_X: self.GENE1_STRAND,
                                                       self.STRAND_Y: self.GENE2_STRAND})

        # Start position of gene1
        fus_plt_ed_g1_s = self.create_data_frame_position(fus_plt_ed, self.START_SIDE, 1)
        # Get start position
        fus_plt_ed_g1_s_pre_k = self.create_data_frame_chromosome(
            fus_plt_ed, gene_chromosome, self.START_SIDE, 1)
        # Get end position
        fus_plt_ed_g1_s_pre_s = self.create_data_frame_chromosome(
            fus_plt_ed, gene_chromosome, self.END_SIDE, 1
        )
        # End position of gene1
        fus_plt_ed_g1_e = self.create_data_frame_position(fus_plt_ed, self.END_SIDE, 1)
        # Start position of gene2
        fus_plt_ed_g2_s = self.create_data_frame_position(fus_plt_ed, self.START_SIDE, 2)
        # Get start position
        fus_plt_ed_g2_s_pre_k = self.create_data_frame_chromosome(
            fus_plt_ed, gene_chromosome, self.START_SIDE, 2)
        # Get end position
        fus_plt_ed_g2_s_pre_s = self.create_data_frame_chromosome(
            fus_plt_ed, gene_chromosome, self.END_SIDE, 2)
        # End position of gene2
        fus_plt_ed_g2_e = self.create_data_frame_position(fus_plt_ed, self.END_SIDE, 2)

        # UNION-Combine start location and the number by base
        fus_plt_ed2 = self.combine_data_frame(
            fus_plt_ed_g1_s, fus_plt_ed_g1_s_pre_k, fus_plt_ed_g1_s_pre_s, fus_plt_ed_g1_e)
        # UNION-Combine start location and the number by base
        fus_plt_ed3 = self.combine_data_frame(
            fus_plt_ed_g2_s, fus_plt_ed_g2_s_pre_k, fus_plt_ed_g2_s_pre_s, fus_plt_ed_g2_e)

        fus_plt_ed2 = self.add_number(fus_plt_ed2)
        fus_plt_ed3 = self.add_number(fus_plt_ed3)

        # UNION-Combine gene1 and gene2
        # Change the number type
        fus_plt_ed4 = fus_plt_ed2.append(fus_plt_ed3, ignore_index=True)
        fus_plt_ed4 = fus_plt_ed4.astype({self.POSITION: 'int64', self.SIZE: 'int64'})
        fus_plt_ed4 = fus_plt_ed4.rename(columns={self.SIZE: self.FREQ})

        fus_plt_ed5 = fus_plt_ed4.rank(method=self.DENSE)
        fus_plt_ed5.set_axis([self.MARKER_NAME_NUMBER, self.MUT_TYPE_NUMBER,
                              self.MARKER_NUMBER, self.STRAND_NUMBER, self.CHROMOSOME_NUMBER,
                              self.POSITION_NUMBER, self.FREQ_NUMBER, self.RANK_NUMBER,
                              self.GENE_NUMBER, self.START_AND_END_NUMBER,
                              self.CHR_NUM_NUMBER],
                             axis='columns', inplace=True)

        # Combine
        gene_fusion_plot = pd.concat([fus_plt_ed4, fus_plt_ed5], axis=1)
        gene_fusion_plot = gene_fusion_plot.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.MARKER, self.STRAND_DIRECTION,
                                                    self.CHROMOSOME, self.POSITION, self.FREQ, self.RANK,
                                                    self.FUS_GENE_NUM, self.START_AND_END, self.MARKER_NAME_NUMBER]]
        gene_fusion_plot = gene_fusion_plot.sort_values(
            by=[self.MARKER_NAME, self.MARKER_NAME_NUMBER], ascending=[True, False])

        return gene_fusion_plot
