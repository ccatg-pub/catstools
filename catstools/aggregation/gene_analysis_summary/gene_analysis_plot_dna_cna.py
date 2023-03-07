#!/usr/bin/env python
# coding: utf-8
import gc
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .gene_analysis_plot import GeneAnalysisPlot


class GeneAnalysisPlotDnaCna(GeneAnalysisPlot):

    TYPE_AMP = 'Type=="Amp"'
    TYPE_LOSS = 'Type=="Loss"'
    MUT_AMPLIFICATION = 'Mut_type=="Amplification"'
    MUT_LOSS = 'Mut_type=="Loss"'
    AMPERSAND = '&'
    SANDE_START = 'SandE=="Start"'
    SANDE_END = 'SandE=="End"'
    SANDE_SANDE = 'SandE=="SandE"'
    OUT_FILENAME = 'cna_geneplot_'

    def draw_gene_analysis_cna(self, input_rawdata, output_for_gene_analysis):
        # Exclude if the end letter is Nan
        transcript_df = pd.read_csv(self.get_input_tid(), sep=self.SEP_TAB)
        # Rename a column (F1_gene_name→Marker_Name)
        transcript_df_r = transcript_df.dropna(subset=[self.AMINO_ACIDS])

        # PANEL
        panel_df_l = pd.read_csv(input_rawdata, sep=self.SEP_TAB)

        # Check if the number of data frame length is not zero
        if len(panel_df_l) == 0:
            print('No variant data in input data')
            return None

        cna_plt_base = panel_df_l.query(self.CONNECT_OR.join([self.TYPE_AMP, self.TYPE_LOSS]))
        target_cna_plt = self.exclude_non_target_gene(cna_plt_base, transcript_df_r)
        cna_plt_pre = self.set_gene_default_position(target_cna_plt, transcript_df_r)

        # Check if the number of cna data frame length is not zero
        if len(cna_plt_pre) == 0:
            print('No cna variant data in input data')
            return None

        cna_plt = self.get_cna_plt(self.get_cna_plt_pre(cna_plt_pre), transcript_df_r)
        gene_cna_plt = self.get_gene_cna_plt(cna_plt)
        marker_name_num = gene_cna_plt[self.MARKER_NAME_NUMBER]

        for index in range(1, (int(max(marker_name_num)) + 1)):
            self.create_graph(gene_cna_plt[marker_name_num == index], output_for_gene_analysis)

    def create_graph(self, gene_cna_graph, output_for_gene_analysis):
        mut_type_pre = self.get_mut_type_pre(gene_cna_graph)

        # Rename columns and select genes
        # gene : Amplification_Start
        gene_cna_amp_s = self.create_frame(gene_cna_graph, self.MUT_AMPLIFICATION, self.START_SIDE)
        # gene : Amplification_End
        gene_cna_amp_e = self.create_frame(gene_cna_graph, self.MUT_AMPLIFICATION, self.END_SIDE)
        # gene : Amplification_Start-END
        gene_cna_amp_se = self.create_frame(gene_cna_graph, self.MUT_AMPLIFICATION, self.START_AND_END)

        # gene：Loss_Start
        gene_cna_loss_s = self.create_frame(gene_cna_graph, self.MUT_LOSS, self.START_SIDE)
        # gene：Loss_End
        gene_cna_loss_e = self.create_frame(gene_cna_graph, self.MUT_LOSS, self.END_SIDE)
        # gene：Loss_Start-END
        gene_cna_loss_se = self.create_frame(gene_cna_graph, self.MUT_LOSS, self.START_AND_END)

        # Drawing graphs
        fig = plt.figure(figsize=(30, 15), dpi=400)
        plt.rcParams[self.FONT_SIZE] = 12
        plt.rcParams[self.LEGEND_FONTSIZE] = 30
        plt.rcParams[self.LEGEND_HANDLELENGTH] = 3

        ax1 = fig.add_axes((0.5, 0.5, 0.5, 0.6))
        ax2 = fig.add_axes((0.5, 0.5, 0.5, 0.6))

        # Determine name to Y-axis
        plt.ylabel(self.FREQUENCY, fontsize=30)
        plt.tick_params(labelsize=30)

        # Don't display the value on the scale of ax1 X-axis and Y-axis
        ax1.set_xticks([])
        ax1.set_yticks([])

        # Align labels
        fig.align_labels()

        # Display Gene1
        self.display_gene_graph(
            ax1, ax2, mut_type_pre, gene_cna_graph, gene_cna_amp_s, gene_cna_amp_e, gene_cna_amp_se,
            gene_cna_loss_s, gene_cna_loss_e, gene_cna_loss_se, output_for_gene_analysis
        )

        # Open memory
        fig.clf()
        # Close a saved graph
        plt.close(fig)

        # Delete figure instance
        del fig
        # Open memory by garbage collection
        gc.collect()

    def mut_type_sort(self, type):
        """
        Sorting function of Mut_type
        """
        if type == self.AMPLIFICATION:
            return 1
        elif type == self.LOSS:
            return 2
        else:
            return 3

    def mut_type_def(self, mut_type_list, marker):
        """
        Add a legend
        """
        if len(mut_type_list) == 2:
            plt.legend(
                fontsize=14, loc=self.BEST,
                labels=[self.AMPLIFICATION_START_POSITION, self.AMPLIFICATION_END_POSITION,
                        self.LOSS_START_POSITION, self.LOSS_END_POSITION],
                title=f'$\it{marker}$ breakpoint',
                labelcolor=[self.ORANGE_RED, self.ORANGE, self.DARKBLUE, self.BLUE]
            )
        else:
            if mut_type_list == [self.AMPLIFICATION]:
                plt.legend(
                    fontsize=14, loc=self.BEST,
                    labels=[self.AMPLIFICATION_START_POSITION,
                            self.AMPLIFICATION_END_POSITION],
                    title=f'$\it{marker}$ breakpoint',
                    labelcolor=[self.ORANGE_RED, self.ORANGE]
                )
            elif mut_type_list == [self.LOSS]:
                plt.legend(
                    fontsize=14, loc=self.BEST,
                    labels=[self.LOSS_START_POSITION, self.LOSS_END_POSITION],
                    title=f'$\it{marker}$ breakpoint',
                    labelcolor=[self.DARKBLUE, self.BLUE]
                )
            else:
                pass

    def get_cna_plt_pre(self, in_cna_plt_pre):
        cna_plt_pre = in_cna_plt_pre.loc[:, [self.MARKER_NAME, self.TYPE, self.CHROMOSOME1,
                                             self.START_POS1, self.END_POS1]]

        cna_plt_pre[self.MUT_TYPE] = cna_plt_pre.apply(
            lambda c: self.AMPLIFICATION if c[self.TYPE] == self.AMP
            else (self.LOSS if c[self.TYPE] == self.LOSS
                  else c[self.TYPE]
                  ), axis=1)
        cna_plt_pre[self.CHROMOSOME1] = cna_plt_pre[self.CHROMOSOME1].str.strip()
        cna_plt_pre[self.START_POS1] = cna_plt_pre[self.START_POS1].astype(str).str.strip().astype(int)

        # Amplification or Loss
        # Case of base
        cna_plt_pre[self.END_POS1] = cna_plt_pre[self.END_POS1].astype(str).str.strip().astype(int)

        return cna_plt_pre

    def get_cna_plt(self, cna_plt_pre, transcript_df_r):
        # Get TranscriptIDVersion and frequency of gene
        cna_plt = pd.merge(cna_plt_pre, transcript_df_r, on=[self.MARKER_NAME], how=self.LEFT)
        cna_plt = cna_plt.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.CHROMOSOME1, self.START_POS1, self.END_POS1,
                                  self.START_P, self.END_P, self.LENGTH, self.STRAND_DIRECTION,
                                  self.TRANSCRIPT_ID_VERSION]]
        cna_plt = cna_plt.groupby(
            [self.MARKER_NAME, self.MUT_TYPE, self.CHROMOSOME1, self.START_POS1, self.END_POS1, self.START_P,
             self.END_P, self.LENGTH, self.STRAND_DIRECTION, self.TRANSCRIPT_ID_VERSION],
            as_index=False
        )
        cna_plt = cna_plt.size()

        # Ranked by gene
        cna_plt[self.START1_END1] = (cna_plt[self.CHROMOSOME1] + ":"
                                     + cna_plt[self.START_POS1].astype(str) + "-"
                                     + cna_plt[self.END_POS1].astype(str))
        cna_plt[self.START1_END1] = cna_plt[self.START1_END1].str.strip()
        cna_plt_rank = cna_plt.groupby([self.MARKER_NAME, self.MUT_TYPE])[self.START1_END1]
        cna_plt[self.RANK] = cna_plt_rank.rank(ascending=True)
        cna_plt[self.START_P_PDIFF] = self.get_diff(cna_plt, self.START_P, self.START_P)
        cna_plt[self.END_P_PDIFF] = self.get_diff(cna_plt, self.END_P, self.START_P)
        cna_plt[self.START1_PDIFF] = self.get_diff(cna_plt, self.START_POS1, self.START_P)
        cna_plt[self.END1_PDIFF] = self.get_diff(cna_plt, self.END_POS1, self.START_P)
        return cna_plt

    def get_gene_cna_plt(self, in_cna_plt):
        cna_plt = in_cna_plt
        # Start location of first gene
        cna_plt_gene_s = self.create_cna_position(cna_plt, self.START_SIDE)
        # Get start location
        cna_plt_gene_s_pre_k = self.create_gene_position(cna_plt, self.START_SIDE)
        # Get end location
        cna_plt_gene_s_pre_s = self.create_gene_position(cna_plt, self.END_SIDE)
        # End location of first gene
        cna_plt_gene_e = self.create_cna_position(cna_plt, self.END_SIDE)

        # UNION-Combine start location and the number by base
        cna_plt = cna_plt_gene_s.append(cna_plt_gene_s_pre_k, ignore_index=True)
        cna_plt = cna_plt.append(cna_plt_gene_s_pre_s, ignore_index=True)
        cna_plt = cna_plt.append(cna_plt_gene_e, ignore_index=True)
        cna_plt[self.MARKER_NAME_NUMBER] = cna_plt[self.MARKER_NAME].rank(method=self.DENSE)

        return cna_plt.sort_values(by=[self.MARKER_NAME, self.POSITION])

    def normalize_decimal(self, f_val):
        """
        Hide after the decimal point
        """
        return str(f_val).rstrip('0').rstrip('.')

    def get_diff(self, cna_df, base_point, sub_point):
        """
        get position difference value from dataframe
        """
        return self.convert_int(cna_df[base_point]) - self.convert_int(cna_df[sub_point])

    def convert_int(self, val):
        return val.astype(str).str.replace(',', '').astype(float).astype(int)

    def create_cna_position(self, cna_df, start_end):
        """
        Create data frames based on gene (start/end) position information
        """
        # Determining the column name to use,
        # depending on start_end_side
        if start_end == self.START_SIDE:
            pos_col = self.START_POS1
            pos_diff_col = self.START1_PDIFF
        else:
            pos_col = self.END_POS1
            pos_diff_col = self.END1_PDIFF

        conv_cna_df = cna_df.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.STRAND_DIRECTION, self.TRANSCRIPT_ID_VERSION,
                                     self.CHROMOSOME1, pos_col, self.SIZE, self.RANK, self.START_P_PDIFF,
                                     self.END_P_PDIFF, self.START1_PDIFF, self.END1_PDIFF]]
        conv_cna_df = conv_cna_df.rename(columns={pos_diff_col: self.POSITION,
                                         self.CHROMOSOME1: self.CHROMOSOME, self.SIZE: self.FREQ})
        conv_cna_df[self.START_AND_END] = start_end

        return conv_cna_df

    def create_gene_position(self, cna_df, start_end):
        """
        Create data frames based on gene position information
        """
        # Determining the column name to use,
        # depending on start_end_side
        if start_end == self.START_SIDE:
            pos_col = self.START_P
            pos_diff_col = self.START_P_PDIFF
        else:
            pos_col = self.END_P
            pos_diff_col = self.END_P_PDIFF

        conv_cna_df = cna_df.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.STRAND_DIRECTION,
                                     self.TRANSCRIPT_ID_VERSION, self.CHROMOSOME1, pos_col,
                                     self.START_P_PDIFF, self.END_P_PDIFF, self.START1_PDIFF,
                                     self.END1_PDIFF]]
        conv_cna_df = conv_cna_df.drop_duplicates()
        conv_cna_df = conv_cna_df.rename(
            columns={pos_diff_col: self.POSITION, self.CHROMOSOME1: self.CHROMOSOME}
        )
        conv_cna_df[self.START_AND_END] = self.START_AND_END
        conv_cna_df[self.FREQ] = 0

        return conv_cna_df

    def create_frame(self, cna_df, mut_type, start_end_type):
        """
        create input data frame for graph plotting
        """

        cna_query_df = cna_df.query(self.get_search_terms(mut_type, start_end_type))

        column_list = [self.MARKER_NAME, self.MUT_TYPE, self.STRAND_DIRECTION, self.CHROMOSOME]

        # start side
        if start_end_type == self.START_SIDE:
            column_list.extend([self.RANK, self.POSITION, self.FREQ, self.START_POS1])
        # end side
        elif start_end_type == self.END_SIDE:
            column_list.extend([self.RANK, self.POSITION, self.FREQ, self.END_POS1])
        # start and end
        elif start_end_type == self.START_AND_END:
            column_list.extend([self.POSITION, self.FREQ, self.START_POS1, self.START_P, self.END_P])
        else:
            raise ValueError(f"Error: start_end_type invalid. val={start_end_type}")

        return cna_query_df.loc[:, column_list]

    def get_search_terms(self, mut_type, start_end):

        if self.START_SIDE == start_end:
            sande = self.SANDE_START
        elif self.END_SIDE == start_end:
            sande = self.SANDE_END
        else:
            sande = self.SANDE_SANDE

        return f"{mut_type} {self.AMPERSAND} {sande}"

    def annotate_cna_graph(self, cna_pos, cna_freq, cna_rank,
                           left_right, start_end_type, plot_color):
        """
        annotate text to cna graph plot
        """
        y_coordinate = self.get_y_coordinate(cna_rank, start_end_type)
        x_coordinate = self.get_x_coordinate(start_end_type)

        # annotate text to cna graph plot
        plt.annotate('{:.0f}'.format(cna_pos), xy=(cna_pos, cna_freq),
                     xytext=(x_coordinate, y_coordinate),
                     textcoords=self.OFFSET_POINTS, ha=left_right, va=self.BOTTOM,
                     fontsize=10, color=plot_color)

    def get_y_coordinate(self, cna_rank, start_end_type):

        # define coordinate to annotate text
        # Determine y-coordinate
        if 1 <= cna_rank <= 6:
            y_coordinate = (cna_rank * 17) - 12
        elif 7 <= cna_rank <= 10:
            y_coordinate = (cna_rank * -17) + 107
        else:
            if start_end_type == self.START_SIDE:
                y_coordinate = 107
            else:
                y_coordinate = 90
        return y_coordinate

    def get_x_coordinate(self, start_end_type):

        # start side
        if start_end_type == self.START_SIDE:
            x_coordinate = -7
        else:
            x_coordinate = 5
        return x_coordinate

    def display_gene_graph(self, axis1, axis2, mut_type_df, cna_df, amp_start_df,
                           amp_end_df, amp_start_and_end_df, loss_start_df,
                           loss_end_df, loss_start_and_end_df, output_dir):
        """
        Display Gene
        """
        y_max_amp = 0
        # In the case of Amp
        if len(amp_start_df) != 0:
            # Specify the gene name of a gene
            # (# Specify marker name to X-axis name)
            marker_name = self.get_last_loc(amp_start_df, self.MARKER_NAME)
            # Specify the gene name pf the chromosome
            chromosome = self.get_last_loc(amp_start_df, self.CHROMOSOME)
            # Specify the orientation of gene
            strand = self.get_last_loc(amp_start_df, self.STRAND_DIRECTION)

            self.draw_gene_amp(amp_start_df, amp_end_df, amp_start_and_end_df)

            # Maximum value of yAmp
            y_max_amp = int(max(amp_start_df.freq)) + 2

            # Display gene name
            print(f"{self.STRAND_GENE} {marker_name} {self.AMPLIFICATION}")

        # In case of Loss
        y_max_loss = 0
        if len(loss_start_df) != 0:
            marker_name = self.get_last_loc(loss_start_df, self.MARKER_NAME)
            chromosome = self.get_last_loc(loss_start_df, self.CHROMOSOME)
            strand = self.get_last_loc(loss_start_df, self.STRAND_DIRECTION)

            self.draw_gene_loss(loss_start_df, loss_end_df, loss_start_and_end_df)

            # Maximum value of yLoss
            y_max_loss = int(max(loss_start_df.freq)) + 2

            # Display gene name
            print(f"{self.STRAND_GENE} {marker_name} {self.LOSS}")

        # Sort Mut_type
        mut_type = sorted((list(mut_type_df)), key=lambda x: self.mut_type_sort(x))

        # Display a legend
        self.mut_type_def(mut_type, marker_name)

        # Determine name to X-axis of axi1
        if strand == '+':
            arrows = '\n →  →  →  →  →  →  →  →  →  →'
            axis_color = self.RED
        else:
            arrows = '\n ←  ←  ←  ←  ←  ←  ←  ←  ←  ←'
            axis_color = self.BLUE
        axis1.set_xlabel(arrows, fontsize=30, color=axis_color, linespacing=0.5)
        axis2.set_xlabel(
            f'\n Chromosome{chromosome} \n $\it{marker_name}$ CNA\n'
            f'(If there was no description of the chromosome or position, '
            f'the $\it{marker_name}$ gene sequence was used.)',
            fontsize=18
        )

        # Determine name to Y-axis of plt
        plt.ylabel(self.FREQUENCY, fontsize=18)

        x_min_p = int(min(cna_df.position))
        x_max_p = int(max(cna_df.position))
        arg_ary = [self.START_P, self.END_P, self.START_POS1, self.END_POS1]
        x_min_rp = cna_df[arg_ary].astype(float).min().min()
        x_max_rp = cna_df[arg_ary].astype(float).max().max()

        y_max = self.get_y_max(amp_start_df, loss_start_df, y_max_amp, y_max_loss)

        # Specify baseline
        baseline = plt.stem([x_min_p, x_max_p], [0, 0], linefmt=None, markerfmt=' ', basefmt=None)
        plt.setp(baseline, color=self.RED)

        if y_max < 6:
            plt.yticks(np.arange(0, y_max, 1), fontsize=18)
            plt.grid(which=self.MAJOR, axis="y", alpha=0.8, linestyle="--", linewidth=1)
        else:
            plt.yticks(np.arange(0, y_max, 5), fontsize=18)
            plt.minorticks_on()
            plt.grid(color=self.BLACK, axis="y", linestyle="--")
            plt.grid(which=self.MINOR, axis="y", alpha=0.8, linestyle="--", linewidth=1)
        plt.xticks((x_min_p, x_max_p), (self.normalize_decimal(x_min_rp), self.normalize_decimal(x_max_rp)),
                   fontsize=18)

        # Specify the save location
        savefig_pass = os.path.join(output_dir, f"{self.OUT_FILENAME}{marker_name}{self.OUT_EXTENSION}")
        plt.savefig(savefig_pass, bbox_inches=self.TIGHT, dpi=400)

    def draw_gene_amp(self, amp_start_df, amp_end_df, amp_start_and_end_df):
        # Draw Gene-AMP_START
        markerline, stemlines, baseline = plt.stem(amp_start_df.position, amp_start_df.freq,
                                                   linefmt=self.ORANGE_RED, markerfmt='o')
        markerline.set_color(self.ORANGE_RED)
        for pos, freq, rank in zip(amp_start_df[self.POSITION], amp_start_df[self.FREQ], amp_start_df[self.RANK]):
            self.annotate_cna_graph(pos, freq, rank, self.RIGHT, self.START_SIDE, self.ORANGE_RED)

        # Draw Gene-AMP_END
        markerline, stemlines, baseline = plt.stem(amp_end_df.position, amp_end_df.freq,
                                                   linefmt=self.ORANGE, markerfmt='o')
        markerline.set_color(self.ORANGE)
        markerline.set_markerfacecolor(self.NONE)
        for pos, freq, rank in zip(amp_end_df[self.POSITION], amp_end_df[self.FREQ], amp_start_df[self.RANK]):
            self.annotate_cna_graph(pos, freq, rank, self.LEFT, self.END_SIDE, self.ORANGE)

        # Draw Gene-Amp_STARTEND
        for pos, rank in zip(amp_start_and_end_df[self.POSITION], amp_start_and_end_df[self.FREQ]):
            plt.annotate('{:.0f}'.format(pos), xy=(pos, rank), xytext=(0, 5), textcoords=self.OFFSET_POINTS,
                         ha=self.CENTER, fontsize=10)

    def draw_gene_loss(self, loss_start_df, loss_end_df, loss_start_and_end_df):
        # Draw Gene-Loss_START
        markerline, stemlines, baseline = plt.stem(loss_start_df.position, loss_start_df.freq,
                                                   linefmt=self.DARKBLUE, markerfmt='o')
        markerline.set_color(self.DARKBLUE)
        for pos, freq, rank in zip(loss_start_df[self.POSITION], loss_start_df[self.FREQ],
                                   loss_start_df[self.RANK]):
            self.annotate_cna_graph(pos, freq, rank, self.RIGHT, self.START_SIDE, self.DARKBLUE)

        # Draw Gene-Loss_END
        markerline, stemlines, baseline = plt.stem(
            loss_end_df.position, loss_end_df.freq, linefmt=self.BLUE, markerfmt='o')
        markerline.set_color(self.BLUE)
        markerline.set_markerfacecolor(self.NONE)
        for pos, freq, rank in zip(loss_end_df[self.POSITION], loss_end_df[self.FREQ], loss_end_df[self.RANK]):
            self.annotate_cna_graph(pos, freq, rank, self.LEFT, self.END_SIDE, self.BLUE)

        # Draw Gene-Loss_STARTEND
        for pos, freq in zip(loss_start_and_end_df[self.POSITION], loss_start_and_end_df[self.FREQ]):
            plt.annotate('{:.0f}'.format(pos), xy=(pos, freq), xytext=(0, 5),
                         textcoords=self.OFFSET_POINTS, ha=self.CENTER, fontsize=10)

    def get_y_max(self, amp_start_df, loss_start_df, y_max_amp, y_max_loss):
        # Amp
        if len(amp_start_df) > 0:
            y_max = y_max_amp
        # Loss
        elif len(loss_start_df) > 0:
            y_max = y_max_loss
        else:
            y_max = max(y_max_amp, y_max_loss)
        return y_max
