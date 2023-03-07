#!/usr/bin/env python
# coding: utf-8
import gc
import os

import matplotlib.pyplot as plt
import pandas as pd

from .gene_analysis_plot import GeneAnalysisPlot


class GeneAnalysisPlotExceptFusion(GeneAnalysisPlot):
    OUT_FILENAME = 'geneplot_'
    draw_mut_list = [GeneAnalysisPlot.DELETION,
                     GeneAnalysisPlot.INSERTION,
                     GeneAnalysisPlot.DELINS,
                     GeneAnalysisPlot.AMPLIFICATION,
                     GeneAnalysisPlot.LOSS]
    sort_mut_list = [GeneAnalysisPlot.SNV] + draw_mut_list

    def draw_gene_analysis_except_fusion(self, input_rawdata, output_for_gene_analysis):
        # Importing data
        transcript_df = pd.read_csv(self.get_input_tid(), sep=self.SEP_TAB)
        panel_df = pd.read_csv(input_rawdata, sep=self.SEP_TAB)

        # Check if the number of data frame length is not zero
        if len(panel_df) == 0:
            print('No variant data in input data')
            return None

        gene_df = panel_df.query(self.CONNECT_AND.join(
            [self.NOT_TYPE_MSI, self.NOT_TYPE_TMB, self.NOT_TYPE_FUSION]))
        target_gene_df = self.exclude_non_target_gene(gene_df, transcript_df)
        panel_df = self.set_gene_default_position(target_gene_df, transcript_df)

        # Check if the number of gene data frame length is not zero
        if len(panel_df) == 0:
            print('No gene variant data in input data')
            return None

        gene_plt = self.get_gene_plt(panel_df)
        gene_plt_no_fus = self.get_gene_plt_no_fus(gene_plt, transcript_df)
        marker_name_num = gene_plt_no_fus[self.MARKER_NAME_NUMBER]

        for index in range(1, (int(max(marker_name_num)) + 1)):
            self.create_graph(gene_plt_no_fus[marker_name_num == index], output_for_gene_analysis)

    def create_graph(self, gene_plt_graph, output_for_gene_analysis):
        mut_type_pre = self.get_mut_type_pre(gene_plt_graph)

        # Sort Mut_type
        mut_type = sorted((list(mut_type_pre)), key=lambda x: self.sort_mutation_type(x))

        # Add a legend
        # SNV
        gene_plt_snv = self.create_data_frame(gene_plt_graph, self.SNV)
        # Deletion
        gene_plt_del = self.create_data_frame(gene_plt_graph, self.DELETION)
        # Insertion
        gene_plt_ins = self.create_data_frame(gene_plt_graph, self.INSERTION)
        # Delins
        gene_plt_delins = self.create_data_frame(gene_plt_graph, self.DELINS)
        # Amplification
        gene_plt_amp = self.create_data_frame(gene_plt_graph, self.AMPLIFICATION)
        # Loss
        gene_plt_loss = self.create_data_frame(gene_plt_graph, self.LOSS)

        # Drawing graphs
        fig = plt.figure(figsize=(40, 15), dpi=400)
        plt.rcParams[self.FONT_SIZE] = 12
        plt.rcParams[self.LEGEND_FONTSIZE] = 16
        plt.rcParams[self.LEGEND_HANDLELENGTH] = 3
        # Determine marker name to X-axis name
        marker_name = self.get_last_loc(gene_plt_graph, self.MARKER_NAME)
        transcript_id_version = self.get_last_loc(gene_plt_graph, self.TRANSCRIPT_ID_VERSION)

        plt.xlabel(
            f'$\it{marker_name}$ Gene (bp) \n {transcript_id_version}\n'
            f'(If there was no description of the CNA chromosome or position, '
            f'the $\it{marker_name}$ gene sequence was used.)',
            fontsize=30
        )
        # Determine name to Y-axis
        plt.ylabel(self.FREQUENCY, fontsize=30)
        plt.tick_params(labelsize=30)

        # Set scale to Y-axis
        plt.minorticks_on()
        plt.grid(color=self.BLACK, axis="y", linestyle="--")
        plt.grid(which=self.MINOR, axis="y", alpha=0.8, linestyle="--", linewidth=1)
        # Display gene name
        print(f"gene {marker_name}")
        save_fig_path = os.path.join(output_for_gene_analysis, f"{self.OUT_FILENAME}{marker_name}{self.OUT_EXTENSION}")

        # Mode frequency for each mutation type
        y_max = self.calculate_frequencies([gene_plt_snv, gene_plt_del, gene_plt_ins,
                                            gene_plt_delins, gene_plt_amp, gene_plt_loss])

        # Plot values in graph
        self.draw_graphical_stem(gene_plt_snv, self.SNV)
        self.draw_graphical_stem(gene_plt_del, self.DELETION)
        self.draw_graphical_stem(gene_plt_ins, self.INSERTION)
        self.draw_graphical_stem(gene_plt_delins, self.DELINS)
        self.draw_graphical_stem(gene_plt_amp, self.AMPLIFICATION)
        self.draw_graphical_stem(gene_plt_loss, self.LOSS)
        # Set scale to Y-axis
        self.draw_graphical_framework(y_max)

        # Add a legend
        self.draw_mutation_type_legend(mut_type)
        self.annotate_gene_graphs([gene_plt_snv, gene_plt_del, gene_plt_ins,
                                   gene_plt_delins, gene_plt_amp, gene_plt_loss])

        # Save the graph image
        plt.savefig(save_fig_path, bbox_inches=self.TIGHT, dpi=400)
        # Open memory
        fig.clf()
        # Close a saved graph
        plt.close(fig)

        # Delete figure instance
        del fig
        # Open memory bi garbage collection
        gc.collect()

    def get_gene_plt_no_fus(self, gene_plt, transcript_df):

        # gene_plt
        gene_plt_m = pd.merge(gene_plt, transcript_df, on=[self.MARKER_NAME], how=self.LEFT)
        gene_plt_m = gene_plt_m.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.CDS_CHANGE_NUM, self.START_POS1,
                                        self.END_POS1, self.START_P, self.END_P, self.LENGTH, self.STRAND_DIRECTION,
                                        self.TRANSCRIPT_ID_VERSION]]
        gene_plt_m[self.P_DIFF] = gene_plt_m.apply(
            lambda x: self.calculate_position_diff(x), axis=1
        )

        # Aggregate the number of cases by base
        gene_plt2_locs = [self.MARKER_NAME, self.MUT_TYPE, self.P_DIFF, self.TRANSCRIPT_ID_VERSION]
        gene_plt2 = gene_plt_m.loc[:, gene_plt2_locs].groupby(gene_plt2_locs, as_index=False).size()

        # Get start location
        gene_plt_st = gene_plt_m.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.TRANSCRIPT_ID_VERSION]]
        gene_plt_st = gene_plt_st.drop_duplicates()
        gene_plt_st[self.P_DIFF] = 0
        gene_plt_st[self.SIZE] = 0
        gene_plt_st.sort_values([self.MARKER_NAME, self.MUT_TYPE, self.P_DIFF])

        # UNION-Combine start location and the number by base
        gene_plt3 = gene_plt2.append(gene_plt_st, ignore_index=True)
        gene_plt3 = gene_plt3.sort_values(by=[self.MARKER_NAME, self.P_DIFF, self.MUT_TYPE])

        # Get end location
        gene_plt_end = gene_plt_m.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.TRANSCRIPT_ID_VERSION, self.LENGTH]]
        gene_plt_end = gene_plt_end.drop_duplicates()
        gene_plt_end[self.SIZE] = 0
        gene_plt_end = gene_plt_end.rename(columns={self.LENGTH: self.P_DIFF})

        # UNION-Combine end location and the number by base
        gene_plt_ed = gene_plt3.append(gene_plt_end, ignore_index=True)
        gene_plt_ed = gene_plt_ed.sort_values(by=[self.MARKER_NAME, self.P_DIFF, self.MUT_TYPE])
        gene_plt_ed = gene_plt_ed.loc[:, [self.MARKER_NAME, self.MUT_TYPE,
                                          self.P_DIFF, self.TRANSCRIPT_ID_VERSION, self.SIZE]]

        # Determine rename a column and Marker_Name_No (prepare to loop)
        gene_plt_pre = gene_plt_ed.rank(method=self.DENSE)
        gene_plt_pre.set_axis(
            [self.MARKER_NAME_NUMBER, self.MUT_TYPE_NUMBER, self.P_DIFF_NO, self.SIZE_NO,
             self.TRANSCRIPT_ID_VERSION_NUMBER], axis=self.COLUMNS, inplace=True
        )

        # Combine
        gene_plt_no_fus = pd.concat([gene_plt_ed, gene_plt_pre], axis=1)
        gene_plt_no_fus = gene_plt_no_fus.loc[:, [self.MARKER_NAME, self.TRANSCRIPT_ID_VERSION, self.MUT_TYPE,
                                                  self.P_DIFF, self.SIZE, self.MARKER_NAME_NUMBER]]
        gene_plt_no_fus = gene_plt_no_fus.sort_values(
            by=[self.MARKER_NAME, self.MARKER_NAME_NUMBER, self.MUT_TYPE, self.P_DIFF],
            ascending=[True, False, True, True]
        )
        return gene_plt_no_fus

    def get_gene_plt(self, panel_df):
        df_type = panel_df.loc[:, [self.MARKER_NAME, self.CDS_CHANGE, self.REFERENCE_ALLELE, self.ALTERNATE_ALLELE,
                                   self.SAMPLE_ID, self.TYPE, self.CHROMOSOME1, self.START_POS1, self.END_POS1]]
        df_type = self.set_allele_info(df_type)

        # Get only base from column
        cds_change_num = df_type.Cds_Change
        df_type[self.CDS_CHANGE_NUM] = cds_change_num.str.extract(self.REG_NUMBER_REP).fillna(0).astype(int)

        # Get base and position from column
        gene_plt = df_type.loc[:, [self.MARKER_NAME, self.MUT_TYPE,
                                   self.CDS_CHANGE_NUM, self.START_POS1, self.END_POS1]]
        gene_plt = gene_plt.sort_values([self.MARKER_NAME, self.MUT_TYPE, self.CDS_CHANGE_NUM])
        return gene_plt

    def calculate_position_diff(self, no_fus_df):
        """
        Subtract mutation position to gene start position
        """
        if no_fus_df.strand == '-':
            return int(no_fus_df.end_P) - int(no_fus_df.Start1)
        else:
            return int(no_fus_df.Start1) - int(no_fus_df.start_P)

    def sort_mutation_type(self, mut_type):
        return super().sort_mutation_type(mut_type, self.sort_mut_list)

    def create_data_frame(self, non_fus_df, mut_type):
        return super().create_data_frame(non_fus_df, mut_type, self.P_DIFF)

    def draw_graphical_stem(self, mut_type_df, mut_type):
        return super().draw_graphical_stem(mut_type_df, mut_type, self.draw_mut_list)

    def draw_mutation_type_legend(self, mut_type_list):
        return super().draw_mutation_type_legend(mut_type_list, 6)
