#!/usr/bin/env python
# coding: utf-8
import gc
import os

import matplotlib.pyplot as plt
import pandas as pd

from .gene_analysis_plot import GeneAnalysisPlot


class GeneAnalysisPlotAminoAcids(GeneAnalysisPlot):
    NOT_TYPE_AMP = 'not Type=="Amp"'
    NOT_TYPE_LOSS = 'not Type=="Loss"'
    MUT_SNV = 'Mut_type=="SNV"'
    MUT_DELETION = 'Mut_type=="Deletion"'
    MUT_INSERTION = 'Mut_type=="Insertion"'
    MUT_DELINS = 'Mut_type=="Delins"'
    OUT_FILENAME = 'aminoplot_'

    draw_mut_list = [GeneAnalysisPlot.DELETION,
                     GeneAnalysisPlot.INSERTION,
                     GeneAnalysisPlot.DELINS]

    sort_mut_list = [GeneAnalysisPlot.SNV] + draw_mut_list

    def draw_gene_analysis_amino(self, input_rawdata, output_for_gene_analysis):
        transcript_df = pd.read_csv(self.get_input_tid(), sep=self.SEP_TAB)
        # Exclude if the end letter is Nan
        transcript_df_r = transcript_df.dropna(subset=[self.AMINO_ACIDS])

        panel_df_pre = pd.read_csv(input_rawdata, sep=self.SEP_TAB)

        # Check if the number of data frame length is not zero
        if len(panel_df_pre) == 0:
            print('No variant data in input data')
            return None

        panel_df = panel_df_pre.query(self.CONNECT_AND.join([self.NOT_TYPE_MSI, self.NOT_TYPE_TMB, self.NOT_TYPE_FUSION,
                                                             self.NOT_TYPE_AMP, self.NOT_TYPE_LOSS]))
        panel_df = self.exclude_non_target_gene(panel_df, transcript_df_r)

        # Check if the number of short-variant data frame length is not zero
        if len(panel_df) == 0:
            print('No short-variant data in input data')
            return None

        aa_plt_graph_pre = self.get_aa_plt_graph_pre(panel_df, transcript_df_r)
        marker_name_number = aa_plt_graph_pre[self.MARKER_NAME_NUMBER]

        for index in range(1, (int(max(marker_name_number)) + 1)):
            self.create_graph(aa_plt_graph_pre[marker_name_number == index], output_for_gene_analysis)

    def create_graph(self, aa_plt_graph, output_for_gene_analysis):
        mut_type_pre = self.get_mut_type_pre(aa_plt_graph)
        # Sort Mut_type
        mut_type = sorted((list(mut_type_pre)), key=lambda x: self.sort_mutation_type(x))

        # Rename columns and select genes
        # SNV
        aa_plt_snv = self.create_data_frame(aa_plt_graph, self.SNV)
        # Deletion
        aa_plt_del = self.create_data_frame(aa_plt_graph, self.DELETION)
        # Insertion
        aa_plt_ins = self.create_data_frame(aa_plt_graph, self.INSERTION)
        # Delins
        aa_plt_delins = self.create_data_frame(aa_plt_graph, self.DELINS)

        fig = plt.figure(figsize=(40, 15), dpi=400)
        plt.rcParams[self.FONT_SIZE] = 12
        plt.rcParams[self.LEGEND_FONTSIZE] = 16
        plt.rcParams[self.LEGEND_HANDLELENGTH] = 3

        # Determine marker name to X-axis name
        marker_name = self.get_last_loc(aa_plt_graph, self.MARKER_NAME)
        protein_id = self.get_last_loc(aa_plt_graph, self.PROTEIN_ID)

        # Change the notes
        plt.xlabel(
            f'{marker_name}({protein_id}) AminoAcids (aa) \n CNAs and '
            f'rearrangements are not aggregated and therefore not shown '
            f'in the graph.', fontsize=30
        )

        # Save graphs
        print(f"Amino Acid {marker_name}")
        save_fig_aa_pig = os.path.join(output_for_gene_analysis,
                                       f"{self.OUT_FILENAME}{marker_name}{self.OUT_EXTENSION}")

        # Determine name to Y-axis
        plt.ylabel(self.FREQUENCY, fontsize=30)
        plt.tick_params(labelsize=30)

        # Mode frequency for each mutation type
        y_max = self.calculate_frequencies([aa_plt_snv, aa_plt_del, aa_plt_ins, aa_plt_delins])

        # Plot values in graph
        self.draw_graphical_stem(aa_plt_snv, self.SNV)
        self.draw_graphical_stem(aa_plt_del, self.DELETION)
        self.draw_graphical_stem(aa_plt_ins, self.INSERTION)
        self.draw_graphical_stem(aa_plt_delins, self.DELINS)
        # Set scale to Y-axis
        self.draw_graphical_framework(y_max)

        # Add a legend
        self.draw_mutation_type_legend(mut_type)
        self.annotate_gene_graphs([aa_plt_snv, aa_plt_del, aa_plt_ins, aa_plt_delins])

        # Save the graph image
        plt.savefig(save_fig_aa_pig, bbox_inches=self.TIGHT, dpi=400)
        # Open memory
        fig.clf()
        # Close a saved graph
        plt.close(fig)

        # Delete figure instance
        del fig
        # Open memory by garbage collection
        gc.collect()

    def get_aa_plt_graph_pre(self, panel_df, transcript_df_r):
        panel_df = self.set_allele_info(panel_df)
        aa_plt = panel_df.loc[:, [self.SAMPLE_ID, self.TYPE, self.VARIANT_TYPE_STATE, self.AMINO_ACIDS_CHANGE,
                                  self.CDS_CHANGE, self.MUT_TYPE, self.MARKER_NAME, self.START_POS1, self.END_POS1]]
        aa_plt = aa_plt.query(self.CONNECT_OR.join([self.MUT_SNV, self.MUT_DELETION, self.MUT_INSERTION,
                                                    self.MUT_DELINS]))

        # Exclude words containing splice site and promoter
        aa_plt = aa_plt.loc[~aa_plt[self.AMINO_ACIDS_CHANGE].str.startswith(self.PROMOTER)]
        aa_plt = aa_plt.loc[~aa_plt[self.AMINO_ACIDS_CHANGE].str.startswith(self.SPLICE_SITE)]

        # Get only amino acids from column
        aa_plt[self.AMINO_ACIDS_NUM] = aa_plt.Amino_Acids_Change.str.extract(self.REG_NUMBER_REP)
        aa_plt[self.AMINO_ACIDS_NUM] = aa_plt[self.AMINO_ACIDS_NUM].fillna(0).astype(int)
        aa_plt2 = aa_plt.loc[:, [self.MARKER_NAME, self.MUT_TYPE, self.AMINO_ACIDS_NUM]]
        aa_plt2 = aa_plt2.groupby([self.MARKER_NAME, self.MUT_TYPE, self.AMINO_ACIDS_NUM], as_index=False)
        aa_plt2 = aa_plt2.size()
        aa_plt_id = pd.merge(aa_plt2, transcript_df_r, on=[self.MARKER_NAME], how=self.LEFT)
        aa_plt_id = aa_plt_id.loc[:, [self.MARKER_NAME, self.PROTEIN_ID,
                                      self.MUT_TYPE, self.AMINO_ACIDS_NUM, self.AMINO_ACIDS, self.SIZE]]
        # Get start location
        aa_plt_id_uq = aa_plt_id.loc[:, [self.MARKER_NAME, self.PROTEIN_ID, self.MUT_TYPE]]
        aa_plt_id_uq = aa_plt_id_uq.drop_duplicates()
        aa_plt_id_uq[self.AMINO_ACIDS_NUM] = 0
        aa_plt_id_uq[self.SIZE] = 0

        # UNION-Combine start location and the number by amino acid
        aa_plt_m = aa_plt_id.loc[:, [self.MARKER_NAME, self.PROTEIN_ID, self.MUT_TYPE, self.AMINO_ACIDS_NUM, self.SIZE]]
        aa_plt_m = aa_plt_m.append(aa_plt_id_uq, ignore_index=True)
        aa_plt_m = aa_plt_m.sort_values(by=[self.MARKER_NAME, self.AMINO_ACIDS_NUM, self.MUT_TYPE])

        # Get end location
        aa_plt_aa = aa_plt_id.loc[:, [self.MARKER_NAME, self.PROTEIN_ID, self.MUT_TYPE, self.AMINO_ACIDS]]
        aa_plt_aa = aa_plt_aa.rename(columns={self.AMINO_ACIDS: self.AMINO_ACIDS_NUM})
        aa_plt_aa = aa_plt_aa.drop_duplicates()
        aa_plt_aa[self.SIZE] = 0

        # UNION-Combine end location and the number by amino acid
        aa_plt_ed_test = aa_plt_m.append(aa_plt_aa, ignore_index=True)
        aa_plt_ed_test = aa_plt_ed_test.sort_values(by=[self.MARKER_NAME, self.AMINO_ACIDS_NUM, self.MUT_TYPE])
        aa_plt_ed_test[self.AMINO_ACIDS_NUM] = (aa_plt_ed_test[self.AMINO_ACIDS_NUM].astype(int))

        # Determine rename a column and Marker_Name_No (prepare to loop)
        aa_plt_ed_test2 = aa_plt_ed_test.rank(method=self.DENSE)
        aa_plt_ed_test2.set_axis(
            [self.MARKER_NAME_NUMBER, self.MUT_TYPE_NUMBER, self.PROTEIN_ID_NUMBER, self.SIZE_NUMBER,
             self.AMINO_ACIDS_NUM_NUMBER], axis=self.COLUMNS, inplace=True
        )

        # combine
        aa_plt_graph_pre = pd.concat([aa_plt_ed_test, aa_plt_ed_test2], axis=1)
        aa_plt_graph_pre = aa_plt_graph_pre.loc[:, [self.MARKER_NAME, self.PROTEIN_ID, self.MUT_TYPE,
                                                    self.AMINO_ACIDS_NUM, self.SIZE, self.MARKER_NAME_NUMBER]]
        aa_plt_graph_pre = aa_plt_graph_pre.sort_values(
            by=[self.MARKER_NAME, self.MARKER_NAME_NUMBER, self.MUT_TYPE, self.AMINO_ACIDS_NUM],
            ascending=[True, False, True, True]
        )
        return aa_plt_graph_pre

    def assign_mutation_type(self, amino_df):
        # If not short-variant, no further processing is performed
        if amino_df[self.TYPE] != self.SNV:
            return None
        return super().assign_mutation_type(amino_df)

    def sort_mutation_type(self, mut_type):
        return super().sort_mutation_type(mut_type, self.sort_mut_list)

    def create_data_frame(self, amino_df, mut_type):
        return super().create_data_frame(amino_df, mut_type, self.AMINO_ACIDS_NUM)

    def draw_graphical_stem(self, mut_type_df, mut_type):
        return super().draw_graphical_stem(mut_type_df, mut_type, self.draw_mut_list)

    def draw_mutation_type_legend(self, mut_type_list):
        return super().draw_mutation_type_legend(mut_type_list, 4)
