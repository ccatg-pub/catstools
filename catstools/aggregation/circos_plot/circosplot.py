#!/usr/bin/env python
# coding: utf-8

import os

import numpy as np
import pandas as pd
from ..aggregation import Aggregation


class CircosPlot(Aggregation):
    # various consts
    # data frame column name(for input data)

    SNP_SNV = 'SNP_SNV'
    CNV_AMP = 'CNV_Amp'
    CNV_LOSS = 'CNV_Loss'
    LINK_FUSION = 'LINK_Fusion'
    TXT_EXTENSION = '.txt'
    JS_EXTENSION = '.js'

    VALUE = 'Value'
    VARIANT_TYPE = 'Type'
    SNV_CHANGE = 'SNV_change'
    DESCRIPTION = 'description'
    FUS_GENE1 = 'Gene1'
    FUS_GENE2 = 'Gene2'
    NOT_CHR = 'not Chromosome1=="-"'

    def draw_circos_plot(self, input_data_rawdata, output_for_circos):
        # Import data
        panel_df = pd.read_csv(input_data_rawdata, sep=self.SEP_TAB)
        transcript_df = pd.read_csv(self.get_input_tid(), sep=self.SEP_TAB)

        # Check if the number of data frame length is not zero
        if len(panel_df) == 0:
            print('No variant data in input data')
            return None

        # Extraction of sample IDs containing SNV, Amp, Loss and Fusion
        cir_plt_base = panel_df.query(self.ORIGIN_SOMATIC).query(self.CONNECT_OR.join(
            [self.TYPE_SNV, self.TYPE_AMP, self.TYPE_LOSS, self.TYPE_FUSION])).sort_values([self.SAMPLE_ID])

        target_cir_plt = self.exclude_non_target_gene(cir_plt_base, transcript_df)
        cir_plt = self.set_gene_default_position(target_cir_plt, transcript_df)

        cir_plt[self.CHROMOSOME1] = cir_plt[self.CHROMOSOME1].astype(str).str.replace(' ', '')
        cir_plt[self.START_POS1] = cir_plt[self.START_POS1].astype(str).str.replace(' ', '')
        cir_plt[self.END_POS1] = cir_plt[self.END_POS1].astype(str).str.replace(' ', '')
        cir_plt[self.VALUE] = cir_plt[self.VALUE].astype(str).str.replace(' ', '')
        cir_plt[self.CHROMOSOME2] = cir_plt[self.CHROMOSOME2].astype(str).str.replace(' ', '')
        cir_plt[self.START_POS2] = cir_plt[self.START_POS2].astype(str).str.replace(' ', '')
        cir_plt[self.END_POS2] = cir_plt[self.END_POS2].astype(str).str.replace(' ', '')

        # SNV
        self.create_snv(output_for_circos, cir_plt)
        # Amp
        self.create_amp(output_for_circos, cir_plt)
        # Loss
        self.create_loss(output_for_circos, cir_plt)
        # Fusion(Link)
        self.create_fusion(output_for_circos, cir_plt)

    def get_input_tid(self) -> str:
        # get main.py current path
        main_py_path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(main_py_path, "input_data", 'transcriptionID.tsv')

    def create_snv(self, output_for_circos, cir_plt):
        snp_snv = os.path.join(output_for_circos, self.SNP_SNV + self.TXT_EXTENSION)
        snv_circos_pre = cir_plt.query(self.TYPE_SNV).loc[:, [
            self.MARKER_NAME, self.VARIANT_TYPE, self.CHROMOSOME1, self.START_POS1, self.VALUE, self.AMINO_ACIDS_CHANGE,
            self.CDS_CHANGE]]

        if len(snv_circos_pre) == 0:
            with open(snp_snv, mode="w", encoding="utf-8") as snv_writer:
                snv_writer.write("")

            return None

        snv_circos_pre[self.SNV_CHANGE] = snv_circos_pre.apply(lambda x: self.choose_snv_change(x), axis=1)
        snv_circos_pre = snv_circos_pre.dropna(how='any')
        snv_circos_pre[self.DESCRIPTION] = snv_circos_pre[self.MARKER_NAME] + \
            ' ' + snv_circos_pre[self.SNV_CHANGE]
        # Renaming
        snv_circos = snv_circos_pre.loc[:, [self.CHROMOSOME1, self.START_POS1, self.VALUE,
                                            self.DESCRIPTION]].sort_values([self.CHROMOSOME1])
        # Saving
        snv_circos.to_csv(snp_snv, sep=self.SEP_TAB, header=None, index=False)

    def create_amp(self, output_for_circos, cir_plt):
        cnv_amp = os.path.join(output_for_circos, self.CNV_AMP + self.TXT_EXTENSION)
        amp_circos_pre = cir_plt.query(self.TYPE_AMP).query(self.NOT_CHR).loc[:, [
            self.MARKER_NAME, self.VARIANT_TYPE, self.CHROMOSOME1, self.START_POS1, self.END_POS1, self.VALUE]]
        # Renaming
        amp_circos = amp_circos_pre.loc[:, [self.CHROMOSOME1, self.START_POS1, self.END_POS1,
                                            self.VALUE]].sort_values([self.CHROMOSOME1])
        # Saving
        amp_circos.to_csv(cnv_amp, sep=self.SEP_TAB, header=None, index=False)

    def create_loss(self, output_for_circos, cir_plt):
        cnv_loss = os.path.join(output_for_circos, self.CNV_LOSS + self.TXT_EXTENSION)
        loss_circos_pre = cir_plt.query(self.TYPE_LOSS).query(self.NOT_CHR).loc[:, [
            self.MARKER_NAME, self.VARIANT_TYPE, self.CHROMOSOME1, self.START_POS1, self.END_POS1, self.VALUE]]
        # Renaming
        loss_circos = loss_circos_pre.loc[:, [self.CHROMOSOME1, self.START_POS1, self.END_POS1,
                                              self.VALUE]].sort_values([self.CHROMOSOME1])
        # Saving
        loss_circos.to_csv(cnv_loss, sep=self.SEP_TAB, header=None, index=False)

    def create_fusion(self, output_for_circos, cir_plt):
        link_fusion = os.path.join(output_for_circos, self.LINK_FUSION + self.TXT_EXTENSION)
        fusion_circos_pre = cir_plt.query(self.TYPE_FUSION).loc[:, [
            self.MARKER_NAME, self.VARIANT_TYPE, self.CHROMOSOME1, self.START_POS1, self.END_POS1, self.CHROMOSOME2,
            self.START_POS2, self.END_POS2]]
        fusion_circos_pre[self.FUS_GENE1] = fusion_circos_pre.Marker_Name.str.extract('(.+)-', expand=True)
        fusion_circos_pre[self.FUS_GENE2] = fusion_circos_pre.Marker_Name.str.extract('-(.+)', expand=True)
        fusion_circos_pre[self.CHROMOSOME1] = fusion_circos_pre[self.CHROMOSOME1].apply(self.change_chromosome)
        fusion_circos_pre[self.CHROMOSOME2] = fusion_circos_pre[self.CHROMOSOME2].apply(self.change_chromosome)
        fusion_circos_pre = fusion_circos_pre.dropna(how='any')
        # Renaming
        fusion_circos = fusion_circos_pre.loc[:, [self.MARKER_NAME, self.CHROMOSOME1, self.START_POS1, self.END_POS1,
                                                  self.FUS_GENE1, self.CHROMOSOME2, self.START_POS2, self.END_POS2,
                                                  self.FUS_GENE2]].sort_values([self.CHROMOSOME1])
        # Saving
        fusion_circos.to_csv(link_fusion, sep=self.SEP_TAB, header=None, index=False)

    def choose_snv_change(self, circos_df):
        if circos_df.Amino_Acids_Change == '-' and circos_df.Cds_Change != '-':
            return circos_df.Cds_Change
        if circos_df.Amino_Acids_Change != '-':
            return circos_df.Amino_Acids_Change
        else:
            return np.NaN

    def change_chromosome(self, chromosome_num):
        # Changing chromosome
        if chromosome_num in {'1', '2', '3', '4', '5', '5', '6', '7', '8', '9',
                              '10', '11', '12', '13', '14', '15', '16', '17',
                              '18', '19', '20', '21', '22', 'X'}:
            return chromosome_num
        else:
            return np.nan
