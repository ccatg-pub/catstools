#!/usr/bin/env python
# coding: utf-8
import os
import sys

# import convert_cats_data
from . import convert_cats_data

from .circos_plot.circosplot import CircosPlot
from .gene_analysis_summary.gene_analysis_plot_amino_acids import GeneAnalysisPlotAminoAcids
from .gene_analysis_summary.gene_analysis_plot_dna_cna import GeneAnalysisPlotDnaCna
from .gene_analysis_summary.gene_analysis_plot_dna_except_fusion import GeneAnalysisPlotExceptFusion
from .gene_analysis_summary.gene_analysis_plot_dna_fusion import GeneAnalysisPlotDnaFusion
from .onco_plot.oncoplot import OncoPlot

# Argument length
MAIN_LENGTH = 4
# Error messages
SELECT_MODE_MESSAGE = "Please select graph mode from 1 to 6.\n" \
    "graph mode [1]: onco_plot\n" \
    "graph mode [2]: circos_plot\n" \
    "graph mode [3]: gene_analysis_plot_amino\n" \
    "graph mode [4]: gene_analysis_plot_dna_except_fusion\n" \
    "graph mode [5]: gene_analysis_plot_cna\n" \
    "graph mode [6]: gene_analysis_plot_fusion"
ARGUMENT_1_MESSAGE = "Argument 1 requires a number that represents the graph aggregation type."
ARGUMENT_2_MESSAGE = "Argument 2 requires a search pattern for the CATS format path to be used for graph aggregation."
ARGUMENT_3_MESSAGE = "Argument 3 requires the path of the output data directory to save the aggregate graph."


def run(graph_mode, input_cats_pattern, output_dir):

    # argument integrity flag
    arg_check_result = True

    # Argument consistency check
    # Graph mode consistency check
    if not graph_mode.isdigit() or int(graph_mode) not in range(1, 7):
        print(f"Unexpect graph mode: {graph_mode}")
        print(SELECT_MODE_MESSAGE)
        arg_check_result = False

    # Output directory path consistency check
    if not os.path.isdir(output_dir):
        print(ARGUMENT_3_MESSAGE)
        arg_check_result = False

    # If arguments are inconsistent, terminate processing
    if not arg_check_result:
        sys.exit(-1)

    # Generate TSV files based on the CATS format.
    tsv_file = convert_cats_data.create_tsv_file(input_cats_pattern, output_dir)

    try:
        # Draw onco_plot
        if graph_mode == "1":
            oncoplot = OncoPlot()
            oncoplot.draw_oncoplot(tsv_file, output_dir)
        # Draw circos_plot
        elif graph_mode == "2":
            circosplot = CircosPlot()
            circosplot.draw_circos_plot(tsv_file, output_dir)
        # Draw aggregate graph by gene(protein)
        elif graph_mode == "3":
            gene_amino = GeneAnalysisPlotAminoAcids()
            gene_amino.draw_gene_analysis_amino(tsv_file, output_dir)
        # Draw aggregate graph by gene(dna except fusion)
        elif graph_mode == "4":
            gene_except_fusion = GeneAnalysisPlotExceptFusion()
            gene_except_fusion.draw_gene_analysis_except_fusion(tsv_file, output_dir)
        # Draw aggregate graph by gene(cna)
        elif graph_mode == "5":
            gene_cna = GeneAnalysisPlotDnaCna()
            gene_cna.draw_gene_analysis_cna(tsv_file, output_dir)
        # Draw aggregate graph by gene(fusion)
        elif graph_mode == "6":
            gene_fusion = GeneAnalysisPlotDnaFusion()
            gene_fusion.draw_gene_analysis_fusion(tsv_file, output_dir)

    except Exception:
        print("The process was aborted due to an error.")
        raise
