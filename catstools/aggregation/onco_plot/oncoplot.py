# this documentation makes use of pandas, numpy, and palettable
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from comut import comut

from catstools.aggregation.aggregation import Aggregation


class OncoPlot(Aggregation):

    # various consts
    # data frame column name(for input data)
    MUTATION_TYPE = 'Mutation type'
    SAMPLE = 'sample'
    SEX = 'sex'
    VALUE = 'value'
    CATEGORY = 'category'
    AGE = 'age'
    CANCER_TYPE = 'cancer_type'
    TMB = 'TMB'
    DISEASE = 'Disease'
    GENDER = 'Gender'
    MALE = '男'
    FEMALE = '女'
    GENDER_MALE = 'Male'
    GENDER_FEMALE = 'Female'
    GENDER_OTHER = '-'
    ID_RANK = "ID rank"
    MUTATION_COUNT = "Mutation Count"

    # variant type
    ALL = 'All'
    MULTIPLE = 'Multiple'

    # various option value
    WIDTH = 'width'
    EDGE_COLOR = 'edgecolor'
    MUT_FREQ = 'Mutfreq'
    DEL_LEGEND_NAME = 'Mutationtype'

    # disease name
    ADRENAL_GLAND = 'Adrenal Gland'
    AMPULLA_OF_VATER = 'Ampulla of Vater'
    BILIARY_TRACT = 'Biliary Tract'
    BLADDER_URINARY_TRACT = 'Bladder/Urinary Tract'
    BONE = 'Bone'
    BOWEL = 'Bowel'
    BREAST = 'Breast'
    CNS_BRAIN = 'CNS/Brain'
    CERVIX = 'Cervix'
    ESOPHAGUS_STOMACH = 'Esophagus/Stomach'
    EYE = 'Eye'
    HEAD_AND_NECK = 'Head and Neck'
    KIDNEY = 'Kidney'
    LIVER = 'Liver'
    LUNG = 'Lung'
    OTHER = 'Other'
    OVARY_FALLOPIAN_TUBE = 'Ovary/Fallopian Tube'
    PANCREAS = 'Pancreas'
    PENIS = 'Penis'
    PERIPHERAL_NERVOUS_SYSTEM = 'Peripheral Nervous System'
    PERITONEUM = 'Peritoneum'
    PLEURA = 'Pleura'
    PROSTATE = 'Prostate'
    SKIN = 'Skin'
    SOFT_TISSUE = 'Soft Tissue'
    TESTIS = 'Testis'
    THYMUS = 'Thymus'
    THYROID = 'Thyroid'
    UTERUS = 'Uterus'
    VULVA_VAGINA = 'Vulva/Vagina'

    # graph color
    GREEN = 'green'
    BROWN = 'brown'
    PURPLE = 'purple'
    YELLOW = 'yellow'
    MEDIUM_AQUA_MARINE = 'mediumaquamarine'
    SLATE_BLUE = 'slateblue'
    DARK_ORANGE = 'darkorange'
    YELLOW_GREEN = 'yellowgreen'
    DARK_MAGENTA = 'darkmagenta'
    WHITE_SMOKE = 'whitesmoke'
    SADDLE_BROWN = 'saddlebrown'
    HOT_PINK = 'hotpink'
    GRAY = 'gray'
    TEAL = 'teal'
    LIGHT_SKY_BLUE = 'lightskyblue'
    DARK_GREEN = 'darkgreen'
    DARK_RED = 'darkred'
    MEDIUM_SEA_GREEN = 'mediumseagreen'
    LIGHT_BLUE = 'lightblue'
    INDIGO = 'indigo'
    GREY = 'grey'
    FOREST_GREEN = 'forestgreen'
    CYAN = 'cyan'
    DIM_GREY = 'dimgrey'
    LIGHT_YELLOW = 'lightyellow'
    DARK_SLATE_BLUE = 'darkslateblue'
    DARK_CYAN = 'darkcyan'
    PEACH_PUFF = 'peachpuff'
    NAVY = 'navy'
    INNER = 'inner'
    HORIZONTAL = 'horizontal'
    MUTS_MB = 'Muts/Mb'
    MUTATION_BURDEN = 'Mutation burden'

    # query
    CATEGORY_NOT_ALL = 'category != "All"'
    SAMPLE_NOT_ALL = 'sample != "All"'

    # output file name
    AGD_SVG = 'AGD.svg'
    TMB_SVG = 'TMB.svg'
    ONCOGENE_X_SVG = 'oncogene_x.svg'
    ONCOGENE_XY_SVG = 'oncogene_xy.svg'
    ONCOGENE_X2_SVG = 'oncogene_x2.svg'

    def draw_oncoplot(self, input_data_rawdata, output_for_oncoplot):
        # Data frames for aggregation
        # (top 100 with the highest number of mutations)
        target_mut_df = self.create_target_df(input_data_rawdata)

        # Check if the number of data frame length is not zero
        if len(target_mut_df) == 0:
            print('No variant data in input data')
            return None

        # TMB Data Frame
        tmb_df = self.get_tmb_df(target_mut_df)
        # Frequency
        mut_df_pre = self.get_mut_df_pre(target_mut_df)
        mut_freq = self.get_mut_freq(mut_df_pre, self.CATEGORY, self.CATEGORY_NOT_ALL)

        val_types = [self.SNV, self.DELETION, self.INSERTION, self.DELINS, self.AMPLIFICATION,
                     self.LOSS, self.FUSION]
        mut_freq_ed = mut_freq.loc[:, [self.CATEGORY] + val_types]

        # mut_freq
        mut_df_x2 = self.get_mut_freq(mut_df_pre, self.SAMPLE, self.SAMPLE_NOT_ALL)
        mut_df_x3 = mut_df_x2.loc[:, [self.SAMPLE] + val_types]
        mut_df = self.get_mut_df(mut_df_pre)

        # Age, Disease, Gender
        self.draw_patient_information_graph(mut_df_pre, mut_df, output_for_oncoplot)
        # TMB
        self.draw_tmb_graph(tmb_df, output_for_oncoplot)

        # define mapping, shrink
        bar_mapping = self.get_bar_mapping()

        # oncogene-x
        self.draw_oncogene_x_graph(mut_df, mut_freq_ed, mut_df_x3, bar_mapping, output_for_oncoplot)
        # oncogene-xy
        mut_df_x_sort, mut_freq_y_sort = \
            self.draw_oncogene_xy_graph(mut_df, mut_freq, mut_freq_ed, mut_df_x2,
                                        mut_df_x3, bar_mapping, output_for_oncoplot)
        # oncogene-x2
        self.draw_oncogene_x2_graph(mut_freq_ed, mut_df_x_sort, mut_freq_y_sort,
                                    mut_df_x3, bar_mapping, output_for_oncoplot)

    def create_target_df(self, input_data_rawdata):
        """
        Generate data frames for the top 100 most mutated cases.
        If the number of mutations is the same,
        priority is given to the one with the younger ID number.
        TMB and MSI shall not be included in the count of mutations.
        """
        # import data
        panel_df = pd.read_csv(input_data_rawdata, sep=self.SEP_TAB)

        # Exclude TMB and MSI information
        mut_df = panel_df.query(self.CONNECT_AND.join([self.NOT_TYPE_MSI, self.NOT_TYPE_TMB]))
        # Data frames prioritized by ID
        mut_df[self.ID_RANK] = mut_df[self.SAMPLE_ID].rank()
        id_rank_df = mut_df.loc[:, [self.SAMPLE_ID, self.ID_RANK]]
        id_rank_df = id_rank_df.drop_duplicates()
        # Data frame containing the number of mutations
        base_mut_cnt_ser = mut_df[self.SAMPLE_ID].value_counts()
        # Series to Data frame
        base_mut_cnt_df = base_mut_cnt_ser.rename_axis(self.SAMPLE_ID)
        base_mut_cnt_df = base_mut_cnt_df.reset_index(name=self.MUTATION_COUNT)
        # Combine the data frame storing the number of mutations
        # and the data frame storing the priority of IDs.
        mut_cnt_df = pd.merge(base_mut_cnt_df, id_rank_df, on=self.SAMPLE_ID, how=self.INNER)
        mut_cnt_df = mut_cnt_df.sort_values([self.MUTATION_COUNT, self.ID_RANK], ascending=[False, True])
        mut_cnt_df = mut_cnt_df.head(100)
        # Generate data frames with only the information
        # from the top 100 most mutated cases.
        target_df = pd.merge(panel_df, mut_cnt_df, on=self.SAMPLE_ID, how=self.INNER)

        return target_df

    def get_tmb_df(self, target_df):
        # import data
        base_tmb_df = target_df.query(self.TYPE_TMB)
        # Extract the top 100 sample IDs
        excerpt_tmb_df = base_tmb_df[self.SAMPLE_ID].drop_duplicates().head(100)
        # Focus on 100 sample IDs
        tmb_df = pd.merge(base_tmb_df, excerpt_tmb_df, on=self.SAMPLE_ID, how=self.INNER)
        # Convert TMB values to float type
        # to correctly reflect y-axis values in graphs
        tmb_df[self.VALUE_COL] = tmb_df[self.VALUE_COL].astype(float)
        # column renaming
        rename_columns = {self.SAMPLE_ID: self.SAMPLE}

        return tmb_df.rename(columns=rename_columns)

    def get_mut_df(self, mut_df_pre):

        drop_columns = [self.PANEL_NAME, self.ORIGIN_TYPE, self.VARIANT_TYPE_STATE, self.REFERENCE_ALLELE,
                        self.ALTERNATE_ALLELE, self.CHROMOSOME1, self.START_POS1, self.END_POS1, self.TYPE,
                        self.CHROMOSOME2, self.START_POS2, self.END_POS2, self.REFERENCE_ALLELE_LEN,
                        self.ALTERNATE_ALLELE_LEN, self.REFERENCE_ALLELE_1STL, self.ALTERNATE_ALLELE_1STL,
                        self.RA_TF]

        return mut_df_pre.drop(columns=drop_columns)

    def get_mut_freq(self, mut_df_pre, in_key, query_key):
        # Frequency
        mut_df_sample = mut_df_pre.loc[:, [in_key, self.VALUE]]
        mut_freq = pd.pivot_table(mut_df_sample, index=[in_key], columns=[self.VALUE], aggfunc=len,
                                  fill_value=0, margins=True)
        mut_freq = mut_freq.query(query_key).reset_index()

        # Actions that KeyError may occur
        return self.check_mutation_type(mut_freq)

    def get_mut_df_pre(self, target_df):
        # Extraction of sample IDs containing SNV, Amp, Loss and Fusion
        panel_df_l = target_df.query(self.CONNECT_AND.join([self.NOT_TYPE_MSI, self.NOT_TYPE_TMB]))
        # SNVs to SNV, Deletion,Insertion, Delins
        df_type = self.set_allele_info(panel_df_l)
        rename_columns = {self.SAMPLE_ID: self.SAMPLE,
                          self.MUT_TYPE: self.VALUE,
                          self.MARKER_NAME: self.CATEGORY}
        return df_type.rename(columns=rename_columns)

    def get_bar_mapping(self):
        bar_mapping = {
            self.SNV: self.MEDIUM_AQUA_MARINE,
            self.INSERTION: self.SLATE_BLUE,
            self.DELETION: self.DARK_ORANGE,
            self.DELINS: self.YELLOW_GREEN,
            self.AMPLIFICATION: self.RED,
            self.LOSS: self.BLUE,
            self.FUSION: self.GREEN,
            self.MULTIPLE: self.BROWN
        }
        return bar_mapping

    def check_mutation_type(self, in_df):
        df_map = in_df
        val_types = [self.SNV, self.DELETION, self.INSERTION, self.DELINS, self.AMPLIFICATION,
                     self.LOSS, self.FUSION]

        for val_type in val_types:
            df_map[val_type] = self.check_mutation_type_exists(df_map, val_type)

        return df_map

    def draw_patient_information_graph(self, mut_df_pre, mut_df, output_dir):
        """
        Draw a graph about patient information(age, disease, gender).
        """
        if any([self.AGE not in mut_df_pre.columns,
                self.CANCER_TYPE not in mut_df_pre.columns,
                self.SEX not in mut_df_pre.columns]):
            print('no Age, Gender, Disease data')
        else:
            # age
            mut_df_age_ed = self.create_patient_df_age(mut_df)

            # Types of Cancer
            mut_df_disease_ed = self.create_patient_df_disease(mut_df)
            # Determine color by cancer type
            disease_color = {
                self.ADRENAL_GLAND: self.PURPLE,
                self.AMPULLA_OF_VATER: self.DARK_MAGENTA,
                self.BILIARY_TRACT: self.GREEN,
                self.BLADDER_URINARY_TRACT: self.YELLOW,
                self.BONE: self.WHITE_SMOKE,
                self.BOWEL: self.SADDLE_BROWN,
                self.BREAST: self.HOT_PINK,
                self.CNS_BRAIN: self.GRAY,
                self.CERVIX: self.TEAL,
                self.ESOPHAGUS_STOMACH: self.LIGHT_SKY_BLUE,
                self.EYE: self.DARK_GREEN,
                self.HEAD_AND_NECK: self.DARK_RED,
                self.KIDNEY: self.ORANGE,
                self.LIVER: self.MEDIUM_SEA_GREEN,
                self.LUNG: self.GAINSBORO,
                self.OTHER: self.BLACK,
                self.OVARY_FALLOPIAN_TUBE: self.LIGHT_BLUE,
                self.PANCREAS: self.INDIGO,
                self.PENIS: self.BLUE,
                self.PERIPHERAL_NERVOUS_SYSTEM: self.GREY,
                self.PERITONEUM: self.FOREST_GREEN,
                self.PLEURA: self.MEDIUM_BLUE,
                self.PROSTATE: self.CYAN,
                self.SKIN: self.DIM_GREY,
                self.SOFT_TISSUE: self.LIGHT_YELLOW,
                self.TESTIS: self.RED,
                self.THYMUS: self.DARK_SLATE_BLUE,
                self.THYROID: self.DARK_CYAN,
                self.UTERUS: self.PEACH_PUFF,
                self.VULVA_VAGINA: self.NAVY
            }

            # gender
            mut_df_sex_ed = self.create_patient_df_gender(mut_df)

            # add other data
            # make CoMut
            cat_comut = comut.CoMut()

            # Types of Cancer(write a graph third from the top)
            cat_comut.add_categorical_data(mut_df_disease_ed, mapping=disease_color, name=self.DISEASE)
            # gender(Write a graph second from the top)
            # Determining color of gender
            gender_color = {self.GENDER_MALE: self.BLUE, self.GENDER_FEMALE: self.RED}
            cat_comut.add_categorical_data(mut_df_sex_ed, mapping=gender_color, name=self.GENDER)

            # age(write the graph at the top)
            # Determining color
            cmap = mpl.cm.jet
            cat_comut.add_continuous_data(mut_df_age_ed, name=self.AGE, mapping=cmap, value_range=(0, 1))

            # drawing
            cat_comut.plot_comut(figsize=(20, 1), x_padding=0.04, y_padding=0.04,
                                 tri_padding=0.03, hspace=0.08, subplot_hspace=0.01)

            # Add a legend
            # age
            # color bars must be added manually based on figure coordinates
            # - [left, bottom, width, height]
            age_ax = cat_comut.figure.add_axes([1, 1.2, 0.08, 0.2])
            # Age ranges from 0 to 1
            norm = mpl.colors.Normalize(vmin=0, vmax=1)
            # create the colorbar with colormap used
            # when the continuous data was added (purp_7)
            age_colorbar = cat_comut.figure.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=age_ax, orientation=self.HORIZONTAL
            )
            # remove tick marks and move tick labels slightly inwards.
            # Also remove black border
            age_colorbar.ax.tick_params(size=0)
            # Decide where to display
            age_colorbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])
            # Specify the value to display
            age_colorbar.set_ticklabels([0, 20, 40, 60, 80, 100])
            age_colorbar.outline.set_visible(False)
            # set title of color bar to line up with other legend elements
            # Decide where to display the legend (displayed above)
            age_colorbar.set_label(self.AGE, labelpad=-50, x=0, fontsize=14)

            # Determine the order of the 30 types of cancer
            disease_order = (
                self.ADRENAL_GLAND, self.AMPULLA_OF_VATER, self.BILIARY_TRACT, self.BLADDER_URINARY_TRACT, self.BONE,
                self.BOWEL, self.BREAST, self.CNS_BRAIN, self.CERVIX, self.ESOPHAGUS_STOMACH, self.EYE,
                self.HEAD_AND_NECK, self.KIDNEY, self.LIVER, self.LUNG, self.OTHER, self.OVARY_FALLOPIAN_TUBE,
                self.PANCREAS, self.PENIS, self.PERIPHERAL_NERVOUS_SYSTEM, self.PERITONEUM, self.PLEURA, self.PROSTATE,
                self.SKIN, self.SOFT_TISSUE, self.TESTIS, self.THYMUS, self.THYROID, self.UTERUS, self.VULVA_VAGINA
            )
            # Order of gender
            gender_order = self.GENDER_MALE, self.GENDER_FEMALE

            # gender
            cat_comut.add_axis_legend(name=self.GENDER, bbox_to_anchor=(1.1, 6), order=gender_order,
                                      title=self.GENDER, ncol=1)
            # Types of Cancer
            cat_comut.add_axis_legend(name=self.DISEASE, bbox_to_anchor=(1, 3), order=disease_order,
                                      title=self.DISEASE, ncol=2)

            # save Age,Gender and Disease
            agd_save = os.path.join(output_dir, self.AGD_SVG)
            cat_comut.figure.savefig(agd_save, bbox_inches=self.TIGHT, dpi=600)
            print('Age, Gender, Disease plot')

            # Open memory
            plt.clf()
            # Close a saved graph
            plt.close()

    def draw_tmb_graph(self, tmb_df, output_dir):
        """
        Draw a graph about TMB information.
        """
        # Check if the number of TMB data frame length is not zero
        if len(tmb_df) == 0:
            print('No TMB data in input data')
            return None
        else:
            # TMB
            mut_df_tmb_ed = tmb_df.loc[:, [self.SAMPLE, self.VALUE_COL]]
            # define mapping, shrink
            bar_mapping = {self.VALUE_COL: self.PURPLE}
            bar_kwargs = {self.WIDTH: 0.8, self.EDGE_COLOR: self.BLACK}

            # make CoMut
            cat_comut = comut.CoMut()

            # add other data
            # TMB
            cat_comut.add_bar_data(mut_df_tmb_ed, name=self.MUTATION_BURDEN, mapping=bar_mapping,
                                   stacked=True, bar_kwargs=bar_kwargs, ylabel=self.MUTS_MB)
            # drawing
            cat_comut.plot_comut(figsize=(30, 10), x_padding=0.04, y_padding=0.04,
                                 tri_padding=0.03, hspace=0.08, subplot_hspace=0.01)
            cat_comut.add_unified_legend()

            # Save TMB
            tmb_save = os.path.join(output_dir, self.TMB_SVG)
            cat_comut.figure.savefig(tmb_save, bbox_inches=self.TIGHT, dpi=600)
            print('TMB plot')

            # Open memory
            plt.clf()
            # Close a saved graph
            plt.close()

    def draw_oncogene_x_graph(self, mut_df, mut_freq_df, gene_df, color_map, output_dir):

        # graph save
        save_fig = os.path.join(output_dir, self.ONCOGENE_X_SVG)
        print('Now outputting the oncogene plot (1/3).')

        # draw once graph
        self.draw_onco_graph(mut_df, None, mut_freq_df, gene_df, color_map, save_fig)
        print('Output of oncogene plot1 completed (1/3).')

    def draw_oncogene_xy_graph(self, mut_df, mut_freq_df, mut_freq_ed_df, mut_df_x2,
                               mut_df_x3, color_map, output_dir):

        # graph save
        save_fig = os.path.join(output_dir, self.ONCOGENE_XY_SVG)

        # Cancer mutations
        # Sort the Y-axis(alphabetical order)
        mut_freq_y_sort_alphabet = self.sort_loc_values(mut_freq_df, [self.CATEGORY], [self.CATEGORY], [True])

        # Prepare before it draw (decide on sorting)
        print('Now outputting the oncogene plot (2/3).')

        # draw onco graph
        self.draw_onco_graph(mut_df, list(reversed(mut_freq_y_sort_alphabet)),
                             mut_freq_ed_df, mut_df_x3, color_map, save_fig)

        sort_values = [self.ALL, self.SNV, self.DELETION, self.INSERTION, self.DELINS,
                       self.AMPLIFICATION, self.LOSS]
        ascendings = [False, False, False, False, False, False, False, True]

        # graph save
        save_fig = os.path.join(output_dir, self.ONCOGENE_XY_SVG)

        # Cancer mutations
        # Sort the Y-axis(alphabetical order)
        mut_freq_y_sort = self.sort_loc_values(mut_freq_df, sort_values + [self.CATEGORY], [self.CATEGORY], ascendings)

        # draw onco graph
        self.draw_onco_graph(mut_df, list(reversed(mut_freq_y_sort)), mut_freq_ed_df, mut_df_x3, color_map, save_fig)

        # When X-axis is sorted
        # Prepare for sorting of number of X-axis
        # Make sample IDs lists
        mut_freq_x_sort = self.sort_loc_values(mut_df_x2, sort_values + [self.SAMPLE], [self.SAMPLE], ascendings)

        # Reorder the data (order of sorted sample IDs)
        mut_df[self.SAMPLE] = pd.Categorical(mut_df[self.SAMPLE], categories=mut_freq_x_sort)
        mut_df_x_sort = mut_df.sort_values(by=[self.SAMPLE])

        # graph save
        oncogene_xy_save = os.path.join(output_dir, self.ONCOGENE_XY_SVG)

        # draw onco graph
        self.draw_onco_graph(mut_df_x_sort, list(reversed(mut_freq_y_sort_alphabet)),  mut_freq_ed_df, mut_df_x3,
                             color_map, oncogene_xy_save)

        print('Output of oncogene plot completed (2/3).')

        return mut_df_x_sort, mut_freq_y_sort

    def draw_onco_graph(self, cat_data, cat_order, side_bar_data, bar_data, color_map, save_fig):

        # make CoMut
        cat_comut = comut.CoMut()

        # Prepare before it draw (decide on sorting)
        cat_comut.add_categorical_data(cat_data, name=self.MUTATION_TYPE, mapping=color_map,
                                       category_order=cat_order)

        # Y-axis Frequency
        # sort
        cat_comut.add_side_bar_data(side_bar_data, paired_name=self.MUTATION_TYPE, name=self.MUT_FREQ,
                                    mapping=color_map, xlabel=self.FREQUENCY, stacked=True, position=self.LEFT)
        # Types of mutations of Axis X Frequency
        # sort
        cat_comut.add_bar_data(bar_data, stacked=True, name=self.DEL_LEGEND_NAME,
                               ylabel=self.FREQUENCY, mapping=color_map)
        # draw graphs and organize legends
        self.organize_oncoplot(cat_comut)
        # graph save
        cat_comut.figure.savefig(save_fig, bbox_inches=self.TIGHT, dpi=600)

        # Open memory
        plt.clf()
        # Close a saved graph
        plt.close()

    def sort_loc_values(self, in_df, sort_values, loc_values, ascendings):
        ret_df = in_df.sort_values(sort_values, ascending=ascendings)
        ret_df = ret_df.loc[:, loc_values].values.tolist()
        return sum(ret_df, [])

    def draw_oncogene_x2_graph(self, mut_freq_df, x_sort_df, y_sort_df, gene_df,
                               color_map, output_dir):
        # graph save
        save_fig = os.path.join(output_dir, self.ONCOGENE_X2_SVG)
        print('Now outputting the oncogene plot (3/3).')

        # draw onco graph
        self.draw_onco_graph(x_sort_df, list(reversed(y_sort_df)), mut_freq_df, gene_df, color_map, save_fig)
        print('Output of oncogene plot completed (3/3).')

    def convert_gender(self, gender):
        gender_dict = {self.MALE: self.GENDER_MALE, self.FEMALE: self.GENDER_FEMALE}

        return gender_dict.get(gender, self.GENDER_OTHER)

    def check_mutation_type_exists(self, mut_df, df_column):
        """
        Check for the presence of a column for each mutation type.
        """
        # If the mutation type column does not exist
        if df_column not in mut_df.columns:
            return 0
        else:
            return mut_df[df_column]

    def create_patient_df_age(self, mut_df):
        """
        Create a data frame regarding patient information(age).
        """
        patient_df = mut_df.loc[:, [self.SAMPLE, self.AGE]]
        patient_df[self.CATEGORY] = self.AGE
        patient_df[self.VALUE] = patient_df[self.AGE] * 0.01

        return patient_df.loc[:, [self.SAMPLE, self.CATEGORY, self.VALUE]].drop_duplicates()

    def create_patient_df_disease(self, mut_df):
        """
        Create a data frame regarding patient information(disease).
        """
        patient_df = mut_df.loc[:, [self.SAMPLE, self.CANCER_TYPE]]
        patient_df[self.CATEGORY] = self.DISEASE
        patient_df[self.VALUE] = patient_df[self.CANCER_TYPE]

        return patient_df.loc[:, [self.SAMPLE, self.CATEGORY, self.VALUE]].drop_duplicates()

    def create_patient_df_gender(self, mut_df):
        """
        Create a data frame regarding patient information(gender).
        """
        patient_df = mut_df.loc[:, [self.SAMPLE, self.SEX]]
        patient_df[self.CATEGORY] = self.GENDER
        patient_df[self.VALUE] = patient_df[self.SEX].apply(self.convert_gender)

        return patient_df.loc[:, [self.SAMPLE, self.CATEGORY, self.VALUE]].drop_duplicates()

    def organize_oncoplot(self, comut_obj):
        """
        Draw graphs and organize legends.
        """
        mut_order = [self.SNV, self.INSERTION, self.DELETION, self.DELINS,
                     self.AMPLIFICATION, self.LOSS, self.FUSION, self.MULTIPLE]
        # drawing
        comut_obj.plot_comut(figsize=(30, 50), x_padding=0.04, y_padding=0.04,
                             tri_padding=0.03, hspace=0.08, subplot_hspace=0.01)
        # remove unified legend
        comut_obj.add_unified_legend()
        if comut_obj.axes[self.DEL_LEGEND_NAME].get_legend() is not None:
            comut_obj.axes[self.DEL_LEGEND_NAME].get_legend().remove()
        # Draw legends again
        comut_obj.add_axis_legend(name=self.MUTATION_TYPE, bbox_to_anchor=(1, 0.85),
                                  title=self.MUTATION_TYPE, order=mut_order, ncol=1)
