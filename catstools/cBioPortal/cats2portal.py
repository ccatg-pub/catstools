import os
from typing import Dict

from catstools.cBioPortal.variant import (
    filter_fusion, iter_cna, iter_expressions, iter_other_biomarkers, iter_snv, iter_sv)
from catstools.common import CatsMeta


def convert_to_portal(meta: CatsMeta) -> Dict:
    """Conversion from CATS format to cBioPortal input data"""
    output_data: Dict = dict()
    output_data.update(get_meta_study(meta))
    output_data.update(get_case_lists(meta))
    output_data.update(get_clinical_patient(meta))
    output_data.update(get_clinical_sample(meta))
    output_data.update(get_mutations_extended(meta))
    output_data.update(get_cna_discrete(meta))
    output_data.update(get_sv(meta))
    output_data.update(get_fusion(meta))
    output_data.update(get_expression(meta))
    return output_data


def get_meta_study(meta: CatsMeta) -> Dict:
    """Creates a mutation summary file value and returns the file name as a key"""
    grc_release = 'hg38' if 'hg38' in meta.json_cats["metaData"]["referenceGenome"]["grcRelease"] else 'hg19'
    meta_study = ['type_of_cancer: type_of_cancer',
                  f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                  f'name: {meta.json_cats["testInfo"]["testId"]}',
                  'description: sample description',
                  f'reference_genome: {grc_release}']
    return {'meta_study.txt': meta_study}


def get_case_lists(meta: CatsMeta) -> Dict:
    """Create values for various target specimen files and return file names as keys"""
    cases_all = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                 f'stable_id: {meta.json_cats["testInfo"]["testId"]}_all',
                 'case_list_name: All samples',
                 'case_list_description: All samples description',
                 'case_list_category: all_cases_in_study',
                 f'case_list_ids: {meta.json_cats["testInfo"]["testId"]}'
                 ]
    cases_cna = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                 f'stable_id: {meta.json_cats["testInfo"]["testId"]}_cna',
                 'case_list_name: CNA samples',
                 'case_list_description: CNA samples description',
                 'case_list_category: all_cases_in_study',
                 f'case_list_ids: {meta.json_cats["testInfo"]["testId"]}'
                 ]
    cases_cnaseq = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                    f'stable_id: {meta.json_cats["testInfo"]["testId"]}_cnaseq',
                    'case_list_name: CNASEQ samples',
                    'case_list_description: CNASEQ samples description',
                    f'case_list_ids: {meta.json_cats["testInfo"]["testId"]}'
                    ]
    cases_sequenced = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                       f'stable_id: {meta.json_cats["testInfo"]["testId"]}_sequenced',
                       'case_list_name: Sequenced samples',
                       'case_list_description: Sequenced samples description',
                       'case_list_category: all_cases_with_mutation_data',
                       f'case_list_ids: {meta.json_cats["testInfo"]["testId"]}'
                       ]
    return {'case_lists/cases_all.txt': cases_all,
            'case_lists/cases_cna.txt': cases_cna,
            'case_lists/cases_cnaseq.txt': cases_cnaseq,
            'case_lists/cases_sequenced.txt': cases_sequenced
            }


def get_clinical_patient(meta: CatsMeta) -> Dict:
    """Create values for various patient files and return file names as keys"""
    meta_clinical_patient = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                             'genetic_alteration_type: CLINICAL',
                             'datatype: PATIENT_ATTRIBUTES',
                             'data_filename: data_clinical_patient.txt'
                             ]
    data_clinical_patient = ['#Patient Identifier\tSex\tDiagnosis Age\tOverall Survival (Months)\tOverall Survival Status\tDisease Free Status\tDisease Free (Months)',
                             '#Identifier to uniquely specify a patient.\tSex\tAge at which a condition or disease was first diagnosed.\tOverall survival in months since initial diagonosis.\tOverall patient survival status.\tDisease free status\tDisease Free (Months)',
                             '#STRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\tNUMBER',
                             '#1\t1\t1\t1\t9\t1\t1',
                             'PATIENT_ID\tSEX\tAGE\tOS_MONTHS\tOS_STATUS\tDFS_STATUS\tDFS_MONTHS'
                             ]
    return {'meta_clinical_patient.txt': meta_clinical_patient,
            'data_clinical_patient.txt': data_clinical_patient
            }


def get_clinical_sample(meta: CatsMeta) -> Dict:
    """Create values for various specimen files and return file names as keys"""
    meta_clinical_sample = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                            'genetic_alteration_type: CLINICAL',
                            'datatype: SAMPLE_ATTRIBUTES',
                            'data_filename: data_clinical_sample.txt'
                            ]
    data_clinical_sample = ['#Patient Identifier\tSample Identifier\tSample Type\tOncotree Code\tNeoplasm Histologic Grade\tMGMT Status\tIDH1 Mutation\t1p/19q Status\tNon-silent mutations in TP53, ATRX, CIC, FUBP1\tNote\tCancer Type\tCancer Type Detailed\tSomatic Status\tTMB (nonsynonymous)\tMSI (sensor score)',
                            '#Identifier to uniquely specify a patient.\tA unique sample identifier.\tThe type of sample (i.e., normal, primary, met, recurrence).\tOncotree Code\tNumeric value to express the degree of abnormality of cancer cells, a measure of differentiation and aggressiveness.\tMGMT Status\tIDH1 mutation is present\t1p/19q Status\tNon-silent mutations in TP53, ATRX, CIC, FUBP1\tNote.\tCancer Type\tCancer Type Detailed\tSomatic Status\tTMB (nonsynonymous)\tMSI (sensor score)',
                            '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\tNUMBER',
                            '#1\t1\t9\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1',
                            'PATIENT_ID\tSAMPLE_ID\tSAMPLE_TYPE\tONCOTREE_CODE\tGRADE\tMGMT_STATUS\tIDH1_MUTATION\tSTATUS_1P_19Q\tNON_SILENT_MUT_TP53_ATRX_CIC_FUBP1\tNOTE\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tSOMATIC_STATUS\tTMB_NONSYNONYMOUS\tMSI_SENSOR_SCORE'
                            ]
    tmb_value = ''
    msi_value = ''
    output_dir: Dict = dict()
    for other_biomarker in iter_other_biomarkers(meta, filters=[]):
        value = other_biomarker.metrics[0].get('value')
        if other_biomarker.biomarker_type == 'TMB':
            tmb_value = str(value)
        elif other_biomarker.biomarker_type == 'MSI':
            msi_value = str(value)
    if tmb_value or msi_value:
        data_clinical_sample.append('\t'.join(['', meta.json_cats["testInfo"]["testId"], '', '', '', '',
                                               '', '', '', '', '', '', '', tmb_value, msi_value]))
        output_dir = {'meta_clinical_sample.txt': meta_clinical_sample,
                      'data_clinical_sample.txt': data_clinical_sample
                      }

    return output_dir


def get_mutations_extended(meta: CatsMeta) -> Dict:
    """Create values for various snv files and return file names as keys"""
    meta_mutations_extended = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                               'genetic_alteration_type: MUTATION_EXTENDED',
                               'stable_id: mutations',
                               'datatype: MAF',
                               'show_profile_in_analysis_tab: true',
                               'profile_description: mutations_extended',
                               'profile_name: Mutations',
                               'data_filename: data_mutations_extended.maf',
                               'swissprot_identifier: name'
                               ]
    data_mutations_extended = ['##fileformat=VCFv4.3',
                               '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
                               ]
    output_dir: Dict = dict()
    for snv in iter_snv(meta, filters=[]):
        data_mutations_extended.append('\t'.join([f"chr{snv.chrom}", str(snv.pos), snv.id,
                                                  snv.ref, snv.alt, snv.qual, snv.filter, snv.info]))

    if len(data_mutations_extended) > 2:
        output_dir = {'meta_mutations_extended.txt': meta_mutations_extended,
                      'data_snv.vcf': data_mutations_extended
                      }

    return output_dir


def get_cna_discrete(meta: CatsMeta) -> Dict:
    """Create values for various cna files and return file names as keys"""
    meta_cna_discrete = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                         'genetic_alteration_type: COPY_NUMBER_ALTERATION',
                         'datatype: DISCRETE',
                         'stable_id: gistic',
                         'show_profile_in_analysis_tab: true',
                         'profile_description: Test cna data. Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.',
                         'profile_name: Putative copy-number alterations from GISTIC',
                         'data_filename: data_cna_discrete.txt'
                         ]
    data_cna_discrete = [f'Hugo_Symbol\tEntrez_Gene_Id\t{meta.json_cats["testInfo"]["testId"]}']
    output_dir: Dict = dict()
    for cna in iter_cna(meta, filters=[]):
        data_cna_discrete.append('\t'.join([cna.gene, '', cna.cna_type]))

    if len(data_cna_discrete) > 1:
        output_dir = {'meta_cna_discrete.txt': meta_cna_discrete,
                      'data_cna_discrete.txt': data_cna_discrete
                      }

    return output_dir


def get_sv(meta: CatsMeta) -> Dict:
    """Create values for various sv files and return file names as keys"""
    meta_sv = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
               'genetic_alteration_type: STRUCTURAL_VARIANT',
               'datatype: SV',
               'data_filename: data_sv.txt',
               'stable_id: structural_variants',
               'profile_name: sv data',
               'profile_description: SVs',
               'show_profile_in_analysis_tab: true'
               ]
    data_sv = ['Sample_ID\tSite1_Entrez_Gene_Id\tSite1_Hugo_Symbol\tSite1_Ensembl_Transcript_Id\tSite1_Exon\tSite1_Chromosome\tSite1_Position\tSite1_Description\tSite2_Entrez_Gene_Id\tSite2_Hugo_Symbol\tSite2_Ensembl_Transcript_Id\tSite2_Exon\tSite2_Chromosome\tSite2_Position\tSite2_Description\tSite2_Effect_On_Frame\tNCBI_Build\tDNA_Support\tRNA_Support\tNormal_Read_Count\tTumor_Read_Count\tNormal_Variant_Count\tTumor_Variant_Count\tNormal_Paired_End_Read_Count\tTumor_Paired_End_Read_Count\tNormal_Split_Read_Count\tTumor_Split_Read_Count\tAnnotation\tBreakpoint_Type\tCenter\tConnection_Type\tEvent_Info\tClass\tLength\tComments\tExternal_Annotation\tcbp_driver\tcbp_driver_annotation\tcbp_driver_tiers\tcbp_driver_tiers_annotation']
    output_dir: Dict = dict()
    for sv in iter_sv(meta, filters=[]):
        grc_release = 'GRCh38' if 'GRCh38' in meta.json_cats["metaData"]["referenceGenome"]["grcRelease"] else 'GRCh37'
        data_sv.append('\t'.join([meta.json_cats["testInfo"]["testId"], 'NA',
                                  sv.breakends[0].gene, sv.breakends[0].transcript_id,
                                  '', str(sv.breakends[0].chrom), str(sv.breakends[0].start), '',
                                  'NA', sv.breakends[1].gene, sv.breakends[1].transcript_id, '',
                                  str(sv.breakends[1].chrom), str(sv.breakends[1].start),
                                  '', 'NA', grc_release,
                                  sv.dna, sv.rna, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                                  f'{sv.breakends[0].gene}-{sv.breakends[1].gene}',
                                  'NA', 'NA', 'NA', 'NA', sv.rearrangement_type, 'NA',
                                  sv.rearrangement_name[0], 'NA', 'NA', 'NA', 'NA', 'NA']))

    if len(data_sv) > 1:
        output_dir = {'meta_sv.txt': meta_sv,
                      'data_sv.txt': data_sv
                      }

    return output_dir


def get_fusion(meta: CatsMeta) -> Dict:
    """Create values for various fusion files and return file names as keys"""
    meta_fusion = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                   'genetic_alteration_type: FUSION',
                   'datatype: FUSION',
                   'stable_id: fusion',
                   'profile_description: Fusions.',
                   'show_profile_in_analysis_tab: true',
                   'profile_name: Fusions',
                   'data_filename: data_fusion.txt'
                   ]
    data_fusion = [
        'Hugo_Symbol\tEntrez_Gene_Id\tCenter\tTumor_Sample_Barcode\tFusion\tDNA_support\t'
        'RNA_support\tMethod\tFrame\tFusion_Status'
    ]
    output_dir: Dict = dict()
    for fusion in iter_sv(meta, filters=[filter_fusion]):
        data_fusion.append('\t'.join([meta.json_cats["testInfo"]["testId"], '', 'NA',
                                      meta.json_cats["testInfo"]["testId"],
                                      f'{fusion.breakends[0].gene}-{fusion.breakends[0].gene} fusion',
                                      fusion.dna, fusion.rna, 'unknown', 'unknown', '']))

    if len(data_fusion) > 1:
        output_dir = {'meta_fusion.txt': meta_fusion,
                      'data_fusion.txt': data_fusion
                      }

    return output_dir


def get_expression(meta: CatsMeta) -> Dict:
    """Create values for various expression files and return file names as keys"""
    meta_expression = [f'cancer_study_identifier: {meta.json_cats["testInfo"]["testId"]}',
                       'genetic_alteration_type: MRNA_EXPRESSION',
                       'datatype: CONTINUOUS',
                       'stable_id: rna_seq_mrna',
                       'show_profile_in_analysis_tab: true',
                       'profile_name: mRNA expression ',
                       'profile_description: Expression levels ',
                       'data_filename: data_expression.txt'
                       ]
    data_expression = [f'Hugo_Symbol\tEntrez_Gene_Id\t{meta.json_cats["testInfo"]["testId"]}']
    output_dir: Dict = dict()
    for expression in iter_expressions(meta, filters=[]):
        value = expression.expression_level_metrics[0].get('value', '')
        data_expression.append('\t'.join([expression.gene, '', '', str(value)]))

    if len(data_expression) > 1:
        output_dir = {'meta_expression.txt': meta_expression,
                      'data_expression.txt': data_expression
                      }

    return output_dir


def write_dataset(output_dir, output_data: Dict[str, list]):
    """Output data for cBioPortal input"""
    os.makedirs(os.path.join(output_dir, 'case_lists'), exist_ok=True)
    for key, value in output_data.items():
        with open(os.path.join(output_dir, key), "w", encoding="utf-8") as f:
            f.write('\n'.join(value))
