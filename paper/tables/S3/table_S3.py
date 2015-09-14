##Code for creating table S3

##load necessary modules
import os


BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))



##create a dictionary for each msigdb gene set
msigdb_dict={}
f=open(os.path.join(BASE_DIR,'figures','figure_4','msigdb.v5.0.entrez.gmt'))
for i in f:
    msigdb_dict[i.split('\t')[0]]=i.strip().split('\t')[2:]


##additional gene sets from PerouLab BMC Med Genomics. 2011 | PMID: 21214954
msigdb_dict['HS_Red7']=['10085', '1009', '10154', '10184', '10234', '10398', '10409', '10417', '10418', '10457', '10468', '10484', '10536', '10544', '10609', '10630', '10631', '10699', '10970', '10981', '11015', '11027', '11031', '11075', '11082', '11096', '11098', '11117', '11167', '11183', '112399', '11259', '113', '11346', '114899', '1176', '1272', '1277', '1278', '1281', '1282', '1289', '1290', '1292', '1293', '1300', '1301', '1306', '1307', '1311', '1318', '1404', '1462', '1490', '150', '1508', '1513', '1519', '1601', '1604', '1634', '165', '1727', '1805', '1809', '1833', '1842', '1893', '1909', '2009', '2012', '2014', '2018', '2028', '2054', '206358', '2121', '2132', '2162', '2191', '2192', '2199', '2200', '2201', '2212', '2274', '22795', '22846', '22891', '23022', '23090', '231', '2316', '23176', '23213', '2331', '23414', '23452', '23516', '23597', '23601', '23627', '23643', '23670', '23705', '23768', '24141', '2512', '2535', '25878', '2589', '25900', '25903', '25907', '25932', '25987', '26010', '26011', '26064', '2619', '2621', '2632', '26585', '2697', '27122', '27286', '27295', '27306', '283537', '284119', '285', '2882', '28962', '28984', '28986', '2922', '2983', '29931', '29940', '29951', '29967', '29995', '30008', '307', '308', '30851', '3339', '3357', '347902', '3491', '355', '3611', '3624', '3671', '3678', '3685', '3693', '3730', '3751', '3912', '3956', '4015', '4016', '4017', '4035', '4052', '4053', '4060', '4092', '4124', '4148', '4237', '4267', '4312', '4313', '4314', '4320', '4322', '4323', '4330', '4616', '4627', '4642', '4675', '4681', '4692', '4811', '4815', '4837', '4907', '4920', '4921', '4958', '4969', '4973', '5033', '50507', '50509', '5054', '50863', '5099', '5118', '51226', '51280', '51313', '51339', '51363', '51375', '51454', '5157', '5159', '51715', '5176', '5218', '5311', '5327', '5328', '5329', '5350', '5352', '5358', '5376', '54209', '54210', '54453', '54587', '5475', '54796', '5480', '54829', '54985', '55033', '55068', '5507', '55076', '55083', '55084', '5530', '55357', '55454', '55568', '55714', '55740', '558', '55816', '55830', '55841', '56034', '5627', '5654', '56925', '56935', '56937', '56944', '57088', '57125', '57333', '57493', '5789', '579', '5793', '58189', '5919', '5954', '5961', '6004', '60681', '6091', '6275', '6277', '6281', '6310', '6320', '633', '6356', '6382', '6383', '6387', '64175', '64236', '6424', '64359', '6444', '6447', '6474', '64754', '64759', '649', '6525', '6526', '6591', '6624', '6678', '6695', '6696', '6764', '683', '6876', '6925', '7041', '7043', '7045', '7057', '7058', '7059', '7060', '7070', '7077', '7078', '7107', '7130', '7132', '716', '7171', '7280', '7292', '7345', '7358', '7414', '7423', '7431', '7472', '7474', '753', '7533', '754', '7837', '79070', '79586', '79614', '79652', '79776', '79783', '79899', '79953', '800', '80176', '80206', '8038', '8076', '80781', '8082', '813', '8321', '83468', '83700', '8406', '84168', '8434', '84561', '8460', '84617', '8483', '84909', '8532', '857', '8630', '871', '8754', '8829', '8839', '8840', '8877', '8933', '8974', '9060', '9189', '9242', '9260', '9270', '9315', '9358', '9448', '9509', '9532', '9540', '9590', '9638', '9645', '9732', '9781', '9843', '9886', '9891', '9902', '9945', '9955']
msigdb_dict['Mclaudin_cluster']=['83481', '5792', '27134', '3815', '7571', '340371', '10809', '2068', '6565', '114569', '163732', '23254', '113452', '85462', '1159', '780', '3898', '1525', '23286', '79977', '6663', '11187', '3664', '7022', '152189', '55286', '219738', '85415', '2001', '54682', '29841', '349633', '3775', '146439', '149428', '91862', '5652', '153562', '4950', '3875', '286077', '254427', '10653', '27237', '2886', '200634', '51599', '1365', '92359', '126695', '10053', '1366', '6768', '4072', '64787', '57662', '79755', '286676', '132014', '6712', '54566', '376267', '22996', '54836', '54894', '26298', '3934', '149466', '6692', '999', '90060', '79784', '1397', '93664', '10207', '55620', '2591', '10279', '26751', '5364', '1299', '9414', '100000000', '245973', '8495', '57475', '55930', '10040', '7274', '200010', '1297', '55107', '230', '5169', '161291', '23541']
msigdb_dict['HS_Green17']=['10067', '10092', '10158', '10228', '10262', '103', '10456', '10472', '10623', '1063', '10638', '10654', '10712', '10753', '10765', '10899', '11243', '11266', '11279', '116832', '1196', '123', '1314', '1382', '142', '159', '1660', '1942', '1945', '2029', '2058', '2155', '2224', '223', '2271', '22796', '22874', '22926', '23219', '23248', '23385', '23528', '23623', '2444', '25782', '25874', '25896', '259266', '25936', '26154', '26750', '27042', '27089', '27097', '27101', '27173', '27246', '27252', '27285', '2745', '28956', '29097', '29104', '29108', '29956', '3068', '3140', '3664', '3707', '375035', '3775', '3790', '3964', '4170', '4486', '4580', '4582', '4720', '4751', '4931', '50848', '51018', '51022', '51027', '51029', '51094', '51097', '51107', '51182', '51205', '5121', '51377', '51430', '51506', '51603', '5194', '5202', '54344', '54499', '54530', '54583', '54953', '54996', '55034', '55105', '55127', '55157', '55249', '55526', '55585', '55699', '55732', '55740', '55746', '55758', '55974', '5664', '5692', '56950', '57111', '57147', '57402', '57645', '5796', '5824', '5877', '5993', '6016', '6045', '6051', '607', '6282', '6286', '63931', '64216', '64222', '64746', '64757', '65123', '6635', '6675', '6726', '6746', '6905', '6944', '7100', '7159', '7203', '7257', '7286', '7678', '7818', '79005', '79017', '79577', '79590', '80232', '81609', '81875', '824', '84288', '8443', '8560', '8564', '8703', '874', '8751', '8799', '8804', '8934', '8991', '9015', '9019', '90806', '9129', '92342', '9453', '9557', '9588', '9673', '9816', '9826', '9857', '9869', '9887', '9898', '9917', '9926']

###generate gene programs based on Hoadley et al. 2014
Proliferation_DNA_Repair={}
for i in msigdb_dict['PUJANA_CHEK2_PCC_NETWORK']:
    Proliferation_DNA_Repair[i]=''
for i in msigdb_dict['REACTOME_CELL_CYCLE_MITOTIC']:
    Proliferation_DNA_Repair[i]=''
for i in msigdb_dict['KEGG_CELL_CYCLE']:
    Proliferation_DNA_Repair[i]=''
for i in msigdb_dict['REACTOME_DNA_REPAIR']:
    Proliferation_DNA_Repair[i]=''

Immunecell={}
for i in msigdb_dict['REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM']:
    Immunecell[i]=''
for i in msigdb_dict['KEGG_HEMATOPOIETIC_CELL_LINEAGE']:
    Immunecell[i]=''
for i in msigdb_dict['REACTOME_TCR_SIGNALING']:
    Immunecell[i]=''
for i in msigdb_dict['KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY']:
    Immunecell[i]=''

Tumor_suppressing_miRNA_targets={}
for i in msigdb_dict['DACOSTA_UV_RESPONSE_VIA_ERCC3_DN']:
    Tumor_suppressing_miRNA_targets[i]=''
for i in msigdb_dict['GTTTGTT,MIR-495']:
    Tumor_suppressing_miRNA_targets[i]=''
for i in msigdb_dict['TGCTTTG,MIR-330']:
    Tumor_suppressing_miRNA_targets[i]=''
for i in msigdb_dict['ACATTCC,MIR-1,MIR-206']:
    Tumor_suppressing_miRNA_targets[i]=''

MES_ECM={}
for i in msigdb_dict['HS_Red7']:
    MES_ECM[i]=''
for i in msigdb_dict['KEGG_ECM_RECEPTOR_INTERACTION']:
    MES_ECM[i]=''
for i in msigdb_dict['ORGAN_DEVELOPMENT']:
    MES_ECM[i]=''
for i in msigdb_dict['KEGG_FOCAL_ADHESION']:
    MES_ECM[i]=''

MYC_targets_TERT={}
for i in msigdb_dict['DANG_BOUND_BY_MYC']:
    MYC_targets_TERT[i]=''
for i in msigdb_dict['DAIRKEE_TERT_TARGETS_UP']:
    MYC_targets_TERT[i]=''

Squamous_differentiation_development={}
for i in msigdb_dict['RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_POORLY_DN']:
    Squamous_differentiation_development[i]=''

Estrogen_signaling={}
for i in msigdb_dict['SMID_BREAST_CANCER_BASAL_DN']:
    Estrogen_signaling[i]=''
for i in msigdb_dict['VANTVEER_BREAST_CANCER_ESR1_UP']:
    Estrogen_signaling[i]=''

FOXO_stemness={}
for i in msigdb_dict['TTGTTT_V$FOXO4_01']:
    FOXO_stemness[i]=''
for i in msigdb_dict['YTATTTTNR_V$MEF2_02']:
    FOXO_stemness[i]=''

Cellcell_adhesion={}
for i in msigdb_dict['Mclaudin_cluster']:
    Cellcell_adhesion[i]=''

Fatty_acid_oxidation={}
for i in msigdb_dict['CARBOXYLIC_ACID_METABOLIC_PROCESS']:
    Fatty_acid_oxidation[i]=''
for i in msigdb_dict['OXIDOREDUCTASE_ACTIVITY']:
    Fatty_acid_oxidation[i]=''

Interferon={}
for i in msigdb_dict['ZHANG_INTERFERON_RESPONSE']:
    Interferon[i]=''
for i in msigdb_dict['BROWNE_INTERFERON_RESPONSIVE_GENES']:
    Interferon[i]=''

Hypoxia_glycolysis={}
for i in msigdb_dict['SEMENZA_HIF1_TARGETS']:
    Hypoxia_glycolysis[i]=''
for i in msigdb_dict['REACTOME_GLYCOLYSIS']:
    Hypoxia_glycolysis[i]=''

Neural_signaling={}
for i in msigdb_dict['MODULE_100']:
    Neural_signaling[i]=''

Plasma_membrane_cellcell_signaling={}
for i in msigdb_dict['MORF_CNTN1']:
    Plasma_membrane_cellcell_signaling[i]=''
for i in msigdb_dict['MORF_CASP2']:
    Plasma_membrane_cellcell_signaling[i]=''
for i in msigdb_dict['MORF_PML']:
    Plasma_membrane_cellcell_signaling[i]=''
for i in msigdb_dict['MORF_DDX11']:
    Plasma_membrane_cellcell_signaling[i]=''

EGF_signaling={}
for i in msigdb_dict['NAGASHIMA_EGF_SIGNALING_UP']:
    EGF_signaling[i]=''

MAPKs={}
for i in msigdb_dict['INTRACELLULAR_SIGNALING_CASCADE']:
    MAPKs[i]=''

Basal_signaling={}
for i in msigdb_dict['SMID_BREAST_CANCER_BASAL_UP']:
    Basal_signaling[i]=''
for i in msigdb_dict['DOANE_BREAST_CANCER_ESR1_DN']:
    Basal_signaling[i]=''

Vesicle={}
for i in msigdb_dict['MEMBRANE_COAT']:
    Vesicle[i]=''

amplicon_1Q={}
for i in msigdb_dict['HS_Green17']:
    amplicon_1Q[i]=''

TAL1={}
for i in msigdb_dict['GNF2_TAL1']:
    TAL1[i]=''
for i in msigdb_dict['GNF2_ANK1']:
    TAL1[i]=''
for i in msigdb_dict['GNF2_SPTB']:
    TAL1[i]=''
for i in msigdb_dict['GNF2_SPTA1']:
    TAL1[i]=''
for i in msigdb_dict['GNF2_MAP2K3']:
    TAL1[i]=''

Antiapoptosis={}
for i in msigdb_dict['MORF_MT4']:
    Antiapoptosis[i]=''

amplicon_16Q22={}
for i in msigdb_dict['chr16q22']:
    amplicon_16Q22[i]=''


pathways=[Proliferation_DNA_Repair,Immunecell,Tumor_suppressing_miRNA_targets,MES_ECM,MYC_targets_TERT,Squamous_differentiation_development,\
          Estrogen_signaling,FOXO_stemness,Cellcell_adhesion,Fatty_acid_oxidation,Interferon,Hypoxia_glycolysis,Neural_signaling,\
          Plasma_membrane_cellcell_signaling,EGF_signaling,MAPKs,Basal_signaling,Vesicle,amplicon_1Q,TAL1,Antiapoptosis,amplicon_16Q22]



pathway_names=['Proliferation_DNA_Repair', 'Immunecell', 'Tumor_suppressing_miRNA_targets', 'MES_ECM', 'MYC_targets_TERT', 'Squamous_differentiation_development',\
               'Estrogen_signaling', 'FOXO_stemness', 'Cellcell_adhesion', 'Fatty_acid_oxidation', 'Interferon', 'Hypoxia_glycolysis', 'Neural_signaling',\
               'Plasma_membrane_cellcell_signaling', 'EGF_signaling', 'MAPKs', 'Basal_signaling', 'Vesicle', 'amplicon_1Q', 'TAL1', 'Antiapoptosis', 'amplicon_16Q22']


##print how many genes are in each pathway

for i,j in zip(pathway_names,pathways):
    print i, '\t', len(j)



