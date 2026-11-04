### Config file with some information multiple analysis scripts need to access

# Marker groups used for broad spinal cord analysis
custom_SC_keys = [
    "dI1",
    "dI2",
    "dI3",
    "dI4",
    "dI5",
    "dI6",
    "V0",
    "V1",
    "V2a",
    "V2b",
    "V2c",
    "V3",
    "MN",
    "inhibitory SC INs",
    "neurons",
    "neural progenitors",
    "cholinergic neurons",
    "glycinergic neurons",
    "sensory neurons",
    "dILA/B",
    "FP",
    "RP",
    "differentiating neurons",
    "temporal TF",
    "excitatory neurons",
    "Astrocytes",
    "Pericytes / smooth muscle cells",
    "Endothelial cells",
    "Microglia / perivascular macrophages",
    "Oligodendrocytes",
    "Cadherins"
]


# Non-neuronal groups excluded from the neuron-focused analysis
non_neuron_keys = [
    "Astrocytes",
    "Pericytes / smooth muscle cells",
    "Endothelial cells",
    "Microglia / perivascular macrophages",
    "Oligodendrocytes",
]


# Marker groups used for neuron-focused analysis
SC_neuron_keys = [
    key
    for key in custom_SC_keys
    if key not in non_neuron_keys
]


cell_type_colors = {
    # dI1
    "dI1": "#ABDDA4",

    # dI2
    "dI2": "#17B684",

    # dI3
    "dI3": "#7BC87C",

    # dI4 - 9 true-blue shades
    "dI4-1": "#08306B",
    "dI4-2": "#08519C",
    "dI4-3": "#2171B5",
    "dI4-4": "#D0E1ED",
    "dI4-5": "#6BAED6",
    "dI4-6": "#9ECAE1",
    "dI4-7": "#BDD7E7",
    "dI4-8": "#4292C6",

    # dI5 - 6 cyan / teal-blue shades
    "dI5-1": "#005B5F",
    "dI5-2": "#007C83",
    "dI5-3": "#B9E3E4",
    "dI5-4": "#45B7BD",
    "dI5-5": "#83CDD0",
    "dI5-6": "#159CA4",

    # dI6
    "dI6": "#E6F598",

    # V0 - more distinct yellows
    "V0-1": "#D99A00",
    "V0-2": "#FEE08B",

    # V1
    "V1": "#FDCB6E",

    # V2a
    "V2a": "#FDAE61",

    # V2b
    "V2b": "#F46D43",

    # V3
    "V3": "#D53E4F",

    # Motor neurons
    "MN": "#8B1E3F",
}

### CAREFUL ABOUT TOUCHING THESE CAUSE IT MAY CHANGE THE NEURON CLUSTERING
select_genes = [
    'Slc17a6',
    'Slc6a5',
    'Slc32a1',
    'Gad2',
    'Pax2',
    'Lhx5',    
    'Slc18a3',
    'Slc10a4',  
    
    'Barhl1',
    'Barhl2',
    'Lhx2',
    'Lhx9',

    'Hmx2',
    'Hmx3',
    'Foxd3',
    'Nr4a2',

    'Isl1',
    'Otp',

    'Reln',
    'Zic1',
    'Lbx1',
    'Skor1',
    'Tfap2b',   
    'Lhx1',
    'Prrxl1',
    'Gbx1',
    'Rorb',
    'Kcnip2',
    'Pdyn',
    'Bhlhe22',
    'Mafb',
    'Nppa',
    'Pax8',
    'Cdh23',
    'Zic2',
    'Zic3',
    'Zic4',
    'Zic5',
    'Lmx1b',
    'Tlx1',
    'Phox2a',
    'Sox21',
    'Tlx3',
    'Grp',
    'Pth2r',
    'Nhlh1',
    
    
    'Irx1',
    'Irx2',
    'Dmrt3',
    'Wt1',

    'Evx1',
    'Evx2',
    'Pitx2',

    'En1',
    'Foxp2',

    'Vsx2',
    'Sox14',
    'Shox2',
    'Lhx3',
    'Tal1',

    'Gata2',
    'Gata3',
    'Mnx1',

    'Isl2',
    'Mecom',
    'Chat',
    'Olig2',
    'Foxp1',
    'Aldh1a2',
    'Arhgap36',
    'Shh',

    'Nkx2-2',
    'Sim1',
    'Olig3',
    
    

    'Pou4f1',
    'Zfhx3',
    'Zfhx4',
    'Nfia',
    'Nfib',
    'Neurod2',
    'Neurod6',
]

# Zics, TTFs, dI4/5 subtype markers, go in supp maybe in "subtype figure"
select_genes_Lora = [
    'Slc17a6',
    'Slc6a5',
    'Slc32a1',
    'Pax2',
    'Slc18a3',
    'Slc10a4',  
    
    'Barhl1',
    'Barhl2',

    'Hmx2',
    'Hmx3',
    'Foxd3',
    'Nr4a2',

    'Isl1',
    'Otp',

    'Lbx1',
    'Tfap2b',  
    'Gbx1',
    'Rorb',
    'Pdyn',
    'Mafb',
    'Sox21',
    'Pax8',
    'Cdh23',
    'Lmx1b',
    'Tlx1',
    'Phox2a',
    'Tlx3',
    'Grp',
    'Nhlh1',
    
    'Irx2',
    'Dmrt3',
    'Wt1',

    'Evx1',
    'Evx2',
    'Pitx2',

    'En1',

    'Vsx2',
    'Sox14',
    'Shox2',
    'Lhx3',
    'Tal1',

    'Gata2',
    'Gata3',
    'Mnx1',

    'Isl2',
    'Chat',

    'Nkx2-2',
    'Sim1',
    'Olig3'
]

MADM_GENES = [
    'Slc17a6',
    'Slc6a5',
    'Slc32a1',
    'Pax2',

    'Lbx1',
    'Tfap2b',  
    #'Rorb',
    #'Pdyn',
    'Mafb',
    #'Sox21',
    'Pax8',
    #'Cdh23',
    
    'Irx2',
    'Dmrt3',

    'En1',
    'Foxd3',
    'Otp',

    'Vsx2',
    'Sox14',
    #'Shox2',
    'Lhx3',
    'Tal1',
    #'Gata2',
    'Gata3',
    #'Mnx1',
]

CLONE1_GENES = [
    'Zfhx3',
    'Zfhx4',
    'Nfia',
    'Nfib',
    'Neurod2',
    'Neurod6',

    "Zic1",
    "Zic2",
    "Zic3",
    "Zic4",
    "Zic5",
    
    'Cdh4',
    'Cdh6',
    'Cdh9',
    'Cdh13',
    'Cdh20',
    'Cdh23',
    
    "Lbx1",
    "Tfap2b",
    "Igfbp5",
    "Reln",
    "Lhx5",
    "Nrp2",
    "Skor1",
    "Mafb",
    "Penk",
    "Neurog2",
    "Mafa",
    "Rorb",
    "Kcnip2",
    "Pax3",
    "Lhx1",
    "Pdyn",
    "Prrxl1",      
]

CLONE4_GENES = [
    'Dmrt3',
    'Pax8',
    'Irx1',
    'En1',
    'Foxd3',
    'Otp',
    'Foxp2',
    'Reln',
    'Adamtsl1',
    'Calb2',
    'Chodl',
    'Vsx2',
    'Sox14',
    'Shox2',
    'Mafa',
    'Lhx3',
    'Ebf3',
    'Prox1',
    'Gata3',
    'Tal1',
    'Mafb',
    'Nr4a2',

    'Zfhx3',
    'Zfhx4',
    'Nfia',
    'Nfib',
    'Neurod2',
    'Neurod6',
    
    'Cdh4',
    'Cdh6',
    'Cdh9',
    'Cdh13',
    'Cdh20',
    'Cdh23',
]

TTF_CDH_ZIC_GENES = [
    'Zfhx3',
    'Zfhx4',
    'Nfia',
    'Nfib',
    'Neurod2',
    'Neurod6',
    
    'Cdh4',
    'Cdh6',
    'Cdh9',
    'Cdh13',
    'Cdh20',
    'Cdh23',
    
    "Zic1",
    "Zic2",
    "Zic3",
    "Zic4",
    "Zic5",
]

V2a_GENES = [
    'Zfhx3',
    'Zfhx4',
    'Nfia',
    'Nfib',
    'Neurod2',
    'Neurod6',
    
    'Cdh4',
    'Cdh6',
    'Cdh9',
    'Cdh13',
    'Cdh20',
    'Cdh23',

    'Vsx2', 
    'Shox2',
    'Sox14', 
    'Lhx3', 
    'Reln', 
    'Chodl',
    'Prox1',
    'Ebf3']

dI6_GENES = [
    'Zfhx3',
    'Zfhx4',
    'Nfia',
    'Nfib',
    'Neurod2',
    'Neurod6',
    
    'Cdh4',
    'Cdh6',
    'Cdh9',
    'Cdh13',
    'Cdh20',
    'Cdh23',

    'Dmrt3',
    'Irx2',
    'Wt1',
    'Ebf3',
    'Nr4a2',
    'Gbx1',
    'Bhlhe22',
    'Penk']

V2b_GENES = [
    'Zfhx3',
    'Zfhx4',
    'Nfia',
    'Nfib',
    'Neurod2',
    'Neurod6',
    
    'Cdh4',
    'Cdh6',
    'Cdh9',
    'Cdh13',
    'Cdh20',
    'Cdh23',
 
    'Gata2',
    'Gata3',
    'Prox1',
    'Bhlhe22',
    'Mafa',
    'Mafb',
    'Reln',
    'Nr4a2']