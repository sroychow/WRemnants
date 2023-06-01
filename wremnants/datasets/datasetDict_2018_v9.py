import copy

BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
xsec_Zmm = 2001.9
xsec_Wpmunu = 11765.9
xsec_Wmmunu = 8703.87
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)

dataDictV9_2018 = {
    'Wplusmunu' : { 
                      'filepaths' : 
                      ["{BASE_PATH}/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_{NANO_PROD_TAG}/*/*/*.root"],
                      'xsec' : xsec_Wpmunu,
                      'group': "Wmunu",
    },
    'Wminusmunu' : { 
                      'filepaths' : 
                      ["{BASE_PATH}/WminusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_{NANO_PROD_TAG}/*/*/*.root"],
                      'xsec' : xsec_Wmmunu,
                      'group': "Wmunu",
    },
    'Zmumu' : { 
                   'filepaths' :
                    ["{BASE_PATH}/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_{NANO_PROD_TAG}/*/*/*.root"],
                   'xsec' : xsec_Zmm,
                   'group': "Zmumu",
    },
}

'''
    'data' : { 
                      'filepaths' : ["{BASE_PATH}/SingleMuon/NanoV9Run2016FData_{NANO_PROD_TAG}/*/*/*.root",
                                     "{BASE_PATH}/SingleMuon/NanoV9Run2016GData_{NANO_PROD_TAG}/*/*/*.root",
                                     "{BASE_PATH}/SingleMuon/NanoV9Run2016HData_{NANO_PROD_TAG}/*/*/*.root",
                                     ],
                      'group': "Data",

    },
    'Ztautau' : { 
                   'filepaths' : 
                   ["{BASE_PATH}/DYJetsToTauTau_M-50_AtLeastOneEorMuDecay_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC_{NANO_PROD_TAG}/*/*/*.root"],
                   # At least one tau->e or mu decay, so everything that's not all other decays
                   'xsec' : xsec_Zmm*Z_TAU_TO_LEP_RATIO,
                   'group': "Ztautau",
    },


    'Wplustaunu' : { 
                         'filepaths' : 
                         ["{BASE_PATH}/WplusJetsToTauNu_TauToMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC_{NANO_PROD_TAG}/*/*/*.root"],
                         'xsec' : BR_TAUToMU*xsec_Wpmunu,
                         'group': "Wtaunu",
    },
    
    'Wminustaunu' : { 
                         'filepaths' : 
                         ["{BASE_PATH}/WminusJetsToTauNu_TauToMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC_{NANO_PROD_TAG}/*/*/*.root"],
                         'xsec' : BR_TAUToMU*xsec_Wmmunu,
                         'group': "Wtaunu",
    },
    'TTLeptonic' : { 
                          'filepaths' : 
                          ["{BASE_PATH}/BKGV9/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*.root"],
                          'xsec' : 88.29,
                          'group' : "Top",
    },

    'TTSemileptonic' : { 
                         'filepaths' : 
                         ["{BASE_PATH}/BKGV9/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*.root"],
                         'xsec' : 366.34,
                         'group' : "Top",
    },
    'SingleTschanLepDecays' : { 
                                          'filepaths' : 
                                          ["{BASE_PATH}/BKGV9/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/*.root"],
                                          'xsec' : 3.609,
                                          'group' : "Top",
    },
    'SingleTtWAntitop' : { 
                                     'filepaths' : 
                                     ["{BASE_PATH}/BKGV9/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                     'xsec' : 19.55, # 35.85 * (1.0-((1-0.1086*3)*(1-0.1086*3))) = 19.5 pb
                                     'group' : "Top",
    },
    'SingleTtWTop' : { 
                                     'filepaths' : 
                                     ["{BASE_PATH}/BKGV9/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                     'xsec' : 19.55,
                                     'group' : "Top",
    },
    'SingleTtchanAntitop' : { 
                                        'filepaths' : 
                                        ["{BASE_PATH}/BKGV9/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                        'xsec' : 80.0,
                                        'group' : "Top",
    },
    'SingleTtchanTop' : { 
                                    'filepaths' : 
                                    ["{BASE_PATH}/BKGV9/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*.root"],
                                    'xsec' : 134.2,
                                    'group' : "Top",
    },   
    # inclusive samples, keep for reference
    # 'WW' : { 
    #                 'filepaths' : 
    #                 ["{BASE_PATH}/BKGV9/WW_TuneCP5_13TeV-pythia8/*.root"],
    #                 'xsec' : 118.7,
    #                 'group' : "Diboson",
    # },
    # 'WZ' : { 
    #                 'filepaths' : 
    #                 ["{BASE_PATH}/BKGV9/WZ_TuneCP5_13TeV-pythia8/*.root"],
    #                 'xsec' : 47.026760,  # to check, taken from WZTo1L1Nu2Q dividing by BR: 10.71/(3*0.1086)/(1-3*0.033658-0.2)
    #                 'group' : "Diboson",
    # },
    ##
    'WWTo2L2Nu' : { 
        'filepaths' : 
        ["{BASE_PATH}/BKGV9/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*/*.root"],
        'xsec' : 12.6, # 118.7*0.1086*0.1086*9
        'group' : "Diboson",
    },
    'WWTo1L1Nu2Q' : { 
        'filepaths' : 
        ["{BASE_PATH}/BKGV9/WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/*/*.root"],
        'xsec' : 52.146, # 118.7*[(3*0.1086)*(1-3*0.1086)]*2 (2 is because one W or the other can go to Q)
        'group' : "Diboson",
    },
    'WZTo3LNu' : { 
        'filepaths' : 
        ["{BASE_PATH}/BKGV9/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/*/*.root"],
        'xsec' : 4.91, # 4.42965*1.109, 1.109 is the NLO to NNLO kfactor, for this one would need to make sure about the NLO XS, depends a lot on the dilepton mass cut
        'group' : "Diboson",
    },
    'WZTo2Q2L' : { 
        'filepaths' : 
        ["{BASE_PATH}/BKGV9/WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/*/*.root"],
        'xsec' : 5.4341, # 4.9*1.109
        'group' : "Diboson",
    },
    'WZTo1L1Nu2Q' : { 
        'filepaths' : 
        ["{BASE_PATH}/BKGV9/WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/*/*.root"],
        'xsec' : 11.781, # 10.71*1.10
        'group' : "Diboson",
    },
    'ZZTo2L2Nu' : { 
        'filepaths' : 
        ["{BASE_PATH}/BKGV9/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/*.root"],
        'xsec' : 0.60,
        'group' : "Diboson",
    },
    'ZZTo2Q2L' : { 
        'filepaths' : 
        ["{BASE_PATH}/BKGV9/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/*/*.root"],
        'xsec' : 5.1,
        'group' : "Diboson",
    },
    'QCDmuEnrichPt15' : { 
        'filepaths' : 
        ["{BASE_PATH}/BKGV9/QCD_Pt-20_MuEnrichedPt15_TuneCP5_13TeV-pythia8/*/*.root"],
        'xsec' : 238800,
        'group' : "QCD",
    }
}
'''
