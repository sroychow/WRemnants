from utilities import boostHistHelpers as hh, common, logging, differential
from utilities.io_tools import output_tools

parser,initargs = common.common_parser(True)

import ROOT
import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_validation, muon_calibration, muon_selections, unfolding_tools
from wremnants.histmaker_tools import scale_to_data, aggregate_groups
from wremnants.datasets.dataset_tools import getDatasets
from wremnants.datasets.datagroups import Datagroups
import hist
import lz4.frame
import math
import time
import os

parser.add_argument("--csVarsHist", action='store_true', help="Add CS variables to dilepton hist")
parser.add_argument("--axes", type=str, nargs="*", default=["mll", "etaPlus", "ptPlus", "etaMinus", "ptMinus"], help="")
parser.add_argument("--finePtBinning", action='store_true', help="Use fine binning for ptll")
parser.add_argument("--useDileptonTriggerSelection", action='store_true', help="Use dilepton trigger selection (default uses the Wlike one, with one triggering muon and odd/even event selection to define its charge, staying agnostic to the other)")
parser.add_argument("--noAuxiliaryHistograms", action="store_true", help="Remove auxiliary histograms to save memory (removed by default with --unfolding or --theoryAgnostic)")

parser = common.set_parser_default(parser, "genVars", ["ptVGen", "absYVGen"])
parser = common.set_parser_default(parser, "pt", [34,26.,60.])
parser = common.set_parser_default(parser, "eta", [48,-2.4,2.4])
parser = common.set_parser_default(parser, "aggregateGroups", ["Diboson", "Top", "Wtaunu", "Wmunu"])
parser = common.set_parser_default(parser, "theoryCorr", ["scetlib_dyturbo", "virtual_ew", "horaceqedew_FSR", "horacelophotosmecoffew_FSR"])

args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

thisAnalysis = ROOT.wrem.AnalysisType.Dilepton if args.useDileptonTriggerSelection else ROOT.wrem.AnalysisType.Wlike
era = args.era

############NOTE
####For mow this is hardcoded to run on Zmumu only
filterProcess=['Zmumu']
datasets = getDatasets(maxFiles=args.maxFiles,
                       filt=filterProcess,
                       excl=args.excludeProcs, 
                       nanoVersion="v9",
                       base_path=args.dataPath,
                       extended = "msht20an3lo" not in args.pdfs,
                       era = era)

# dilepton invariant mass cuts
mass_min = 60
mass_max = 120

ewMassBins = theory_tools.make_ew_binning(mass = 91.1535, width = 2.4932, initialStep=0.010)
dilepton_ptV_binning = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 20, 23, 27, 32, 40, 54, 100] if not args.finePtBinning else range(60)
# available axes for dilepton validation plots
all_axes = {
    "mll": hist.axis.Regular(60, 60., 120., name = "mll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    # "mll": hist.axis.Variable([60,70,75,78,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,100,102,105,110,120], name = "mll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "yll": hist.axis.Regular(20, -2.5, 2.5, name = "yll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "absYll": hist.axis.Regular(10, 0., 2.5, name = "absYll", underflow=False, overflow=not args.excludeFlow),
    "ptll": hist.axis.Variable(dilepton_ptV_binning, name = "ptll", underflow=False, overflow=not args.excludeFlow),
    "etaPlus": hist.axis.Variable([-2.4,-1.2,-0.3,0.3,1.2,2.4], name = "etaPlus"),
    "etaMinus": hist.axis.Variable([-2.4,-1.2,-0.3,0.3,1.2,2.4], name = "etaMinus"),
    "absEtaPlus": hist.axis.Regular(8, 0, 2.4, name = "absEtaPlus"),
    "absEtaMinus": hist.axis.Regular(8, 0, 2.4, name = "absEtaMinus"),
    "etaSum": hist.axis.Regular(12, -4.8, 4.8, name = "etaSum"),
    "etaDiff": hist.axis.Variable([-4.8, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 4.8], name = "etaDiff"),
    "ptPlus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "ptPlus"),
    "ptMinus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "ptMinus"),
    "cosThetaStarll": hist.axis.Regular(20, -1., 1., name = "cosThetaStarll", underflow=False, overflow=False),
    "phiStarll": hist.axis.Regular(20, -math.pi, math.pi, circular = True, name = "phiStarll"),
    #"charge": hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge") # categorical axes in python bindings always have an overflow bin, so use a regular
    "massVgen": hist.axis.Variable(ewMassBins, name = "massVgen", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewMll": hist.axis.Variable(ewMassBins, name = "ewMll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewMlly": hist.axis.Variable(ewMassBins, name = "ewMlly", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewLogDeltaM": hist.axis.Regular(100, -10, 4, name = "ewLogDeltaM", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
}

auxiliary_gen_axes = ["massVgen", # preFSR variables
    "ewMll", "ewMlly", "ewLogDeltaM" # ew variables
    ]

for a in args.axes:
    if a not in all_axes.keys():
        logger.error(f" {a} is not a known axes! Supported axes choices are {list(axes.keys())}")

nominal_cols = args.axes

if args.csVarsHist:
    nominal_cols += ["cosThetaStarll", "phiStarll"]

nominal_axes = [all_axes[a] for a in nominal_cols] 

gen_axes = {
    "ptVGen": hist.axis.Variable(dilepton_ptV_binning, name = "ptVGen", underflow=False, overflow=args.poiAsNoi),
    "absYVGen": hist.axis.Regular(10, 0, 2.5, name = "absYVGen", underflow=False, overflow=args.poiAsNoi),  
}


# define helpers
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

pileup_helper = wremnants.make_pileup_helper(era = era)
vertex_helper = wremnants.make_vertex_helper(era = era)

calib_filepaths = common.calib_filepaths
closure_filepaths = common.closure_filepaths
diff_weights_helper = ROOT.wrem.SplinesDifferentialWeightsHelper(calib_filepaths['tflite_file']) if (args.muonScaleVariation == 'smearingWeightsSplines' or args.validationHists) else None
mc_jpsi_crctn_helper, data_jpsi_crctn_helper, mc_jpsi_crctn_unc_helper, data_jpsi_crctn_unc_helper = muon_calibration.make_jpsi_crctn_helpers(args, calib_filepaths, make_uncertainty_helper=True)
z_non_closure_parametrized_helper, z_non_closure_binned_helper = muon_calibration.make_Z_non_closure_helpers(args, calib_filepaths, closure_filepaths)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = muon_calibration.make_muon_calibration_helpers(args,era=era)

smearing_helper, smearing_uncertainty_helper = (None, None) if args.noSmearing else muon_calibration.make_muon_smearing_helpers()

bias_helper = muon_calibration.make_muon_bias_helpers(args) 

corr_helpers = theory_corrections.load_corr_helpers([d.name for d in datasets if d.name in common.vprocs], args.theoryCorr)

def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isWorZ = isW or isZ

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper
    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    axes = nominal_axes
    cols = nominal_cols
    #This is for HLT
    df = df.Filter(muon_selections.hlt_string(era))
    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)

    df = muon_calibration.define_corrected_muons(df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper)

    df = muon_selections.select_veto_muons(df, nMuons=2)
    df = muon_selections.select_good_muons(df, args.pt[1], args.pt[2], dataset.group, nMuons=2, use_isolation=True, isoDefinition=args.isolationDefinition, condition=">=")

    df = df.Define("rGaus4", "wrem::gausInSlot(rdfslot_, event)")
    #for debugging Random generation
    #df = df.Define("ff", "wrem::plotInSlot(rdfslot_, event, rGaus4)")
    #df = df.Filter('ff')
    
    df = df.Define("pairs", "wrem::muPair(rdfslot_,Muon_correctedPt, Muon_correctedCharge, Muon_correctedEta, Muon_correctedPhi, Muon_dxy, Muon_dz, GenPart_status, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, rGaus4)")
    df = df.Define("mll_reco","return get<2>(pairs);");
    df = df.Filter("mll_reco>1.0");


    df = df.Define("posTrackPt","return get<3>(pairs);")\
           .Define("negTrackPt","return get<4>(pairs);")\
           .Define("mll_diff_reco","return get<5>(pairs);")\
           .Define("jacobian_weight_mll_diff_reco","return get<5>(pairs)*std::copysign(1.0, genWeight);")\
           .Define("mll_smear","return get<6>(pairs);")\
           .Define("posPtSmear","return get<7>(pairs);")\
           .Define("negPtSmear","return get<8>(pairs);")\
           .Define("mll_diff_smear","return get<9>(pairs);")\
           .Define("jacobian_weight_mll_diff_smear", "return get<9>(pairs)*std::copysign(1.0, genWeight);")\
           .Define("mll_gen","return get<10>(pairs);")\
           .Define("posPtGen","return get<11>(pairs);")\
           .Define("negPtGen","return get<12>(pairs);")\
           .Define("jacobian_weight_mll_diff_squared_smear","return get<13>(pairs)*std::copysign(1.0, genWeight);")\
           .Define("smear_beta_weight","return get<14>(pairs)*std::copysign(1.0, genWeight);")\
           .Define("posPtSmearBetaVal","return get<15>(pairs);")\
           .Define("negPtSmearBetaVal","return get<16>(pairs);")\
           .Define("jacobian_weight_mll_diff_squared_reco","return get<17>(pairs)*std::copysign(1.0, genWeight);")\
           .Define("posTrackEta","return Muon_correctedEta[get<0>(pairs)];")\
           .Define("negTrackEta","return Muon_correctedEta[get<1>(pairs)];")\


    #if dataset.is_data:
    #    #
    #    results.append(df.HistoBoost("mll", [all_axes['mll']], ['mll_reco']))
    #else:
    df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
    df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])
    df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])
    
    if era == "2016PostVFP":
        weight_expr = "weight_pu*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
    else:
        weight_expr = "weight_pu*L1PreFiringWeight_Muon_Nom*L1PreFiringWeight_ECAL_Nom"
        
        
    if not args.noVertexWeight:
        weight_expr += "*weight_vtx"            
            
    logger.debug(f"Experimental weight defined: {weight_expr}")
    df = df.Define("exp_weight", weight_expr)
    df = theory_tools.define_theory_weights_and_corrs(df, dataset.name, corr_helpers, args)
    results.append(df.HistoBoost("weight", [hist.axis.Regular(100, -2, 2)], ["nominal_weight"], storage=hist.storage.Double()))
    results.append(df.HistoBoost("mll", [all_axes['mll']], ['mll_reco', 'nominal_weight']))

    
    return results, weightsum

logger.debug(f"Datasets are {[d.name for d in datasets]}")
resultdict = narf.build_and_run(datasets, build_graph)

if not args.noScaleToData:
    scale_to_data(resultdict)
    aggregate_groups(datasets, resultdict, args.aggregateGroups)

output_tools.write_analysis_output(resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args)
