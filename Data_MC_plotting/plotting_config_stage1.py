# Default config for plotting from the stage 1 ntuples, grouping uds = light jets

import os
import sys
import Zqq_plots 
import Zqq_processes 
import Zqq_selections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import run_list

class PlottingConfig:

    # Meta-data about the data plotted -> some of it used to define filenames as well
    year = "1994"
    sel_tag = "Selected events"
    lumi = run_list.luminosity_pb(year) # in pb-1, of the runs kept by stage1 (data/lumi/run_list_<year>.csv)
    ecm = 91. # in GeV
    norm_file = "normalisation.json" # json file with the cross-section information, following FCCAnalyses format
    mc_type = "zqq"
    ana_stage = "stage1"
    ntuple_tag ="v13_test_dEdx_dataMC"
    tree_name="events"

    #Run options
    n_threads = 32 #Put to "None" or other empty var to disable multithreading


    # Dictionary of all variables to plot: info about binning, range, labels etc - edit or add in Zqq_plots.py
    plots_dict = Zqq_plots.Zqq_data_MC_vars_dEdx

    # selection to apply
    selection_name="dedx_data_mc_quality_hi_p_cut"
    selection=Zqq_selections.selection_dictionary[selection_name]

    # need to mask branches to be able to apply selection on jet constituents
    branches_to_mask=Zqq_plots.masked_branches_dict["Zqq_data_MC_vars_dEdx"] #TODO unify with plots dict get both of them!

    # Dictionary of processes: info about samples, colours, labels - edit or add in Zqq_processes.py
    data = Zqq_processes.zqq_data["data"]
    mc_processes = Zqq_processes.MC_group_light_jets


    # Plotting style settings:
    do_log_y = True
    add_overflow = False
    ratio_range = (0., 2.)
    weighted = False # whether to apply generator level weights or not, currently not needed keep on false

    # Output file settings:
    out_format = ".png"
    store_root_file= False #if true a ROOT file of the plotted histograms is written that can be used by e.g. combine for fitting

    # Input and output paths 
    inputs_path_data = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/{year}/{ana_stage}/{ntuple_tag}"
    inputs_path_mc = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/{year}/{mc_type}/{ana_stage}/{ntuple_tag}"
    outputs_path = f"/eos/user/b/bistapf/ALEPH_plots/{ntuple_tag}/{selection_name}"