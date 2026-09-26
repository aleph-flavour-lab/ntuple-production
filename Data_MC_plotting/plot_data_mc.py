# Script to plot data and MC from ALEPH

import ROOT
import os 
import sys
from collections import namedtuple
from array import array
import copy
import numpy
from argparse import ArgumentParser
import importlib
import json
import glob

# from plotting_config_stage1 import PlottingConfig

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

#Helper inline C++ for the RDF plotting

ROOT.gInterpreter.Declare("""
ROOT::RVec<double> filterByFlag(const ROOT::RVec<double>& values, const ROOT::RVec<int>& flags) {
   ROOT::RVec<size_t> idx;
   for(size_t i = 0; i < flags.size(); i++) {
      if(flags[i] == 1) idx.push_back(i);
   }
   return ROOT::VecOps::Take(values, idx);
}

double filterByFlagLead(const ROOT::RVec<double>& values, const ROOT::RVec<int>& flags) {
   ROOT::RVec<size_t> idx;
   for(size_t i = 0; i < flags.size(); i++) {
      if(flags[i] == 1) idx.push_back(i);
   }
   return idx.size() > 0 ? values[idx[0]] : -999.;
}
""")


# helper function to apply selection on jet constituent level
ROOT.gInterpreter.Declare("""
    ROOT::RVec<double> maskPAndDedx(const ROOT::RVec<double>& values, 
                                    const ROOT::RVec<double>& p, 
                                    const ROOT::RVec<double>& dedx,
                                    double p_min=0., double dedx_min = 0., double dedx_max=0., double default_val=-9.) {
       ROOT::RVec<double> result = values;
       
       if (p_min == 0. && dedx_min && dedx_max == 0.){
          return result;
       }

       for(size_t i = 0; i < p.size(); i++) {
          if(p[i] <= p_min || dedx[i] <= dedx_min ||  dedx[i] >= dedx_max  ) {
             result[i] = default_val;
          }
       }
       return result;
    }
    """)

# helper function to select jet constituents passing some momentum and dedx cuts:
# def mask_branches_jet_constituents(branches_to_mask, p_min, dedx_max ):
#     for branch


def addOverflowToLastBin(hist):
    # print("Adding overflow to last bin for histogram:", hist.GetName())
    overflow_index = hist.GetNbinsX()+1
    if hist.IsBinOverflow(overflow_index):
        # print ("Found overflow bin with index", overflow_index, "merging content with previous bin!")
        overflow = hist.GetBinContent(overflow_index)
        lastbin = hist.GetBinContent(overflow_index-1)
        hist.SetBinContent(overflow_index-1, overflow+lastbin)
        # print ("Added together overflow =", overflow, "and last bin content", lastbin, "in histogram of name", hist.GetName())



def usable_chunks(input_filepath, tree_name="events"):
    # completed chunk files with a non-empty tree; a chunk whose input had no selected event has none,
    # and a job that stopped with an error leaves no eventsSelected
    usable = []
    for f in glob.glob(os.path.join(input_filepath, "chunk*.root")):
        tf = ROOT.TFile.Open(f)
        if not tf or tf.IsZombie():
            print(f"File {f} is not a valid ROOT file. Skipping.")
            continue
        tree = tf.Get(tree_name)
        if not tree or tree.GetEntries() == 0:
            print(f"Tree '{tree_name}' not found or empty in file {f}. Skipping.")
            continue
        if not tf.GetListOfKeys().Contains("eventsSelected"):
            print(f"File {f} is incomplete (no eventsSelected). Skipping.")
            continue
        usable.append(f)
    return usable

def is_valid_rdf(input_filepath, tree_name="events"):

    # Check if input is a directory and if yes, if it has a usable chunk
    if os.path.isdir(input_filepath):
        if not usable_chunks(input_filepath, tree_name):
            print(f"No usable files found in {input_filepath}. Skipping.")
            return False
    # Handle the case where the path is a single file
    else:
        if not os.path.isfile(input_filepath):
            print("Filepath ", input_filepath, "not found")
            return False 

        tf = ROOT.TFile.Open(input_filepath)

        if not tf or tf.IsZombie():
            print(f"File {input_filepath} is not a valid ROOT file. Skipping.")
            return False 

        if not tf.Get(tree_name):
            print(f"Tree '{tree_name}' not found in file {input_filepath}. Skipping.")
            return False 
    
    return True

def get_rdf(input_filepath, tree_name="events"):

    rdf = None

    # Detect if a single ROOT file (<input>.root) exists.
    # If yes, use that; otherwise, assume the input is a directory with chunks.

    # Check if a single ROOT file exists
    single_file = input_filepath + ".root"

    if os.path.isfile(single_file):
        print(f"Detected single ROOT file: {single_file}")

        if not is_valid_rdf(single_file, tree_name):
            print("Process ", single_file, " is empty. Skipping.")
            return None
        
        rdf = ROOT.RDataFrame(tree_name, single_file)

    # assume dir with chunks otherwise
    else:
        print(f"Assuming directory with chunk files: {input_filepath}")

        if not is_valid_rdf(input_filepath, tree_name):
            print("Process ", input_filepath, " is empty or broken. Skipping.")
            return None

        # rdf = ROOT.RDataFrame("events", input_filepath+"/*/chunk*") #allow for top level dir.

        # Get list of chunks allowing for inconvenient FCCAna extra top level dir
        chunk_files = usable_chunks(input_filepath, tree_name)
        if not chunk_files:
            raise FileNotFoundError(f"No chunk files found under {input_filepath}")

        print(f"Found {len(chunk_files)} chunk files")
        rdf = ROOT.RDataFrame(tree_name, chunk_files)  # Pass the list directly

        # check rdf from chunks is ok
        if not rdf:
            print("Empty file for:", input_filepath, " Skipping.")
            return None

    # print(rdf.GetColumnNames())

    return rdf

def get_param_value_from_file(param_name, filepath):
    # Extract the the parameter requested either from single file or chunks case 

    param_value = 1.0  # default fallback
    single_file = filepath + ".root"

    # Case 1: single ROOT file
    if os.path.isfile(single_file):
        resolved_path = single_file
        print(f"Detected single ROOT file: {resolved_path}")

        inputfile = ROOT.TFile.Open(resolved_path, "READ")
        if not inputfile or inputfile.IsZombie():
            print(f"File {resolved_path} is not a valid ROOT file. Cannot read '{param_name}'.")
            return param_value

        if not inputfile.GetListOfKeys().Contains(param_name):
            raise KeyError(f"TParameter '{param_name}' not found in file '{resolved_path}'")

        param_value = inputfile.Get(param_name).GetVal()
        inputfile.Close()

    # Case 2: directory → handle chunk files
    else:
        print(f"Detected directory with chunks: {filepath}")

        
        chunk_files = glob.glob(os.path.join(filepath, "chunk*.root"))
        if not chunk_files:
            raise FileNotFoundError(f"No chunk ROOT files found in directory: {filepath}")

        param_total = 0.0
        for chunk in chunk_files:
            f = ROOT.TFile.Open(chunk, "READ")
            if not f or f.IsZombie():
                print(f"Warning: Skipping invalid chunk file {chunk}")
                continue

            if not f.GetListOfKeys().Contains(param_name):
                print(f"Warning: TParameter '{param_name}' not found in {chunk}, skipping.")
                f.Close()
                continue

            val = f.Get(param_name).GetVal()
            param_total += val
            f.Close()

        param_value = param_total
        print(f"Summed '{param_name}' from {len(chunk_files)} chunks: {param_value}")

    return param_value


# def get_param_value_from_file(param_name, filepath):
#     param_value = 1.

#     inputfile = ROOT.TFile.Open(filepath, "READ")
#     if not inputfile or inputfile.IsZombie():
#         print(f"File {filepath} is not a valid ROOT file. Cannot read paramValue for normalisation. Skipping.")
#         return param_value 
    
#     if not inputfile.GetListOfKeys().Contains(param_name):
#         raise KeyError(f"TParameter '{param_name}' not found in file '{filepath}'")
    
#     param_value = inputfile.Get(param_name).GetVal()

#     inputfile.Close()
#     return param_value

def get_norm_factor(sample, filepath, norm_file, weighted=False, lumi=1., print_debug=False):

    norm_factor = -1.

    # Read the cross-section and its correction factors from the .json file
    if '.json' in norm_file:
        with open(norm_file, 'r') as dict_file:
            proct_dict = json.load(dict_file)
    else:
        raise Exception("Error in get_norm_factor: Only support .json files for the cross-section values currently!")

    # For the total number of events generated or sum of weights, read them from the file (input needs to be ntuple produced with FCCAnalyses!)
    if weighted:
        total_sow = get_param_value_from_file("SumOfWeights", filepath)
        norm_factor = lumi*proct_dict[sample]["crossSection"]*proct_dict[sample]["kfactor"]*proct_dict[sample]["matchingEfficiency"]/total_sow
    else:
        total_nevts = get_param_value_from_file("eventsProcessed", filepath)
        norm_factor = lumi*proct_dict[sample]["crossSection"]*proct_dict[sample]["kfactor"]*proct_dict[sample]["matchingEfficiency"]/total_nevts

    if print_debug:
        print(f"{sample} - built normfactor from:")
        if weighted:
            print("sum of weights = ", total_sow)
        else:
            print("total nevts = ", total_nevts)

        print("x-sec = ", proct_dict[sample]["crossSection"])
        print("k-factor = ", proct_dict[sample]["kfactor"])
        print("match-eff = ", proct_dict[sample]["matchingEfficiency"])
        print("lumi = ", lumi)
        print("norm_factor = ", norm_factor)

    if norm_factor < 0:
        raise Exception("Error in get_norm_factor: Probably missing a factor in the normalisation ...")

    return norm_factor


def get_hist_from_tree(proc_name, input_filepath, plot_specs, norm_file=None,
                       lumi=1., is_data= False, add_overflow=False, weighted=False, tree_name="events",
                       selection=None, branches_to_mask=None):

    if is_data and weighted:
        raise Exception("Incompatible options in get_hist_from_tree() : is_data and weighted (MC weights) both true.")
    
    if not norm_file and not is_data:
        raise Exception("Error in get_hist_from_tree() on MC: no json file with normalisation info provided!")

    #pass the full plot_specs and then break down
    hist_var = plot_specs.name
    hist_xmin = plot_specs.xmin
    hist_xmax = plot_specs.xmax
    hist_nbins = plot_specs.nbins

    #going to need a TH1Model to fill:
    has_variable_binning = False
    if not isinstance(hist_nbins, int):
        has_variable_binning = True
        hist_binEdges = array("d", hist_nbins)
        hist_nBins = len(hist_nbins)-1
        #init the histogram with variable bin widths:
        hist_model = ROOT.RDF.TH1DModel(f"hist_{proc_name}_{hist_var}", f"hist_{proc_name}_{hist_var}", hist_nBins, hist_binEdges)

    else:
        hist_model = ROOT.RDF.TH1DModel(f"hist_{proc_name}_{hist_var}", f"hist_{proc_name}_{hist_var}", hist_nbins, hist_xmin, hist_xmax)

    # get the dataframe and fill the histogram model
    rdf = get_rdf(input_filepath, tree_name)
    if rdf is None:
        # nothing usable in this sample: an empty histogram, which the caller skips
        return hist_model.GetHistogram()

    # mask branches for jet constituent selection if requested:
    if selection and branches_to_mask:
        for branch in branches_to_mask:
            new_branch_names = (f"{branch}_jet1_masked", f"{branch}_jet2_masked")
            print(f"Going to mask branch {branch} for jet one and jet two with names {new_branch_names[0]}, {new_branch_names[1]}")
            rdf = rdf.Define(new_branch_names[0], f"maskPAndDedx({branch}[0], pfcand_p[0], pfcand_dEdx_wires_value[0], {selection.p_min}, {selection.dedx_min}, {selection.dedx_max}, -9.)")
            rdf = rdf.Define(new_branch_names[1], f"maskPAndDedx({branch}[1], pfcand_p[1], pfcand_dEdx_wires_value[1], {selection.p_min}, {selection.dedx_min}, {selection.dedx_max}, -9.)")

    #if needed, define a new column so we can plot it:
    if plot_specs.redefine:
        print(f"Defining new column for requested variable {plot_specs.name} as {plot_specs.redefine}")
        rdf = rdf.Define(f"{plot_specs.name}", f"{plot_specs.redefine}")


    tmp_hist = rdf.Histo1D(hist_model, hist_var).GetValue() # DO NOT USE WEIGHTS

    if add_overflow:
        addOverflowToLastBin(tmp_hist)

    # scale MC with xsec and lumi
    if not is_data:
        norm_factor = get_norm_factor(proc_name, input_filepath, norm_file, weighted=weighted, lumi=lumi )
        tmp_hist.Scale(norm_factor)

    print("{} \t {:.2f}".format(proc_name, tmp_hist.Integral()))

    hist_to_add = copy.deepcopy(tmp_hist)
    hist_to_add.SetTitle(proc_name)
    return hist_to_add

# function to draw all the things #todo: add back support for applying extra cut?
def make_plot(plot_name, plot, input_dir_data, input_dir_mc, data_proc, mc_processes, out_dir_base, 
              year="", sel_tag ="", lumi=1., ecm=1., norm_file=None,
              do_log_y=True, add_overflow=False, fix_ratio_range=(), weighted=False,
              out_format = ".png", store_root_file=False, tree_name="events", 
              selection=None, branches_to_mask=None, 
              ):

    print("Plotting", plot.name)

    if weighted:
        print("Applying generator event weights to MC")
    
    #get data histogram
    data_hist = None
    for data_file in data_proc.sample_list:
        data_filepath = os.path.join(input_dir_data, f"{data_file}")
        print(f"Looking for data file: {data_file} in {data_filepath}")
        data_tmp_hist = get_hist_from_tree(data_file, data_filepath, plot, 
                                           is_data= True, add_overflow=add_overflow, tree_name=tree_name, selection=selection, branches_to_mask=branches_to_mask) #dont need a norm file, nor weights for data

        #skip if nothing passes selection:
        if not data_tmp_hist.GetEntries():
            print("Empty histogram for:", data_file, " Skipping.")
            continue

        #set plotting properties:
        data_tmp_hist.SetTitle(data_proc.title)
        data_tmp_hist.SetLineColor(data_proc.colour_key)
        data_tmp_hist.SetMarkerStyle(ROOT.kFullCircle)

        # make summed up histogram of all samples in the process
        if not data_hist: 
            data_hist = copy.deepcopy(data_tmp_hist)
            data_hist.SetName("data_hist")
        else:
            data_hist.Add(data_tmp_hist)

    # get MC prediction
    mc_hists_list = []

    for proc in mc_processes:
        #loop over files in our mc process and get each histogram, then add them together
        mc_proc_hist = None
        for sample in mc_processes[proc].sample_list:
            filepath = os.path.join(input_dir_mc, f"{sample}")
            print(f"Looking for MC sample: {sample} in {filepath}")
            
            tmp_hist = get_hist_from_tree(sample, filepath, plot, norm_file=norm_file,
                                          lumi=lumi, is_data= False, add_overflow=add_overflow, weighted=weighted, tree_name=tree_name, selection=selection, branches_to_mask=branches_to_mask)

            #skip if nothing passes selection:
            if not tmp_hist.GetEntries():
                print("Empty histogram for:", sample, " Skipping.")
                continue

            #set plotting properties:
            tmp_hist.SetTitle(mc_processes[proc].title)
            tmp_hist.SetLineColor(mc_processes[proc].colour_key)
            tmp_hist.SetFillColor(mc_processes[proc].colour_key)

            # make summed up histogram of all samples in the process
            if not mc_proc_hist: 
                mc_proc_hist = copy.deepcopy(tmp_hist)
                mc_proc_hist.SetName(f"{proc}_hist")
            else:
                mc_proc_hist.Add(tmp_hist)
            
        # add the total histogram for this process to the MC stack 
        mc_hists_list.append(mc_proc_hist)
    

    #sort the MC by size:
    mc_hists_list = sorted(mc_hists_list, key=lambda hist: hist.Integral())
    mc_stack = ROOT.THStack("mc_stack","")
    for mc_hist in mc_hists_list:
        mc_stack.Add(mc_hist)

    #draw stuff:
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800) 
    #add ratio:
    pad_up = ROOT.TPad("pad_up", "pad_up", 0., 0., 1., 1.)
    pad_up.SetFillStyle(0)
    pad_up.SetBottomMargin(0.32)
    pad_up.SetTopMargin(0.03)
    pad_up.SetLeftMargin(0.13)
    pad_up.SetRightMargin(0.05)
    pad_up.SetLogy(do_log_y)
    pad_up.Draw()

    pad_low = ROOT.TPad("pad_low", "pad_low", 0., 0., 1., 1.);
    pad_low.SetFillStyle(0)
    pad_low.SetBottomMargin(0.12)
    pad_low.SetTopMargin(0.72)
    pad_low.SetLeftMargin(0.13)
    pad_low.SetRightMargin(0.05)
    pad_low.SetGrid()
    pad_low.Draw()

    #
    pad_up.cd()

    mc_stack.Draw("hist")
    data_hist.Draw("E0P same")

    #axes:
    mc_stack.SetMinimum(1.e0)
    mc_stack.SetMaximum(1.e10) #MAKE INPUT ARGUMENT
    # mc_stack.SetMaximum(1.e7)
    mc_stack.GetYaxis().SetTitle("Events")
    mc_stack.GetXaxis().SetLabelSize(0)
    data_hist.GetXaxis().SetTitle(plot.label)

    #labels:
    Text = ROOT.TLatex()

    Text.SetNDC() 

    Text.SetTextAlign(12);
    Text.SetNDC(ROOT.kTRUE) 
    Text.SetTextSize(0.025) 
    Text.DrawLatex(0.17, 0.95, "#it{ALEPH data and simulation}") 
    lumi_tag = "#bf{{ #sqrt{{s}} = {:.0f} GeV, L = {:.2f} pb^{{-1}} }}".format(ecm, lumi)
    Text.DrawLatex(0.17, 0.9, lumi_tag)
    ana_tag = "#bf{{ {} , {} }}".format(year, sel_tag)
    Text.DrawLatex(0.17, 0.85, ana_tag)
    if selection:
        const_sel_tag = "#bf{{ {} }}".format(selection.label)
        Text.DrawLatex(0.17, 0.80, const_sel_tag)

    #legend
    legsize = 0.05*((len(mc_hists_list)+1)/2)
    # leg = ROOT.TLegend(0.58,0.86 - legsize,0.86,0.88)
    # leg.Draw()
    leg = ROOT.TLegend(0.6, 0.95 - legsize, 0.95, 0.95)
    for mc_hist in reversed(mc_hists_list):
        leg.AddEntry(mc_hist, mc_hist.GetTitle(), "f")
    leg.AddEntry(data_hist, data_hist.GetTitle(), "ep")

    leg.SetFillStyle( 0 )
    leg.SetBorderSize( 0 )
    leg.SetTextFont( 43 )
    leg.SetTextSize( 22 )
    leg.SetNColumns( 2 )
    leg.SetColumnSeparation(-0.05)
    leg.Draw()

    #ratio panel
    hist_ratio = data_hist.Clone()
    hist_ratio.Divide( mc_stack.GetStack().Last() )
    hist_ratio.GetYaxis().SetTitle("data/MC")
    hist_ratio.GetYaxis().SetTitleOffset(1.95)
    hist_ratio.GetYaxis().SetNdivisions(6)

    #set some fixed ratio: 
    if fix_ratio_range:
        hist_ratio.GetYaxis().SetRangeUser(fix_ratio_range[0], fix_ratio_range[1])

    pad_low.cd()
    pad_low.Update()
    hist_ratio.Draw("E0P")
    pad_low.RedrawAxis()

    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()	

    if not os.path.exists(out_dir_base):
        os.makedirs(out_dir_base)

    if weighted:
        filename = f"hists_{plot_name}_weighted"
    else:
        filename = f"hists_{plot_name}"
    
    if do_log_y:
        filename+="_logY"
    
    if add_overflow:
        filename+="_addedOverflow"
        
    fileout = os.path.join(out_dir_base, f"{filename}{out_format}")

    canvas.SaveAs(fileout)

    #write out a root file for use in e.g. WS building if requested:
    if store_root_file:
        out_filepath = os.path.join(out_dir_base, f"histograms_{plot_name}.root")

        fileout = ROOT.TFile(out_filepath, "RECREATE")
        fileout.cd()

        data_hist.Write()

        for mc_hist in mc_hists_list:
            mc_hist.Write()
        
        # add the metadata as TParameters for book-keeping
        lumi_param = ROOT.TParameter(type(lumi).__name__)("Luminosity", lumi)
        lumi_param.Write()

        ecm_param = ROOT.TParameter(type(ecm).__name__)("Energy", ecm)
        ecm_param.Write()


        print("Wrote root file with histograms to:", out_filepath)

        fileout.Close()
    
def load_config_module(module_name: str):
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError:
        print(f"ERROR: Config module '{module_name}' not found.")
        sys.exit(1)
     

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-c", "--config", type=str, required=True, help="Name of the plotting config module (e.g. 'plotting_config_stage1')")
    args = parser.parse_args()

    #check the config file has a PlottingConfig class, and if yes import
    config_module = load_config_module(args.config)

    if not hasattr(config_module, "PlottingConfig"):
        print(f"ERROR: Module '{config_module.__name__}' does not contain a 'PlottingConfig' class.")
        sys.exit(1)

    PlottingConfig = config_module.PlottingConfig

    #enable multithreading if requested in config:
    if PlottingConfig.n_threads:
        print(f"Enabling Multithreading in RDF with {PlottingConfig.n_threads} cores")
        ROOT.EnableImplicitMT(PlottingConfig.n_threads)

    # loop over all plots : check the config imported for details
    for plot_name, plot_specs in PlottingConfig.plots_dict.items():
        make_plot(plot_name, plot_specs, PlottingConfig.inputs_path_data, PlottingConfig.inputs_path_mc, PlottingConfig.data, PlottingConfig.mc_processes, PlottingConfig.outputs_path, 
              year=PlottingConfig.year, sel_tag =PlottingConfig.sel_tag, lumi=PlottingConfig.lumi, ecm=PlottingConfig.ecm, norm_file=PlottingConfig.norm_file,
              do_log_y=PlottingConfig.do_log_y, add_overflow=PlottingConfig.add_overflow, fix_ratio_range=PlottingConfig.ratio_range, 
              weighted=PlottingConfig.weighted, out_format=PlottingConfig.out_format, store_root_file=PlottingConfig.store_root_file,
              tree_name=PlottingConfig.tree_name, selection=PlottingConfig.selection, branches_to_mask=PlottingConfig.branches_to_mask,
           )
    

