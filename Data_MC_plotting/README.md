# Data & MC plotting for ALEPH 

Flexible script to make plots with data and MC for the revived ALEPH datasets, all plotting details are set via config files.

## Setup and input preparation

Follow the instructions of the main repo to set up the environment (`source setup.sh` from the repo root: key4hep stack, FCCAnalyses, and `src/` on `PYTHONPATH` for the shared `run_list` module that the configs import). 

Common input files can be found on eos in `/eos/user/h/hfatehi/D0fliped-good/` currently. Otherwise, use the scripts in this repo to produce ` stage 1` or `inference` type of input files. 

## Run 

With a given `plotting_config.py`, run the following command:

```
python plot_data_mc.py -c plotting_config
```

## Existing config files

### Stage 1 

The config file `plotting_config_stage1.py` is used to plot variables from the `stage1` ntuples that are used in the tagger. 

### Inference 

No config for the ntuples written by the tagger inference is included. Their plots are defined in `Zqq_plots.Zqq_data_MC_inference`; to draw them, write a config like `plotting_config_stage1.py` with `plots_dict = Zqq_plots.Zqq_data_MC_inference`, `ana_stage` and the input paths pointing to the inference ntuples, and a `selection` / `branches_to_mask` suited to them (the stage1 config applies a dE/dx constituent selection).

## Config variables

The following variables are set in the config files:

- `inputs_path`: Path to the input ntuples.  
- `outputs_path` : Path where to store the plots. 
- `plots_dict` : This is a dictionary of named tuples specifying the exact plot settings. It is loaded from the helper file `Zqq_plots.py`. Add new dictionaries there as needed. See below for explanation. 
- `data` and `mc_processes`: These define the filenames/samples for data and MC processes plotted. Again they use named dictionaries and are defined in the helper file `Zqq_processes.py`. See below for explanation. 
- `year`: The data taking year, will be added as label on the plots.
- `sel_tag`: A string that specifies the analysis or selection level of events, will be added as label on the plots.
- `lumi` : The luminosity to normalise MC to,  will be added as label on the plots. `plotting_config_stage1.py` takes it from `src/run_list.py`, i.e. the luminosity of the runs kept by stage1 (`data/lumi/run_list_<year>.csv`).
- `ecm` : The center of mass energy,  will be added as label on the plots.
- `norm_file` : This is a `.json` file which contains the normalisation info for each process, so cross-section, k-factor, matching efficiency. Follows the same standard as used by the `EventProducer` and `FCCAnalyses` approach. Note that the values for `numberOfEvents` and `sumOfWeights` are placeholders, since these will be read from the input ntuples, as they are recalculated during each production to account for failed jobs or otherwise missing files. 
- `do_log_y` : Whether to set the y-axis to logarithmic.
- `add_overflow` : If true, the histogram overflow is added to the last bin.
- `ratio_range` : The y-axis range for the data/MC ratio panel. 
- `weighted` : Whether to use event weights or not.
- `out_format` : The output format in which to store the plots.
- `store_root_file` : If a ROOT file of the plotted histograms is written into the `outputs_path` that can be used by e.g. combine for fitting. 


### PlotSpecs config

### ProcessSpecs config

