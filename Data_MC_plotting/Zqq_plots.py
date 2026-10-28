from collections import namedtuple
import ROOT

# some plots need masking of branches, to enable selection of jet constituents
masked_branches_dict ={
	"Zqq_data_MC_vars_dEdx":[
		"pfcand_p",
		"pfcand_dEdx_wires_value",
        "pfcand_dEdx_wires_pull_ele",
        "pfcand_dEdx_wires_pull_mu",
        "pfcand_dEdx_wires_pull_pi",
        "pfcand_dEdx_wires_pull_kaon",
        "pfcand_dEdx_wires_pull_proton"
	]
}

#variables to plot:
default_binning = 25 
PlotSpecs = namedtuple('PlotSpecs', ['name', 'xmin', 'xmax', 'label', 'nbins', 'redefine'])

#binning for dE/dx plots:
xmin_dedx = 0.
xmax_dedx = 3.
bin_width_dedx = 0.1
nbins_dedx = int((xmax_dedx-xmin_dedx)/bin_width_dedx)



Zqq_data_MC_vars_dEdx = {

	# check modelling of jet constituent momentum
	"pfcand_p_all":PlotSpecs(name="pfcand_p_all", xmin=0., xmax=50, label="p in GeV all pfCands in jet1 ", nbins=100, 
										redefine="pfcand_p_jet1_masked" ), 

	# dE/dx wires values 

	# # all pfcandidates
	# "dEdx_wires_val_jet1_all":PlotSpecs(name="dEdx_wires_val_jet1_all", xmin=xmin_dedx, xmax=xmax_dedx, label="dE/dx all pfCands in jet1 (wires) ", nbins=nbins_dedx, 
	# 									redefine="pfcand_dEdx_wires_value_jet1_masked" ),
	# "dEdx_wires_val_jet2_all":PlotSpecs(name="dEdx_wires_val_jet2_all", xmin=xmin_dedx, xmax=xmax_dedx, label="dE/dx all pfCands in jet2 (wires) ", nbins=nbins_dedx, 
	# 									redefine="pfcand_dEdx_wires_value_jet2_masked" ),
	
	# # # leading pfcandidate
	# "dEdx_wires_val_jet1_leading":PlotSpecs(name="dEdx_wires_val_jet1_leading_pfCand", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx, 
	# 										label="dE/dx leading pfCand in jet1 (wires) ", 
	# 										redefine="pfcand_dEdx_wires_value_jet1_masked[0]" ), 
	
	# "dEdx_wires_val_jet2_leading":PlotSpecs(name="dEdx_wires_val_jet2_leading_pfCand", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx, 
	# 										label="dE/dx leading pfCand in jet2 (wires) ", 
	# 										redefine="pfcand_dEdx_wires_value_jet2_masked[0]" ), 
	
	# # # all electrons
	# "dEdx_wires_val_jet1_electrons":PlotSpecs(name="dEdx_wires_val_jet1_electrons", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
	# 										  label="dE/dx all electrons in jet1 (wires) ",  
	# 										  redefine="filterByFlag(pfcand_dEdx_wires_value_jet1_masked, pfcand_isEl[0])" ),
	
	# "dEdx_wires_val_jet2_electrons":PlotSpecs(name="dEdx_wires_val_jet2_electrons", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
	# 										  label="dE/dx all electrons in jet2 (wires) ",  
	# 										  redefine="filterByFlag(pfcand_dEdx_wires_value_jet2_masked, pfcand_isEl[1])" ),
	
	# # leading electron
	# "dEdx_wires_val_jet1_electron_lead":PlotSpecs(name="dEdx_wires_val_jet1_electron_lead", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
	# 											 label="dE/dx leading electron in jet1 (wires) ",  
	# 											 redefine="filterByFlagLead(pfcand_dEdx_wires_value_jet1_masked, pfcand_isEl[0])" ),

	# "dEdx_wires_val_jet2_electron_lead":PlotSpecs(name="dEdx_wires_val_jet2_electron_lead", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
	# 											 label="dE/dx leading electron in jet2 (wires) ",  
	# 											 redefine="filterByFlagLead(pfcand_dEdx_wires_value_jet2_masked, pfcand_isEl[1])" ),

	# # # all muons
	# "dEdx_wires_val_jet1_muons": PlotSpecs(name="dEdx_wires_val_jet1_muons", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
    #                                    label="dE/dx all muons in jet1 (wires)",  
    #                                    redefine="filterByFlag(pfcand_dEdx_wires_value_jet1_masked, pfcand_isMu[0])"),

	# "dEdx_wires_val_jet2_muons": PlotSpecs(name="dEdx_wires_val_jet2_muons", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
	# 									label="dE/dx all muons in jet2 (wires)",  
	# 									redefine="filterByFlag(pfcand_dEdx_wires_value_jet2_masked, pfcand_isMu[1])"),
	
	# # leading muon
	# "dEdx_wires_val_jet1_muon_lead": PlotSpecs(name="dEdx_wires_val_jet1_muon_lead", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
    #                                        label="dE/dx leading muon in jet1 (wires)",  
    #                                        redefine="filterByFlagLead(pfcand_dEdx_wires_value_jet1_masked, pfcand_isMu[0])"),

	# "dEdx_wires_val_jet2_muon_lead": PlotSpecs(name="dEdx_wires_val_jet2_muon_lead", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
	# 										label="dE/dx leading muon in jet2 (wires)",  
	# 										redefine="filterByFlagLead(pfcand_dEdx_wires_value_jet2_masked, pfcand_isMu[1])"),

	
	# # # all charged hadrons 
	# "dEdx_wires_val_jet1_chargedHad": PlotSpecs(name="dEdx_wires_val_jet1_chargedHad", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
    #                                         label="dE/dx all charged hadrons in jet1 (wires)",  
    #                                         redefine="filterByFlag(pfcand_dEdx_wires_value_jet1_masked, pfcand_isChargedHad[0])"),

	# "dEdx_wires_val_jet2_chargedHad": PlotSpecs(name="dEdx_wires_val_jet2_chargedHad", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
	# 											label="dE/dx all charged hadrons in jet2 (wires)",  
	# 											redefine="filterByFlag(pfcand_dEdx_wires_value_jet2_masked, pfcand_isChargedHad[1])"),

	# # leading charged hadron
	# "dEdx_wires_val_jet1_chargedHad_lead": PlotSpecs(name="dEdx_wires_val_jet1_chargedHad_lead", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
    #                                              label="dE/dx leading charged hadron in jet1 (wires)",  
    #                                              redefine="filterByFlagLead(pfcand_dEdx_wires_value_jet1_masked, pfcand_isChargedHad[0])"),

	# "dEdx_wires_val_jet2_chargedHad_lead": PlotSpecs(name="dEdx_wires_val_jet2_chargedHad_lead", xmin=xmin_dedx, xmax=xmax_dedx, nbins=nbins_dedx,
	# 												label="dE/dx leading charged hadron in jet2 (wires)",  
	# 												redefine="filterByFlagLead(pfcand_dEdx_wires_value_jet2_masked, pfcand_isChargedHad[1])"),

	# # # one tester plot on neutral hadrons to be sure the plotting logic works 

	# "DEBUG_dEdx_wires_val_jet1_neutralHad": PlotSpecs(name="dEdx_wires_val_jet1_neutralHad", xmin=-10, xmax=15, nbins=nbins_dedx,
    #                                         label="dE/dx all neutral hadrons in jet1 (wires)",  
    #                                         redefine="filterByFlag(pfcand_dEdx_wires_value_jet1_masked, pfcand_isNeutralHad[0])"),


	# dE/dx wires errors # TODO !!!!!

	# #all pfcandidates
	# "dEdx_wires_err_jet1_all":PlotSpecs(name="dEdx_wires_err_jet1_all", xmin=0, xmax=15, label="#Delta dE/dx all pfCands in jet1 (wires) ", nbins=default_binning, redefine="pfcand_dEdx_wires_error[0]" ),
	# "dEdx_wires_err_jet2_all":PlotSpecs(name="dEdx_wires_err_jet2_all", xmin=0, xmax=15, label="#Delta dE/dx all pfCands in jet2 (wires) ", nbins=default_binning, redefine="pfcand_dEdx_wires_error[1]" ),
	
	# # leading pfcandidate
	# "dEdx_wires_err_jet1_leading":PlotSpecs(name="dEdx_wires_err_jet1_leading_pfCand", xmin=0, xmax=15, nbins=default_binning, 
	# 										label="#Delta dE/dx leading pfCand in jet1 (wires) ", 
	# 										redefine="pfcand_dEdx_wires_error[0][0]" ), 
	
	# "dEdx_wires_err_jet2_leading":PlotSpecs(name="dEdx_wires_err_jet2_leading_pfCand", xmin=0, xmax=15, nbins=default_binning, 
	# 										label="dE/dx leading pfCand in jet2 (wires) ", 
	# 										redefine="pfcand_dEdx_wires_error[1][0]" ), 


	# PID p-values from Bethe Bloch fits

	#all pfcandidates
	# "PID_pval_electron_wires_jet1_all":PlotSpecs(name="pfcand_PID_pval_wires_ele_jet1_all", xmin=-1.1, xmax=1.1, label="p-val electron all pfCands in jet1 (wires) ", nbins=default_binning, 
	# 									redefine="pfcand_PID_pval_wires_ele[0]" ),
	# "PID_pval_electron_wires_jet2_all":PlotSpecs(name="PID_pval_electron_wires_jet2_all", xmin=-1.1, xmax=1.1, label="p-val electron all pfCands in jet2 (wires) ", nbins=default_binning, 
	# 									redefine="pfcand_PID_pval_wires_ele[1]" ),
	
	# # leading pfcandidate
	# "PID_pval_electron_wires_jet1_leading":PlotSpecs(name="pfcand_PID_pval_wires_ele_jet1_leading", xmin=-1.1, xmax=1.1, nbins=default_binning, 
	# 										label="p-val electron leading pfCand in jet1 (wires) ", 
	# 										redefine="pfcand_PID_pval_wires_ele[0][0]" ), 
	
	# "PID_pval_electron_wires_jet2_leading":PlotSpecs(name="PID_pval_electron_wires_jet2_leading", xmin=-1.1, xmax=1.1, nbins=default_binning, 
	# 										label="p-val electron leading pfCand in jet2 (wires) ", 
	# 										redefine="pfcand_PID_pval_wires_ele[1][0]" ), 

	# # # all electrons

	# "PID_pval_electron_wires_jet1_electrons": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet1_electrons", xmin=-1.1, xmax=1.1, nbins=default_binning,
    #                                                 label="p-val electron electrons in jet1 (wires)",  
    #                                                 redefine="filterByFlag(pfcand_PID_pval_wires_ele[0], pfcand_isEl[0])"),

													

	# "PID_pval_electron_wires_jet2_electrons": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet2_electrons", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 													label="p-val electron electrons in jet2 (wires)",  
	# 													redefine="filterByFlag(pfcand_PID_pval_wires_ele[1], pfcand_isEl[1])"),


	# # # leading electron
	# "PID_pval_electron_wires_jet1_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet1_electron_leading", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 														label="p-val electron leading electron in jet1 (wires)",  
	# 														redefine="filterByFlagLead(pfcand_PID_pval_wires_ele[0], pfcand_isEl[0])"),

	# "PID_pval_electron_wires_jet2_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet2_electron_leading", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 														label="p-val electron leading electron in jet2 (wires)",  
	# 														redefine="filterByFlagLead(pfcand_PID_pval_wires_ele[0], pfcand_isEl[0])"),

	# # # all muons
	# "PID_pval_electron_wires_jet1_muons": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet1_muons", xmin=-1.1, xmax=1.1, nbins=default_binning,
    #                                             label="p-val electron muons in jet1 (wires)",  
    #                                             redefine="filterByFlag(pfcand_PID_pval_wires_ele[0], pfcand_isMu[0])"),

	# "PID_pval_electron_wires_jet2_muons": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet2_muons", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 												label="p-val electron muons in jet2 (wires)",  
	# 												redefine="filterByFlag(pfcand_PID_pval_wires_ele[1], pfcand_isMu[1])"),
	# # # leading muon
	# "PID_pval_electron_wires_jet1_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet1_muon_leading", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 													label="p-val electron leading muon in jet1 (wires)",  
	# 													redefine="filterByFlagLead(pfcand_PID_pval_wires_ele[0], pfcand_isMu[0])"),

	# "PID_pval_electron_wires_jet2_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet2_muon_leading", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 													label="p-val electron leading muon in jet2 (wires)",  
	# 													redefine="filterByFlagLead(pfcand_PID_pval_wires_ele[1], pfcand_isMu[1])"),

	# # # all charged hadrons
	# "PID_pval_electron_wires_jet1_chargedHad": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet1_chargedHad", xmin=-1.1, xmax=1.1, nbins=default_binning,
    #                                                  label="p-val electron charged hadrons in jet1 (wires)",  
    #                                                  redefine="filterByFlag(pfcand_PID_pval_wires_ele[0], pfcand_isChargedHad[0])"),

	# "PID_pval_electron_wires_jet2_chargedHad": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet2_chargedHad", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 													label="p-val electron charged hadrons in jet2 (wires)",  
	# 													redefine="filterByFlag(pfcand_PID_pval_wires_ele[1], pfcand_isChargedHad[1])"),
	# # leading charged hadron
	# "PID_pval_electron_wires_jet1_chargedHad_leading": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet1_chargedHad_leading", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 															label="p-val electron leading charged hadron in jet1 (wires)",  
	# 															redefine="filterByFlagLead(pfcand_PID_pval_wires_ele[0], pfcand_isChargedHad[0])"),

	# "PID_pval_electron_wires_jet2_chargedHad_leading": PlotSpecs(name="pfcand_PID_pval_wires_ele_jet2_chargedHad_leading", xmin=-1.1, xmax=1.1, nbins=default_binning,
	# 															label="p-val electron leading charged hadron in jet2 (wires)",  
	# 															redefine="filterByFlagLead(pfcand_PID_pval_wires_ele[1], pfcand_isChargedHad[1])"),

	# p-value for muons

	# all pfcandidates
	# "PID_pval_muon_wires_jet1_all": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet1_all", xmin=-1.1, xmax=1.1, 
	# 	label="p-val muon all pfCands in jet1 (wires)", nbins=default_binning, 
	# 	redefine="pfcand_PID_pval_wires_mu[0]"),

	# "PID_pval_muon_wires_jet2_all": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet2_all", xmin=-1.1, xmax=1.1, 
	# 	label="p-val muon all pfCands in jet2 (wires)", nbins=default_binning, 
	# 	redefine="pfcand_PID_pval_wires_mu[1]"),

	# # leading pfcandidate
	# "PID_pval_muon_wires_jet1_leading": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet1_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon leading pfCand in jet1 (wires)", 
	# 	redefine="pfcand_PID_pval_wires_mu[0][0]"),

	# "PID_pval_muon_wires_jet2_leading": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet2_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon leading pfCand in jet2 (wires)", 
	# 	redefine="pfcand_PID_pval_wires_mu[1][0]"),

	# # all electrons
	# "PID_pval_muon_wires_jet1_electrons": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet1_electrons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon electrons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_mu[0], pfcand_isEl[0])"),

	# "PID_pval_muon_wires_jet2_electrons": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet2_electrons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon electrons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_mu[1], pfcand_isEl[1])"),

	# # leading electron
	# "PID_pval_muon_wires_jet1_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet1_electron_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val muon leading electron in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_mu[0], pfcand_isEl[0])"),

	# "PID_pval_muon_wires_jet2_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet2_electron_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val muon leading electron in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_mu[1], pfcand_isEl[1])"),

	# # all muons
	# "PID_pval_muon_wires_jet1_muons": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet1_muons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon muons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_mu[0], pfcand_isMu[0])"),

	# "PID_pval_muon_wires_jet2_muons": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet2_muons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon muons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_mu[1], pfcand_isMu[1])"),

	# # leading muon
	# "PID_pval_muon_wires_jet1_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet1_muon_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon leading muon in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_mu[0], pfcand_isMu[0])"),

	# "PID_pval_muon_wires_jet2_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet2_muon_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon leading muon in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_mu[1], pfcand_isMu[1])"),

	# # all charged hadrons
	# "PID_pval_muon_wires_jet1_chargedHad": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet1_chargedHad", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon charged hadrons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_mu[0], pfcand_isChargedHad[0])"),

	# "PID_pval_muon_wires_jet2_chargedHad": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet2_chargedHad", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val muon charged hadrons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_mu[1], pfcand_isChargedHad[1])"),

	# # leading charged hadron
	# "PID_pval_muon_wires_jet1_chargedHad_leading": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet1_chargedHad_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val muon leading charged hadron in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_mu[0], pfcand_isChargedHad[0])"),

	# "PID_pval_muon_wires_jet2_chargedHad_leading": PlotSpecs(name="pfcand_PID_pval_wires_mu_jet2_chargedHad_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val muon leading charged hadron in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_mu[1], pfcand_isChargedHad[1])"),
	
	# # pions p-values

	# # all pfcandidates
	# "PID_pval_pion_wires_jet1_all": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet1_all", xmin=-1.1, xmax=1.1, 
	# 	label="p-val pion all pfCands in jet1 (wires)", nbins=default_binning, 
	# 	redefine="pfcand_PID_pval_wires_pi[0]"),

	# "PID_pval_pion_wires_jet2_all": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet2_all", xmin=-1.1, xmax=1.1, 
	# 	label="p-val pion all pfCands in jet2 (wires)", nbins=default_binning, 
	# 	redefine="pfcand_PID_pval_wires_pi[1]"),

	# # leading pfcandidate
	# "PID_pval_pion_wires_jet1_leading": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet1_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion leading pfCand in jet1 (wires)", 
	# 	redefine="pfcand_PID_pval_wires_pi[0][0]"),

	# "PID_pval_pion_wires_jet2_leading": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet2_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion leading pfCand in jet2 (wires)", 
	# 	redefine="pfcand_PID_pval_wires_pi[1][0]"),

	# # all electrons
	# "PID_pval_pion_wires_jet1_electrons": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet1_electrons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion electrons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_pi[0], pfcand_isEl[0])"),

	# "PID_pval_pion_wires_jet2_electrons": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet2_electrons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion electrons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_pi[1], pfcand_isEl[1])"),

	# # leading electron
	# "PID_pval_pion_wires_jet1_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet1_electron_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val pion leading electron in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_pi[0], pfcand_isEl[0])"),

	# "PID_pval_pion_wires_jet2_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet2_electron_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val pion leading electron in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_pi[1], pfcand_isEl[1])"),

	# # all muons
	# "PID_pval_pion_wires_jet1_muons": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet1_muons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion muons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_pi[0], pfcand_isMu[0])"),

	# "PID_pval_pion_wires_jet2_muons": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet2_muons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion muons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_pi[1], pfcand_isMu[1])"),

	# # leading muon
	# "PID_pval_pion_wires_jet1_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet1_muon_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion leading muon in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_pi[0], pfcand_isMu[0])"),

	# "PID_pval_pion_wires_jet2_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet2_muon_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion leading muon in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_pi[1], pfcand_isMu[1])"),

	# # all charged hadrons
	# "PID_pval_pion_wires_jet1_chargedHad": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet1_chargedHad", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion charged hadrons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_pi[0], pfcand_isChargedHad[0])"),

	# "PID_pval_pion_wires_jet2_chargedHad": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet2_chargedHad", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val pion charged hadrons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_pi[1], pfcand_isChargedHad[1])"),

	# # leading charged hadron
	# "PID_pval_pion_wires_jet1_chargedHad_leading": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet1_chargedHad_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val pion leading charged hadron in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_pi[0], pfcand_isChargedHad[0])"),

	# "PID_pval_pion_wires_jet2_chargedHad_leading": PlotSpecs(name="pfcand_PID_pval_wires_pi_jet2_chargedHad_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val pion leading charged hadron in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_pi[1], pfcand_isChargedHad[1])"),

	# ## kaons p-values
	# # all pfcandidates
	# "PID_pval_kaon_wires_jet1_all": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet1_all", xmin=-1.1, xmax=1.1, 
	# 	label="p-val kaon all pfCands in jet1 (wires)", nbins=default_binning, 
	# 	redefine="pfcand_PID_pval_wires_kaon[0]"),

	# "PID_pval_kaon_wires_jet2_all": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet2_all", xmin=-1.1, xmax=1.1, 
	# 	label="p-val kaon all pfCands in jet2 (wires)", nbins=default_binning, 
	# 	redefine="pfcand_PID_pval_wires_kaon[1]"),

	# # leading pfcandidate
	# "PID_pval_kaon_wires_jet1_leading": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet1_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val kaon leading pfCand in jet1 (wires)", 
	# 	redefine="pfcand_PID_pval_wires_kaon[0][0]"),

	# "PID_pval_kaon_wires_jet2_leading": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet2_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val kaon leading pfCand in jet2 (wires)", 
	# 	redefine="pfcand_PID_pval_wires_kaon[1][0]"),

	# # all electrons
	# "PID_pval_kaon_wires_jet1_electrons": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet1_electrons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val kaon electrons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_kaon[0], pfcand_isEl[0])"),

	# "PID_pval_kaon_wires_jet2_electrons": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet2_electrons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val kaon electrons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_kaon[1], pfcand_isEl[1])"),

	# # leading electron
	# "PID_pval_kaon_wires_jet1_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet1_electron_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val kaon leading electron in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_kaon[0], pfcand_isEl[0])"),

	# "PID_pval_kaon_wires_jet2_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet2_electron_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val kaon leading electron in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_kaon[1], pfcand_isEl[1])"),

	# # all muons
	# "PID_pval_kaon_wires_jet1_muons": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet1_muons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val kaon muons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_kaon[0], pfcand_isMu[0])"),

	# "PID_pval_kaon_wires_jet2_muons": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet2_muons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val kaon muons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_kaon[1], pfcand_isMu[1])"),

	# # leading muon
	# "PID_pval_kaon_wires_jet1_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet1_muon_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val kaon leading muon in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_kaon[0], pfcand_isMu[0])"),

	# "PID_pval_kaon_wires_jet2_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_kaon_jet2_muon_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val kaon leading muon in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_kaon[1], pfcand_isMu[1])"),

	# all charged hadrons
	"dEdx_pull_kaon_wires_jet1_chargedHad": PlotSpecs(name="pfcand_dEdx_wires_pull_kaon_jet1_chargedHad", xmin=-5., xmax=5.,
		nbins=default_binning, label="kaon dE/dx pull charged hadrons in jet1 (wires)",
		redefine="filterByFlag(pfcand_dEdx_wires_pull_kaon[0], pfcand_isChargedHad[0])"),

	"dEdx_pull_kaon_wires_jet2_chargedHad": PlotSpecs(name="pfcand_dEdx_wires_pull_kaon_jet2_chargedHad", xmin=-5., xmax=5.,
		nbins=default_binning, label="kaon dE/dx pull charged hadrons in jet2 (wires)",
		redefine="filterByFlag(pfcand_dEdx_wires_pull_kaon[1], pfcand_isChargedHad[1])"),

	# leading charged hadron
	"dEdx_pull_kaon_wires_jet1_chargedHad_leading": PlotSpecs(name="pfcand_dEdx_wires_pull_kaon_jet1_chargedHad_leading",
		xmin=-5., xmax=5., nbins=default_binning, label="kaon dE/dx pull leading charged hadron in jet1 (wires)",
		redefine="filterByFlagLead(pfcand_dEdx_wires_pull_kaon[0], pfcand_isChargedHad[0])"),

	"dEdx_pull_kaon_wires_jet2_chargedHad_leading": PlotSpecs(name="pfcand_dEdx_wires_pull_kaon_jet2_chargedHad_leading",
		xmin=-5., xmax=5., nbins=default_binning, label="kaon dE/dx pull leading charged hadron in jet2 (wires)",
		redefine="filterByFlagLead(pfcand_dEdx_wires_pull_kaon[1], pfcand_isChargedHad[1])"),
	
	# ## proton p-values
	# # all pfcandidates
	# "PID_pval_proton_wires_jet1_all": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet1_all", xmin=-1.1, xmax=1.1, 
	# 	label="p-val proton all pfCands in jet1 (wires)", nbins=default_binning, 
	# 	redefine="pfcand_PID_pval_wires_proton[0]"),

	# "PID_pval_proton_wires_jet2_all": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet2_all", xmin=-1.1, xmax=1.1, 
	# 	label="p-val proton all pfCands in jet2 (wires)", nbins=default_binning, 
	# 	redefine="pfcand_PID_pval_wires_proton[1]"),

	# # leading pfcandidate
	# "PID_pval_proton_wires_jet1_leading": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet1_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton leading pfCand in jet1 (wires)", 
	# 	redefine="pfcand_PID_pval_wires_proton[0][0]"),

	# "PID_pval_proton_wires_jet2_leading": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet2_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton leading pfCand in jet2 (wires)", 
	# 	redefine="pfcand_PID_pval_wires_proton[1][0]"),

	# # all electrons
	# "PID_pval_proton_wires_jet1_electrons": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet1_electrons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton electrons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_proton[0], pfcand_isEl[0])"),

	# "PID_pval_proton_wires_jet2_electrons": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet2_electrons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton electrons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_proton[1], pfcand_isEl[1])"),

	# # leading electron
	# "PID_pval_proton_wires_jet1_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet1_electron_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val proton leading electron in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_proton[0], pfcand_isEl[0])"),

	# "PID_pval_proton_wires_jet2_electron_leading": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet2_electron_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val proton leading electron in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_proton[1], pfcand_isEl[1])"),

	# # all muons
	# "PID_pval_proton_wires_jet1_muons": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet1_muons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton muons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_proton[0], pfcand_isMu[0])"),

	# "PID_pval_proton_wires_jet2_muons": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet2_muons", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton muons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_proton[1], pfcand_isMu[1])"),

	# # leading muon
	# "PID_pval_proton_wires_jet1_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet1_muon_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton leading muon in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_proton[0], pfcand_isMu[0])"),

	# "PID_pval_proton_wires_jet2_muon_leading": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet2_muon_leading", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton leading muon in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_proton[1], pfcand_isMu[1])"),

	# # all charged hadrons
	# "PID_pval_proton_wires_jet1_chargedHad": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet1_chargedHad", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton charged hadrons in jet1 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_proton[0], pfcand_isChargedHad[0])"),

	# "PID_pval_proton_wires_jet2_chargedHad": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet2_chargedHad", xmin=-1.1, xmax=1.1, 
	# 	nbins=default_binning, label="p-val proton charged hadrons in jet2 (wires)",  
	# 	redefine="filterByFlag(pfcand_PID_pval_wires_proton[1], pfcand_isChargedHad[1])"),

	# # leading charged hadron
	# "PID_pval_proton_wires_jet1_chargedHad_leading": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet1_chargedHad_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val proton leading charged hadron in jet1 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_proton[0], pfcand_isChargedHad[0])"),

	# "PID_pval_proton_wires_jet2_chargedHad_leading": PlotSpecs(name="pfcand_PID_pval_wires_proton_jet2_chargedHad_leading", 
	# 	xmin=-1.1, xmax=1.1, nbins=default_binning, label="p-val proton leading charged hadron in jet2 (wires)",  
	# 	redefine="filterByFlagLead(pfcand_PID_pval_wires_proton[1], pfcand_isChargedHad[1])"),


#




	

	# event multiplicities & kinematics
	# "n_jets":PlotSpecs(name="event_njet", xmin=0, xmax=5, label="Number of jets", nbins=5),
	# "m_jj":PlotSpecs(name="event_invariant_mass", xmin=50, xmax=120, label="m_{jj} in GeV", nbins=70),

	# # jet properties
	# "jet_mass":PlotSpecs(name="jet_mass", xmin=0, xmax=60, label="Jet mass in GeV", nbins=60),
	# "jet_pT":PlotSpecs(name="jet_pT", xmin=0, xmax=60, label="p_{T} jet in GeV", nbins=60),
	# # "jet_eta":PlotSpecs(name="jet_eta", xmin=-2, xmax=2, label="#eta jet", nbins=30),
	# "jet_phi":PlotSpecs(name="jet_phi", xmin=0, xmax=7, label="#phi jet", nbins=30),

	# # jet constituents multiplicities
	# "jet_nconst":PlotSpecs(name="jet_nconst", xmin=0, xmax=50, label="Number of constituents in jet", nbins=50),
	# "jet_nmu":PlotSpecs(name="jet_nmu", xmin=0, xmax=5, label="Number of muons in jet", nbins=5),
	# "jet_nel":PlotSpecs(name="jet_nel", xmin=0, xmax=5, label="Number of electrons in jet", nbins=5),
	# "jet_ngamma":PlotSpecs(name="jet_ngamma", xmin=0, xmax=20, label="Number of photons in jet", nbins=20),
	# "jet_nchad":PlotSpecs(name="jet_nchad", xmin=0, xmax=30, label="Number of charged hadrons in jet", nbins=30),
	# "jet_nnhad":PlotSpecs(name="jet_nnhad", xmin=0, xmax=20, label="Number of neutral hadrons in jet", nbins=20),

	# all input features of the tagger

	
}

Zqq_data_MC_vars_stage2 = {

	# check dE/dx and pvals
	# "dEdx_wires_val_jet1_leading":PlotSpecs(name="pfcand_dEdx_wires_value_jet1_masked[0]", xmin=0, xmax=15, label="dE/dx leading pfCand in jet1 (wires) ", nbins=20),
	# "dEdx_wires_err_allcands_jet1_leading":PlotSpecs(name="pfcand_dEdx_wires_error[0][0", xmin=0, xmax=15, label="Error dE/dx leading pfCand in jet1 (wires) ", nbins=20),

	# event multiplicities & kinematics
	# "n_jets":PlotSpecs(name="event_njet", xmin=0, xmax=5, label="Number of jets", nbins=5),
	# "m_jj":PlotSpecs(name="event_invariant_mass", xmin=50, xmax=120, label="m_{jj} in GeV", nbins=70),

	# # jet properties
	# "jet_mass":PlotSpecs(name="jet_mass", xmin=0, xmax=60, label="Jet mass in GeV", nbins=60),
	# "jet_pT":PlotSpecs(name="jet_pT", xmin=0, xmax=60, label="p_{T} jet in GeV", nbins=60),
	# # "jet_eta":PlotSpecs(name="jet_eta", xmin=-2, xmax=2, label="#eta jet", nbins=30),
	# "jet_phi":PlotSpecs(name="jet_phi", xmin=0, xmax=7, label="#phi jet", nbins=30),

	# # jet constituents multiplicities
	# "jet_nconst":PlotSpecs(name="jet_nconst", xmin=0, xmax=50, label="Number of constituents in jet", nbins=50),
	# "jet_nmu":PlotSpecs(name="jet_nmu", xmin=0, xmax=5, label="Number of muons in jet", nbins=5),
	# "jet_nel":PlotSpecs(name="jet_nel", xmin=0, xmax=5, label="Number of electrons in jet", nbins=5),
	# "jet_ngamma":PlotSpecs(name="jet_ngamma", xmin=0, xmax=20, label="Number of photons in jet", nbins=20),
	# "jet_nchad":PlotSpecs(name="jet_nchad", xmin=0, xmax=30, label="Number of charged hadrons in jet", nbins=30),
	# "jet_nnhad":PlotSpecs(name="jet_nnhad", xmin=0, xmax=20, label="Number of neutral hadrons in jet", nbins=20),

	# all input features of the tagger

	
}

Zqq_data_MC_inference = {
	"dijet_score":PlotSpecs(name="dijet_score", xmin=0, xmax=1, label="B-jet score product", nbins=50, redefine=""),
	"dijet_score_coarse":PlotSpecs(name="dijet_score", xmin=0, xmax=1, label="B-jet score product", nbins=25, redefine=""),
	"dijet_score_coarsest":PlotSpecs(name="dijet_score", xmin=0, xmax=1, label="B-jet score product", nbins=10, redefine=""),
	"dijet_score_luka":PlotSpecs(name="dijet_score", xmin=0, xmax=1, label="B-jet score product", nbins=20, redefine=""),
	
}
