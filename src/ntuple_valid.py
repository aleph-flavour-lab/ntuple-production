#!/usr/bin/env python3

import sys
import ROOT

# helper function to read branches

def get_value(obj):
    if hasattr(obj, "__len__") and not isinstance(obj, (int, float)):
        if len(obj) == 0:
            return None
        return obj[0]
    return obj


# helper function which returns the unique run & event_number combinations in a file as python set
def load_events(filename, tree_name, run_branch, event_branch, do_flavour_filter = False, flavour_val = 5):
    file = ROOT.TFile.Open(filename)
    if not file or file.IsZombie():
        raise RuntimeError(f"Cannot open file: {filename}")

    tree = file.Get(tree_name)
    if not tree:
        raise RuntimeError(f"Cannot find tree '{tree_name}' in {filename}")

    events = set()

    for entry in tree:
        run_vec = getattr(entry, run_branch)
        evt_vec = getattr(entry, event_branch)

        # Safety checks
        if len(run_vec) == 0 or len(evt_vec) == 0:
            continue

        if len(run_vec) != 1 or len(evt_vec) != 1:
            raise RuntimeError(
                f"Expected exactly 1 element per event, "
                f"got {len(run_vec)} (run) and {len(evt_vec)} (event)"
            )
        
        # luka's files aren't filtered by flavour, ours are, add the filter here:
        if do_flavour_filter:
            evt_type = getattr(entry, "genEventType")

            if evt_type != flavour_val:
                continue

        run = int(run_vec[0])
        evt = int(evt_vec[0])

        events.add((int(run), int(evt)))

    file.Close()
    return events

# function that compares the common events branch by branch
# translation of branch names needs to be passed 

def compare_events(filepath1, filepath2,
                   common_events,
                   branch_map,
                   max_print=10, tree_name = "events"):

    #debug:
    print("file1:", filepath1)
    print("file2:", filepath2)

    #open files again
    file1 = ROOT.TFile.Open(filepath1)
    file2 = ROOT.TFile.Open(filepath2)

    # get the tree in each file
    tree1 = file1.Get(tree_name)
    tree2 = file2.Get(tree_name)

    # Build ROOT internal indices
    print("Building indices...")
    tree1.BuildIndex(branch_map["run"][0], branch_map["event"][0])
    tree2.BuildIndex(branch_map["run"][1], branch_map["event"][1])

    n_diff_events = 0
    list_of_branches_with_diffs = []

    for run_number, event_number in common_events:
        entry_number_file1 = tree1.GetEntryNumberWithIndex(run_number, event_number)
        tree1.GetEntry(entry_number_file1)

        entry_number_file2 = tree2.GetEntryNumberWithIndex(run_number, event_number)
        tree2.GetEntry(entry_number_file2)

        event_has_diff = False

        for branch_name, (branch_1, branch_2) in branch_map.items():

            val1 = get_value(getattr(tree1, branch_1))
            val2 = get_value(getattr(tree2, branch_2))

            # print(branch_name, val1, val2) #DEBUG REMOVE LATER

            if val1 != val2:
                if not event_has_diff:
                    if n_diff_events < max_print:
                        print(f"\nDifference in event {run_number} {event_number}")
                    event_has_diff = True
                    n_diff_events += 1

                print(f"  Branch mismatch: {branch_1} vs {branch_2}")
                print(f"    file1 (Luka): {val1}")
                print(f"    file2 (ours): {val2}")

                if not branch_name in list_of_branches_with_diffs:
                    list_of_branches_with_diffs.append(branch_name)
            
            # if n_diff_events >= max_print:
            #     print("\nReached print limit.")
            #     break
    
    print(f"\nTotal events with differences: {n_diff_events}")
    print("Branches with differences found:")
    print(list_of_branches_with_diffs)


def main():

    tree_name = "events"

    # First we get the sets of run and event number from each file and find the overlap (=common events)
    file1 = "/eos/user/l/llambrec/aleph-data/ntuples-withks/eventlevel/mc/output_qqb_1.root"
    run_branch_1 = "runNumber"
    event_branch_1 = "eventNumber"

    print("Checking events in file {}".format(file1))
    # Note: Luka's file is not filtered by flavour, ours is, so apply the flavour filter here
    events1 = load_events(file1, tree_name, run_branch_1, event_branch_1, do_flavour_filter = True, flavour_val = 5.)
    print(f"File1: {len(events1)} events")

    # file2 = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/v09_ntuple_valid/ntuple_valid_tester_5.root" 
    file2 = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/v15_SV_test/Zbb.root" 
    run_branch_2 = "run_number"
    event_branch_2 = "event_number"

    print("Checking events in file {}".format(file2))
    events2 = load_events(file2, tree_name, run_branch_2, event_branch_2, do_flavour_filter = False, flavour_val = 5.)
    print(f"File2: {len(events2)} events")

    only_in_1 = events1 - events2
    only_in_2 = events2 - events1
    common = events1 & events2

    print("\nComparison results:")
    print(f"Common events: {len(common)}")
    print(f"Only in file1: {len(only_in_1)}") # This should be 0 !!!
    print(f"Only in file2: {len(only_in_2)}") # Note: out file processed larger fraction of input, so will be  > 0


    # print(common)

    # now check the common events branch by branch

    # names are not unified, so need translation map:
    branch_translation = {
        #"name":"(name_file1, name_file2)",
        "run":("runNumber", "run_number"),
        "event":("eventNumber", "event_number"),
        #input to the primary vertex fit:
        # "n_selected_tracks":("Event_nSelectedTracks", "n_tracks_sel"),
        # # "n_selected_tracks_vertex":("", "n_trackstates_sel"), #luka doesnt store this?
        # # output of primary vertex fit
        # "n_primary_tracks":("Event_nPrimaryTracks", "n_primary_tracks"),
        # "n_secondary_tracks":("Event_nSecondaryTracks", "n_secondary_tracks"),
        # # vertex position
        # "vertex_x":("PV_x", "Vertex_refit_x"),
        # "vertex_y":("PV_y", "Vertex_refit_y"),
        # "vertex_z":("PV_z", "Vertex_refit_z"),
        # #truth vertex
        # "gen_vertex_x":("GenPV_x", "gen_vertex_x"),
        # "gen_vertex_y":("GenPV_y", "gen_vertex_y"),
        # "gen_vertex_z":("GenPV_z", "gen_vertex_z"),
        #secondary vertices:
        "n_sv_event":("Event_nSV", "n_sv_event"),
        "n_sv_jets":("Jets_nSV", "n_sv_jets"),
        "sv_ntracks":("SecondaryVertices_nTracks", "sv_ntracks"),
        "n_v0_event":("Event_nV0Candidates", "n_v0_event"),



    }

    compare_events(file1, file2, common, branch_translation)
    print(f"Reminder: Common events: {len(common)}")



if __name__ == "__main__":
    main()
