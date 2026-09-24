"""Per-run event-class counts of the converted ALEPH data files, input of build_run_list.py.

    python count_event_classes.py --out class_counts_raw_1994.json [--inputs GLOB] [--threads N]
"""
import argparse
import glob
import json
import ROOT

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--out", required=True, help="output JSON (keep it outside the repository)")
ap.add_argument("--inputs", default="/eos/experiment/fcc/ee/analyses/case-studies/aleph/LEP1_DATA/1994/data_*.root")
ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--run-min", type=int, default=25500, help="lowest run number of the year (one histogram bin per run)")
ap.add_argument("--run-max", type=int, default=30500, help="one above the highest run number")
a = ap.parse_args()
NRUN = a.run_max - a.run_min

ROOT.EnableImplicitMT(a.threads)
ROOT.gInterpreter.Declare('''
ROOT::VecOps::RVec<int> bits(const ROOT::VecOps::RVec<uint32_t>& cb){ ROOT::VecOps::RVec<int> o; if(cb.empty()) return o; for(int i=0;i<32;i++) if((cb[0]>>i)&1) o.push_back(i+1); return o;}
''')
files = sorted(glob.glob(a.inputs))
print("files", len(files), flush=True)
df = ROOT.RDataFrame("events", files).Filter("EventHeader.runNumber.size() == 1").Define("run", "(int)EventHeader.runNumber[0]").Define("cls", "bits(ClassBitset)").Define("ncls", "(int)cls.size()")
h = df.Histo2D(("h", "h", NRUN, a.run_min, a.run_max, 33, 0, 33), "run", "cls")
hr = df.Histo1D(("hr", "hr", NRUN, a.run_min, a.run_max), "run")
h0 = df.Filter("ncls==0").Histo1D(("h0", "h0", NRUN, a.run_min, a.run_max), "run")
n = df.Count()
N = n.GetValue()
outside = hr.GetBinContent(0) + hr.GetBinContent(NRUN + 1)
if outside > 0:
    raise SystemExit(f"{int(outside)} events have a run number outside [{a.run_min}, {a.run_max}); widen --run-min/--run-max")
out = {"n_events": N, "per_run": {}}
for i in range(1, hr.GetNbinsX() + 1):
    t = hr.GetBinContent(i)
    if t <= 0:
        continue
    r = int(hr.GetBinLowEdge(i))
    out["per_run"][r] = {"total": int(t), "noclass": int(h0.GetBinContent(i)),
                         "cls": {k: int(h.GetBinContent(i, k + 1)) for k in range(1, 33) if h.GetBinContent(i, k + 1) > 0}}
json.dump(out, open(a.out, "w"))
print("events", N, "runs", len(out["per_run"]), "done")
