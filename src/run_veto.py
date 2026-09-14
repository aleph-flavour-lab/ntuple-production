"""Run veto for the 1994 ALEPH data: stage1.py applies it in data mode
before any selection (--noRunVeto disables it, --excludeRuns adds runs).

Interim list derived from our own per-run checks, to be replaced by the
official ALEPH good-run list once available."""

VETO_RUNS_1994 = {
    # unusable for vertexing: VDET off (15 runs; 26097 also carries a +-9.5 cm
    # z offset per TPC drift half), or sigma(d0) doubled with ITC hit loss (27776)
    "vertexing": (
        25521, 25522, 25527, 25528, 25530, 25531, 25907, 26097,
        26725, 26726, 26745, 26747, 27776, 29321, 29322, 29323,
    ),
    # run-quality criteria: barrel TPC hit deficit (bootstrap z <= -4) or
    # sigma(d0) degradation with an ITC/VDET hit-loss correlate
    "tracking_quality": (
        25520, 25692, 25712, 25781, 25908, 25934, 25953, 26084, 26088,
        26746, 26854, 27518, 28067, 28206, 28530, 28941, 28990, 29116,
        29539, 29545, 29601, 29977, 30079, 30080, 30086, 30095, 30159,
        30160, 30177, 30188, 30322, 30340, 30347,
    ),
    # solenoid off-plateau: 1.5176 T (+1.2 %) from the current word and from the
    # track curvature; compensation-coil words at 0.4 % of their usual values
    "field": (28125,),
}


def vetoed_runs():
    """Sorted list of all vetoed run numbers."""
    return sorted(set().union(*VETO_RUNS_1994.values()))
