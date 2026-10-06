"""fccanalysis analysis script of the stage1 CI.

Runs the repository's src/stage1.py unchanged; only the input directory, the output
directory and the number of threads are redirected, from the environment:
  CI_INPUT_DIR   directory holding 1994/*.root (data) and QQB/*.root (MC)
  CI_OUTPUT_DIR  where the stage1 output is written
  CI_NTHREADS    number of threads
Usage: fccanalysis run .github/ci/stage1_ci.py -- <stage1 arguments>
"""
import importlib.util
import os

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src")
_spec = importlib.util.spec_from_file_location("stage1_under_test", os.path.join(SRC, "stage1.py"))
_stage1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_stage1)


class Analysis(_stage1.Analysis):
    def __init__(self, cmdline_args):
        super().__init__(cmdline_args)
        self.input_dir = os.environ["CI_INPUT_DIR"]
        self.output_dir = os.environ["CI_OUTPUT_DIR"]
        self.n_threads = int(os.environ["CI_NTHREADS"])
        # no copy of the output to the production area on EOS
        self.__dict__.pop("output_dir_eos", None)
        # the header paths are relative to src/, not to this script
        self.include_paths = [p if os.path.isabs(p) else os.path.join(SRC, p) for p in self.include_paths]
