# Sources the key4hep stack pinned in FCCAnalyses/.fccana/stack_pin, or the
# latest key4hep release if that file does not exist (it is not part of the
# repository, see README.md), unless a key4hep stack is already set up in the
# shell, and sets up all FCCAnalyses environment variables.
# To set or update the pin, edit FCCAnalyses/.fccana/stack_pin.
# Also puts src/ on PYTHONPATH for the shared python modules (run_list).
export PYTHONPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/src${PYTHONPATH:+:${PYTHONPATH}}"
source $(dirname "${BASH_SOURCE[0]}")/FCCAnalyses/setup.sh
