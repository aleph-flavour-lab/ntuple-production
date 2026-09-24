# Sources the pinned key4hep stack (recorded in FCCAnalyses/.fccana/stack_pin)
# and sets up all FCCAnalyses environment variables.
# To update the pin, edit FCCAnalyses/.fccana/stack_pin.
# Also puts src/ on PYTHONPATH for the shared python modules (run_list).
export PYTHONPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/src${PYTHONPATH:+:${PYTHONPATH}}"
source $(dirname "${BASH_SOURCE[0]}")/FCCAnalyses/setup.sh
