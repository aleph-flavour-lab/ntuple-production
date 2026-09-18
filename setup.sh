# Sources the pinned key4hep stack (recorded in FCCAnalyses/.fccana/stack_pin)
# and sets up all FCCAnalyses environment variables.
# To update the pin, edit FCCAnalyses/.fccana/stack_pin.
# Also puts src/ on PYTHONPATH so that the analysis and plotting scripts can import
# the shared modules there (run_list) from any directory.
export PYTHONPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/src${PYTHONPATH:+:${PYTHONPATH}}"
source $(dirname "${BASH_SOURCE[0]}")/FCCAnalyses/setup.sh
