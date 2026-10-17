#!/usr/bin/env bash
set -euo pipefail

# All sample, input and output locations are supplied by the caller.
exec hybcapseq_quant "$@"
