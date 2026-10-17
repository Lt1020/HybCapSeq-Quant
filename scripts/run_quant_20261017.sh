#!/usr/bin/env bash
set -euo pipefail

# Legacy filename retained; this entry point now takes the same generic arguments.
exec hybcapseq_quant "$@"
