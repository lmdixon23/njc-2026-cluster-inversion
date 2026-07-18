#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python code/certify_open_family.py config/witness_common.json
python code/certify_open_family.py config/witness_anisotropic.json
