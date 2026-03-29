#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"

read -r -p "Input .geo file name: " input_geo

if [[ -z "${input_geo}" ]]; then
    echo "Error: .geo file name is required."
    exit 1
fi

if [[ "${input_geo}" != *.geo ]]; then
    input_geo="${input_geo}.geo"
fi

if [[ ! -f "${SCRIPT_DIR}/${input_geo}" ]]; then
    echo "Error: File '${input_geo}' not found."
    exit 1
fi

output_msh="${input_geo%.geo}.msh"
output_vtu="${input_geo%.geo}.vtu"
output_csv="${input_geo%.geo}.csv"

if ! command -v gmsh >/dev/null 2>&1; then
    echo "Error: gmsh is not installed or not in PATH."
    exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda is not installed or not in PATH."
    exit 1
fi

echo "Generating quadratic T6 mesh '${output_msh}' from '${input_geo}'..."
gmsh "${SCRIPT_DIR}/${input_geo}" -2 -order 2 -format msh4 -o "${SCRIPT_DIR}/${output_msh}"

echo "Running tri2poly.py with '${output_msh}'..."
(
    cd "${SCRIPT_DIR}"
    conda run -n python_env python tri2poly.py "${output_msh}"
)

mkdir -p "${DATA_DIR}"
if [[ -f "${SCRIPT_DIR}/${output_csv}" ]]; then
    cp "${SCRIPT_DIR}/${output_csv}" "${DATA_DIR}/${output_csv}"
    echo "Copied '${output_csv}' to '${DATA_DIR}'."
else
    echo "Error: Expected polygon CSV '${output_csv}' was not generated."
    exit 1
fi

echo "Generated '${output_msh}' and '${output_vtu}' remain in '${SCRIPT_DIR}'."
