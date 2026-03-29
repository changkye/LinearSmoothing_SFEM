#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

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

if ! command -v gmsh >/dev/null 2>&1; then
    echo "Error: gmsh is not installed or not in PATH."
    exit 1
fi

echo "Generating '${output_msh}' from '${input_geo}'..."
gmsh "${SCRIPT_DIR}/${input_geo}" -2 -order 2 -format msh4 -o "${SCRIPT_DIR}/${output_msh}"

echo "Generating '${output_vtu}' for ParaView mesh inspection..."
python3 "${SCRIPT_DIR}/msh_to_vtu.py" "${SCRIPT_DIR}/${output_msh}" "${SCRIPT_DIR}/${output_vtu}"

echo "Generated '${output_msh}' and '${output_vtu}' in '${SCRIPT_DIR}'."
echo "The C++ Cook problem reads this .msh file directly."
