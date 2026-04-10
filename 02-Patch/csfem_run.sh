#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_FILE="${ROOT_DIR}/bin/csfem"
COOK_MSH="${ROOT_DIR}/data/Cook.msh"
COOK_GEO="${ROOT_DIR}/data/Cook.geo"
RUN_MESH_SCRIPT="${ROOT_DIR}/data/run_mesh.sh"

print_help() {
    cat <<'EOF'
Usage:
  bash csfem_run.sh [linear_patch|cantilever|cook|cook1|cook2|cook3|cook4|cook5] [options]

Examples:
  bash csfem_run.sh linear_patch
  bash csfem_run.sh cantilever
  bash csfem_run.sh cook
  bash csfem_run.sh cook1
EOF
}

ensure_cook_mesh() {
    local variant="${1:-Cook}"
    local msh="${ROOT_DIR}/data/${variant}.msh"
    local geo="${ROOT_DIR}/data/${variant}.geo"
    if [[ ! -f "${msh}" || "${geo}" -nt "${msh}" || "${RUN_MESH_SCRIPT}" -nt "${msh}" ]]; then
        printf '%s\n' "${variant}" | "${RUN_MESH_SCRIPT}"
    fi
}

if (($# == 0)); then
    set -- linear_patch
fi

if [[ "${1:-}" == "--build-clean" ]]; then
    bash "${ROOT_DIR}/csfem_build.sh" --build-clean
    exit 0
fi

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    print_help
    exit 0
fi

CASE_NAME="${1:-linear_patch}"
shift
CASE_NAME="$(printf '%s' "${CASE_NAME}" | tr '[:upper:]' '[:lower:]')"

has_solver_option=0
has_maxiter_option=0
for arg in "$@"; do
    if [[ "${arg}" == "--solver" ]]; then
        has_solver_option=1
    fi
    if [[ "${arg}" == "--maxiter" ]]; then
        has_maxiter_option=1
    fi
done

case "${CASE_NAME}" in
    linear_patch|high_order_patch|cantilever|cook|cook1|cook2|cook3|cook4|cook5)
        bash "${ROOT_DIR}/csfem_build.sh" --ensure
        if [[ "${CASE_NAME}" == cook* ]]; then
            # Capitalise first letter to match .geo/.msh file names (cook1 -> Cook1)
            variant="$(tr '[:lower:]' '[:upper:]' <<< "${CASE_NAME:0:1}")${CASE_NAME:1}"
            ensure_cook_mesh "${variant}"
        fi
        "${BIN_FILE}" "${CASE_NAME}" "$@"
        ;;
    *)
        echo "Unsupported case: ${CASE_NAME}" >&2
        print_help >&2
        exit 1
        ;;
esac
