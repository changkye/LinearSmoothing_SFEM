#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_FILE="${ROOT_DIR}/bin/esfem"
COOK_MSH="${ROOT_DIR}/data/Cook.msh"
COOK_GEO="${ROOT_DIR}/data/Cook.geo"
RUN_MESH_SCRIPT="${ROOT_DIR}/data/run_mesh.sh"

print_help() {
    cat <<'EOF'
Usage:
  bash esfem_run.sh [linear_patch|nonlinear_patch|cantilever|bending_block|cook] [options]

Examples:
  bash esfem_run.sh linear_patch
  bash esfem_run.sh nonlinear_patch
  bash esfem_run.sh cook
  bash esfem_run.sh linear_patch --num-els 4
  bash esfem_run.sh nonlinear_patch --num-els 4 --nstep 20
EOF
}

ensure_cook_mesh() {
    if [[ ! -f "${COOK_MSH}" || "${COOK_GEO}" -nt "${COOK_MSH}" || "${RUN_MESH_SCRIPT}" -nt "${COOK_MSH}" ]]; then
        printf 'Cook\n' | "${RUN_MESH_SCRIPT}"
    fi
}

if (($# == 0)); then
    set -- linear_patch
fi

if [[ "${1:-}" == "--build-clean" ]]; then
    bash "${ROOT_DIR}/esfem_build.sh" --build-clean
    exit 0
fi

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    print_help
    exit 0
fi

CASE_NAME="${1:-linear_patch}"
shift
CASE_NAME="$(printf '%s' "${CASE_NAME}" | tr '[:upper:]' '[:lower:]')"

case "${CASE_NAME}" in
    linear_patch|nonlinear_patch|cantilever|bending_block|cook)
        bash "${ROOT_DIR}/esfem_build.sh" --ensure
        if [[ "${CASE_NAME}" == "cook" ]]; then
            ensure_cook_mesh
        fi
        "${BIN_FILE}" "${CASE_NAME}" "$@"
        ;;
    *)
        echo "Unsupported case: ${CASE_NAME}" >&2
        print_help >&2
        exit 1
        ;;
esac
