#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_FILE="${ROOT_DIR}/bin/csfem"

print_help() {
    cat <<'EOF'
Usage:
  bash csfem_run.sh [linear_patch|nonlinear_patch|cantilever|bending_block] [options]

Examples:
  bash csfem_run.sh linear_patch
  bash csfem_run.sh nonlinear_patch
  bash csfem_run.sh nonlinear_patch --num-els 4 --nstep 20
EOF
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
    linear_patch|nonlinear_patch|cantilever|bending_block)
        bash "${ROOT_DIR}/csfem_build.sh" --ensure
        if [[ "${CASE_NAME}" == "bending_block" && "${has_solver_option}" -eq 0 && "${has_maxiter_option}" -eq 0 ]]; then
            "${BIN_FILE}" bending_block --solver sparselu --maxiter 80 "$@"
        elif [[ "${CASE_NAME}" == "bending_block" && "${has_solver_option}" -eq 0 ]]; then
            "${BIN_FILE}" bending_block --solver sparselu "$@"
        elif [[ "${CASE_NAME}" == "nonlinear_patch" && "${has_solver_option}" -eq 0 ]]; then
            "${BIN_FILE}" nonlinear_patch --solver sparselu "$@"
        else
            "${BIN_FILE}" "${CASE_NAME}" "$@"
        fi
        ;;
    *)
        echo "Unsupported case: ${CASE_NAME}" >&2
        print_help >&2
        exit 1
        ;;
esac
