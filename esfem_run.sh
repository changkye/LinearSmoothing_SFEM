#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_FILE="${ROOT_DIR}/bin/esfem"

print_help() {
    cat <<'EOF'
Usage:
  bash esfem_run.sh [linear_patch|nonlinear_patch|bending_block] [options]

Examples:
  bash esfem_run.sh linear_patch
  bash esfem_run.sh nonlinear_patch
  bash esfem_run.sh linear_patch --num-els 4
  bash esfem_run.sh nonlinear_patch --num-els 4 --nstep 20
EOF
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

case "${CASE_NAME}" in
    linear_patch|nonlinear_patch|bending_block)
        bash "${ROOT_DIR}/esfem_build.sh" --ensure
        "${BIN_FILE}" "${CASE_NAME}" "$@"
        ;;
    *)
        echo "Unsupported case: ${CASE_NAME}" >&2
        print_help >&2
        exit 1
        ;;
esac
