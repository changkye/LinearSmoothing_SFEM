#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${ROOT_DIR}/build"
BIN_FILE="${ROOT_DIR}/bin/csfem"
REBUILD=0
CLEAN=0
BUILD_CLEAN=0
CLEAN_ALL=0
INTERNAL_ENSURE=0

print_removed_path() {
    local target="$1"
    if [[ -e "${target}" ]]; then
        rm -rf "${target}"
        echo "Removed: ${target}"
    else
        echo "Not found, skipped: ${target}"
    fi
}

print_help() {
    cat <<'EOF'
Usage:
  bash csfem_build.sh
  bash csfem_build.sh --rebuild
  bash csfem_build.sh --build-clean
  bash csfem_build.sh --clean
  bash csfem_build.sh --clean-all
  bash csfem_build.sh --ensure

Options:
  --rebuild    Remove build before rebuilding the CS-FEM executable
  --build-clean Remove build and bin before rebuilding the CS-FEM executable
  --clean      Remove build and the CS-FEM binary
  --clean-all  Remove build, bin, and res
  --ensure     Internal mode. Build only if needed
  --help       Show this help
EOF
}

while (($# > 0)); do
    case "$1" in
        --rebuild)
            REBUILD=1
            shift
            ;;
        --build-clean)
            BUILD_CLEAN=1
            shift
            ;;
        --clean)
            CLEAN=1
            shift
            ;;
        --clean-all)
            CLEAN_ALL=1
            shift
            ;;
        --ensure)
            INTERNAL_ENSURE=1
            shift
            ;;
        --help|-h)
            print_help
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            print_help >&2
            exit 1
            ;;
    esac
done

if [[ "${BUILD_CLEAN}" -eq 1 ]]; then
    print_removed_path "${BUILD_DIR}"
    print_removed_path "${ROOT_DIR}/bin"
    exit 0
fi

if [[ "${CLEAN}" -eq 1 ]]; then
    print_removed_path "${BUILD_DIR}"
    print_removed_path "${BIN_FILE}"
    exit 0
fi

if [[ "${CLEAN_ALL}" -eq 1 ]]; then
    print_removed_path "${BUILD_DIR}"
    print_removed_path "${ROOT_DIR}/bin"
    print_removed_path "${ROOT_DIR}/res"
    exit 0
fi

build_dir_matches_root() {
    local cache_file cached_root
    cache_file="${BUILD_DIR}/CMakeCache.txt"
    if [[ ! -f "${cache_file}" ]]; then
        return 0
    fi

    cached_root="$(sed -n 's#^CMAKE_HOME_DIRECTORY:INTERNAL=##p' "${cache_file}" | head -n 1)"
    [[ -z "${cached_root}" || "${cached_root}" == "${ROOT_DIR}" ]]
}

needs_rebuild() {
    if [[ "${REBUILD}" -eq 1 ]]; then
        return 0
    fi

    if ! build_dir_matches_root; then
        return 0
    fi

    if [[ ! -d "${BUILD_DIR}" || ! -f "${BIN_FILE}" ]]; then
        return 0
    fi

    if find "${ROOT_DIR}/src" "${ROOT_DIR}/include" \
        -type f \( -name '*.cpp' -o -name '*.hpp' \) \
        -newer "${BIN_FILE}" | grep -q .; then
        return 0
    fi

    if [[ "${ROOT_DIR}/CMakeLists.txt" -nt "${BIN_FILE}" ||
          "${ROOT_DIR}/main.cpp" -nt "${BIN_FILE}" ||
          "${ROOT_DIR}/csfem_build.sh" -nt "${BIN_FILE}" ]]; then
        return 0
    fi

    return 1
}

build_project() {
    if [[ "${REBUILD}" -eq 1 ]] || ! build_dir_matches_root; then
        rm -rf "${BUILD_DIR}"
    fi

    cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}"
    cmake --build "${BUILD_DIR}" --target csfem -j
}

if [[ "${INTERNAL_ENSURE}" -eq 1 ]]; then
    if needs_rebuild; then
        build_project
    else
        echo "Using existing executable: ${BIN_FILE}"
    fi
    exit 0
fi

build_project
