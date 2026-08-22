#!/usr/bin/env bash
set -euo pipefail

if [ -z "${1:-}" ]; then
    echo "Usage: $0 <output-file>" >&2
    exit 1
fi

OUTPUT_FILE="$1"
[[ "$OUTPUT_FILE" != *.pdf ]] && OUTPUT_FILE="${OUTPUT_FILE}.pdf"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIGURES_DIR="${REPO_ROOT}/results/figures"

if [ ! -d "$FIGURES_DIR" ]; then
    echo "Error: Directory '${FIGURES_DIR}' does not exist." >&2
    exit 1
fi

shopt -s nullglob
files=()
for f in "${FIGURES_DIR}"/*.pdf; do
    [ -f "$f" ] || continue
    # Skip destination output file if it already exists in the figures directory
    if [ -f "$OUTPUT_FILE" ] && [ "$(realpath "$f")" = "$(realpath "$OUTPUT_FILE")" ]; then
        continue
    fi
    files+=("$f")
done

if [ ${#files[@]} -eq 0 ]; then
    echo "No PDF figures found to merge in ${FIGURES_DIR}." >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"

echo "Merging ${#files[@]} figures into ${OUTPUT_FILE}..."
qpdf --empty --pages "${files[@]}" -- "$OUTPUT_FILE"

echo "Removing old figure files..."
rm -f "${files[@]}"

echo "Successfully merged figures into ${OUTPUT_FILE}."
