#!/usr/bin/env bash
# conda_publish.sh — Build and upload a conda package to your personal
# anaconda.org channel (no bioconda review required).
#
# Prerequisites:
#   conda activate bioconda-build   # needs conda-build + anaconda-client
#   anaconda login                  # one-time, saved in ~/.anaconda/config.yaml
#
# Usage:
#   ./scripts/publishing/conda_publish.sh
#
# This builds from conda/meta.yaml using the PyPI sdist as the source, so
# make sure post_release.sh has already patched the sha256.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

die()  { echo -e "${RED}Error:${NC} $1" >&2; exit 1; }
info() { echo -e "${GREEN}✓${NC} $1"; }
warn() { echo -e "${YELLOW}⚠${NC} $1"; }

cd "$(git rev-parse --show-toplevel)" || die "Not in a git repository"

# ── Pre-flight ──────────────────────────────────────────────────────────
command -v conda-build >/dev/null 2>&1 || \
    die "conda-build not found. Run: conda activate bioconda-build"
command -v anaconda >/dev/null 2>&1 || \
    die "anaconda-client not found. Run: conda activate bioconda-build"

anaconda whoami >/dev/null 2>&1 || die "Not logged in to anaconda.org. Run: anaconda login"

ANACONDA_USER=$(anaconda whoami 2>&1 | awk '/Username:/{print $2}')
[[ -n "$ANACONDA_USER" ]] || die "Could not determine anaconda.org username"
info "Logged in to anaconda.org as: ${ANACONDA_USER}"

# Extract version from conda/meta.yaml.
VERSION=$(grep -E '\{% set version' conda/meta.yaml | sed -E 's/.*"([^"]+)".*/\1/')
[[ -n "$VERSION" ]] || die "Could not parse version from conda/meta.yaml"

# Sanity check: meta.yaml must have a real sha256 (not a placeholder).
if grep -Eq 'sha256:[[:space:]]*(UPDATE_WITH_ACTUAL_HASH|PLACEHOLDER|TBD)' conda/meta.yaml; then
    die "conda/meta.yaml has a placeholder sha256.
Run: ./scripts/publishing/post_release.sh ${VERSION}"
fi

info "Building rigel ${VERSION}"

# ── Build ───────────────────────────────────────────────────────────────
OUTPUT_DIR="${PWD}/conda-build-artifacts"
mkdir -p "$OUTPUT_DIR"

echo ""
echo "Building conda package from conda/meta.yaml…"
echo "────────────────────────────────────────────"

conda-build conda/ \
    --channel conda-forge \
    --channel bioconda \
    --python 3.12 \
    --numpy 1.26 \
    --output-folder "$OUTPUT_DIR" \
    --no-anaconda-upload

PKG_PATH=$(find "$OUTPUT_DIR" -type f \
    \( -name "rigel-${VERSION}-*.conda" -o -name "rigel-${VERSION}-*.tar.bz2" \) \
    | head -n 1)

[[ -f "$PKG_PATH" ]] || die "Built package not found in $OUTPUT_DIR"
info "Built: ${PKG_PATH}"

# ── Upload ──────────────────────────────────────────────────────────────
echo ""
echo "Uploading to anaconda.org/${ANACONDA_USER}…"
echo "────────────────────────────────────────────"

anaconda upload --user "$ANACONDA_USER" "$PKG_PATH"
info "Uploaded"

cat <<EOF

${GREEN}Done.${NC} Install with:

  conda install -c ${ANACONDA_USER} -c conda-forge -c bioconda rigel==${VERSION}

EOF
