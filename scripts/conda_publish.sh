# conda_publish.sh — Build and upload a conda package to your personal
# anaconda.org channel (no bioconda review required).
#
# Prerequisites:
#   conda activate bioconda-build   (needs conda-build + anaconda-client)
#   anaconda login                  (one-time login to anaconda.org)
#
# Usage (source it so conda is available):
#   conda activate bioconda-build
#   source scripts/conda_publish.sh
#
# This builds from conda/meta.yaml using the PyPI sdist as the source.
# Make sure the PyPI release exists first (run release.sh + create GitHub Release).

# Use return instead of exit so sourcing doesn't kill the shell
_conda_publish() {
    local RED='\033[0;31m'
    local GREEN='\033[0;32m'
    local NC='\033[0m'

    _die() { echo -e "${RED}Error: $1${NC}" >&2; return 1; }
    _info() { echo -e "${GREEN}✓${NC} $1"; }

    cd "$(git rev-parse --show-toplevel)" || { _die "Not in a git repository"; return 1; }

    # ── Pre-flight checks ────────────────────────────────────────────────
    command -v conda-build >/dev/null 2>&1 || \
        { _die "conda-build not found. Run: conda activate bioconda-build"; return 1; }
    command -v anaconda >/dev/null 2>&1 || \
        { _die "anaconda-client not found. Run: conda activate bioconda-build"; return 1; }

    if ! anaconda whoami >/dev/null 2>&1; then
        _die "Not logged in to anaconda.org. Run: anaconda login"; return 1
    fi
    local ANACONDA_USER
    ANACONDA_USER=$(anaconda whoami 2>/dev/null | grep 'Username:' | awk '{print $2}')
    _info "Logged in to anaconda.org as: ${ANACONDA_USER}"

    if grep -q 'UPDATE_WITH_ACTUAL_HASH' conda/meta.yaml; then
        _die "conda/meta.yaml still has the placeholder SHA256.
Run: source scripts/bioconda_hash.sh <version>"; return 1
    fi

    local VERSION
    VERSION=$(grep '{% set version' conda/meta.yaml | sed 's/.*"\(.*\)".*/\1/')
    _info "Building rigel ${VERSION}"

    # ── Build ────────────────────────────────────────────────────────────
    echo ""
    echo "Building conda package from conda/meta.yaml..."
    echo "─────────────────────────────────────────────────"

    conda-build conda/ \
        --channel conda-forge \
        --channel bioconda \
        --no-anaconda-upload || { _die "conda-build failed"; return 1; }

    local PKG_PATH
    PKG_PATH=$(conda-build conda/ \
        --channel conda-forge \
        --channel bioconda \
        --output 2>/dev/null)

    [[ -f "$PKG_PATH" ]] || { _die "Built package not found at: $PKG_PATH"; return 1; }
    _info "Package built: ${PKG_PATH}"

    # ── Upload ───────────────────────────────────────────────────────────
    echo ""
    echo "Uploading to anaconda.org/${ANACONDA_USER}..."
    echo "─────────────────────────────────────────────────"

    anaconda upload --user "${ANACONDA_USER}" "$PKG_PATH" || { _die "Upload failed"; return 1; }

    _info "Uploaded!"
    echo ""
    echo -e "${GREEN}Done!${NC} Install with:"
    echo ""
    echo "  conda install -c ${ANACONDA_USER} -c conda-forge -c bioconda rigel==${VERSION}"
}

# Run the function, then clean up
_conda_publish
unset -f _conda_publish _die _info
