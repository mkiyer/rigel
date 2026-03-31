#!/usr/bin/env bash
# bioconda_hash.sh — Fetch the PyPI sdist SHA256 and update conda/meta.yaml.
#
# Usage:
#   ./scripts/bioconda_hash.sh 0.3.0
#
# Run this AFTER the PyPI publish workflow completes successfully.
# It fetches the SHA256 from the PyPI JSON API and patches conda/meta.yaml.
#
# For the FIRST bioconda submission:
#   1. Fork https://github.com/bioconda/bioconda-recipes
#   2. Copy conda/meta.yaml → recipes/rigel/meta.yaml in your fork
#   3. Open a PR against bioconda/bioconda-recipes
#
# For SUBSEQUENT releases:
#   Bioconda's auto-bump bot watches PyPI and opens PRs automatically.
#   If the bot doesn't fire, manually update recipes/rigel/meta.yaml.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

die() { echo -e "${RED}Error: $1${NC}" >&2; exit 1; }
info() { echo -e "${GREEN}✓${NC} $1"; }

[[ $# -eq 1 ]] || die "Usage: $0 <version>  (e.g. $0 0.3.0)"
VERSION="$1"

cd "$(git rev-parse --show-toplevel)" || die "Not in a git repository"

# ── Fetch SHA256 from PyPI ──────────────────────────────────────────────
echo "Fetching SHA256 for rigel-rnaseq==${VERSION} from PyPI..."

PYPI_JSON=$(curl -sf "https://pypi.org/pypi/rigel-rnaseq/${VERSION}/json") || \
    die "Could not fetch PyPI metadata. Is version ${VERSION} published?"

SHA256=$(echo "$PYPI_JSON" | python3 -c '
import json, sys
data = json.load(sys.stdin)
for item in data["urls"]:
    if item["packagetype"] == "sdist":
        print(item["digests"]["sha256"])
        sys.exit(0)
print("NOT_FOUND", file=sys.stderr)
sys.exit(1)
') || die "Could not find sdist SHA256 in PyPI response"

info "SHA256: ${SHA256}"

# ── Update conda/meta.yaml ──────────────────────────────────────────────
META="conda/meta.yaml"
[[ -f "$META" ]] || die "conda/meta.yaml not found"

# Replace the sha256 comment/placeholder with the real hash
if grep -q '# sha256:' "$META"; then
    # First time: replace the placeholder comment
    sed -i '' "s|# sha256:.*|sha256: ${SHA256}|" "$META"
elif grep -q 'sha256:' "$META"; then
    # Subsequent: replace the existing hash
    sed -i '' "s|sha256:.*|sha256: ${SHA256}|" "$META"
else
    die "Could not find sha256 line in $META"
fi

info "Updated $META with SHA256"

echo ""
echo "Updated conda/meta.yaml:"
echo "─────────────────────────"
head -12 "$META"
echo ""

# ── Check if this is the first submission ────────────────────────────────
echo -e "${GREEN}Next steps for bioconda:${NC}"
echo ""
echo "  FIRST submission (one-time):"
echo "    1. Fork https://github.com/bioconda/bioconda-recipes"
echo "    2. Clone your fork and create a branch:"
echo "       git clone https://github.com/mkiyer/bioconda-recipes.git"
echo "       cd bioconda-recipes"
echo "       git checkout -b add-rigel-${VERSION}"
echo "    3. Copy the recipe:"
echo "       mkdir -p recipes/rigel"
echo "       cp /path/to/rigel/conda/meta.yaml recipes/rigel/meta.yaml"
echo "    4. Validate locally (optional):"
echo "       bioconda-utils lint --packages rigel"
echo "    5. Commit and push:"
echo "       git add recipes/rigel/meta.yaml"
echo "       git commit -m 'Add rigel ${VERSION}'"
echo "       git push origin add-rigel-${VERSION}"
echo "    6. Open a PR at https://github.com/bioconda/bioconda-recipes/pulls"
echo ""
echo "  SUBSEQUENT releases:"
echo "    Bioconda's auto-bump bot usually detects PyPI updates within 24h."
echo "    If not, repeat the steps above with an 'update-rigel-${VERSION}' branch."
