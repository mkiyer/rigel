#!/usr/bin/env bash
# post_release.sh — Finish a Rigel release after PyPI publish succeeds.
#
# Usage:
#   ./scripts/publishing/post_release.sh 0.4.0
#
# What it does (stage 2 of the release pipeline):
#   1. Polls https://pypi.org/pypi/rigel-rnaseq/<version>/json until the sdist
#      is live (timeout ~30 min).
#   2. Fetches the sdist SHA256 and patches conda/meta.yaml.
#   3. Commits the recipe update and pushes to origin/main (no new tag).
#   4. If `--conda-publish` is passed (or you answer yes at the prompt), also
#      runs conda_publish.sh to build and upload to your personal anaconda.org
#      channel.
#
# Bioconda auto-bump usually opens a PR against bioconda-recipes within 24h
# of the PyPI upload. If it doesn't, use the snippet printed at the end.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

die()  { echo -e "${RED}Error:${NC} $1" >&2; exit 1; }
info() { echo -e "${GREEN}✓${NC} $1"; }
warn() { echo -e "${YELLOW}⚠${NC} $1"; }
note() { echo -e "${BLUE}→${NC} $1"; }

sed_inplace() {
    if sed --version >/dev/null 2>&1; then
        sed -i "$@"
    else
        sed -i '' "$@"
    fi
}

# ── Args ─────────────────────────────────────────────────────────────────
AUTO_CONDA_PUBLISH=0
POSITIONAL=()
for arg in "$@"; do
    case "$arg" in
        --conda-publish) AUTO_CONDA_PUBLISH=1 ;;
        --skip-conda-publish) AUTO_CONDA_PUBLISH=-1 ;;
        -h|--help)
            sed -n '2,20p' "$0"; exit 0 ;;
        *) POSITIONAL+=("$arg") ;;
    esac
done

[[ ${#POSITIONAL[@]} -eq 1 ]] || die "Usage: $0 <version> [--conda-publish|--skip-conda-publish]"
VERSION="${POSITIONAL[0]}"
TAG="v${VERSION}"

[[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || die "Version must be X.Y.Z (got: $VERSION)"

cd "$(git rev-parse --show-toplevel)" || die "Not in a git repository"

git tag -l "$TAG" | grep -q . || warn "Tag $TAG does not exist locally — did release.sh run?"

# Working tree must be clean (we're going to commit).
if [[ -n "$(git status --porcelain)" ]]; then
    die "Working tree is dirty. Commit or stash before running post_release."
fi

# ── Poll PyPI ────────────────────────────────────────────────────────────
PYPI_URL="https://pypi.org/pypi/rigel-rnaseq/${VERSION}/json"
note "Waiting for rigel-rnaseq==${VERSION} on PyPI …"

MAX_ATTEMPTS=60   # 60 * 30s = 30 min
SLEEP_SECS=30
PYPI_JSON=""
for ((i = 1; i <= MAX_ATTEMPTS; i++)); do
    if PYPI_JSON=$(curl -sf "$PYPI_URL" 2>/dev/null); then
        info "Found on PyPI (attempt ${i})"
        break
    fi
    printf "  [%2d/%2d] not yet live, sleeping %ds…\r" "$i" "$MAX_ATTEMPTS" "$SLEEP_SECS"
    sleep "$SLEEP_SECS"
done
echo ""
[[ -n "$PYPI_JSON" ]] || die "Timed out after $((MAX_ATTEMPTS * SLEEP_SECS / 60)) min waiting for PyPI.
Check the workflow:  https://github.com/mkiyer/rigel/actions"

# ── Extract sdist sha256 ─────────────────────────────────────────────────
SHA256=$(echo "$PYPI_JSON" | python3 -c '
import json, sys
data = json.load(sys.stdin)
for item in data.get("urls", []):
    if item.get("packagetype") == "sdist":
        print(item["digests"]["sha256"])
        sys.exit(0)
sys.exit(1)
') || die "Could not find sdist SHA256 in PyPI response"
info "sha256: ${SHA256}"

# ── Patch conda/meta.yaml ────────────────────────────────────────────────
META="conda/meta.yaml"
[[ -f "$META" ]] || die "$META not found"

# Verify version matches (release.sh should already have bumped it).
grep -q "{% set version = \"${VERSION}\" %}" "$META" || \
    die "$META version is not ${VERSION}. Did release.sh run?"

# Replace any existing sha256 line (including placeholder comment).
if grep -Eq '^[[:space:]]*#?[[:space:]]*sha256:' "$META"; then
    sed_inplace -E "s|^([[:space:]]*)#?[[:space:]]*sha256:.*|\1sha256: ${SHA256}|" "$META"
else
    die "Could not locate a sha256 line in $META"
fi
grep -q "sha256: ${SHA256}" "$META" || die "Failed to patch $META"
info "Patched $META"

echo ""
echo "Recipe update:"
echo "──────────────"
git --no-pager diff "$META"
echo ""

# ── Commit + push ────────────────────────────────────────────────────────
read -rp "Commit and push conda recipe update to origin/main? (y/N) " answer
if ! [[ "$answer" =~ ^[Yy]$ ]]; then
    warn "Skipping git commit. Recipe is patched but uncommitted."
else
    git add "$META"
    git commit -m "Release ${TAG}: update conda sha256"
    git push origin main
    info "Pushed recipe update"
fi

# ── Optional: personal anaconda.org upload ───────────────────────────────
RUN_CONDA_PUB=0
if [[ $AUTO_CONDA_PUBLISH -eq 1 ]]; then
    RUN_CONDA_PUB=1
elif [[ $AUTO_CONDA_PUBLISH -eq 0 ]]; then
    echo ""
    read -rp "Build and upload to your personal anaconda.org channel now? (y/N) " answer
    [[ "$answer" =~ ^[Yy]$ ]] && RUN_CONDA_PUB=1
fi

if [[ $RUN_CONDA_PUB -eq 1 ]]; then
    echo ""
    note "Running conda_publish.sh (requires bioconda-build env + anaconda login)…"
    "$(dirname "$0")/conda_publish.sh"
fi

# ── Next steps ───────────────────────────────────────────────────────────
cat <<EOF

${GREEN}PyPI release ${VERSION} finalized.${NC}

  PyPI:       https://pypi.org/project/rigel-rnaseq/${VERSION}/
  Install:    pip install rigel-rnaseq==${VERSION}

${BLUE}Bioconda:${NC} the auto-bump bot usually opens a PR within 24h. To track or
force an update manually:

  https://github.com/bioconda/bioconda-recipes/pulls?q=rigel

If the bot doesn't fire, run this to open a manual PR:

  git clone https://github.com/mkiyer/bioconda-recipes.git
  cd bioconda-recipes
  git checkout -b update-rigel-${VERSION}
  mkdir -p recipes/rigel
  cp ../rigel/conda/meta.yaml recipes/rigel/meta.yaml
  git add recipes/rigel/meta.yaml
  git commit -m "Update rigel to ${VERSION}"
  git push origin update-rigel-${VERSION}
  # then open a PR at https://github.com/bioconda/bioconda-recipes/pulls

EOF
