#!/usr/bin/env bash
# release.sh — Bump version, commit, tag, and push to trigger PyPI publishing.
#
# Usage:
#   ./scripts/release.sh 0.3.0
#
# What it does:
#   1. Validates the version argument
#   2. Checks that CHANGELOG.md has an entry for the new version
#   3. Bumps the version in pyproject.toml and conda/meta.yaml
#   4. Commits the version bump
#   5. Tags the commit as vX.Y.Z
#   6. Pushes the commit and tag to origin
#
# After this script completes, go to GitHub and create a Release from the tag.
# That triggers publish.yml → sdist + wheels → PyPI upload.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

die() { echo -e "${RED}Error: $1${NC}" >&2; exit 1; }
info() { echo -e "${GREEN}✓${NC} $1"; }
warn() { echo -e "${YELLOW}⚠${NC} $1"; }

# ── Validate arguments ──────────────────────────────────────────────────
[[ $# -eq 1 ]] || die "Usage: $0 <version>  (e.g. $0 0.3.0)"
VERSION="$1"
TAG="v${VERSION}"

# Validate version format (semver-ish: X.Y.Z)
[[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || \
    die "Version must be in X.Y.Z format, got: $VERSION"

# ── Pre-flight checks ───────────────────────────────────────────────────
cd "$(git rev-parse --show-toplevel)" || die "Not in a git repository"

# Clean working tree?
if [[ -n "$(git status --porcelain)" ]]; then
    warn "Working tree has uncommitted changes."
    echo "  Files with changes:"
    git status --short
    echo ""
    read -rp "Continue anyway? (y/N) " answer
    [[ "$answer" =~ ^[Yy]$ ]] || exit 0
fi

# On main branch?
BRANCH=$(git branch --show-current)
[[ "$BRANCH" == "main" ]] || die "Must be on 'main' branch (currently on '$BRANCH')"

# Tag doesn't already exist?
git tag -l "$TAG" | grep -q . && die "Tag $TAG already exists"

# CHANGELOG.md has an entry for this version?
if ! grep -q "^\#\# \[${VERSION}\]" CHANGELOG.md; then
    die "CHANGELOG.md has no entry for [${VERSION}].
Add a section like:

  ## [${VERSION}] - $(date +%Y-%m-%d)

  ### Added
  - ...

Then re-run this script."
fi

# ── Read current version ────────────────────────────────────────────────
OLD_VERSION=$(grep '^version' pyproject.toml | head -1 | sed 's/.*"\(.*\)"/\1/')
info "Current version: $OLD_VERSION → New version: $VERSION"

# ── Bump versions ───────────────────────────────────────────────────────
# pyproject.toml
sed -i '' "s/^version = \"${OLD_VERSION}\"/version = \"${VERSION}\"/" pyproject.toml
info "Updated pyproject.toml"

# conda/meta.yaml
sed -i '' "s/{% set version = \"${OLD_VERSION}\" %}/{% set version = \"${VERSION}\" %}/" conda/meta.yaml
info "Updated conda/meta.yaml"

# ── Show diff for confirmation ──────────────────────────────────────────
echo ""
echo "Changes to be committed:"
echo "─────────────────────────"
git diff --stat
echo ""
git diff pyproject.toml conda/meta.yaml
echo ""

read -rp "Commit, tag ${TAG}, and push? (y/N) " answer
[[ "$answer" =~ ^[Yy]$ ]] || { warn "Aborted."; exit 0; }

# ── Commit, tag, push ───────────────────────────────────────────────────
git add pyproject.toml conda/meta.yaml CHANGELOG.md
git commit -m "Release ${TAG}"
info "Committed"

git tag -a "$TAG" -m "Release ${TAG}"
info "Tagged ${TAG}"

git push origin main --tags
info "Pushed to origin"

echo ""
echo -e "${GREEN}Done!${NC} Next steps:"
echo ""
echo "  1. Go to: https://github.com/mkiyer/rigel/releases/new?tag=${TAG}"
echo "  2. Set title: Rigel ${VERSION}"
echo "  3. Paste the CHANGELOG entry as the release body"
echo "  4. Click 'Publish release' → triggers publish.yml → PyPI upload"
echo ""
echo "  After PyPI upload succeeds:"
echo "  5. Run: ./scripts/bioconda_hash.sh ${VERSION}"
