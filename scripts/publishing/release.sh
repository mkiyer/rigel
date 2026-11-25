#!/usr/bin/env bash
# release.sh — Cut a new Rigel release.
#
# Usage:
#   ./scripts/publishing/release.sh 0.4.0
#
# What it does (stage 1 of the release pipeline):
#   1. Validates the version and repo state (no uncommitted changes to tracked files, on main, tag free).
#   2. Auto-updates CHANGELOG.md:
#        - replaces `## [Unreleased]` with `## [X.Y.Z] - YYYY-MM-DD` (if present)
#        - verifies a `## [X.Y.Z]` section exists
#        - appends the `[X.Y.Z]: ...compare/vPREV...vX.Y.Z` link at the bottom
#   3. Bumps the version string in pyproject.toml and conda/meta.yaml.
#   4. Shows a diff and asks for confirmation.
#   5. Commits ("Release vX.Y.Z"), tags `vX.Y.Z`, and pushes both to origin in one atomic push; --force
#      replaces an existing tag only after the confirmation, and PyPI must answer "not found" for the version.
#
# Pushing the tag triggers .github/workflows/publish.yml, which builds sdist +
# wheels and uploads to PyPI via OIDC; the tag push is the only trigger.
#
# After this script finishes, wait for the workflow to go green, then run:
#
#   ./scripts/publishing/post_release.sh 0.4.0
#
# which patches the conda recipe with the real PyPI sha256 and (optionally)
# uploads to your personal anaconda.org channel.

set -euo pipefail

RED=$'\033[0;31m'
GREEN=$'\033[0;32m'
YELLOW=$'\033[1;33m'
BLUE=$'\033[0;34m'
NC=$'\033[0m'

die()  { echo -e "${RED}Error:${NC} $1" >&2; exit 1; }
info() { echo -e "${GREEN}✓${NC} $1"; }
warn() { echo -e "${YELLOW}⚠${NC} $1"; }
note() { echo -e "${BLUE}→${NC} $1"; }

# Portable in-place sed (GNU vs BSD).
sed_inplace() {
    if sed --version >/dev/null 2>&1; then
        sed -i "$@"          # GNU
    else
        sed -i '' "$@"       # BSD / macOS
    fi
}

# ── Args ─────────────────────────────────────────────────────────────────
FORCE=0
POSITIONAL=()
for arg in "$@"; do
    case "$arg" in
        -f|--force) FORCE=1 ;;
        -h|--help) sed -n '2,/^$/p' "$0"; exit 0 ;;
        -*) die "Unknown flag: $arg" ;;
        *) POSITIONAL+=("$arg") ;;
    esac
done
[[ ${#POSITIONAL[@]} -eq 1 ]] || die "Usage: $0 <version> [--force]  (e.g. $0 0.4.0)"
VERSION="${POSITIONAL[0]}"
TAG="v${VERSION}"
TODAY=$(date +%Y-%m-%d)

[[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || \
    die "Version must be X.Y.Z (got: $VERSION)"

# ── Repo state ───────────────────────────────────────────────────────────
cd "$(git rev-parse --show-toplevel)" || die "Not in a git repository"

BRANCH=$(git branch --show-current)
[[ "$BRANCH" == "main" ]] || die "Must be on 'main' (currently on '$BRANCH')"

# The release commit must be exactly the tagged tree: an uncommitted change to a tracked file would be
# left out of it. Untracked files (a sandbox note) do not block a release.
git diff --quiet && git diff --cached --quiet || \
    die "Tracked files have uncommitted changes. Commit or stash them before releasing."

# PyPI refuses re-uploads, even of a deleted file, so a version PyPI knows is never reused. Only PyPI's own "not
# found" clears the version: any other answer, or none, stops the release.
PYPI_STATUS=$(curl -s -o /dev/null -w '%{http_code}' "https://pypi.org/pypi/rigel-rnaseq/${VERSION}/json") || PYPI_STATUS="no answer"
case "$PYPI_STATUS" in
    404) info "PyPI has no rigel-rnaseq==${VERSION}" ;;
    200) die "PyPI already has rigel-rnaseq==${VERSION}. Pick a new version (PyPI forbids re-uploads)." ;;
    *)   die "Could not ask PyPI whether rigel-rnaseq==${VERSION} exists (HTTP ${PYPI_STATUS}). Retry when it answers." ;;
esac

# A pre-existing tag. Nothing is deleted here: a --force release replaces the tag only after the confirmation below,
# and the remote one only if it still points where it pointed now.
TAG_EXISTS_LOCAL=0
git tag -l "$TAG" | grep -q . && TAG_EXISTS_LOCAL=1
REMOTE_REFS=$(git ls-remote origin "refs/tags/$TAG") || die "Could not read origin's tags"
REMOTE_TAG_SHA=$(awk -v ref="refs/tags/$TAG" '$2 == ref {print $1}' <<< "$REMOTE_REFS")

if [[ $TAG_EXISTS_LOCAL -eq 1 || -n "$REMOTE_TAG_SHA" ]]; then
    [[ $FORCE -eq 1 ]] || die "Tag $TAG already exists (local=$TAG_EXISTS_LOCAL remote=${REMOTE_TAG_SHA:-none}).
Re-run with --force to replace it (PyPI has no ${VERSION}, checked above)."
    warn "Tag $TAG exists (local=$TAG_EXISTS_LOCAL remote=${REMOTE_TAG_SHA:-none}); it will be replaced on confirmation (--force)"
fi

# Previous version = latest v-tag that isn't the one we're about to (re)create.
PREV_TAG=$(git tag -l 'v*' | grep -v "^${TAG}$" | sort -V | tail -1)
PREV_VERSION="${PREV_TAG#v}"

# Current version in pyproject.toml (source of truth).
OLD_VERSION=$(grep -E '^version = ' pyproject.toml | head -1 | sed -E 's/.*"([^"]+)".*/\1/')

info "Current:  ${OLD_VERSION}  (last tag: ${PREV_TAG:-none})"
info "New:      ${VERSION}  (tag: ${TAG})"
# Equal versions are OK under --force (re-release of an existing commit).
if [[ "$VERSION" == "$OLD_VERSION" && $FORCE -eq 0 ]]; then
    die "New version equals current version (use --force to re-release)"
fi

# ── CHANGELOG.md auto-update ─────────────────────────────────────────────
CHANGELOG="CHANGELOG.md"
[[ -f "$CHANGELOG" ]] || die "$CHANGELOG not found"

if grep -q "^## \[${VERSION}\]" "$CHANGELOG"; then
    info "CHANGELOG already has a [${VERSION}] section"
elif grep -q "^## \[Unreleased\]" "$CHANGELOG"; then
    sed_inplace "s|^## \[Unreleased\].*$|## [${VERSION}] - ${TODAY}|" "$CHANGELOG"
    info "CHANGELOG: [Unreleased] → [${VERSION}] - ${TODAY}"
else
    die "CHANGELOG.md has neither [${VERSION}] nor [Unreleased] heading.
Add a section like:

  ## [Unreleased]

  ### Added
  - ...

or

  ## [${VERSION}] - ${TODAY}

  ### Added
  - ..."
fi

# Append comparison link at the bottom if it isn't already there.
if ! grep -q "^\[${VERSION}\]: " "$CHANGELOG"; then
    if [[ -n "$PREV_TAG" ]]; then
        LINK="[${VERSION}]: https://github.com/mkiyer/rigel/compare/${PREV_TAG}...${TAG}"
    else
        LINK="[${VERSION}]: https://github.com/mkiyer/rigel/releases/tag/${TAG}"
    fi
    # Ensure file ends with a newline before appending.
    [[ -z "$(tail -c 1 "$CHANGELOG")" ]] || echo >> "$CHANGELOG"
    echo "$LINK" >> "$CHANGELOG"
    info "CHANGELOG: appended link $LINK"
fi

# ── Bump versions (idempotent) ───────────────────────────────────────────
if [[ "$OLD_VERSION" != "$VERSION" ]]; then
    sed_inplace "s|^version = \"${OLD_VERSION}\"|version = \"${VERSION}\"|" pyproject.toml
    info "Bumped pyproject.toml: ${OLD_VERSION} → ${VERSION}"
fi
grep -q "^version = \"${VERSION}\"" pyproject.toml || die "pyproject.toml version is not ${VERSION}"

META_OLD=$(grep -E '\{% set version' conda/meta.yaml | sed -E 's/.*"([^"]+)".*/\1/')
if [[ "$META_OLD" != "$VERSION" ]]; then
    sed_inplace "s|{% set version = \"${META_OLD}\" %}|{% set version = \"${VERSION}\" %}|" conda/meta.yaml
    info "Bumped conda/meta.yaml: ${META_OLD} → ${VERSION}"
fi
grep -q "{% set version = \"${VERSION}\" %}" conda/meta.yaml || die "conda/meta.yaml version is not ${VERSION}"

# ── Show & confirm ───────────────────────────────────────────────────────
echo ""
echo "Pending changes:"
echo "─────────────────"
git --no-pager diff --stat pyproject.toml conda/meta.yaml "$CHANGELOG"
echo ""
git --no-pager diff pyproject.toml conda/meta.yaml
echo ""

read -rp "Commit, tag ${TAG}, and push to origin? (y/N) " answer
if ! [[ "$answer" =~ ^[Yy]$ ]]; then
    warn "Aborted. Your working tree still has the version/changelog edits;"
    warn "revert with: git checkout -- pyproject.toml conda/meta.yaml ${CHANGELOG}"
    exit 0
fi

# ── Commit / tag / push ──────────────────────────────────────────────────
if git diff --cached --quiet pyproject.toml conda/meta.yaml "$CHANGELOG" && \
   git diff --quiet pyproject.toml conda/meta.yaml "$CHANGELOG"; then
    info "No version/changelog changes to commit (tagging HEAD)"
else
    git add pyproject.toml conda/meta.yaml "$CHANGELOG"
    git commit -m "Release ${TAG}"
    info "Committed"
fi

[[ $TAG_EXISTS_LOCAL -eq 1 ]] && git tag -d "$TAG" >/dev/null
git tag -a "$TAG" -m "Release ${TAG}"
info "Tagged ${TAG}"

# main and the tag land together or not at all. A replaced remote tag is overwritten only if it still points where it
# did when this script read it.
PUSH_ARGS=(--atomic)
[[ -n "$REMOTE_TAG_SHA" ]] && PUSH_ARGS+=("--force-with-lease=refs/tags/${TAG}:${REMOTE_TAG_SHA}")
git push "${PUSH_ARGS[@]}" origin main "refs/tags/${TAG}"
info "Pushed main + ${TAG} to origin"

# ── Done ─────────────────────────────────────────────────────────────────
cat <<EOF

${GREEN}Release ${TAG} initiated.${NC}

The tag push will trigger .github/workflows/publish.yml, which builds the
sdist and wheels and uploads them to PyPI via OIDC trusted publishing.

  Monitor the build:
    https://github.com/mkiyer/rigel/actions/workflows/publish.yml

Once the workflow is green (≈15–30 min), finish the release with:

  ${BLUE}./scripts/publishing/post_release.sh ${VERSION}${NC}

That script will:
  • poll PyPI for rigel-rnaseq==${VERSION}
  • patch conda/meta.yaml with the real sdist sha256
  • commit + push the recipe update
  • optionally build and upload to your personal anaconda.org channel

EOF
