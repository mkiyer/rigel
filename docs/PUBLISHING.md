# Publishing Rigel

Cutting a new release is now two commands, plus one wait.

> **Naming**
> - PyPI distribution: `rigel-rnaseq` (`pip install rigel-rnaseq`)
> - Bioconda package: `rigel` (`conda install -c bioconda rigel`)
> - Python import / CLI: `rigel`

---

## TL;DR

```bash
# 1. Add a `## [Unreleased]` section to CHANGELOG.md with your notes,
#    then cut the release:
./scripts/publishing/release.sh 0.4.0

# 2. Wait for the GitHub Actions publish workflow to go green
#    (~15–30 min): https://github.com/mkiyer/rigel/actions

# 3. Finalize: patch conda recipe with PyPI sha256 and (optionally)
#    upload to your personal anaconda.org channel:
./scripts/publishing/post_release.sh 0.4.0
```

That's it. No GitHub UI clicks, no manual sha256 lookups, no
`sed`-by-hand of YAML files.

---

## How it works

### `release.sh X.Y.Z` — stage 1 (fast)

1. Validates version format, clean tree, branch is `main`, tag free.
2. Auto-finishes the CHANGELOG:
   - renames `## [Unreleased]` → `## [X.Y.Z] - YYYY-MM-DD` (or verifies an
     explicit `## [X.Y.Z]` section already exists);
   - appends the `[X.Y.Z]: …compare/vPREV…vX.Y.Z` link at the bottom.
3. Bumps the version string in `pyproject.toml` **and** `conda/meta.yaml`.
4. Shows the diff and asks for confirmation.
5. `git commit`, `git tag vX.Y.Z`, `git push origin main --tags`.

**Pushing the tag is the trigger.** `.github/workflows/publish.yml` runs on
any `v*` tag push: it builds the sdist + Linux x86_64/aarch64 + macOS arm64
wheels and uploads to PyPI via OIDC trusted publishing. You no longer need
to open the GitHub Releases UI.

### `post_release.sh X.Y.Z` — stage 2 (run after the workflow is green)

1. Polls `https://pypi.org/pypi/rigel-rnaseq/X.Y.Z/json` every 30 s (up to
   30 min) until the sdist is live.
2. Extracts the sdist `sha256` from the PyPI API.
3. Patches `conda/meta.yaml` (portable Linux / macOS sed).
4. Commits `Release vX.Y.Z: update conda sha256` and pushes to `main`.
5. Optionally runs `conda_publish.sh` to build and upload to your personal
   anaconda.org channel (pass `--conda-publish` to skip the prompt, or
   `--skip-conda-publish` to suppress it entirely).

### `conda_publish.sh` — personal channel upload (optional)

Builds the conda recipe from the PyPI sdist and uploads it to your own
`anaconda.org/<user>` channel. Requires a `bioconda-build` env plus
`anaconda login`. Useful when you want the package installable right away
while bioconda is still reviewing.

```bash
conda activate bioconda-build
./scripts/publishing/conda_publish.sh
```

---

## Full walkthrough (first time)

### 1. Write changelog notes

Under the top of `CHANGELOG.md`, add an Unreleased section with your
release notes:

```markdown
## [Unreleased]

### Added
- ...

### Fixed
- ...
```

You do **not** need to fill in the version number or the date — stage 1
will do that for you. You also don't need to append the comparison link at
the bottom; stage 1 adds it.

### 2. Cut the release

```bash
./scripts/publishing/release.sh 0.4.0
```

The script prints the diff and asks for confirmation. If you answer "N",
your working tree keeps the staged edits and you can revert with
`git checkout -- pyproject.toml conda/meta.yaml CHANGELOG.md`.

### 3. Watch the PyPI build

<https://github.com/mkiyer/rigel/actions/workflows/publish.yml>

The `publish` job uses PyPI OIDC trusted publishing, so no API token is
needed.

Verify once green:

```bash
pip install rigel-rnaseq==0.4.0
rigel --version
```

### 4. Finalize

```bash
./scripts/publishing/post_release.sh 0.4.0
```

It will block until PyPI has the sdist, then patch/commit/push the conda
recipe and (on prompt) upload to your personal channel. The next bioconda
auto-bump bot run (within ~24h) picks up the new PyPI release and opens a
PR against `bioconda-recipes`.

If the auto-bump bot doesn't fire, `post_release.sh` prints the exact git
commands to open a manual PR.

---

## Recovery scenarios

### I answered "N" at the confirmation prompt

Nothing was committed, tagged, or pushed. Revert the staged edits and
start over:

```bash
git checkout -- pyproject.toml conda/meta.yaml CHANGELOG.md
./scripts/publishing/release.sh 0.4.0
```

### The tag pushed but `publish.yml` failed before uploading to PyPI

You can re-run the workflow against the same tag:

```bash
# From the Actions tab → publish.yml run → "Re-run failed jobs"
```

Or delete the tag and retry:

```bash
git tag -d v0.4.0
git push origin --delete v0.4.0
# fix the underlying problem, commit, then:
./scripts/publishing/release.sh 0.4.0
```

### PyPI upload succeeded but the release is broken

PyPI refuses re-uploads of the same version. Bump to a patch:

```bash
# Add a ## [Unreleased] section to CHANGELOG.md with the fix,
# then:
./scripts/publishing/release.sh 0.4.1
```

### `release.sh` says the tag already exists

You've already run stage 1. Either continue with
`./scripts/publishing/post_release.sh <same-version>`, or delete the tag
locally + remotely and retry (see above).

---

## One-time setup (already configured)

### PyPI trusted publisher

Registered at
<https://pypi.org/manage/project/rigel-rnaseq/settings/publishing/>:

| Field | Value |
|-------|-------|
| Owner | `mkiyer` |
| Repository | `rigel` |
| Workflow | `publish.yml` |
| Environment | `pypi` |

A matching GitHub environment named `pypi` exists in repo settings.

### Workflows

| Workflow | Trigger | What it does |
|----------|---------|--------------|
| `ci.yml` | Push / PR to `main` | Tests on Ubuntu + macOS, Python 3.12 + 3.13 |
| `publish.yml` | **Tag push `v*`** / GitHub Release / manual dispatch | Builds sdist + wheels → PyPI |

### Wheel matrix

| Platform | Arch | Runner / image |
|----------|------|----------------|
| Linux | x86_64 | `manylinux_2_28` |
| Linux | aarch64 | `manylinux_2_28` (QEMU) |
| macOS | arm64 | `macos-latest` |

Intel Mac users install from the sdist. Linux glibc ≥ 2.28 required.

---

## Release checklist

- [ ] `CHANGELOG.md` has a `## [Unreleased]` section with real notes
- [ ] `./scripts/publishing/release.sh X.Y.Z`
- [ ] `publish.yml` green on <https://github.com/mkiyer/rigel/actions>
- [ ] `pip install rigel-rnaseq==X.Y.Z` works in a clean env
- [ ] `./scripts/publishing/post_release.sh X.Y.Z`
- [ ] Bioconda PR opened (manually, or by the auto-bump bot within 24h)
