# Publishing Rigel

This is the release procedure: how a version of Rigel reaches PyPI and Bioconda, the two
scripts that do it, what the GitHub Actions workflow builds, and how to recover when a step
fails. It does not describe how to build or test the code (`README.md` and `CLAUDE.md`) or
what a release contains (the changelog and git). Cutting a release is two commands plus one wait.

> **Naming**
> - PyPI distribution: `rigel-rnaseq` (`pip install rigel-rnaseq`)
> - Bioconda package: `rigel` (`conda install -c bioconda rigel`)
> - Python import / CLI: `rigel`

---

## TL;DR

```bash
# 1. Add a `## [Unreleased]` section to CHANGELOG.md with your notes,
#    then cut the release:
./scripts/publishing/release.sh X.Y.Z

# 2. Wait for the GitHub Actions publish workflow to go green
#    (~15–30 min): https://github.com/mkiyer/rigel/actions

# 3. Finalize: patch conda recipe with PyPI sha256 and (optionally)
#    upload to your personal anaconda.org channel:
./scripts/publishing/post_release.sh X.Y.Z
```


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
any `v*` tag push (and on manual dispatch, never on a GitHub Release — that
would double-fire and the second upload would fail as "file already exists"):
it builds the sdist + Linux x86_64/aarch64 + macOS arm64 wheels and uploads
to PyPI via OIDC trusted publishing.

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

Builds the conda recipe from the PyPI sdist and uploads it to your own `anaconda.org/<user>`
channel (`conda activate bioconda-build && ./scripts/publishing/conda_publish.sh`; needs
`anaconda login`). Useful while bioconda is still reviewing.

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

Stage 1 fills in the version, the date and the comparison link.

### 2. Cut the release

```bash
./scripts/publishing/release.sh X.Y.Z
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
pip install rigel-rnaseq==X.Y.Z
rigel --version
```

### 4. Finalize

```bash
./scripts/publishing/post_release.sh X.Y.Z
```

It blocks until PyPI has the sdist, then patches, commits and pushes the conda recipe and (on
prompt) uploads to your personal channel. The bioconda auto-bump bot (within ~24h) opens a PR
against `bioconda-recipes`; if it does not fire, `post_release.sh` prints the git commands for a
manual PR.

---

## Recovery scenarios

- **Answered "N" at the confirmation prompt.** Nothing was committed, tagged or pushed:
  `git checkout -- pyproject.toml conda/meta.yaml CHANGELOG.md`, then run `release.sh` again.
- **The tag pushed but `publish.yml` failed before uploading.** Re-run the failed jobs from the
  Actions tab, or delete the tag and retry after fixing the problem:

  ```bash
  git tag -d vX.Y.Z && git push origin --delete vX.Y.Z
  ./scripts/publishing/release.sh X.Y.Z
  ```

- **PyPI upload succeeded but the release is broken.** PyPI refuses re-uploads of a version:
  add a `## [Unreleased]` section with the fix and release the next patch version.
- **`release.sh` says the tag already exists.** Stage 1 already ran: continue with
  `./scripts/publishing/post_release.sh <same-version>`, or delete the tag as above and retry.

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
| `publish.yml` | **Tag push `v*`** / manual dispatch (deliberately not GitHub Release, which would double-fire) | Builds sdist + wheels → PyPI |

### Wheel matrix

| Platform | Arch | Runner / image |
|----------|------|----------------|
| Linux | x86_64 | `manylinux_2_28` |
| Linux | aarch64 | `manylinux_2_28` (QEMU) |
| macOS | arm64 | `macos-latest` |

`CIBW_BUILD: "cp312-*"` (with `*-musllinux_*` skipped) builds CPython 3.12 wheels only, so
Python 3.13 — which CI tests and the package classifies — installs from the sdist on every
platform, as do Intel Macs. Linux glibc ≥ 2.28 required.

---

## Release checklist

- [ ] `CHANGELOG.md` has a `## [Unreleased]` section with real notes
- [ ] `./scripts/publishing/release.sh X.Y.Z`
- [ ] `publish.yml` green on <https://github.com/mkiyer/rigel/actions>
- [ ] `pip install rigel-rnaseq==X.Y.Z` works in a clean env
- [ ] `./scripts/publishing/post_release.sh X.Y.Z`
- [ ] Bioconda PR opened (manually, or by the auto-bump bot within 24h)
