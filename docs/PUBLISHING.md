# Publishing Rigel

Step-by-step release process for PyPI and Bioconda.

> **Naming:**
> - PyPI distribution: `rigel-rnaseq` (`pip install rigel-rnaseq`)
> - Bioconda package: `rigel` (`conda install -c bioconda rigel`)
> - Python import / CLI: `rigel`

---

## Quick reference

```bash
# 1. Edit CHANGELOG.md with the new version's notes
# 2. Run the release script:
./scripts/release.sh 0.3.0
# 3. Create a GitHub Release from the tag (triggers PyPI publish)
# 4. After PyPI succeeds, fetch the hash for bioconda:
./scripts/bioconda_hash.sh 0.3.0
```

---

## Full walkthrough

### Step 1: Update the changelog

Add an entry to `CHANGELOG.md` for the new version **before** running the
release script. The script checks for this entry and will refuse to proceed
without it.

```markdown
## [0.3.0] - 2026-04-15

### Added
- ...

### Fixed
- ...
```

Also add the comparison link at the bottom of the file:

```markdown
[0.3.0]: https://github.com/mkiyer/rigel/compare/v0.2.0...v0.3.0
```

### Step 2: Run the release script

```bash
./scripts/release.sh 0.3.0
```

This script:
1. Validates that `CHANGELOG.md` has the version entry
2. Bumps the version in `pyproject.toml` and `conda/meta.yaml`
3. Shows a diff for confirmation
4. Commits, tags `v0.3.0`, and pushes to `origin main`

### Step 3: Create a GitHub Release

Go to the link printed by the script (or navigate to **Releases → Draft a new
release** on GitHub). Select the tag, paste the changelog entry as the body,
and click **Publish release**.

Publishing the release triggers `.github/workflows/publish.yml`, which:
1. Builds the source distribution (sdist)
2. Builds binary wheels for Linux x86_64, Linux aarch64, and macOS arm64
3. Uploads everything to PyPI via OIDC trusted publishing (no API token needed)

Monitor the workflow at: `https://github.com/mkiyer/rigel/actions`

### Step 4: Verify the PyPI release

```bash
pip install rigel-rnaseq==0.3.0
rigel --version
```

Check: https://pypi.org/project/rigel-rnaseq/

### Step 5: Publish to your personal conda channel (fast path)

This uploads to your own `anaconda.org/mkiyer` channel — no review required.
Useful for immediate availability while waiting for bioconda.

```bash
# One-time login (saved in ~/.anaconda/config.yaml):
conda activate bioconda-build
anaconda login

# Build and upload:
./scripts/conda_publish.sh
```

The script builds the package from `conda/meta.yaml` (using the PyPI sdist as
source), then uploads to your personal channel. Install with:

```bash
conda install -c mkiyer -c conda-forge -c bioconda rigel==0.3.0
```

### Step 6: Update bioconda

After the PyPI upload succeeds, fetch the sdist SHA256:

```bash
./scripts/bioconda_hash.sh 0.3.0
```

This patches `conda/meta.yaml` with the real hash from PyPI.

**First submission** (one-time setup):

```bash
# Fork https://github.com/bioconda/bioconda-recipes on GitHub, then:
git clone https://github.com/mkiyer/bioconda-recipes.git
cd bioconda-recipes
git checkout -b add-rigel
mkdir -p recipes/rigel
cp /path/to/rigel/conda/meta.yaml recipes/rigel/meta.yaml
git add recipes/rigel/meta.yaml
git commit -m "Add rigel version"
git push origin add-rigel
# Open a PR at https://github.com/bioconda/bioconda-recipes/pulls
```

**Subsequent releases**: Bioconda's auto-bump bot watches PyPI and usually
opens a PR within 24 hours. If it doesn't, manually update
`recipes/rigel/meta.yaml` with the new version and hash.

### Step 7: Verify bioconda

After the bioconda PR merges (allow ~24h for the package to build):

```bash
conda install -c conda-forge -c bioconda rigel==0.3.0
rigel --version
```

---

## Recovery: re-releasing the same version

If CI fails after `release.sh` ran (e.g., a missing dependency), you need to
fix the code and re-tag. The exact steps depend on how far the release got.

### Case A: Tag was pushed but no GitHub Release yet

The tag exists on GitHub but publish.yml hasn't run (no Release created).

```bash
# 1. Fix the bug on main
git add <files>
git commit -m "Fix: description of fix"

# 2. Delete the old tag (local + remote)
git tag -d v0.3.0
git push origin --delete v0.3.0

# 3. Re-run the release script (it will re-tag the new HEAD)
./scripts/release.sh 0.3.0
```

### Case B: GitHub Release was created and publish.yml failed

The Release exists on GitHub but the PyPI upload didn't happen (or partially
failed).

```bash
# 1. Delete the GitHub Release
#    Go to https://github.com/mkiyer/rigel/releases
#    Click on the release → Delete release

# 2. Fix the bug on main
git add <files>
git commit -m "Fix: description of fix"

# 3. Delete the old tag (local + remote)
git tag -d v0.3.0
git push origin --delete v0.3.0

# 4. Re-run the release script
./scripts/release.sh 0.3.0

# 5. Create a new GitHub Release from the new tag
```

### Case C: PyPI upload succeeded but the release is broken

PyPI does **not** allow re-uploading the same version. You must bump to a
patch release.

```bash
# 1. Fix the bug on main
# 2. Update CHANGELOG.md with a 0.3.1 entry
# 3. Run release for the new version:
./scripts/release.sh 0.3.1
```

### Case D: release.sh failed before creating the tag

The script died mid-way (e.g., CHANGELOG check failed, or you answered "N"
at the confirmation prompt). No tag was created.

```bash
# 1. Fix whatever caused the failure
# 2. Commit the fix
# 3. Re-run — the script is idempotent when no tag exists:
./scripts/release.sh 0.3.0
```

---

## One-time setup (already done)

### PyPI trusted publisher

Configured at https://pypi.org/manage/project/rigel-rnaseq/settings/publishing/:

| Field | Value |
|-------|-------|
| Owner | `mkiyer` |
| Repository | `rigel` |
| Workflow | `publish.yml` |
| Environment | `pypi` |

A matching GitHub environment named `pypi` exists in the repository settings.

### CI/CD workflows

| Workflow | Trigger | What it does |
|----------|---------|--------------|
| `ci.yml` | Push/PR to `main` | Tests on Ubuntu + macOS, Python 3.12 + 3.13 |
| `publish.yml` | GitHub Release or manual dispatch | Builds sdist + wheels, uploads to PyPI |

### Wheel targets

| Platform | Architecture | Image / Runner |
|----------|-------------|----------------|
| Linux | x86_64 | `manylinux_2_28` (AlmaLinux 8) |
| Linux | aarch64 | `manylinux_2_28` (QEMU emulated) |
| macOS | arm64 | `macos-latest` (macOS 15) |

Intel Mac users install from the sdist (`pip install rigel-rnaseq` compiles
locally). Linux glibc >= 2.28 is required for the binary wheel (RHEL 8+,
Ubuntu 20.04+, Debian 10+).

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| Trusted publisher mismatch | Verify project name, repo, workflow name, and environment match exactly |
| Wheel build fails on Linux | Check htslib compilation in `CIBW_BEFORE_ALL_LINUX` |
| pyarrow install fails in test | Ensure `CIBW_MANYLINUX_*_IMAGE` is `manylinux_2_28` (pyarrow dropped manylinux2014) |
| Missing dependency in CI | Add it to `pyproject.toml` dependencies, `conda/meta.yaml` run deps, and `publish.yml` `CIBW_TEST_REQUIRES` |
| Bioconda hash mismatch | Re-run `./scripts/bioconda_hash.sh` after PyPI upload completes |
| sdist missing native sources | Check `MANIFEST.in` includes `src/rigel/native/*.cpp` and `*.h` |
| `release.sh` says tag exists | See "Recovery: re-releasing the same version" above |

---

## Release checklist

- [ ] `CHANGELOG.md` updated with new version entry
- [ ] `./scripts/release.sh X.Y.Z` — bumps versions, commits, tags, pushes
- [ ] GitHub Release created from tag -> `publish.yml` runs
- [ ] `publish.yml` succeeds (all wheel + sdist jobs green)
- [ ] `pip install rigel-rnaseq==X.Y.Z` works in a clean environment
- [ ] `./scripts/bioconda_hash.sh X.Y.Z` — patches conda recipe with SHA256
- [ ] `./scripts/conda_publish.sh` — builds and uploads to personal anaconda.org channel
- [ ] Bioconda recipe PR opened (first time) or auto-bump confirmed
