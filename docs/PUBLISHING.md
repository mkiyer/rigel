# Publishing Rigel to PyPI and Bioconda

This guide covers the complete process for publishing Rigel releases.

> **Note on naming:** The GitHub repository, Python import name, Bioconda package,
> and CLI command are all **`rigel`**. Only the **PyPI distribution name** is
> **`rigel-rnaseq`** (because `rigel` was already taken on PyPI).
>
> - Install from PyPI: `pip install rigel-rnaseq`
> - Install from Bioconda: `conda install -c bioconda rigel`
> - Import in Python: `import rigel`
> - Run CLI: `rigel --version`

---

## Prerequisites

1. **GitHub repository** at `https://github.com/mkiyer/rigel` (public)
2. **PyPI account** at [pypi.org](https://pypi.org/account/register/)
3. **GitHub account** with fork of [bioconda-recipes](https://github.com/bioconda/bioconda-recipes)

---

## Part 1: Publishing to PyPI

### One-Time Setup (Trusted Publisher)

PyPI supports "trusted publishing" — no API tokens needed. GitHub Actions
authenticates directly via OpenID Connect.

1. Go to [pypi.org/manage/account/publishing](https://pypi.org/manage/account/publishing/)
2. Add a new **pending publisher**:
   - **PyPI project name:** `rigel-rnaseq`
   - **Owner:** `mkiyer`
   - **Repository:** `rigel`
   - **Workflow name:** `publish.yml`
   - **Environment name:** `pypi`
3. In your GitHub repo, go to **Settings → Environments** and create an
   environment called `pypi`. Optionally add a deployment protection rule
   requiring manual approval.

### Release Workflow

```bash
# 1. Update version in pyproject.toml
#    version = "0.1.0" → "0.2.0"

# 2. Update CHANGELOG.md

# 3. Commit and tag
git add pyproject.toml CHANGELOG.md
git commit -m "Release v0.2.0"
git tag v0.2.0
git push origin main --tags

# 4. Create a GitHub Release from the tag
#    Go to: https://github.com/mkiyer/rigel/releases/new
#    Select tag: v0.2.0
#    Title: "Rigel v0.2.0"
#    Paste changelog notes
#    Click "Publish release"
```

Publishing the GitHub Release triggers the `publish.yml` workflow which:
1. Builds an sdist (source distribution)
2. Builds binary wheels for Linux (x86_64, aarch64) and macOS (x86_64, arm64)
3. Uploads everything to PyPI via trusted publishing

### Manual Testing Before First Release

```bash
# Build locally and inspect
pip install build
python -m build

# Test the sdist builds correctly
pip install dist/rigel_rnaseq-0.1.0.tar.gz

# Test the wheel installs correctly
pip install dist/rigel_rnaseq-0.1.0-cp312-abi3-*.whl
rigel --help
```

### Verifying the PyPI Upload

After the workflow completes:
- Check [pypi.org/project/rigel-rnaseq](https://pypi.org/project/rigel-rnaseq/)
- Verify `pip install rigel-rnaseq` works in a clean environment

---

## Part 2: Publishing to Bioconda

Bioconda requires that the package is **already on PyPI** (or has a stable
source tarball URL). The typical flow is: PyPI first, then Bioconda.

### Step 1: Get the sdist SHA256

After the PyPI upload, get the hash of the source tarball:

```bash
pip download --no-binary :all: --no-deps rigel-rnaseq==0.1.0 -d /tmp/rigel-dl
sha256sum /tmp/rigel-dl/rigel_rnaseq-0.1.0.tar.gz
```

Or from the PyPI JSON API:
```bash
curl -s https://pypi.org/pypi/rigel-rnaseq/0.1.0/json | python -c "
import sys, json
data = json.load(sys.stdin)
for f in data['urls']:
    if f['packagetype'] == 'sdist':
        print(f['digests']['sha256'])
"
```

### Step 2: Fork and Clone bioconda-recipes

```bash
# One-time setup
git clone https://github.com/mkiyer/bioconda-recipes.git
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git
```

### Step 3: Create the Recipe Branch

```bash
git fetch upstream
git checkout -b add-rigel upstream/master

# Create the recipe directory
mkdir -p recipes/rigel
```

### Step 4: Write the Recipe

Copy the template from `conda/meta.yaml` in this repo to
`recipes/rigel/meta.yaml` in the bioconda-recipes fork:

```bash
cp /path/to/rigel/conda/meta.yaml recipes/rigel/meta.yaml
```

**Update the sha256 hash** with the actual value from Step 1:

```yaml
source:
  url: https://pypi.org/packages/source/r/rigel-rnaseq/rigel_rnaseq-0.1.0.tar.gz
  sha256: ACTUAL_HASH_HERE
```

### Step 5: Test Locally with bioconda-utils

```bash
# Install bioconda-utils
conda install -c conda-forge -c bioconda bioconda-utils

# Lint the recipe
bioconda-utils lint --packages rigel

# Build locally (requires Docker on Linux)
bioconda-utils build --packages rigel --docker
```

### Step 6: Submit the Pull Request

```bash
git add recipes/rigel/meta.yaml
git commit -m "Add rigel 0.1.0"
git push origin add-rigel
```

Then open a PR to `bioconda/bioconda-recipes`. The PR template will appear
with a checklist. Key requirements:

- [x] Recipe follows [bioconda guidelines](https://bioconda.github.io/contributor/guidelines.html)
- [x] Package is already on PyPI (or has a stable URL)
- [x] License is open source (GPL-3.0-or-later ✓)
- [x] `test` section verifies the install (import + CLI)
- [x] `about.license_file` points to LICENSE
- [x] `extra.recipe-maintainers` lists your GitHub username

### Step 7: Review Process

Bioconda uses automated CI bots:
1. **BiocondaBot** checks recipe formatting and runs linting
2. **Azure Pipelines** builds the package on Linux and macOS
3. A bioconda team member reviews and merges

Common issues to watch for:
- Missing `run_exports` for C libraries
- Build failures from missing `pkg-config` or compiler packages
- Test failures if `rigel --help` has non-zero exit or unexpected output
- Hash mismatch if you updated the PyPI package after getting the hash

After merge, the package appears on bioconda within ~30 minutes:
```bash
conda install -c conda-forge -c bioconda rigel
```

### Updating Bioconda for New Releases

For subsequent releases, the bioconda bot can auto-bump:

```bash
cd bioconda-recipes
git fetch upstream
git checkout -b update-rigel-0.2.0 upstream/master

# Edit recipes/rigel/meta.yaml:
#   version: "0.2.0"
#   sha256: NEW_HASH
#   number: 0

git add recipes/rigel/meta.yaml
git commit -m "Update rigel to 0.2.0"
git push origin update-rigel-0.2.0
# Open PR
```

Alternatively, if you've set up the **autobump bot**, it will detect new
PyPI releases and open PRs automatically.

---

## Part 3: Checklist Summary

### Before First Release

- [ ] Verify `pyproject.toml` metadata (name, version, description, URLs, license)
- [ ] Ensure `MANIFEST.in` includes all C/C++ sources and headers
- [ ] Ensure `CHANGELOG.md` has release notes
- [ ] Test local build: `python -m build` succeeds
- [ ] Test local install: `pip install dist/rigel-*.tar.gz` works
- [ ] Set up PyPI trusted publisher (one-time)
- [ ] Create GitHub `pypi` environment (one-time)

### For Each Release

1. [ ] Bump version in `pyproject.toml`
2. [ ] Update `CHANGELOG.md`
3. [ ] Commit, tag, push
4. [ ] Create GitHub Release → triggers PyPI publish
5. [ ] Verify on PyPI: `pip install rigel-rnaseq==X.Y.Z`
6. [ ] Get sdist SHA256 from PyPI
7. [ ] Update bioconda recipe with new version + hash
8. [ ] Open PR to bioconda-recipes
9. [ ] Verify on Bioconda: `conda install -c bioconda rigel==X.Y.Z`

### Files Created/Modified

| File | Purpose |
|------|---------|
| `.github/workflows/ci.yml` | CI testing on push/PR |
| `.github/workflows/publish.yml` | Build wheels + publish to PyPI on release |
| `conda/meta.yaml` | Bioconda recipe template |
| `MANIFEST.in` | Ensures sdist includes C/C++ sources |
| `pyproject.toml` | Added project URLs, verified metadata |
| `CMakeLists.txt` | Added `RIGEL_PORTABLE` option for distributable wheels |
