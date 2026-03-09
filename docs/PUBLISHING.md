# Publishing Rigel

This guide covers the release workflow for PyPI and Bioconda.

> Naming summary:
>
> - PyPI distribution: `rigel-rnaseq`
> - Bioconda package: `rigel`
> - Python import: `rigel`
> - CLI command: `rigel`

---

## 1. Prerequisites

You need:

1. A public GitHub repository at `mkiyer/rigel`
2. A PyPI account
3. A fork of `bioconda/bioconda-recipes`

---

## 2. PyPI publishing

### 2.1 One-time trusted publisher setup

Rigel is configured for trusted publishing from GitHub Actions.

Create a pending publisher on PyPI with:

- project name: `rigel-rnaseq`
- owner: `mkiyer`
- repository: `rigel`
- workflow: `publish.yml`
- environment: `pypi`

Then create a GitHub environment named `pypi` in the repository settings.

### 2.2 Release preparation

Before a release:

1. Bump `version` in `pyproject.toml`
2. Update `CHANGELOG.md`
3. Ensure the docs are current
4. Verify that local builds still work

Local validation:

```bash
pip install build
python -m build

pip install dist/rigel_rnaseq-X.Y.Z.tar.gz
rigel --help

pip install dist/rigel_rnaseq-X.Y.Z-cp312-abi3-*.whl
rigel --version
```

### 2.3 Create the release

```bash
git add pyproject.toml CHANGELOG.md README.md docs/
git commit -m "Release vX.Y.Z"
git tag vX.Y.Z
git push origin main --tags
```

Create a GitHub Release from `vX.Y.Z`. Publishing the release triggers the
`publish.yml` workflow, which should:

1. build the sdist
2. build wheels for supported Linux and macOS targets
3. upload all artifacts to PyPI via trusted publishing

### 2.4 Post-release checks

Verify:

```bash
pip install rigel-rnaseq==X.Y.Z
rigel --version
```

Also check the project page:

- `https://pypi.org/project/rigel-rnaseq/`

---

## 3. Bioconda publishing

Bioconda should be done after PyPI, because the recipe source URL points to the
PyPI sdist.

### 3.1 Fetch the sdist hash

```bash
pip download --no-binary :all: --no-deps rigel-rnaseq==X.Y.Z -d /tmp/rigel-dl
shasum -a 256 /tmp/rigel-dl/rigel_rnaseq-X.Y.Z.tar.gz
```

Or via the PyPI JSON API:

```bash
curl -s https://pypi.org/pypi/rigel-rnaseq/X.Y.Z/json | python -c '
import json, sys
data = json.load(sys.stdin)
for item in data["urls"]:
      if item["packagetype"] == "sdist":
            print(item["digests"]["sha256"])
'
```

### 3.2 Update the Bioconda recipe

Use the recipe template from this repository:

```bash
cp /path/to/rigel/conda/meta.yaml recipes/rigel/meta.yaml
```

Update at least:

- version
- source URL
- source SHA256
- build number if needed

The source stanza should look like:

```yaml
source:
   url: https://pypi.org/packages/source/r/rigel-rnaseq/rigel_rnaseq-X.Y.Z.tar.gz
   sha256: ACTUAL_SHA256_HERE
```

### 3.3 Validate locally

```bash
conda install -c conda-forge -c bioconda bioconda-utils
bioconda-utils lint --packages rigel
bioconda-utils build --packages rigel --docker
```

### 3.4 Submit the PR

```bash
git fetch upstream
git checkout -b update-rigel-X.Y.Z upstream/master
git add recipes/rigel/meta.yaml
git commit -m "Update rigel to X.Y.Z"
git push origin update-rigel-X.Y.Z
```

Open a pull request against `bioconda/bioconda-recipes`.

After merge, verify:

```bash
conda install -c conda-forge -c bioconda rigel==X.Y.Z
rigel --version
```

---

## 4. Release checklist

### Before tagging

- `pyproject.toml` version bumped
- `CHANGELOG.md` updated
- `README.md`, `docs/MANUAL.md`, `docs/METHODS.md`, and `docs/PUBLISHING.md` reviewed
- `MANIFEST.in` still includes native sources and docs needed for the sdist
- `python -m build` succeeds locally

### After the GitHub Release

- PyPI workflow completed successfully
- `pip install rigel-rnaseq==X.Y.Z` works in a clean environment
- Bioconda recipe updated with the new hash
- Bioconda PR opened or autobump verified

---

## 5. Common failure modes

- PyPI trusted publisher mismatch: verify project name, repo name, workflow name, and environment
- Source build failure: verify `MANIFEST.in` still includes native source files
- Wheel portability issue: verify CI builds use the portable wheel configuration
- Bioconda hash mismatch: regenerate the SHA256 from the final PyPI sdist
- CLI smoke-test failure: ensure `rigel --help` and `rigel --version` both work in the built artifact
