# Installing STAR Aligner on macOS ARM (Apple Silicon M1/M2/M3)

## Background

The bioconda STAR package does not work on Apple Silicon (ARM) because
it ships x86-64 binaries. You must compile from source. The good news:
STAR now uses **SIMDe** (SIMD Everywhere), which transparently translates
x86 SSE/AVX intrinsics to ARM NEON, so native ARM compilation is possible.

> **Note:** STAR's README lists "x86-64 compatible processors" as a hardware
> requirement, but this is outdated — SIMDe support was added specifically
> to enable non-x86 architectures.

---

## Prerequisites

### 1. Install Homebrew (if not installed)

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### 2. Install GCC via Homebrew

Apple's Clang does not support OpenMP, which STAR requires.
You need **real GCC** from Homebrew:

```bash
brew install gcc
```

This installs `g++-14` (or similar version) under Homebrew's prefix.
Verify:

```bash
# Find the actual g++ binary
ls $(brew --prefix gcc)/bin/g++-*
# Example output: /opt/homebrew/Cellar/gcc/14.2.0/bin/g++-14
```

Note the full path — you will need it below. On Apple Silicon Macs,
Homebrew installs to `/opt/homebrew/` (not `/usr/local/`).

---

## Compile STAR from Source

### 1. Get STAR source code

```bash
# Option A: Download release tarball
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xzf 2.7.11b.tar.gz
cd STAR-2.7.11b

# Option B: Clone the repository
git clone https://github.com/alexdobin/STAR.git
cd STAR
```

### 2. Build for macOS ARM

The default Makefile sets `-mavx2` (an x86 SIMD flag) which will fail
on ARM. Override it with an empty `CXXFLAGS_SIMD` and point to real gcc:

```bash
cd source

# Find your gcc version
GCC_BIN=$(brew --prefix gcc)/bin/g++-14  # adjust version number if needed

# Clean any previous build artifacts
make clean

# Build for macOS with ARM-native optimization
make STARforMacStatic \
    CXX=${GCC_BIN} \
    CXXFLAGS_SIMD="" \
    CXXFLAGSextra="-mcpu=native"
```

**What the flags do:**
- `STARforMacStatic` — macOS build target with correct linker flags
- `CXX=${GCC_BIN}` — use Homebrew gcc instead of Apple Clang
- `CXXFLAGS_SIMD=""` — disable x86-specific `-mavx2` flag
  (SIMDe handles SIMD translation automatically)
- `CXXFLAGSextra="-mcpu=native"` — optimize for your specific ARM chip
  (note: ARM uses `-mcpu` not `-march`)

### 3. Install the binary

```bash
# Copy to a location on your PATH
sudo cp STAR /usr/local/bin/
# Or for user-local install:
cp STAR ~/bin/  # make sure ~/bin is on your PATH
```

### 4. Verify installation

```bash
STAR --version
# Should print: 2.7.11b
```

---

## Troubleshooting

### "error: unrecognized command-line option '-mavx2'"

You forgot to set `CXXFLAGS_SIMD=""`. ARM gcc does not recognize
x86 SIMD flags. Re-run make with the flag override shown above.

### "fatal error: 'omp.h' file not found" or OpenMP errors

You are using Apple Clang instead of Homebrew gcc. Verify:

```bash
${GCC_BIN} --version
# Should say "g++ (Homebrew GCC ...)" NOT "Apple clang"
```

### htslib compilation errors

If htslib fails to compile, clean and rebuild:

```bash
make CLEAN
make STARforMacStatic CXX=${GCC_BIN} CXXFLAGS_SIMD="" CXXFLAGSextra="-mcpu=native"
```

### Shared memory errors at runtime

STAR uses shared memory for genome loading. On macOS, you may need to
increase shared memory limits:

```bash
sudo sysctl -w kern.sysv.shmmax=36000000000
sudo sysctl -w kern.sysv.shmall=36000000000
```

Or use `--genomeLoad NoSharedMemory` (the default) to avoid this entirely.

### Link-time errors: "ld: library not found for -lz"

Install zlib:

```bash
brew install zlib
export LDFLAGS="-L$(brew --prefix zlib)/lib"
export CPPFLAGS="-I$(brew --prefix zlib)/include"
```

Then re-run the `make` command.

---

## Alternative: Rosetta 2 (fallback)

If native ARM compilation fails for any reason, you can run the x86-64
pre-compiled binary under Rosetta 2:

```bash
# 1. Install Rosetta 2 (if not already installed)
softwareupdate --install-rosetta --agree-to-license

# 2. Download the pre-compiled Linux binary (won't work)
#    OR use the x86 Homebrew to install gcc and compile:

# Start an x86 shell
arch -x86_64 /bin/zsh

# Install x86 Homebrew (installs to /usr/local/)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install x86 gcc
/usr/local/bin/brew install gcc

# Build STAR using x86 gcc
cd STAR/source
make clean
make STARforMacStatic CXX=/usr/local/bin/g++-14

# The resulting binary will run under Rosetta 2 translation
sudo cp STAR /usr/local/bin/STAR
```

> **Performance note:** Rosetta 2 adds ~20-30% overhead. Native ARM
> compilation (the method above) is strongly preferred.

---

## Quick Reference: STAR Alignment Command

Once installed, a typical paired-end alignment with STAR:

```bash
# Step 1: Generate genome index (once)
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles /path/to/genome.fa \
     --sjdbGTFfile /path/to/genes.gtf \
     --sjdbOverhang 100 \
     --runThreadN 8

# Step 2: Align reads
STAR --runMode alignReads \
     --genomeDir /path/to/star_index \
     --readFilesIn reads_R1.fastq.gz reads_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS NM MD XS \
     --outMultimapperOrder Random \
     --runThreadN 8

# Key flags for hulkrna compatibility:
#   NH — number of reported alignments (used for multimap detection)
#   HI — hit index (groups mate pairs across multimappers)
#   XS — strand for splice junctions
```

---

## Memory Requirements

STAR genome indexing and alignment require substantial RAM:

| Genome       | Index Generation | Alignment |
|-------------|-----------------|-----------|
| Human (hg38) | ~32 GB         | ~32 GB    |
| Mouse (mm10) | ~32 GB         | ~32 GB    |

The M3 MacBook Pro with 36 GB or more unified memory should work.
For 16-18 GB machines, reduce `--genomeSAsparseD 2` during index
generation to trade speed for memory (~16 GB).
