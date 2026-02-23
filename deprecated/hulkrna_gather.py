import re
import os
import sys
import argparse
import glob
import collections
import logging
import json
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

# import local modules
SCRIPT_PATH = Path(__file__).resolve().parent
HULKRNA_LIB_PATH = SCRIPT_PATH / '..' / 'lib'
if HULKRNA_LIB_PATH.exists() and str(HULKRNA_LIB_PATH) not in sys.path:
    sys.path.insert(0, str(HULKRNA_LIB_PATH))

# import from local modules
from hulkrna.txref import TxIndex
from hulkrna.core import CountType


def _series_find_max_length(series):
    """
    Returns the max string length in a pandas Series (excluding NAs).
    """
    # If the column has non-string objects, cast to string and get the length
    return series.dropna().astype(str).map(len).max()

def _series_to_fixed_width_bytes(series):
    """
    Converts a pandas Series of object/string dtype to a fixed-width bytes ndarray.
    Pads/truncates as needed. Skips NA handling for brevity.
    """
    maxlen = _series_find_max_length(series)
    dtype = f"S{maxlen}"
    # Convert: first to str, then encode
    return series.fillna('').astype(str).str.encode('utf-8').to_numpy(dtype=dtype)


class HulkData:
    COUNT_COMPRESSION = 'gzip'
    COUNT_DTYPE = 'f8'
    NDIM = 3

    def __init__(self):
        self._h5f = None

    def close(self):
        if self._h5f is not None:
            self._h5f.close()
            self._h5f = None
 
    @staticmethod
    def create(filename, nrows, ncols):
        f = h5py.File(filename, 'w')
        # count dataset
        f.create_dataset('counts', 
                         shape=(nrows, ncols, len(CountType)), 
                         dtype=HulkData.COUNT_DTYPE, 
                         compression=HulkData.COUNT_COMPRESSION)
        # empty metadata
        for i in range(HulkData.NDIM):
            f.create_group(f'meta/axis{i}', track_order=True)
        self = HulkData()
        self._h5f = f
        # add the CountType enum to metadata
        df = pd.DataFrame(CountType.columns(), columns=['count_type'])
        self.add_metadata(df, axis=2)
        return self

    def add_metadata(self, df, axis):
        # ensure metadata matches axis length
        axis_len = self._h5f['counts'].shape[axis]
        assert len(df) == axis_len, f"Metadata length {len(df)} does not match axis {axis} length {axis_len}"
        
        group = self._h5f[f'meta/axis{axis}']
        for col_name, col_dtype in df.dtypes.items():
            # R does not support int64, so convert to int32
            if pd.api.types.is_integer_dtype(col_dtype):
                # convert all integer types to 32-bit
                col_data = df[col_name].to_numpy(dtype=np.int32)            
            if col_dtype == 'object':
                col_data = _series_to_fixed_width_bytes(df[col_name])
            else:
                col_data = df[col_name].to_numpy(dtype=col_dtype)
            group.create_dataset(col_name, data=col_data, dtype=col_data.dtype)

    def add_counts(self, i, df):
        # convert DataFrame to numpy array, ensuring that 
        # count columns in DataFrame exist in correct order
        arr = df[CountType.columns()].to_numpy(dtype=HulkData.COUNT_DTYPE)
        # write to the hdf5 dataset
        dset = self._h5f['counts']
        dset[:, i, :] = arr


# output file constants
FASTP_JSON_FILE = 'fastp.json'
K2_RRNA_REPORT_FILE = 'k2_rrna.report.tsv.gz'
SAMTOOLS_MARKDUP_STATS_FILE = 'samtools_markdup_stats.json'
FIXUMI_STATS_FILE = 'fixumi_stats.json'
HULKRNA_COUNT_SUMMARY_FILE = 'summary.json'
STAR_LOG_FILE = 'Log.final.out'
HULKRNA_COUNT_GENE_FILE = 'gene_counts.feather'

# output file names
OUTPUT_H5_FILE = 'counts.h5'



def _flatten_dict(d, parent_key='', sep='_'):
    '''
    flatten nested dictionary to a single level
    '''
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(_flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def parse_hulkrna_count_summary(lib_dir):
    # process json file
    filename = lib_dir / 'hulkrna' / 'summary.json'
    with open(filename) as f:
        data = json.load(f)
    flat_data = _flatten_dict(data)
    d = {}
    for k, v in flat_data.items():
        k = f'hulkrna_{k}'
        d[k] = v
    return d


def parse_hulkrna_bam_stats(lib_dir):
    # process json file
    filename = lib_dir / 'hulkrna' / 'bam_stats.json'
    with open(filename) as f:
        data = json.load(f)
    flat_data = _flatten_dict(data)
    d = {}
    for k, v in flat_data.items():
        k = f'bam_{k}'
        d[k] = v
    return d


def _clean_star_log_key(key):
    key = key.strip().lower()
    # Replace "number" as a word with "num"
    key = re.sub(r'\bnumber\b', "num", key)
    # Replace "%" with "pct"
    key = key.replace('%', 'pct')
    # Remove "of" as a word (not part of another word)
    key = re.sub(r'\bof\b', '', key)
    # Replace one or more whitespace chars with a single underscore
    key = re.sub(r'\s+', '_', key)
    # Remove any character that's not a-z, 0-9, or underscore
    key = re.sub(r'[^a-z0-9_]', '', key)
    # Collapse multiple underscores (which may appear after cleaning)
    key = re.sub(r'_+', '_', key)
    # Strip leading/trailing underscores
    key = key.strip('_')
    return key

def parse_star_log(lib_dir):
    filename = lib_dir / 'star' / STAR_LOG_FILE
    data = {}
    with open(filename) as fh:
        for line in fh:
            if "|" not in line:
                continue
            key, value = line.split("|", 1)
            key = _clean_star_log_key(key)
            key = f'star_{key}'
            # clean up value
            value = value.strip().rstrip('%')
            data[key] = value
    return data


def parse_fixumi_stats(lib_dir):
    d = {
        'fixumi_num_reads': 0,
        'fixumi_umis_found': 0,
        'fixumi_frac_umis_found': 0.0
    }
    filename = lib_dir / FIXUMI_STATS_FILE
    if filename.exists():
        with open(filename, 'r') as f:
            data = json.load(f)
            d['fixumi_num_reads'] = int(data['num_reads'])
            d['fixumi_umis_found'] = int(data['num_umis_found'])
            d['fixumi_frac_umis_found'] = float(data['frac_umis_found'])
    return d


def parse_samtools_markdup_stats(lib_dir):
    filename = lib_dir / SAMTOOLS_MARKDUP_STATS_FILE
    d = {}
    with open(filename, 'r') as f:
        data = json.load(f) 
        d['markdup_total_reads'] = int(data['READ'])
        d['markdup_excluded'] = int(data['EXCLUDED'])
        d['markdup_duplicates'] = int(data['DUPLICATE TOTAL'])
        d['markdup_duplicate_pct'] = round(100.0 * float(d['markdup_duplicates']) / d['markdup_total_reads'], 2)
        d['markdup_est_lib_size'] = int(data['ESTIMATED_LIBRARY_SIZE'])
    return d


def parse_kraken2_rrna(lib_dir):
    '''
    Kraken 2's standard sample report format is tab-delimited with one line per taxon. 
    The fields of the output, from left-to-right, are as follows:
    - Percentage of fragments covered by the clade rooted at this taxon
    - Number of fragments covered by the clade rooted at this taxon
    - Number of fragments assigned directly to this taxon
    - A rank code, indicating (U)nclassified, (R)oot, (D)omain, 
      (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. 
      Taxa that are not at any of these 10 ranks have a rank code that is 
      formed by using the rank code of the closest ancestor rank with a 
      number indicating the distance from that rank. E.g., "G2" is a 
      rank code indicating a taxon is between genus and species and the 
      grandparent taxon is at the genus rank. 
    - taxonomic ID number
    - indented scientific name
    '''
    filename = lib_dir / 'kraken2' / K2_RRNA_REPORT_FILE
    # header_fields = ['pct', 'frags', 'frags_taxon', 'rank', 'taxid', 'taxname']
    # field_pct = header_fields.index('pct')
    # field_frags = header_fields.index('frags')
    # field_taxid = header_fields.index('taxid')
    K2_REPORT_PCT_FIELD = 0
    K2_REPORT_FRAGS_FIELD = 1
    K2_REPORT_TAXID_FIELD = 4
    UNCLASSIFIED_TAXID = 0
    # init dictionary with zeros
    d = {
        'k2_rrna_pct_unclassified': 0.0,
        'k2_rrna_frags_unclassified': 0
    }
    # read the file
    with gzip.open(filename, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            taxid = int(fields[K2_REPORT_TAXID_FIELD].strip())
            if taxid == UNCLASSIFIED_TAXID:
                d['k2_rrna_pct_unclassified'] = float(fields[K2_REPORT_PCT_FIELD])
                d['k2_rrna_frags_unclassified'] = int(fields[K2_REPORT_FRAGS_FIELD])
                break
    return d


def parse_fastp_metrics(lib_dir):
    filename = lib_dir / 'fastp' / FASTP_JSON_FILE
    d = {}
    with open(filename, 'r') as f:
        data = json.load(f) 
        d['fastp_total_reads_before_filtering'] = int(data['summary']['before_filtering']['total_reads'])
        d['fastp_total_megabases_before_filtering'] = float(data['summary']['before_filtering']['total_bases']) / 1e6
        d['fastp_total_reads_after_filtering'] = int(data['summary']['after_filtering']['total_reads'])
        d['fastp_total_megabases_after_filtering'] = float(data['summary']['after_filtering']['total_bases']) / 1e6
        d['fastp_passed_filter_reads'] = int(data['filtering_result']['passed_filter_reads'])
        d['fastp_insert_size_peak'] = int(data['frag_length']['peak'])
        d['fastp_adapter_trimmed_reads'] = int(data['adapter_cutting']['adapter_trimmed_reads'])
        d['fastp_adapter_trimmed_megabases'] = float(data['adapter_cutting']['adapter_trimmed_bases']) / 1e6
    return d


def parse_metadata(lib_dir):
    lib_dir = Path(lib_dir)
    d = {}
    d.update(parse_fastp_metrics(lib_dir))
    d.update(parse_kraken2_rrna(lib_dir))
    d.update(parse_samtools_markdup_stats(lib_dir))
    d.update(parse_fixumi_stats(lib_dir))
    d.update(parse_star_log(lib_dir))
    d.update(parse_hulkrna_bam_stats(lib_dir))
    d.update(parse_hulkrna_count_summary(lib_dir))
    return d


def _find_results(search_dirs=None, libs_file=None):
    '''
    Finds and aggregates library result directories based on the provided
    search directories and libraries file.
    '''
    # read libs file
    search_libs = set()
    if libs_file is not None:
        logging.info(f'Reading libraries file {libs_file}')
        with open(libs_file) as f:
            search_libs.update(x.strip() for x in f)
        logging.info(f'Found libraries: {len(search_libs)}')
    else:
        logging.info(f'No libraries file specified, using all libraries')

    # search directories
    log_msg = 'Searching directories'
    logging.info(f'[START] {log_msg}')
    lib_dirs = set()
    if search_dirs:
        for search_dir in search_dirs:
            logging.debug(f'Searching {search_dir}')
            n = 0
            total = 0
            for x in glob.glob(os.path.join(search_dir, '*')):
                if os.path.exists(os.path.join(x, 'done.txt')):
                    if len(search_libs) == 0 or (os.path.basename(x) in search_libs):
                        lib_dirs.add(x)
                        n += 1
                total += 1
            logging.info(f'Found {n} matches from {total} library results in {search_dir}')
    logging.info(f'Found a total of {len(lib_dirs)} library result dirs')
    logging.info(f'[DONE] {log_msg}')
    lib_dirs = sorted(lib_dirs)
    return lib_dirs


def _create_data_object(h5_file, g_meta_df, lib_meta_df):
    # create data object
    logging.info(f'Creating HulkData object at {h5_file}')
    hdata = HulkData.create(h5_file, len(g_meta_df), len(lib_meta_df))
    # add gene metadata
    logging.info(f'Adding gene metadata...')
    hdata.add_metadata(g_meta_df, axis=0)
    # add library metadata
    hdata.add_metadata(lib_meta_df, axis=1)
    return hdata


def hulkrna_gather(output_file, index_dir, search_dirs=None, libs_file=None):
    # search for library result directories
    lib_dirs = _find_results(search_dirs=search_dirs, libs_file=libs_file)
    libs = [os.path.basename(x) for x in lib_dirs]
    lib_df = pd.DataFrame(libs, columns=['library'])

    # load transcript metadata from index dir
    # logging.info(f'Loading transcript metadata from {index_dir}')
    # t_df = pd.read_feather(index_dir / 'transcripts.feather')

    # load gene metadata from index dir
    logging.info(f'Loading gene metadata from {index_dir}')
    g_df = pd.read_feather(index_dir / 'genes.feather')
    g_df.drop(columns=['t_id', 't_index'], inplace=True)

    # setup output directory
    hdata = _create_data_object(output_file, g_df, lib_df)

    # add libraries to data object
    logging.info(f'Parsing result directories')
    meta = collections.OrderedDict()
    for i, lib_dir in enumerate(lib_dirs):
        lib_dir = Path(lib_dir)
        lib_name = libs[i]
        logging.debug(f'Parsing metadata {lib_name} ({i + 1}/{len(lib_dirs)})')
        # parse metadata
        meta[lib_name] = parse_metadata(lib_dir)
        # parse gene counts
        logging.debug(f'Writing gene counts {lib_name} ({i + 1}/{len(lib_dirs)})')
        filename = lib_dir / 'hulkrna' / HULKRNA_COUNT_GENE_FILE
        g_counts_df = pd.read_feather(filename)
        hdata.add_counts(i, g_counts_df)

    logging.info('Adding metadata to data object')
    df = pd.DataFrame.from_dict(meta, orient='index')
    hdata.add_metadata(df, axis=1)
    hdata.close()
    return 0


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', dest='num_cpus', type=int, default=1)
    parser.add_argument('--mem', dest='mem_gb', type=int, default=16)
    parser.add_argument('--index', dest='index_dir', required=True, 
                        help='directory with index files')
    parser.add_argument('-d', '--search-dir', dest='search_dirs', nargs='+', action='extend',
                        help='directories to search for library results')
    parser.add_argument('--libs', dest='libs_file', default=None,
                        help='text file with list of libraries to aggregate (one library per line)')
    parser.add_argument('-o', dest='output_file', required=True,
                        help='HDF5 output file')
    args = parser.parse_args()

    # check args
    for d in args.search_dirs:
        if not Path(d).exists():
            parser.error(f'search dir {d} not found')
    if args.libs_file is not None:
        if not Path(args.libs_file).exists():
            parser.error(f'libraries file {args.libs_file} not found')
    if not Path(args.index_dir).exists():
        parser.error(f'index directory {args.index_dir} not found')

    # setup output files
    return(hulkrna_gather(
        output_file = Path(args.output_file), 
        index_dir = Path(args.index_dir), 
        search_dirs = args.search_dirs, 
        libs_file = args.libs_file
    ))


if __name__ == '__main__':
    sys.exit(main())
