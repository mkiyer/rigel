import argparse
import collections
import sys

from hulkrna.transcript import Transcript

def main():
    parser = argparse.ArgumentParser(description="Enumerate genes by gene_type from a GTF file.")
    parser.add_argument("gtf_file", help="Path to the input GTF file")
    args = parser.parse_args()

    print(f"Reading {args.gtf_file}...")

    # 1. Read all transcripts using the provided Transcript class
    # This parses exons and constructs Transcript objects
    try:
        transcripts = Transcript.read_gtf(args.gtf_file)
    except Exception as e:
        print(f"Error reading GTF: {e}")
        sys.exit(1)

    # 2. Deduplicate genes
    # Transcript.read_gtf returns a list of transcripts. Multiple transcripts 
    # share the same gene_id. We map gene_id -> gene_type to get unique genes.
    genes = {}
    
    print("Processing transcripts...")
    for t in transcripts:
        # Access gene_id and gene_type attributes defined in Transcript class
        g_id = t.g_id
        g_type = t.g_type

        if g_id:
            # We assume all transcripts of the same gene share the same gene_type.
            genes[g_id] = g_type
        else:
            # Handle cases where gene_id might be missing (rare in standard GTFs)
            pass

    # 3. Count frequencies
    # Count the occurrences of each gene_type
    type_counts = collections.Counter(genes.values())

    # 4. Print results
    print(f"\nFound {len(genes)} unique genes.\n")
    print(f"{'GENE TYPE':<30} | {'COUNT':<10}")
    print("-" * 45)
    
    # Sort by count descending
    for g_type, count in type_counts.most_common():
        # Handle None values for gene_type just in case
        label = g_type if g_type else "N/A"
        print(f"{label:<30} | {count:<10}")

if __name__ == "__main__":
    main()