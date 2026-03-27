import argparse
import gzip
from collections import Counter
from typing import Optional

RMSK_HEADER_ROWS = 3
RMSK_FIELD_NAMES = [
    'sw_score', 'perc_div', 'perc_del', 'perc_ins', 'qname', 'qstart',
    'qend', 'qleft', 'strand', 'repname', 'repclassfamily', 'repstart',
    'repend', 'repleft', 'id'
]


def _parse_rmsk_line(line: str) -> Optional[dict]:
    parts = line.strip().split()
    if len(parts) < len(RMSK_FIELD_NAMES):
        return None
    return dict(zip(RMSK_FIELD_NAMES, parts[:len(RMSK_FIELD_NAMES)]))


def count_classes(rmsk_path: str):
    class_counts = Counter()
    family_counts = Counter()

    open_func = gzip.open if rmsk_path.endswith('.gz') else open
    with open_func(rmsk_path, 'rt') as fh:
        for _ in range(RMSK_HEADER_ROWS):
            if not fh.readline():
                break

        for line in fh:
            row = _parse_rmsk_line(line)
            if row is None:
                continue
            repclass, repfamily = (row['repclassfamily'].split('/') + [None])[:2]
            if repclass:
                class_counts[repclass] += 1
            if repfamily:
                family_counts[repfamily] += 1

    return class_counts, family_counts


def main():
    parser = argparse.ArgumentParser(description="Count RepeatMasker classes and families")
    parser.add_argument('rmsk', help='RepeatMasker output file (plain or gzipped)')
    args = parser.parse_args()

    classes, families = count_classes(args.rmsk)

    print('RepeatMasker class counts:')
    for name, count in classes.most_common():
        print(f'{name}\t{count}')
    print('\nRepeatMasker family counts:')
    for name, count in families.most_common():
        print(f'{name}\t{count}')


if __name__ == '__main__':
    main()
