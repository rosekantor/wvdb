#!/usr/bin/env python3
"""
parse_clusters_v2.py — extract contig pairs and candidate IDs from vclust cluster output.

Modes:
  rank12:   longest vs second-longest (original branch A behaviour)
  rank23:   second-longest vs third-longest, clusters with >=3 contigs only
            (original branch B behaviour; supports --restrict-reps)
  rank123:  all three ranks in one pass — produces pairs12, pairs13, pairs23,
            singletons, and a candidate_ids file for seqkit grep + minimap2.
            This is the recommended mode for the minimap2-based trimming pipeline.

rank12 and rank23 modes are retained for backward compatibility.
"""

import argparse
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Shared helpers (unchanged from original)
# ---------------------------------------------------------------------------

def read_id_set(path):
    """Read newline-delimited IDs into a set (ignores blank lines)."""
    ids = set()
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                ids.add(s)
    return ids


def get_lengths(length_file):
    """Read tab-delimited contig/length table produced by seqkit fx2tab -nl."""
    contig_lengths = {}
    with open(length_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            contig, length = line.split('\t')
            if contig == 'name':   # drop header if present
                continue
            contig_lengths[contig] = int(length)
    return contig_lengths


def get_clusters(cluster_file, contig_lengths, restrict_reps=None):
    """
    Read vclust clusters.tsv into a dict:
      clusters[rep] = ([name, ...], [length, ...])
    The rep is the vclust-chosen representative (longest contig) and is
    also present as one of the member names.

    If restrict_reps is provided, only clusters whose rep is in that set
    are retained.
    """
    clusters = defaultdict(lambda: [[], []])
    with open(cluster_file) as f:
        for line in f:
            if line.startswith('object'):
                continue
            line = line.strip()
            if not line:
                continue

            target, rep = line.split('\t')

            if restrict_reps is not None and rep not in restrict_reps:
                continue

            if target not in contig_lengths:
                raise KeyError(
                    f"Contig '{target}' found in cluster file but missing from "
                    f"length file. Check that the same FASTA was used for both."
                )

            clusters[rep][0].append(target)
            clusters[rep][1].append(contig_lengths[target])

    return clusters


def get_sorted_by_length(names, lengths):
    """Return [(name, length), ...] sorted by length descending."""
    return sorted(zip(names, lengths), key=lambda x: x[1], reverse=True)


# ---------------------------------------------------------------------------
# Original single-mode logic (rank12 / rank23) — unchanged
# ---------------------------------------------------------------------------

def make_contig_list(clusters, mode):
    """
    Build singletons list and pairs list for rank12 or rank23 mode.
    Returns (singletons, pairs) where pairs is a list of [name_a, name_b].
    """
    if mode == 'rank12':
        ranks = [0, 1]
        min_size = 2
    elif mode == 'rank23':
        ranks = [1, 2]
        min_size = 3
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    singletons = []
    pairs = []

    for rep, (names, lengths) in clusters.items():
        n = len(names)
        if n == 1:
            singletons.append(rep)
        if n < min_size:
            continue
        sorted_contigs = get_sorted_by_length(names, lengths)
        if max(ranks) >= len(sorted_contigs):
            continue
        pairs.append([sorted_contigs[r][0] for r in ranks])

    return singletons, pairs


# ---------------------------------------------------------------------------
# New rank123 logic
# ---------------------------------------------------------------------------

def make_rank123_lists(clusters):
    """
    Single-pass extraction of all three ranks per cluster.

    Returns:
      singletons    list[str]         — clusters with exactly 1 member
      pairs12       list[[r1, r2]]    — all clusters with >=2 members
      pairs13       list[[r1, r3]]    — clusters with >=3 members
      pairs23       list[[r2, r3]]    — clusters with >=3 members
      candidate_ids set[str]          — union of all rank1/2/3 IDs

    Notes:
      - pairs12, pairs13, pairs23 all use the same [query, target] column order
        so trim_genomes.py can treat them uniformly.
      - candidate_ids is used for the seqkit grep step that precedes minimap2.
        It is intentionally a set (no duplicates) even though the same ID may
        appear in multiple pair lists (e.g. rank2 is query in pairs23 but
        target in pairs12).
    """
    singletons = []
    pairs12 = []
    pairs13 = []
    pairs23 = []
    candidate_ids = set()

    for rep, (names, lengths) in clusters.items():
        n = len(names)

        if n == 1:
            singletons.append(rep)
            # singletons go to branch C; not included in minimap2 candidates
            continue

        sorted_contigs = get_sorted_by_length(names, lengths)
        r1 = sorted_contigs[0][0]
        r2 = sorted_contigs[1][0]

        # All clusters with >=2 members get a 1v2 pair
        pairs12.append([r1, r2])
        candidate_ids.update([r1, r2])

        if n >= 3:
            r3 = sorted_contigs[2][0]
            pairs13.append([r1, r3])
            pairs23.append([r2, r3])
            candidate_ids.add(r3)

    return singletons, pairs12, pairs13, pairs23, candidate_ids


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def write_pairs(path, pairs):
    with open(path, 'w') as f:
        for pair in pairs:
            f.write('\t'.join(pair) + '\n')


def write_ids(path, ids):
    with open(path, 'w') as f:
        for id_ in sorted(ids):
            f.write(id_ + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Extract contig pairs from vclust cluster output for trimming."
    )
    parser.add_argument(
        '-c', '--cluster_file', required=True,
        help='vclust clusters TSV (columns: object<TAB>rep)'
    )
    parser.add_argument(
        '-l', '--length_file', required=True,
        help='TSV of contig names and lengths (seqkit fx2tab -nl output)'
    )
    parser.add_argument(
        '-m', '--pair_mode', required=True,
        choices=['rank12', 'rank23', 'rank123'],
        help=(
            'rank12: longest vs second-longest (branch A); '
            'rank23: second vs third-longest (branch B); '
            'rank123: all three ranks in one pass (recommended for minimap2 pipeline)'
        )
    )

    # --- rank12 / rank23 mode arguments (original interface, unchanged) ---
    parser.add_argument(
        '-o', '--out_pairs',
        help='[rank12/rank23] output file for contig pairs'
    )
    parser.add_argument(
        '-s', '--out_singletons',
        help='output file for singleton contig IDs'
    )
    parser.add_argument(
        '--restrict-reps',
        help='[rank23] file of rep IDs; only clusters whose rep is in this list are kept'
    )

    # --- rank123 mode arguments ---
    parser.add_argument(
        '--out-pairs12',
        help='[rank123] output TSV for rank1 vs rank2 pairs'
    )
    parser.add_argument(
        '--out-pairs13',
        help='[rank123] output TSV for rank1 vs rank3 pairs'
    )
    parser.add_argument(
        '--out-pairs23',
        help='[rank123] output TSV for rank2 vs rank3 pairs'
    )
    parser.add_argument(
        '--out-candidate-ids',
        help='[rank123] output file of all rank1/2/3 IDs for seqkit grep'
    )

    args = parser.parse_args()

    # --- Argument validation ---
    if args.pair_mode in ('rank12', 'rank23'):
        if not args.out_pairs:
            parser.error(f'-o/--out_pairs is required for --pair_mode {args.pair_mode}')
        if not args.out_singletons:
            parser.error(f'-s/--out_singletons is required for --pair_mode {args.pair_mode}')
        if args.pair_mode != 'rank23' and args.restrict_reps:
            parser.error('--restrict-reps is only valid with --pair_mode rank23')

    if args.pair_mode == 'rank123':
        missing = [
            flag for flag, val in [
                ('--out-pairs12',      args.out_pairs12),
                ('--out-pairs13',      args.out_pairs13),
                ('--out-pairs23',      args.out_pairs23),
                ('--out-candidate-ids', args.out_candidate_ids),
                ('-s/--out_singletons', args.out_singletons),
            ] if not val
        ]
        if missing:
            parser.error(f'rank123 mode requires: {", ".join(missing)}')
        if args.restrict_reps:
            parser.error('--restrict-reps is not used in rank123 mode')

    # --- Load data ---
    contig_lengths = get_lengths(args.length_file)

    restrict_reps = None
    if args.pair_mode == 'rank23' and args.restrict_reps:
        restrict_reps = read_id_set(args.restrict_reps)

    clusters = get_clusters(args.cluster_file, contig_lengths, restrict_reps=restrict_reps)

    # --- Run selected mode ---
    if args.pair_mode in ('rank12', 'rank23'):
        singletons, pairs = make_contig_list(clusters, args.pair_mode)

        print(f'total clusters: {len(clusters)}')
        print(f'singletons: {len(singletons)}')
        print(f'pairs ({args.pair_mode}): {len(pairs)}')
        if args.pair_mode == 'rank23' and restrict_reps is not None:
            print(f'restrict_reps applied: {len(restrict_reps)} IDs')

        with open(args.out_singletons, 'w') as f:
            if singletons:
                f.write('\n'.join(singletons) + '\n')
        write_pairs(args.out_pairs, pairs)

    else:  # rank123
        singletons, pairs12, pairs13, pairs23, candidate_ids = make_rank123_lists(clusters)

        n_2member = len(pairs12) - len(pairs13)   # clusters with exactly 2 members
        print(f'total clusters: {len(clusters)}', file=sys.stderr)
        print(f'singletons: {len(singletons)}', file=sys.stderr)
        print(f'2-member clusters (pairs12 only): {n_2member}', file=sys.stderr)
        print(f'3+-member clusters: {len(pairs13)}', file=sys.stderr)
        print(f'pairs12: {len(pairs12)}', file=sys.stderr)
        print(f'pairs13: {len(pairs13)}', file=sys.stderr)
        print(f'pairs23: {len(pairs23)}', file=sys.stderr)
        print(f'minimap2 candidate sequences: {len(candidate_ids)}', file=sys.stderr)

        with open(args.out_singletons, 'w') as f:
            if singletons:
                f.write('\n'.join(singletons) + '\n')

        write_pairs(args.out_pairs12, pairs12)
        write_pairs(args.out_pairs13, pairs13)
        write_pairs(args.out_pairs23, pairs23)
        write_ids(args.out_candidate_ids, candidate_ids)


if __name__ == '__main__':
    main()
