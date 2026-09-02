#!/usr/bin/env python3
"""
detect_genome_duplication.py — detect and correct whole-genome tandem
duplications flagged by CheckV's kmer_freq metric.

CheckV computes kmer_freq as a measure of internal sequence repetition;
values near 2.0 indicate the assembly contains two tandem copies of the
same genome (a common assembly artifact, especially for small/circular
viral genomes reproduced identically across multiple samples).

For each flagged sequence, this script:
  1. Self-aligns the sequence against itself with nucmer
  2. Excludes the trivial full-length self-hit
  3. Looks for the best non-trivial self-alignment block consistent with
     a clean tandem duplication (starts near position 1, ends near the
     full length, contiguous copies, high identity, block length close
     to total_length / round(kmer_freq))
  4. If validated, extracts a single copy using the empirical alignment
     offset as the repeat period (not just length/2)
  5. If not validated, leaves the sequence untouched and flags it for
     manual review

Outputs:
  deduplicated_all.fasta    — all input sequences; corrected ones replaced
                              by their single copy, everything else unchanged
  corrected_only.fasta      — only the successfully corrected sequences
  flagged_for_review.fasta  — sequences that looked duplicated but did not
                              cleanly validate as a simple tandem repeat
  dedup_report.tsv          — per-sequence metrics and decision
  plots/<seq_id>.png        — self-alignment dot plot per flagged sequence,
                              with the selected repeat block highlighted

Usage:
  detect_genome_duplication.py \\
      --fasta       votus.fasta \\
      --checkv      quality_summary.tsv \\
      --outdir      dedup_out \\
      --threshold   1.2 \\
      --min-identity 95 \\
      --max-period-deviation 0.05
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Probe-vs-full alignment: find repeated copies without guessing boundaries
# ---------------------------------------------------------------------------

def align_probe_vs_full(probe_seq, full_seq, tmpdir):
    """
    Align a short probe (subsequence from the start of the contig) against
    the full original sequence. Because the probe is genuinely shorter than
    and distinct in length from the full sequence, nucmer reports each
    repeated occurrence as a separate hit — this avoids the trivial
    full-length self-match that occurs when a sequence is aligned against
    an identical copy of itself.
    """
    tmp = Path(tmpdir)
    probe_path = tmp / "probe.fasta"
    full_path  = tmp / "full.fasta"
    probe_path.write_text(f">probe\n{probe_seq}\n")
    full_path.write_text(f">full\n{full_seq}\n")

    prefix = tmp / "probe_vs_full"
    subprocess.run(
        ["nucmer", "--maxmatch", "-p", str(prefix), str(probe_path), str(full_path)],
        check=True, capture_output=True, text=True
    )
    coords_txt = subprocess.run(
        ["show-coords", "-r", "-T", f"{prefix}.delta"],
        check=True, capture_output=True, text=True
    ).stdout

    rows = []
    for line in coords_txt.splitlines():
        parts = line.split('\t')
        if len(parts) < 9:
            continue
        try:
            s1, e1, s2, e2, len1, len2, idy = (
                int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]),
                int(parts[4]), int(parts[5]), float(parts[6])
            )
        except ValueError:
            continue
        rows.append(dict(S1=s1, E1=e1, S2=s2, E2=e2, LEN1=len1, LEN2=len2, IDY=idy))

    return pd.DataFrame(rows)


def evaluate_tandem_hypothesis(seq_str, expected_copies, tmpdir,
                                min_identity, min_coverage):
    """
    Probe-based detection of clean whole-genome duplication, including
    inverted (reverse-complement) duplications.

    Takes a probe from the start of the sequence — sized at half the
    expected single-copy length, so it safely fits within one copy even
    if the period estimate is imprecise — and aligns it against the full
    original sequence with nucmer, which searches both strands by default.

    nucmer's show-coords convention: a hit on the plus strand of the
    reference has S2 < E2; a hit on the minus strand (i.e. the probe
    matches a reverse-complemented copy) has S2 > E2. Both are treated
    as valid evidence of a duplicate copy.

    If nucmer finds `expected_copies` well-separated, high-identity,
    high-coverage hits (one per copy, each covering nearly the entire
    probe, evenly spaced), the sequence is a clean duplication — tandem
    or inverted — and the copy boundary is read directly from the hit
    coordinates, not guessed. Only the first (forward, reference) copy
    is kept; any inverted or repeated copies are discarded as redundant.

    A rotated duplication (e.g. a circular genome linearized at different
    points in each copy) produces fragmented, partial-coverage hits and
    is correctly rejected rather than mis-corrected.

    Returns a dict describing the decision and supporting metrics.
    """
    total_len = len(seq_str)
    period_estimate = total_len // expected_copies
    probe_len = max(50, period_estimate // 2)

    if probe_len < 50:
        return dict(decision='flagged_for_review',
                    reason=f'estimated period too short ({period_estimate} bp) to probe reliably')

    probe = seq_str[:probe_len]
    coords = align_probe_vs_full(probe, seq_str, tmpdir)

    if coords.empty:
        return dict(decision='no_repeat_found',
                    reason='probe did not align anywhere (unexpected)', hits=[])

    # Normalise strand: for reverse hits S2 > E2. Record strand and use
    # min/max for position bookkeeping regardless of orientation.
    coords['strand']    = np.where(coords.S2 <= coords.E2, '+', '-')
    coords['ref_start']  = coords[['S2', 'E2']].min(axis=1)
    coords['ref_end']    = coords[['S2', 'E2']].max(axis=1)
    coords['coverage']  = coords['LEN1'] / probe_len

    passing = coords[
        (coords.IDY >= min_identity) & (coords.coverage >= min_coverage)
    ].copy()
    passing = passing.sort_values('ref_start').reset_index(drop=True)

    hits_summary = [
        dict(ref_start=int(r.ref_start), ref_end=int(r.ref_end), strand=r.strand,
            identity=float(r.IDY), coverage=float(r.coverage))
        for _, r in coords.iterrows()
    ]

    if len(passing) < expected_copies:
        return dict(
            decision='flagged_for_review',
            reason=(f"only {len(passing)}/{expected_copies} full-coverage, "
                    f"high-identity probe occurrences found — repeat structure "
                    f"is not a clean duplication"),
            hits=hits_summary,
        )

    # Derive the true copy period from each additional-copy hit's position.
    # This must be strand-aware:
    #   '+' strand hit (copy i, 1-indexed, i=2..N): probe re-occurs at the
    #     START of copy i, so period = (ref_start - 1) / (i - 1)
    #   '-' strand hit (copy i is reverse-complemented relative to copy 1):
    #     the forward probe matches copy i's reverse complement, which
    #     places the match at the END of copy i's span, not its start.
    #     Working through the coordinate algebra: ref_end = i * period,
    #     so period = ref_end / i
    # All copies are assumed equal length; the derived periods from each
    # hit should agree closely if this is a clean, uniform duplication.
    # The hit with the smallest ref_start is always the probe matching its
    # own source position (copy 1 itself) — this is expected and not one of
    # the "extra" duplicate copies we need to validate.
    passing = passing.sort_values('ref_start').reset_index(drop=True)
    extra_hits = passing.iloc[1:expected_copies].reset_index(drop=True)

    candidate_periods = []
    for i, row in enumerate(extra_hits.itertuples(), start=2):  # i = copy index, 1-based
        if row.strand == '+':
            # ref_start is 1-based; copy 1 occupies 0-based [0, ref_start-1),
            # so its length (period) is ref_start - 1
            candidate_periods.append((row.ref_start - 1) / (i - 1))
        else:
            # ref_end (1-based) = i * period exactly for a reverse-strand
            # copy at position i (derived from coordinate algebra above)
            candidate_periods.append(row.ref_end / i)

    period_cv = (np.std(candidate_periods) / np.mean(candidate_periods)) if candidate_periods else 1.0

    if period_cv > 0.02:  # coefficient of variation — derived periods should agree closely
        return dict(
            decision='flagged_for_review',
            reason=(f"probe occurrences found but derived periods are "
                    f"inconsistent ({[round(p) for p in candidate_periods]}, "
                    f"CV={period_cv:.1%}) — possible rotated or irregular "
                    f"repeat structure"),
            hits=hits_summary,
        )

    true_period = round(np.mean(candidate_periods))

    # Final sanity check: total length should be close to N * period
    length_check = abs(total_len - expected_copies * true_period) / total_len
    if length_check > 0.02:
        return dict(
            decision='flagged_for_review',
            reason=(f"derived period ({true_period}) x expected_copies "
                    f"({expected_copies}) = {expected_copies * true_period}, "
                    f"which deviates {length_check:.1%} from total length "
                    f"({total_len})"),
            hits=hits_summary,
        )

    worst_identity = extra_hits['IDY'].min()
    worst_coverage = extra_hits['coverage'].min()
    strands_found  = extra_hits['strand'].tolist()

    return dict(
        decision='corrected',
        period=true_period,
        identity=float(worst_identity),
        coverage=float(worst_coverage),
        strands=strands_found,
        hits=hits_summary,
    )


def plot_probe_hits(seq_id, total_len, probe_len, result, out_path):
    """
    Visualize where the probe aligns along the full sequence: one bar per
    hit, colored by strand (+/-), showing coverage and identity. Confirms
    at a glance whether copies are evenly spaced and full coverage.
    """
    hits = result.get('hits', [])
    fig, ax = plt.subplots(figsize=(8, 2.5))

    ax.axhline(0, color='lightgray', linewidth=1, zorder=0)
    ax.add_patch(plt.Rectangle((0, -0.15), total_len, 0.3,
                               facecolor='whitesmoke', edgecolor='gray', linewidth=0.5))

    for h in hits:
        color = 'steelblue' if h['strand'] == '+' else 'darkorange'
        passed = h['identity'] >= 95 and h['coverage'] >= 0.90
        alpha = 1.0 if passed else 0.4
        ax.add_patch(plt.Rectangle(
            (h['ref_start'], -0.15), h['ref_end'] - h['ref_start'], 0.3,
            facecolor=color, alpha=alpha, edgecolor='black', linewidth=0.5
        ))
        ax.text((h['ref_start'] + h['ref_end']) / 2, 0.22,
                f"{h['strand']} {h['identity']:.1f}%\ncov={h['coverage']:.0%}",
                ha='center', fontsize=7)

    ax.set_xlim(-total_len * 0.02, total_len * 1.02)
    ax.set_ylim(-0.5, 0.6)
    ax.set_yticks([])
    ax.set_xlabel('Position (bp)')
    ax.set_title(
        f"{seq_id}  (length={total_len}, probe_len={probe_len})\n"
        f"decision={result['decision']}",
        fontsize=9
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Detect and correct whole-genome duplications (tandem or "
                    "inverted) via probe-vs-full-sequence alignment."
    )
    p.add_argument('--fasta',    required=True)
    p.add_argument('--checkv',   required=True, help='CheckV quality_summary.tsv')
    p.add_argument('--outdir',   required=True)
    p.add_argument('--threshold', type=float, default=1.2,
                   help='kmer_freq threshold to flag a sequence (default: 1.2)')
    p.add_argument('--min-identity', type=float, default=95.0,
                   help='Minimum %% identity for a probe hit to count as a copy (default: 95)')
    p.add_argument('--min-coverage', type=float, default=0.90,
                   help='Minimum fractional coverage of the probe for a hit to count (default: 0.90)')
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plots_dir = outdir / 'plots'
    plots_dir.mkdir(exist_ok=True)

    checkv_df = pd.read_csv(args.checkv, sep='\t')
    kmer_freq = dict(zip(checkv_df['contig_id'].astype(str),
                        pd.to_numeric(checkv_df['kmer_freq'], errors='coerce')))

    records = {r.id: r for r in SeqIO.parse(args.fasta, 'fasta')}

    flagged = [sid for sid, kf in kmer_freq.items()
              if sid in records and pd.notna(kf) and kf > args.threshold]

    print(f"[detect_genome_duplication] {len(records)} total sequences, "
          f"{len(flagged)} flagged (kmer_freq > {args.threshold})", file=sys.stderr)

    report_rows = []
    corrected_records   = []
    flagged_review_recs = []
    final_records       = dict(records)  # will be mutated for corrected seqs

    with tempfile.TemporaryDirectory() as tmpdir:
        for seq_id in flagged:
            rec = records[seq_id]
            seq_str = str(rec.seq)
            total_len = len(seq_str)
            kf = kmer_freq[seq_id]
            expected_copies = max(2, round(kf))
            probe_len = max(50, (total_len // expected_copies) // 2)

            result = evaluate_tandem_hypothesis(
                seq_str, expected_copies, tmpdir,
                args.min_identity, args.min_coverage
            )

            plot_probe_hits(seq_id, total_len, probe_len, result,
                            plots_dir / f"{seq_id}.png")

            row = dict(
                seq_id=seq_id, original_length=total_len, kmer_freq=round(kf, 3),
                expected_copies=expected_copies, decision=result['decision'],
            )

            if result['decision'] == 'corrected':
                period = result['period']
                corrected_seq = seq_str[:period]
                corrected_rec = rec.__class__(
                    seq=rec.seq.__class__(corrected_seq), id=seq_id, description=''
                )
                corrected_records.append(corrected_rec)
                final_records[seq_id] = corrected_rec

                strands = result.get('strands', [])
                orientation = (
                    'inverted' if any(s == '-' for s in strands) else 'tandem'
                )

                row.update(
                    corrected_length=period,
                    identity_pct=round(result['identity'], 2),
                    coverage_pct=round(result['coverage'] * 100, 2),
                    orientation=orientation,
                    notes='auto-corrected: single copy extracted',
                )
            else:
                flagged_review_recs.append(rec)
                row.update(
                    corrected_length=None, identity_pct=None,
                    coverage_pct=None, orientation=None,
                    notes=result.get('reason', result['decision']),
                )

            report_rows.append(row)

    # --- Write outputs ---
    with open(outdir / 'deduplicated_all.fasta', 'w') as f:
        SeqIO.write(final_records.values(), f, 'fasta')

    with open(outdir / 'corrected_only.fasta', 'w') as f:
        SeqIO.write(corrected_records, f, 'fasta')

    with open(outdir / 'flagged_for_review.fasta', 'w') as f:
        SeqIO.write(flagged_review_recs, f, 'fasta')

    report_df = pd.DataFrame(report_rows, columns=[
        'seq_id', 'original_length', 'kmer_freq', 'expected_copies', 'decision',
        'corrected_length', 'orientation', 'identity_pct', 'coverage_pct', 'notes'
    ])
    report_df.to_csv(outdir / 'dedup_report.tsv', sep='\t', index=False)

    n_corrected = len(corrected_records)
    n_flagged   = len(flagged_review_recs)
    print(
        f"[detect_genome_duplication] {n_corrected} corrected, "
        f"{n_flagged} flagged for manual review",
        file=sys.stderr
    )
    print(f"[detect_genome_duplication] outputs written to {outdir}/", file=sys.stderr)


if __name__ == '__main__':
    main()
