#!/usr/bin/env python3
"""
categorize_host_terms.py — LLM-based categorization of free-text host and
isolation-source strings pulled from core_nt BLAST hits (via Entrez in
merge_annotations.py).

This is a separate, cacheable step from guess_host.py's decision-tree
logic — the LLM categorization is the slow, non-deterministic, API-cost-
incurring part; the decision tree is cheap, deterministic pandas logic
that should be safe to rerun freely without re-querying an LLM. See
guess_host.py for the second step.

Caching: if --host-dict / --isolation-dict already exist (e.g. from a
previous run, or hand-edited to correct LLM mistakes), only genuinely new
unseen values are sent to the LLM. The updated dict is written back out
in full, so corrections you make by hand persist across reruns as long
as you keep passing the same dict file back in.

Categories:
  hosts_ntBlastHit_simple:            bacteria | vertebrates | plants |
                                       invertebrates | fungi | unknown
  isolation_source_ntBlastHit_simple: fecal associated | digestive tract
                                       non-fecal | plant-associated | other

Usage:
  categorize_host_terms.py \\
      --merged            merged_annotations.tsv \\
      --env-file           /path/to/.env \\
      --provider           anthropic \\
      --model              claude-sonnet-5 \\
      --host-dict          host_dict.tsv \\
      --isolation-dict     isolation_dict.tsv \\
      --out-host-dict      host_dict_updated.tsv \\
      --out-isolation-dict isolation_dict_updated.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from llm_utils import load_api_key, call_llm, parse_json_response


HOST_CATEGORIES = ['bacteria', 'vertebrates', 'plants', 'invertebrates', 'fungi', 'unknown']
ISOLATION_CATEGORIES = ['fecal associated', 'digestive tract non-fecal', 'plant-associated', 'other']

HOST_SYSTEM_PROMPT = f"""You are categorizing host organism names from GenBank
BLAST hit metadata into broad taxonomic groups for a viral genome database.

Valid categories (use exactly these strings): {HOST_CATEGORIES}

Rules:
- A host name may be a common name, scientific name, or a mix (e.g. "human",
  "Homo sapiens", "Sus scrofa", "rice", "Enterobacteriaceae").
- Bacterial genus/family names (e.g. "Enterobacteriaceae", "Pseudomonas
  putida") → bacteria.
- Diatoms, protists, and other unicellular eukaryotes that are not plants,
  fungi, or animals → invertebrates (as the closest broad non-vertebrate
  eukaryote bucket used by this pipeline).
- If you cannot confidently determine the category, use "unknown" — do not
  guess.

Respond with ONLY a JSON object mapping each input term (exactly as given)
to one of the valid categories. No other text, no markdown code fences.
"""

ISOLATION_SYSTEM_PROMPT = f"""You are categorizing isolation source
descriptions from GenBank metadata for a viral genome database, to flag
samples that may reflect dietary/ingested content rather than true host
association.

Valid categories (use exactly these strings): {ISOLATION_CATEGORIES}

Rules:
- Fecal/stool samples of any kind → "fecal associated"
- Other digestive tract samples (gut, intestine, rumen, oral, etc. that are
  NOT specifically fecal/stool) → "digestive tract non-fecal"
- Plant tissue, soil directly associated with a plant, rhizosphere, etc.
  → "plant-associated"
- Anything else, including missing/null/uninformative descriptions →
  "other"

Respond with ONLY a JSON object mapping each input term (exactly as given)
to one of the valid categories. No other text, no markdown code fences.
"""


def load_existing_dict(path):
    """Load a previously-saved category dict TSV. Empty dict if missing."""
    if not path:
        return {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    df = pd.read_csv(p, sep='\t')
    if df.shape[1] < 2:
        return {}
    key_col, val_col = df.columns[0], df.columns[1]
    return dict(zip(df[key_col].astype(str), df[val_col].astype(str)))


def validate_response(result, batch, valid_categories):
    """
    Strictly validate an LLM categorization response against the batch of
    terms that were sent and the list of allowed categories.

    A term is only accepted into the returned dict if:
      - it appears in the response (no missing keys)
      - the LLM did not invent extra keys not in the batch (those are
        ignored, not merged in)
      - its value is an exact match (case-sensitive) to one of
        valid_categories — no fuzzy matching, no fallback substitution
      - its value is not blank/None/whitespace-only

    Terms that fail any check are left OUT of the returned dict entirely
    — not defaulted to any category — so they remain "uncategorized" and
    will be retried on the next run rather than silently cached with a
    guessed value.

    Returns (accepted: dict, rejected: list of (term, reason) tuples).
    """
    accepted = {}
    rejected = []

    extra_keys = set(result.keys()) - set(batch)
    if extra_keys:
        print(f"    WARNING: LLM response included {len(extra_keys)} key(s) "
              f"not in the requested batch — ignoring: {sorted(extra_keys)}",
              file=sys.stderr)

    for term in batch:
        if term not in result:
            rejected.append((term, 'missing from response'))
            continue

        value = result[term]

        if value is None or (isinstance(value, str) and not value.strip()):
            rejected.append((term, 'blank/null value'))
            continue

        if not isinstance(value, str) or value not in valid_categories:
            rejected.append((term, f'invalid category: {value!r}'))
            continue

        accepted[term] = value

    return accepted, rejected


def categorize_new_terms(terms, existing_dict, system_prompt, valid_categories,
                          provider, model, api_key, base_url=None, batch_size=200):
    """
    Categorize only the terms not already in existing_dict. Returns the
    full updated dict (existing + newly categorized). Terms whose LLM
    response fails validation (wrong/missing category, blank value) are
    left out of the returned dict — not defaulted to any category — so
    they remain uncategorized and are retried on the next run rather
    than silently cached with a guessed or fallback value.
    """
    new_terms = sorted(set(terms) - set(existing_dict.keys()))
    updated = dict(existing_dict)
    all_rejected = []

    if not new_terms:
        print(f"  All {len(terms)} unique term(s) already categorized "
              f"(cached) — no LLM calls needed", file=sys.stderr)
        return updated, all_rejected

    print(f"  {len(new_terms)} new term(s) to categorize "
          f"({len(existing_dict)} already cached)", file=sys.stderr)

    for i in range(0, len(new_terms), batch_size):
        batch = new_terms[i:i + batch_size]
        print(f"  LLM batch {i // batch_size + 1}: {len(batch)} term(s)",
              file=sys.stderr)

        user_prompt = "Categorize these terms:\n" + "\n".join(
            f"- {t}" for t in batch
        )

        response_text = call_llm(
            provider=provider, model=model, api_key=api_key,
            system_prompt=system_prompt, user_prompt=user_prompt,
            base_url=base_url,
        )

        try:
            result = parse_json_response(response_text)
        except Exception as e:
            print(f"  WARNING: failed to parse LLM response as JSON "
                  f"({e}) — skipping this batch, terms will remain "
                  f"uncategorized: {batch}", file=sys.stderr)
            all_rejected.extend((t, 'unparseable batch response') for t in batch)
            continue

        accepted, rejected = validate_response(result, batch, valid_categories)
        updated.update(accepted)
        all_rejected.extend(rejected)

        if rejected:
            print(f"    {len(rejected)}/{len(batch)} term(s) in this batch "
                  f"failed validation and remain uncategorized:",
                  file=sys.stderr)
            for term, reason in rejected:
                print(f"      {term!r}: {reason}", file=sys.stderr)

    return updated, all_rejected


def main():
    p = argparse.ArgumentParser(
        description="LLM-categorize host and isolation-source free text "
                    "from core_nt BLAST hits, with caching."
    )
    p.add_argument('--merged',           required=True,
                   help='merged_annotations.tsv from merge_annotations.py')
    p.add_argument('--env-file',         required=True,
                   help='Path to .env file containing API key(s)')
    p.add_argument('--provider',         required=True,
                   choices=['anthropic', 'openai', 'openai_compatible'])
    p.add_argument('--model',            required=True,
                   help='Model name — must be specified explicitly, no default')
    p.add_argument('--llm-base-url',     default=None,
                   help='Base URL for --provider openai_compatible (e.g. an '
                        'internal enterprise LLM gateway). Required only for '
                        'that provider; ignored otherwise.')
    p.add_argument('--llm-api-key-env-var', default=None,
                   help='Override which .env variable holds the API key '
                        '(default: ANTHROPIC_API_KEY / OPENAI_API_KEY / '
                        'OPENAI_COMPATIBLE_API_KEY depending on --provider). '
                        'Use this if your organization names it differently.')
    p.add_argument('--host-dict',        default=None,
                   help='Existing host category dict TSV (cache), if any')
    p.add_argument('--isolation-dict',   default=None,
                   help='Existing isolation-source category dict TSV (cache), if any')
    p.add_argument('--out-host-dict',    required=True)
    p.add_argument('--out-isolation-dict', required=True)
    p.add_argument('--out-rejected-terms', default=None,
                   help='Optional: TSV of terms that failed validation and '
                        'remain uncategorized (missing/blank/invalid response)')
    args = p.parse_args()

    df = pd.read_csv(args.merged, sep='\t')

    if 'hosts_ntBlastHit' not in df.columns:
        print("Warning: 'hosts_ntBlastHit' column not found in merged "
              "annotations — core_nt BLASTn may not have been run. "
              "Writing empty category dicts.", file=sys.stderr)
        pd.DataFrame(columns=['term', 'category']).to_csv(
            args.out_host_dict, sep='\t', index=False)
        pd.DataFrame(columns=['term', 'category']).to_csv(
            args.out_isolation_dict, sep='\t', index=False)
        if args.out_rejected_terms:
            pd.DataFrame(columns=['term', 'category_type', 'reason']).to_csv(
                args.out_rejected_terms, sep='\t', index=False)
        return

    has_blast = df['hosts_ntBlastHit'].notna()
    host_terms = df.loc[has_blast, 'hosts_ntBlastHit'].dropna().unique().tolist()
    isolation_terms = (
        df.loc[has_blast, 'isolation_source_ntBlastHit'].dropna().unique().tolist()
        if 'isolation_source_ntBlastHit' in df.columns else []
    )

    if args.provider == 'openai_compatible' and not args.llm_base_url:
        print("ERROR: --llm-base-url is required when --provider openai_compatible",
              file=sys.stderr)
        sys.exit(1)

    api_key = load_api_key(args.env_file, args.provider, args.llm_api_key_env_var)

    print(f"Categorizing {len(host_terms)} unique host term(s)...", file=sys.stderr)
    existing_host_dict = load_existing_dict(args.host_dict)
    host_dict, host_rejected = categorize_new_terms(
        host_terms, existing_host_dict, HOST_SYSTEM_PROMPT, HOST_CATEGORIES,
        args.provider, args.model, api_key, base_url=args.llm_base_url,
    )

    print(f"Categorizing {len(isolation_terms)} unique isolation source term(s)...",
          file=sys.stderr)
    existing_isolation_dict = load_existing_dict(args.isolation_dict)
    isolation_dict, isolation_rejected = categorize_new_terms(
        isolation_terms, existing_isolation_dict, ISOLATION_SYSTEM_PROMPT,
        ISOLATION_CATEGORIES, args.provider, args.model, api_key,
        base_url=args.llm_base_url,
    )

    pd.DataFrame(
        list(host_dict.items()), columns=['hosts_ntBlastHit', 'hosts_ntBlastHit_simple']
    ).to_csv(args.out_host_dict, sep='\t', index=False)

    pd.DataFrame(
        list(isolation_dict.items()),
        columns=['isolation_source_ntBlastHit', 'isolation_source_ntBlastHit_simple']
    ).to_csv(args.out_isolation_dict, sep='\t', index=False)

    if args.out_rejected_terms:
        rejected_rows = (
            [{'term': t, 'category_type': 'host', 'reason': r} for t, r in host_rejected] +
            [{'term': t, 'category_type': 'isolation_source', 'reason': r} for t, r in isolation_rejected]
        )
        pd.DataFrame(rejected_rows, columns=['term', 'category_type', 'reason']).to_csv(
            args.out_rejected_terms, sep='\t', index=False)

    n_rejected = len(host_rejected) + len(isolation_rejected)
    print(f"[categorize_host_terms] {len(host_dict)} host terms → {args.out_host_dict}",
          file=sys.stderr)
    print(f"[categorize_host_terms] {len(isolation_dict)} isolation terms → "
          f"{args.out_isolation_dict}", file=sys.stderr)
    if n_rejected:
        print(f"[categorize_host_terms] {n_rejected} term(s) failed validation "
              f"and remain uncategorized (will retry next run)"
              + (f" — see {args.out_rejected_terms}" if args.out_rejected_terms else ""),
              file=sys.stderr)


if __name__ == '__main__':
    main()
