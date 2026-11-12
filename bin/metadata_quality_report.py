#!/usr/bin/env python3
"""Generate an actionable metadata quality report for Pathoplexus submissions.

Checks Pathoplexus metadata against the INRB linelist and produces a
plain-language report that can be sent directly to data generators.

  PART 1 — NEEDS YOUR INPUT
    - Sequences not in the linelist
    - Unresolved duplicates
    - Missing specimenCollectorSampleId (when linked via submissionId)
    - Date/location mismatches between Pathoplexus and linelist
    - Location issues (health area vs zone, unknown, missing)
    - Incomplete or missing collection dates

  PART 2 — FOR YOUR RECORDS (already resolved)
    Resolved duplicates with KEEP/REVOKE recommendations.

Non-consortium sequences are written to a separate file.

Usage:
    python bin/metadata_quality_report.py \
        -m data/20261029/metadata.tsv \
        -a data/aliases.csv \
        -l data/CONFIDENTIAL-INRB-BDBV_with_pathoplexus_ID-24Oct2026.xlsx \
        -s data/health_area_zone_province_slim.csv \
        -o data/20261029/metadata_quality_report.csv \
        --exclude-out data/20261029/exclude_auto.txt \
        -f data/20261029/aligned.fasta \
        --non-consortium-out data/20261029/non_consortium_sequences.csv
"""

import argparse
import sys
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


# ── Helpers ─────────────────────────────────────────────────────────

def acc_base(accession):
    s = str(accession)
    return s.rsplit('.', 1)[0] if '.' in s else s


def build_alias_map(aliases_df):
    amap = {}
    for _, row in aliases_df.iterrows():
        obs = str(row['observed_name']).strip()
        canon = str(row['canonical_nom']).strip()
        reason = str(row.get('mismatch_reason', '')) if pd.notna(row.get('mismatch_reason')) else ''
        notes = str(row.get('notes', '')) if pd.notna(row.get('notes')) else ''
        amap[obs.lower()] = (canon, reason, notes)
    return amap


def normalize_location(val, alias_map):
    if pd.isna(val) or str(val).strip() == '':
        return ''
    lookup = alias_map.get(str(val).strip().lower())
    if lookup:
        return lookup[0]
    return str(val).strip()


# ── Linelist linking ────────────────────────────────────────────────

def link_to_linelist(meta, linelist_path):
    """Link each accession to its linelist row. Returns a dict: accessionVersion -> {status, ll_sample_id, ll_row}."""
    ll = pd.read_excel(linelist_path)
    ll = ll.dropna(subset=['Pathoplexus ID'])

    ll_by_sample_id = {}
    ll_by_pp_base = {}
    for idx, row in ll.iterrows():
        sid = str(row['Sample ID']).strip()
        pp_base = acc_base(str(row['Pathoplexus ID']).strip())
        ll_by_sample_id[sid] = row
        ll_by_pp_base[pp_base] = row

    links = {}
    for _, r in meta.iterrows():
        acc = r['accessionVersion']
        specimen = str(r['specimenCollectorSampleId']).strip() if pd.notna(r['specimenCollectorSampleId']) else ''
        submission = str(r['submissionId']).strip() if pd.notna(r['submissionId']) else ''
        ab = acc_base(acc)

        if specimen and specimen in ll_by_sample_id:
            links[acc] = {'status': 'linked', 'match_key': 'specimenCollectorSampleId',
                          'll_sample_id': specimen, 'll_row': ll_by_sample_id[specimen]}
        elif submission and submission in ll_by_sample_id:
            links[acc] = {'status': 'linked_via_submission', 'match_key': 'submissionId',
                          'll_sample_id': submission, 'll_row': ll_by_sample_id[submission]}
        elif ab in ll_by_pp_base:
            links[acc] = {'status': 'linked', 'match_key': 'accessionVersion',
                          'll_sample_id': str(ll_by_pp_base[ab]['Sample ID']), 'll_row': ll_by_pp_base[ab]}
        else:
            links[acc] = {'status': 'not_linked', 'match_key': None,
                          'll_sample_id': None, 'll_row': None}

    n_linked = sum(1 for v in links.values() if v['status'] in ('linked', 'linked_via_submission'))
    n_not = sum(1 for v in links.values() if v['status'] == 'not_linked')
    print(f'Linelist linking: {n_linked} linked, {n_not} not linked (of {len(meta)})')

    return links, ll_by_sample_id


# ── Checkers ─────────────────────────────────────────────────────────

def check_duplicates(meta, links, ll_by_sample_id):
    issues = []
    seen_pairs = set()

    def _resolve_pair(accs, sample_id, group):
        # Tier 1: the linelist explicitly names one accession as its Pathoplexus ID
        ll_row = ll_by_sample_id.get(sample_id) if ll_by_sample_id else None
        if ll_row is not None:
            ll_pp_id = str(ll_row.get('Pathoplexus ID', '')).strip()
            ll_pp_base = acc_base(ll_pp_id) if ll_pp_id else ''
            listed = [a for a in accs if acc_base(a) == ll_pp_base]
            if len(listed) == 1:
                return 'resolved', listed[0], 'linelist match'

        # Tier 2: one accession is linked to any linelist entry, the other isn't
        in_ll = [a for a in accs if links.get(a, {}).get('status') in ('linked', 'linked_via_submission')]
        not_in_ll = [a for a in accs if a not in in_ll]
        if len(in_ll) == 1 and len(not_in_ll) >= 1:
            return 'resolved', in_ll[0], 'linelist match'

        # Both or neither linked — cannot resolve automatically.
        # Temporarily keep the higher-completeness one, but flag for team.
        rows = group[group['accessionVersion'].isin(accs)]
        temp_keep = None
        if 'completeness' in rows.columns and rows['completeness'].notna().all() and len(rows) >= 2:
            sorted_rows = rows.sort_values(
                ['completeness', 'totalAmbiguousNucs'],
                ascending=[False, True],
            )
            best = sorted_rows.iloc[0]
            second = sorted_rows.iloc[1]
            comps = {r['accessionVersion']: f"{r['completeness']:.3f}" for _, r in rows.iterrows()}
            reason = f"temporarily keeping {best['accessionVersion']} (completeness: {', '.join(f'{a}={c}' for a, c in comps.items())})"
            temp_keep = best['accessionVersion']
        else:
            reason = 'no linelist match — please advise which to keep'
            temp_keep = accs[0]

        return 'unresolved', temp_keep, reason

    def _process_group(sample_id, group_df, field_name):
        accs = group_df['accessionVersion'].tolist()
        pair_key = tuple(sorted(accs))
        if pair_key in seen_pairs:
            return
        seen_pairs.add(pair_key)

        status, best_acc, reason = _resolve_pair(accs, sample_id, group_df)

        if status == 'resolved':
            for acc in accs:
                keep = (acc == best_acc)
                issues.append({
                    'accessionVersion': acc,
                    'issue': 'duplicate_resolved',
                    'field': field_name,
                    'current_value': sample_id,
                    'suggested_fix': 'KEEP' if keep else 'REVOKE',
                    'details': f'keep {best_acc} ({reason})' if keep
                               else f'revoke (duplicate of {best_acc}, {reason})',
                })
        else:
            for acc in accs:
                keep = (acc == best_acc)
                issues.append({
                    'accessionVersion': acc,
                    'issue': 'duplicate_unresolved',
                    'field': field_name,
                    'current_value': sample_id,
                    'suggested_fix': 'TEMP_KEEP' if keep else 'TEMP_REVOKE',
                    'details': reason,
                })

    # Pass 1: group by specimenCollectorSampleId
    has_specimen = meta.dropna(subset=['specimenCollectorSampleId'])
    dup_mask = has_specimen['specimenCollectorSampleId'].duplicated(keep=False)
    dups1 = has_specimen[dup_mask]
    for sid, group in dups1.groupby('specimenCollectorSampleId'):
        _process_group(sid, group, 'specimenCollectorSampleId')

    # Pass 2: group by submissionId where at least one has blank specimenCollectorSampleId
    has_submission = meta.dropna(subset=['submissionId'])
    dup_sub_mask = has_submission['submissionId'].duplicated(keep=False)
    dups2 = has_submission[dup_sub_mask]
    for sid, group in dups2.groupby('submissionId'):
        if group['specimenCollectorSampleId'].isna().any():
            _process_group(sid, group, 'submissionId')

    return issues


def check_not_in_linelist(meta, links):
    issues = []
    for _, row in meta.iterrows():
        acc = row['accessionVersion']
        link = links.get(acc, {})
        if link.get('status') == 'not_linked':
            specimen = str(row['specimenCollectorSampleId']).strip() if pd.notna(row['specimenCollectorSampleId']) else ''
            submission = str(row['submissionId']).strip() if pd.notna(row['submissionId']) else ''
            has_id = bool(specimen) or bool(submission)
            if has_id:
                id_used = specimen if specimen else submission
                issues.append({
                    'accessionVersion': acc,
                    'issue': 'not_in_linelist',
                    'field': 'specimenCollectorSampleId' if specimen else 'submissionId',
                    'current_value': id_used,
                    'suggested_fix': 'Please add this sample to the linelist, or check for typos',
                    'details': '',
                })
    return issues


def check_sample_ids(meta, links):
    issues = []
    missing = meta[meta['specimenCollectorSampleId'].isna()]
    for _, row in missing.iterrows():
        acc = row['accessionVersion']
        link = links.get(acc, {})
        if link.get('status') == 'linked_via_submission':
            ll_sid = link['ll_sample_id']
            details = f'matched via submissionId to Sample ID {ll_sid}'
            fix = 'Please add specimenCollectorSampleId'
        else:
            details = ''
            fix = 'Please add the sample ID (e.g. BIA-XXXX) used during collection'
        issues.append({
            'accessionVersion': acc,
            'issue': 'missing_sample_id',
            'field': 'specimenCollectorSampleId',
            'current_value': '',
            'suggested_fix': fix,
            'details': details,
        })
    return issues


def check_linelist_mismatches(meta, links, alias_map):
    issues = []
    for _, row in meta.iterrows():
        acc = row['accessionVersion']
        link = links.get(acc, {})
        if link.get('status') not in ('linked', 'linked_via_submission'):
            continue
        ll_row = link['ll_row']
        if ll_row is None:
            continue

        # Date mismatch
        pp_date = str(row['sampleCollectionDate']).strip() if pd.notna(row['sampleCollectionDate']) else ''
        ll_date_raw = ll_row.get('Sampling date')
        if pd.notna(ll_date_raw):
            if hasattr(ll_date_raw, 'strftime'):
                ll_date = ll_date_raw.strftime('%Y-%m-%d')
            else:
                ll_date = str(ll_date_raw).strip()[:10]

            if pp_date and ll_date and len(pp_date) >= 10 and pp_date[:10] != ll_date[:10]:
                issues.append({
                    'accessionVersion': acc,
                    'issue': 'date_mismatch',
                    'field': 'sampleCollectionDate',
                    'current_value': pp_date,
                    'suggested_fix': f'Linelist says {ll_date} — please check which is correct',
                    'details': f'Pathoplexus: {pp_date}, linelist: {ll_date}',
                })

        # Location mismatch
        pp_loc = normalize_location(row['geoLocAdmin2'], alias_map)
        ll_zone = str(ll_row.get('Health Zone', '')).strip() if pd.notna(ll_row.get('Health Zone')) else ''
        if pp_loc and ll_zone and pp_loc.lower() != ll_zone.lower():
            # Also check alias of linelist zone
            ll_zone_normalized = normalize_location(ll_zone, alias_map)
            if pp_loc.lower() != ll_zone_normalized.lower():
                issues.append({
                    'accessionVersion': acc,
                    'issue': 'location_mismatch',
                    'field': 'geoLocAdmin2',
                    'current_value': str(row['geoLocAdmin2']).strip() if pd.notna(row['geoLocAdmin2']) else '',
                    'suggested_fix': f'Linelist says {ll_zone} — please check which is correct',
                    'details': f'Pathoplexus: {pp_loc}, linelist: {ll_zone}',
                })

    return issues


def check_admin2(meta, alias_map, province_names=None):
    canonical_values = {v[0] for v in alias_map.values()}
    lower_canonicals = {c.lower() for c in canonical_values}
    lower_provinces = {p.lower(): p for p in (province_names or [])}
    issues = []

    for _, row in meta.iterrows():
        acc = row['accessionVersion']
        val = row['geoLocAdmin2']

        if pd.isna(val) or str(val).strip() == '':
            issues.append({
                'accessionVersion': acc,
                'issue': 'missing_location',
                'field': 'geoLocAdmin2',
                'current_value': '',
                'suggested_fix': 'Please add the health zone (zone de sante)',
                'details': '',
            })
            continue

        val_stripped = str(val).strip()
        lookup = alias_map.get(val_stripped.lower())

        if lookup:
            canon, reason, notes = lookup
            if reason == 'health_area':
                issues.append({
                    'accessionVersion': acc,
                    'issue': 'health_area_not_zone',
                    'field': 'geoLocAdmin2',
                    'current_value': val_stripped,
                    'suggested_fix': canon,
                    'details': notes,
                })
        elif val_stripped.lower() in lower_provinces:
            province = lower_provinces[val_stripped.lower()]
            issues.append({
                'accessionVersion': acc,
                'issue': 'province_not_zone',
                'field': 'geoLocAdmin2',
                'current_value': val_stripped,
                'suggested_fix': f'Please provide the specific health zone within {province}',
                'details': f'"{val_stripped}" is a province, not a health zone',
            })
        elif val_stripped.lower() not in lower_canonicals:
            issues.append({
                'accessionVersion': acc,
                'issue': 'unknown_location',
                'field': 'geoLocAdmin2',
                'current_value': val_stripped,
                'suggested_fix': 'Please check — this does not match any known health zone',
                'details': '',
            })

    return issues


def check_dates(meta):
    issues = []
    for _, row in meta.iterrows():
        acc = row['accessionVersion']
        d = str(row['sampleCollectionDate']) if pd.notna(row['sampleCollectionDate']) else ''

        if not d or d == 'nan':
            issues.append({
                'accessionVersion': acc,
                'issue': 'missing_date',
                'field': 'sampleCollectionDate',
                'current_value': '',
                'suggested_fix': 'Please add the collection date (YYYY-MM-DD)',
                'details': '',
            })
        elif len(d) < 10:
            issues.append({
                'accessionVersion': acc,
                'issue': 'incomplete_date',
                'field': 'sampleCollectionDate',
                'current_value': d,
                'suggested_fix': 'Please provide the full date (YYYY-MM-DD)',
                'details': '',
            })

    return issues


# ── Sequence quality ────────────────────────────────────────────────

ADAR_PAIRS = {('A', 'G'), ('T', 'C')}


def read_fasta(path):
    seqs = {}
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if name is not None:
        seqs[name] = ''.join(parts)
    return seqs


def compute_consensus(seqs):
    if not seqs:
        return ''
    length = len(next(iter(seqs.values())))
    consensus = []
    for i in range(length):
        counts = Counter()
        for seq in seqs.values():
            base = seq[i] if i < len(seq) else 'N'
            if base not in ('N', '-'):
                counts[base] += 1
        if counts:
            consensus.append(counts.most_common(1)[0][0])
        else:
            consensus.append('N')
    return ''.join(consensus)


def check_sequence_quality(meta, fasta_path, completeness_threshold=0.70,
                           mismatch_threshold=0.10):
    seqs = read_fasta(fasta_path)
    consensus = compute_consensus(seqs)
    seq_len = len(consensus)

    issues = []
    for _, row in meta.iterrows():
        acc = row['accessionVersion']
        seq = seqs.get(acc)
        if seq is None:
            continue

        # Completeness check (from metadata field — more reliable than counting)
        completeness = row.get('completeness')
        if pd.notna(completeness) and float(completeness) < completeness_threshold:
            issues.append({
                'accessionVersion': acc,
                'issue': 'low_completeness',
                'field': 'completeness',
                'current_value': f'{float(completeness):.3f}',
                'suggested_fix': 'AUTO-EXCLUDED',
                'details': f'completeness {float(completeness):.1%} < {completeness_threshold:.0%} threshold',
            })

        # Non-ADAR mismatch check vs consensus
        aligned_positions = 0
        total_mismatches = 0
        adar_mismatches = 0
        for i in range(min(len(seq), seq_len)):
            s, c = seq[i], consensus[i]
            if s in ('N', '-') or c in ('N', '-'):
                continue
            aligned_positions += 1
            if s != c:
                total_mismatches += 1
                if (c, s) in ADAR_PAIRS:
                    adar_mismatches += 1

        non_adar = total_mismatches - adar_mismatches
        if aligned_positions > 0:
            non_adar_frac = non_adar / aligned_positions
            if non_adar_frac > mismatch_threshold:
                issues.append({
                    'accessionVersion': acc,
                    'issue': 'high_non_adar_mismatches',
                    'field': 'alignment',
                    'current_value': f'{non_adar}/{aligned_positions}',
                    'suggested_fix': 'AUTO-EXCLUDED',
                    'details': (f'{non_adar_frac:.1%} non-ADAR mismatches vs consensus '
                                f'(total: {total_mismatches}, ADAR: {adar_mismatches}, '
                                f'non-ADAR: {non_adar})'),
                })

    return issues


# ── Report formatting ────────────────────────────────────────────────

def _row_line(meta, acc, fields=('submissionId', 'sampleCollectionDate', 'geoLocAdmin2')):
    r = meta[meta['accessionVersion'] == acc]
    if r.empty:
        return acc
    r = r.iloc[0]
    parts = [acc]
    labels = {'submissionId': 'submission', 'sampleCollectionDate': 'date', 'geoLocAdmin2': 'location'}
    for f in fields:
        val = r.get(f, '')
        if pd.notna(val) and str(val).strip():
            parts.append(f'{labels.get(f, f)}: {val}')
    return '  |  '.join(parts)


def summarize(issues_df, meta, out_file=None):
    today = date.today().strftime('%Y-%m-%d')
    n_total = len(meta)

    lines = []

    def p(s=''):
        lines.append(s)

    p()
    p('=' * 72)
    p(f'  PATHOPLEXUS METADATA QUALITY REPORT')
    p(f'  Generated: {today}')
    p(f'  Accessions reviewed: {n_total}')
    p('=' * 72)

    if issues_df.empty:
        p('\n  No issues found. All metadata looks good!\n')
        p('=' * 72)
        text = '\n'.join(lines)
        print(text)
        if out_file:
            Path(out_file).write_text(text + '\n')
        return

    unresolved_types = [
        'not_in_linelist', 'duplicate_unresolved', 'missing_sample_id',
        'date_mismatch', 'location_mismatch',
        'health_area_not_zone', 'province_not_zone', 'unknown_location',
        'missing_location', 'incomplete_date', 'missing_date',
    ]
    auto_excluded_types = ['low_completeness', 'high_non_adar_mismatches']
    resolved_types = ['duplicate_resolved']

    part1 = issues_df[issues_df['issue'].isin(unresolved_types)]
    part2 = issues_df[issues_df['issue'].isin(resolved_types)]

    # ── PART 1 ──
    if not part1.empty:
        n1 = part1['accessionVersion'].nunique()
        p()
        p('─' * 72)
        p('  PART 1: NEEDS YOUR INPUT')
        p('  These issues require action from the data generating team.')
        p(f'  ({n1} accessions affected)')
        p('─' * 72)

        section = 0

        # 1. Not in linelist
        not_ll = part1[part1['issue'] == 'not_in_linelist']
        if not not_ll.empty:
            section += 1
            p(f'\n  {section}. NOT IN LINELIST ({len(not_ll)} accessions)')
            p('     These accessions have a sample/submission ID but we could')
            p('     not find a matching entry in the linelist.')
            p('     Please add them to the linelist or check for typos.\n')
            for _, row in not_ll.iterrows():
                id_val = row['current_value']
                p(f'     - {_row_line(meta, row["accessionVersion"])}  (ID: {id_val})')
            p()

        # 2. Unresolved duplicates
        dups_unresolved = issues_df[issues_df['issue'] == 'duplicate_unresolved']
        if not dups_unresolved.empty:
            n_groups = dups_unresolved['current_value'].nunique()
            section += 1
            p(f'\n  {section}. DUPLICATE SAMPLES — NEED CLARIFICATION ({n_groups} pairs)')
            p('     The following sample IDs appear on more than one accession,')
            p('     and neither is in the linelist. Please tell us which to keep.')
            p('     We have temporarily kept the higher-completeness one.\n')
            for sid in dups_unresolved['current_value'].unique():
                group = dups_unresolved[dups_unresolved['current_value'] == sid]
                p(f'     Sample ID: {sid}')
                for _, row in group.iterrows():
                    acc = row['accessionVersion']
                    label = 'TEMP KEEP' if row['suggested_fix'] == 'TEMP_KEEP' else 'TEMP REVOKE'
                    meta_row = meta[meta['accessionVersion'] == acc]
                    comp = f"{meta_row.iloc[0]['completeness']:.3f}" if not meta_row.empty and 'completeness' in meta_row.columns and pd.notna(meta_row.iloc[0]['completeness']) else '?'
                    p(f'       {label}:  {_row_line(meta, acc)}  (completeness: {comp})')
                p()

        # 3. Missing sample IDs
        missing_ids = part1[part1['issue'] == 'missing_sample_id']
        if not missing_ids.empty:
            section += 1
            linked_via_sub = missing_ids[missing_ids['details'].str.startswith('matched via', na=False)]
            truly_missing = missing_ids[~missing_ids['details'].str.startswith('matched via', na=False)]

            p(f'\n  {section}. MISSING SAMPLE IDS ({len(missing_ids)} accessions)')
            p('     These accessions do not have a specimenCollectorSampleId')
            p('     in their Pathoplexus metadata. Please add it.\n')

            if not linked_via_sub.empty:
                p(f'     Matched via submissionId ({len(linked_via_sub)} accessions):')
                p('     We identified these via their submissionId.\n')
                for _, row in linked_via_sub.iterrows():
                    p(f'     - {_row_line(meta, row["accessionVersion"])}')
                p()

            if not truly_missing.empty:
                p(f'     No match found ({len(truly_missing)} accessions):')
                p('     Without a sample ID, we cannot link them to the linelist.\n')
                for _, row in truly_missing.iterrows():
                    p(f'     - {_row_line(meta, row["accessionVersion"])}')
                p()

        # 4. Linelist mismatches
        date_mm = part1[part1['issue'] == 'date_mismatch']
        loc_mm = part1[part1['issue'] == 'location_mismatch']
        if not date_mm.empty or not loc_mm.empty:
            n_mm = len(date_mm) + len(loc_mm)
            section += 1
            p(f'\n  {section}. LINELIST MISMATCHES ({n_mm} discrepancies)')
            p('     The Pathoplexus metadata differs from the linelist.')
            p('     Please check which value is correct.\n')

            if not date_mm.empty:
                p(f'     Date mismatches ({len(date_mm)}):')
                for _, row in date_mm.iterrows():
                    p(f'     - {_row_line(meta, row["accessionVersion"])}')
                    p(f'       {row["details"]}')
                p()

            if not loc_mm.empty:
                p(f'     Location mismatches ({len(loc_mm)}):')
                for _, row in loc_mm.iterrows():
                    p(f'     - {_row_line(meta, row["accessionVersion"])}')
                    p(f'       {row["details"]}')
                p()

        # 5. Location issues
        ha_issues = part1[part1['issue'] == 'health_area_not_zone']
        prov_issues = part1[part1['issue'] == 'province_not_zone']
        unk_loc = part1[part1['issue'] == 'unknown_location']
        miss_loc = part1[part1['issue'] == 'missing_location']
        loc_total = len(ha_issues) + len(prov_issues) + len(unk_loc) + len(miss_loc)

        if loc_total > 0:
            section += 1
            p(f'\n  {section}. LOCATION ISSUES ({loc_total} accessions)')

            if not ha_issues.empty:
                by_fix = ha_issues.groupby(['current_value', 'suggested_fix']).agg(
                    count=('accessionVersion', 'size'),
                    accessions=('accessionVersion', list),
                ).reset_index()
                p(f'\n     Health area used instead of health zone ({len(ha_issues)} accessions):')
                p('     The location field should contain the health zone (zone de sante),')
                p('     not the health area (aire de sante).\n')
                for _, r in by_fix.iterrows():
                    p(f'     "{r["current_value"]}" should be "{r["suggested_fix"]}" ({r["count"]} accessions)')
                    for a in r['accessions'][:10]:
                        p(f'       - {_row_line(meta, a)}')
                    if r['count'] > 10:
                        p(f'       ... and {r["count"] - 10} more')
                    p()

            if not prov_issues.empty:
                by_val = prov_issues.groupby('current_value').agg(
                    count=('accessionVersion', 'size'),
                    accessions=('accessionVersion', list),
                ).reset_index()
                p(f'\n     Province name used instead of health zone ({len(prov_issues)} accessions):')
                p('     Please provide the specific health zone.\n')
                for _, r in by_val.iterrows():
                    p(f'     "{r["current_value"]}" — please provide the health zone ({r["count"]} accessions)')
                    for a in r['accessions'][:10]:
                        p(f'       - {_row_line(meta, a)}')
                    if r['count'] > 10:
                        p(f'       ... and {r["count"] - 10} more')
                    p()

            if not unk_loc.empty:
                by_val = unk_loc.groupby('current_value').agg(
                    count=('accessionVersion', 'size'),
                    accessions=('accessionVersion', list),
                ).reset_index()
                p(f'\n     Unrecognized location ({len(unk_loc)} accessions):')
                p('     Please check the spelling or let us know which health zone it belongs to.\n')
                for _, r in by_val.iterrows():
                    p(f'     "{r["current_value"]}" ({r["count"]} accessions)')
                    for a in r['accessions'][:10]:
                        p(f'       - {_row_line(meta, a)}')
                    if r['count'] > 10:
                        p(f'       ... and {r["count"] - 10} more')
                    p()

            if not miss_loc.empty:
                p(f'\n     Missing location ({len(miss_loc)} accessions):')
                p('     Please add the health zone (zone de sante).\n')
                for _, row in miss_loc.iterrows():
                    p(f'     - {_row_line(meta, row["accessionVersion"])}')
                p()

        # 6. Date issues
        inc_date = part1[part1['issue'] == 'incomplete_date']
        miss_date = part1[part1['issue'] == 'missing_date']
        date_total = len(inc_date) + len(miss_date)

        if date_total > 0:
            section += 1
            p(f'\n  {section}. DATE ISSUES ({date_total} accessions)')

            if not inc_date.empty:
                p(f'\n     Incomplete collection dates ({len(inc_date)} accessions):')
                p('     We need the exact day (YYYY-MM-DD) for our analyses.\n')
                for _, row in inc_date.iterrows():
                    p(f'     - {_row_line(meta, row["accessionVersion"])}  (current date: {row["current_value"]})')
                p()

            if not miss_date.empty:
                p(f'\n     Missing collection dates ({len(miss_date)} accessions):')
                p('     Please add the collection date (YYYY-MM-DD).\n')
                for _, row in miss_date.iterrows():
                    p(f'     - {_row_line(meta, row["accessionVersion"])}')
                p()

    # ── AUTO-EXCLUDED (quality) ──
    quality_excluded = issues_df[issues_df['issue'].isin(auto_excluded_types)]
    if not quality_excluded.empty:
        p()
        p('─' * 72)
        p('  AUTO-EXCLUDED: SEQUENCE QUALITY')
        p('  These sequences were automatically excluded from analysis.')
        p(f'  ({len(quality_excluded)} accessions)')
        p('─' * 72)

        low_comp = quality_excluded[quality_excluded['issue'] == 'low_completeness']
        high_mm = quality_excluded[quality_excluded['issue'] == 'high_non_adar_mismatches']

        if not low_comp.empty:
            p(f'\n     Low completeness ({len(low_comp)} accessions):')
            p('     Genome completeness below 70%.\n')
            for _, row in low_comp.iterrows():
                p(f'     - {_row_line(meta, row["accessionVersion"])}')
                p(f'       {row["details"]}')
            p()

        if not high_mm.empty:
            p(f'\n     High non-ADAR mismatches ({len(high_mm)} accessions):')
            p('     More than 10% of aligned positions differ from consensus')
            p('     (after excluding ADAR edits).\n')
            for _, row in high_mm.iterrows():
                p(f'     - {_row_line(meta, row["accessionVersion"])}')
                p(f'       {row["details"]}')
            p()

    # ── PART 2 ──
    if not part2.empty:
        keep_rows = part2[part2['suggested_fix'] == 'KEEP']
        revoke_rows = part2[part2['suggested_fix'] == 'REVOKE']
        n_pairs = revoke_rows.shape[0]

        p()
        p('─' * 72)
        p('  PART 2: FOR YOUR RECORDS (already resolved)')
        p('  These are duplicate submissions that we have already resolved.')
        p('  No urgent action is needed, but when you have time, please')
        p('  revoke the extra accessions on Pathoplexus.')
        p(f'  ({n_pairs} accessions to revoke)')
        p('─' * 72)

        for sid in revoke_rows['current_value'].unique():
            sid_keep = keep_rows[keep_rows['current_value'] == sid]
            sid_revoke = revoke_rows[revoke_rows['current_value'] == sid]

            p(f'\n     Sample ID: {sid}')
            if not sid_keep.empty:
                acc_keep = sid_keep.iloc[0]['accessionVersion']
                reason = sid_keep.iloc[0]['details']
                p(f'       KEEP:   {_row_line(meta, acc_keep)}')
                p(f'               ({reason})')
            for _, row in sid_revoke.iterrows():
                p(f'       REVOKE: {_row_line(meta, row["accessionVersion"])}')

        p()

    p('=' * 72)
    n_issues = issues_df['accessionVersion'].nunique()
    p(f'  Total accessions with issues: {n_issues} / {n_total}')
    if not part1.empty:
        p(f'  Issues requiring input: {part1["accessionVersion"].nunique()}')
    if not quality_excluded.empty:
        p(f'  Auto-excluded (quality): {quality_excluded["accessionVersion"].nunique()}')
    if not part2.empty:
        p(f'  Resolved duplicates (revoke when convenient): {part2[part2["suggested_fix"] == "REVOKE"].shape[0]}')
    p('=' * 72)
    p()

    text = '\n'.join(lines)
    print(text)
    if out_file:
        Path(out_file).write_text(text + '\n')
        print(f'Wrote report to {out_file}')


def main():
    parser = argparse.ArgumentParser(
        description='Generate actionable metadata quality report for Pathoplexus submissions.')
    parser.add_argument('--metadata', '-m', required=True,
                        help='Pathoplexus metadata TSV')
    parser.add_argument('--aliases', '-a', required=True,
                        help='Aliases CSV (observed_name -> canonical_nom)')
    parser.add_argument('--linelist', '-l', required=True,
                        help='INRB linelist Excel with Pathoplexus ID and Sample ID columns')
    parser.add_argument('--hz-lookup', '-s',
                        help='Health area/zone/province CSV (for province name detection)')
    parser.add_argument('-o', '--output', default='metadata_quality_report.csv',
                        help='Output CSV path for issue details')
    parser.add_argument('--report-txt',
                        help='Output path for human-readable report (default: same as -o but .txt)')
    parser.add_argument('--fasta', '-f',
                        help='Aligned FASTA for sequence quality checks (completeness, mismatch rate)')
    parser.add_argument('--exclude-out',
                        help='Output exclude list (resolved duplicates + quality failures)')
    parser.add_argument('--non-consortium-out',
                        help='Output CSV for non-consortium sequences')
    parser.add_argument('--consortium-group', default='Pathogen Genomic Lab',
                        help='Group name for consortium sequences (default: Pathogen Genomic Lab)')
    args = parser.parse_args()

    meta = pd.read_csv(args.metadata, sep='\t')
    aliases = pd.read_csv(args.aliases)
    alias_map = build_alias_map(aliases)

    province_names = set()
    if args.hz_lookup:
        hz_df = pd.read_csv(args.hz_lookup)
        province_names = set(hz_df['province'].unique())
        print(f'Loaded {len(province_names)} province names from {args.hz_lookup}')

    print(f'Loaded {len(meta)} accessions from metadata')
    print(f'Loaded {len(aliases)} alias mappings')

    # Split consortium vs non-consortium
    consortium = meta[meta['groupName'] == args.consortium_group].copy()
    non_consortium = meta[meta['groupName'] != args.consortium_group].copy()
    print(f'Consortium ({args.consortium_group}): {len(consortium)}')
    print(f'Non-consortium: {len(non_consortium)} ({non_consortium["groupName"].nunique()} groups)')

    if args.non_consortium_out and not non_consortium.empty:
        cols = ['accessionVersion', 'submissionId', 'groupName',
                'sampleCollectionDate', 'geoLocCountry', 'geoLocAdmin2', 'completeness']
        out_cols = [c for c in cols if c in non_consortium.columns]
        non_consortium[out_cols].to_csv(args.non_consortium_out, index=False)
        print(f'Wrote {len(non_consortium)} non-consortium sequences to {args.non_consortium_out}')

    # Link to linelist
    links, ll_by_sample_id = link_to_linelist(consortium, args.linelist)

    # Run checks on consortium sequences only
    all_issues = []
    dup_issues = check_duplicates(consortium, links, ll_by_sample_id)
    all_issues.extend(dup_issues)

    # Skip accessions handled as duplicates (REVOKE or TEMP_REVOKE) from other checks
    revoked_accs = {i['accessionVersion'] for i in dup_issues
                    if i['suggested_fix'] in ('REVOKE', 'TEMP_REVOKE')}
    not_ll_issues = check_not_in_linelist(consortium, links)
    all_issues.extend([i for i in not_ll_issues if i['accessionVersion'] not in revoked_accs])

    all_issues.extend([i for i in check_sample_ids(consortium, links) if i['accessionVersion'] not in revoked_accs])
    all_issues.extend([i for i in check_linelist_mismatches(consortium, links, alias_map) if i['accessionVersion'] not in revoked_accs])
    all_issues.extend([i for i in check_admin2(consortium, alias_map, province_names) if i['accessionVersion'] not in revoked_accs])
    all_issues.extend([i for i in check_dates(consortium) if i['accessionVersion'] not in revoked_accs])

    if args.fasta:
        quality_issues = check_sequence_quality(consortium, args.fasta)
        all_issues.extend(quality_issues)
        n_low = sum(1 for i in quality_issues if i['issue'] == 'low_completeness')
        n_mm = sum(1 for i in quality_issues if i['issue'] == 'high_non_adar_mismatches')
        print(f'Sequence quality: {n_low} low completeness, {n_mm} high non-ADAR mismatches')

    issues_df = pd.DataFrame(all_issues)

    if not issues_df.empty:
        extra_cols = ['submissionId', 'sampleCollectionDate', 'geoLocAdmin2']
        for col in extra_cols:
            if col in consortium.columns:
                col_map = consortium.set_index('accessionVersion')[col]
                issues_df[col] = issues_df['accessionVersion'].map(col_map)

        issues_df = issues_df.sort_values(['issue', 'accessionVersion'])
        issues_df.to_csv(args.output, index=False)
        print(f'\nWrote {len(issues_df)} issues to {args.output}')
    else:
        pd.DataFrame(columns=['accessionVersion', 'issue', 'field', 'current_value',
                               'suggested_fix', 'details']).to_csv(args.output, index=False)
        print('\nNo issues found.')

    # Report
    txt_path = args.report_txt
    if not txt_path:
        txt_path = str(Path(args.output).with_suffix('.txt'))
    summarize(issues_df, consortium, out_file=txt_path)

    # Write exclude list (duplicates + temp revokes + quality failures)
    if args.exclude_out and not issues_df.empty:
        dup_revoke = issues_df[(issues_df['issue'] == 'duplicate_resolved') &
                               (issues_df['suggested_fix'] == 'REVOKE')]
        temp_revoke = issues_df[(issues_df['issue'] == 'duplicate_unresolved') &
                                (issues_df['suggested_fix'] == 'TEMP_REVOKE')]
        quality_fail = issues_df[issues_df['issue'].isin(
            ['low_completeness', 'high_non_adar_mismatches'])]

        with open(args.exclude_out, 'w') as f:
            f.write(f'# Auto-generated exclude list — {date.today().strftime("%Y-%m-%d")}\n')
            if not dup_revoke.empty:
                f.write(f'# Resolved duplicates ({len(dup_revoke)})\n')
                for _, row in dup_revoke.iterrows():
                    f.write(f'{row["accessionVersion"]}\t# {row["current_value"]}: {row["details"]}\n')
            if not temp_revoke.empty:
                f.write(f'# Unresolved duplicates — temporarily excluded ({len(temp_revoke)})\n')
                for _, row in temp_revoke.iterrows():
                    f.write(f'{row["accessionVersion"]}\t# {row["current_value"]}: {row["details"]}\n')
            if not quality_fail.empty:
                f.write(f'# Sequence quality failures ({len(quality_fail)})\n')
                for _, row in quality_fail.iterrows():
                    f.write(f'{row["accessionVersion"]}\t# {row["issue"]}: {row["details"]}\n')

        total_excluded = len(dup_revoke) + len(temp_revoke) + len(quality_fail)
        print(f'Wrote {total_excluded} accessions to {args.exclude_out} '
              f'({len(dup_revoke)} resolved duplicates, {len(temp_revoke)} temp revokes, {len(quality_fail)} quality)')


if __name__ == '__main__':
    main()
