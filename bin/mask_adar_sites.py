#!/usr/bin/env python3
import argparse
import csv
import treeswift
from collections import defaultdict
from pathlib import Path


GP_START_DEFAULT = 5875
GP_END_DEFAULT = 8279


def parse_alignment(fasta_path):
    seqs = {}
    current = None
    chunks = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current is not None:
                    seqs[current] = list(''.join(chunks))
                current = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        seqs[current] = list(''.join(chunks))
    return seqs


def write_alignment(seqs, path):
    with open(path, 'w') as f:
        for header, seq_list in seqs.items():
            f.write(f'>{header}\n')
            f.write(''.join(seq_list) + '\n')


def parse_state_file(state_path):
    states = {}
    with open(state_path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('Node\t'):
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                node, site, state = parts[0], int(parts[1]), parts[2]
                states[(node, site)] = state
    return states


def get_state(node_name, site, asr_states, alignment_seqs):
    key = (node_name, site)
    if key in asr_states:
        return asr_states[key]
    for header, seq in alignment_seqs.items():
        if header == node_name or header.startswith(node_name + '|'):
            idx = site - 1
            if 0 <= idx < len(seq):
                return seq[idx]
    return '?'


def parse_raccoon_csv(csv_path):
    """Return list of (site, parent, child) for ADAR entries."""
    entries = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['mutation_type'] != 'adar':
                continue
            site = int(row['site'])
            for branch_str in row['present_in'].split(';'):
                parent, child = branch_str.strip().split('_', 1)
                entries.append((site, parent, child))
    return entries


def resolve_internal_nodes(tree_path, child_nodes_set):
    """Traverse tree once, return {internal_node_label: [leaf_labels]} for
    internal nodes in child_nodes_set."""
    tree = treeswift.read_tree_newick(tree_path)
    internal_to_leaves = {}
    for n in tree.traverse_preorder():
        if n.label in child_nodes_set and not n.is_leaf():
            leaves = [leaf.label for leaf in n.traverse_leaves()]
            internal_to_leaves[n.label] = leaves
    return internal_to_leaves


def classify_substitution(parent_state, child_state):
    """Return substitution type string, or None if not an editing signature."""
    p, c = parent_state.upper(), child_state.upper()
    if p == 'T' and c == 'C':
        return 'T>C'
    if p == 'A' and c == 'G':
        return 'A>G'
    if p == 'G' and c == 'A':
        return 'G>A'
    return None


def get_branch_label(child_name):
    """Short label for a branch child (first field before |)."""
    return child_name.split('|')[0]


def mode_mask(alignment_seqs, branch_site_subs, gp_start, gp_end):
    """Apply masking: GP whole-column N for T>C sites, non-GP per-branch
    keep-lowest mask-rest for T>C sites. Returns (modified_seqs, report_lines)."""
    headers = list(alignment_seqs.keys())
    report = []

    # Collect all T>C sites in GP for whole-column masking
    gp_columns_to_mask = set()
    for branch_child, site_subs in branch_site_subs.items():
        for site, sub_type in site_subs:
            if sub_type == 'T>C' and gp_start <= site <= gp_end:
                gp_columns_to_mask.add(site)

    # Mask GP columns
    if gp_columns_to_mask:
        report.append(f"GP region ({gp_start}-{gp_end}): masking {len(gp_columns_to_mask)} columns across all genomes")
        for site in sorted(gp_columns_to_mask):
            idx = site - 1
            for header in headers:
                alignment_seqs[header][idx] = 'N'
            report.append(f"  site {site}: entire column -> N")

    # Non-GP: per branch, keep lowest T>C site, mask rest in affected genomes
    non_gp_masked = 0
    for branch_child, site_subs in branch_site_subs.items():
        tc_sites = sorted([s for s, st in site_subs if st == 'T>C' and not (gp_start <= s <= gp_end)])
        if len(tc_sites) <= 1:
            if tc_sites:
                report.append(f"  {get_branch_label(branch_child)}: 1 non-GP T>C site ({tc_sites[0]}), kept")
            continue

        kept = tc_sites[0]
        to_mask = tc_sites[1:]
        genomes = branch_site_subs[branch_child]  # this is site_subs, need genomes
        report.append(f"  {get_branch_label(branch_child)}: {len(tc_sites)} non-GP T>C sites, kept {kept}, masking {to_mask}")

        for site in to_mask:
            idx = site - 1
            for genome in branch_child if isinstance(branch_child, list) else [branch_child]:
                for header in headers:
                    if header == genome or header.startswith(genome.split('|')[0] + '|'):
                        alignment_seqs[header][idx] = 'N'
                        non_gp_masked += 1
                        break

    report.append(f"\nTotal: {len(gp_columns_to_mask)} GP columns masked, {non_gp_masked} non-GP genome-site pairs masked")
    return alignment_seqs, report


def main():
    parser = argparse.ArgumentParser(description="Handle ADAR-edited sites: mask or exclude")
    parser.add_argument("-a", "--alignment", required=True, help="Input aligned FASTA")
    parser.add_argument("-r", "--raccoon", required=True, help="Raccoon phylo_flags CSV")
    parser.add_argument("-t", "--tree", required=True, help="Newick tree file")
    parser.add_argument("-s", "--state", required=True, help="IQ-TREE .state file")
    parser.add_argument("--mode", required=True, choices=['mask', 'exclude_tc', 'exclude_all'],
                        help="mask: site masking; exclude_tc: remove T>C genomes; exclude_all: remove T>C/A>G/G>A genomes")
    parser.add_argument("-p", "--prefix", required=True, help="Output prefix")
    parser.add_argument("--gp-start", type=int, default=GP_START_DEFAULT,
                        help=f"GP region start in alignment coords (default: {GP_START_DEFAULT})")
    parser.add_argument("--gp-end", type=int, default=GP_END_DEFAULT,
                        help=f"GP region end in alignment coords (default: {GP_END_DEFAULT})")
    args = parser.parse_args()

    alignment_seqs = parse_alignment(args.alignment)
    asr_states = parse_state_file(args.state)
    raccoon_entries = parse_raccoon_csv(args.raccoon)

    # Identify internal nodes in raccoon output
    child_nodes_in_raccoon = set()
    for site, parent, child in raccoon_entries:
        child_nodes_in_raccoon.add(child)

    internal_to_leaves = resolve_internal_nodes(args.tree, child_nodes_in_raccoon)

    # Build per-branch-child → [(site, sub_type)] and per-branch-child → [genome_headers]
    # "branch child" here means the resolved tip genomes, grouped by the original branch
    #
    # We need two structures:
    #   branch_site_subs: {branch_id: [(site, sub_type), ...]}
    #   branch_genomes:   {branch_id: [genome_header, ...]}
    # where branch_id is the original child node label from raccoon

    branch_site_subs = defaultdict(list)
    branch_genomes = {}

    for site, parent, child in raccoon_entries:
        parent_state = get_state(parent, site, asr_states, alignment_seqs)
        child_state = get_state(child, site, asr_states, alignment_seqs)
        sub_type = classify_substitution(parent_state, child_state)

        if child in internal_to_leaves:
            if child not in branch_genomes:
                branch_genomes[child] = internal_to_leaves[child]
            # Use parent→internal-node substitution (from ASR); fall back to
            # checking a leaf if the internal-node state is ambiguous
            if sub_type is not None:
                branch_site_subs[child].append((site, sub_type))
            else:
                for leaf in internal_to_leaves[child]:
                    leaf_state = get_state(leaf, site, asr_states, alignment_seqs)
                    leaf_sub = classify_substitution(parent_state, leaf_state)
                    if leaf_sub is not None:
                        branch_site_subs[child].append((site, leaf_sub))
                        break
        else:
            if child not in branch_genomes:
                matching = [h for h in alignment_seqs if h == child or h.startswith(child.split('|')[0] + '|')]
                branch_genomes[child] = matching
            if sub_type is not None:
                branch_site_subs[child].append((site, sub_type))

    # Deduplicate sites per branch
    for branch in branch_site_subs:
        branch_site_subs[branch] = sorted(set(branch_site_subs[branch]), key=lambda x: x[0])

    report_lines = []
    report_lines.append(f"ADAR handling report")
    report_lines.append(f"Mode: {args.mode}")
    report_lines.append(f"Input: {len(alignment_seqs)} sequences")
    report_lines.append(f"Raccoon entries: {len(raccoon_entries)} (site, branch) pairs")
    report_lines.append(f"Branches with editing signatures: {len(branch_site_subs)}")
    report_lines.append(f"GP region: {args.gp_start}-{args.gp_end}")
    report_lines.append("")

    # Per-branch summary
    report_lines.append("Per-branch details:")
    for branch, site_subs in sorted(branch_site_subs.items(), key=lambda x: -len(x[1])):
        label = get_branch_label(branch)
        genomes = branch_genomes.get(branch, [])
        genome_labels = [get_branch_label(g) for g in genomes]
        tc = [(s, st) for s, st in site_subs if st == 'T>C']
        ag = [(s, st) for s, st in site_subs if st == 'A>G']
        ga = [(s, st) for s, st in site_subs if st == 'G>A']
        parts = []
        if tc:
            parts.append(f"{len(tc)} T>C")
        if ag:
            parts.append(f"{len(ag)} A>G")
        if ga:
            parts.append(f"{len(ga)} G>A")
        sub_summary = ', '.join(parts) if parts else 'none confirmed'
        genome_str = ', '.join(genome_labels) if len(genome_labels) <= 5 else f"{len(genome_labels)} genomes"
        sites_str = ', '.join(f"{s}({st})" for s, st in site_subs)
        report_lines.append(f"  {label} ({sub_summary}): {sites_str} -> {genome_str}")
    report_lines.append("")

    if args.mode == 'mask':
        # GP: whole column N for T>C sites
        gp_tc_sites = set()
        for branch, site_subs in branch_site_subs.items():
            for site, sub_type in site_subs:
                if sub_type == 'T>C' and args.gp_start <= site <= args.gp_end:
                    gp_tc_sites.add(site)

        headers = list(alignment_seqs.keys())
        if gp_tc_sites:
            report_lines.append(f"GP masking: {len(gp_tc_sites)} columns -> N (all genomes)")
            for site in sorted(gp_tc_sites):
                idx = site - 1
                for header in headers:
                    alignment_seqs[header][idx] = 'N'
                report_lines.append(f"  site {site}: entire column -> N")
        else:
            report_lines.append("GP masking: no T>C sites in GP region")
        report_lines.append("")

        # Non-GP: per branch, keep lowest T>C, mask rest in affected genomes only
        non_gp_masked_count = 0
        report_lines.append("Non-GP masking (per branch, keep lowest T>C site):")
        for branch, site_subs in sorted(branch_site_subs.items(), key=lambda x: get_branch_label(x[0])):
            tc_non_gp = sorted([s for s, st in site_subs if st == 'T>C' and not (args.gp_start <= s <= args.gp_end)])
            if not tc_non_gp:
                continue
            genomes = branch_genomes.get(branch, [])
            if len(tc_non_gp) == 1:
                report_lines.append(f"  {get_branch_label(branch)}: 1 non-GP T>C site ({tc_non_gp[0]}), kept")
                continue

            kept = tc_non_gp[0]
            to_mask = tc_non_gp[1:]
            report_lines.append(f"  {get_branch_label(branch)}: {len(tc_non_gp)} non-GP T>C sites, kept {kept}, masking {to_mask}")
            for site in to_mask:
                idx = site - 1
                for genome_header in genomes:
                    for header in headers:
                        if header == genome_header or header.startswith(genome_header.split('|')[0] + '|'):
                            alignment_seqs[header][idx] = 'N'
                            non_gp_masked_count += 1
                            break

        report_lines.append(f"\nSummary: {len(gp_tc_sites)} GP columns masked, {non_gp_masked_count} non-GP genome-site pairs masked")
        report_lines.append(f"Output: {len(alignment_seqs)} sequences (all retained)")

    elif args.mode in ('exclude_tc', 'exclude_all'):
        # Build genome → sites mapping (propagates branch sites to descendant genomes)
        genome_site_subs = defaultdict(list)
        exclude_genomes = set()
        for branch, site_subs in branch_site_subs.items():
            dominated_types = set(st for _, st in site_subs)
            if args.mode == 'exclude_tc':
                should_exclude = 'T>C' in dominated_types
            else:
                should_exclude = bool(dominated_types & {'T>C', 'A>G', 'G>A'})

            if should_exclude:
                for genome in branch_genomes.get(branch, []):
                    exclude_genomes.add(genome)
                    genome_site_subs[genome].extend(site_subs)

        # Deduplicate and sort sites per genome
        for genome in genome_site_subs:
            genome_site_subs[genome] = sorted(set(genome_site_subs[genome]), key=lambda x: x[0])

        report_lines.append(f"Excluding {len(exclude_genomes)} genomes:")
        for g in sorted(exclude_genomes, key=get_branch_label):
            subs = genome_site_subs[g]
            sub_str = ', '.join(f"{s}({st})" for s, st in subs)
            report_lines.append(f"  {get_branch_label(g)} -- {len(subs)} sites: {sub_str}")

        # Remove excluded genomes
        filtered = {}
        for header, seq in alignment_seqs.items():
            if header not in exclude_genomes:
                filtered[header] = seq
        alignment_seqs = filtered

        report_lines.append(f"\nOutput: {len(alignment_seqs)} sequences")

        # Write exclusion list
        excl_path = f"{args.prefix}.adar_exclusions.txt"
        with open(excl_path, 'w') as f:
            f.write("# ADAR-excluded accessions\n")
            f.write(f"# Mode: {args.mode}\n")
            for genome in sorted(exclude_genomes, key=get_branch_label):
                label = get_branch_label(genome)
                acc = genome.split('|')[1].split('.')[0] if '|' in genome else label
                subs = genome_site_subs[genome]
                sub_str = ', '.join(f"{s}({st})" for s, st in subs)
                f.write(f"{acc}  # {label} -- {len(subs)} sites: {sub_str}\n")
        print(f"Wrote exclusion list: {excl_path}")

    # Write output alignment
    out_aln = f"{args.prefix}.aln"
    write_alignment(alignment_seqs, out_aln)
    print(f"Wrote alignment: {out_aln} ({len(alignment_seqs)} sequences)")

    # Write report
    report_path = f"{args.prefix}.adar_report.txt"
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines) + '\n')
    print(f"Wrote report: {report_path}")


if __name__ == '__main__':
    main()
