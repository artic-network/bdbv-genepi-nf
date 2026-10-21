#!/usr/bin/env python3
"""Find the best-fit temporal root of a phylogenetic tree and identify RTT outliers.

Reimplements PearTree's temporalRootGraph algorithm: searches every branch for
the root position that minimises the RMSE of an OLS root-to-tip regression,
enforcing a positive clock rate. O(N) via subtree aggregate sums.

Outputs:
  {prefix}.bestFitRoot.tree         — rerooted Newick
  {prefix}.bestFitRoot.rootToTip.csv — RTT data (PearTree-compatible format)
  {prefix}.outliers.txt             — tip labels with |z| > threshold
  {prefix}.rtt_summary.txt          — clock rate, R², root date, RMSE
  {prefix}.rtt_plot.pdf             — RTT scatter (optional, if matplotlib available)

Usage:
  python temporal_root.py -t tree.treefile
  python temporal_root.py -t tree.treefile -z 3.0 -o output_dir --prefix my_analysis
"""
import argparse
import csv
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import toYearFraction, toDatestring


# ---------------------------------------------------------------------------
# Tree data structure (minimal, to avoid treeswift reroot bugs)
# ---------------------------------------------------------------------------

class Node:
    __slots__ = ('idx', 'label', 'children', 'parent', 'edge_length',
                 'adjacents', 'lengths')

    def __init__(self, idx, label=None, edge_length=None):
        self.idx = idx
        self.label = label
        self.children = []
        self.parent = None
        self.edge_length = edge_length
        self.adjacents = []
        self.lengths = []

    def is_leaf(self):
        return len(self.children) == 0


def parse_newick(path):
    """Parse a Newick file into a list of Node objects. Returns (nodes, root_idx)."""
    with open(path) as f:
        nwk = f.read().strip()
    if nwk.endswith(';'):
        nwk = nwk[:-1]

    nodes = []
    idx_counter = [0]

    def new_node(label=None, edge_length=None):
        n = Node(idx_counter[0], label=label, edge_length=edge_length)
        nodes.append(n)
        idx_counter[0] += 1
        return n

    root = new_node()
    current = root
    i = 0
    while i < len(nwk):
        c = nwk[i]
        if c == '(':
            child = new_node()
            child.parent = current
            current.children.append(child)
            current = child
            i += 1
        elif c == ')':
            current = current.parent
            i += 1
            # read label and/or branch length after closing paren
            label_chars = []
            while i < len(nwk) and nwk[i] not in '(),;:':
                label_chars.append(nwk[i])
                i += 1
            if label_chars:
                current.label = ''.join(label_chars)
        elif c == ',':
            current = current.parent
            child = new_node()
            child.parent = current
            current.children.append(child)
            current = child
            i += 1
        elif c == ':':
            i += 1
            num_chars = []
            while i < len(nwk) and nwk[i] not in '(),;:':
                num_chars.append(nwk[i])
                i += 1
            current.edge_length = float(''.join(num_chars))
        else:
            label_chars = []
            while i < len(nwk) and nwk[i] not in '(),;:':
                label_chars.append(nwk[i])
                i += 1
            current.label = ''.join(label_chars)

    return nodes, root.idx


def get_leaves(nodes, root_idx):
    """Return list of leaf node indices."""
    return [n.idx for n in nodes if n.is_leaf()]


def build_adjacency(nodes, root_idx):
    """Build undirected adjacency lists from the rooted tree."""
    for n in nodes:
        n.adjacents = []
        n.lengths = []
    for n in nodes:
        if n.parent is not None:
            p = n.parent
            el = n.edge_length if n.edge_length is not None else 0.0
            n.adjacents.append(p.idx)
            n.lengths.append(el)
            p.adjacents.append(n.idx)
            p.lengths.append(el)


def root_to_tip_distances(nodes, root_idx):
    """Compute root-to-tip distances via BFS from root."""
    dist = {root_idx: 0.0}
    stack = [root_idx]
    while stack:
        cur = stack.pop()
        n = nodes[cur]
        for child in n.children:
            d = child.edge_length if child.edge_length is not None else 0.0
            dist[child.idx] = dist[cur] + d
            stack.append(child.idx)
    return dist


# ---------------------------------------------------------------------------
# Reroot
# ---------------------------------------------------------------------------

def reroot_on_edge(nodes, root_idx, child_idx, dist_from_parent):
    """Reroot the tree on the edge leading to child_idx, at dist_from_parent
    from the parent of child_idx. Returns (new_nodes, new_root_idx)."""
    import copy

    child_node = nodes[child_idx]
    parent_node = child_node.parent
    full_edge = child_node.edge_length if child_node.edge_length is not None else 0.0

    # Create new root node
    new_root = Node(len(nodes), label=None, edge_length=None)

    # Deep copy all existing nodes
    new_nodes = [None] * (len(nodes) + 1)
    for n in nodes:
        nn = Node(n.idx, label=n.label, edge_length=n.edge_length)
        new_nodes[n.idx] = nn
    new_nodes[new_root.idx] = new_root

    # Rebuild parent/children relationships
    for n in nodes:
        nn = new_nodes[n.idx]
        nn.children = [new_nodes[c.idx] for c in n.children]
        nn.parent = new_nodes[n.parent.idx] if n.parent is not None else None

    # Now modify the tree: insert new root on the edge parent->child
    cn = new_nodes[child_idx]
    pn = cn.parent

    # Remove child from parent's children
    pn.children = [c for c in pn.children if c.idx != child_idx]

    # Set up new root with two children
    cn.edge_length = dist_from_parent
    cn.parent = new_root
    new_root.children.append(cn)

    # Now reverse the path from parent up to old root
    remaining_edge = full_edge - dist_from_parent
    prev = new_root
    cur = pn

    while cur is not None:
        next_up = cur.parent
        # Remove cur from its parent's children
        if next_up is not None:
            next_up.children = [c for c in next_up.children if c.idx != cur.idx]

        # Reverse: cur becomes child of prev
        cur.parent = prev
        prev.children.append(cur)
        cur.edge_length = remaining_edge

        if next_up is not None:
            remaining_edge = cur.edge_length if next_up.parent is not None else cur.edge_length
            # The edge from next_up to cur (which we are reversing)
            # In the original tree, cur.edge_length was the edge from cur to its parent
            # After reversal, we need the edge from next_up to cur
            # which was next_up's edge_length in the original (but cur was already set)
            remaining_edge = next_up.edge_length if next_up.edge_length is not None else 0.0

        prev = cur
        cur = next_up

    # Handle the old root — if it was the top of the path, it needs special handling
    # The old root may have had multiple children; the ones NOT on the reversal path
    # stay as children. If old root has exactly one remaining child after reversal,
    # we can suppress it.

    return new_nodes, new_root.idx


def reroot_simple(nodes, root_idx, child_idx, dist_from_parent):
    """Simpler reroot: write the tree to Newick, use a fresh parse approach.
    Actually, let's do it by rebuilding from adjacency."""

    build_adjacency(nodes, root_idx)
    n_nodes = len(nodes)

    child_node = nodes[child_idx]
    parent_idx = child_node.parent.idx
    full_edge = child_node.edge_length if child_node.edge_length is not None else 0.0

    # Build new tree from adjacency with new root inserted on the target edge
    new_root_idx = n_nodes

    # Create adjacency for new root
    adj = [[] for _ in range(n_nodes + 1)]
    elen = [[] for _ in range(n_nodes + 1)]

    for n in nodes:
        for a, l in zip(n.adjacents, n.lengths):
            if (n.idx == child_idx and a == parent_idx) or \
               (n.idx == parent_idx and a == child_idx):
                continue
            adj[n.idx].append(a)
            elen[n.idx].append(l)

    # Connect new root to child and parent
    dist_to_child = full_edge - dist_from_parent
    adj[new_root_idx].append(child_idx)
    elen[new_root_idx].append(dist_to_child)
    adj[child_idx].append(new_root_idx)
    elen[child_idx].append(dist_to_child)

    adj[new_root_idx].append(parent_idx)
    elen[new_root_idx].append(dist_from_parent)
    adj[parent_idx].append(new_root_idx)
    elen[parent_idx].append(dist_from_parent)

    # BFS from new root to build rooted tree
    new_nodes = []
    for i in range(n_nodes + 1):
        label = nodes[i].label if i < n_nodes else None
        nn = Node(i, label=label)
        new_nodes.append(nn)

    visited = set()
    stack = [new_root_idx]
    visited.add(new_root_idx)
    while stack:
        cur = stack.pop()
        for a, l in zip(adj[cur], elen[cur]):
            if a not in visited:
                visited.add(a)
                new_nodes[a].parent = new_nodes[cur]
                new_nodes[a].edge_length = l
                new_nodes[cur].children.append(new_nodes[a])
                stack.append(a)

    return new_nodes, new_root_idx


# ---------------------------------------------------------------------------
# Newick output
# ---------------------------------------------------------------------------

def to_newick(nodes, root_idx):
    """Convert rooted tree to Newick string (iterative to avoid stack overflow)."""
    result = []
    # Iterative post-order traversal
    stack = [(root_idx, False)]
    while stack:
        idx, processed = stack.pop()
        node = nodes[idx]
        if processed or node.is_leaf():
            if node.is_leaf():
                result.append(node.label if node.label else '')
            else:
                result.append(')')
                if node.label:
                    result.append(node.label)
            if node.edge_length is not None and node.parent is not None:
                result.append(f':{node.edge_length}')
            if stack:
                next_idx, next_proc = stack[-1]
                if not next_proc and not nodes[next_idx].is_leaf():
                    pass
                # add comma if next sibling
            continue

        # Push closing marker, then children in reverse
        stack.append((idx, True))
        for i, child in enumerate(reversed(node.children)):
            if i > 0:
                stack.append((',', None))
            stack.append((child.idx, False))
        result.append('(')

    # This approach is getting complicated. Use simple recursion with increased limit.
    return _newick_recursive(nodes, root_idx) + ';'


def _newick_recursive(nodes, idx):
    node = nodes[idx]
    if node.is_leaf():
        s = node.label if node.label else ''
    else:
        children_str = ','.join(_newick_recursive(nodes, c.idx) for c in node.children)
        s = f'({children_str})'
        if node.label:
            s += node.label
    if node.edge_length is not None and node.parent is not None:
        s += f':{node.edge_length}'
    return s


# ---------------------------------------------------------------------------
# Temporal rooting algorithm (translated from PearTree's phylograph.js)
# ---------------------------------------------------------------------------

def build_root_opt_state(nodes, root_idx, tip_dates):
    """Build subtree aggregate sums for O(1) per-branch evaluation.

    This is a direct translation of PearTree's _buildRootOptState.
    The tree must already be in adjacency form (call build_adjacency first).

    Args:
        nodes: list of Node objects with adjacents/lengths populated
        root_idx: index of current root node
        tip_dates: dict mapping node index -> decimal year date

    Returns:
        dict with all intermediate sums needed by eval_branch
    """
    n_nodes = len(nodes)
    anchor = root_idx

    # DFS to collect traversal order and identify tips
    order = []
    parent_of = [-1] * n_nodes
    tip_indices = []
    tip_pos = {}

    stack = [(anchor, -1)]
    while stack:
        cur, frm = stack.pop()
        order.append(cur)
        parent_of[cur] = frm
        children = [a for a in nodes[cur].adjacents if a != frm]
        if not children and cur in tip_dates:
            tip_pos[cur] = len(tip_indices)
            tip_indices.append(cur)
        elif children:
            for c in reversed(children):
                stack.append((c, cur))

    N = len(tip_indices)
    t = np.array([tip_dates.get(i, 0.0) for i in tip_indices])
    sum_t = t.sum()
    sum_tt = (t * t).sum()
    t_bar = sum_t / N
    C = sum_tt - sum_t * sum_t / N

    # Bottom-up: subtree sums
    sub_n = np.zeros(n_nodes, dtype=np.int64)
    sub_t = np.zeros(n_nodes)
    sub_y = np.zeros(n_nodes)
    sub_yy = np.zeros(n_nodes)
    sub_ty = np.zeros(n_nodes)

    for oi in range(len(order) - 1, -1, -1):
        cur = order[oi]
        frm = parent_of[cur]
        if cur in tip_pos:
            ai = tip_pos[cur]
            sub_n[cur] = 1
            sub_t[cur] = t[ai]
        else:
            for i, adj in enumerate(nodes[cur].adjacents):
                if adj == frm:
                    continue
                e = nodes[cur].lengths[i]
                nc = sub_n[adj]
                sub_n[cur] += nc
                sub_t[cur] += sub_t[adj]
                sub_y[cur] += sub_y[adj] + nc * e
                sub_yy[cur] += sub_yy[adj] + 2 * e * sub_y[adj] + nc * e * e
                sub_ty[cur] += sub_ty[adj] + e * sub_t[adj]

    # Top-down: total-from-node sums
    total_y = np.zeros(n_nodes)
    total_yy = np.zeros(n_nodes)
    total_ty = np.zeros(n_nodes)

    total_y[anchor] = sub_y[anchor]
    total_yy[anchor] = sub_yy[anchor]
    total_ty[anchor] = sub_ty[anchor]

    for oi in range(len(order)):
        p = order[oi]
        for i, c in enumerate(nodes[p].adjacents):
            if c == parent_of[p]:
                continue
            e = nodes[p].lengths[i]
            nc = sub_n[c]
            sub_y_from_p = sub_y[c] + nc * e
            total_y[c] = total_y[p] + (N - 2 * nc) * e
            total_yy[c] = total_yy[p] + 2 * e * (total_y[p] - 2 * sub_y_from_p) + N * e * e
            total_ty[c] = total_ty[p] + e * (sum_t - 2 * sub_t[c])

    return {
        'nodes': nodes, 'anchor': anchor, 'order': order, 'parent_of': parent_of,
        'N': N, 'Nd': N, 't': t, 't_bar': t_bar, 'C': C, 'sum_t': sum_t,
        'sub_n': sub_n, 'sub_t': sub_t, 'sub_y': sub_y, 'sub_yy': sub_yy, 'sub_ty': sub_ty,
        'total_y': total_y, 'total_yy': total_yy, 'total_ty': total_ty,
        'tip_indices': tip_indices, 'tip_pos': tip_pos,
    }


def eval_branch(child_idx, parent_idx, state):
    """Evaluate the optimal root position on a single branch.

    Direct translation of PearTree's _evalBranch with forcePositiveRate=True.

    Returns:
        (child_idx, dist_from_parent, score) or None if branch is invalid
    """
    nodes = state['nodes']
    N = state['N']
    Nd = state['Nd']
    sum_t = state['sum_t']
    t_bar = state['t_bar']
    C = state['C']
    sub_n = state['sub_n']
    sub_t = state['sub_t']
    sub_y = state['sub_y']
    sub_yy = state['sub_yy']
    sub_ty = state['sub_ty']
    total_y = state['total_y']
    total_yy = state['total_yy']
    total_ty = state['total_ty']

    # Find edge length from child to parent in adjacency
    child_node = nodes[child_idx]
    L = None
    for a, l in zip(child_node.adjacents, child_node.lengths):
        if a == parent_idx:
            L = l
            break
    if L is None or L <= 0:
        return None

    nd = sub_n[child_idx]
    np_ = Nd - nd
    if nd == 0 or np_ == 0:
        return None

    # B-side sums (below child, measured from child)
    sum_rB = sub_y[child_idx]
    sum_rrB = sub_yy[child_idx]
    sum_trB = sub_ty[child_idx]
    sum_tB = sub_t[child_idx]

    # P-side sums
    sum_tP = sum_t - sum_tB
    sum_rP = total_y[child_idx] - sum_rB - np_ * L
    sum_rrP = total_yy[child_idx] - sum_rrB - 2 * L * sum_rP - np_ * L * L
    sum_trP = total_ty[child_idx] - sum_trB - L * sum_tP

    # Heights at d=0 (root at parent)
    sum_hB0 = sum_rB + nd * L
    sum_dBL2 = sum_rrB + 2 * L * sum_rB + nd * L * L
    sum_dP2 = sum_rrP
    M0 = (sum_rP + sum_hB0) / Nd

    # Polynomial coefficients for N*Var(heights)
    B2 = 4 * nd * np_ / Nd
    if B2 * L <= 1e-20:
        return None
    B0 = (sum_dP2 - 2 * M0 * sum_rP + np_ * M0 * M0) + \
         (sum_dBL2 - 2 * M0 * sum_hB0 + nd * M0 * M0)
    sumV_B = sum_hB0 - nd * M0
    sumV_P = sum_rP - np_ * M0
    B1 = 2 * (sumV_P - sumV_B)

    if C > 1e-20:
        # Heterochronous: minimise OLS regression residual
        sum_yAdj = sum_rP + sum_hB0
        sum_ty0 = sum_trP + sum_trB + L * sum_tB
        A0 = sum_ty0 - t_bar * sum_yAdj
        A1 = (sum_tP - np_ * t_bar) - (sum_tB - nd * t_bar)

        # Enforce positive rate
        d_lo = 0.0
        d_hi = L
        if A1 > 1e-20:
            d_lo = max(d_lo, -A0 / A1)
        elif A1 < -1e-20:
            d_hi = min(d_hi, -A0 / A1)
        elif A0 <= 0:
            return None
        if d_lo >= d_hi:
            return None

        Q2 = B2 - A1 * A1 / C
        Q1 = B1 - 2 * A0 * A1 / C

        if Q2 > 1e-30:
            d = max(d_lo, min(d_hi, -Q1 / (2 * Q2)))
        else:
            flo = B0 + B1 * d_lo + B2 * d_lo * d_lo - (A0 + A1 * d_lo) ** 2 / C
            fhi = B0 + B1 * d_hi + B2 * d_hi * d_hi - (A0 + A1 * d_hi) ** 2 / C
            d = d_lo if flo <= fhi else d_hi

        ssxy_new = A0 + A1 * d
        ssyy_new = B0 + B1 * d + B2 * d * d
        score = (ssyy_new - ssxy_new * ssxy_new / C) / Nd
    else:
        # Homochronous: minimise variance
        d = max(0, min(L, -B1 / (2 * B2)))
        score = (B0 + B1 * d + B2 * d * d) / Nd

    return (child_idx, d, score)


def find_temporal_root(nodes, root_idx, tip_dates):
    """Find the best temporal root across all branches.

    Returns:
        (child_idx, dist_from_parent) for the optimal root position
    """
    build_adjacency(nodes, root_idx)
    state = build_root_opt_state(nodes, root_idx, tip_dates)
    anchor = state['anchor']
    parent_of = state['parent_of']

    best_score = float('inf')
    best_child = None
    best_dist = 0.0

    # Evaluate all non-root edges
    for node in nodes:
        p = parent_of[node.idx]
        if p == -1:
            continue
        result = eval_branch(node.idx, p, state)
        if result and result[2] < best_score:
            best_score = result[2]
            best_child = result[0]
            best_dist = result[1]

    if best_child is None:
        # Fallback: midpoint of first edge
        for node in nodes:
            if node.parent is not None and node.edge_length:
                return (node.idx, node.edge_length / 2)
        return (nodes[1].idx if len(nodes) > 1 else 0, 0.0)

    return (best_child, best_dist)


# ---------------------------------------------------------------------------
# OLS regression
# ---------------------------------------------------------------------------

def ols_regression(x, y):
    """Simple OLS: y = slope * x + intercept. Returns dict with stats."""
    n = len(x)
    x_bar = x.mean()
    y_bar = y.mean()
    ssxy = ((x - x_bar) * (y - y_bar)).sum()
    ssxx = ((x - x_bar) ** 2).sum()
    ssyy = ((y - y_bar) ** 2).sum()

    slope = ssxy / ssxx if ssxx > 0 else 0.0
    intercept = y_bar - slope * x_bar
    predicted = slope * x + intercept
    residuals = y - predicted
    ss_res = (residuals ** 2).sum()
    r_squared = 1.0 - ss_res / ssyy if ssyy > 0 else 0.0
    rmse = np.sqrt(ss_res / n) if n > 0 else 0.0
    x_intercept = -intercept / slope if abs(slope) > 1e-30 else float('nan')

    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_squared,
        'rmse': rmse,
        'x_intercept': x_intercept,
        'predicted': predicted,
        'residuals': residuals,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def extract_date(label, delimiter='|', field=-1):
    """Extract date string from a tip label."""
    return label.split(delimiter)[field]


def _run_single_pass(tree_path, tip_dates_by_label, delimiter, date_field,
                     exclude_labels=None):
    """Parse tree, find temporal root (excluding labels), reroot, compute
    regression on ALL tips. Returns (new_nodes, new_root_idx, names,
    dates_arr, divs_arr, reg, zscores, child_idx, dist)."""
    exclude_labels = exclude_labels or set()

    nodes, root_idx = parse_newick(tree_path)
    leaves = get_leaves(nodes, root_idx)

    all_tip_dates = {}
    filtered_tip_dates = {}
    for li in leaves:
        label = nodes[li].label
        if not label or label not in tip_dates_by_label:
            continue
        all_tip_dates[li] = tip_dates_by_label[label]
        if label not in exclude_labels:
            filtered_tip_dates[li] = tip_dates_by_label[label]

    child_idx, dist = find_temporal_root(nodes, root_idx, filtered_tip_dates)
    new_nodes, new_root_idx = reroot_simple(nodes, root_idx, child_idx, dist)
    dists = root_to_tip_distances(new_nodes, new_root_idx)

    dates_arr, divs_arr, names = [], [], []
    for n in new_nodes:
        if not (n.is_leaf() and n.label and n.label in tip_dates_by_label):
            continue
        dates_arr.append(tip_dates_by_label[n.label])
        divs_arr.append(dists[n.idx])
        names.append(n.label)

    dates_arr = np.array(dates_arr)
    divs_arr = np.array(divs_arr)
    reg = ols_regression(dates_arr, divs_arr)
    zscores = (reg['residuals'] / reg['rmse']
               if reg['rmse'] > 0 else np.zeros_like(reg['residuals']))

    return (new_nodes, new_root_idx, names, dates_arr, divs_arr,
            reg, zscores, child_idx, dist)


def iterative_temporal_root(tree_path, tip_dates_by_label, delimiter, date_field,
                            zscore_threshold=2.0, max_iter=10):
    """Iteratively refine the temporal root by excluding outliers from the
    rooting regression, then re-evaluating all tips."""
    exclude = set()
    prev_outliers = None

    for iteration in range(1, max_iter + 1):
        result = _run_single_pass(tree_path, tip_dates_by_label, delimiter,
                                  date_field, exclude_labels=exclude)
        (new_nodes, new_root_idx, names, dates_arr, divs_arr,
         reg, zscores, child_idx, dist) = result

        outlier_labels = {names[i] for i in range(len(names))
                          if abs(zscores[i]) > zscore_threshold}
        n_out = len(outlier_labels)
        root_dec = reg['x_intercept']
        root_str = toDatestring(root_dec) if not np.isnan(root_dec) else 'N/A'

        print(f"  Iteration {iteration}: root edge={child_idx}, "
              f"d={dist:.8e}, root date={root_str}, "
              f"outliers={n_out}, excluded={len(exclude)}")

        if outlier_labels == prev_outliers:
            print(f"  Converged after {iteration} iterations.")
            break

        prev_outliers = outlier_labels
        exclude = outlier_labels
    else:
        print(f"  Did not converge after {max_iter} iterations.")

    return (new_nodes, new_root_idx, names, dates_arr, divs_arr,
            reg, zscores, iteration)


def main():
    parser = argparse.ArgumentParser(
        description='Find best-fit temporal root and identify RTT outliers')
    parser.add_argument('-t', '--tree', required=True, help='Newick tree file')
    parser.add_argument('-d', '--delimiter', default='|',
                        help='Character separating fields in tip labels (default: |)')
    parser.add_argument('-f', '--date-field', type=int, default=-1,
                        help='0-indexed position of the date field in tip labels, '
                             'or negative to count from end (default: -1, i.e. last field)')
    parser.add_argument('-z', '--zscore-threshold', type=float, default=2.0,
                        help='Z-score threshold for outlier flagging (default: 2.0)')
    parser.add_argument('-o', '--outdir', default='.', help='Output directory')
    parser.add_argument('--prefix', default=None,
                        help='Output prefix (default: tree filename)')
    parser.add_argument('--iterate', action='store_true',
                        help='Iteratively refine root by excluding outliers '
                             'from the rooting regression')
    parser.add_argument('--max-iter', type=int, default=10,
                        help='Maximum iterations for --iterate (default: 10)')
    parser.add_argument('--plot', action='store_true',
                        help='Generate RTT scatter plot')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.prefix or os.path.splitext(os.path.basename(args.tree))[0]
    out_prefix = os.path.join(args.outdir, prefix)

    # Parse tree
    sys.setrecursionlimit(10000)
    nodes, root_idx = parse_newick(args.tree)
    leaves = get_leaves(nodes, root_idx)
    print(f"Loaded tree: {len(nodes)} nodes, {len(leaves)} tips")

    # Extract tip dates
    tip_dates = {}
    for li in leaves:
        label = nodes[li].label
        if not label:
            print(f"WARNING: tip {li} has no label, skipping")
            continue
        try:
            date_str = extract_date(label, args.delimiter, args.date_field)
            tip_dates[li] = toYearFraction(date_str)
        except Exception as e:
            print(f"WARNING: cannot parse date for {label}: {e}")

    if len(tip_dates) < 3:
        print("ERROR: need at least 3 tips with dates")
        sys.exit(1)

    print(f"Parsed dates for {len(tip_dates)}/{len(leaves)} tips")
    date_range = max(tip_dates.values()) - min(tip_dates.values())
    print(f"Date range: {toDatestring(min(tip_dates.values()))} to "
          f"{toDatestring(max(tip_dates.values()))} ({date_range:.4f} years)")

    # Build label -> date mapping for all tips
    tip_dates_by_label = {}
    for li in tip_dates:
        tip_dates_by_label[nodes[li].label] = tip_dates[li]

    if args.iterate:
        print(f"\nIterative temporal rooting (max {args.max_iter} iterations)...")
        (new_nodes, new_root_idx, names, dates_arr, divs_arr,
         reg, zscores, n_iter) = iterative_temporal_root(
            args.tree, tip_dates_by_label, args.delimiter, args.date_field,
            args.zscore_threshold, args.max_iter)
    else:
        print("\nSearching for optimal temporal root...")
        (new_nodes, new_root_idx, names, dates_arr, divs_arr,
         reg, zscores, _, _) = _run_single_pass(
            args.tree, tip_dates_by_label, args.delimiter, args.date_field)
        n_iter = 1

    # Write rerooted tree
    tree_out = f"{out_prefix}.bestFitRoot.tree"
    nwk = _newick_recursive(new_nodes, new_root_idx) + ';'
    with open(tree_out, 'w') as f:
        f.write(nwk + '\n')
    print(f"Wrote rerooted tree: {tree_out}")

    # Summary
    root_date_dec = reg['x_intercept']
    root_date_str = toDatestring(root_date_dec) if not np.isnan(root_date_dec) else 'N/A'
    clock_rate = reg['slope']
    n_outliers = int(np.sum(np.abs(zscores) > args.zscore_threshold))

    print(f"\nRegression results:")
    print(f"  Clock rate:  {clock_rate:.6e} subs/site/year")
    print(f"  R²:          {reg['r_squared']:.6f}")
    print(f"  RMSE:        {reg['rmse']:.6e}")
    print(f"  Root date:   {root_date_str} ({root_date_dec:.4f})")
    print(f"  Outliers:    {n_outliers} (|z| > {args.zscore_threshold})")

    # Write RTT CSV
    csv_path = f"{out_prefix}.bestFitRoot.rootToTip.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['name', 'date', 'date (decimal year)', 'divergence',
                         'regression', 'residual', 'residual z-score'])
        for i in range(len(names)):
            writer.writerow([
                names[i],
                extract_date(names[i], args.delimiter, args.date_field),
                dates_arr[i],
                divs_arr[i],
                reg['predicted'][i],
                reg['residuals'][i],
                zscores[i],
            ])
    print(f"Wrote RTT CSV: {csv_path}")

    # Write outliers
    outlier_path = f"{out_prefix}.outliers.txt"
    outlier_names = []
    with open(outlier_path, 'w') as f:
        f.write(f"# Temporal outliers (|z-score| > {args.zscore_threshold})\n")
        for i in range(len(names)):
            if abs(zscores[i]) > args.zscore_threshold:
                f.write(f"{names[i]}\tz={zscores[i]:+.2f}\n")
                outlier_names.append((names[i], zscores[i]))
    print(f"Wrote outliers: {outlier_path} ({n_outliers} sequences)")
    for name, z in sorted(outlier_names, key=lambda x: -abs(x[1])):
        sign = '+' if z > 0 else '-'
        print(f"  {name.split('|')[0]}  z={z:+.2f}  ({sign}ve)")

    # Write summary
    summary_path = f"{out_prefix}.rtt_summary.txt"
    with open(summary_path, 'w') as f:
        f.write(f"Clock rate: {clock_rate:.6e} subs/site/year\n")
        f.write(f"R-squared: {reg['r_squared']:.6f}\n")
        f.write(f"RMSE: {reg['rmse']:.6e}\n")
        f.write(f"Root date: {root_date_str} ({root_date_dec:.4f})\n")
        f.write(f"Tips: {len(names)}\n")
        f.write(f"Outliers (|z| > {args.zscore_threshold}): {n_outliers}\n")
        if args.iterate:
            f.write(f"Iterative refinement: {n_iter} iterations\n")
    print(f"Wrote summary: {summary_path}")

    # Optional plot
    if args.plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import matplotlib.dates as mdates
            from matplotlib.colors import to_rgba
            from datetime import datetime, timedelta

            def _fracyear_to_mdate(fy):
                year = int(fy)
                rem = fy - year
                days = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
                dt = datetime(year, 1, 1) + timedelta(days=rem * days)
                return mdates.date2num(dt)

            dates_md = np.array([_fracyear_to_mdate(d) for d in dates_arr])
            outlier_mask = np.abs(zscores) > args.zscore_threshold

            fig, ax = plt.subplots(1, 1, figsize=(7, 4.5), facecolor='white')
            ax.set_facecolor('white')

            tip_color = '#6699BB'
            outlier_color = '#CC6655'

            ax.scatter(dates_md[~outlier_mask], divs_arr[~outlier_mask],
                       s=12, color=to_rgba(tip_color, 0.5),
                       edgecolors=to_rgba(tip_color, 0.8), linewidths=0.3,
                       zorder=5, label='Tips', clip_on=False)
            ax.scatter(dates_md[outlier_mask], divs_arr[outlier_mask],
                       s=18, color=to_rgba(outlier_color, 0.55),
                       edgecolors=to_rgba(outlier_color, 0.9), linewidths=0.3,
                       zorder=6, label=f'Outliers (|z| > {args.zscore_threshold})',
                       clip_on=False)

            x_pad_fy = (dates_arr.max() - dates_arr.min()) * 0.03
            x_line_fy = np.linspace(dates_arr.min() - x_pad_fy,
                                    dates_arr.max() + x_pad_fy, 100)
            x_line = np.array([_fracyear_to_mdate(d) for d in x_line_fy])
            y_line = reg['slope'] * x_line_fy + reg['intercept']
            band = args.zscore_threshold * reg['rmse']
            ax.plot(x_line, y_line, '-', color='0.15', linewidth=1.2, zorder=4)
            ax.plot(x_line, y_line + band, '--', color='0.45', linewidth=0.8,
                    zorder=3)
            ax.plot(x_line, y_line - band, '--', color='0.45', linewidth=0.8,
                    zorder=3)

            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
            ax.set_xlabel('Sampling date', fontsize=10)
            ax.set_ylabel('Root-to-tip divergence', fontsize=10)
            ax.tick_params(axis='both', labelsize=9, direction='out',
                           length=3, width=0.6)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(0.6)
            ax.spines['bottom'].set_linewidth(0.6)
            ax.yaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f'{v:.1e}' if v != 0 else '0'))

            res_mean_sq = np.mean(reg['residuals'] ** 2)
            cv = reg['rmse'] / np.mean(divs_arr) if np.mean(divs_arr) > 0 else 0.0
            stats_text = (
                f'n = {len(names)}\n'
                f'Rate = {clock_rate:.2e} s/s/y\n'
                f'R² = {reg["r_squared"]:.4f}\n'
                f'Root = {root_date_str}\n'
                f'Res. mean sq. = {res_mean_sq:.2e}\n'
                f'CV = {cv:.4f}'
            )
            ax.text(0.02, 0.97, stats_text, transform=ax.transAxes,
                    fontsize=8, fontfamily='monospace', verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                              edgecolor='0.7', linewidth=0.5, alpha=0.9))

            ax.legend(fontsize=8, frameon=True, loc='lower right',
                      markerscale=1.2, handletextpad=0.4,
                      fancybox=True, edgecolor='0.7',
                      facecolor='white', framealpha=0.9)

            plot_path = f"{out_prefix}.rtt_plot.pdf"
            plt.tight_layout()
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.savefig(plot_path.replace('.pdf', '.png'), dpi=300,
                        bbox_inches='tight')
            print(f"Wrote plot: {plot_path}")
        except ImportError:
            print("matplotlib not available, skipping plot")


if __name__ == '__main__':
    main()
