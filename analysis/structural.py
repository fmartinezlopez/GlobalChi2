"""
Build the structural bin-combination matrix A from a binning configuration file.
Decomposes variables into shared (used by 2+ blocks) and private (used by 1 block),
enumerates only over shared variables, and expands private-variable bin choices.

Supports selection indices: blocks with different selection indices represent
mutually exclusive event populations and are treated independently.

Usage:
    python structural.py <config_file>
"""

import numpy as np
from itertools import product
import time


# ── Parsing bin config ─────────────────────────────────────────────────────

def parse_binning_config(filename):
    """Parse binning configuration from a text file.

    Supports optional 'sel <idx>' at end of block header to specify
    mutually exclusive event selections.  Default selection is 0.
    """
    blocks = []

    with open(filename, "r") as f:
        lines = [
            l.strip()
            for l in f
            if l.strip() and not l.strip().startswith("#")
        ]

    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if parts[0] != "block":
            raise ValueError(f"Expected 'block' keyword, got: {lines[i]}")

        block_idx = int(parts[1])

        # Check for 'sel <idx>' at end of header
        sel_idx = 0
        if len(parts) >= 2 and "sel" in parts:
            sel_pos = parts.index("sel")
            sel_idx = int(parts[sel_pos + 1])
            # Remove sel and its argument for dimension detection
            parts = parts[:sel_pos]

        n_var_parts = len(parts) - 2  # subtract 'block' and block_idx

        if n_var_parts == 1:
            # 1D block
            var_idx = int(parts[2])
            i += 1
            edges = list(map(float, lines[i].split()))
            blocks.append({
                "idx": block_idx,
                "dim": 1,
                "var": var_idx,
                "sel": sel_idx,
                "edges": edges,
                "n_bins": len(edges) - 1,
            })
            i += 1

        elif n_var_parts == 2:
            # 2D block
            slice_var_idx = int(parts[2])
            bin_var_idx = int(parts[3])
            i += 1
            slices = []
            while i < len(lines) and not lines[i].startswith("block"):
                sp = lines[i].split()
                s_idx = int(sp[0])
                s_low = float(sp[1])
                s_high = float(sp[2])
                i += 1
                bin_edges = list(map(float, lines[i].split()))
                slices.append({
                    "idx": s_idx,
                    "low": s_low,
                    "high": s_high,
                    "bin_edges": bin_edges,
                    "n_bins": len(bin_edges) - 1,
                })
                i += 1

            all_slice_edges = set()
            for s in slices:
                all_slice_edges.add(s["low"])
                all_slice_edges.add(s["high"])

            blocks.append({
                "idx": block_idx,
                "dim": 2,
                "slice_var": slice_var_idx,
                "bin_var": bin_var_idx,
                "sel": sel_idx,
                "slice_edges": sorted(all_slice_edges),
                "slices": slices,
                "n_bins": sum(s["n_bins"] for s in slices),
            })
        else:
            raise ValueError(f"Unexpected block header: {lines[i]}")

    return blocks


# ── Helpers ──────────────────────────────────────────────────────────────────

def _find_bin(value, edges):
    """Return the bin index for `value` in sorted edges. None if out of range."""
    if value < edges[0] or value > edges[-1]:
        return None
    for j in range(len(edges) - 1):
        if value < edges[j + 1] or j == len(edges) - 2:
            return j
    return None


def _find_slice(value, slices):
    """Return the slice dict whose range contains `value`."""
    max_high = max(s["high"] for s in slices)
    for s in slices:
        if s["low"] <= value < s["high"]:
            return s
        if value == s["high"] == max_high:
            return s
    return None


def _bin_offset_in_block(target_slice, slices):
    """Cumulative bin count before `target_slice`."""
    offset = 0
    for s in slices:
        if s["idx"] == target_slice["idx"]:
            return offset
        offset += s["n_bins"]
    return offset


def _block_global_bin(block, slice_mid, bin_mid):
    """Given midpoints in slice and bin variables, return the global bin index
    for a 2D block, or None if out of range."""
    orig_slice = _find_slice(slice_mid, block["slices"])
    if orig_slice is None:
        return None
    orig_bin = _find_bin(bin_mid, orig_slice["bin_edges"])
    if orig_bin is None:
        return None
    offset = _bin_offset_in_block(orig_slice, block["slices"])
    return block["_offset"] + offset + orig_bin


# ── Core ─────────────────────────────────────────────────────────────────────

def _get_block_bins_for_shared_cell(block, shared_mids, shared_vars):
    """For a given shared-variable cell (specified by midpoints), return the
    list of possible global bin indices for this block.

    Returns None if the shared cell is out of range for this block.
    Returns a list of global bin indices (length >= 1) if valid.
    """

    if block["dim"] == 1:
        v = block["var"]
        if v in shared_vars:
            orig_bin = _find_bin(shared_mids[v], block["edges"])
            if orig_bin is None:
                return None
            return [block["_offset"] + orig_bin]
        else:
            return list(range(block["_offset"],
                              block["_offset"] + block["n_bins"]))

    # 2D block
    vs = block["slice_var"]
    vb = block["bin_var"]
    stype = block["_shared_type"]

    if stype == "both_shared":
        gbin = _block_global_bin(block, shared_mids[vs], shared_mids[vb])
        if gbin is None:
            return None
        return [gbin]

    elif stype == "slice_shared":
        orig_slice = _find_slice(shared_mids[vs], block["slices"])
        if orig_slice is None:
            return None
        offset = block["_offset"] + _bin_offset_in_block(
            orig_slice, block["slices"])
        return list(range(offset, offset + orig_slice["n_bins"]))

    elif stype == "bin_shared":
        possible = set()
        for s in block["slices"]:
            orig_bin = _find_bin(shared_mids[vb], s["bin_edges"])
            if orig_bin is not None:
                offset = block["_offset"] + _bin_offset_in_block(
                    s, block["slices"])
                possible.add(offset + orig_bin)
        if not possible:
            return None
        return sorted(possible)

    elif stype == "both_private":
        return list(range(block["_offset"],
                          block["_offset"] + block["n_bins"]))

    return None


def _expand_representative(block_bin_options):
    """Given a list of lists (one per block) of possible global bin indices,
    return a representative set of indicator tuples that spans the same
    column space as the full Cartesian product.
    """
    combos = set()

    # Base combo: first option for every block
    base = tuple(opts[0] for opts in block_bin_options)
    combos.add(base)

    # For each block, vary its choice while others stay at base
    for k, opts in enumerate(block_bin_options):
        for j in range(1, len(opts)):
            row = list(base)
            row[k] = opts[j]
            combos.add(tuple(row))

    return combos


# Threshold: if the full Cartesian product exceeds this, use representative
_MAX_FULL_EXPAND = 500_000


def _build_A_for_selection(sel_blocks, n_bins_total, verbose=True,
                           sel_label=""):
    """Build the structural A matrix for a single selection group.

    The blocks already have their global _offset set (pointing into
    the full bin vector), so the resulting rows of A have the correct
    column positions even though we only enumerate over this selection's
    blocks.
    """

    prefix = f"[sel {sel_label}] " if sel_label else ""

    # ── Classify variables (within this selection only) ───────────
    var_usage = {}
    for b in sel_blocks:
        if b["dim"] == 1:
            var_usage.setdefault(b["var"], set()).add(b["idx"])
        elif b["dim"] == 2:
            var_usage.setdefault(b["slice_var"], set()).add(b["idx"])
            var_usage.setdefault(b["bin_var"], set()).add(b["idx"])

    shared_vars = {v for v, users in var_usage.items() if len(users) >= 2}
    private_vars = {v for v, users in var_usage.items() if len(users) == 1}

    if verbose:
        print(f"{prefix}Shared variables:  {sorted(shared_vars)}")
        print(f"{prefix}Private variables: {sorted(private_vars)}")

    # ── Refined edges for shared variables only ──────────────────────
    shared_var_edges = {}
    for block in sel_blocks:
        if block["dim"] == 1:
            v = block["var"]
            if v in shared_vars:
                shared_var_edges.setdefault(v, set()).update(block["edges"])
        elif block["dim"] == 2:
            vs = block["slice_var"]
            vb = block["bin_var"]
            if vs in shared_vars:
                shared_var_edges.setdefault(vs, set()).update(
                    block["slice_edges"])
            if vb in shared_vars:
                for s in block["slices"]:
                    shared_var_edges.setdefault(vb, set()).update(
                        s["bin_edges"])

    shared_refined = {v: sorted(edges)
                      for v, edges in shared_var_edges.items()}
    shared_n_refined = {v: len(e) - 1 for v, e in shared_refined.items()}
    shared_vars_sorted = sorted(shared_refined.keys())

    if verbose:
        for v in shared_vars_sorted:
            print(f"{prefix}  var {v}: {shared_n_refined[v]} refined bins")
        n_shared_cells = 1
        for v in shared_vars_sorted:
            n_shared_cells *= shared_n_refined[v]
        print(f"{prefix}Shared-variable cells to enumerate: "
              f"{n_shared_cells}")

    # ── Classify each block's variable usage ─────────────────────────
    for block in sel_blocks:
        if block["dim"] == 1:
            block["_shared_type"] = ("shared" if block["var"] in shared_vars
                                     else "private")
        elif block["dim"] == 2:
            s_shared = block["slice_var"] in shared_vars
            b_shared = block["bin_var"] in shared_vars
            if s_shared and b_shared:
                block["_shared_type"] = "both_shared"
            elif s_shared and not b_shared:
                block["_shared_type"] = "slice_shared"
            elif not s_shared and b_shared:
                block["_shared_type"] = "bin_shared"
            else:
                block["_shared_type"] = "both_private"

    # ── Enumerate shared cells and expand ────────────────────────────
    shared_ranges = [range(shared_n_refined[v]) for v in shared_vars_sorted]
    unique_combinations = set()
    n_shared_enum = 0
    n_total_combos = 0
    n_representative_used = 0

    t0 = time.time()

    for shared_cell in product(*shared_ranges):
        n_shared_enum += 1

        shared_mids = {}
        for i, v in enumerate(shared_vars_sorted):
            edges = shared_refined[v]
            shared_mids[v] = 0.5 * (edges[shared_cell[i]]
                                     + edges[shared_cell[i] + 1])

        block_bin_options = []
        valid = True

        for block in sel_blocks:
            opts = _get_block_bins_for_shared_cell(
                block, shared_mids, shared_vars)
            if opts is None:
                valid = False
                break
            block_bin_options.append(opts)

        if not valid:
            continue

        full_size = 1
        for opts in block_bin_options:
            full_size *= len(opts)
            if full_size > _MAX_FULL_EXPAND:
                break

        if full_size <= _MAX_FULL_EXPAND:
            for combo in product(*block_bin_options):
                n_total_combos += 1
                unique_combinations.add(combo)
        else:
            n_representative_used += 1
            rep = _expand_representative(block_bin_options)
            n_total_combos += len(rep)
            unique_combinations.update(rep)

    elapsed = time.time() - t0

    # ── Build matrix ─────────────────────────────────────────────────
    sorted_combos = sorted(unique_combinations)
    n_rows = len(sorted_combos)
    A = np.zeros((n_rows, n_bins_total), dtype=int)
    for i, combo in enumerate(sorted_combos):
        for global_bin in combo:
            A[i, global_bin] = 1

    if verbose:
        print(f"{prefix}Rows: {n_rows}, "
              f"combos expanded: {n_total_combos}, "
              f"representative used: {n_representative_used}, "
              f"time: {elapsed:.2f}s")

    return A


def build_structural_A(blocks, verbose=True):
    """Build the structural bin-combination matrix, handling multiple
    mutually exclusive selections.

    Blocks are grouped by their 'sel' index. Each selection group is
    processed independently (no event shares across selections), and
    the resulting matrices are stacked vertically.
    """

    # ── Global bin offsets (across ALL blocks, all selections) ────
    n_bins_total = 0
    for block in blocks:
        block["_offset"] = n_bins_total
        n_bins_total += block["n_bins"]

    if verbose:
        print(f"Total bins: {n_bins_total}")
        for b in blocks:
            print(f"  Block {b['idx']} (sel {b['sel']}): "
                  f"offset={b['_offset']}, n_bins={b['n_bins']}")

    # ── Group blocks by selection ────────────────────────────────
    sel_groups = {}
    for b in blocks:
        sel_groups.setdefault(b["sel"], []).append(b)

    sel_indices = sorted(sel_groups.keys())

    if verbose:
        print(f"\nSelections: {sel_indices}")
        for s in sel_indices:
            block_ids = [b["idx"] for b in sel_groups[s]]
            print(f"  sel {s}: blocks {block_ids}")
        print()

    # ── Build A per selection and stack ───────────────────────────
    A_parts = []
    for s in sel_indices:
        if verbose:
            print(f"── Selection {s} "
                  f"({len(sel_groups[s])} blocks) ──")
        A_sel = _build_A_for_selection(
            sel_groups[s], n_bins_total, verbose=verbose,
            sel_label=str(s))
        A_parts.append(A_sel)
        if verbose:
            print()

    A = np.vstack(A_parts)

    if verbose:
        print(f"Combined A: {A.shape[0]} rows × {A.shape[1]} columns")

    return A, n_bins_total


# ── Analysis ─────────────────────────────────────────────────────────────────

def analyze_structural_A(A, n_bins, blocks=None, verbose=True):
    """Compute and print rank and null-space dimension of structural A."""
    rank = np.linalg.matrix_rank(A)
    n_null = n_bins - rank

    if verbose:
        print(f"\n{'='*50}")
        print(f"  N_bins (total)          : {n_bins}")
        print(f"  N_rows (structural A)   : {A.shape[0]}")
        print(f"  rank(A_struct)          : {rank}")
        print(f"  N_null (structural)     : {n_null}")
        print(f"{'='*50}")

    if n_null > 0:
        if A.shape[0] > 10000:
            ATA = A.T.astype(float) @ A.astype(float)
            eigvals, eigvecs = np.linalg.eigh(ATA)
            tol = max(eigvals) * 1e-10
            null_mask = eigvals < tol
            null_vecs = eigvecs[:, null_mask]
        else:
            U, s, Vt = np.linalg.svd(A.astype(float), full_matrices=True)
            null_vecs = Vt[rank:].T

        if verbose and blocks is not None:
            print(f"\nNull-space basis vectors:")
            for col in range(null_vecs.shape[1]):
                vec = null_vecs[:, col]
                print(f"\n  Vector {col}:")
                for b in blocks:
                    o = b["_offset"]
                    n = b["n_bins"]
                    vals = vec[o:o + n]
                    nz = np.count_nonzero(np.abs(vals) > 1e-10)
                    if nz > 0:
                        sel_str = f", sel={b['sel']}" if "sel" in b else ""
                        print(f"    Block {b['idx']}{sel_str} "
                              f"(bins {o}-{o+n-1}): "
                              f"{np.round(vals, 4)}")
                    else:
                        sel_str = f", sel={b['sel']}" if "sel" in b else ""
                        print(f"    Block {b['idx']}{sel_str} "
                              f"(bins {o}-{o+n-1}): [all zero]")
        elif verbose:
            print(f"\nNull-space basis vectors (columns):")
            print(np.round(null_vecs, 4))

    return rank, n_null


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python structural.py <config_file>")
        sys.exit(1)

    blocks = parse_binning_config(sys.argv[1])
    print("Parsed blocks:")
    for b in blocks:
        sel_str = f", sel={b['sel']}"
        if b["dim"] == 1:
            print(f"  Block {b['idx']}: 1D, var={b['var']}{sel_str}, "
                  f"n_bins={b['n_bins']}")
        else:
            print(f"  Block {b['idx']}: 2D, slice_var={b['slice_var']}, "
                  f"bin_var={b['bin_var']}{sel_str}, "
                  f"n_bins={b['n_bins']}, n_slices={len(b['slices'])}")
    print()

    A, n_bins = build_structural_A(blocks)
    analyze_structural_A(A, n_bins, blocks=blocks)
    np.savetxt("structure_" + sys.argv[1], A, fmt="%d")