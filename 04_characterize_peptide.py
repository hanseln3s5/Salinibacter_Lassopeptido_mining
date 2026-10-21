#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_characterize_peptide.py - Physicochemical, surface, geometric and
topological characterisation of a modelled lasso peptide, from its PDB file.

WHY A PURPOSE-BUILT SCRIPT AND NOT PROTPARAM
--------------------------------------------
In a lasso peptide the alpha-amino group of residue 1 and the side-chain
carboxyl of the acceptor Asp/Glu are NOT free: they form the isopeptide bond
that closes the macrolactam. A general-purpose tool such as ProtParam counts
both as titratable groups and therefore returns an incorrect isoelectric point
and net charge. The error is not marginal - removing one positive and one
negative titratable group changes the shape of the whole titration curve, and
for a peptide with few ionisable residues it can move the pI by several units.

This script excludes both groups from the charge calculation and subtracts the
water molecule lost when the bond forms. Declared disulfide bonds are handled
the same way: bridged cysteine thiols no longer titrate, and each bridge costs
two hydrogens.

BOTH the naive and the corrected values are always reported. The difference is
itself a result worth documenting in the methods section, and reporting only
the corrected figure would leave a reader unable to compare with values
obtained from standard tools.

WHAT IT COMPUTES
----------------
  1. Topology    : isopeptide bond detection and length, ring size, and
                   verification of tail threading by ray-triangle intersection
                   (Moller-Trumbore). Repeating the check performed in step 03
                   is deliberate: the two scripts handle residue naming
                   differently, so agreement between them is a free cross-check.
  2. Physicochem : molecular weight, pI, net charge versus pH, aliphatic index,
                   GRAVY hydropathy, extinction coefficient, composition.
  3. Surface     : total, polar and apolar solvent-accessible surface area by
                   Shrake-Rupley; per-residue SASA and relative accessibility;
                   the list of exposed residues.
  4. Geometry    : radius of gyration, maximum dimension, centre of mass.

The exposed-residue list feeds directly into the docking step: those residues
are the defensible choice for HADDOCK passive residues on the ligand.

WHAT IT DELIBERATELY DOES NOT COMPUTE
-------------------------------------
The Guruprasad instability index. It requires a 400-entry dipeptide weight
table which is not reproduced here, because transcribing it from memory risks
silent numerical errors in a thesis. Take that value from ProtParam if needed.
Its meaning for a 15-residue constrained peptide is questionable in any case.

A NOTE ON THE ISOELECTRIC POINT
-------------------------------
When a peptide has no basic residues and its alpha-amino group is blocked by
the isopeptide bond, the net charge never changes sign between pH 0 and 14:
there is no accessible isoelectric point. The script detects this and reports
it explicitly rather than returning a meaningless pI of 0.00. That outcome is
a genuine physicochemical property of such peptides, not a failure of the
calculation.

USAGE
-----
  python3 04_characterize_peptide.py model_LINK.pdb
  python3 04_characterize_peptide.py model_LINK.pdb --isopeptide 1:9
  python3 04_characterize_peptide.py model_LINK.pdb --disulfide 6:15 --pymol

Use the LINK-bearing PDB written by 03_validate_lasso_fold.py, not the raw
predictor output.

OUTPUT (in the directory given by -o, default 'characterization/')
  report.txt              human-readable report
  summary.tsv             one row per quantity, for the annex
  per_residue.tsv         SASA, relative accessibility and exposure per residue
  charge_vs_pH.tsv        titration curve
  surface_apbs.pml        PyMOL script for the electrostatic surface

No dependencies: Python 3.8+ standard library only.
Licence: MIT
"""

import argparse
import math
import os
import sys
from collections import Counter, OrderedDict, defaultdict

# --------------------------------------------------------------- constants
# Average residue masses (Da). Peptide mass = sum of residues + one water.
AA_MASS = {'G': 57.0519, 'A': 71.0788, 'S': 87.0782, 'P': 97.1167, 'V': 99.1326,
           'T': 101.1051, 'C': 103.1388, 'L': 113.1594, 'I': 113.1594,
           'N': 114.1038, 'D': 115.0886, 'Q': 128.1307, 'K': 128.1741,
           'E': 129.1155, 'M': 131.1926, 'H': 137.1411, 'F': 147.1766,
           'R': 156.1875, 'Y': 163.1760, 'W': 186.2132}
WATER = 18.0153

# Kyte and Doolittle (1982) hydropathy scale.
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

# pKa values, EMBOSS/Expasy set. Stated explicitly because different sets give
# slightly different isoelectric points and the choice must be reproducible.
PKA_NTERM, PKA_CTERM = 8.6, 3.6
PKA_SIDECHAIN = {'D': 3.9, 'E': 4.1, 'C': 8.5, 'Y': 10.1,
                 'H': 6.5, 'K': 10.8, 'R': 12.5}
ACIDIC, BASIC = set('DECY'), set('HKR')

# Maximum accessible surface area per residue (Tien et al. 2013, theoretical
# values, A^2). Used as the denominator of relative solvent accessibility.
ASA_MAX = {'A': 129, 'R': 274, 'N': 195, 'D': 193, 'C': 167, 'E': 223, 'Q': 225,
           'G': 104, 'H': 224, 'I': 197, 'L': 201, 'K': 236, 'M': 224, 'F': 240,
           'P': 159, 'S': 155, 'T': 172, 'W': 285, 'Y': 263, 'V': 174}

THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

VDW_RADIUS = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'H': 1.20}
ATOM_MASS = {'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06, 'H': 1.008}

# Side-chain carboxyl carbon that can form the isopeptide bond. Asn and Gln are
# absent on purpose: their side chains are amides and cannot close the ring.
CARBOXYL = {'ASP': 'CG', 'GLU': 'CD'}


# ------------------------------------------------------------------ PDB input
class Atom:
    __slots__ = ('name', 'resname', 'chain', 'resseq', 'x', 'y', 'z',
                 'element', 'radius', 'mass')

    def __init__(self, line):
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21].strip() or 'A'
        self.resseq = int(line[22:26])
        self.x, self.y, self.z = (float(line[30:38]), float(line[38:46]),
                                  float(line[46:54]))
        e = line[76:78].strip().upper()
        if not e:
            e = ''.join(c for c in self.name if c.isalpha())[:1].upper()
        self.element = e
        self.radius = VDW_RADIUS.get(e, 1.70)
        self.mass = ATOM_MASS.get(e, 12.0)


def read_pdb(path, keep_hydrogens=False):
    """Parse the first model only, dropping waters, hydrogens and alternates.

    Only the first MODEL is read: a multi-model NMR-style file would otherwise
    superimpose several conformers and inflate every surface and distance.
    """
    atoms, links = [], []
    with open(path, errors='replace') as fh:
        for line in fh:
            if line.startswith('LINK'):
                links.append(line.rstrip('\n'))
            if line.startswith('ENDMDL'):
                break
            if not line.startswith(('ATOM', 'HETATM')) or len(line) < 54:
                continue
            altloc = line[16].strip()
            if altloc and altloc not in ('A', '1'):
                continue
            if line[17:20].strip() in ('HOH', 'WAT'):
                continue
            try:
                a = Atom(line)
            except ValueError:
                continue
            if not keep_hydrogens and a.element == 'H':
                continue
            atoms.append(a)
    return atoms, links


def group_residues(atoms):
    """OrderedDict {resseq: {'resname', 'one_letter', 'atoms'}} in file order."""
    d = OrderedDict()
    for a in atoms:
        if a.resseq not in d:
            d[a.resseq] = {'resname': a.resname,
                           'one_letter': THREE_TO_ONE.get(a.resname, 'X'),
                           'atoms': []}
        d[a.resseq]['atoms'].append(a)
    return d


def atom_named(residue, name):
    for a in residue['atoms']:
        if a.name == name:
            return a
    return None


# ------------------------------------------------------------------- geometry
def distance(a, b):
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def moller_trumbore(origin, direction, v0, v1, v2, eps=1e-9):
    """Segment-triangle intersection. Returns t in [0,1] if it crosses."""
    e1 = [v1[i] - v0[i] for i in range(3)]
    e2 = [v2[i] - v0[i] for i in range(3)]
    p = [direction[1] * e2[2] - direction[2] * e2[1],
         direction[2] * e2[0] - direction[0] * e2[2],
         direction[0] * e2[1] - direction[1] * e2[0]]
    det = sum(e1[i] * p[i] for i in range(3))
    if abs(det) < eps:
        return None
    inv = 1.0 / det
    s = [origin[i] - v0[i] for i in range(3)]
    u = inv * sum(s[i] * p[i] for i in range(3))
    if u < 0 or u > 1:
        return None
    q = [s[1] * e1[2] - s[2] * e1[1],
         s[2] * e1[0] - s[0] * e1[2],
         s[0] * e1[1] - s[1] * e1[0]]
    v = inv * sum(direction[i] * q[i] for i in range(3))
    if v < 0 or u + v > 1:
        return None
    t = inv * sum(e2[i] * q[i] for i in range(3))
    return t if 0.0 <= t <= 1.0 else None


def detect_isopeptide(residues, forced=None, cutoff=2.0):
    """Locate the isopeptide bond. Returns (res1, acceptor_res, distance).

    The cutoff here is tighter than in step 03 (2.0 vs 4.0 A) because this
    script runs on a model that has already passed validation, so a genuine
    bond is expected at ~1.33 A and a looser threshold would only invite false
    assignments when several acidic residues sit near the N-terminus.
    """
    keys = list(residues.keys())
    if not keys:
        return None
    first = keys[0]
    n_term = atom_named(residues[first], 'N')
    if n_term is None:
        return None
    if forced:
        r1, racc = forced
        if racc not in residues:
            return None
        acc = atom_named(residues[racc],
                         CARBOXYL.get(residues[racc]['resname'], 'CG'))
        return (r1, racc, distance(n_term, acc) if acc else float('nan'))
    best = None
    for ri in keys[1:]:
        name = CARBOXYL.get(residues[ri]['resname'])
        if not name:
            continue
        acc = atom_named(residues[ri], name)
        if acc is None:
            continue
        d = distance(n_term, acc)
        if d <= cutoff and (best is None or d < best[2]):
            best = (first, ri, d)
    return best


def verify_threading(residues, r1, racc):
    """Count how many times the tail crosses the capped macrocycle surface."""
    keys = list(residues.keys())
    ring_keys = [k for k in keys if r1 <= k <= racc]
    tail_keys = [k for k in keys if k > racc]
    ca = {k: atom_named(residues[k], 'CA') for k in keys}
    ring_keys = [k for k in ring_keys if ca[k]]
    tail_keys = [k for k in tail_keys if ca[k]]
    if len(ring_keys) < 3 or len(tail_keys) < 2:
        return {'crossings': 0, 'threaded': False, 'ring': len(ring_keys),
                'tail': len(tail_keys), 'detail': []}
    points = [(ca[k].x, ca[k].y, ca[k].z) for k in ring_keys]
    centre = [sum(p[i] for p in points) / len(points) for i in range(3)]
    triangles = [(centre, points[i], points[(i + 1) % len(points)])
                 for i in range(len(points))]
    crossings, detail = 0, []
    for a, b in zip(tail_keys, tail_keys[1:]):
        origin = (ca[a].x, ca[a].y, ca[a].z)
        direction = (ca[b].x - origin[0], ca[b].y - origin[1], ca[b].z - origin[2])
        for v0, v1, v2 in triangles:
            t = moller_trumbore(origin, direction, v0, v1, v2)
            if t is not None:
                crossings += 1
                detail.append((a, b, round(t, 3)))
                break      # one crossing per segment; see 03 for the reason
    return {'crossings': crossings, 'threaded': crossings % 2 == 1,
            'ring': len(ring_keys), 'tail': len(tail_keys), 'detail': detail}


# ----------------------------------------------------------------------- SASA
def sphere_points(n):
    """Quasi-uniform points on the unit sphere (Fibonacci spiral)."""
    pts, phi = [], math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = 1 - (2.0 * i + 1.0) / n
        r = math.sqrt(max(0.0, 1 - y * y))
        t = phi * i
        pts.append((math.cos(t) * r, y, math.sin(t) * r))
    return pts


def solvent_accessible_area(atoms, probe=1.4, n_points=256):
    """Shrake-Rupley numerical SASA. Returns (per_atom, per_residue).

    The probe radius of 1.4 A is the conventional water radius. Accuracy is set
    by n_points; 256 is a good compromise for a small peptide, where the whole
    calculation takes a second or two.
    """
    pts = sphere_points(n_points)
    per_atom, per_residue = {}, defaultdict(float)
    for i, a in enumerate(atoms):
        ra = a.radius + probe
        neighbours = []
        for b in atoms:
            if b is a:
                continue
            rb = b.radius + probe
            limit = ra + rb
            d2 = (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
            if d2 < limit * limit:
                neighbours.append((b, rb * rb))
        exposed = 0
        for px, py, pz in pts:
            sx, sy, sz = a.x + ra * px, a.y + ra * py, a.z + ra * pz
            if not any((sx - b.x) ** 2 + (sy - b.y) ** 2 + (sz - b.z) ** 2 < rb2
                       for b, rb2 in neighbours):
                exposed += 1
        area = 4.0 * math.pi * ra * ra * exposed / n_points
        per_atom[i] = area
        per_residue[a.resseq] += area
    return per_atom, dict(per_residue)


# ------------------------------------------------------------- physicochemistry
def net_charge(seq, pH, nterm_free=True, cterm_free=True, excluded=None):
    """Henderson-Hasselbalch net charge.

    excluded: list of (1-based index, residue letter) whose side chain does NOT
    titrate, i.e. the macrolactam acceptor and any bridged cysteines.
    """
    excl = set(excluded or [])
    q = 0.0
    if nterm_free:
        q += 1.0 / (1.0 + 10 ** (pH - PKA_NTERM))
    if cterm_free:
        q -= 1.0 / (1.0 + 10 ** (PKA_CTERM - pH))
    for i, aa in enumerate(seq, 1):
        if (i, aa) in excl or aa not in PKA_SIDECHAIN:
            continue
        pk = PKA_SIDECHAIN[aa]
        if aa in BASIC:
            q += 1.0 / (1.0 + 10 ** (pH - pk))
        else:
            q -= 1.0 / (1.0 + 10 ** (pk - pH))
    return q


def isoelectric_point(seq, **kw):
    """pI by bisection, or None when the charge never changes sign.

    Returning None rather than a number matters: a peptide with no basic
    residues and a blocked alpha-amino group has no accessible pI, and
    reporting 0.00 would misrepresent that as a measured value.
    """
    q0, q14 = net_charge(seq, 0.0, **kw), net_charge(seq, 14.0, **kw)
    if q0 * q14 > 0:
        return None
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if net_charge(seq, mid, **kw) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def format_pi(value, seq, **kw):
    if value is not None:
        return '%12.2f' % value
    return '      %s' % ('   < 0' if net_charge(seq, 0.0, **kw) < 0 else '  > 14')


def molecular_weight(seq, isopeptide_bonds=0, disulfides=0):
    """Each isopeptide bond releases one water; each disulfide, two hydrogens."""
    m = sum(AA_MASS.get(a, 0.0) for a in seq) + WATER
    return m - isopeptide_bonds * WATER - disulfides * 2.016


def gravy(seq):
    values = [KD[a] for a in seq if a in KD]
    return sum(values) / len(values) if values else float('nan')


def aliphatic_index(seq):
    """Ikai (1980). Relative volume occupied by aliphatic side chains."""
    n = len(seq)
    if not n:
        return float('nan')
    pct = lambda a: 100.0 * seq.count(a) / n
    return pct('A') + 2.9 * pct('V') + 3.9 * (pct('I') + pct('L'))


def extinction_coefficient(seq, cystines=0):
    """Molar extinction at 280 nm (Gill and von Hippel, 1989)."""
    return 1490 * seq.count('Y') + 5500 * seq.count('W') + 125 * cystines


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description='Physicochemical, surface, geometric and topological '
                    'characterisation of a modelled lasso peptide.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n'
               '  python3 04_characterize_peptide.py model_LINK.pdb '
               '--isopeptide 1:9 --disulfide 6:15\n')
    ap.add_argument('pdb', help='PDB model, preferably the LINK-bearing file '
                                'written by 03_validate_lasso_fold.py')
    ap.add_argument('-o', '--outdir', default='characterization')
    ap.add_argument('--isopeptide', metavar='R1:RACC',
                    help='force the isopeptide bond between residue R1 and '
                         'acceptor RACC when auto-detection fails')
    ap.add_argument('--disulfide', metavar='RA:RB', action='append', default=[],
                    help='disulfide bridge between two Cys; repeatable. Their '
                         'thiols stop being counted as titratable')
    ap.add_argument('--exposed-rsa', type=float, default=25.0,
                    help='relative accessibility above which a residue counts '
                         'as exposed (25%%)')
    ap.add_argument('--points', type=int, default=256,
                    help='sphere points per atom in Shrake-Rupley (256)')
    ap.add_argument('--pymol', action='store_true',
                    help='also write the PyMOL script for the APBS surface')
    args = ap.parse_args()

    if not os.path.isfile(args.pdb):
        sys.exit('Not found: %s' % args.pdb)
    os.makedirs(args.outdir, exist_ok=True)

    atoms, links = read_pdb(args.pdb)
    if not atoms:
        sys.exit('No atoms were read from %s.' % args.pdb)
    residues = group_residues(atoms)
    keys = list(residues.keys())
    seq = ''.join(residues[k]['one_letter'] for k in keys)
    n = len(seq)

    out = []
    W = out.append
    W('=' * 72)
    W('LASSO PEPTIDE MODEL CHARACTERISATION')
    W('=' * 72)
    W('File      : %s' % os.path.abspath(args.pdb))
    W('Residues  : %d' % n)
    W('Atoms     : %d (hydrogens excluded)' % len(atoms))
    W('Sequence  : %s' % seq)
    W('')

    # ------------------------------------------------------------- 1. topology
    W('-' * 72)
    W('1. TOPOLOGY')
    W('-' * 72)
    forced = None
    if args.isopeptide:
        a, b = args.isopeptide.split(':')
        forced = (int(a), int(b))
    iso = detect_isopeptide(residues, forced)
    threading = None
    if iso is None:
        W('  [!] No isopeptide bond detected (N-terminus within 2.0 A of an')
        W('      Asp/Glu side-chain carboxyl).')
        W('      If the model has one, declare it with --isopeptide R1:RACC.')
        W('      WARNING: without a closed macrolactam this is not a lasso peptide.')
    else:
        r1, racc, d = iso
        W('  Isopeptide bond  : residue %d (N-terminus) -- residue %d (%s)'
          % (r1, racc, residues[racc]['resname']))
        W('  N--C distance    : %.3f A  [%s]'
          % (d, 'consistent with a covalent bond' if d < 4.0 else 'TOO LONG'))
        W('  Ring size        : %d residues' % (racc - r1 + 1))
        threading = verify_threading(residues, r1, racc)
        W('  Tail residues    : %d' % threading['tail'])
        W('  Tail/ring crossings: %d  (Moller-Trumbore over the triangle fan'
          % threading['crossings'])
        W('                       of the ring CA atoms)')
        W('')
        if threading['threaded']:
            W('  >>> LASSO TOPOLOGY CONFIRMED (odd number of crossings)')
        elif threading['crossings'] == 0:
            W('  >>> NOT THREADED: the tail does not pass through the ring.')
            W('      This is a branched macrocycle, not a lasso peptide.')
        else:
            W('  >>> EVEN NUMBER OF CROSSINGS (%d): the tail enters and exits.'
              % threading['crossings'])
            W('      This is not a lasso topology.')
        if threading['detail']:
            W('      Crossing segments: ' +
              ', '.join('%d-%d' % (a, b) for a, b, _ in threading['detail']))
        W('')
        W('  Confirm independently with LassoProt. Check that the PDB carries')
        W('  the LINK record for the isopeptide bond; without it the server')
        W('  reads the chain as linear.')
        W('  LINK records in this file: %d' % len(links))
    W('')

    # ------------------------------------------------------ 2. physicochemistry
    W('-' * 72)
    W('2. PHYSICOCHEMICAL PROPERTIES')
    W('-' * 72)

    excluded = []
    iso_note = ''
    if iso:
        r1, racc, _ = iso
        acc_index = keys.index(racc) + 1
        excluded.append((acc_index, residues[racc]['one_letter']))
        iso_note = ('the alpha-amino group of residue %d and the side-chain '
                    'carboxyl of residue %d (%s) are engaged in the isopeptide '
                    'bond' % (r1, racc, residues[racc]['resname']))
    bridges = []
    for spec in args.disulfide:
        a, b = (int(x) for x in spec.split(':'))
        for r in (a, b):
            if r in keys:
                excluded.append((keys.index(r) + 1, residues[r]['one_letter']))
        bridges.append((a, b))

    mw_naive = molecular_weight(seq, 0, 0)
    mw_corrected = molecular_weight(seq, 1 if iso else 0, len(bridges))
    corrected_kw = dict(nterm_free=not bool(iso), excluded=excluded)
    pi_naive = isoelectric_point(seq)
    pi_corrected = isoelectric_point(seq, **corrected_kw)
    q_naive = net_charge(seq, 7.0)
    q_corrected = net_charge(seq, 7.0, **corrected_kw)

    W('                                    uncorrected      corrected')
    W('  Molecular weight (Da)          %12.2f  %12.2f' % (mw_naive, mw_corrected))
    W('  Isoelectric point (pI)         %s  %s'
      % (format_pi(pi_naive, seq), format_pi(pi_corrected, seq, **corrected_kw)))
    W('  Net charge at pH 7.0           %12.2f  %12.2f' % (q_naive, q_corrected))
    W('')
    if pi_corrected is None:
        W('  [!] The net charge does not change sign between pH 0 and 14: the')
        W('      molecule has no accessible isoelectric point and stays')
        W('      negatively charged across the whole range. This follows from')
        W('      having no basic residues and an alpha-amino group blocked by')
        W('      the isopeptide bond. It is a result, not a calculation error.')
        W('')
    if iso:
        W('  Correction applied: %s,' % iso_note)
        W('  so neither group titrates. The weight also subtracts the water')
        W('  molecule released when the bond forms.')
    if bridges:
        W('  Declared disulfide bridges: %s'
          % ', '.join('Cys%d-Cys%d' % p for p in bridges))
        W('  Their thiols are likewise excluded from the titration.')
    W('')
    W('  Aliphatic index                %12.2f' % aliphatic_index(seq))
    W('  GRAVY (Kyte-Doolittle)         %12.3f  [%s]'
      % (gravy(seq), 'hydrophilic' if gravy(seq) < 0 else 'hydrophobic'))
    W('  Extinction coeff. 280 nm       %12d  M-1 cm-1'
      % extinction_coefficient(seq, len(bridges)))
    W('')
    composition = Counter(seq)
    n_acidic = sum(composition[a] for a in 'DE')
    n_basic = sum(composition[a] for a in 'KR') + composition['H']
    W('  Acidic residues (D+E)          %12d  (%.1f %%)'
      % (n_acidic, 100.0 * n_acidic / n))
    W('  Basic residues (K+R+H)         %12d  (%.1f %%)'
      % (n_basic, 100.0 * n_basic / n))
    W('  Acid-base balance              %12d' % (n_acidic - n_basic))
    W('  Composition: ' + '  '.join('%s%d' % (a, c)
                                    for a, c in sorted(composition.items())))
    W('')
    W('  The Guruprasad instability index is not computed here; it needs a')
    W('  400-entry dipeptide table. Take it from ProtParam if required.')
    W('')

    # -------------------------------------------------------------- 3. surface
    W('-' * 72)
    W('3. SOLVENT-ACCESSIBLE SURFACE')
    W('-' * 72)
    per_atom, per_residue = solvent_accessible_area(atoms, n_points=args.points)
    total = sum(per_atom.values())
    polar = sum(a for i, a in per_atom.items() if atoms[i].element in ('N', 'O'))
    apolar = total - polar
    W('  Total SASA                     %12.1f A^2' % total)
    W('  Polar SASA   (N, O)            %12.1f A^2  (%.1f %%)'
      % (polar, 100.0 * polar / total if total else 0))
    W('  Apolar SASA  (C, S)            %12.1f A^2  (%.1f %%)'
      % (apolar, 100.0 * apolar / total if total else 0))
    W('')
    exposed, buried, rows = [], [], []
    for i, k in enumerate(keys, 1):
        one = residues[k]['one_letter']
        area = per_residue.get(k, 0.0)
        maximum = ASA_MAX.get(one)
        rsa = 100.0 * area / maximum if maximum else float('nan')
        zone = ('ring' if iso and iso[0] <= k <= iso[1]
                else ('tail' if iso else '-'))
        (exposed if rsa >= args.exposed_rsa else buried).append('%s%d' % (one, k))
        rows.append((i, k, one, residues[k]['resname'], zone, area, rsa))
    W('  Exposed residues (RSA >= %.0f %%): %s'
      % (args.exposed_rsa, ', '.join(exposed) or 'none'))
    W('  Buried residues                : %s' % (', '.join(buried) or 'none'))
    W('')
    W('  The exposed set is the defensible choice for HADDOCK passive residues')
    W('  on the ligand; full detail in per_residue.tsv')
    W('')

    # ------------------------------------------------------------- 4. geometry
    W('-' * 72)
    W('4. GEOMETRY')
    W('-' * 72)
    mass = sum(a.mass for a in atoms)
    cx = sum(a.mass * a.x for a in atoms) / mass
    cy = sum(a.mass * a.y for a in atoms) / mass
    cz = sum(a.mass * a.z for a in atoms) / mass
    rg = math.sqrt(sum(a.mass * ((a.x - cx) ** 2 + (a.y - cy) ** 2
                                 + (a.z - cz) ** 2) for a in atoms) / mass)
    dmax, pair = 0.0, None
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            d = distance(atoms[i], atoms[j])
            if d > dmax:
                dmax, pair = d, (atoms[i], atoms[j])
    W('  Radius of gyration             %12.2f A' % rg)
    W('  Maximum dimension              %12.2f A  (%s%d %s -- %s%d %s)'
      % (dmax, pair[0].resname, pair[0].resseq, pair[0].name,
         pair[1].resname, pair[1].resseq, pair[1].name))
    W('  Centre of mass                 (%.2f, %.2f, %.2f)' % (cx, cy, cz))
    W('')
    W('=' * 72)

    report = '\n'.join(out)
    print(report)
    with open(os.path.join(args.outdir, 'report.txt'), 'w') as fh:
        fh.write(report + '\n')

    # ------------------------------------------------------------------ tables
    with open(os.path.join(args.outdir, 'per_residue.tsv'), 'w') as fh:
        fh.write('position\tpdb_resseq\taa\tresname\tzone\tSASA_A2\tRSA_pct\texposed\n')
        for i, k, one, resname, zone, area, rsa in rows:
            fh.write('%d\t%d\t%s\t%s\t%s\t%.2f\t%.1f\t%s\n'
                     % (i, k, one, resname, zone, area, rsa,
                        'yes' if rsa >= args.exposed_rsa else 'no'))

    with open(os.path.join(args.outdir, 'summary.tsv'), 'w') as fh:
        fh.write('quantity\tvalue\tunit\n')
        R = lambda k, v, u='': fh.write('%s\t%s\t%s\n' % (k, v, u))
        R('sequence', seq)
        R('n_residues', n)
        R('mw_uncorrected', '%.2f' % mw_naive, 'Da')
        R('mw_corrected', '%.2f' % mw_corrected, 'Da')
        R('pI_uncorrected', '%.2f' % pi_naive if pi_naive is not None
          else 'no pI in 0-14')
        R('pI_corrected', '%.2f' % pi_corrected if pi_corrected is not None
          else 'no pI in 0-14')
        R('charge_pH7_uncorrected', '%.2f' % q_naive)
        R('charge_pH7_corrected', '%.2f' % q_corrected)
        R('aliphatic_index', '%.2f' % aliphatic_index(seq))
        R('GRAVY', '%.3f' % gravy(seq))
        R('extinction_280', extinction_coefficient(seq, len(bridges)), 'M-1 cm-1')
        R('acidic_residues', n_acidic)
        R('basic_residues', n_basic)
        R('SASA_total', '%.1f' % total, 'A^2')
        R('SASA_polar', '%.1f' % polar, 'A^2')
        R('SASA_apolar', '%.1f' % apolar, 'A^2')
        R('radius_of_gyration', '%.2f' % rg, 'A')
        R('maximum_dimension', '%.2f' % dmax, 'A')
        if iso:
            R('isopeptide_residues', '%d-%d' % (iso[0], iso[1]))
            R('isopeptide_distance', '%.3f' % iso[2], 'A')
            R('ring_size', iso[1] - iso[0] + 1, 'residues')
        if threading:
            R('tail_ring_crossings', threading['crossings'])
            R('lasso_topology', 'yes' if threading['threaded'] else 'no')

    with open(os.path.join(args.outdir, 'charge_vs_pH.tsv'), 'w') as fh:
        fh.write('pH\tcharge_uncorrected\tcharge_corrected\n')
        for i in range(0, 141):
            ph = i / 10.0
            fh.write('%.1f\t%.3f\t%.3f\n'
                     % (ph, net_charge(seq, ph),
                        net_charge(seq, ph, **corrected_kw)))

    if args.pymol:
        pml = os.path.join(args.outdir, 'surface_apbs.pml')
        with open(pml, 'w') as fh:
            fh.write(
                '# Electrostatic surface potential.\n'
                '# First: pdb2pqr --ff=AMBER --with-ph=7.0 %s model.pqr\n'
                '#        apbs apbs.in           (produces model.dx)\n'
                '# Then:  pymol %s\n\n'
                'load model.pqr\n'
                'load model.dx, potential\n'
                'set surface_quality, 1\n'
                'set transparency, 0.0\n'
                'hide everything\n'
                'show surface\n'
                'ramp_new ramp, potential, [-5, 0, 5], [red, white, blue]\n'
                'set surface_color, ramp\n'
                'bg_color white\n'
                'set ray_opaque_background, 0\n'
                'orient\n'
                'ray 2000, 1600\n'
                'png electrostatic_surface.png, dpi=300\n'
                % (os.path.basename(args.pdb), os.path.basename(pml)))

    print('\nWritten to %s/: report.txt, summary.tsv, per_residue.tsv, '
          'charge_vs_pH.tsv%s'
          % (args.outdir, ', surface_apbs.pml' if args.pymol else ''))


if __name__ == '__main__':
    main()
