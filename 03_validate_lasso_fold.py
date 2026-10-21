#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_validate_lasso_fold.py - Local geometric and topological validation of a
modelled lasso peptide.

WHY THIS STEP EXISTS
--------------------
A lasso peptide is defined by its threaded topology, not by its sequence. A
structure predictor can return a model with excellent local geometry that is
nevertheless NOT a lasso peptide, because the C-terminal tail does not actually
pass through the macrolactam ring. Reporting such a model as a lasso peptide
would be a category error, so the topology must be verified explicitly before
any downstream characterisation or docking is performed.

Two independent criteria are checked here:

  1. GEOMETRIC - is the macrolactam closed? The isopeptide bond links the
     alpha-amino group of residue 1 to the side-chain carboxyl of an Asp or Glu
     at position 7, 8 or 9. We measure that N...C distance. A real covalent
     bond is ~1.33 A; we accept < 4.0 A as "formed" to tolerate the coordinate
     error of predicted models, and report the actual value so the reader can
     judge.

  2. TOPOLOGICAL - does the tail thread the ring? The ring CA atoms are
     triangulated as a fan from their centroid, and every consecutive CA-CA
     segment of the tail is tested for intersection with that surface using the
     Moller-Trumbore ray-triangle algorithm. An ODD number of crossings means
     the tail enters and does not come back out, i.e. it is threaded. An even
     number means it enters and exits again, which is not a lasso fold.

     Choosing the parity of the crossing count, rather than merely "does it
     cross at all", is what distinguishes a true lasso (class L1 in LassoProt
     terminology) from a tail that merely brushes through the ring plane twice.

The script also writes the PDB LINK record for the isopeptide bond. This is not
cosmetic: LassoProt, CHARMM-GUI and most molecular-dynamics setup tools do not
infer non-canonical covalent bonds from geometry alone. Without the LINK record
the macrolactam is read as absent, LassoProt classifies the chain as linear, and
an MD topology built from the file will let the ring fall open during the
simulation. Emitting the record here is therefore a prerequisite for every
subsequent step of the pipeline.

The validation implemented here is deliberately independent of LassoProt so
that the result can be audited without relying on an external web service.
Agreement between the two is reported in the thesis as a cross-check.

USAGE
-----
    python3 03_validate_lasso_fold.py model.pdb
    python3 03_validate_lasso_fold.py model.pdb --acceptor 9
    python3 03_validate_lasso_fold.py model.pdb --list
    python3 03_validate_lasso_fold.py model.pdb --out model_LINK.pdb

EXIT STATUS
-----------
    0  the model is compatible with a lasso fold
    1  the model fails one of the two criteria

This makes the script usable as a gate inside a larger shell pipeline.

The parser understands the AMBER residue naming found in minimised PDB files
(GLH, ASH, GLX, ASX, HIE/HID/HIP, and the N-/C-terminal prefixes such as NGLU).
When a residue name is not recognised, the side-chain carboxyl is identified
from the atom names actually present, which keeps the script working on files
written by tools with idiosyncratic naming.

Requirements: numpy only.
Licence: MIT
"""

import argparse
import sys

import numpy as np

# Maximum N...C distance accepted as evidence of a formed isopeptide bond.
# A genuine amide C-N bond is ~1.33 A. The threshold is deliberately generous
# because predicted models carry coordinate error; the measured value is always
# printed so the reader is not forced to trust the threshold.
BOND_CUTOFF = 4.0

# Side-chain carboxyl carbon, by residue name. Includes the AMBER variants:
# GLH/ASH (protonated), GLX/ASX (ambiguous). Asn and Gln are listed explicitly
# with None because their side chains are amides, not carboxyls, and so cannot
# form the macrolactam - listing them prevents a silent false positive.
SIDE_CHAIN_CARBOXYL = {'ASP': 'CG', 'ASH': 'CG', 'ASX': 'CG', 'ASN': None,
                       'GLU': 'CD', 'GLH': 'CD', 'GLX': 'CD', 'GLN': None}


def side_chain_carboxyl(atoms, resname):
    """Return the side-chain carboxyl carbon of a residue, or None.

    The residue name is tried first, then the four-letter terminal variants
    (NGLU -> GLU). If the name is still unknown, the carboxyl is inferred from
    the atoms present, which is what rescues non-standard naming.
    """
    base = resname.lstrip('NC') if len(resname) == 4 else resname   # NGLU -> GLU
    for name in (resname, base):
        if name in SIDE_CHAIN_CARBOXYL:
            c = SIDE_CHAIN_CARBOXYL[name]
            return c if (c and c in atoms) else None
    if 'CG' in atoms and {'OD1', 'OD2'} & set(atoms):               # aspartate-like
        return 'CG'
    if 'CD' in atoms and {'OE1', 'OE2'} & set(atoms):               # glutamate-like
        return 'CD'
    return None


# --------------------------------------------------------------- PDB parsing
def read_pdb(path, chain=None):
    """Return ({(resseq, resname): {atom_name: xyz}}, ordered_keys).

    Only the first alternate conformer is kept: a model carrying alternates
    would otherwise contribute two coordinates for the same atom and corrupt
    both the distance measurement and the ring triangulation.
    """
    residues = {}
    order = []
    names = {}           # residue number -> name, fixed by its first atom
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            altloc = line[16]
            if altloc not in (' ', 'A'):      # skip alternate conformers
                continue
            ch = line[21]
            if chain and ch != chain:
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            key = (resseq, names.setdefault(resseq, resname))
            if key not in residues:
                residues[key] = {}
                order.append(key)
            residues[key][atom_name] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not residues:
        sys.exit('No atoms were read from %s%s'
                 % (path, ' in chain %s' % chain if chain else ''))
    return residues, order


# ------------------------------------------------------------- knot geometry
def intersects(p0, p1, v0, v1, v2):
    """Moller-Trumbore: does segment p0->p1 pass through triangle v0,v1,v2?

    The test is restricted to 0 < t < 1 so that only the segment itself counts,
    not its infinite extension. Without that restriction a tail pointing away
    from the ring would still register a crossing.
    """
    d = p1 - p0
    e1, e2 = v1 - v0, v2 - v0
    h = np.cross(d, e2)
    a = np.dot(e1, h)
    if abs(a) < 1e-9:
        return False                      # segment parallel to the triangle
    f = 1.0 / a
    s = p0 - v0
    u = f * np.dot(s, h)
    if u < 0.0 or u > 1.0:
        return False
    q = np.cross(s, e1)
    v = f * np.dot(d, q)
    if v < 0.0 or u + v > 1.0:
        return False
    t = f * np.dot(e2, q)
    return 0.0 < t < 1.0                  # inside the segment, not beyond it


def ring_crossings(ring, tail):
    """Number of times the tail passes through the capped macrocycle surface.

    The ring has no intrinsic surface, so one is constructed: a triangle fan
    from the centroid of the ring CA atoms to each consecutive pair. This is
    the standard way to cap a closed polygon for intersection testing, and it
    is exact for a convex ring and a good approximation otherwise.

    At most one crossing is counted per tail segment. A segment that passes
    exactly through a shared triangle edge would otherwise be counted twice and
    flip the parity, turning a genuine lasso into an apparent non-lasso.
    """
    centre = ring.mean(axis=0)
    n = len(ring)
    total = 0
    for i in range(len(tail) - 1):
        p0, p1 = tail[i], tail[i + 1]
        for j in range(n):                # triangle fan from the centroid
            if intersects(p0, p1, centre, ring[j], ring[(j + 1) % n]):
                total += 1
                break                     # one crossing per segment, see above
    return total


# ---------------------------------------------------------------- LINK record
def link_record(at1, res1, ch1, seq1, at2, res2, ch2, seq2, dist):
    """Build a PDB LINK record in fixed-column format.

    The PDB format is column-based, not whitespace-delimited. Tools that read
    LINK records parse them by column, so the fields are placed by absolute
    position rather than assembled with string formatting that merely looks
    right. '1555' is the identity symmetry operator for both partners, i.e. the
    bond is within the asymmetric unit.
    """
    buf = list(' ' * 80)

    def put(start, text):                 # start is a 1-based column number
        buf[start - 1:start - 1 + len(text)] = list(text)

    put(1, 'LINK')
    put(13, '%s' % at1.center(4) if len(at1) < 4 else at1)
    put(18, '%3s' % res1)
    put(22, ch1)
    put(23, '%4d' % seq1)
    put(43, '%s' % at2.center(4) if len(at2) < 4 else at2)
    put(48, '%3s' % res2)
    put(52, ch2)
    put(53, '%4d' % seq2)
    put(60, '%6s' % '1555')
    put(67, '%6s' % '1555')
    put(74, '%5.2f' % dist)
    return ''.join(buf).rstrip()


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(
        description='Local validation of lasso-peptide topology: isopeptide '
                    'bond geometry, tail threading, and PDB LINK record.')
    ap.add_argument('pdb', help='modelled structure, e.g. LassoPred output')
    ap.add_argument('--acceptor', type=int, default=None,
                    help='residue whose side-chain carboxyl closes the ring; '
                         'auto-detected if omitted')
    ap.add_argument('--list', action='store_true', dest='list_residues',
                    help='print residues with their atoms and exit; use this '
                         'when auto-detection fails to see the actual naming')
    ap.add_argument('--chain', default=None,
                    help='chain to analyse (default: all chains read)')
    ap.add_argument('--out', default=None,
                    help='output PDB carrying the LINK record '
                         '(default: <input>_LINK.pdb)')
    args = ap.parse_args()

    residues, order = read_pdb(args.pdb, args.chain)
    first = order[0]
    chain_id = args.chain or 'A'

    if args.list_residues:
        print('%d residues in %s\n' % (len(order), args.pdb))
        for key in order:
            atoms = ' '.join(sorted(residues[key]))
            flag = ('  <-- side-chain carboxyl'
                    if side_chain_carboxyl(residues[key], key[1]) else '')
            print('  %4s%-4d %s%s' % (key[1], key[0], atoms, flag))
        return 0

    # --- 1. geometric criterion -------------------------------------------
    if 'N' not in residues[first]:
        sys.exit('Residue %s%d has no N atom; cannot locate the N-terminus.'
                 % (first[1], first[0]))
    n_term = residues[first]['N']

    # Candidates: every residue carrying a side-chain carboxyl. Ranking them by
    # distance and reporting the runners-up matters, because a model with two
    # acidic residues close to the N-terminus is ambiguous and the reader
    # should see that ambiguity rather than a single confident answer.
    candidates = []
    for key in order:
        c = side_chain_carboxyl(residues[key], key[1])
        if c:
            candidates.append(
                (key, c, float(np.linalg.norm(n_term - residues[key][c]))))
    if not candidates:
        sys.exit('No residue with a side-chain carboxyl was found.\n'
                 'Run with --list to inspect residue and atom naming.')
    candidates.sort(key=lambda x: x[2])

    if args.acceptor is None:
        acceptor, acc_atom, _ = candidates[0]
        print('[acceptor auto-detected: %s%d]' % (acceptor[1], acceptor[0]))
    else:
        chosen = [c for c in candidates if c[0][0] == args.acceptor]
        if not chosen:
            present = next((k for k in order if k[0] == args.acceptor), None)
            msg = ('Residue %d is %s and has no side-chain carboxyl.'
                   % (args.acceptor, present[1]) if present
                   else 'Residue %d is not present in the file.' % args.acceptor)
            near = ', '.join('%s%d (%.1f A)' % (k[1], k[0], d)
                             for k, _, d in candidates[:4])
            sys.exit('%s\nResidues with a side-chain carboxyl: %s\n'
                     'Use --list to inspect the file.' % (msg, near))
        acceptor, acc_atom, _ = chosen[0]
    acc_resname = acceptor[1]
    acc_coord = residues[acceptor][acc_atom]

    if len(candidates) > 1:
        print('  other residues with a side-chain carboxyl: '
              + ', '.join('%s%d %.1f A' % (k[1], k[0], d)
                          for k, _, d in candidates[1:5]))

    distance = float(np.linalg.norm(n_term - acc_coord))
    geometry_ok = distance < BOND_CUTOFF

    print('=' * 64)
    print('  File        : %s' % args.pdb)
    print('  Macrocycle  : %s%d (N) -- %s%d (%s)'
          % (first[1], first[0], acc_resname, acceptor[0], acc_atom))
    print('=' * 64)
    print('\n1. GEOMETRIC CRITERION')
    print('   N(term)...%s distance = %.2f A   (threshold < %.1f A)'
          % (acc_atom, distance, BOND_CUTOFF))
    print('   -> %s' % ('PASS' if geometry_ok
                        else 'FAIL: the isopeptide bond is not formed'))

    # --- 2. topological criterion ------------------------------------------
    # The ring is everything from the N-terminus up to and including the
    # acceptor; the tail is everything after it. This follows directly from the
    # definition of the macrolactam and requires no further assumption.
    acc_index = order.index(acceptor)
    ring_keys = order[:acc_index + 1]
    tail_keys = order[acc_index + 1:]
    ring = np.array([residues[k]['CA'] for k in ring_keys if 'CA' in residues[k]])
    tail = np.array([residues[k]['CA'] for k in tail_keys if 'CA' in residues[k]])

    print('\n2. TAIL THREADING')
    print('   ring: %d residues   |   tail: %d residues' % (len(ring), len(tail)))
    if len(tail) < 2:
        print('   -> tail too short to define a segment; cannot be evaluated')
        crossings, topology_ok = 0, False
    else:
        crossings = ring_crossings(ring, tail)
        topology_ok = crossings % 2 == 1
        print('   tail crossings through the macrocycle: %d' % crossings)
        print('   -> %s' % ('THREADED (consistent with class L1)' if topology_ok
                            else 'NOT THREADED (consistent with class L0)'))

    # --- 3. LINK record for downstream tools -------------------------------
    link = link_record('N', first[1], chain_id, first[0],
                       acc_atom, acc_resname, chain_id, acceptor[0], distance)
    print('\n3. LINK RECORD FOR LassoProt / CHARMM-GUI / MD SETUP')
    print('   %s' % link)

    destination = args.out or args.pdb.replace('.pdb', '_LINK.pdb')
    with open(args.pdb, encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()
    # The LINK record must precede the coordinate section, per the PDB format.
    cut = next((i for i, x in enumerate(lines)
                if x.startswith(('ATOM  ', 'HETATM'))), 0)
    lines.insert(cut, link + '\n')
    with open(destination, 'w', encoding='utf-8') as fh:
        fh.writelines(lines)
    print('   written to: %s' % destination)
    print('   Use THIS file, not the original, for every downstream step.')

    # --- verdict -----------------------------------------------------------
    print('\n' + '=' * 64)
    if geometry_ok and topology_ok:
        print('  VERDICT: consistent with a lasso (threaded) topology.')
        print('  Confirm independently by submitting the LINK-bearing PDB to')
        print('  LassoProt; agreement between the two is reported as a check.')
    elif geometry_ok and not topology_ok:
        print('  VERDICT: the macrolactam is closed but the tail does NOT pass')
        print('  through it. This is a branched macrocycle, not a lasso peptide.')
        print('  Re-examine the precursor annotation and the predicted model.')
    else:
        print('  VERDICT: the isopeptide bond is not formed. The model does not')
        print('  represent a lasso peptide and must not be used downstream.')
    print('=' * 64)
    return 0 if (geometry_ok and topology_ok) else 1


if __name__ == '__main__':
    sys.exit(main())
