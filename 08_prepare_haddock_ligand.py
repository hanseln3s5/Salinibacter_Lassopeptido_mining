#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_prepare_haddock_ligand.py - Turn a predicted lasso peptide into a HADDOCK-
ready ligand: standard residue naming, a fixed segid, an ensemble PDB, and the
unambiguous restraint that keeps the macrolactam closed during docking.

WHY THIS STEP EXISTS
--------------------
Four separate things block a predicted lasso peptide from being docked as-is.

  1. AMBER residue naming. Structures that have been minimised come back with
     names such as GLH, ASH, HIE/HID/HIP, CYX and the terminal forms NGLU or
     CALA. HADDOCK's engine is CNS, which does not recognise them and fails
     during topology generation. The failure message points at the topology
     file rather than the residue name, so the cause is not obvious.

  2. Chain and segid. HADDOCK identifies molecules by segid. If the ligand does
     not carry the segid used in the restraint file, the restraint silently
     refers to nothing and the macrolactam is unrestrained.

  3. Ensemble input. HADDOCK can sample several starting conformers of the same
     molecule if they are supplied as MODEL/ENDMDL blocks in one file. For a
     peptide whose tail is the least reliable part of the prediction, giving
     the docking more than one starting conformation is the honest way to
     handle that uncertainty rather than committing to a single model.

  4. The isopeptide bond. CNS builds topology from residue connectivity and
     does not know about the non-canonical macrolactam. Without an explicit
     distance restraint, the ring opens during the simulated annealing and the
     molecule being docked is no longer a lasso peptide. The restraint written
     here pins the bond at 1.33 A, the characteristic length of an amide bond.

WHY THE RESTRAINT IS "UNAMBIGUOUS"
----------------------------------
HADDOCK distinguishes ambiguous interaction restraints, which drive the two
molecules together and are allowed to be partially violated, from unambiguous
restraints, which encode known covalent or experimental geometry and are
enforced throughout. The macrolactam is a covalent bond, so it belongs in the
second category and goes in the "Unambiguous distance restraints" field, never
in the AIR field.

USAGE
-----
  python3 08_prepare_haddock_ligand.py min1.pdb min2.pdb min3.pdb
  python3 08_prepare_haddock_ligand.py min*.pdb --segid B \
        --out peptide.pdb --tbl isopeptide.tbl
  python3 08_prepare_haddock_ligand.py min1.pdb --keep-hydrogens

Hydrogens are stripped by default because HADDOCK regenerates them with its own
topology; keeping predictor-assigned hydrogens can conflict with the protonation
state CNS assigns.

No dependencies: Python 3.8+ standard library only.
Licence: MIT
"""

import argparse
import os
import sys

# AMBER naming -> standard name understood by CNS.
TRANSLATE = {
    'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',
    'CYX': 'CYS', 'CYM': 'CYS',
    'ASH': 'ASP', 'ASX': 'ASP',
    'GLH': 'GLU', 'GLX': 'GLU',
    'LYN': 'LYS', 'ARN': 'ARG',
}
STANDARD = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL'}
# Side-chain carboxyl carbon that can close the macrolactam.
CARBOXYL = {'ASP': 'CG', 'GLU': 'CD'}


def normalise(resname):
    """Return (standard_name, was_changed)."""
    r = resname.strip().upper()
    if r in TRANSLATE:
        return TRANSLATE[r], True
    # Four-letter terminal forms such as NGLU or CALA.
    if len(r) == 4 and r[0] in 'NC' and r[1:] in STANDARD:
        return r[1:], True
    return r, False


def process(path, segid, strip_hydrogens):
    """Return (atom_lines, renamings, residues) for one PDB."""
    lines, renamings, residues = [], {}, {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            if not raw.startswith(('ATOM  ', 'HETATM')):
                continue
            # Pad to 80 columns: many predictors truncate the line, and the
            # segid field at columns 73-76 would otherwise not exist to write to.
            raw = raw.rstrip('\n').ljust(80)
            element = raw[76:78].strip() or raw[12:16].strip()[:1]
            if strip_hydrogens and element.upper().startswith('H'):
                continue
            resname = raw[17:20]
            new_name, changed = normalise(resname)
            if changed:
                renamings[resname.strip()] = new_name
            try:
                resseq = int(raw[22:26])
            except ValueError:
                continue
            residues.setdefault(resseq, {'resname': new_name, 'atoms': set()})
            residues[resseq]['atoms'].add(raw[12:16].strip())
            # Rewrite residue name, chain identifier and segid in place, by
            # column, because the PDB format is column-based.
            buf = list(raw)
            buf[17:20] = list('%3s' % new_name)
            buf[21] = segid
            buf[72:76] = list('%-4s' % segid)
            lines.append(''.join(buf).rstrip())
    if not lines:
        sys.exit('No atoms were read from %s' % path)
    return lines, renamings, residues


def find_acceptors(residues):
    """Asp/Glu residues whose side-chain carboxyl is present, in sequence order."""
    candidates = []
    for resseq in sorted(residues):
        resname = residues[resseq]['resname']
        atom = CARBOXYL.get(resname)
        if atom and atom in residues[resseq]['atoms']:
            candidates.append((resseq, resname, atom))
    return candidates


def main():
    ap = argparse.ArgumentParser(
        description='Prepare a lasso peptide as a HADDOCK ligand: naming, '
                    'segid, ensemble PDB and the isopeptide restraint.')
    ap.add_argument('pdbs', nargs='+',
                    help='one or more conformers (min1.pdb min2.pdb ...)')
    ap.add_argument('--segid', default='B',
                    help='chain and segid to assign (default: B, the second '
                         'molecule of the docking)')
    ap.add_argument('--out', default='peptide_haddock.pdb',
                    help='output ensemble PDB')
    ap.add_argument('--tbl', default='isopeptide.tbl',
                    help='output restraint file')
    ap.add_argument('--acceptor', type=int, default=None,
                    help='force the macrolactam acceptor residue')
    ap.add_argument('--keep-hydrogens', action='store_true',
                    help='do not strip hydrogens (they are removed by default '
                         'because HADDOCK regenerates them)')
    args = ap.parse_args()

    if len(args.segid) != 1:
        sys.exit('The segid must be a single character (A, B, ...).')

    all_renamings, reference_residues, blocks = {}, None, []
    for path in args.pdbs:
        if not os.path.isfile(path):
            sys.exit('File not found: %s' % path)
        lines, renamings, residues = process(path, args.segid,
                                             not args.keep_hydrogens)
        all_renamings.update(renamings)
        if reference_residues is None:
            reference_residues = residues
        elif sorted(reference_residues) != sorted(residues):
            # HADDOCK requires every member of an ensemble to have identical
            # composition and numbering; otherwise restraints resolve to
            # different atoms in different conformers.
            print('  WARNING: %s does not share the residue numbering of %s.'
                  % (path, args.pdbs[0]))
            print('           HADDOCK requires all ensemble members to match.')
        blocks.append((path, lines))

    # ---- ensemble PDB
    with open(args.out, 'w', encoding='utf-8') as fh:
        for i, (_, lines) in enumerate(blocks, 1):
            fh.write('MODEL     %4d\n' % i)
            fh.write('\n'.join(lines) + '\n')
            fh.write('ENDMDL\n')
        fh.write('END\n')

    print('=' * 68)
    print('  Ensemble PDB   : %s  (%d model(s))' % (args.out, len(blocks)))
    print('  chain / segid  : %s' % args.segid)
    print('  hydrogens      : %s'
          % ('kept' if args.keep_hydrogens
             else 'stripped (HADDOCK regenerates them)'))
    if all_renamings:
        print('  residues renamed:')
        for k, v in sorted(all_renamings.items()):
            print('      %s -> %s' % (k, v))
    else:
        print('  residues renamed: none (already standard)')

    # ---- acceptor and restraint
    candidates = find_acceptors(reference_residues)
    if not candidates:
        sys.exit('\nNo Asp/Glu with a side-chain carboxyl was found. '
                 'Check the input file.')
    if args.acceptor is not None:
        chosen = [c for c in candidates if c[0] == args.acceptor]
        if not chosen:
            listing = ', '.join('%s%d' % (r, i) for i, r, _ in candidates)
            sys.exit('\nResidue %d is not a valid acceptor. Candidates: %s'
                     % (args.acceptor, listing))
        resseq, resname, atom = chosen[0]
    else:
        resseq, resname, atom = candidates[0]
        if len(candidates) > 1:
            others = ', '.join('%s%d' % (r, i) for i, r, _ in candidates[1:])
            print('  WARNING: more than one acceptor candidate. %s%d was chosen;'
                  % (resname, resseq))
            print('           the others are %s.' % others)
            print('           If that is not the right one, re-run with '
                  '--acceptor N.')

    first = min(reference_residues)
    if 'N' not in reference_residues[first]['atoms']:
        sys.exit('\nResidue %d has no N atom; the bond cannot be defined.' % first)

    tbl = ('! Lasso peptide isopeptide bond\n'
           '! N-terminus (%s%d) -- %s of %s%d\n'
           '! Target distance 1.33 A, the characteristic length of an amide bond\n'
           '! This is a COVALENT bond, so it belongs in the HADDOCK field\n'
           '! "Unambiguous distance restraints", not in the AIR field.\n'
           'assign (segid %s and resid %d and name N) '
           '(segid %s and resid %d and name %s) 1.33 0.1 0.1\n'
           % (reference_residues[first]['resname'], first, atom, resname, resseq,
              args.segid, first, args.segid, resseq, atom))
    with open(args.tbl, 'w', encoding='utf-8') as fh:
        fh.write(tbl)

    print('-' * 68)
    print('  Restraint      : %s' % args.tbl)
    print()
    print(tbl.rstrip())
    print('=' * 68)
    print('\nNext step: upload the PDB as molecule 2 and the .tbl in the')
    print('"Unambiguous distance restraints" field. Verify that molecule 2')
    print('receives segid %s; if HADDOCK assigns a different one, re-run with'
          % args.segid)
    print('--segid to match, or the restraint will refer to nothing.')


if __name__ == '__main__':
    main()
