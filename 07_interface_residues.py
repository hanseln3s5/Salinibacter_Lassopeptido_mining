#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07_interface_residues.py - Extract interface residues from a reference complex
and print them in the format the HADDOCK web server expects.

WHY THIS STEP EXISTS
--------------------
HADDOCK drives docking with ambiguous interaction restraints built from user-
supplied "active" residues. Those must be the residues of the REAL binding
site, numbered AS THEY APPEAR IN THE FILE THAT WAS UPLOADED. Two things go
wrong routinely:

  1. Residue lists taken from the literature use the numbering of the original
     publication. A PDB that has been renumbered, had chains merged, or been
     processed by a preparation pipeline no longer matches, and the restraints
     then point at the wrong residues without any error being raised.

  2. The binding site is identified by eye from a figure rather than measured
     from coordinates, which is neither reproducible nor auditable.

This script reads the residue list from the structure itself. For the lasso
peptide case the reference complex is 6N60, RNA polymerase bound to microcin
J25: the polymerase residues contacting the microcin are the experimentally
defined binding site, and those are what belong in the active list.

TWO TRAPS IN THE HADDOCK FORM THIS SCRIPT ADDRESSES
---------------------------------------------------
  "Remove buried active residues", default RSA threshold 15%.
  The RNA polymerase secondary channel is a tunnel; its lining residues are
  buried by definition. Leaving that option enabled with a 15% threshold
  deletes precisely the residues that define the site. Run with --sasa and the
  script reports which residues would be removed, so the decision is informed
  rather than accidental.

  Multi-chain receptors.
  The secondary channel is formed between subunits, so contacts fall on more
  than one chain. HADDOCK takes one list per molecule, so the script says so
  explicitly and prints the list per chain rather than silently merging them
  into a single meaningless list.

A NOTE ON THE LIGAND SIDE
-------------------------
Active residues should normally NOT be declared on the peptide. Which face of
the peptide contacts the target is unknown; imposing active residues there
dictates the orientation in advance and the docking then confirms an assumption
rather than testing it. Declaring the peptide's surface residues as passive and
letting HADDOCK resolve the orientation is the defensible setup. The script
prints this reminder with its output.

USAGE
-----
  python3 07_interface_residues.py 6n60.pdb --list
  python3 07_interface_residues.py 6n60.pdb --ligand X --cutoff 5.0
  python3 07_interface_residues.py 6n60.pdb --ligand X --receptor C D --sasa

No dependencies: Python 3.8+ standard library only.
Licence: MIT
"""

import argparse
import math
import sys
from collections import OrderedDict, defaultdict

# Van der Waals radii (A). Used both for contact detection and for SASA.
VDW_RADIUS = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
              'P': 1.80, 'H': 1.20, 'F': 1.47, 'SE': 1.90}
VDW_DEFAULT = 1.70

# Maximum accessible surface area per residue (Tien et al. 2013, theoretical).
ASA_MAX = {'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP': 193, 'CYS': 167,
           'GLU': 223, 'GLN': 225, 'GLY': 104, 'HIS': 224, 'ILE': 197,
           'LEU': 201, 'LYS': 236, 'MET': 224, 'PHE': 240, 'PRO': 159,
           'SER': 155, 'THR': 172, 'TRP': 285, 'TYR': 263, 'VAL': 174}

THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

WATERS = {'HOH', 'WAT', 'DOD'}


# --------------------------------------------------------------------- input
class Atom:
    __slots__ = ('chain', 'resseq', 'icode', 'resname', 'name', 'element',
                 'x', 'y', 'z', 'radius', 'hetatm')

    def __init__(self, line):
        self.hetatm = line.startswith('HETATM')
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21].strip() or '_'
        self.resseq = int(line[22:26])
        self.icode = line[26].strip()
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        e = line[76:78].strip().upper()
        if not e:
            e = ''.join(c for c in self.name if c.isalpha())[:1].upper()
        self.element = e
        self.radius = VDW_RADIUS.get(e, VDW_DEFAULT)

    @property
    def key(self):
        return (self.chain, self.resseq, self.icode)


def read_pdb(path, keep_hydrogens=False, keep_waters=False):
    """Parse the first model only, one conformer per atom."""
    atoms = []
    model = 0
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            if line.startswith('MODEL'):
                model += 1
                if model > 1:
                    break                       # first model only
            if line.startswith('ENDMDL'):
                break
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            if len(line) < 54:
                continue
            altloc = line[16].strip()
            if altloc and altloc not in ('A', '1'):
                continue                        # a single alternate conformer
            try:
                a = Atom(line)
            except ValueError:
                continue
            if not keep_waters and a.resname in WATERS:
                continue
            if not keep_hydrogens and a.element == 'H':
                continue
            atoms.append(a)
    return atoms


def by_chain(atoms):
    d = OrderedDict()
    for a in atoms:
        d.setdefault(a.chain, []).append(a)
    return d


# ---------------------------------------------------------------- neighbour grid
class Grid:
    """Uniform cell list, so contacts are found without an all-versus-all scan.

    A ribosome- or polymerase-sized receptor has tens of thousands of atoms;
    the quadratic comparison would take minutes, the grid takes a second.
    """

    def __init__(self, atoms, cell):
        self.cell = cell
        self.cells = defaultdict(list)
        for a in atoms:
            self.cells[self._index(a.x, a.y, a.z)].append(a)

    def _index(self, x, y, z):
        c = self.cell
        return (int(math.floor(x / c)), int(math.floor(y / c)),
                int(math.floor(z / c)))

    def near(self, x, y, z):
        i, j, k = self._index(x, y, z)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for dk in (-1, 0, 1):
                    for a in self.cells.get((i + di, j + dj, k + dk), ()):
                        yield a


# -------------------------------------------------------------------- contacts
def contacting_residues(receptor, ligand, cutoff):
    """Receptor residues with any atom within `cutoff` of any ligand atom.

    Returns {(chain, resseq, icode): (resname, minimum distance)}. The minimum
    distance is kept so the report can be ordered by proximity, which helps
    when a borderline cutoff has to be justified.
    """
    grid = Grid(receptor, max(cutoff, 4.0))
    cutoff2 = cutoff * cutoff
    found = {}
    for b in ligand:
        for a in grid.near(b.x, b.y, b.z):
            d2 = (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
            if d2 <= cutoff2:
                d = math.sqrt(d2)
                key = a.key
                if key not in found or d < found[key][1]:
                    found[key] = (a.resname, d)
    return found


# ------------------------------------------------------------------------ SASA
def sphere_points(n):
    """Quasi-uniform points on the unit sphere (Fibonacci spiral)."""
    points = []
    phi = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = 1 - (2.0 * i + 1.0) / n
        r = math.sqrt(max(0.0, 1 - y * y))
        t = phi * i
        points.append((math.cos(t) * r, y, math.sin(t) * r))
    return points


def sasa_per_residue(atoms, probe=1.4, n_points=92):
    """Shrake-Rupley SASA on the receptor ALONE.

    Computing it on the receptor without the ligand is deliberate: the question
    is whether a residue is accessible in the unbound receptor, which is the
    state HADDOCK evaluates when it decides whether an active residue is
    buried. Including the ligand would mask the very residues of interest.

    n_points is lower than in step 04 because receptors here are large; 92 is
    enough to classify residues as exposed or buried.
    """
    points = sphere_points(n_points)
    max_radius = max(a.radius for a in atoms) + probe
    grid = Grid(atoms, 2 * max_radius)
    area = defaultdict(float)
    resnames = {}
    for a in atoms:
        ra = a.radius + probe
        neighbours = []
        for b in grid.near(a.x, a.y, a.z):
            if b is a:
                continue
            rb = b.radius + probe
            limit = ra + rb
            d2 = (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
            if d2 < limit * limit:
                neighbours.append((b, rb * rb))
        exposed = 0
        for px, py, pz in points:
            sx, sy, sz = a.x + ra * px, a.y + ra * py, a.z + ra * pz
            covered = False
            for b, rb2 in neighbours:
                if (sx - b.x) ** 2 + (sy - b.y) ** 2 + (sz - b.z) ** 2 < rb2:
                    covered = True
                    break
            if not covered:
                exposed += 1
        area[a.key] += 4.0 * math.pi * ra * ra * exposed / n_points
        resnames[a.key] = a.resname
    return {k: (resnames[k], v) for k, v in area.items()}


def format_list(keys):
    """Comma-separated residue numbers, ordered, as the HADDOCK form expects."""
    return ','.join(str(r) for _, r, _ in sorted(keys, key=lambda k: k[1]))


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description='Interface residues in HADDOCK format, read from your own PDB.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n'
               '  python3 07_interface_residues.py 6n60.pdb --list\n'
               '  python3 07_interface_residues.py 6n60.pdb --ligand X --sasa\n')
    ap.add_argument('pdb', help='reference complex, e.g. 6N60')
    ap.add_argument('--list', action='store_true', dest='list_chains',
                    help='print the chains in the file and exit')
    ap.add_argument('--ligand', '-l', metavar='CHAIN',
                    help='ligand chain, e.g. microcin J25 in 6N60')
    ap.add_argument('--receptor', '-r', nargs='+', metavar='CHAIN',
                    help='receptor chains (default: everything but the ligand)')
    ap.add_argument('--cutoff', '-c', type=float, default=5.0,
                    help='contact distance in angstrom (5.0). HADDOCK '
                         'documentation commonly uses 3.9-5.0')
    ap.add_argument('--sasa', action='store_true',
                    help='compute accessibility and flag the residues HADDOCK '
                         'would discard as buried (slow on large receptors)')
    ap.add_argument('--rsa-min', type=float, default=15.0,
                    help="HADDOCK's relative accessibility threshold, %% (15.0)")
    ap.add_argument('--keep-waters', action='store_true')
    args = ap.parse_args()

    atoms = read_pdb(args.pdb, keep_waters=args.keep_waters)
    if not atoms:
        sys.exit('No atoms were read from %s. Check that it is a valid PDB.'
                 % args.pdb)
    chains = by_chain(atoms)

    print('\nFile   : %s' % args.pdb)
    print('Atoms  : %d   Chains: %d\n' % (len(atoms), len(chains)))
    print('  %-8s %8s %8s   %s' % ('chain', 'atoms', 'residues', 'first residues'))
    for c, group in chains.items():
        residues = OrderedDict()
        for a in group:
            residues[a.key] = a.resname
        sample = ' '.join('%s%d' % (THREE_TO_ONE.get(rn, 'X'), k[1])
                          for k, rn in list(residues.items())[:6])
        print('  %-8s %8d %8d   %s' % (c, len(group), len(residues), sample))

    if args.list_chains or not args.ligand:
        if not args.list_chains:
            print('\n[!] --ligand is missing. Look at the table above and pick')
            print('    the small peptide chain; it is usually the one with the')
            print('    fewest residues.')
        return

    if args.ligand not in chains:
        sys.exit('\nChain "%s" is not present in the file.' % args.ligand)

    ligand = chains[args.ligand]
    receptor_chains = args.receptor or [c for c in chains if c != args.ligand]
    missing = [c for c in receptor_chains if c not in chains]
    if missing:
        sys.exit('\nReceptor chains not found: %s' % ', '.join(missing))
    receptor = [a for c in receptor_chains for a in chains[c]]

    print('\n' + '=' * 70)
    print('Ligand   : chain %s (%d atoms)' % (args.ligand, len(ligand)))
    print('Receptor : chains %s (%d atoms)'
          % (', '.join(receptor_chains), len(receptor)))
    print('Cutoff   : %.1f A' % args.cutoff)
    print('=' * 70)

    contacts = contacting_residues(receptor, ligand, args.cutoff)
    if not contacts:
        print('\nNo contacts at %.1f A. Try a larger cutoff.' % args.cutoff)
        return

    rsa = {}
    if args.sasa:
        print('\nComputing accessibility (Shrake-Rupley, receptor alone)...')
        sys.stdout.flush()
        areas = sasa_per_residue(receptor)
        for key, (resname, area) in areas.items():
            maximum = ASA_MAX.get(resname)
            if maximum:
                rsa[key] = 100.0 * area / maximum

    grouped = defaultdict(list)
    for key, (resname, d) in contacts.items():
        grouped[key[0]].append((key, resname, d))

    print('\n--- DETAIL BY CHAIN ---')
    for c in sorted(grouped):
        entries = sorted(grouped[c], key=lambda t: t[0][1])
        print('\nChain %s  (%d contacting residues)' % (c, len(entries)))
        for key, resname, d in entries:
            flag = ''
            if key in rsa:
                flag = ('  RSA %5.1f%%%s'
                        % (rsa[key], '  <-- HADDOCK would discard this'
                           if rsa[key] < args.rsa_min else ''))
            print('    %s%-5d %-4s  d=%.2f A%s'
                  % (THREE_TO_ONE.get(resname, 'X'), key[1], resname, d, flag))

    print('\n\n' + '=' * 70)
    print('TO PASTE INTO THE HADDOCK FORM')
    print('=' * 70)
    print('\nMolecule 1 (receptor) -> "Active residues":\n')
    if len(grouped) == 1:
        c = list(grouped)[0]
        print('  %s' % format_list([k for k, _, _ in grouped[c]]))
    else:
        print('  [NOTE] the interface spans several chains. HADDOCK takes ONE')
        print('  list per molecule, so either upload each chain as a separate')
        print('  molecule, or merge the chains into one renumbered file before')
        print('  uploading, and renumber this list to match.\n')
        for c in sorted(grouped):
            print('  chain %s:  %s' % (c, format_list([k for k, _, _ in grouped[c]])))

    if rsa:
        buried = [k for k, _, _ in (t for lst in grouped.values() for t in lst)
                  if k in rsa and rsa[k] < args.rsa_min]
        if buried:
            print('\n  %d of these fall below RSA %.0f%%.' % (len(buried), args.rsa_min))
            print('  For a deep channel or pocket that is NORMAL, so UNTICK')
            print('  "Remove buried active residues" in the form, or the binding')
            print('  site itself will be deleted from the restraints.')
            print('  Accessible-only list (RSA >= %.0f%%), if you prefer it:'
                  % args.rsa_min)
            accessible = [k for k, _, _ in
                          (t for lst in grouped.values() for t in lst)
                          if k not in rsa or rsa[k] >= args.rsa_min]
            print('    %s' % (format_list(accessible) or '(none)'))

    print('\nMolecule 1 -> "Passive residues":')
    print('  leave empty and tick "Automatically define passive residues')
    print('  around the active residues"')
    print('\nMolecule 2 (your peptide) -> "Active residues":')
    print('  LEAVE EMPTY. Tick "Automatically define surface residues as')
    print('  passive". You do not know which face of the peptide binds; if you')
    print('  impose active residues there you are dictating the orientation to')
    print('  the docking instead of testing it.')
    print('=' * 70 + '\n')


if __name__ == '__main__':
    main()
