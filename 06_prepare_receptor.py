#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_prepare_receptor.py - Split a reference complex into receptor and peptide,
renumber it for HADDOCK, and extract the experimentally defined contact residues.

Written for 6N60 (Escherichia coli RNA polymerase bound to microcin J25) but it
works on any protein-peptide complex.

WHY THIS STEP EXISTS
--------------------
Docking a predicted lasso peptide against RNA polymerase only means something
if the receptor is prepared the same way every time and the binding site is
read from coordinates rather than from a figure. Four practical obstacles have
to be cleared first, and each of them silently produces a wrong result rather
than an error message.

  1. The deposited complex contains far more than the receptor. 6N60 carries
     nucleic acid (the transcription bubble), ions, waters and other
     heteroatoms. Leaving them in changes the surface the docking sees and the
     contacts the script would report.

  2. HADDOCK requires residue numbering to be UNIQUE ACROSS CHAINS. A
     multi-subunit polymerase numbers each subunit from 1, so the same residue
     number occurs in several chains. HADDOCK rejects the file, or worse,
     resolves a restraint to the wrong subunit. Two renumbering strategies are
     offered, and the mapping back to the original numbering is always written
     to disk so results can be reported in the numbering of the publication.

  3. The web server imposes an atom-count limit. A full polymerase exceeds it,
     so an optional truncation keeps only residues within a given radius of the
     peptide.

  4. The binding site must be measured. The receptor residues in contact with
     the microcin in the experimental structure ARE the binding site; that is
     the defensible source for HADDOCK active residues.

A NOTE ON TRUNCATION
--------------------
Truncating the receptor creates artificial chain termini and removes long-range
electrostatics. It is a legitimate practical compromise but it must be declared
among the limitations of the model. The script prints that reminder whenever
truncation is used, and reports exactly how many residues and atoms were
removed so the figure can be quoted in the methods section.

RELATIONSHIP TO 07_interface_residues.py
----------------------------------------
This script prepares files AND reports contacts, in the renumbered frame it
just produced. Script 07 only reports contacts, works on any complex, and adds
solvent-accessibility filtering. Use this one when preparing the receptor; use
07 when you only need the residue list or want the accessibility analysis.

USAGE
-----
  python3 06_prepare_receptor.py 6n60.pdb --list
  python3 06_prepare_receptor.py 6n60.pdb --peptide X
  python3 06_prepare_receptor.py 6n60.pdb --peptide X --truncate 15 --chains C,D

Requirements: numpy.
Licence: MIT
"""

import argparse
import sys

import numpy as np

# Residue names counted as amino acids, including selenomethionine and the
# AMBER protonation variants that appear in minimised structures.
AMINO_ACIDS = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
               'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
               'TYR', 'VAL', 'MSE', 'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH',
               'GLX', 'ASX'}
NUCLEOTIDES = {'DA', 'DT', 'DG', 'DC', 'DU', 'A', 'U', 'G', 'C', 'T'}
WATERS = {'HOH', 'WAT', 'DOD'}


def read_pdb(path):
    """Return atoms as dictionaries, preserving the original line.

    The raw line is kept so that the output files are byte-faithful to the
    input except for the fields deliberately rewritten. Rebuilding PDB lines
    from parsed values is a common source of subtle formatting corruption.
    """
    atoms = []
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('ENDMDL'):
                break                       # first model only
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            altloc = line[16]
            if altloc not in (' ', 'A'):
                continue
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]),
                                float(line[46:54])])
                resseq = int(line[22:26])
            except ValueError:
                continue
            atoms.append({
                'line': line.rstrip('\n'),
                'name': line[12:16].strip(),
                'resname': line[17:20].strip(),
                'chain': line[21],
                'resseq': resseq,
                'xyz': xyz,
                'hetatm': line.startswith('HETATM'),
            })
    if not atoms:
        sys.exit('No atoms were read from %s' % path)
    return atoms


def summarise_chains(atoms):
    """Classify each chain as protein, nucleic acid or heteroatoms.

    The classification is by majority rather than by any single residue,
    because deposited structures routinely carry a few modified residues or
    ligands inside an otherwise protein chain.
    """
    info = {}
    for a in atoms:
        c = info.setdefault(a['chain'],
                            {'residues': {}, 'aa': 0, 'nt': 0, 'other': 0})
        if a['resseq'] not in c['residues']:
            c['residues'][a['resseq']] = a['resname']
            if a['resname'] in AMINO_ACIDS:
                c['aa'] += 1
            elif a['resname'] in NUCLEOTIDES:
                c['nt'] += 1
            elif a['resname'] not in WATERS:
                c['other'] += 1
    for _, d in info.items():
        if d['aa'] >= max(1, d['nt'] + d['other']):
            d['type'] = 'protein'
        elif d['nt'] > 0:
            d['type'] = 'nucleic acid'
        else:
            d['type'] = 'heteroatoms / ions'
        d['n_residues'] = len(d['residues'])
    return info


def write_pdb(path, atoms):
    with open(path, 'w', encoding='utf-8') as fh:
        for a in atoms:
            fh.write(a['line'] + '\n')
        fh.write('END\n')


def main():
    ap = argparse.ArgumentParser(
        description='Split a reference complex, renumber the receptor for '
                    'HADDOCK and extract the contact residues.')
    ap.add_argument('pdb')
    ap.add_argument('--peptide', default=None,
                    help='peptide chain; defaults to the shortest protein chain')
    ap.add_argument('--cutoff', type=float, default=5.0,
                    help='contact distance in angstrom (5.0)')
    ap.add_argument('--list', action='store_true', dest='list_chains',
                    help='describe the chains and exit')
    ap.add_argument('--truncate', type=float, default=None, metavar='RADIUS',
                    help='keep only residues with an atom within RADIUS A of '
                         'the peptide (suggested 12-15). Off by default')
    ap.add_argument('--chains', default=None,
                    help='keep only these receptor chains, e.g. C,D')
    ap.add_argument('--renumber', choices=['no', 'offset', 'sequential'],
                    default='offset',
                    help='HADDOCK requires unique numbering across chains. '
                         '"offset" (default) adds a multiple of 1000 per chain, '
                         'preserving the original number in the last digits; '
                         '"sequential" renumbers 1..N; "no" leaves it alone')
    ap.add_argument('--receptor-out', default='receptor.pdb')
    ap.add_argument('--reference-out', default='reference_peptide.pdb')
    args = ap.parse_args()

    atoms = read_pdb(args.pdb)
    info = summarise_chains(atoms)

    print('=' * 70)
    print('  %s: %d chains' % (args.pdb, len(info)))
    print('=' * 70)
    for c in sorted(info):
        d = info[c]
        print('   chain %s  |  %5d residues  |  %s' % (c, d['n_residues'], d['type']))
    if args.list_chains:
        return 0

    # --- choose the peptide chain
    protein_chains = {c: d for c, d in info.items() if d['type'] == 'protein'}
    if not protein_chains:
        sys.exit('\nNo protein chain was identified.')
    if args.peptide:
        peptide_chain = args.peptide
        if peptide_chain not in info:
            sys.exit('\nChain %s does not exist. Chains: %s'
                     % (peptide_chain, ', '.join(sorted(info))))
    else:
        # The shortest protein chain is almost always the bound peptide in a
        # protein-peptide complex, but the guess is stated so it can be checked.
        peptide_chain = min(protein_chains,
                            key=lambda c: protein_chains[c]['n_residues'])
        print('\n  [peptide chain detected: %s, the shortest protein chain '
              'with %d residues]'
              % (peptide_chain, protein_chains[peptide_chain]['n_residues']))
        print('  If that is not correct, re-run with --peptide <chain>.')

    # --- split
    peptide_atoms = [x for x in atoms if x['chain'] == peptide_chain]
    receptor_atoms = [x for x in atoms
                      if x['chain'] != peptide_chain
                      and info[x['chain']]['type'] == 'protein'
                      and x['resname'] in AMINO_ACIDS
                      and not x['hetatm']]
    if args.chains:
        wanted = {c.strip() for c in args.chains.split(',') if c.strip()}
        missing = wanted - {x['chain'] for x in receptor_atoms}
        if missing:
            sys.exit('\nThese chains are not protein or do not exist: %s'
                     % ', '.join(sorted(missing)))
        receptor_atoms = [x for x in receptor_atoms if x['chain'] in wanted]
        print('\n  [receptor restricted to chains: %s]' % ', '.join(sorted(wanted)))
    if not receptor_atoms:
        sys.exit('\nThe receptor is empty. Inspect the chains with --list.')

    peptide_xyz = np.array([x['xyz'] for x in peptide_atoms])

    if args.truncate:
        atoms_before = len(receptor_atoms)
        residues_before = len({(x['chain'], x['resseq']) for x in receptor_atoms})
        near = set()
        receptor_xyz = np.array([x['xyz'] for x in receptor_atoms])
        # Blocked distance computation: the full pairwise matrix for a
        # polymerase-sized receptor would not fit in memory.
        for i in range(0, len(receptor_xyz), 4000):
            block = receptor_xyz[i:i + 4000]
            dmin = np.linalg.norm(block[:, None, :] - peptide_xyz[None, :, :],
                                  axis=2).min(axis=1)
            for j in np.where(dmin < args.truncate)[0]:
                at = receptor_atoms[i + j]
                near.add((at['chain'], at['resseq']))
        receptor_atoms = [x for x in receptor_atoms
                          if (x['chain'], x['resseq']) in near]
        if not receptor_atoms:
            sys.exit('\nTruncation at %s A left no residues.' % args.truncate)
        print('\n  [truncated at %s A from the peptide: %d -> %d residues, '
              '%d -> %d atoms]'
              % (args.truncate, residues_before, len(near),
                 atoms_before, len(receptor_atoms)))
        print('  Declare this truncation in the thesis: it creates artificial')
        print('  chain termini and belongs among the limitations of the model.')

    # --- unique numbering across chains (a HADDOCK requirement)
    mapping = {}
    receptor_chains = sorted({x['chain'] for x in receptor_atoms})
    has_duplicates = (len({(x['chain'], x['resseq']) for x in receptor_atoms})
                      != len({x['resseq'] for x in receptor_atoms}))
    if args.renumber != 'no' and len(receptor_chains) > 1 and has_duplicates:
        if args.renumber == 'offset':
            # Offsetting by 1000 per chain keeps the original number readable
            # in the last three digits, which makes manual cross-checking easy.
            for k, c in enumerate(receptor_chains):
                offset = k * 1000
                for resseq in sorted({x['resseq'] for x in receptor_atoms
                                      if x['chain'] == c}):
                    mapping[(c, resseq)] = resseq + offset
        else:
            n = 0
            for c in receptor_chains:
                for resseq in sorted({x['resseq'] for x in receptor_atoms
                                      if x['chain'] == c}):
                    n += 1
                    mapping[(c, resseq)] = n
        highest = max(mapping.values())
        if highest > 9999:
            # The PDB residue-number field is four columns wide; a larger value
            # would overflow into the neighbouring field and corrupt the file.
            sys.exit('\nRenumbering reaches %d, which does not fit the PDB '
                     'format.\nUse --truncate to reduce the receptor, or '
                     '--renumber sequential.' % highest)
        for x in receptor_atoms:
            new = mapping[(x['chain'], x['resseq'])]
            buf = list(x['line'].ljust(80))
            buf[22:26] = list('%4d' % new)
            x['line'] = ''.join(buf).rstrip()
            x['new_resseq'] = new
        with open('receptor_numbering.tsv', 'w', encoding='utf-8') as fh:
            fh.write('chain\toriginal_residue\thaddock_residue\n')
            for (c, old), new in sorted(mapping.items()):
                fh.write('%s\t%d\t%d\n' % (c, old, new))
        print('\n  [renumbered for HADDOCK: mode "%s", %d residues, max %d]'
              % (args.renumber, len(mapping), highest))
        print('  Mapping to the original numbering -> receptor_numbering.tsv')
        print('  Report results in the ORIGINAL numbering using that table.')
    else:
        for x in receptor_atoms:
            x['new_resseq'] = x['resseq']
        if len(receptor_chains) > 1 and has_duplicates:
            print('\n  WARNING: residue numbers are duplicated across chains. '
                  'HADDOCK will reject the file.')

    write_pdb(args.reference_out, peptide_atoms)
    write_pdb(args.receptor_out, receptor_atoms)
    print('\n  receptor   -> %s   (%d atoms, %d residues)'
          % (args.receptor_out, len(receptor_atoms),
             len({(x['chain'], x['resseq']) for x in receptor_atoms})))
    print('  reference  -> %s   (%d atoms, %d residues)'
          % (args.reference_out, len(peptide_atoms),
             len({x['resseq'] for x in peptide_atoms})))
    print('  Nucleic acids, waters, ions and heteroatoms were discarded.')
    if len(receptor_atoms) > 50000:
        print('\n  WARNING: %d atoms exceed the 50,000-atom limit of the web '
              'server.' % len(receptor_atoms))
        print('           Use --truncate 15 or --chains to reduce it.')
    else:
        print('  Receptor size: %d atoms (server limit: 50,000).'
              % len(receptor_atoms))

    # --- contacts
    receptor_xyz = np.array([x['xyz'] for x in receptor_atoms])
    contacts = {}
    for i in range(0, len(receptor_xyz), 4000):        # blocked, see above
        block = receptor_xyz[i:i + 4000]
        d = np.linalg.norm(block[:, None, :] - peptide_xyz[None, :, :], axis=2)
        for j in np.where(d.min(axis=1) < args.cutoff)[0]:
            at = receptor_atoms[i + j]
            key = (at['chain'], at['resseq'])
            dmin = float(d[j].min())
            if key not in contacts or dmin < contacts[key][2]:
                contacts[key] = (at['resname'], at['new_resseq'], dmin)

    print('\n' + '=' * 70)
    print('  RECEPTOR ACTIVE RESIDUES  (within %.1f A of chain %s)'
          % (args.cutoff, peptide_chain))
    print('=' * 70)
    if not contacts:
        print('  None. Check that the peptide chain is the right one.')
        return 1

    by_chain = {}
    for (c, resseq), (resname, new, dmin) in sorted(contacts.items()):
        by_chain.setdefault(c, []).append((resseq, new, resname, dmin))

    renumbered = any(k[1] != v[1] for k, v in contacts.items())
    for c in sorted(by_chain):
        entries = by_chain[c]
        print('\n  Chain %s  -  %d residues' % (c, len(entries)))
        header = ('      residue     original   haddock   distance' if renumbered
                  else '      residue     original   distance')
        print(header)
        for old, new, resname, dmin in entries:
            if renumbered:
                print('      %-4s       %6d    %6d    %5.2f A'
                      % (resname, old, new, dmin))
            else:
                print('      %-4s       %6d    %5.2f A' % (resname, old, dmin))

    print('\n' + '-' * 70)
    print('  ACTIVE RESIDUES TO PASTE INTO THE FORM')
    print('  (HADDOCK numbering, all chains in a single list)')
    print('-' * 70)
    print('      ' + ','.join(str(v[1]) for v in sorted(contacts.values(),
                                                        key=lambda v: v[1])))

    print('\n' + '=' * 70)
    print('  Review the list before using it: a 5 A cutoff can include')
    print('  marginal contacts. For a stricter definition use --cutoff 4.0.')
    print('  Do NOT declare active residues on the peptide side; see')
    print('  07_interface_residues.py for the reasoning.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
