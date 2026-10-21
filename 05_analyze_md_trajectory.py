#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_analyze_md_trajectory.py - Molecular dynamics analysis of a lasso peptide,
including the one check no general-purpose tool performs: whether the threaded
topology survives the simulation.

WHY THIS STEP EXISTS
--------------------
RMSD and radius of gyration report whether a structure moves. They do not
report whether it is still a lasso peptide. The two are independent:

  - A model can show a very low RMSD and have lost its topology, because the
    tail can slip out of the ring with little overall coordinate change in a
    small peptide.
  - A model can keep its topology while showing a high RMSD, because the tail
    of a lasso peptide is intrinsically flexible.

So the topology has to be measured directly. This script recomputes, frame by
frame, how many times the tail crosses the plane of the ring, using the same
Moller-Trumbore ray-triangle test as steps 03 and 04. An odd count means the
fold is intact.

WHAT THE RESULT MEANS
---------------------
If the threaded fraction stays at essentially 100% over the trajectory, the
steric plug retains the tail under the force field and the topology produced by
the structure predictor is not an artefact of the prediction. That is a
substantially stronger statement than a low RMSD.

If the fraction decays over time, the plug does not hold and the model was
unstable. That is also a reportable result, and the script prints the time of
the first unthreaded frame so the event can be located in the trajectory.

CONTROLS WORTH RUNNING ALONGSIDE
--------------------------------
  1. The LINEAR peptide of the same sequence, without the macrolactam, run with
     --linear. If threading contributes rigidity it should show up as a
     difference in per-residue RMSF.

  2. The model WITH and WITHOUT a disulfide bridge, when the core carries
     compatible cysteines. This tests directly whether the bridge contributes
     to stability, which is what separates classes I, III and IV from class II.

A CAVEAT THAT BELONGS IN THE METHODS SECTION
--------------------------------------------
Molecular dynamics tests whether a model is STABLE under a force field. It does
not test whether the model is CORRECT. An incorrect structure can be perfectly
stable. The simulation supports the plausibility of the model; it does not
validate the prediction.

A PRACTICAL WARNING ABOUT THE TOPOLOGY FILE
-------------------------------------------
The isopeptide bond is chemically an ordinary amide and standard force fields
describe it well, but GROMACS does not detect it automatically: it must be
declared when the system is built. CHARMM-GUI handles this crosslink correctly.
If the bond is missing from the topology the macrolactam simply opens during
the simulation and every subsequent number is meaningless. This script tracks
the bond length for exactly that reason and warns when it exceeds 2.0 A.

USAGE
-----
  python3 05_analyze_md_trajectory.py topology.gro trajectory.xtc --ring 1:9
  python3 05_analyze_md_trajectory.py md.pdb md.dcd --ring 1:9 --stride 5
  python3 05_analyze_md_trajectory.py md.gro md.xtc --ring 1:9 \
        --isopeptide N:1:CG:9

OUTPUT (in the directory given by -o, default 'md_analysis/')
  time_series.tsv       per frame: time, RMSD, Rg, crossings, threaded,
                        isopeptide bond length
  rmsf_per_residue.tsv  mean fluctuation per residue, with ring/tail zone
  summary.txt           statistics and the verdict on topology

Requirements: MDAnalysis and NumPy, which will already be in an MD environment.
  conda install -c conda-forge mdanalysis

Licence: MIT
"""

import argparse
import os
import sys

try:
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
except ImportError as exc:
    sys.exit('Missing dependency (%s).\n'
             'Install with:  conda install -c conda-forge mdanalysis' % exc)


# -------------------------------------------------------- vectorised geometry
def ring_crossings(ring_ca, tail_ca):
    """Count intersections of consecutive tail segments with the ring surface.

    The ring has no intrinsic surface, so one is built as a triangle fan from
    the centroid of the ring CA atoms. This is the same construction used in
    steps 03 and 04, so the three scripts agree by design.

    At most one crossing is counted per segment: a segment passing exactly
    through a shared triangle edge would otherwise be counted twice and flip
    the parity, turning a threaded frame into an apparently unthreaded one.
    """
    if len(ring_ca) < 3 or len(tail_ca) < 2:
        return 0
    centre = ring_ca.mean(axis=0)
    triangles = [(centre, ring_ca[i], ring_ca[(i + 1) % len(ring_ca)])
                 for i in range(len(ring_ca))]
    count = 0
    for origin, nxt in zip(tail_ca[:-1], tail_ca[1:]):
        d = nxt - origin
        for v0, v1, v2 in triangles:
            e1, e2 = v1 - v0, v2 - v0
            h = np.cross(d, e2)
            det = np.dot(e1, h)
            if abs(det) < 1e-9:
                continue                       # segment parallel to triangle
            inv = 1.0 / det
            s = origin - v0
            u = inv * np.dot(s, h)
            if u < 0.0 or u > 1.0:
                continue
            q = np.cross(s, e1)
            v = inv * np.dot(d, q)
            if v < 0.0 or u + v > 1.0:
                continue
            t = inv * np.dot(e2, q)
            if 0.0 <= t <= 1.0:                # within the segment itself
                count += 1
                break
    return count


def main():
    ap = argparse.ArgumentParser(
        description='MD trajectory analysis of a lasso peptide, with frame-by-'
                    'frame verification of the threaded topology.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n'
               '  python3 05_analyze_md_trajectory.py md.gro md.xtc --ring 1:9\n')
    ap.add_argument('topology', help='gro, pdb, psf...')
    ap.add_argument('trajectory', help='xtc, dcd, trr...')
    ap.add_argument('--ring', required=True, metavar='R1:RACC',
                    help='first and last residue of the macrolactam, e.g. 1:9')
    ap.add_argument('--selection', default='protein',
                    help='MDAnalysis selection for the peptide (protein)')
    ap.add_argument('--isopeptide', metavar='AT1:R1:AT2:R2',
                    help='atoms of the isopeptide bond to monitor, e.g. N:1:CG:9')
    ap.add_argument('--stride', type=int, default=1,
                    help='analyse one frame in every N (1)')
    ap.add_argument('--linear', action='store_true',
                    help='linear control: skip the topology check')
    ap.add_argument('-o', '--outdir', default='md_analysis')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    u = mda.Universe(args.topology, args.trajectory)
    peptide = u.select_atoms(args.selection)
    if not len(peptide):
        sys.exit('Selection "%s" returned no atoms.' % args.selection)

    r1, racc = (int(x) for x in args.ring.split(':'))
    ring_sel = peptide.select_atoms('name CA and resid %d:%d' % (r1, racc))
    tail_sel = peptide.select_atoms('name CA and resid %d:%d'
                                    % (racc + 1, peptide.residues.resids.max()))
    backbone = peptide.select_atoms('backbone')
    if not len(backbone):
        backbone = peptide.select_atoms('name CA')

    isopeptide_pair = None
    if args.isopeptide:
        a1, i1, a2, i2 = args.isopeptide.split(':')
        g1 = peptide.select_atoms('name %s and resid %s' % (a1, i1))
        g2 = peptide.select_atoms('name %s and resid %s' % (a2, i2))
        if len(g1) == 1 and len(g2) == 1:
            isopeptide_pair = (g1, g2)
        else:
            print('[warning] could not locate the isopeptide bond atoms; '
                  'bond length will not be monitored.')

    print('Peptide : %d atoms, %d residues' % (len(peptide), len(peptide.residues)))
    print('Ring    : residues %d-%d (%d CA)' % (r1, racc, len(ring_sel)))
    print('Tail    : %d CA' % len(tail_sel))
    print('Frames  : %d (analysing one in every %d)\n'
          % (len(u.trajectory), args.stride))

    # The first frame is the RMSD reference. Using the starting model rather
    # than an average structure makes the curve directly interpretable as
    # "how far has it drifted from the prediction".
    u.trajectory[0]
    reference = backbone.positions.copy()
    reference_centred = reference - backbone.center_of_mass()

    rows = []
    threaded_frames = 0
    for ts in u.trajectory[::args.stride]:
        positions = backbone.positions.copy()
        centred = positions - backbone.center_of_mass()
        _, rmsd = align.rotation_matrix(centred, reference_centred)
        rg = peptide.radius_of_gyration()

        if args.linear:
            crossings, threaded = -1, ''
        else:
            crossings = ring_crossings(ring_sel.positions.copy(),
                                       tail_sel.positions.copy())
            threaded = 'yes' if crossings % 2 == 1 else 'no'
            if crossings % 2 == 1:
                threaded_frames += 1

        bond_length = ''
        if isopeptide_pair:
            bond_length = '%.3f' % np.linalg.norm(
                isopeptide_pair[0].positions[0] - isopeptide_pair[1].positions[0])
        rows.append((ts.time, rmsd, rg, crossings, threaded, bond_length))

    # ------------------------------------------------------------------- RMSF
    # The trajectory is aligned on the backbone first, so that the fluctuation
    # measured is internal flexibility rather than overall tumbling.
    u.trajectory[0]
    align.AlignTraj(u, u, select='backbone', in_memory=True).run()
    ca = peptide.select_atoms('name CA')
    coordinates = np.array([ca.positions.copy() for _ in u.trajectory[::args.stride]])
    mean_structure = coordinates.mean(axis=0)
    rmsf = np.sqrt(((coordinates - mean_structure) ** 2).sum(axis=2).mean(axis=0))

    # ---------------------------------------------------------------- writing
    with open(os.path.join(args.outdir, 'time_series.tsv'), 'w') as fh:
        fh.write('time_ps\tRMSD_A\tRg_A\tcrossings\tthreaded\t'
                 'isopeptide_distance_A\n')
        for t, rmsd, rg, crossings, threaded, bond in rows:
            fh.write('%.1f\t%.3f\t%.3f\t%s\t%s\t%s\n'
                     % (t, rmsd, rg, '' if crossings < 0 else crossings,
                        threaded, bond))

    with open(os.path.join(args.outdir, 'rmsf_per_residue.tsv'), 'w') as fh:
        fh.write('resid\tresname\tzone\tRMSF_A\n')
        for residue, value in zip(ca.residues, rmsf):
            zone = 'ring' if r1 <= residue.resid <= racc else 'tail'
            fh.write('%d\t%s\t%s\t%.3f\n'
                     % (residue.resid, residue.resname, zone, value))

    n = len(rows)
    out = []
    W = out.append
    W('=' * 64)
    W('MD TRAJECTORY ANALYSIS - LASSO PEPTIDE')
    W('=' * 64)
    W('Topology   : %s' % args.topology)
    W('Trajectory : %s' % args.trajectory)
    W('Frames analysed: %d' % n)
    W('')
    rmsds = [r[1] for r in rows]
    rgs = [r[2] for r in rows]
    W('Backbone RMSD   mean %.2f  sd %.2f  max %.2f  A'
      % (np.mean(rmsds), np.std(rmsds), np.max(rmsds)))
    W('Radius of gyr.  mean %.2f  sd %.2f  A' % (np.mean(rgs), np.std(rgs)))
    ring_rmsf = [v for res, v in zip(ca.residues, rmsf)
                 if r1 <= res.resid <= racc]
    tail_rmsf = [v for res, v in zip(ca.residues, rmsf) if res.resid > racc]
    if ring_rmsf:
        W('Mean RMSF, ring %.2f A' % np.mean(ring_rmsf))
    if tail_rmsf:
        W('Mean RMSF, tail %.2f A' % np.mean(tail_rmsf))
    if ring_rmsf and tail_rmsf:
        W('  The tail is expected to fluctuate more than the ring. This is')
        W('  consistent with it being the region least well reproduced by')
        W('  structure predictors, and it is why conclusions that depend on')
        W('  the exact tail orientation must be qualified.')
    W('')
    if isopeptide_pair:
        lengths = [float(r[5]) for r in rows if r[5]]
        W('Isopeptide bond length: mean %.3f  max %.3f A'
          % (np.mean(lengths), np.max(lengths)))
        if np.max(lengths) > 2.0:
            W('  [!] The bond stretches beyond 2.0 A, which almost certainly')
            W('      means it is NOT defined in the simulation topology.')
            W('      Rebuild the system with the crosslink declared (CHARMM-GUI)')
            W('      before interpreting anything else in this report.')
        W('')
    if not args.linear:
        fraction = 100.0 * threaded_frames / n if n else 0
        W('TOPOLOGY')
        W('  Threaded frames: %d of %d  (%.1f %%)' % (threaded_frames, n, fraction))
        if fraction >= 99:
            W('  >>> The lasso topology is maintained throughout the')
            W('      simulation. The threading is not an artefact of the')
            W('      structure prediction.')
        elif fraction >= 80:
            W('  >>> Largely maintained, with isolated losses. Check whether')
            W('      they cluster at the start (equilibration) or increase with')
            W('      time (progressive unthreading).')
        elif fraction > 0:
            W('  >>> The topology is lost during the simulation. The steric')
            W('      plug does not retain the tail under this force field.')
        else:
            W('  >>> No threading detected in any frame. Verify that --ring is')
            W('      correct and that the starting model was threaded.')
        W('')
        first_loss = next((i for i, r in enumerate(rows) if r[4] == 'no'), None)
        if first_loss is not None:
            W('  First unthreaded frame: t = %.1f ps' % rows[first_loss][0])
    W('=' * 64)
    W('')
    W('Remember: molecular dynamics tests whether the model is STABLE under a')
    W('force field, not whether it is CORRECT. An incorrect structure can be')
    W('perfectly stable.')

    text = '\n'.join(out)
    print(text)
    with open(os.path.join(args.outdir, 'summary.txt'), 'w') as fh:
        fh.write(text + '\n')
    print('\nWritten to %s/' % args.outdir)


if __name__ == '__main__':
    main()
