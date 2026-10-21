#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aux_osmoadaptation_profile.py - Classify the osmoadaptation strategy of each
genome or MAG from its proteome.

WHY THIS SCRIPT IS PART OF THIS PROJECT
---------------------------------------
Salinibacter ruber uses the "salt-in" strategy: it accumulates molar potassium
chloride rather than synthesising compatible solutes. That has a consequence
which matters directly for the interpretation of any peptide it encodes. To
remain soluble and folded in an intracellular brine, the ENTIRE proteome must
be acidified, and that acidification is expected to extend to the products of
secondary metabolism, including RiPP core sequences.

So this script provides the physiological context in which a predicted lasso
peptide's composition is read. A strongly acidic core peptide from such an
organism is not an anomaly to be explained away; it is what the physiology
predicts. It also sets expectations for tools trained on proteomes from other
lineages, which will not have seen this compositional bias.

TWO INDEPENDENT LINES OF EVIDENCE
---------------------------------
  1. PROTEOME ACIDITY. The salt-in strategy leaves an unmistakable signature:
     the distribution of isoelectric points shifts towards acidic values and
     loses the basic peak that non-halophiles show. This is a physiological
     state legible directly from genome sequence, with no need to measure
     intracellular solutes.

  2. MARKER GENES. The alternative "salt-out" strategy depends on synthesising
     or importing compatible solutes (ectoine, glycine betaine, trehalose), on
     potassium uptake and on sodium export. Those systems are searched for by
     gene name and product description in the annotation.

Using both matters: acidity alone cannot distinguish a moderately halophilic
proteome from a salt-in one at the boundary, and marker genes alone are easily
over-called because many transporters are promiscuously annotated.

CLASSIFICATION THRESHOLDS
-------------------------
Derived from the observation that known salt-in halophiles (Halobacteria,
Salinibacter) show median proteome pI around 4.2-5.0, against 6.0-7.0 in
non-halophiles. They are stated as module constants so they can be adjusted and
the adjustment reported; they are heuristic boundaries, not measurements.

The output is a PHYSIOLOGICAL HYPOTHESIS, not a measurement. Confirming it
requires growth at controlled salinity, and the script says so when it finishes.

USAGE
-----
    python3 aux_osmoadaptation_profile.py --faa proteomes/
    python3 aux_osmoadaptation_profile.py --faa proteomes/ --annotation annot/
    python3 aux_osmoadaptation_profile.py --faa proteomes/ --figure profile.png

  --faa         folder holding one .faa per genome (Prodigal, Bakta or Prokka)
  --annotation  folder holding the annotation (.tsv from Bakta or .gff from
                Prokka), with the same base name as the .faa

Requirements: numpy. matplotlib only if --figure is requested.
Licence: MIT
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

# pKa values, EMBOSS set. Stated explicitly because different sets shift the
# computed pI by a few tenths, and the choice must be reproducible.
PKA_POSITIVE = {'K': 10.8, 'R': 12.5, 'H': 6.5}
PKA_NEGATIVE = {'D': 3.9, 'E': 4.1, 'C': 8.5, 'Y': 10.1}
PKA_NTERM, PKA_CTERM = 8.6, 3.6

# Classification thresholds; see the module docstring for their origin.
PI_SALT_IN = 5.2
PI_INTERMEDIATE = 6.0
ACIDIC_FRACTION_SALT_IN = 0.70

# Marker gene families, matched by gene name or product description. Patterns
# are deliberately loose because annotation vocabulary varies between tools;
# the counts are indicative, which is why they never decide the classification
# on their own.
MARKERS = {
    'Ectoine (synthesis)': [r'\bect[ABCD]\b', r'ectoine synthase', r'diaminobutyr',
                            r'ectoine hydroxylase'],
    'Glycine betaine (synthesis)': [r'\bbet[ABIT]\b', r'choline dehydrogenase',
                                    r'betaine-aldehyde dehydrogenase',
                                    r'glycine betaine', r'\bgbsA\b'],
    'Glycine betaine (transport)': [r'\bopu[ABCD]', r'\bpro[UVWX]\b', r'\bproVWX\b',
                                    r'betaine.*transport',
                                    r'glycine betaine/proline', r'\bbetT\b',
                                    r'\bbetP\b'],
    'Trehalose': [r'\bots[AB]\b', r'\btre[SYZ]\b', r'trehalose-6-phosphate',
                  r'trehalose synthase', r'\bglgE\b'],
    'Glucosylglycerate': [r'glucosylglycerate', r'\bggpS\b', r'glucosylglycerol'],
    'Potassium uptake': [r'\btrk[AGH]\b', r'\bkdp[ABCDEF]\b', r'\bkup\b',
                         r'\bktrAB?\b', r'potassium.*transport',
                         r'potassium uptake'],
    'Na+/H+ antiport': [r'\bnha[ABCDEKP]\b', r'\bmrp[ABCDEFG]\b',
                        r'\bmnh[ABCDEFG]\b', r'sodium/(proton|hydrogen) antiport',
                        r'Na\(\+\)/H\(\+\) antiport'],
    'Mechanosensitive channels': [r'\bmsc[LSK]\b', r'mechanosensitive'],
    'Rhodopsins': [r'rhodopsin', r'\bbop\b', r'\bhop\b', r'bacteriorhodopsin',
                   r'halorhodopsin', r'\bblh\b', r'\bbrp\b'],
    'Bacterioruberin': [r'\bcrt[BIDY]\b', r'\blyeJ\b', r'\bcruF\b', r'phytoene',
                        r'lycopene elongase'],
}


def read_fasta(path):
    """Yield (header, sequence) pairs."""
    header, sequence = None, []
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence)
                header, sequence = line[1:], []
            elif line:
                sequence.append(line)
    if header is not None:
        yield header, ''.join(sequence)


def charge_function(sequence):
    """Return a function giving net charge as a function of pH.

    Residue counts are taken once and closed over, rather than recounted at
    every pH, because the bisection evaluates the function ~100 times per
    protein and a proteome holds thousands of proteins.
    """
    positive_counts = {a: sequence.count(a) for a in PKA_POSITIVE}
    negative_counts = {a: sequence.count(a) for a in PKA_NEGATIVE}

    def q(ph):
        total = 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))          # amino terminus
        for a, pka in PKA_POSITIVE.items():
            if positive_counts[a]:
                total += positive_counts[a] / (1.0 + 10 ** (ph - pka))
        total -= 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))         # carboxy terminus
        for a, pka in PKA_NEGATIVE.items():
            if negative_counts[a]:
                total -= negative_counts[a] / (1.0 + 10 ** (pka - ph))
        return total
    return q


def isoelectric_point(sequence, lo=0.0, hi=14.0, tol=1e-4):
    """pI by bisection on the charge curve, clamped to the 0-14 interval."""
    q = charge_function(sequence)
    if q(lo) < 0:
        return lo
    if q(hi) > 0:
        return hi
    for _ in range(100):
        mid = (lo + hi) / 2.0
        value = q(mid)
        if abs(value) < tol or (hi - lo) < tol:
            return mid
        if value > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def analyse_proteome(path, min_length=40):
    """Isoelectric points of every protein above the length cutoff.

    Short predicted ORFs are excluded because their pI is dominated by the
    terminal charges and would blur the distribution the method relies on.
    Non-standard characters are stripped so that X, B, Z or stop symbols do not
    silently enter the counts.
    """
    points = []
    for _, sequence in read_fasta(path):
        sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
        if len(sequence) < min_length:
            continue
        points.append(isoelectric_point(sequence))
    return np.array(points)


def classify(median_pi, acidic_fraction):
    if median_pi < PI_SALT_IN and acidic_fraction >= ACIDIC_FRACTION_SALT_IN:
        return 'salt-in (acidified proteome)'
    if median_pi < PI_INTERMEDIATE:
        return 'intermediate / moderate halophile'
    return 'salt-out (non-acidified proteome)'


def count_markers(annotation_path):
    """Count matches per marker category in a .tsv or .gff annotation file."""
    try:
        text = open(annotation_path, encoding='utf-8', errors='replace').read()
    except OSError:
        return None
    out = {}
    for category, patterns in MARKERS.items():
        n = 0
        for pattern in patterns:
            n += len(re.findall(pattern, text, flags=re.IGNORECASE))
        out[category] = n
    return out


def main():
    ap = argparse.ArgumentParser(
        description='Osmoadaptation profile of hypersaline genomes and MAGs')
    ap.add_argument('--faa', required=True, help='folder holding the .faa proteomes')
    ap.add_argument('--annotation', default=None,
                    help='folder holding the annotation (.tsv or .gff)')
    ap.add_argument('--out', default='osmoadaptation_profile.tsv')
    ap.add_argument('--figure', default=None, help='PNG of the pI distributions')
    ap.add_argument('--min-length', type=int, default=40,
                    help='minimum protein length in residues (40)')
    args = ap.parse_args()

    proteomes = sorted(glob.glob(os.path.join(args.faa, '*.faa')) +
                       glob.glob(os.path.join(args.faa, '*.fasta')) +
                       glob.glob(os.path.join(args.faa, '*.fa')))
    if not proteomes:
        sys.exit('No proteomes found in %s' % args.faa)

    rows, distributions = [], {}
    categories = list(MARKERS)
    print('%d proteomes\n' % len(proteomes))
    print('%-28s%8s%9s%10s  strategy' % ('genome', 'n prot', 'median pI', '% acidic'))
    print('-' * 90)

    for path in proteomes:
        name = os.path.splitext(os.path.basename(path))[0]
        points = analyse_proteome(path, args.min_length)
        if len(points) == 0:
            print('%-28s  no usable proteins' % name)
            continue
        median_pi = float(np.median(points))
        # Fraction below pI 5.0: a stricter cut than the classification
        # threshold, chosen to capture the bulk of the acidic peak rather than
        # the boundary region.
        acidic_fraction = float((points < 5.0).mean())
        strategy = classify(median_pi, acidic_fraction)
        distributions[name] = points
        row = {'genome': name, 'n_proteins': len(points),
               'median_pI': round(median_pi, 2),
               'mean_pI': round(float(points.mean()), 2),
               'fraction_pI_below_5': round(acidic_fraction, 3),
               'strategy': strategy}

        if args.annotation:
            annotation_file = None
            for ext in ('.tsv', '.gff3', '.gff', '.txt'):
                candidate = os.path.join(args.annotation, name + ext)
                if os.path.isfile(candidate):
                    annotation_file = candidate
                    break
            markers = count_markers(annotation_file) if annotation_file else None
            for c in categories:
                row[c] = markers[c] if markers else 'NA'
        rows.append(row)
        print('%-28s%8d%9.2f%9.1f%%  %s'
              % (name, len(points), median_pi, acidic_fraction * 100, strategy))

    if not rows:
        sys.exit('\nNo proteome could be analysed.')

    columns = ['genome', 'n_proteins', 'median_pI', 'mean_pI',
               'fraction_pI_below_5', 'strategy']
    if args.annotation:
        columns += categories
    with open(args.out, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(columns) + '\n')
        for row in rows:
            fh.write('\t'.join(str(row.get(c, 'NA')) for c in columns) + '\n')
    print('\nTable -> %s' % args.out)

    if args.figure:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        n = len(distributions)
        n_cols = min(3, n)
        n_rows = (n + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(4.6 * n_cols, 3.1 * n_rows),
                                 squeeze=False, facecolor='#FCFCFB')
        for ax, (name, points) in zip(axes.ravel(), distributions.items()):
            median_pi = np.median(points)
            # Colour encodes the classification, so the panel is readable
            # without consulting the table.
            colour = ('#2a78d6' if median_pi < PI_SALT_IN
                      else ('#eda100' if median_pi < PI_INTERMEDIATE else '#eb6834'))
            ax.hist(points, bins=np.arange(3, 13.2, 0.2), color=colour,
                    edgecolor='none')
            ax.axvline(median_pi, color='#0b0b0b', lw=1.4, ls='--')
            ax.set_title('%s\nmedian pI = %.2f' % (name, median_pi), fontsize=9.5)
            ax.set_xlabel('isoelectric point', fontsize=8.5)
            ax.set_facecolor('#FCFCFB')
            for side in ('top', 'right'):
                ax.spines[side].set_visible(False)
            ax.tick_params(labelsize=8)
        for ax in axes.ravel()[n:]:
            ax.axis('off')
        fig.tight_layout()
        fig.savefig(args.figure, dpi=190, facecolor='#FCFCFB')
        print('Figure -> %s' % args.figure)

    print('\nReading the result: a proteome with median pI below 5.2 and more')
    print('than 70% acidic proteins indicates salt-in. Non-halophiles show a')
    print('bimodal distribution with medians between 6 and 7. This is a')
    print('physiological HYPOTHESIS, not a measurement; confirm it with growth')
    print('at controlled salinity.')


if __name__ == '__main__':
    main()
