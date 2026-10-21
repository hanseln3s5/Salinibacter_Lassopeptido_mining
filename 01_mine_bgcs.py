#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_mine_bgcs.py - Targeted, antiSMASH-independent search for lasso peptide
biosynthetic gene clusters across a set of annotated genomes.

WHY THIS STEP EXISTS
--------------------
It answers the question a reviewer will ask about any prevalence figure: "are
the negative genomes really negative, or did antiSMASH miss the cluster?"

Reporting that a cluster occurs in N of M genomes is a claim about biology only
if the M-N negatives are true absences. Two technical artefacts produce false
negatives, and both are testable:

  - Fragmented assemblies. A cluster split across two contigs may not satisfy
    antiSMASH's rule-based detection.
  - Contig edges. A cluster truncated by the end of a contig loses the genes
    the rule needs.

So the script does three things, independently of antiSMASH.

  1. ASSEMBLY QUALITY. Contig count, N50, total length and GC content per
     genome. If the negatives are systematically more fragmented than the
     positives, the observed prevalence is a LOWER BOUND, not a biological
     result, and must be reported as such.

  2. HOMOLOGY SEARCH. hmmsearch of the marker domains against the proteomes,
     bypassing antiSMASH's rules entirely:

       PF13471 (Transglut_core3)  lasso peptide B protein: RiPP recognition
                                  element plus leader peptidase.
                                  THIS IS THE DISCRIMINATING MARKER.

       PF00733 (Asn_synthase)     the C cyclase, homologous to asparagine
                                  synthetase B.
                                  NOT specific on its own: housekeeping
                                  asparagine synthetase is present in virtually
                                  every genome, so counting PF00733 hits as
                                  evidence of a lasso peptide cluster would
                                  make almost every genome a false positive.
                                  It counts ONLY when CO-LOCALISED with PF13471.

     Hence the primary criterion is PF13471, and PF00733 serves to confirm that
     the operon's cyclase is present alongside it.

  3. CONTIG EDGES. For each hit, the distance to the nearest end of its contig.
     A cluster sitting against an edge may be truncated, which is the second
     route by which a real BGC is missed.

The result is cross-tabulated against the antiSMASH positives and probable
false negatives are flagged.

USAGE
-----
    # assembly-quality report only (does not need HMMER)
    python3 01_mine_bgcs.py --prokka prokka_out/ --no-hmm

    # full search
    python3 01_mine_bgcs.py --prokka prokka_out/ --hmm hmms/ \
            --antismash antismash_positives.txt --outdir results/

Expected layout under --prokka: one subdirectory per genome, each containing
the .faa and .gff files Prokka produces. The subdirectory name is used as the
genome identifier.

Fetch the HMMs from InterPro/Pfam with:
    python3 01_mine_bgcs.py --download-hmm hmms/

Requirements: HMMER (hmmsearch) on PATH for the search. Nothing else; Python 3
standard library only.

Licence: MIT
"""

import argparse
import gzip
import os
import subprocess
import sys
import urllib.request

DOMAINS = {
    'PF13471': 'B protein (RRE + leader peptidase) - DISCRIMINATING MARKER',
    'PF00733': 'C cyclase (Asn_synthase) - valid only when co-localised with PF13471',
}
PFAM_URL = 'https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{}?annotation=hmm'


# ------------------------------------------------------------------- utilities
def download_hmms(destination):
    """Fetch the required Pfam HMMs."""
    os.makedirs(destination, exist_ok=True)
    for pf in DOMAINS:
        path = os.path.join(destination, '%s.hmm' % pf)
        if os.path.isfile(path):
            print('  %s: already present, skipped' % pf)
            continue
        print('  downloading %s ...' % pf, end=' ', flush=True)
        try:
            with urllib.request.urlopen(PFAM_URL.format(pf), timeout=60) as r:
                data = r.read()
            if data[:2] == b'\x1f\x8b':          # served gzipped
                data = gzip.decompress(data)
            with open(path, 'wb') as fh:
                fh.write(data)
            print('ok')
        except Exception as exc:
            print('FAILED (%s)' % exc)
            print('    Download it manually from: %s' % PFAM_URL.format(pf))
    return 0


def n50(lengths):
    """N50 of a list of contig lengths: the length at which half the assembly
    is contained in contigs of that size or longer."""
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2.0
    cumulative = 0
    for L in ordered:
        cumulative += L
        if cumulative >= half:
            return L
    return ordered[-1]


def fasta_stats(path):
    """(n_contigs, total_length, N50, GC%) of a FASTA, without dependencies."""
    lengths, gc_count, base_count, current = [], 0, 0, 0
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('>'):
                if current:
                    lengths.append(current)
                current = 0
            else:
                s = line.strip()
                current += len(s)
                upper = s.upper()
                gc_count += upper.count('G') + upper.count('C')
                base_count += len(s)
    if current:
        lengths.append(current)
    gc = (100.0 * gc_count / base_count) if base_count else 0.0
    return len(lengths), sum(lengths), n50(lengths), gc


def read_gff(path):
    """Return (coordinates, contig_lengths).

    coordinates    : {protein_id: (contig, start, end)}
    contig_lengths : {contig: length}

    Parsing stops at Prokka's trailing ##FASTA block, and contig lengths are
    taken from the ##sequence-region directives; without them the distance to
    the contig edge cannot be computed.
    """
    coordinates, contig_lengths = {}, {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                break
            if line.startswith('##sequence-region'):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        contig_lengths[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
                continue
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            protein_id = None
            for attribute in fields[8].split(';'):
                if attribute.startswith('ID='):
                    protein_id = attribute[3:].strip()
                    break
            if not protein_id:
                continue
            try:
                coordinates[protein_id] = (fields[0], int(fields[3]), int(fields[4]))
            except ValueError:
                continue
    return coordinates, contig_lengths


def run_hmmsearch(hmm, faa, evalue, tmp):
    """Return {protein_id: best_evalue} for one HMM against one proteome."""
    cmd = ['hmmsearch', '--tblout', tmp, '-E', str(evalue), '--cpu', '1', hmm, faa]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE)
    except FileNotFoundError:
        sys.exit('hmmsearch not found. Install it:  '
                 'conda install -c bioconda hmmer')
    except subprocess.CalledProcessError as exc:
        print('    WARNING: hmmsearch failed on %s: %s'
              % (faa, exc.stderr.decode()[:200]))
        return {}
    hits = {}
    with open(tmp, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            columns = line.split()
            if len(columns) < 5:
                continue
            protein_id, ev = columns[0], float(columns[4])
            if protein_id not in hits or ev < hits[protein_id]:
                hits[protein_id] = ev
    return hits


def locate_files(directory):
    """Find the .faa / .gff / .fna inside one genome directory."""
    faa = gff = fna = None
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if not os.path.isfile(path):
            continue
        if name.endswith('.faa') and faa is None:
            faa = path
        elif name.endswith('.gff') and gff is None:
            gff = path
        elif name.endswith(('.fna', '.fa', '.fasta')) and fna is None:
            fna = path
    return faa, gff, fna


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description='Targeted, antiSMASH-independent search for lasso peptide BGCs')
    ap.add_argument('--prokka',
                    help='directory with one subdirectory per genome (Prokka output)')
    ap.add_argument('--hmm', help='directory holding PF13471.hmm and PF00733.hmm')
    ap.add_argument('--download-hmm', metavar='DIR',
                    help='download the Pfam HMMs and exit')
    ap.add_argument('--antismash',
                    help='file with one genome ID per line: the antiSMASH positives')
    ap.add_argument('--outdir', default='lasso_results', help='output directory')
    ap.add_argument('--evalue', type=float, default=1e-5,
                    help='E-value threshold (1e-5)')
    ap.add_argument('--window', type=int, default=10000,
                    help='maximum distance in bp between B and C for them to '
                         'count as co-localised (10000). Lasso operons are '
                         'compact, so this is generous')
    ap.add_argument('--edge', type=int, default=5000,
                    help='distance to the contig end below which truncation is '
                         'flagged (5000)')
    ap.add_argument('--no-hmm', action='store_true',
                    help='assembly-quality report only')
    args = ap.parse_args()

    if args.download_hmm:
        print('Downloading Pfam HMMs:')
        return download_hmms(args.download_hmm)

    if not args.prokka:
        sys.exit('--prokka is required (or use --download-hmm).')
    if not os.path.isdir(args.prokka):
        sys.exit('Directory not found: %s' % args.prokka)
    if not args.no_hmm and not args.hmm:
        sys.exit('--hmm is required. Fetch them with --download-hmm hmms/ '
                 '(or use --no-hmm)')

    genomes = sorted(d for d in os.listdir(args.prokka)
                     if os.path.isdir(os.path.join(args.prokka, d)))
    if not genomes:
        sys.exit('No genome subdirectories found in %s' % args.prokka)
    os.makedirs(args.outdir, exist_ok=True)

    antismash_positives = set()
    if args.antismash:
        with open(args.antismash, encoding='utf-8') as fh:
            antismash_positives = {l.strip() for l in fh
                                   if l.strip() and not l.startswith('#')}
        print('antiSMASH positives declared: %d' % len(antismash_positives))

    print('Genomes found: %d\n' % len(genomes))

    quality_rows, hit_rows, summary_rows = [], [], []
    tmp = os.path.join(args.outdir, '_hmm.tmp')

    for i, genome in enumerate(genomes, 1):
        directory = os.path.join(args.prokka, genome)
        faa, gff, fna = locate_files(directory)
        print('[%3d/%d] %s' % (i, len(genomes), genome), end='')

        # --- assembly quality
        n_contigs = total_length = assembly_n50 = 0
        gc = 0.0
        if fna:
            n_contigs, total_length, assembly_n50, gc = fasta_stats(fna)
        quality_rows.append((genome, n_contigs, total_length, assembly_n50,
                             '%.2f' % gc))

        if args.no_hmm:
            print('   contigs=%d N50=%d' % (n_contigs, assembly_n50))
            continue

        if not faa or not gff:
            print('   NO .faa/.gff - search skipped')
            summary_rows.append((genome, 'NA', 'NA', 'NA', 'NA',
                                 'yes' if genome in antismash_positives else 'no',
                                 'no_annotation'))
            continue

        coordinates, contig_lengths = read_gff(gff)

        # --- hmmsearch
        hits = {}
        for pf in DOMAINS:
            hmm_path = os.path.join(args.hmm, '%s.hmm' % pf)
            if not os.path.isfile(hmm_path):
                sys.exit('%s not found. Use --download-hmm.' % hmm_path)
            hits[pf] = run_hmmsearch(hmm_path, faa, args.evalue, tmp)

        # --- contig-edge distance and hit records
        def edge_distance(protein_id):
            if protein_id not in coordinates:
                return None
            contig, start, end = coordinates[protein_id]
            length = contig_lengths.get(contig)
            if not length:
                return None
            return min(start, length - end)

        for pf, pf_hits in hits.items():
            for protein_id, ev in sorted(pf_hits.items(), key=lambda x: x[1]):
                contig, start, end = coordinates.get(protein_id, ('NA', 0, 0))
                edge = edge_distance(protein_id)
                hit_rows.append((genome, pf, protein_id, contig, start, end,
                                 '%.2g' % ev, edge if edge is not None else 'NA'))

        # --- criterion: PF13471 discriminates; PF00733 only if co-localised
        b_hits = hits['PF13471']
        c_hits = hits['PF00733']
        has_b = len(b_hits) > 0
        co_localised = False
        edge_dist = 'NA'
        if has_b:
            edges = []
            for b_id in b_hits:
                if b_id not in coordinates:
                    continue
                b_contig, b_start, b_end = coordinates[b_id]
                b_edge = edge_distance(b_id)
                if b_edge is not None:
                    edges.append(b_edge)
                for c_id in c_hits:
                    if c_id not in coordinates:
                        continue
                    c_contig, c_start, c_end = coordinates[c_id]
                    # Same contig and within the window, measured between the
                    # facing ends of the two genes.
                    if (c_contig == b_contig
                            and min(abs(c_start - b_end),
                                    abs(b_start - c_end)) <= args.window):
                        co_localised = True
            if edges:
                edge_dist = min(edges)

        in_antismash = genome in antismash_positives
        if has_b and co_localised:
            verdict = 'confirmed' if in_antismash else 'PROBABLE_FALSE_NEGATIVE'
        elif has_b:
            verdict = ('partially_confirmed' if in_antismash
                       else 'B_without_colocalised_cyclase')
        else:
            verdict = 'absent'
        if has_b and edge_dist != 'NA' and edge_dist < args.edge:
            verdict += '_CONTIG_EDGE'

        summary_rows.append((genome, len(b_hits), len(c_hits),
                             'yes' if co_localised else 'no', edge_dist,
                             'yes' if in_antismash else 'no', verdict))
        print('   PF13471=%d  PF00733=%d  -> %s'
              % (len(b_hits), len(c_hits), verdict))

    if os.path.isfile(tmp):
        os.remove(tmp)

    # --------------------------------------------------------------- writing
    def write_tsv(name, header, rows):
        path = os.path.join(args.outdir, name)
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write('\t'.join(header) + '\n')
            for row in rows:
                fh.write('\t'.join(str(x) for x in row) + '\n')
        return path

    p1 = write_tsv('assembly_quality.tsv',
                   ['genome', 'n_contigs', 'total_length', 'N50', 'GC_pct'],
                   quality_rows)
    print('\n  assembly quality   -> %s' % p1)

    if not args.no_hmm:
        p2 = write_tsv('hmm_hits.tsv',
                       ['genome', 'domain', 'protein', 'contig', 'start', 'end',
                        'evalue', 'contig_edge_distance'], hit_rows)
        p3 = write_tsv('lasso_summary.tsv',
                       ['genome', 'n_PF13471', 'n_PF00733', 'B_and_C_colocalised',
                        'contig_edge_distance', 'antismash_positive', 'verdict'],
                       summary_rows)
        print('  domain hits        -> %s' % p2)
        print('  per-genome summary -> %s' % p3)

        n_with_b = sum(1 for r in summary_rows if r[1] not in ('NA', 0))
        n_colocalised = sum(1 for r in summary_rows if r[3] == 'yes')
        n_false_neg = sum(1 for r in summary_rows
                          if 'FALSE_NEGATIVE' in str(r[6]))
        n_edge = sum(1 for r in summary_rows if 'CONTIG_EDGE' in str(r[6]))

        print('\n' + '=' * 68)
        print('  SUMMARY')
        print('=' * 68)
        print('  Genomes analysed                        : %d' % len(genomes))
        print('  With PF13471 (discriminating marker)    : %d' % n_with_b)
        print('  With B and cyclase co-localised         : %d' % n_colocalised)
        if antismash_positives:
            print('  Declared positive by antiSMASH          : %d'
                  % len(antismash_positives))
            print('  PROBABLE FALSE NEGATIVES                : %d' % n_false_neg)
        print('  With the marker next to a contig edge   : %d' % n_edge)
        print('=' * 68)
        if n_false_neg:
            print('\n  Some genomes carry the operon machinery that antiSMASH did')
            print('  not flag. Inspect those regions manually, looking for the')
            print('  precursor gene (40-60 codons), which Prodigal frequently')
            print('  misses because it falls below its minimum ORF length.')
        if n_edge:
            print('\n  At least one cluster sits against a contig edge and may be')
            print('  truncated. Declare this in the limitations: the prevalence')
            print('  you report would be a lower bound.')
        print('\n  Now compare n_contigs and N50 between positives and negatives')
        print('  in assembly_quality.tsv. If the negatives are more fragmented,')
        print('  the difference may be technical rather than biological.')

    return 0


if __name__ == '__main__':
    sys.exit(main())
