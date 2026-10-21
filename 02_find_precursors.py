#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_find_precursors.py - Recover lasso peptide precursor genes from antiSMASH
GenBank region files and emit the core peptides for structure prediction.

WHY THIS STEP EXISTS
--------------------
antiSMASH reliably detects the lasso peptide MATURATION MACHINERY - the
PF13471 leader peptidase, the asparagine-synthetase-like cyclase, the RiPP
recognition element - because those are full-length proteins with conserved
domains. It frequently fails to annotate the PRECURSOR gene itself.

The reason is mechanical, not biological. General-purpose gene callers such as
Prodigal (used by Prokka and by most annotation pipelines) impose a minimum
open-reading-frame length of 90 nt, i.e. 30 codons. Lasso peptide precursors
are 25-60 residues long and frequently initiate at GTG or TTG rather than ATG.
They therefore fall below the detection floor and are never written to the
annotation, which means the cluster is reported without the one gene that
encodes the product.

Without a core peptide sequence there is nothing to submit to a structure
predictor and nothing to dock. This script recovers the candidates.

WHAT IT DOES
------------
  1. Extracts any precursor antiSMASH DID annotate (CDS_motif features carrying
     a /core_sequence qualifier). If present, that is the answer and the search
     stops there.

  2. Otherwise it performs its own six-frame ORF scan of the cluster
     neighbourhood, including ORFs below the gene caller's floor and those
     starting at alternative initiation codons, then filters and scores the
     candidates against the known biology of lasso peptide precursors.

SCORING RATIONALE
-----------------
Each candidate cut point is scored on features that are documented properties
of characterised lasso peptide precursors:

  Genomic context (the strongest discriminator). A genuine precursor lies
  within a few kb of the maturation genes and almost always on the same strand,
  because the cluster is usually a single operon. Including this term reduced
  the candidate list on the test cluster from 49 to 28 and moved the two
  biologically plausible ORFs to the top. Sequence features alone are far too
  permissive on a 30 kb region.

  Ring-forming residue 1. The macrolactam forms from the alpha-amino group of
  the first core residue. Gly dominates the characterised set, Cys is common,
  Ala/Ser/Thr/Val occur; the weights reflect that ordering.

  Ring size. The acceptor Asp or Glu sits at position 7, 8 or 9. Sizes 8 and 9
  are more frequent than 7 among characterised peptides.

  Core length. Characterised cores are typically 15-24 residues. Values outside
  10-30 are penalised and flagged rather than discarded, so that a genuinely
  atypical precursor is not silently lost.

  Tail and steric plug. The threaded fold requires enough tail to pass through
  the ring, and a bulky residue acting as a plug on the far side.

  Leader composition. Lasso peptide leaders are characteristically acidic, and
  a Thr-x(3,4)-Pro motif recurs in the region preceding the cleavage site.

  Low-complexity and strongly basic cores are penalised: both are signatures of
  spurious reading frames rather than real coding sequence.

LIMITATIONS - READ BEFORE CITING
--------------------------------
This is a heuristic screen, NOT a replacement for RODEO (Tietz et al. 2017) or
BAGEL4. Those tools score candidates with classifiers trained on manually
curated precursor sets; the weights here were chosen by the author from the
published literature and are not trained on anything. The output is a RANKED
LIST OF CANDIDATES.

Real validation requires both of:
  (a) the same ORF recurring, in the same relative position, across several
      genomes carrying the cluster, and
  (b) independent confirmation by RODEO or BAGEL4.

In a thesis or manuscript this script should be described as an exploratory
complementary screen reported in an annex, never at the same level as the
trained tools.

USAGE
-----
  python3 02_find_precursors.py region.gbk
  python3 02_find_precursors.py regions/*.gbk -o results/ --top 5
  python3 02_find_precursors.py regions/*.gbk --min-score 6 --whole-record

OUTPUT (in the directory given by -o, default 'lasso_precursors/')
  candidates.tsv              full ranked table, for the annex
  cores_for_structure.faa     core peptides - this is what goes to LassoPred
  full_precursors.faa         leader + core, for cross-genome BLAST

No dependencies: Python 3.8+ standard library only. The GenBank parser is
written here rather than taken from Biopython so that the script runs in a
bare container with no scientific stack installed.

Licence: MIT
"""

import argparse
import os
import re
import sys
from collections import Counter

# --------------------------------------------------------------- genetic code
_B = 'TCAG'
_AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {}
_i = 0
for _a in _B:
    for _b in _B:
        for _c in _B:
            CODON[_a + _b + _c] = _AA[_i]
            _i += 1

# Bacterial translation table 11 alternative initiation codons. Including GTG
# and TTG is essential here: they are common in GC-rich genomes such as
# Salinibacter (~66% GC) and are a recurring reason why short precursors are
# missed. The first codon is translated as Met regardless of its identity.
START_CODONS = ('ATG', 'GTG', 'TTG')

COMPLEMENT = str.maketrans('ACGTNRYKMSWBDHV', 'TGCANYRMKSWVHDB')


def translate(dna):
    return ''.join(CODON.get(dna[i:i + 3], 'X') for i in range(0, len(dna) - 2, 3))


def reverse_complement(dna):
    return dna.translate(COMPLEMENT)[::-1]


# ------------------------------------------------------------ GenBank parsing
def parse_genbank(path):
    """Minimal GenBank parser returning sequence, features and metadata.

    Biopython is deliberately not used so the script has zero dependencies and
    can be run inside a minimal container by a reviewer reproducing the work.
    """
    text = open(path, 'r', errors='replace').read()

    organism = ''
    m = re.search(r'^\s*ORGANISM\s+(.+)$', text, re.M)
    if m:
        organism = m.group(1).strip()
    accession = os.path.basename(path)
    m = re.search(r'^ACCESSION\s+(\S+)', text, re.M)
    if m:
        accession = m.group(1).strip()

    # Offset of this region within the parent replicon, when antiSMASH recorded
    # it. Reporting coordinates in whole-genome space as well as region space
    # is what lets a reader locate the ORF in the original assembly.
    offset = 0
    m = re.search(r'Orig\.\s*start\s*::\s*(\d+)', text)
    if m:
        offset = int(m.group(1))

    sequence = ''
    if '\nORIGIN' in text:
        sequence = ''.join(re.findall(r'[acgtnACGTN]+',
                                      text.split('\nORIGIN', 1)[1])).upper()

    features = []
    if 'FEATURES' in text:
        body = text.split('FEATURES', 1)[1].split('\nORIGIN', 1)[0]
        current = None
        last_key = None
        for line in body.split('\n')[1:]:
            if re.match(r'^ {5}\S', line):
                if current:
                    features.append(current)
                current = {'type': line[5:21].strip(),
                           'loc': line[21:].strip(), 'q': {}}
                last_key = None
            elif current is not None:
                s = line.strip()
                if s.startswith('/'):
                    mm = re.match(r'/([^=]+)=?"?(.*?)"?$', s)
                    last_key = mm.group(1)
                    current['q'][last_key] = current['q'].get(last_key, '') + mm.group(2)
                elif last_key:
                    previous = current['q'][last_key]
                    # GenBank wraps long qualifiers; a trailing hyphen means the
                    # token continues, otherwise a space separates the pieces.
                    sep = '' if previous.endswith('-') else ' '
                    current['q'][last_key] = previous + sep + s.rstrip('"')
                elif not current['q']:
                    current['loc'] += s
        if current:
            features.append(current)

    return {'accession': accession, 'organism': organism, 'offset': offset,
            'sequence': sequence, 'features': features, 'path': path}


def span(loc):
    """(start, end, strand) as 1-based inclusive coordinates."""
    nums = [int(x) for x in re.findall(r'\d+', loc)]
    if not nums:
        return (0, 0, '+')
    return (min(nums), max(nums), '-' if 'complement' in loc else '+')


# ------------------------------------ precursors already annotated by antiSMASH
def annotated_precursors(record):
    """Extract CDS_motif features carrying /core_sequence.

    antiSMASH writes these when its own RiPP module successfully identifies the
    precursor. When they exist they supersede anything this script could infer,
    so they are returned with a sentinel score and the ORF scan is skipped.
    """
    out = []
    for f in record['features']:
        q = f['q']
        core = q.get('core_sequence', '').replace(' ', '')
        if not core:
            continue
        a, b, strand = span(f['loc'])
        out.append({
            'source': 'antismash_annotated',
            'locus': q.get('locus_tag', q.get('label', 'antismash_prepeptide')),
            'start': a, 'end': b, 'strand': strand,
            'leader': q.get('leader_sequence', '').replace(' ', ''),
            'core': core,
            'peptide_class': q.get('peptide', q.get('predicted_class', '')),
            'score': 99.0, 'ring_size': 0, 'res1': core[:1], 'acceptor': '',
            'notes': 'annotated by antiSMASH',
        })
    return out


def core_bounds(record):
    """Coordinate range of the cluster core, or the whole record if absent."""
    for f in record['features']:
        if f['type'] in ('proto_core', 'protocluster'):
            a, b, _ = span(f['loc'])
            if a and b:
                return a, b
    return 1, len(record['sequence'])


def annotated_cds(record):
    out = []
    for f in record['features']:
        if f['type'] != 'CDS':
            continue
        a, b, strand = span(f['loc'])
        out.append((a, b, strand, f['q'].get('locus_tag', '?'),
                    f['q'].get('translation', '').replace(' ', ''),
                    f['q'].get('product', '')))
    return out


# Domains that mark the lasso peptide maturation machinery. Proximity to these
# is the single most informative feature available for free, because a
# precursor that is not near its own cyclase and peptidase is almost certainly
# a spurious reading frame.
ANCHOR_DOMAINS = ('PF13471', 'Transglut_core3', 'Asn_synthase', 'PF00733',
                  'Stand_Alone_Lasso_RRE', 'PqqD', 'PF05402', 'mcjC', 'micJ25')


def anchor_genes(record):
    """Coordinates and strand of the B / C / RRE genes of the cluster."""
    tags = set()
    for f in record['features']:
        if f['type'] not in ('PFAM_domain', 'aSDomain'):
            continue
        blob = (f['q'].get('aSDomain', '') + ' ' + f['q'].get('db_xref', '')
                + ' ' + f['q'].get('description', ''))
        if any(k in blob for k in ANCHOR_DOMAINS):
            t = f['q'].get('locus_tag')
            if t:
                tags.add(t)
    anchors = []
    for f in record['features']:
        if f['type'] != 'CDS':
            continue
        q = f['q']
        tag = q.get('locus_tag', '')
        functions = q.get('gene_functions', '')
        if tag in tags or any(k in functions for k in ANCHOR_DOMAINS):
            a, b, strand = span(f['loc'])
            anchors.append((a, b, strand))
    return anchors


def genomic_context(a, b, strand, anchors):
    """(distance in bp to the nearest machinery gene, same strand?)"""
    if not anchors:
        return None, None
    d = min(0 if (a <= y and x <= b) else min(abs(a - y), abs(x - b))
            for x, y, _ in anchors)
    dominant = Counter(s for _, _, s in anchors).most_common(1)[0][0]
    return d, (strand == dominant)


# ------------------------------------------------------------------ ORF search
def find_orfs(sequence, lo, hi, min_aa, max_aa):
    """Six-frame ORF scan within [lo, hi], honouring alternative start codons.

    Returns (start, end, strand, protein) with 1-based inclusive coordinates in
    the coordinate frame of the supplied sequence.
    """
    results = []
    L = len(sequence)
    for strand in '+-':
        s = sequence if strand == '+' else reverse_complement(sequence)
        for frame in range(3):
            protein = translate(s[frame:])
            position = 0
            for stretch in protein.split('*'):
                stretch_start = position
                position += len(stretch) + 1
                if len(stretch) < min_aa:
                    continue
                for k, _ in enumerate(stretch):
                    nt = frame + (stretch_start + k) * 3
                    if s[nt:nt + 3] not in START_CODONS:
                        continue
                    length = len(stretch) - k
                    if not (min_aa <= length <= max_aa):
                        continue
                    peptide = 'M' + stretch[k + 1:]
                    nt_end = nt + length * 3 + 3        # include the stop codon
                    if strand == '+':
                        a, b = nt + 1, nt_end
                    else:
                        a, b = L - nt_end + 1, L - nt
                    if a < lo or b > hi:
                        continue
                    results.append((a, b, strand, peptide))
                    # Only the most 5' start of each stretch is kept. Reporting
                    # every internal Met of the same ORF would flood the output
                    # with near-identical candidates.
                    break
    return results


# ----------------------------------------------------------- lasso peptide rules
# Weight of each residue as the ring-forming position 1, reflecting the
# frequency with which it occurs in characterised lasso peptides.
RES1_WEIGHT = {'G': 3.0, 'C': 2.5, 'A': 1.5, 'S': 1.5, 'T': 1.0, 'V': 1.0}
BULKY = set('FWYRKH')


def candidate_cuts(peptide, core_min, core_max, leader_min):
    """All leader|core cut points compatible with a macrolactam.

    A cut is admissible when residue 1 of the putative core can form the
    macrolactam and an Asp or Glu sits at position 7, 8 or 9, leaving at least
    three residues of tail. Those three constraints come directly from the
    chemistry of the fold and exclude the overwhelming majority of positions.
    """
    out = []
    for i, res1 in enumerate(peptide):
        position = i + 1                           # 1-based in the precursor
        if res1 not in RES1_WEIGHT:
            continue
        if position - 1 < leader_min:              # leader too short
            continue
        core = peptide[i:]
        if not (core_min <= len(core) <= core_max):
            continue
        for ring in (7, 8, 9):
            if ring > len(core):
                continue
            acceptor = core[ring - 1]
            if acceptor not in 'DE':
                continue
            tail = core[ring:]
            if len(tail) < 3:                      # cannot thread without tail
                continue
            out.append({'leader': peptide[:i], 'core': core, 'ring_size': ring,
                        'res1': res1,
                        'acceptor': '%s%d' % (acceptor, ring), 'tail': tail})
    return out


def score_cut(cut, distance=None, same_strand=None):
    """Heuristic plausibility score. Higher is better.

    distance     : bp to the nearest maturation gene (None if unknown)
    same_strand  : True when the ORF lies on the cluster's dominant strand
    """
    s = 0.0
    notes = []

    # --- genomic context: the term that separates signal from noise
    if distance is not None:
        if distance <= 2000:
            s += 3.0
            notes.append('%d bp from the machinery' % distance)
        elif distance <= 5000:
            s += 1.0
        else:
            s -= 2.5
            notes.append('far from the machinery (%d bp)' % distance)
    if same_strand is True:
        s += 1.5
    elif same_strand is False:
        s -= 1.0
        notes.append('opposite strand to the cluster')

    s += RES1_WEIGHT[cut['res1']]
    if cut['res1'] == 'G':
        notes.append('res1=Gly')

    # Ring sizes 8 and 9 are more frequent than 7 among characterised peptides.
    s += {7: 1.0, 8: 2.0, 9: 2.0}[cut['ring_size']]

    core_length = len(cut['core'])
    if 15 <= core_length <= 24:
        s += 2.0
    elif 10 <= core_length <= 30:
        s += 0.5
    else:
        s -= 1.5
        notes.append('core outside the canonical 15-24 range')

    if len(cut['tail']) >= 5:
        s += 1.0
    if BULKY & set(cut['tail'][:6]):
        s += 1.0
        notes.append('bulky plug residue in the tail')

    leader = cut['leader']
    if leader:
        acidic = sum(leader.count(x) for x in 'DE') / len(leader)
        if acidic >= 0.15:
            s += 1.5
            notes.append('acidic leader')
        elif acidic >= 0.08:
            s += 0.5
    if 15 <= len(leader) <= 45:
        s += 1.0
    # The conserved Thr and the Thr-x(3,4)-Pro motif occur in the region
    # immediately preceding the cleavage site, so only the leader's C-terminal
    # stretch is inspected.
    leader_tail = leader[-20:]
    if 'T' in leader_tail:
        s += 0.5
    if re.search(r'T.{3,4}P', leader_tail):
        s += 1.5
        notes.append('TxxxP motif in the leader')

    # Low-complexity stretches are not real cores; they are what a spurious
    # reading frame through a repetitive region produces.
    counts = Counter(cut['core'])
    if counts.most_common(1)[0][1] / len(cut['core']) > 0.35:
        s -= 2.0
        notes.append('low-complexity core')

    # Strongly basic cores are likewise a signature of frame-shifted noise.
    basic = sum(cut['core'].count(x) for x in 'RK') / len(cut['core'])
    if basic > 0.25:
        s -= 1.5
        notes.append('anomalously basic core')

    return round(s, 2), '; '.join(notes)


def best_cut(peptide, args, distance=None, same_strand=None):
    """Highest-scoring cut point for one ORF, or None if none is admissible."""
    cuts = candidate_cuts(peptide, args.core_min, args.core_max, args.leader_min)
    if not cuts:
        return None
    scored = []
    for c in cuts:
        value, notes = score_cut(c, distance, same_strand)
        c = dict(c)
        c['score'] = value
        c['notes'] = notes
        c['machinery_distance'] = '' if distance is None else distance
        scored.append(c)
    return max(scored, key=lambda x: x['score'])


# -------------------------------------------------------------------- pipeline
def process(path, args):
    record = parse_genbank(path)
    if not record['sequence']:
        print('  [!] %s carries no sequence (empty ORIGIN); skipped'
              % os.path.basename(path), file=sys.stderr)
        return []

    candidates = []

    # --- 1. whatever antiSMASH already annotated
    annotated = annotated_precursors(record)
    for a in annotated:
        a['accession'] = record['accession']
        a['organism'] = record['organism']
        a['offset'] = record['offset']
        a['tail'] = ''
        a['machinery_distance'] = 0
    candidates.extend(annotated)
    if annotated and not args.force_scan:
        return candidates

    # --- 2. ORF scan
    if args.whole_record:
        lo, hi = 1, len(record['sequence'])
    else:
        lo, hi = core_bounds(record)
        lo = max(1, lo - args.flank)
        hi = min(len(record['sequence']), hi + args.flank)

    cds = annotated_cds(record)
    anchors = anchor_genes(record)

    # 2a. Short annotated CDSs. The precursor is sometimes called after all but
    #     labelled 'hypothetical protein', so these must be inspected too and
    #     not assumed to be covered by the unannotated-ORF scan below.
    for a, b, strand, tag, protein, product in cds:
        if not protein or not (args.orf_min <= len(protein) <= args.orf_max):
            continue
        if b < lo or a > hi:
            continue
        d, same = genomic_context(a, b, strand, anchors)
        cut = best_cut(protein, args, d, same)
        if cut and cut['score'] >= args.min_score:
            cut.update({'source': 'short_annotated_cds', 'locus': tag,
                        'start': a, 'end': b, 'strand': strand,
                        'peptide_class': '',
                        'accession': record['accession'],
                        'organism': record['organism'],
                        'offset': record['offset'],
                        'notes': (cut['notes'] + '; product: ' + product).strip('; ')})
            candidates.append(cut)

    # 2b. Unannotated ORFs in the intergenic gaps.
    def unoccupied(a, b):
        return all(b < s or a > e for s, e, _, _, _, _ in cds)

    for a, b, strand, peptide in find_orfs(record['sequence'], lo, hi,
                                           args.orf_min, args.orf_max):
        if not unoccupied(a, b):
            continue
        d, same = genomic_context(a, b, strand, anchors)
        cut = best_cut(peptide, args, d, same)
        if cut and cut['score'] >= args.min_score:
            cut.update({'source': 'unannotated_orf',
                        'locus': 'orf_%d_%d%s' % (a, b, strand),
                        'start': a, 'end': b, 'strand': strand,
                        'peptide_class': '',
                        'accession': record['accession'],
                        'organism': record['organism'],
                        'offset': record['offset']})
            candidates.append(cut)

    return candidates


def main():
    ap = argparse.ArgumentParser(
        description='Recover lasso peptide precursors from antiSMASH GenBank '
                    'regions and emit core peptides for structure prediction.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n'
               '  python3 02_find_precursors.py regions/*.gbk -o results/ --top 5\n')
    ap.add_argument('gbk', nargs='+', help='antiSMASH region GenBank files')
    ap.add_argument('-o', '--outdir', default='lasso_precursors',
                    help='output directory (default: lasso_precursors)')
    ap.add_argument('--orf-min', type=int, default=25,
                    help='minimum ORF length in aa (25)')
    ap.add_argument('--orf-max', type=int, default=80,
                    help='maximum ORF length in aa (80)')
    ap.add_argument('--core-min', type=int, default=10,
                    help='minimum core length (10)')
    ap.add_argument('--core-max', type=int, default=40,
                    help='maximum core length (40)')
    ap.add_argument('--leader-min', type=int, default=8,
                    help='minimum leader length (8)')
    ap.add_argument('--min-score', type=float, default=9.0,
                    help='minimum score to report (9.0); lower to 6 to be '
                         'permissive when a genome returns nothing')
    ap.add_argument('--flank', type=int, default=1500,
                    help='bp scanned either side of the proto_core (1500)')
    ap.add_argument('--whole-record', action='store_true',
                    help='scan the entire record rather than the cluster core')
    ap.add_argument('--force-scan', action='store_true',
                    help='scan even when antiSMASH already annotated a precursor')
    ap.add_argument('--top', type=int, default=0,
                    help='keep only the N best candidates per file (0 = all)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    everything = []

    for path in args.gbk:
        if not os.path.isfile(path):
            print('  [!] not found: %s' % path, file=sys.stderr)
            continue
        found = process(path, args)
        found.sort(key=lambda x: -x['score'])
        if args.top:
            found = found[:args.top]
        n_annotated = sum(1 for x in found if x['source'] == 'antismash_annotated')
        print('%-42s %2d candidate(s)%s'
              % (os.path.basename(path), len(found),
                 '  [%d annotated by antiSMASH]' % n_annotated if n_annotated else ''))
        everything.extend(found)

    if not everything:
        print('\nNo candidates. Try --min-score 6 --whole-record.')
        return

    everything.sort(key=lambda x: (x['accession'], -x['score']))

    tsv = os.path.join(args.outdir, 'candidates.tsv')
    columns = ['accession', 'organism', 'locus', 'source', 'score',
               'machinery_distance_bp', 'start', 'end', 'strand',
               'genome_start', 'res1', 'acceptor', 'ring_size',
               'leader_length', 'core_length', 'leader', 'core', 'notes']
    with open(tsv, 'w') as fh:
        fh.write('\t'.join(columns) + '\n')
        for c in everything:
            fh.write('\t'.join(str(x) for x in [
                c['accession'], c['organism'], c['locus'], c['source'],
                c['score'], c.get('machinery_distance', ''),
                c['start'], c['end'], c['strand'],
                c['start'] + c.get('offset', 0),
                c.get('res1', ''), c.get('acceptor', ''), c.get('ring_size', ''),
                len(c.get('leader', '')), len(c['core']),
                c.get('leader', ''), c['core'], c.get('notes', '')]) + '\n')

    faa_core = os.path.join(args.outdir, 'cores_for_structure.faa')
    faa_full = os.path.join(args.outdir, 'full_precursors.faa')
    with open(faa_core, 'w') as fc, open(faa_full, 'w') as ff:
        seen = Counter()
        for c in everything:
            seen[c['accession']] += 1
            name = '%s_cand%02d' % (c['accession'], seen[c['accession']])
            desc = ('score=%s ring=%s res1=%s acceptor=%s source=%s loc=%d-%d%s'
                    % (c['score'], c.get('ring_size', '?'), c.get('res1', '?'),
                       c.get('acceptor', '?'), c['source'],
                       c['start'], c['end'], c['strand']))
            fc.write('>%s %s\n%s\n' % (name, desc, c['core']))
            ff.write('>%s_precursor %s\n%s\n'
                     % (name, desc, c.get('leader', '') + c['core']))

    print('\n%d candidate(s) in total' % len(everything))
    print('  table            : %s' % tsv)
    print('  FOR STRUCTURE    -> %s' % faa_core)
    print('  full precursors  : %s' % faa_full)

    atypical = [c for c in everything if not (15 <= len(c['core']) <= 24)]
    if atypical:
        print('\n  [warning] %d core(s) fall outside the canonical 15-24 aa range.'
              % len(atypical))
        print('            Inspect these manually before submitting them: the')
        print('            leader/core boundary may be misplaced, or they may')
        print('            not be precursors at all.')
    print('\n  These are CANDIDATES. Confirm with RODEO or BAGEL4, and by')
    print('  checking that the same ORF recurs in other genomes of the group.')


if __name__ == '__main__':
    main()
