#!/usr/bin/env python
# (C) 2026 Tom Eulenfeld, MIT License

from collections import defaultdict
import json
import os
from subprocess import run
import subprocess
import shutil
from time import perf_counter
import re

from assay_analyze import _preprocess_hits
from assay_blast import run_blast
from sugar import BioSeq, Feature, FeatureList, read, read_fts


PATH = 'data'
TMP = 'data'


def compare_versions(data=None):
    with open(f'{PATH}/queries.json') as f:
        queries = json.load(f)
    if data is None:
        data = {'header': ['organism', 'qlen', 'mismatch', 'version', 'success', 'nhits', 'time', 'evalue'],
                'rows': []}
    for olabel in queries:
        for query_name, queryseq in queries[olabel].items():
            for mm in (0, 2):
                oname = query_name.split('_')[1]
                if (oname == 'human' and len(queryseq) in (10, 20, 50, 100) or
                        oname == 'yeast' and len(queryseq) in (10, 20, 50) or
                        oname == 'bacteria' and len(queryseq) in (10, 20) or
                        oname == 'virus' and len(queryseq) in (10, 20)) and 'mm2' not in query_name:
                    query = f'{{path}}/tmp_query_{oname}_{query_name}.fasta'
                    genome = f'{{path}}/genome_{oname}.fasta'
                    out = f'{PATH}/blast_evaluation_{{version}}_{query_name}_mm{mm}.tsv'
                    BioSeq(queryseq, id=f'primer_{oname}_{query_name}').write(query.format(path=PATH))
                    for version in ('v2', 'v1_e1k', 'v1'):
                        print(f'Run assayBLAST {version} for {oname} with {query_name}...')
                        outv = out.format(version=version)
                        runtime = 0.
                        num_runs = 0
                        while runtime < 10:
                            t0 = perf_counter()
                            if version == 'v2':
                                run_blast(
                                    query.format(path=PATH),
                                    [genome.format(path=PATH)],
                                    outv,
                                    db=f'{TMP}/db_tmp/{oname}', keep_db=True, mismatch=mm)
                            else:
                                # if os.path.exists(f'{TMP}/db_tmp'):
                                #     shutil.rmtree(f'{TMP}/db_tmp')
                                suf = version[2:]
                                path = '../' + PATH
                                call = f"python assayBLAST{suf}.py -q {query.format(path=path)} -db {TMP}/db_tmp/{oname} -m {mm} -g {genome.format(path=path)}"
                                print(call)
                                try:
                                    run(call.split(), check=True, cwd='AssayBLASTv1')
                                    shutil.copyfile('AssayBLASTv1/blast_results.tsv', outv)
                                except subprocess.CalledProcessError as e:
                                    print(f"Error occurred while running {version}: {e}")
                            runtime += perf_counter() - t0
                            num_runs += 1
                        runtime_mean = round(runtime / num_runs, 2)
                        if version == 'v2':
                            fts = read_fts(outv)
                            _preprocess_hits(fts)  # correctly assign number of mismatches
                            fts = fts.select(type_in=('primer', 'probe'), mismatch_le=mm)
                            nhits = len(fts)
                            fts.write(outv + '.filtered.tsv', 'tsv', keys='type seqid name start stop strand mismatch')
                            with open(outv) as f:
                                text = f.read()
                            match = re.search(r'evalue ([0-9.eE+-]+) from option', text)
                            evalue = float(f'{float(match.group(1)):.2g}')
                            success = nhits > 0
                        else:
                            if os.path.exists(outv):
                                with open(outv) as f:
                                    text = f.read()
                                nhits = text.count('pos')
                                success = 'genome' in text
                            else:
                                nhits = -1
                                success = False
                            evalue = None
                        print(f'{version} finished in {runtime_mean} seconds, success: {success}')
                        data['rows'].append((olabel, len(queryseq), mm, version, success, nhits, runtime_mean, evalue))
                        with open(f'{PATH}/evaluation.json', 'w') as f:
                            json.dump(data, f, indent=2)
                    #os.remove(query.format(path=PATH))

def _fmt_evalue(evalue):
    if evalue < 1e-6 or evalue > 1e8:
        return rf'{evalue:.1e}'.replace('+', '')
    elif evalue >= 10:
        return rf'${evalue:_.0f}$'.replace('_', r'\,')
    else:
        return rf'${evalue:.2g}$'


def print_table():
    with open(f'{PATH}/evaluation.json') as f:
        data = json.load(f)
    print('\nevaluation table 4')
    print(r'\begin{tabular}{lrrr rcr rcr rcr}')
    print(r'\toprule')
    print('% ' + ' '.join(data['header']))
    print(r'\midrule')
    results = defaultdict(lambda: defaultdict(list))
    evaluev2 = {}
    nhitsv2 = {}
    for row in data['rows']:
        olabel, qlen, mismatch, version, success, nhits, time, evalue = row
        results[olabel, qlen, mismatch][version] = (success, nhits, time, evalue)
        if evalue:
            evaluev2[olabel, qlen, mismatch] = evalue
            nhitsv2[olabel, qlen, mismatch] = nhits
    # Sort by e-value
    results = dict(sorted(results.items(), key=lambda item: evaluev2.get(item[0])))
    for (olabel, qlen, mismatch), vresults in results.items():
        def tstr(time):
            if time < 10:
                d = 2
            elif time < 100:
                d = 1
            else:
                d = 0
            return rf'${time:.{d}f}\,\mathrm{{s}}$'
        def mark(nhits, version):
            if nhits == nhitsv2[olabel, qlen, mismatch]:
                if olabel == 'Candida albicans' and qlen == 20 and mismatch == 2 and version=='v1':
                    # both versions find 4 hits, but v1 finds one incorrect hit with gaps
                    return r'\xmark'
                return r'\cmark'
            else:
                return r'\xmark'
        asterisk = '$^*$' if olabel == 'Candida albicans' and qlen == 20 and mismatch == 2 else ''
        print(rf'{asterisk}\textit{{{olabel}}} & ${qlen}\,\mathrm{{nt}}$ & {mismatch} & ' +
                _fmt_evalue(evaluev2[olabel, qlen, mismatch]) + ' & ' +
                ' & '.join(rf'{nhits} & {mark(nhits, version)} & {tstr(time)}'
                            for version, (success, nhits, time, evalue) in vresults.items()) +
                r'\\')
    print(r'\bottomrule')
    print(r'\end{tabular}')


def parse_v1_results(fname, filename_as_id=True):
    with open(fname) as f:
        text = f.read()
    fts = []
    query_names = None
    for line in text.splitlines():
        if query_names is None:
            query_names = line.strip().split('\t')
        elif (filename_as_id and '|' in line) or (not filename_as_id and '.fna' in line):
            id_, *hits = line.strip().split('\t')
            if filename_as_id:
                seqid = id_.split('|')[1]
            else:
                seqid = id_.removesuffix('.fna')
            for fts_str, qname in zip(hits, query_names):
                for ft_str in fts_str.split(';'):
                    if ft_str.strip():
                        start, stop = ft_str.removesuffix(')').split()[-1].split('-')
                        strand = '-' if 'rev' in qname else '+'
                        name = qname.rsplit('_', maxsplit=1)[0]
                        type_= name.split('_')[0]
                        fts.append(Feature(type_, meta=dict(seqid=seqid, name=name), start=int(start)-1, stop=int(stop), strand=strand))
    return FeatureList(fts)


def _filter_fts(fts1, fts2, inverse=False):
    locs = {ft.loc for ft in fts2}
    if inverse:
        fts1.data = [ft for ft in fts1 if ft.loc not in locs]
    else:
        fts1.data = [ft for ft in fts1 if ft.loc in locs]
    return fts1


def _add_evalue_v1(fts):
    from assay_blast import OUTFMT7
    out = 'blast_evaluation_v1_query_virus_20_mm2_full.tsv'
    # call BLAST with parameters used in v1
    call = (f'blastn -query {PATH}/tmp_query_yeast_query_yeast_20.fasta -db {PATH}/db_tmp/yeast '
            '-evalue 100000 -gapopen 10 -gapextend 6 -reward 5 -penalty -4 -word_size 7 '
            f"-outfmt '{OUTFMT7}' -dust no -out {TMP}/{out}")
    os.system(call)
    fts2 = read_fts(f'{TMP}/{out}')
    # print(f'all potential v1 hits')
    # print(fts2.tostr(w=0))
    # only return the features that are actually reported in v1
    fts2 = _filter_fts(fts2, fts)
    fts.data = fts2.data
    return fts


def _merge_fts(fts_v2, fts_v1):
    fts = fts_v2.copy()
    fts_v1_locs = {ft.loc for ft in fts_v1}
    for ft in fts:
        ft.meta.v2 = 'yes'
        ft.meta.v2m = r'\cmark'
        ft.meta.v1 = 'yes' if ft.loc in fts_v1_locs else 'no'
        ft.meta.v1m = r'\cmark' if ft.loc in fts_v1_locs else r'\xmark'
    fts_v1_incorrect = _filter_fts(fts_v1, fts_v2, inverse=True)
    for ft in fts_v1_incorrect:
        ft.meta.v2 = 'no'
        ft.meta.v2m = r'\cmark'
        ft.meta.v1 = 'yes'
        ft.meta.v1m = r'\xmark'
    fts.extend(fts_v1_incorrect)
    return fts


def _add_aliscores(fts, seqs, query):
    seqs2 = []
    for ft, seq in zip(fts, seqs):
        if len(seq) == 19:
            seq += '.'
        elif len(seq) == 18:
            seq = seq[:9] + '--' + seq[9:]
        seqs2.append(seq)
        score = sum(5 if a == b else -8 if a == '-' else 0 if a == '.' else -4 for a, b in zip(seq.data, query.data, strict=True))
        ft.meta.aliscore = score
    seqs3 = []
    for seq in seqs2:
        fancy_seq = ''.join(a if a == b else
                            r'\textcolor{red!70!black}{-}' if a == '-' else
                            rf'\textcolor{{white!70!black}}{{{a}}}'
                            for a, b in zip(seq.data, query.data, strict=True))
        seqs3.append(fancy_seq)
    return seqs3

REPORT = 'https://www.ncbi.nlm.nih.gov/nuccore/{id}?report=fasta&from={start}&to={stop}&strand={rc}'

def analyze_yeast_20nt_mm2():
    seqs = read('data/genome_yeast.fasta')
    query = read('data/tmp_query_yeast_query_yeast_20.fasta')[0]
    fts_v2 = read_fts('data/blast_evaluation_v2_query_yeast_20_mm2.tsv.filtered.tsv')
    fts_v2_with_evalue = read_fts('data/blast_evaluation_v2_query_yeast_20_mm2.tsv')
    fts_v2 = _filter_fts(fts_v2_with_evalue, fts_v2)
    fts_v1 = parse_v1_results('data/blast_evaluation_v1_query_yeast_20_mm2.tsv')
    fts_v2.sort()
    fts_v1.sort()
    _add_evalue_v1(fts_v1)
    print('\nv2 hits features')
    print(fts_v2)
    print('v1 hits features')
    print(fts_v1)
    print('query')
    print(query)
    print('v2 hits sequences')
    print(seqs[fts_v2])
    print('v1 hits sequences')
    print(seqs[fts_v1])
    print('\nhits table 4')
    # merge features
    fts = _merge_fts(fts_v2, fts_v1)
    fancy_seqs = _add_aliscores(fts, seqs[fts], query)
    print('table')
    print(rf'query & \texttt{{{query}}}\\ \midline')
    for ft, seq in zip(fts, fancy_seqs):
        id_ = ft.seqid.replace('_', r'\_')
        locstr = f'{id_}:{ft.loc.start+1}-{ft.loc.stop} ({ft.loc.strand})'
        locurl = REPORT.format(id=ft.seqid, start=ft.loc.start+1, stop=ft.loc.stop, rc='true' if ft.loc.strand == '-' else 'false')
        print((rf"\href{{{locurl}}}{{\texttt{{{locstr}}}}} & \texttt{{{seq}}} "
               rf"& {ft.meta.aliscore} & {_fmt_evalue(ft.meta.evalue)} & {ft.meta.v2} & {ft.meta.v2m} & {ft.meta.v1} & {ft.meta.v1m}\\"))


if __name__ == '__main__':
    # run twice to exclude the timing for building the BLAST dbs
    compare_versions()
    print_table()
    analyze_yeast_20nt_mm2()
