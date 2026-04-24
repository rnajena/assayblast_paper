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
from sugar import BioSeq, read_fts


PATH = 'data'
TMP = 'data'


def compare_versions(data=None):
    with open(f'{PATH}/primers.json') as f:
        primers = json.load(f)
    if data is None:
        data = {'header': ['organism', 'plen', 'mismatch', 'version', 'success', 'nhits', 'time', 'evalue'],
                'rows': []}
    for olabel in primers:
        for primer_name, primer in primers[olabel].items():
            for mm in (0, 2):
                oname = primer_name.split('_')[1]
                if (oname == 'human' and len(primer) in (10, 20, 50, 100) or
                        oname == 'yeast' and len(primer) in (10, 20, 50) or
                        oname == 'bacteria' and len(primer) in (10, 20) or
                        oname == 'virus' and len(primer) in (10, 20)) and 'mm2' not in primer_name:
                    query = f'{{path}}/tmp_primer_{oname}_{primer_name}.fasta'
                    genome = f'{{path}}/genome_{oname}.fasta'
                    out = f'{PATH}/blast_evaluation_{{version}}_{primer_name}_mm{mm}.tsv'
                    BioSeq(primer, id=f'{oname}_{primer_name}').write(query.format(path=PATH))
                    for version in ('v2', 'v1_e1k', 'v1'):
                        print(f'Run assayBLAST {version} for {oname} with {primer_name}...')
                        outv = out.format(version=version)
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
                        time = round(perf_counter() - t0, 2)
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
                        print(f'{version} finished in {time} seconds, success: {success}')
                        data['rows'].append((olabel, len(primer), mm, version, success, nhits, time, evalue))
                        with open(f'{PATH}/evaluation.json', 'w') as f:
                            json.dump(data, f, indent=2)
                    os.remove(query.format(path=PATH))


def print_table():
    with open(f'{PATH}/evaluation.json') as f:
        data = json.load(f)
    print(r'\begin{tabular}{lrrr rcr rcr rcr}')
    print(r'\toprule')
    print('% ' + ' '.join(data['header']))
    print(r'\midrule')
    results = defaultdict(lambda: defaultdict(list))
    evaluev2 = {}
    nhitsv2 = {}
    for row in data['rows']:
        olabel, plen, mismatch, version, success, nhits, time, evalue = row
        results[olabel, plen, mismatch][version] = (success, nhits, time, evalue)
        if evalue:
            evaluev2[olabel, plen, mismatch] = evalue
            nhitsv2[olabel, plen, mismatch] = nhits
    # Sort by e-value
    results = dict(sorted(results.items(), key=lambda item: evaluev2.get(item[0])))
    for (olabel, plen, mismatch), vresults in results.items():
        def tstr(time):
            if time < 10:
                d = 2
            elif time < 100:
                d = 1
            else:
                d = 0
            return rf'${time:.{d}f}\,\mathrm{{s}}$'
        def mark(nhits):
            if nhits == nhitsv2[olabel, plen, mismatch]:
                return r'\cmark'
            else:
                return r'\xmark'
        def fmt_evalue(evalue):
            if evalue < 1e-8 or evalue > 1e8:
                return rf'{evalue:.1e}'.replace('+', '')
            elif evalue >= 10:
                return rf'${evalue:_.0f}$'.replace('_', r'\,')
            else:
                return rf'${evalue:.2g}$'

        print(rf'\textit{{{olabel}}} & ${plen}\,\mathrm{{nt}}$ & {mismatch} & ' +
                fmt_evalue(evaluev2[olabel, plen, mismatch]) + ' & ' +
                ' & '.join(rf'{nhits} & {mark(nhits)} & {tstr(time)}'
                            for version, (success, nhits, time, evalue) in vresults.items()) +
                r'\\')
    print(r'\bottomrule')
    print(r'\end{tabular}')


if __name__ == '__main__':
    # run twice to exclude the timing for building the BLAST dbs
    compare_versions()
    print_table()
