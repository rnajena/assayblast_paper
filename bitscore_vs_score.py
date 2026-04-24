#!/usr/bin/env python
# (C) 2026 Tom Eulenfeld, MIT License

import json
from sugar import BioSeq, read_fts
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from assay_blast import run_blast, _score
import seaborn as sns
from scipy.stats import linregress


PATH = 'data'
TMP = '/tmp'


def call_blast():
    with open(f'{PATH}/primers.json') as f:
        primers = json.load(f)
    for olabel in primers:
        for primer_name, primer in primers[olabel].items():
            oname = primer_name.split('_')[1]
            print(f'Run assayBLAST for {olabel} with {primer_name}...')
            BioSeq(primer, id=primer_name).write(f"{TMP}/tmp_{primer_name}.fasta")
            run_blast(
                f'{TMP}/tmp_{primer_name}.fasta',
                [f'{PATH}/genome_{oname}.fasta'],
                f'{PATH}/blast_bitscore_{primer_name}.tsv',
                db=f'{TMP}/db/{oname}', keep_db=True,
                mismatch=2 if 'mm2' in primer_name else 0,
                max_target_seqs=1)
            os.remove(f'{TMP}/tmp_{primer_name}.fasta')


def _adapt_axis(df):
    plt.xlim(0, None)
    plt.ylim(0, None)
    plt.xlabel('score $S$')
    plt.ylabel('bit score $S_2$')
    lin = linregress(df['score'], df['bit score'])
    x = np.array([20, 520])
    plt.plot(x, lin.intercept + lin.slope * x, color='black', alpha=0.2, linestyle='--', zorder=-100)
    plt.annotate(f'$S_2={lin.slope:.4f}\\,S + {lin.intercept:.3f}$\n$R={lin.rvalue:.3f}$', xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top')


def plot_bitscore():
    fts = read_fts(f'{PATH}/blast_bitscore_*.tsv')
    ndata = {}
    for ft in fts:
        if ft.name not in ndata and ft.meta._blast.mismatch == (2 if ft.name.endswith('mm2') else 0):
            ndata[ft.name] = ft
    fts.data = list(ndata.values())
    for ft in fts:
        assert ft.meta._blast.mismatch == (2 if ft.name.endswith('mm2') else 0)
    for ft in fts:
        ft.meta['primer length'] = ft.meta._blast.qlen
        ft.meta.mismatch = ft.meta._blast.mismatch
        ft.meta.organism = ft.name.split('_')[1]
    df = fts.topandas(['organism', 'primer length', 'mismatch', 'score'])
    df['bit score'] = df['score']
    df['score'] = df.apply(lambda row: _score(row['primer length'], row['mismatch']), axis=1)
    sns.relplot(df, x='score', y='bit score', hue='organism', style='organism', height=3, aspect=1.6)
    _adapt_axis(df)
    plt.savefig('fig_bitscore_vs_score_organisms.pdf')
    sns.relplot(df, x='score', y='bit score', hue='primer length', style='mismatch',
                height=3, aspect=1.6, hue_norm=mpl.colors.LogNorm())
    _adapt_axis(df)
    plt.savefig('fig_bitscore_vs_score.pdf', bbox_inches='tight', pad_inches=0.05)


if __name__ == '__main__':
    call_blast()
    plot_bitscore()
