#!/usr/bin/env python
# (C) 2026 Tom Eulenfeld, MIT License

import matplotlib.pyplot as plt
import pandas as pd


def map_probe_names(syn, exp):
    exp_probes = list(exp.columns)
    syn_probes = list(syn.columns)
    m = {}
    # first map probes with the same name
    for p in syn_probes.copy():
        p = p.removeprefix('probe_')
        for p2 in exp_probes.copy():
            if p == p2:
                m[p2] = 'probe_' + p
                exp_probes.remove(p2)
                syn_probes.remove('probe_' + p)
                break
    # map probes via prefixes, and syn probe names longer than 2 chars
    # each syn probe may match several exp probes
    for p in syn_probes.copy():
        p = p.removeprefix('probe_')
        for p2 in exp_probes.copy():
            if len(p) > 2 and p2.startswith(p):
                m[p2] = 'probe_' + p
                print(f'Match probe_{p} to {p2}')
                exp_probes.remove(p2)
                try:
                    syn_probes.remove('probe_' + p)
                except ValueError:
                    print(f'probe_{p} twice in syn_probes')
    # map remaining probes via prefixes
    # each syn probe may match several exp probes
    for p in syn_probes.copy():
        p = p.removeprefix('probe_')
        for p2 in exp_probes.copy():
            if p2.startswith(p):
                m[p2] = 'probe_' + p
                print(f'Match probe_{p} to {p2}')
                exp_probes.remove(p2)
                try:
                    syn_probes.remove('probe_' + p)
                except ValueError:
                    print(f'probe_{p} twice in syn_probes')
    # rename syn probes
    for col in syn.columns:
        for k, v in m.items():
            if v == col:
                syn[k] = syn[col]
        del syn[col]
    # drop a single exp probes
    print(f'Left-over exp probes, remove them: {exp_probes}')
    exp.drop(exp_probes, axis=1, inplace=True)
    return m


def plot_confusion_matrix(syn, exp):
    syn_pos = syn.apply(lambda x: x.str.startswith('lin'))
    exp_pos = exp.apply(lambda x: x > 0.5)

    print('number of values:' , exp.size)
    pp = (exp_pos & syn_pos).values.sum()
    nn = (~exp_pos & ~syn_pos).values.sum()
    pn = (exp_pos & ~syn_pos).values.sum()
    np = (~exp_pos & syn_pos).values.sum()

    print('         syn')
    print('         + -')
    print('exp +', pp, pn)
    print('    -', np, nn)

    plt.figure(figsize=(3, 3))
    plt.imshow([[pp, pn], [np, nn]], cmap='Blues', vmin=0)
    plt.annotate(pp, (0, 0,), ha='center', va='center')
    plt.annotate(pn, (1, 0,), ha='center', va='center')
    plt.annotate(np, (0, 1,), ha='center', va='center')
    plt.annotate(nn, (1, 1,), ha='center', va='center', color='w')
    plt.colorbar(shrink=0.6)
    plt.xticks([0, 1], ['positive', 'negative'])
    plt.yticks([0, 1], ['positive', 'negative'], rotation=90, va='center')
    plt.xlabel('AssayBLAST v2')
    plt.ylabel('Microarray')
    plt.savefig('fig_confusion_matrix.pdf', bbox_inches='tight', pad_inches=0.02)


if __name__ == '__main__':
    exp = pd.read_csv('assay/microarray_results.txt', sep='\t')
    syn = pd.read_csv('assay/blast_results_assay_overview.tsv', sep='\t')

    exp = exp.rename(columns={'probe': 'Genome'}).set_index('Genome')
    syn = syn.set_index('Genome')

    m = map_probe_names(syn, exp)
    exp = exp.sort_index(axis=1)
    syn = syn.sort_index(axis=1)

    # check that only linear or none
    assert syn.apply(lambda x: x.str.startswith('exp') + x.str.startswith('const')).values.sum() == 0

    plot_confusion_matrix(syn, exp)
