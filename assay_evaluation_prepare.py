#!/usr/bin/env python
# (C) 2026 Tom Eulenfeld, MIT License

import pandas as pd
from sugar.web import Entrez
from sugar import BioBasket


def download_genomes():
    client = Entrez()
    assay = pd.read_csv('assay/microarray_results.txt', sep='\t')
    seqs = []
    for sid in assay['probe']:
        print(f'Download {sid}')
        if '-' in sid:
            sid1, sid2 = sid.split('-')
            sid2 = sid1[:-3] + sid2
            seq = client.get_seq(sid1, overwrite=True)
            seq2 = client.get_seq(sid2, overwrite=True)
            print(f'Join seqs with lengths {len(seq)}, {len(seq2)}')
            seq = seq + seq2
            seq.id = sid
        else:
            seq = client.get_seq(sid, overwrite=True)
        seqs.append(seq)
    seqs = BioBasket(seqs)
    print(seqs)
    seqs.write('assay/genomes.fasta')


if __name__ == '__main__':
    download_genomes()