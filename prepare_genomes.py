#!/usr/bin/env python
# (C) 2026 Tom Eulenfeld, MIT License

import json
from sugar import read


PATH = 'data'


def download_genomes_assign_primers():
    with open(f'{PATH}/organisms.json') as f:
        organisms = json.load(f)
    for org in organisms:
        print(f"Download {org['long_label']} genome...")
        seqs = read(org['url'])
        seqs.write(f"{PATH}/genome_{org['name']}.fasta")


def assign_primers():
    with open(f'{PATH}/organisms.json') as f:
        organisms = json.load(f)
    primers = {}
    for org in organisms:
        print(f"Assign primers for {org['long_label']}...")
        oname = org['name']
        olabel = org['label']
        seqs = read(f"{PATH}/genome_{org['name']}.fasta")
        print('Genome length:', sum(len(seq) for seq in seqs))
        primers[olabel] = {}
        for plen in (10, 20, 50, 100):
            primer = seqs[0][10:10+plen]
            assert len(primer) == plen
            n1 = plen // 2
            n2 = plen // 2 + 2
            primer_mm = primer[:n1] + primer[n1:n2].complement() + primer[n2:]
            primers[olabel][f'primer_{oname}_{plen}'] = str(primer)
            primers[olabel][f'primer_{oname}_{plen}_mm2'] = str(primer_mm)
    with open(f'{PATH}/primers.json', 'w') as f:
        json.dump(primers, f, indent=2)


if __name__ == '__main__':
    download_genomes_assign_primers()
    assign_primers()
