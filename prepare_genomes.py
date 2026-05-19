#!/usr/bin/env python
# (C) 2026 Tom Eulenfeld, MIT License

import json
from sugar import read


PATH = 'data'


def download_genomes_assign_queries():
    with open(f'{PATH}/organisms.json') as f:
        organisms = json.load(f)
    for org in organisms:
        print(f"Download {org['long_label']} genome...")
        seqs = read(org['url'])
        seqs.write(f"{PATH}/genome_{org['name']}.fasta")


def assign_queries():
    with open(f'{PATH}/organisms.json') as f:
        organisms = json.load(f)
    queries = {}
    for org in organisms:
        print(f"Assign queries for {org['long_label']}...")
        oname = org['name']
        olabel = org['label']
        seqs = read(f"{PATH}/genome_{org['name']}.fasta")
        print('Genome length:', sum(len(seq) for seq in seqs))
        queries[olabel] = {}
        for plen in (10, 20, 50, 100):
            query = seqs[0][10:10+plen]
            assert len(query) == plen
            n1 = plen // 2
            n2 = plen // 2 + 2
            query_mm = query[:n1] + query[n1:n2].complement() + query[n2:]
            queries[olabel][f'query_{oname}_{plen}'] = str(query)
            queries[olabel][f'query_{oname}_{plen}_mm2'] = str(query_mm)
    with open(f'{PATH}/queries.json', 'w') as f:
        json.dump(queries, f, indent=2)


if __name__ == '__main__':
    download_genomes_assign_queries()
    assign_queries()
