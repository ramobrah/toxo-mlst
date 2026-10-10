"""Add ToxoDB strain sequences for one marker to src/panel.json. Usage: python3 add_marker.py MARKER file.fasta 'note text'"""
import json, sys
from Bio import SeqIO
marker, path, desc = sys.argv[1], sys.argv[2], sys.argv[3]
panel = [p for p in json.load(open('src/panel.json')) if p['marker'] != marker]
order = ['GT1', 'ME49', 'VEG']; key = lambda n: (order.index(n) if n in order else 9, n)
recs = sorted([(r.id[2:].split('_')[0], r) for r in SeqIO.parse(path, 'fasta')], key=lambda x: key(x[0]))
for st, r in recs: panel.append({'marker': marker, 'name': st, 'note': f'ToxoDB {r.id} (strain {st}), {desc}', 'seq': str(r.seq).upper()})
json.dump(panel, open('src/panel.json', 'w'), indent=0)
print(marker, len(recs), 'references added')
