// node check_marker.js MARKER : report region, alleles and simulated GT1/ME49/VEG calls
const M = require('./src/core.js'), fs = require('fs'), m = process.argv[2];
const panel = M.buildPanel(JSON.parse(fs.readFileSync('src/panel.json')), JSON.parse(fs.readFileSync('src/primers.json')));
const mk = panel.markers[m];
console.log(m, 'region', mk.window[1] - mk.window[0] + 1, 'sites', mk.sites.length, 'products', mk.primerSummary.products, 'primers found', mk.primerSummary.found + '/' + mk.primerSummary.total, '\nalleles (' + mk.alleles.length + '):', mk.alleles.map(a => mk.refs.filter(r => r.rep === a).sort((x, y) => x.order - y.order).map(r => r.name).join('=')).join(' | '));
const rd = (n, s) => ({ fileName: n, recName: n, format: 'text', seq: s, qual: null, traces: null, ploc: null }); const reads = [];
for (const name of ['GT1', 'ME49', 'VEG']) { const r = mk.refs.find(x => x.name === name); if (!r || !r.cut) continue; const [a, b] = r.cut, amp = r.seq.slice(Math.max(0, a - 10), b + 10); reads.push(rd('Z_' + name + ' F.ab1', amp), rd('Z_' + name + ' R.ab1', M.revcomp(amp))); }
for (const r of M.analyze(reads, panel).results) console.log(' sim', r.sample, r.marker, '->', r.call, r.status);
