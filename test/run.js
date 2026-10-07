const M = require('../src/core.js'); const fs = require('fs'), path = require('path');
const dir = process.argv[2];
const reads = [];
for (const f of fs.readdirSync(dir)) { const b = fs.readFileSync(path.join(dir, f)); reads.push(...M.parseFile(f, b.buffer.slice(b.byteOffset, b.byteOffset + b.length))); }
const refsRaw = JSON.parse(fs.readFileSync(path.join(__dirname, '../src/panel.json')));
const panel = M.buildPanel(refsRaw);
const mk = panel.markers.PK1; console.log('window', mk.window, 'sites', mk.sites.length, mk.refs.map(r => r.name + ':' + r.win));
const t = Date.now(); const out = M.analyze(reads, panel); console.log('ms', Date.now() - t);
for (const r of out.results) {
  console.log(`\n${r.sample} [${r.marker}] -> ${r.call} (${r.status})  files=${r.files.length} len=${r.stats && r.stats.length} conflicts=${r.stats && r.stats.conflicts} indelc=${r.stats && r.stats.indelConflicts}`);
  if (r.comps) console.log('   ', r.comps.map(c => `${c.ref.name}:${c.events}`).join('  '));
  if (r.best) console.log('    best diffs:', r.best.diffs.map(d => `${d.anchorPos + 1}${d.ref}>${d.cons}(x${d.cov ?? '-'})`).join(' '));
  r.flags.forEach(f => console.log('    [' + f.level + ']', f.text));
}
console.log('problems', out.problems);
console.log('refs detected', out.reads.filter(r => r.isRef).map(r => r.fileName));
