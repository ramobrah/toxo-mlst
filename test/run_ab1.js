const M = require('../src/core.js'); const fs = require('fs'), path = require('path');
const dir = process.argv[2]; const reads = [];
for (const f of fs.readdirSync(dir)) { const b = fs.readFileSync(path.join(dir, f)); reads.push(...M.parseFile(f, b.buffer.slice(b.byteOffset, b.byteOffset + b.length))); }
console.log(reads.map(r => `${r.fileName}: ${r.format} len=${r.seq.length} q=${!!r.qual} tr=${!!r.traces}`).slice(0,3).join('\n'));
const panel = M.buildPanel(JSON.parse(fs.readFileSync(path.join(__dirname, '../src/panel.json'))));
const out = M.analyze(reads, panel);
for (const r of out.results) {
  console.log(`\n${r.sample} -> ${r.call} (${r.status}) len=${r.stats.length} conflicts=${r.stats.conflicts} mixed=${r.stats.mixed}`);
  r.used.forEach(u => console.log('   ', u.fileName, 'trim', u.trimmed, 'peaks', u.peaks.length, 'rev', u.reversed));
  r.flags.forEach(f => console.log('    [' + f.level + ']', f.text));
}
