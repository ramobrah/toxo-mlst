/* MLST allele caller – core analysis (no DOM). Runs in browser and in Node for testing. */
(function (root) {
'use strict';

/* ---------- IUPAC helpers ---------- */
const IUPAC = { A:'A', C:'C', G:'G', T:'T', R:'AG', Y:'CT', S:'CG', W:'AT', K:'GT', M:'AC',
  B:'CGT', D:'AGT', H:'ACT', V:'ACG', N:'ACGT' };
const CODE_OF = {}; for (const [k, v] of Object.entries(IUPAC)) CODE_OF[v] = k;
const COMP = { A:'T', T:'A', C:'G', G:'C', R:'Y', Y:'R', S:'S', W:'W', K:'M', M:'K', B:'V', V:'B', D:'H', H:'D', N:'N', '-':'-' };
const isAmb = b => b in IUPAC && IUPAC[b].length > 1 && b !== 'N';
const setOf = b => IUPAC[b] || '';
const codeFor = bases => { const s = [...new Set(bases.split(''))].sort((a, b) => 'ACGT'.indexOf(a) - 'ACGT'.indexOf(b)).join(''); return CODE_OF[s] || 'N'; };
function compatible(a, b) { const sa = setOf(a), sb = setOf(b); for (const c of sa) if (sb.includes(c)) return true; return false; }
const revcomp = s => s.split('').reverse().map(c => COMP[c] || 'N').join('');
function cleanSeq(s) { return s.toUpperCase().replace(/U/g, 'T').replace(/[^ACGTRYSWKMBDHVN]/g, ''); }

/* ---------- File parsing ---------- */
function parseABIF(buf) {
  const dv = new DataView(buf);
  const u8 = new Uint8Array(buf);
  const str = (o, n) => String.fromCharCode(...u8.subarray(o, o + n));
  if (str(0, 4) !== 'ABIF') throw new Error('Not an ABIF file');
  const dirCount = dv.getInt32(18), dirOff = dv.getInt32(26);
  const tags = {};
  for (let i = 0; i < dirCount; i++) {
    const o = dirOff + i * 28;
    const name = str(o, 4), num = dv.getInt32(o + 4), elType = dv.getInt16(o + 8), elSize = dv.getInt16(o + 10);
    const numEl = dv.getInt32(o + 12), dataSize = dv.getInt32(o + 16);
    const dataOff = dataSize <= 4 ? o + 20 : dv.getInt32(o + 20);
    tags[name + num] = { elType, elSize, numEl, dataSize, dataOff };
  }
  const chars = t => t ? str(t.dataOff, t.numEl) : null;
  const bytes = t => t ? Array.from(u8.subarray(t.dataOff, t.dataOff + t.numEl)) : null;
  const shorts = t => { if (!t) return null; const a = new Array(t.numEl); for (let i = 0; i < t.numEl; i++) a[i] = dv.getInt16(t.dataOff + 2 * i); return a; };
  const pstr = t => { if (!t) return null; const n = u8[t.dataOff]; return str(t.dataOff + 1, n); };
  const seq = chars(tags.PBAS2 || tags.PBAS1);
  if (!seq) throw new Error('No base calls (PBAS) in file');
  const qual = bytes(tags.PCON2 || tags.PCON1);
  const ploc = shorts(tags.PLOC2 || tags.PLOC1);
  const order = chars(tags.FWO_1) || 'GATC';
  const tr = {};
  for (let k = 0; k < 4; k++) { const t = tags['DATA' + (9 + k)]; if (t) tr[order[k]] = shorts(t); }
  const hasTraces = Object.keys(tr).length === 4 && ploc && ploc.length === seq.length;
  return { seq: cleanSeq(seq), qual: qual && qual.length === seq.length ? qual : null,
    traces: hasTraces ? tr : null, ploc: hasTraces ? ploc : null, sampleName: pstr(tags.SMPL1) };
}
function parseFasta(text) {
  const recs = []; let cur = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith('>')) { cur = { header: line.slice(1).trim(), seq: '' }; recs.push(cur); }
    else if (cur) cur.seq += line.trim();
    else if (line.trim()) { cur = { header: '', seq: line.trim() }; recs.push(cur); }
  }
  recs.forEach(r => r.seq = cleanSeq(r.seq));
  return recs.filter(r => r.seq.length);
}
/** Returns array of read objects {fileName, recName, seq, qual, traces, ploc, format} */
function parseFile(fileName, arrayBuffer) {
  const u8 = new Uint8Array(arrayBuffer);
  if (u8.length >= 4 && String.fromCharCode(u8[0], u8[1], u8[2], u8[3]) === 'ABIF') {
    const r = parseABIF(arrayBuffer);
    return [{ fileName, recName: fileName, format: 'ab1', ...r }];
  }
  const text = new TextDecoder().decode(u8);
  const recs = parseFasta(text);
  if (!recs.length) throw new Error('Could not read any sequence from this file');
  return recs.map(r => ({ fileName, recName: r.header || fileName, header: r.header, format: 'text',
    seq: r.seq, qual: null, traces: null, ploc: null }));
}

/* ---------- Quality trimming (Mott's modified algorithm, as in phred/Geneious/sangeranalyseR) ---------- */
function mottTrim(qual, cutoff) {
  let best = 0, bestS = 0, bestE = 0, run = 0, runS = 0;
  for (let i = 0; i < qual.length; i++) {
    run += cutoff - Math.pow(10, -qual[i] / 10);
    if (run < 0) { run = 0; runS = i + 1; }
    if (run > best) { best = run; bestS = runS; bestE = i + 1; }
  }
  return [bestS, bestE];
}

/* ---------- Secondary peak detection (mixed / heterozygous positions) ---------- */
function secondaryPeaks(read, ratio) {
  const out = [];
  if (!read.traces || !read.ploc) return out;
  const B = 'ACGT';
  for (let i = 0; i < read.seq.length; i++) {
    const p = read.ploc[i];
    const h = B.split('').map(b => { const t = read.traces[b]; let m = 0; for (let k = p - 2; k <= p + 2; k++) if (k >= 0 && k < t.length && t[k] > m) m = t[k]; return m; });
    const order = [0, 1, 2, 3].sort((a, b) => h[b] - h[a]);
    if (h[order[0]] <= 0) continue;
    const r = h[order[1]] / h[order[0]];
    if (r >= ratio) out.push({ i, primary: B[order[0]], secondary: B[order[1]], ratio: r, code: codeFor(B[order[0]] + B[order[1]]) });
  }
  return out;
}

/* ---------- k-mer based marker + orientation detection ---------- */
const K = 12;
function kmerSet(s) { const set = new Set(); for (let i = 0; i + K <= s.length; i++) { const k = s.substr(i, K); if (/^[ACGT]+$/.test(k)) set.add(k); } return set; }
function sharedKmers(set, s) { let n = 0; for (let i = 0; i + K <= s.length; i++) if (set.has(s.substr(i, K))) n++; return n; }

/* ---------- Pairwise alignment: affine-gap, free end gaps on both sequences (overlap alignment) ---------- */
function scorePair(a, b) {
  if (a === b && !isAmb(a) && a !== 'N') return 2;
  if (a === 'N' || b === 'N') return 0;
  if (compatible(a, b)) return 1;
  return -3;
}
const BIT = { A: 1, C: 2, G: 4, T: 8 };
const MASK = {}; for (const [k, v] of Object.entries(IUPAC)) MASK[k] = [...v].reduce((m, c) => m | BIT[c], 0);
const SC = new Int8Array(256);
for (let x = 0; x < 16; x++) for (let y = 0; y < 16; y++) {
  const pc = v => (v & 1) + ((v >> 1) & 1) + ((v >> 2) & 1) + ((v >> 3) & 1);
  let s; if (x === y && pc(x) === 1) s = 2; else if (x === 15 || y === 15 || x === 0 || y === 0) s = 0; else if (x & y) s = 1; else s = -3;
  SC[x * 16 + y] = s;
}
const encode = s => { const e = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) e[i] = MASK[s[i]] || 15; return e; };
function diagonal(a, b) {
  const pos = new Map();
  for (let i = 0; i + K <= a.length; i++) { const k = a.substr(i, K); pos.set(k, pos.has(k) ? -1 : i); }
  const cnt = new Map(); let best = null, bc = 0;
  for (let j = 0; j + K <= b.length; j++) { const i = pos.get(b.substr(j, K)); if (i === undefined || i < 0) continue; const d = j - i; const c = (cnt.get(d) || 0) + 1; cnt.set(d, c); if (c > bc) { bc = c; best = d; } }
  return bc >= 5 ? best : null;
}
function align(a, b, opt = {}) {
  const ea = encode(a), eb = encode(b);
  const dg = opt.band === 0 ? null : diagonal(a, b), BAND = opt.band || 80;
  const GO = opt.gapOpen ?? 8, GE = opt.gapExt ?? 2;
  const n = a.length, m = b.length, W = m + 1, NEG = -100000000;
  const M = new Int32Array((n + 1) * W), X = new Int32Array((n + 1) * W), Y = new Int32Array((n + 1) * W);
  const tM = new Uint8Array((n + 1) * W), tX = new Uint8Array((n + 1) * W), tY = new Uint8Array((n + 1) * W);
  // M: a_i aligned to b_j ; X: a_i aligned to gap ; Y: b_j aligned to gap. Leading gaps free.
  M.fill(NEG); X.fill(NEG); Y.fill(NEG);
  for (let i = 0; i <= n; i++) {
   const jlo = dg === null ? 0 : Math.max(0, i + dg - BAND), jhi = dg === null ? m : Math.min(m, i + dg + BAND);
   const cols0 = i === 0 ? m : 0; // row 0 and column 0 always initialised
   for (let jj = 0; jj <= cols0; jj++) { const k = i * W + jj; M[k] = (i === 0 && jj === 0) ? 0 : NEG; X[k] = jj === 0 && i > 0 ? 0 : NEG; Y[k] = i === 0 && jj > 0 ? 0 : NEG; }
   for (let j = Math.max(1, jlo); j <= jhi; j++) {
    if (i === 0) continue;
    const k = i * W + j;
    const d = (i - 1) * W + (j - 1), s = SC[ea[i - 1] * 16 + eb[j - 1]];
    let bm = M[d], bt = 0; if (X[d] > bm) { bm = X[d]; bt = 1; } if (Y[d] > bm) { bm = Y[d]; bt = 2; }
    M[k] = bm + s; tM[k] = bt;
    const u = (i - 1) * W + j;
    let xo = M[u] - GO, xe = X[u] - GE, xy = Y[u] - GO;
    if (xo >= xe && xo >= xy) { X[k] = xo; tX[k] = 0; } else if (xe >= xy) { X[k] = xe; tX[k] = 1; } else { X[k] = xy; tX[k] = 2; }
    const l = i * W + (j - 1);
    let yo = M[l] - GO, ye = Y[l] - GE, yx = X[l] - GO;
    if (yo >= ye && yo >= yx) { Y[k] = yo; tY[k] = 0; } else if (ye >= yx) { Y[k] = ye; tY[k] = 2; } else { Y[k] = yx; tY[k] = 1; }
   }
  }
  // Trailing gaps free: best end on last row or last column
  let best = -Infinity, bi = n, bj = m, bs = 0;
  const consider = (i, j) => { const k = i * W + j; for (const [v, s] of [[M[k], 0], [X[k], 1], [Y[k], 2]]) if (v > best) { best = v; bi = i; bj = j; bs = s; } };
  for (let j = 0; j <= m; j++) consider(n, j);
  for (let i = 0; i <= n; i++) consider(i, m);
  const cols = [];
  for (let j = m - 1; j >= bj; j--) cols.push([-1, j]);
  for (let i = n - 1; i >= bi; i--) cols.push([i, -1]);
  // order: trailing overhangs; put a-overhang after b-overhang is arbitrary; reverse later
  let i = bi, j = bj, s = bs;
  while (i > 0 && j > 0) {
    const k = i * W + j;
    if (s === 0) { cols.push([i - 1, j - 1]); s = tM[k]; i--; j--; }
    else if (s === 1) { cols.push([i - 1, -1]); s = tX[k]; i--; }
    else { cols.push([-1, j - 1]); s = tY[k]; j--; }
  }
  while (i > 0) { cols.push([--i, -1]); }
  while (j > 0) { cols.push([-1, --j]); }
  cols.reverse();
  return { score: best, cols };
}

/* ---------- Reference panel ---------- */
/* ---------- PCR primers ---------- */
/** Find a primer (IUPAC-aware) in seq, either strand. Returns best hit {pos, end, strand, mm} or null. */
function findPrimer(seq, primer, maxMM = 3) {
  primer = cleanSeq(primer);
  let best = null;
  for (const [strand, pr] of [['+', primer], ['-', revcomp(primer)]]) {
    const L = pr.length;
    for (let i = 0; i + L <= seq.length; i++) {
      let mm = 0;
      for (let k = 0; k < L && mm <= maxMM; k++) if (!compatible(seq[i + k], pr[k])) mm++;
      if (mm <= maxMM && (!best || mm < best.mm)) best = { pos: i, end: i + L, strand, mm };
    }
  }
  return best;
}

/** All sites a primer could bind (either strand), keeping only hits close to the best match. */
function findPrimerHits(seq, primer, maxMM = 3) {
  primer = cleanSeq(primer);
  const hits = [];
  for (const [strand, pr] of [['+', primer], ['-', revcomp(primer)]]) {
    const L = pr.length;
    for (let i = 0; i + L <= seq.length; i++) {
      let mm = 0;
      for (let k = 0; k < L && mm <= maxMM; k++) if (!compatible(seq[i + k], pr[k])) mm++;
      if (mm <= maxMM) hits.push({ pos: i, end: i + L, strand, mm });
    }
  }
  if (!hits.length) return [];
  const best = Math.min(...hits.map(h => h.mm));
  return hits.filter(h => h.mm <= best + 1);
}
/**
 * Locate both primers. When a primer can bind more than one site (e.g. inside a tandem repeat), use the
 * innermost pair: the shortest product with the primers facing each other. Returns {fwd, rev, extraSites}.
 */
function locatePrimers(seq, pr) {
  const F = findPrimerHits(seq, pr.fwd), R = findPrimerHits(seq, pr.rev);
  let best = null;
  for (const f of F) for (const r of R) {
    if (f.strand === r.strand) continue;
    const [a, b] = f.strand === '+' ? [f, r] : [r, f]; // a must be upstream on the + strand
    if (a.end > b.pos) continue;
    const len = b.end - a.pos;
    if (!best || len < best.len || (len === best.len && f.mm + r.mm < best.mm)) best = { fwd: f, rev: r, len, mm: f.mm + r.mm };
  }
  if (best) return { fwd: best.fwd, rev: best.rev, extraSites: F.length + R.length - 2, product: best.len };
  const pick = H => H.length ? H.reduce((x, y) => y.mm < x.mm ? y : x) : null;
  return { fwd: pick(F), rev: pick(R), extraSites: Math.max(0, F.length - 1) + Math.max(0, R.length - 1), product: null };
}

/**
 * refs: [{marker, name, seq, note}], primers: {MARKER: {fwd, rev, source}}
 * opts.region: 'auto' (between primers when primers are known) | 'primers' | 'full'
 */
function buildPanel(refs, primers = {}, opts = {}) {
  const region = opts.region || 'auto';
  const markers = {};
  refs.forEach((r, i) => (markers[r.marker] = markers[r.marker] || { name: r.marker, refs: [] }).refs.push({ ...r, seq: cleanSeq(r.seq), order: i })); // order = position in the panel file, used to pick display names
  for (const mk of Object.values(markers)) {
    const pr = primers[mk.name] || null;
    mk.primers = pr;
    mk.region = pr && region !== 'full' ? 'primers' : 'full';
    // anchor: prefer a reference that contains both primer sites, oriented so the forward primer reads forward
    let anchor = mk.refs[0];
    if (pr) {
      for (const r of mk.refs) {
        r.primerHits = locatePrimers(r.seq, pr);
      }
      const both = mk.refs.find(r => r.primerHits.fwd && r.primerHits.rev);
      if (both) {
        anchor = both;
        if (both.primerHits.fwd.strand === '-') { both.seq = revcomp(both.seq); both.primerHits = locatePrimers(both.seq, pr); }
      }
      mk.refs = [anchor, ...mk.refs.filter(r => r !== anchor)];
    }
    mk.anchor = anchor.seq;
    mk.kmers = kmerSet(anchor.seq);
    for (const r of mk.refs) {
      if (r !== anchor && sharedKmers(mk.kmers, revcomp(r.seq)) > sharedKmers(mk.kmers, r.seq)) { r.seq = revcomp(r.seq); r.flippedToPanel = true; }
      if (pr && r !== anchor) r.primerHits = locatePrimers(r.seq, pr);
      const al = r === anchor ? { cols: [...r.seq].map((_, i) => [i, i]) } : align(anchor.seq, r.seq);
      r.toAnchor = new Array(r.seq.length).fill(-1);  // ref pos -> anchor pos
      r.fromAnchor = new Array(anchor.seq.length).fill(-1); // anchor pos -> ref pos
      for (const [ai, ri] of al.cols) if (ai >= 0 && ri >= 0) { r.toAnchor[ri] = ai; r.fromAnchor[ai] = ri; }
      const cov = r.fromAnchor.map((v, i) => v >= 0 ? i : -1).filter(v => v >= 0);
      r.anchorSpan = [cov[0], cov[cov.length - 1]];
      // the typed region lies strictly between the primers (primer sequence reflects the primer, not the template)
      r.cut = null;
      if (mk.region === 'primers' && r.primerHits) {
        const hits = [r.primerHits.fwd, r.primerHits.rev].filter(Boolean);
        let lo = 0, hi = r.seq.length - 1;
        if (hits.length === 2) { const [h1, h2] = hits.sort((x, y) => x.pos - y.pos); lo = h1.end; hi = h2.pos - 1; }
        else if (hits.length === 1) { const h = hits[0]; if (h.pos < r.seq.length / 2) lo = h.end; else hi = h.pos - 1; }
        if (hits.length) {
          r.cut = [lo, hi];
          const near = (p, dir) => { while (p >= 0 && p < r.seq.length && r.toAnchor[p] < 0) p += dir; return p >= 0 && p < r.seq.length ? r.toAnchor[p] : null; };
          const a0 = near(lo, 1), a1 = near(hi, -1);
          if (a0 !== null) r.anchorSpan[0] = Math.max(r.anchorSpan[0], a0);
          if (a1 !== null) r.anchorSpan[1] = Math.min(r.anchorSpan[1], a1);
        }
      }
    }
    mk.window = [Math.max(...mk.refs.map(r => r.anchorSpan[0])), Math.min(...mk.refs.map(r => r.anchorSpan[1]))];
    for (const r of mk.refs) {
      let s = mk.window[0], e = mk.window[1];
      while (r.fromAnchor[s] < 0 && s < e) s++;
      while (r.fromAnchor[e] < 0 && e > s) e--;
      r.win = [r.fromAnchor[s], r.fromAnchor[e]]; // inclusive, in ref coords
    }
    // k-mers of the compared region only, so markers cut from the same gene (alt.SAG2, 3′-SAG2) stay distinct
    mk.kmersAll = new Set(); for (const r of mk.refs) for (const k of kmerSet(r.seq.slice(r.win[0], r.win[1] + 1))) mk.kmersAll.add(k);
    // references that are identical over the compared region share one allele; the first in panel order represents it
    const seen = new Map();
    for (const r of [...mk.refs].sort((a, b) => a.order - b.order)) { const k = r.seq.slice(r.win[0], r.win[1] + 1); if (!seen.has(k)) seen.set(k, r); r.rep = seen.get(k); }
    mk.alleles = [...seen.values()];
    mk.pos = p => p - mk.window[0] + 1; // positions are numbered from the start of the compared region
    mk.primerSummary = pr ? { multi: mk.refs.filter(r => r.primerHits && r.primerHits.extraSites > 0).length, products: [...new Set(mk.refs.map(r => r.primerHits && r.primerHits.product).filter(Boolean))].sort((a, b) => a - b), found: mk.refs.filter(r => r.primerHits && r.primerHits.fwd && r.primerHits.rev).length, oneSide: mk.refs.filter(r => r.primerHits && !!r.primerHits.fwd !== !!r.primerHits.rev).length, total: mk.refs.length } : null;
    mk.sites = [];
    for (let p = mk.window[0]; p <= mk.window[1]; p++) {
      const bases = mk.refs.map(r => r.fromAnchor[p] >= 0 ? r.seq[r.fromAnchor[p]] : '-');
      if (new Set(bases).size > 1) mk.sites.push(p);
    }
  }
  return { markers, primers };
}
/** Parse a panel FASTA. Reference records: >MARKER|Allele|note. Primer records (no sequence lines): >MARKER|primers|FWD|REV */
function panelFromFasta(text, defaultMarker) {
  const refs = [], primers = {};
  for (const line of text.split(/\r?\n/)) {
    const m = line.match(/^>\s*([^|]+)\|\s*primers?\s*\|\s*([A-Za-z]+)\s*\|\s*([A-Za-z]+)/i);
    if (m) primers[m[1].trim()] = { fwd: cleanSeq(m[2]), rev: cleanSeq(m[3]), source: 'loaded file' };
  }
  for (const r of parseFasta(text)) {
    const parts = r.header.split('|').map(s => s.trim());
    if (parts[1] && /^primers?$/i.test(parts[1])) continue;
    if (parts.length >= 2) refs.push({ marker: parts[0], name: parts[1], note: parts.slice(2).join(' | '), seq: r.seq });
    else refs.push({ marker: defaultMarker || 'Marker', name: r.header || 'Reference', seq: r.seq });
  }
  return { refs, primers };
}

/* ---------- File name interpretation ---------- */
const DIR_WORDS = [[/forward|foward|fwd/i, 'F'], [/reverse|rev(?![a-z])/i, 'R']];
function markerRegex(name) {
  // punctuation and spaces in marker names are optional in file names: alt.SAG2 = altSAG2 = alt SAG2; 3′-SAG2 = 3SAG2 = 3'SAG2
  let esc = name.split(/[.'’′\-\s_]+/).map(part => part.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')).join("[.'’′\\-\\s_]*");
  esc = esc.replace(/1/g, '[1I]'); // PK1 is often typed PKI
  return new RegExp('(^|[^A-Za-z0-9])(' + esc + ')([FR]?)(?=$|[^A-Za-z0-9])', 'i');
}
function interpretName(fileName, markerNames) {
  let base = fileName.replace(/\.(ab1|abi|scf|seq|fasta|fa|fas|txt)$/i, '');
  let marker = null, dir = null, sample = base;
  const sorted = [...markerNames].sort((a, b) => b.length - a.length);
  for (const mn of sorted) {
    const m = base.match(markerRegex(mn));
    if (m) {
      marker = mn;
      if (m[3]) dir = m[3].toUpperCase();
      sample = base.slice(0, m.index + m[1].length);
      const rest = base.slice(m.index + m[0].length);
      if (!dir) for (const [re, d] of DIR_WORDS) if (re.test(rest)) { dir = d; break; }
      if (!dir) { const t = rest.match(/(?:^|[^A-Za-z0-9])([FR])$/i); if (t) dir = t[1].toUpperCase(); }
      break;
    }
  }
  if (!marker) {
    for (const [re, d] of DIR_WORDS) if (re.test(base)) { dir = d; sample = base.replace(re, ''); break; }
    if (!dir) { const t = base.match(/^(.*?)[\s_-]+([FR])$/i); if (t) { dir = t[2].toUpperCase(); sample = t[1]; } }
  }
  sample = sample.replace(/[\s_+\-.]+$/, '').replace(/^[\s_+\-.]+/, '').replace(/\s+/g, ' ');
  const isRef = /(^|[^a-z])(ref|reference)([^a-z]|$)/i.test(base);
  return { sample: sample || base, marker, dir, isRef };
}
const sampleKey = s => s.toLowerCase().replace(/[^a-z0-9]/g, '');

/* ---------- Analysis ---------- */
const DEFAULTS = { mottCutoff: 0.05, peakRatio: 0.33, minKmers: 15, minTrimmedLen: 100 };

function assignMarker(read, panel) {
  let best = null;
  for (const mk of Object.values(panel.markers)) {
    const f = sharedKmers(mk.kmersAll, read.seq), r = sharedKmers(mk.kmersAll, revcomp(read.seq));
    const score = Math.max(f, r);
    if (!best || score > best.score) best = { marker: mk.name, score, reversed: r > f };
  }
  return best;
}

/** Prepare a read: orient to panel, trim, flag mixed peaks. Returns working read object. */
function prepareRead(raw, mk, opts) {
  const o = { ...raw, notes: [] };
  const L = raw.seq.length;
  // trimming on original orientation
  let [ts, te] = [0, L];
  if (raw.qual) { [ts, te] = mottTrim(raw.qual, opts.mottCutoff); o.trimmed = [ts, te]; }
  // secondary peaks -> IUPAC
  let seq = raw.seq.split('');
  o.peaks = [];
  if (raw.traces) {
    o.peaks = secondaryPeaks(raw, opts.peakRatio).filter(p => p.i >= ts && p.i < te);
    for (const p of o.peaks) seq[p.i] = p.code;
  }
  // working sequence + index back into original read
  let w = seq.slice(ts, te), idx = [...Array(te - ts).keys()].map(k => k + ts);
  const f = sharedKmers(mk.kmersAll, w.join('')), r = sharedKmers(mk.kmersAll, revcomp(w.join('')));
  o.reversed = r > f;
  if (o.reversed) { w = revcomp(w.join('')).split(''); idx = idx.reverse(); }
  o.work = w.join(''); o.idx = idx; o.trimLen = w.length; o.rawLen = L;
  return o;
}

function buildConsensus(reads) {
  if (reads.length === 1) {
    const r = reads[0];
    return { seq: r.work, cols: [...r.work].map((b, k) => ({ base: b, cov: 1, status: isAmb(b) ? 'ambiguous' : 'single', src: [[0, r.idx[k]]] })) };
  }
  const [a, b] = reads;
  const al = align(a.work, b.work);
  const both = al.cols.map((c, k) => c[0] >= 0 && c[1] >= 0 ? k : -1).filter(k => k >= 0);
  const ovS = both.length ? both[0] : Infinity, ovE = both.length ? both[both.length - 1] : -1;
  const cols = [];
  al.cols.forEach(([i, j], k) => {
    const x = i >= 0 ? a.work[i] : null, y = j >= 0 ? b.work[j] : null;
    const src = []; if (i >= 0) src.push([0, a.idx[i]]); if (j >= 0) src.push([1, b.idx[j]]);
    if (x && y) {
      if (x === y) cols.push({ base: x, cov: 2, status: x === 'N' ? 'unresolved' : (isAmb(x) ? 'mixed' : 'agree'), src });
      else if (compatible(x, y) && (setOf(x).length === 1 || setOf(y).length === 1)) {
        const clean = setOf(x).length === 1 ? x : y;
        cols.push({ base: clean, cov: 2, status: 'resolved', other: clean === x ? y : x, src });
      } else cols.push({ base: 'N', cov: 2, status: 'conflict', pair: x + '/' + y, src });
    } else {
      const base = x || y;
      if (k > ovS && k < ovE) cols.push({ base, cov: 1, status: 'indel-conflict', src });
      else cols.push({ base, cov: 1, status: isAmb(base) ? 'ambiguous' : 'single', src });
    }
  });
  // indel-conflict columns: one read has a base the other lacks (common base-caller error). Kept, but never counted as evidence.
  return { seq: cols.map(c => c.base).join(''), cols, overlap: al };
}

function compareToRef(cons, ref, mk) {
  const al = align(ref.seq, cons.seq);
  const [ws, we] = ref.win;
  const diffs = []; let uncovered = 0, unresolved = [], ambCompat = [], singleCov = 0, covered = 0;
  const uncoveredSet = new Set(); const refToCons = new Array(ref.seq.length).fill(-1); let lastCons = -1;
  const firstC = al.cols.findIndex(c => c[1] >= 0), lastC = al.cols.length - 1 - [...al.cols].reverse().findIndex(c => c[1] >= 0);
  let lastRefPos = -1;
  al.cols.forEach(([ri, ci], k) => {
    if (ri >= 0) lastRefPos = ri;
    if (ci >= 0) lastCons = ci;
    if (ri >= 0 && ci >= 0) refToCons[ri] = ci;
    const inWin = ri >= 0 ? (ri >= ws && ri <= we) : (lastRefPos >= ws && lastRefPos < we);
    if (!inWin) return;
    const anchorPos = ri >= 0 ? ref.toAnchor[ri] : (lastRefPos >= 0 ? ref.toAnchor[lastRefPos] : -1);
    const within = k >= firstC && k <= lastC;
    if (ri >= 0 && ci < 0) {
      if (!within) { uncovered++; uncoveredSet.add(ri); return; }
      diffs.push({ type: 'deletion', refPos: ri, anchorPos, consPos: Math.max(0, lastCons), ref: ref.seq[ri], cons: '-', cov: 2 }); return;
    }
    if (ri < 0 && ci >= 0) { if (cons.cols[ci].status === 'indel-conflict') { unresolved.push({ refPos: lastRefPos, anchorPos, consPos: ci, ref: '-', detail: 'extra ' + cons.seq[ci] + ' in one read' }); return; } diffs.push({ type: 'insertion', refPos: lastRefPos, anchorPos, ref: '-', cons: cons.seq[ci], consPos: ci, cov: cons.cols[ci].cov }); return; }
    const c = cons.seq[ci], r = ref.seq[ri], col = cons.cols[ci];
    covered++; if (col.cov < 2) singleCov++;
    if (c === 'N') { unresolved.push({ refPos: ri, anchorPos, consPos: ci, ref: r, detail: col.pair || 'N' }); return; }
    if (c === r) return;
    if (col.status === 'indel-conflict') { unresolved.push({ refPos: ri, anchorPos, consPos: ci, ref: r, detail: c + ' in one read, missing in other' }); return; }
    if (isAmb(c) && setOf(c).includes(r)) { ambCompat.push({ refPos: ri, anchorPos, consPos: ci, ref: r, cons: c }); return; }
    diffs.push({ type: 'substitution', refPos: ri, anchorPos, consPos: ci, ref: r, cons: c, cov: col.cov });
  });
  // collapse adjacent indel columns into events
  let events = 0, prev = null;
  for (const d of diffs) { if (d.type === 'substitution' || !prev || prev.type !== d.type || Math.abs(d.refPos - prev.refPos) > 1) events++; prev = d; }
  return { ref, diffs, events, uncovered, uncoveredSet, refToCons, unresolved, ambCompat, singleCov, covered, winLen: we - ws + 1 };
}

function callAllele(cons, mk) {
  const comps = mk.refs.map(r => compareToRef(cons, r, mk));
  comps.sort((a, b) => a.events - b.events || a.uncovered - b.uncovered || a.ref.order - b.ref.order);
  const exact = comps.filter(c => c.events === 0);
  const flags = [];
  let call, status, sameAs = [];
  // which anchor positions discriminate among the candidate references?
  const discrim = (refs) => {
    const s = new Set();
    for (let p = mk.window[0]; p <= mk.window[1]; p++) {
      const b = refs.map(r => r.fromAnchor[p] >= 0 ? r.seq[r.fromAnchor[p]] : '-');
      if (new Set(b).size > 1) s.add(p);
    }
    return s;
  };
  const allSites = discrim(mk.refs);
  const best = comps[0];
  // problem positions (unresolved / ambiguous) at sites that separate references
  const unresAtSites = best.unresolved.filter(u => allSites.has(u.anchorPos));
  const ambAtSites = best.ambCompat.filter(u => allSites.has(u.anchorPos));
  if (exact.length === 1) { call = exact[0].ref.name; status = 'ok'; }
  else if (exact.length > 1) {
    const d = discrim(exact.map(c => c.ref));
    call = exact.map(c => c.ref.name).join(' / ');
    if (d.size === 0) {
      // references with the same allele over the compared region: report the first (panel order) and list the rest
      call = exact[0].ref.name; status = 'ok';
      sameAs = [...new Set(exact.slice(1).map(c => c.ref.name))].filter(n => n !== call);
      if (sameAs.length) flags.push({ level: 'info', text: `${call} is identical over the compared region of ${mk.name} to: ${sameAs.join(', ')}.` });
    }
    else { status = 'review'; flags.push({ level: 'warn', text: `Matches ${call} equally; the ${d.size} position(s) that tell them apart are not resolved in this sample.` }); }
  } else {
    call = 'No exact match'; status = 'novel';
    flags.push({ level: 'bad', text: `Closest reference is ${best.ref.name} with ${best.events} difference(s). Check the chromatograms, then BLAST the consensus.` });
  }
  if (best.uncovered > 0.3 * best.winLen && mk.region === 'full') {
    flags.unshift({ level: 'bad', text: `The ${mk.name} references are much longer than this sample (it covers ${Math.round(100 * (1 - best.uncovered / best.winLen))}% of them), so differences at the read ends may be spurious. ` +
      (mk.primers ? `Set "Region compared" to "Between the PCR primers" in Settings.` : `Add the PCR primers for ${mk.name} in the Reference panel, or use references trimmed to the amplicon.`) });
    if (status === 'ok') status = 'review';
  } else if (best.uncovered > 0) {
    const missingSites = [...allSites].filter(p => { const rp = best.ref.fromAnchor[p]; return rp >= 0 && !isCoveredRefPos(cons, best, rp); });
    flags.push({ level: missingSites.length ? 'warn' : 'info', text: `Consensus does not cover ${best.uncovered} of ${best.winLen} bases of the comparison region` + (missingSites.length ? `, including ${missingSites.length} position(s) that distinguish references.` : '.') });
    if (missingSites.length && status === 'ok') status = 'review';
  }
  if (unresAtSites.length) { flags.push({ level: 'warn', text: `${unresAtSites.length} position(s) that distinguish references have forward/reverse disagreement: ${unresAtSites.map(u => 'pos ' + mk.pos(u.anchorPos) + ' (' + u.detail + ')').join(', ')}.` }); if (status === 'ok') status = 'review'; }
  const refBase = (r, p) => r.fromAnchor[p] >= 0 ? r.seq[r.fromAnchor[p]] : '-';
  const mixed = ambAtSites.filter(a => { const bases = new Set(mk.refs.map(r => refBase(r, a.anchorPos))); return [...setOf(a.cons)].filter(x => bases.has(x)).length >= 2; });
  if (mixed.length) {
    // find the pair of references whose combination explains every distinguishing position
    const consAt = new Map(); best.ambCompat.forEach(a => consAt.set(a.anchorPos, a.cons));
    const unresSet = new Set(best.unresolved.map(u => u.anchorPos));
    const pairs = [];
    const reps = mk.alleles || mk.refs;
    for (let x = 0; x < reps.length; x++) for (let y = x + 1; y < reps.length; y++) {
      const A = reps[x], B = reps[y]; let ok = true;
      for (const p of allSites) {
        if (unresSet.has(p)) continue;
        const want = new Set([refBase(A, p), refBase(B, p)]);
        const c = consAt.get(p) || refBase(best.ref, p); // non-ambiguous positions equal the best ref
        const have = new Set(setOf(c) || c);
        if (want.size !== have.size || [...want].some(b => !have.has(b))) { ok = false; break; }
      }
      if (ok) pairs.push([A.name, B.name]);
    }
    const pairTxt = pairs.length ? pairs.map(p => p.join(' + ')).join(' or ') : null;
    flags.push({ level: 'warn', text: `Double peaks at ${mixed.length} position(s) that distinguish references (${mixed.map(a => 'pos ' + mk.pos(a.anchorPos) + ' ' + a.cons).join(', ')}). ` +
      (pairTxt ? `The pattern fits a mixture of ${pairTxt}.` : 'The pattern does not fit a simple mixture of two references.') + ' Possible mixed infection; review the chromatograms.' });
    if (pairTxt) { call = 'Mixed? ' + pairTxt; for (let i = flags.length - 1; i >= 0; i--) if (flags[i].text.startsWith('Matches ')) flags.splice(i, 1); }
    status = 'review';
  } else if (best.ambCompat.length) flags.push({ level: 'info', text: `${best.ambCompat.length} ambiguous base(s) outside reference-distinguishing positions; they do not affect the call.` });
  if (best.unresolved.length - unresAtSites.length > 0) flags.push({ level: 'info', text: `${best.unresolved.length - unresAtSites.length} forward/reverse disagreement(s) at positions that do not affect the call (excluded).` });
  const singleDiffs = best.diffs.filter(d => d.cov === 1);
  if (status === 'novel' && singleDiffs.length && !cons.assembled) {
    if (singleDiffs.length === best.diffs.length) {
      call = best.ref.name; status = 'review';
      flags.unshift({ level: 'warn', text: `Every difference from ${best.ref.name} is seen in only one read, which often means a base-calling error at the end of a read. Check the chromatogram before accepting ${best.ref.name}.` });
      const i = flags.findIndex(f => f.level === 'bad'); if (i >= 0) flags.splice(i, 1);
    } else flags.push({ level: 'warn', text: `${singleDiffs.length} of the ${best.diffs.length} differences are seen in only one read.` });
  }
  return { call, status, flags, sameAs, comps, best, sites: [...allSites].sort((a, b) => a - b) };
}
function isCoveredRefPos(cons, comp, rp) { return !comp.uncoveredSet || !comp.uncoveredSet.has(rp); }

/**
 * Main entry. files: [{fileName, reads:[parsed read objects]}]; panel from buildPanel.
 * overrides: {fileName: {sample, marker}} user edits.
 */
function planGroups(readsIn, panel, opts = {}, overrides = {}) {
  opts = { ...DEFAULTS, ...opts };
  const markerNames = Object.keys(panel.markers);
  const reads = [], problems = [];
  for (const raw of readsIn) {
    const nm = interpretName(raw.fileName, markerNames);
    const ov = overrides[raw.fileName] || {};
    const asg = assignMarker(raw, panel);
    const r = { ...raw, parsed: nm, sample: ov.sample ?? nm.sample, dir: nm.dir, isRef: ov.isRef ?? nm.isRef,
      marker: ov.marker || (asg && asg.score >= opts.minKmers ? asg.marker : null) || nm.marker, kmerScore: asg ? asg.score : 0, notes: [] };
    if (nm.marker && asg && asg.score >= opts.minKmers && nm.marker !== asg.marker && !ov.marker) r.notes.push(`File name says ${nm.marker} but the sequence matches ${asg.marker}; using ${asg.marker}.`);
    if (!asg || asg.score < opts.minKmers) r.notes.push('Sequence does not match any marker in the reference panel.');
    reads.push(r);
  }
  // group sample reads
  const groups = new Map();
  for (const r of reads) {
    if (r.isRef) continue;
    if (!r.marker || !panel.markers[r.marker]) { problems.push({ file: r.fileName, text: r.notes.join(' ') || 'Unknown marker' }); continue; }
    const key = sampleKey(r.sample) + '||' + r.marker;
    if (!groups.has(key)) groups.set(key, { sample: r.sample, marker: r.marker, reads: [] });
    groups.get(key).reads.push(r);
  }
  return { reads, problems, groups: [...groups.values()], opts };
}
function analyzeGroup(g, panel, opts) {
  opts = { ...DEFAULTS, ...opts };
  {
    const mk = panel.markers[g.marker];
    const res = { sample: g.sample, marker: g.marker, files: g.reads.map(r => r.fileName), flags: [], notes: [] };
    let prepared = g.reads.map(r => prepareRead(r, mk, opts));
    for (const p of prepared) {
      if (p.trimmed) p.notes.push(`Quality trimmed to ${p.trimmed[0] + 1}–${p.trimmed[1]} of ${p.rawLen} bases.`);
      if (p.trimLen < opts.minTrimmedLen) res.flags.push({ level: 'bad', text: `${p.fileName}: only ${p.trimLen} good-quality bases. Read probably failed.` });
      res.notes.push(...p.notes.map(n => p.fileName + ': ' + n));
    }
    const good = prepared.filter(p => p.trimLen >= opts.minTrimmedLen);
    // prefer one forward + one reverse
    good.sort((a, b) => (a.dir === 'F' ? 0 : 1) - (b.dir === 'F' ? 0 : 1));
    if (good.length > 2) res.flags.push({ level: 'info', text: `${good.length} reads found; the first forward and reverse reads were used.` });
    let used = good.slice(0, 2);
    if (good.length > 2) { const f = good.find(r => r.dir === 'F'), rv = good.find(r => r.dir === 'R'); if (f && rv) used = [f, rv]; }
    res.used = used;
    if (!used.length) { res.call = 'Failed'; res.status = 'bad'; return res; }
    const assembled = used.length === 1 && !used[0].dir;
    if (assembled) res.flags.push({ level: 'info', text: 'Single sequence with no forward/reverse label, treated as a finished consensus (for example from nanopore or an earlier assembly).' });
    else if (used.length === 1) res.flags.push({ level: 'warn', text: 'Only one usable read; the call rests on single-read coverage.' });
    if (used.length === 2 && used[0].reversed === used[1].reversed) res.flags.push({ level: 'warn', text: 'Both reads point the same direction. Check that one forward and one reverse read were supplied.' });
    const cons = buildConsensus(used);
    cons.assembled = assembled;
    res.consensus = cons;
    const conflicts = cons.cols.filter(c => c.status === 'conflict').length;
    res.stats = { length: cons.seq.length, twoRead: cons.cols.filter(c => c.cov === 2).length, conflicts, indelConflicts: cons.cols.filter(c => c.status === 'indel-conflict').length,
      mixed: cons.cols.filter(c => c.status === 'mixed' || c.status === 'ambiguous').length };
    const c = callAllele(cons, mk);
    Object.assign(res, { call: c.call, status: c.status, sameAs: c.sameAs, comps: c.comps, best: c.best, sites: c.sites });
    res.flags.push(...c.flags);
    return res;
  }
}
const sortResults = rs => rs.sort((a, b) => a.sample.localeCompare(b.sample, undefined, { numeric: true }) || a.marker.localeCompare(b.marker));
function analyze(readsIn, panel, opts = {}, overrides = {}) {
  const plan = planGroups(readsIn, panel, opts, overrides);
  const results = sortResults(plan.groups.map(g => analyzeGroup(g, panel, plan.opts)));
  return { reads: plan.reads, results, problems: plan.problems };
}

const API = { IUPAC, revcomp, cleanSeq, parseABIF, parseFasta, parseFile, mottTrim, secondaryPeaks, align, buildPanel, panelFromFasta,
  interpretName, sampleKey, findPrimer, locatePrimers, analyze, planGroups, analyzeGroup, sortResults, DEFAULTS, encode, isAmb, setOf };
if (typeof module !== 'undefined' && module.exports) module.exports = API; else root.MLST = API;
})(typeof window !== 'undefined' ? window : globalThis);
