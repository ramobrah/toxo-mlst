# Independent Python reference implementation (Biopython) for checking the JS app
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import glob,re,collections
recs={f:str(next(SeqIO.parse(f,'fasta')).seq).upper() for f in glob.glob('*.ab1')}
rc=lambda s:str(Seq(s).reverse_complement())
al=PairwiseAligner(); al.mode='global'; al.match_score=2; al.mismatch_score=-3; al.open_gap_score=-6; al.extend_gap_score=-2
al.target_end_gap_score=0; al.query_end_gap_score=0
IUPAC={'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'CG','W':'AT','K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}
anchor=recs['TypeI_PK1_Ref.ab1']
def orient(s):
    lo=PairwiseAligner(); lo.mode='local'
    return s if lo.score(anchor,s)>=lo.score(anchor,rc(s)) else rc(s)
def proj(seq):
    """project seq onto anchor coords: list of anchor-position -> base ('-' gap, None uncovered); ignore insertions but record them"""
    a=al.align(anchor,seq)[0]
    out=[None]*len(anchor); ins=[]
    for (t0,t1),(q0,q1) in zip(*a.aligned):
        for i in range(t1-t0): out[t0+i]=seq[q0+i]
    # gaps inside aligned block
    c=a.coordinates
    for k in range(c.shape[1]-1):
        t0,t1,q0,q1=c[0][k],c[0][k+1],c[1][k],c[1][k+1]
        if t1>t0 and q1==q0 and 0<q0<len(seq): 
            for i in range(t0,t1): out[i]='-'
        if q1>q0 and t1==t0 and 0<t0<len(anchor): ins.append((t0,seq[q0:q1]))
    return out,ins
refs={}
for f,s in recs.items():
    if 'Ref' in f: refs[f]=proj(orient(s))
for f,(p,ins) in refs.items(): print(f, 'ins',ins, 'covered',sum(x is not None for x in p))
# informative sites among refs
names=list(refs)
win=[i for i in range(len(anchor)) if all(refs[n][0][i] is not None for n in names)]
inf=[i for i in win if len(set(refs[n][0][i] for n in names))>1]
print('window',win[0],win[-1],len(win),'informative',len(inf))
for n in names: print(f"{n:30s}",''.join(refs[n][0][i] for i in inf))
# samples
samp=collections.defaultdict(dict)
for f,s in recs.items():
    if 'Ref' in f: continue
    m=re.match(r'(.*?)[ _]*(?:_?PK1|PKI)',f); name=re.sub(r'[ _]+$','',f.split('PK1')[0]).strip(' _')
    d='F' if re.search(r'(Foward|Forward|PK1F|_F\.ab1|PK1_F)',f) else 'R'
    samp[name][d]=proj(orient(s))
for name,d in samp.items():
    F,_=d['F'];R,_=d['R']
    cons=[]
    for i in range(len(anchor)):
        a,b=F[i],R[i]
        if a is None: cons.append(b); continue
        if b is None: cons.append(a); continue
        if a==b: cons.append(a); continue
        sa,sb=set(IUPAC.get(a,'')),set(IUPAC.get(b,''))
        if sa and sb and (sa<=sb or sb<=sa): cons.append(a if len(sa)<len(sb) else b)  # clean base wins
        else: cons.append('?')
    res={}
    for n in names:
        rp=refs[n][0]; diff=[]
        for i in win:
            c,r=cons[i],rp[i]
            if c is None: diff.append((i,'nocov')); continue
            if c==r: continue
            if c in IUPAC and len(IUPAC[c])>1 and r in IUPAC.get(c,''): diff.append((i,'amb'+c)); continue
            diff.append((i,f'{r}>{c}'))
        res[n]=diff
    print(f"\n{name}: sample at informative sites:",''.join(str(cons[i]) for i in inf))
    for n in names: print(f"   {n:30s} confident={[x for x in res[n] if '?' not in x[1]]} unresolved={[x[0] for x in res[n] if '?' in x[1]]}")
