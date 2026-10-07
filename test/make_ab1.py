"""Write synthetic ABIF (.ab1) chromatograms from FASTA basecalls, for testing the parser/trimmer/peak detector."""
import struct, math, random, sys, os, glob
from Bio import SeqIO
IU={'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'CG','W':'AT','K':'GT','M':'AC','N':'ACGT'}
def write_ab1(path, seq, qual, sample, secondary_ratio=0.6, seed=0):
    rnd=random.Random(seed); sp=12; L=len(seq)*sp+20
    ch={b:[0.0]*L for b in 'GATC'}; ploc=[]
    for i,b in enumerate(seq):
        c=10+i*sp; ploc.append(c)
        bases=IU.get(b,'ACGT'); 
        q=qual[i]; h=900+rnd.randint(-150,150)
        hs=[h]+[h*secondary_ratio]*(len(bases)-1) if len(bases)<=2 else [h*0.3]*len(bases)
        for bb,hh in zip(bases,hs):
            for x in range(max(0,c-10),min(L,c+11)): ch[bb][x]+=hh*math.exp(-((x-c)**2)/(2*2.6**2))
        if q<20:  # noisy background in poor regions
            for bb in 'ACGT':
                if bb not in bases:
                    hh=h*rnd.uniform(0.1,0.45)
                    for x in range(max(0,c-8),min(L,c+9)): ch[bb][x]+=hh*math.exp(-((x-c-rnd.randint(-3,3))**2)/(2*3**2))
    basecalls=''.join(b if b in 'ACGT' else IU.get(b,'ACGT')[0] for b in seq)  # instrument calls primary base
    entries=[]; data=b''; base_off=128
    def add(name,num,eltype,elsize,arr_bytes,numel):
        nonlocal data
        entries.append((name,num,eltype,elsize,numel,len(arr_bytes),arr_bytes))
    order='GATC'
    for k,b in enumerate(order):
        arr=[min(32000,int(v)) for v in ch[b]]; add('DATA',9+k,4,2,struct.pack('>%dh'%len(arr),*arr),len(arr))
    add('FWO_',1,2,1,order.encode(),4)
    add('PBAS',2,2,1,basecalls.encode(),len(basecalls))
    add('PCON',2,2,1,bytes(qual),len(qual))
    add('PLOC',2,4,2,struct.pack('>%dh'%len(ploc),*ploc),len(ploc))
    sm=sample.encode()[:255]; add('SMPL',1,18,1,bytes([len(sm)])+sm,len(sm)+1)
    # lay out data
    off=base_off; dirents=b''; blob=b''
    for name,num,et,es,ne,ds,by in entries:
        if ds<=4: dirents+=struct.pack('>4sihhiii',name.encode(),num,et,es,ne,ds,0)[:20]+by.ljust(4,b'\0')+struct.pack('>i',0)
        else: dirents+=struct.pack('>4sihhiiii',name.encode(),num,et,es,ne,ds,off+len(blob),0); blob+=by
    dir_off=off+len(blob)
    header=b'ABIF'+struct.pack('>h',101)+struct.pack('>4sihhiiii',b'tdir',1,1023,28,len(entries),len(dirents),dir_off,0)
    header=header.ljust(128,b'\0')
    open(path,'wb').write(header+blob+dirents)

def qual_profile(n, rnd, start_bad=25, end_bad=40):
    q=[]
    for i in range(n):
        if i<start_bad: q.append(rnd.randint(4,15) if i<start_bad-8 else rnd.randint(15,30))
        elif i>=n-end_bad: q.append(rnd.randint(3,14))
        else: q.append(rnd.randint(35,55))
    return q

if __name__=='__main__':
    src,dst=sys.argv[1],sys.argv[2]; os.makedirs(dst,exist_ok=True)
    for i,f in enumerate(sorted(glob.glob(src+'/*.ab1'))):
        if 'Ref' in f: continue
        rec=next(SeqIO.parse(f,'fasta')); s=str(rec.seq).upper(); rnd=random.Random(i)
        write_ab1(os.path.join(dst,os.path.basename(f)), s, qual_profile(len(s),rnd), rec.id, seed=i)
