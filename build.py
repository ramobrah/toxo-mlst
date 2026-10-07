import json,glob,os,sys
NOEX='--no-examples' in sys.argv
root=os.path.dirname(os.path.abspath(__file__))
s=open(root+'/src/app.html').read()
core=open(root+'/src/core.js').read()
panel=json.load(open(root+'/src/panel.json'))
ex=[]
for f in ([] if NOEX else sorted(glob.glob(root+'/examples/*.ab1'))):
    if 'Ref' in f: continue
    ex.append({'name':os.path.basename(f),'text':open(f).read()})
s=s.replace('/*CORE*/',core).replace('/*PANEL*/',json.dumps(panel)).replace('/*EXAMPLE*/',json.dumps(ex))

open(root+'/index.html','w').write(s)
print(len(s))
