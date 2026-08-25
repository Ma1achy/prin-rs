"""Parser for prin-rs .prnq quadtree dumps.

Layout, determined by inspection: an ASCII header ending with a  line, then a
12-byte binary prefix (u32 n_quads, u32 pad, u32 n_fields), then n_quads * n_fields
little-endian f64. Integer-valued fields (level, decision, is_leaf, n_hot_*) are stored
as f64 like everything else -- do not reinterpret them as ints.
"""
import numpy as np, struct, glob, os
def load(path):
    d=open(path,'rb').read()
    i=d.find(b'fields='); j=d.find(b'\n',i)
    f=[x.strip() for x in d[i+7:j].decode().strip().split(',')]
    hdr=d[:i].decode(errors='replace')
    body=d[j+1:]
    n,_,nf=struct.unpack('<3I',body[:12])
    a=np.frombuffer(body[12:12+n*nf*8],dtype='<f8').reshape(n,nf)
    return {k:a[:,x] for x,k in enumerate(f)}, hdr
def hdrval(h,k):
    for tok in h.replace('\n',' ').split():
        if tok.startswith(k+'='): return tok[len(k)+1:]
    return None
