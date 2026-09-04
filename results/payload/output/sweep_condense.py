#!/usr/bin/env python3
"""Condensed README tables from the Phase 3 sweep logs. usage: sweep_condense.py <scratchpad dir>"""
import re, sys, os
D = sys.argv[1]
HDR = re.compile(r'^=== (live|march) (\S+) (.*?)\s+(\d\d:\d\d:\d\d)$')
DESC = re.compile(r'descent: (\d+) quads \((\d+) leaves, depth (\d+)\) in ([\d.]+)s, ([\d.e+-]+) substeps(?:, catch-up ([\d.e+-]+) substeps \(([\d.]+)% of the total\))?, stop \[(.*)\]')
MEM = re.compile(r'memory: (\d+) quads ever computed \(mem_all\), resident peak (\d+) and final (\d+) \(mem_resident\), (\d+) children merged')
PAY = re.compile(r'payload/event_class/(indicator|resolvable)/eps=([\de.-]+): (?:final )?error ([\d.]+); dp needs B=(\d+) \(ratio ([\d.]+)x\); uniform needs B=(\d+) \(ratio ([\d.]+)x\); sea_fraction ([\d.]+)')
R = {}
def load(name):
    p = os.path.join(D, name)
    if not os.path.exists(p): return
    cur = None
    for line in open(p):
        m = HDR.match(line.strip())
        if m:
            cur = dict(kind=m.group(1), target=m.group(2), variant=m.group(3), log=name); R[(m.group(1), m.group(2), name, m.group(3))] = cur; continue
        if cur is None: continue
        m = DESC.search(line)
        if m:
            cur.update(quads=int(m.group(1)), leaves=int(m.group(2)), substeps=float(m.group(5)))
            if m.group(6): cur['catchup_pct'] = float(m.group(7))
            cur['stops'] = {k: int(v) for k, v in (kv.split(':') for kv in m.group(8).split())}; continue
        m = MEM.search(line)
        if m: cur.update(mem_all=int(m.group(1)), res_peak=int(m.group(2)), res_final=int(m.group(3)), merged=int(m.group(4))); continue
        m = PAY.search(line)
        if m: cur[m.group(1)] = dict(err=float(m.group(3)), dp_r=float(m.group(5)), uni_r=float(m.group(7)), sea=float(m.group(8)))
for f in ['phase3_alo.txt', 'phase3_noagree.txt', 'phase3_alo0.txt', 'phase3_march2.txt', 'phase3_ladder.txt', 'phase3_nodim.txt', 'phase3_fine.txt', 'phase3_fine2.txt', 'phase3_march3.txt']: load(f)
T = ['near-field', 'deep_interior', 'preset_prho', 'preset_shape', 'config_stability', 'preset_shape_h1']
def g(kind, t, log, var):
    r = R.get((kind, t, log, var)); return r if r and 'indicator' in r else None
def fl(r): return r['stops'].get('floor', 0)
def st(r): return r['stops'].get('stationary', 0)
def e(r, k='indicator'): return f"{r[k]['err']:.4f}"
print("### The area floor, static, all six targets\n")
print("| target | sea | full depth: quads | vs ref | floor 0.2: quads | vs ref | resolvable left | vs optimum | vs uniform | floors | arm off: quads | vs ref | resolvable left | vs uniform | floors |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in T:
    full = g('live', t, 'phase3_alo0.txt', 'stationary=0 alpha_lo=0 (guarded)') or g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0')
    on = g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0.2'); off = g('live', t, 'phase3_noagree.txt', 'stationary=0 alpha_lo=0.2 agreement=0')
    if not (full and on and off): print(f"| {t} | (incomplete) |"); continue
    print(f"| {t} | {on['indicator']['sea']:.4f} | {full['quads']} | {e(full)} | {on['quads']} | {e(on)} | {e(on,'resolvable')} | {on['indicator']['dp_r']:.2f}x | {on['indicator']['uni_r']:.2f}x | {fl(on)} | {off['quads']} | {e(off)} | {e(off,'resolvable')} | {off['indicator']['uni_r']:.2f}x | {fl(off)} |")
print("\n### The stationarity stop, static\n")
print("| target | full depth: fires | vs ref off -> on | floor 0.2: fires | vs ref off -> on |")
print("|---|---|---|---|---|")
for t in T:
    f0 = g('live', t, 'phase3_alo0.txt', 'stationary=0 alpha_lo=0 (guarded)') or g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0')
    f1 = g('live', t, 'phase3_alo0.txt', 'stationary=1 alpha_lo=0 (guarded)') or g('live', t, 'phase3_alo.txt', 'stationary=1 alpha_lo=0')
    a0 = g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0.2'); a1 = g('live', t, 'phase3_alo.txt', 'stationary=1 alpha_lo=0.2')
    if not (f0 and f1 and a0 and a1): print(f"| {t} | (incomplete) |"); continue
    print(f"| {t} | {st(f1)} | {e(f0)} -> {e(f1)} | {st(a1)} | {e(a0)} -> {e(a1)} |")
print("\n### The live march against the static tree, floor 0.2\n")
print("| target | static: quads | vs ref | march: computed | resident final | merged | vs ref | resolvable left | vs uniform | catch-up | pin | arm off march: computed | resident | merged | vs ref |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in T:
    s = g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0.2')
    m = g('march', t, 'phase3_march2.txt', 'expiry-on-structure'); pin = 'expiry'
    if not m: m = g('march', t, 'phase3_alo.txt', 'default (alpha_lo 0.2, merge on)'); pin = 'pre-expiry'
    o = g('march', t, 'phase3_noagree.txt', 'agreement=0')
    if not (s and m and o): print(f"| {t} | (incomplete) |"); continue
    print(f"| {t} | {s['quads']} | {e(s)} | {m['mem_all']} | {m['res_final']} | {m['merged']} | {e(m)} | {e(m,'resolvable')} | {m['indicator']['uni_r']:.2f}x | {m['catchup_pct']:.0f}% | {pin} | {o['mem_all']} | {o['res_final']} | {o['merged']} | {e(o)} |")
print("\n### The alpha_lo ladder, static, arm on\n")
print("| target | alpha_lo | quads | vs ref | resolvable left | vs optimum | vs uniform | floors | quads saved |")
print("|---|---|---|---|---|---|---|---|---|")
RUNGS = [('0', 'phase3_alo0.txt', 'stationary=0 alpha_lo=0 (guarded)'), ('0.001', None, 'stationary=0 alpha_lo=0.001'),
         ('0.005', None, 'stationary=0 alpha_lo=0.005'), ('0.02', None, 'stationary=0 alpha_lo=0.02'),
         ('0.05', 'phase3_ladder.txt', 'stationary=0 alpha_lo=0.05'), ('0.1', 'phase3_ladder.txt', 'stationary=0 alpha_lo=0.1'),
         ('0.2', 'phase3_alo.txt', 'stationary=0 alpha_lo=0.2'), ('0.3', 'phase3_ladder.txt', 'stationary=0 alpha_lo=0.3')]
for t in T:
    full = g('live', t, 'phase3_alo0.txt', 'stationary=0 alpha_lo=0 (guarded)') or g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0')
    for a, log, var in RUNGS:
        r = None
        for lg in ([log] if log else ['phase3_fine.txt', 'phase3_fine2.txt']):
            r = r or g('live', t, lg, var)
        if a == '0' and not r: r = g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0')
        if not r: continue
        saved = f"{100*(1-r['quads']/full['quads']):.0f}%" if full else ""
        print(f"| {t} | {a} | {r['quads']} | {e(r)} | {e(r,'resolvable')} | {r['indicator']['dp_r']:.2f}x | {r['indicator']['uni_r']:.2f}x | {fl(r)} | {saved} |")
print("\n### The live march at alpha_lo 0.005 against 0.2\n")
print("| target | static 0.005: quads | vs ref | march 0.005: computed | resident final | merged | vs ref | resolvable left | vs uniform | catch-up | march 0.2: computed | resident | merged | vs ref |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in ['preset_shape_h1', 'config_stability', 'preset_shape']:
    s5 = g('live', t, 'phase3_fine.txt', 'stationary=0 alpha_lo=0.005') or g('live', t, 'phase3_fine2.txt', 'stationary=0 alpha_lo=0.005')
    m5 = g('march', t, 'phase3_fine2.txt', 'alpha_lo=0.005')
    m2 = g('march', t, 'phase3_march2.txt', 'expiry-on-structure') or g('march', t, 'phase3_alo.txt', 'default (alpha_lo 0.2, merge on)')
    if not (s5 and m5 and m2): print(f"| {t} | (pending) |"); continue
    print(f"| {t} | {s5['quads']} | {e(s5)} | {m5['mem_all']} | {m5['res_final']} | {m5['merged']} | {e(m5)} | {e(m5,'resolvable')} | {m5['indicator']['uni_r']:.2f}x | {m5['catchup_pct']:.0f}% | {m2['mem_all']} | {m2['res_final']} | {m2['merged']} | {e(m2)} |")
print("\n### The dimension floor off, noise stop on, alpha_lo 0.2, arm on\n")
print("| target | sea | full depth: quads | floor on: quads | vs ref | resolvable left | vs uniform | floors | dim off: quads | vs ref | resolvable left | vs optimum | vs uniform | floors | quads saved |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in T:
    full = g('live', t, 'phase3_alo0.txt', 'stationary=0 alpha_lo=0 (guarded)') or g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0')
    on = g('live', t, 'phase3_alo.txt', 'stationary=0 alpha_lo=0.2'); nd = g('live', t, 'phase3_nodim.txt', 'stationary=0 alpha_lo=0.2 dim_floor=0')
    if not (full and on and nd): print(f"| {t} | (pending) |"); continue
    print(f"| {t} | {on['indicator']['sea']:.4f} | {full['quads']} | {on['quads']} | {e(on)} | {e(on,'resolvable')} | {on['indicator']['uni_r']:.2f}x | {fl(on)} | {nd['quads']} | {e(nd)} | {e(nd,'resolvable')} | {nd['indicator']['dp_r']:.2f}x | {nd['indicator']['uni_r']:.2f}x | {fl(nd)} | {100*(1-nd['quads']/full['quads']):.0f}% |")
print("\n### The live march with the dimension floor off\n")
print("| target | static dim off: quads | vs ref | march: computed | resident final | merged | vs ref | resolvable left | vs uniform | catch-up |")
print("|---|---|---|---|---|---|---|---|---|---|")
for t in ['preset_shape_h1', 'config_stability', 'preset_shape']:
    s_ = g('live', t, 'phase3_nodim.txt', 'stationary=0 alpha_lo=0.2 dim_floor=0'); m = g('march', t, 'phase3_nodim.txt', 'dim_floor=0')
    if not (s_ and m): print(f"| {t} | (pending) |"); continue
    print(f"| {t} | {s_['quads']} | {e(s_)} | {m['mem_all']} | {m['res_final']} | {m['merged']} | {e(m)} | {e(m,'resolvable')} | {m['indicator']['uni_r']:.2f}x | {m['catchup_pct']:.0f}% |")

print("\n### The live march at 0.005 after the cap re-decision fix\n")
print("| target | static 0.005: quads | vs ref | march before fix: computed | resident | merged | vs ref | march after fix: computed | resident | merged | vs ref | resolvable left | vs uniform | catch-up |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in ['preset_shape_h1', 'config_stability', 'preset_shape']:
    s5 = g('live', t, 'phase3_fine.txt', 'stationary=0 alpha_lo=0.005') or g('live', t, 'phase3_fine2.txt', 'stationary=0 alpha_lo=0.005')
    m5 = g('march', t, 'phase3_fine2.txt', 'alpha_lo=0.005'); m3 = g('march', t, 'phase3_march3.txt', 'alpha_lo=0.005 capfix')
    if not (s5 and m5 and m3): print(f"| {t} | (pending) |"); continue
    print(f"| {t} | {s5['quads']} | {e(s5)} | {m5['mem_all']} | {m5['res_final']} | {m5['merged']} | {e(m5)} | {m3['mem_all']} | {m3['res_final']} | {m3['merged']} | {e(m3)} | {e(m3,'resolvable')} | {m3['indicator']['uni_r']:.2f}x | {m3['catchup_pct']:.0f}% |")
