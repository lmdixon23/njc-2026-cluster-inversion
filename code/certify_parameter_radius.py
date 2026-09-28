#!/usr/bin/env python3
"""Fail-closed parameter-box certificate for a common l-infinity radius.

For a requested delta >= 0, every entry of A, B, c is replaced by the interval
[theta0-delta, theta0+delta].  The script certifies, uniformly for all parameters
in that box:
  (i) the Cauchy-Binet minor sign pattern has exactly one negative product (2,3),
  (ii) det DF > 0 outside radius R using the same radial domination lemma as the
       point certificate, with all constants evaluated over parameter intervals,
  (iii) det DF > 0 on [-R,R]^2 by interval branch-and-bound.

A PASS is rigorous relative to mpmath.iv's outward interval operations. Any
unresolved leaf, sign ambiguity, or budget exhaustion is a FAIL. The separate
replay scripts provide partition and arithmetic cross checks.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, math, os, sys, time
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path
import mpmath
from mpmath import iv, mp
from report_validation import validate_planar_four, require




def _strip_timing(payload):
    """Return a copy with wall-clock fields removed, recursively.

    Persisted artifacts are content-addressed in SHA256SUMS.txt, so their bytes
    must depend only on the mathematics.  Timings stay on stdout.
    """
    if isinstance(payload, dict):
        return {k: _strip_timing(v) for k, v in payload.items()
                if k != "elapsed_seconds"}
    if isinstance(payload, list):
        return [_strip_timing(v) for v in payload]
    return payload

def load_witness(path: Path):
    raw = path.read_bytes()
    data = json.loads(raw.decode('utf-8'), parse_float=Decimal)
    A = [[F(v) for v in row] for row in data['A']]
    B = [[F(v) for v in row] for row in data['B']]
    c = [F(v) for v in data['c']]
    validate_planar_four(A, B, c)
    return A, B, c, hashlib.sha256(raw).hexdigest()


def q_to_iv(q: F):
    return iv.mpf(q.numerator) / iv.mpf(q.denominator)


def decimal_radius(text: str) -> F:
    return F(Decimal(text))


SPLIT_RULE = 'longer-side; x on ties; child0=lower/left, child1=upper/right'


def canonical_hash(value: dict) -> str:
    payload=json.dumps(value,sort_keys=True,separators=(',',':')).encode('utf-8')
    return hashlib.sha256(payload).hexdigest()


def path_box_exact(path: str, R: F):
    if set(path)-{'0','1'}:
        raise ValueError('state contains a malformed split path')
    x1=y1=-R; x2=y2=R
    for bit in path:
        w=x2-x1; h=y2-y1
        if w>=h:
            midpoint=(x1+x2)/2
            if bit=='0': x2=midpoint
            else: x1=midpoint
        else:
            midpoint=(y1+y2)/2
            if bit=='0': y2=midpoint
            else: y1=midpoint
    return x1,y1,x2,y2


def exact_partition(paths: list[str]) -> bool:
    if not paths or any(not isinstance(p, str) or set(p)-{'0','1'} for p in paths):
        return False
    ordered=sorted(paths)
    if len(set(ordered))!=len(ordered):
        return False
    if any(ordered[i+1].startswith(ordered[i]) for i in range(len(ordered)-1)):
        return False
    depth=max(map(len,ordered))
    return sum(1<<(depth-len(path)) for path in ordered)==1<<depth


def certificate_context(kind,sha,deltaA,deltaB,deltaC,R,grid,hmin,budget,dps):
    return {
      'schema':'parameter-radius-context-v2','kind':kind,'engine':'certify_parameter_radius.py-v2',
      'witness_sha256':sha,'delta_A':str(deltaA),'delta_B':str(deltaB),'delta_c':str(deltaC),
      'R':int(R),'grid':int(grid),'hmin':str(hmin),'budget':int(budget),'precision_dps':int(dps),
      'split_rule':SPLIT_RULE,
    }


def box_q(q: F, delta: F):
    lo, hi = q-delta, q+delta
    L, H = q_to_iv(lo), q_to_iv(hi)
    return iv.mpf([L.a, H.b])


def iv_abs(X):
    if X.a >= 0: return X
    if X.b <= 0: return -X
    return iv.mpf([0, max(-X.a, X.b)])


def abs_upper(X):
    return max(abs(mp.mpf(X.a)), abs(mp.mpf(X.b)))


def sigma_prime_iv(T):
    return 1/(2 + iv.exp(T) + iv.exp(-T))


def parameter_intervals(A,B,c,deltaA,deltaB=None,deltaC=None):
    validate_planar_four(A, B, c)
    deltaB = deltaA if deltaB is None else deltaB
    deltaC = deltaA if deltaC is None else deltaC
    require(all(F(d) >= 0 for d in (deltaA, deltaB, deltaC)), 'parameter radii must be nonnegative')
    AI = [[box_q(v,deltaA) for v in row] for row in A]
    BI = [[box_q(v,deltaB) for v in row] for row in B]
    cI = [box_q(v,deltaC) for v in c]
    return AI, BI, cI


def minor_products_iv(AI,BI):
    validate_planar_four(AI, BI)
    out={}
    for i,j in itertools.combinations(range(4),2):
        dA=AI[0][i]*AI[1][j]-AI[0][j]*AI[1][i]
        dB=BI[i][0]*BI[j][1]-BI[i][1]*BI[j][0]
        out[(i,j)]=dA*dB
    return out


def sign_structure(mI):
    pos=[]; neg=[]; amb=[]
    for I,M in mI.items():
        if M.a>0: pos.append(I)
        elif M.b<0: neg.append(I)
        else: amb.append(I)
    return pos,neg,amb


def det_lower_positive(box, BI, cI, mI):
    x1,y1,x2,y2=box
    X=iv.mpf([q_to_iv(x1).a,q_to_iv(x2).b])
    Y=iv.mpf([q_to_iv(y1).a,q_to_iv(y2).b])
    sp=[sigma_prime_iv(BI[i][0]*X+BI[i][1]*Y+cI[i]) for i in range(4)]
    acc=iv.mpf(0)
    for (i,j),M in mI.items():
        acc += M*sp[i]*sp[j]
    return bool(acc.a>0), mp.mpf(acc.a)


def row_norm_upper(BI,i):
    ax=abs_upper(BI[i][0]); ay=abs_upper(BI[i][1])
    return mp.sqrt(ax*ax+ay*ay)*(1+mp.mpf(2)**(-20))


def tail_setup(BI,cI,mI,pos,neg):
    In=neg[0]
    mneg=-mI[In]  # positive interval
    def rho(i,ct,st):
        return iv_abs(BI[i][0]*ct+BI[i][1]*st)
    K={}
    for J in pos:
        # interval ratio; denominator positive by sign gate
        K[J]=iv.log(16*mneg/mI[J])
        for idx in (J[0],J[1],In[0],In[1]):
            K[J]+=iv_abs(cI[idx])
    L={J: row_norm_upper(BI,In[0])+row_norm_upper(BI,In[1])+row_norm_upper(BI,J[0])+row_norm_upper(BI,J[1]) for J in pos}
    return In,rho,K,L


def tail_certify(R,BI,cI,mI,pos,neg,Mgrid=200000,early=True,state_file='',chunk_seconds=0):
    require(type(R) is int and R > 0 and type(Mgrid) is int and Mgrid > 0, 'R and grid must be positive integers')
    require(len(BI) == len(cI) == 4, 'tail certificate requires four ridges')
    if state_file and Path(state_file).exists():
        raise ValueError('tail resume files are disabled for this release; rerun the tail from direction zero')
    if state_file and chunk_seconds:
        raise ValueError('tail checkpointing is disabled because a compact state cannot prove prior angular coverage')
    In,rho,K,L=tail_setup(BI,cI,mI,pos,neg)
    Rm=mp.mpf(R)
    half_arc=(mp.pi/Mgrid)*(1+mp.mpf(2)**(-30))
    worst=mp.inf; worstk=None; start=0; t0=time.time()
    for k in range(start,Mgrid):
        th=iv.mpf(2)*iv.pi*k/Mgrid
        ct,st=iv.cos(th),iv.sin(th)
        best=-mp.inf
        for J in pos:
            g=(rho(In[0],ct,st)+rho(In[1],ct,st)-rho(J[0],ct,st)-rho(J[1],ct,st))
            glo=mp.mpf(g.a)-L[J]*half_arc
            if glo<=0: continue
            margin=Rm*glo-mp.mpf(K[J].b)
            if margin>best: best=margin
        if best<worst:
            worst=best; worstk=k
        if early and best<=0:
            return False,worst,worstk
        if chunk_seconds and time.time()-t0>=chunk_seconds:
            return None,worst,worstk
    return bool(worst>0),worst,worstk


def interior_certify(R,BI,cI,mI,hmin=F('0.008'),budget=40_000_000,progress=0,leaf_out='',state_file='',chunk_seconds=0,state_context=None):
    require(F(R) > 0 and F(hmin) > 0 and type(budget) is int and budget > 0, 'invalid compact domain or budget')
    require(len(BI) == len(cI) == 4, 'compact certificate requires four ridges')
    Rf=F(R)
    root=(-Rf,-Rf,Rf,Rf)
    stack=[(root,'')]; processed=certified=0; fails=[]; worst=mp.inf; leaves=[]
    if state_file and Path(state_file).exists():
        st=json.loads(Path(state_file).read_text())
        if st.get('schema')!='parameter-radius-interior-state-v2':
            raise ValueError('interior state has an unsupported schema')
        if not state_context or st.get('context')!=state_context or st.get('context_sha256')!=canonical_hash(state_context):
            raise ValueError('interior state is not bound to the selected witness and configuration')
        processed=int(st['processed']); certified=int(st['certified']); leaves=list(st.get('paths',[]))
        stack=[]
        for rec in st['stack']:
            box=tuple(F(z) for z in rec['box']); stack.append((box,rec['path']))
            if box!=path_box_exact(rec['path'],Rf):
                raise ValueError('interior state box does not match its split path')
        if not stack:
            raise ValueError('finished or empty interior state is invalid; completed runs delete their state')
        frontier=leaves+[path for _,path in stack]
        if certified!=len(leaves) or processed!=2*certified+len(stack)-1 or not exact_partition(frontier):
            raise ValueError('interior state does not prove exact prefix-free coverage of the root')
        # A context hash and a complete partition do not prove the stored leaves.
        # Rescore them under this run's exact inputs before trusting a resume.
        for path in leaves:
            ok, lo = det_lower_positive(path_box_exact(path, Rf), BI, cI, mI)
            if not ok:
                raise ValueError('resumed certified leaf failed arithmetic replay: ' + path)
            worst = min(worst, lo)
    t0=time.time()
    while stack and processed<budget:
        box,path=stack.pop(); processed+=1
        ok,lo=det_lower_positive(box,BI,cI,mI)
        worst=min(worst,lo)
        if ok:
            certified+=1
            if leaf_out or state_file: leaves.append(path)
        else:
            x1,y1,x2,y2=box; w=x2-x1; h=y2-y1
            if w<hmin and h<hmin:
                fails.append({'box':[str(z) for z in box],'lower':str(lo)})
                break
            elif w>=h:
                mx=(x1+x2)/2
                stack.append(((x1,y1,mx,y2),path+'0')); stack.append(((mx,y1,x2,y2),path+'1'))
            else:
                my=(y1+y2)/2
                stack.append(((x1,y1,x2,my),path+'0')); stack.append(((x1,my,x2,y2),path+'1'))
        if progress and processed%progress==0:
            print(f'  processed={processed:,} stack={len(stack):,} cert={certified:,} elapsed={time.time()-t0:.1f}s',flush=True)
        if chunk_seconds and time.time()-t0 >= chunk_seconds:
            break
    paused=bool(stack) and not fails and bool(chunk_seconds) and processed<budget
    budget_exhausted=bool(stack) and not fails and processed>=budget
    if state_file:
        if paused:
            payload={'schema':'parameter-radius-interior-state-v2','context':state_context,'context_sha256':canonical_hash(state_context),'processed':processed,'certified':certified,'paths':leaves,'worst':str(worst),'stack':[{'box':[str(z) for z in box],'path':path} for box,path in stack]}
            Path(state_file).write_text(json.dumps(payload,separators=(',',':'))+'\n', newline="\n")
        elif Path(state_file).exists():
            Path(state_file).unlink()
    if leaf_out and not stack and not fails:
        Path(leaf_out).write_text(json.dumps({'R':str(R),'split_rule':SPLIT_RULE,'paths':leaves},separators=(',',':'))+'\n', newline="\n")
    return {
        'passed': not fails and not stack,
        'paused': paused,
        'processed': processed,
        'certified_leaves': certified,
        'unresolved': len(fails),
        'budget_exhausted': budget_exhausted,
        'first_failures': fails,
        'elapsed_seconds': time.time()-t0,
        'smallest_seen_lower': str(worst),
    }


def certify(delta_text,witness,R,grid,hmin,budget,progress,skip_tail=False,skip_interior=False,leaf_out='',state_file='',chunk_seconds=0,delta_a_text='',delta_b_text='',delta_c_text='',tail_state_file='',tail_chunk_seconds=0,dps=40):
    require(not (skip_tail and skip_interior), 'cannot skip both certificate domains')
    require(type(R) is int and R > 0 and type(grid) is int and grid > 0, 'R and grid must be positive integers')
    require(type(dps) is int and dps >= 30, 'precision must be at least 30 decimal digits')
    require(F(hmin) > 0 and type(budget) is int and budget > 0, 'invalid resolution or budget')
    require(math.isfinite(chunk_seconds) and chunk_seconds >= 0 and math.isfinite(tail_chunk_seconds) and tail_chunk_seconds >= 0, 'invalid chunk duration')
    mp.dps = iv.dps = dps
    delta=decimal_radius(delta_text)
    require(delta >= 0, 'parameter radius must be nonnegative')
    deltaA=decimal_radius(delta_a_text) if delta_a_text else delta
    deltaB=decimal_radius(delta_b_text) if delta_b_text else delta
    deltaC=decimal_radius(delta_c_text) if delta_c_text else delta
    A,B,c,sha=load_witness(Path(witness))
    AI,BI,cI=parameter_intervals(A,B,c,deltaA,deltaB,deltaC)
    mI=minor_products_iv(AI,BI)
    pos,neg,amb=sign_structure(mI)
    result={'delta':delta_text,'delta_A':str(deltaA),'delta_B':str(deltaB),'delta_c':str(deltaC),'witness_sha256':sha,'R':R,'grid':grid,'hmin':str(hmin)}
    result['scope'] = 'compact-only' if skip_tail else 'tail-only' if skip_interior else 'compact-and-tail'
    result['minor_intervals']={str(k):[str(v.a),str(v.b)] for k,v in mI.items()}
    result['minor_signs']={'positive':pos,'negative':neg,'ambiguous':amb}
    if amb or neg!=[(2,3)] or len(pos)!=5:
        result['verdict']='FAIL'; result['failure']='minor sign pattern not uniform'
        return result
    if not skip_tail:
        t0=time.time(); ok,margin,k=tail_certify(R,BI,cI,mI,pos,neg,Mgrid=grid,state_file=tail_state_file,chunk_seconds=tail_chunk_seconds)
        result['tail']={'passed':ok,'worst_margin':str(margin),'worst_grid_index':k,'elapsed_seconds':time.time()-t0}
        if ok is None:
            result['verdict']='PAUSED'; return result
        if not ok:
            result['verdict']='FAIL'; result['failure']='tail'
            return result
    if not skip_interior:
        context=certificate_context('interior',sha,deltaA,deltaB,deltaC,R,grid,hmin,budget,dps)
        interior=interior_certify(R,BI,cI,mI,hmin=hmin,budget=budget,progress=progress,leaf_out=leaf_out,state_file=state_file,chunk_seconds=chunk_seconds,state_context=context)
        result['interior']=interior
        if interior.get('paused'):
            result['verdict']='PAUSED'; return result
        if not interior['passed']:
            result['verdict']='FAIL'; result['failure']='interior'
            return result
    result['verdict']='PASS'
    return result


def main():
    ap=argparse.ArgumentParser()
    root=Path(__file__).resolve().parents[1]
    ap.add_argument('--witness',default=str(root/'data'/'witness.json'))
    ap.add_argument('--delta',default='0')
    ap.add_argument('--delta-a',default='')
    ap.add_argument('--delta-b',default='')
    ap.add_argument('--delta-c',default='')
    ap.add_argument('--R',type=int,default=582)
    ap.add_argument('--grid',type=int,default=200000)
    ap.add_argument('--hmin',default='0.008')
    ap.add_argument('--budget',type=int,default=40000000)
    ap.add_argument('--progress',type=int,default=50000)
    ap.add_argument('--dps',type=int,default=40)
    ap.add_argument('--skip-tail',action='store_true')
    ap.add_argument('--skip-interior',action='store_true')
    ap.add_argument('--json-out',default='')
    ap.add_argument('--leaf-out',default='')
    ap.add_argument('--state-file',default='')
    ap.add_argument('--chunk-seconds',type=float,default=0)
    ap.add_argument('--tail-state-file',default='')
    ap.add_argument('--tail-chunk-seconds',type=float,default=0)
    args=ap.parse_args(); mp.dps=args.dps
    try:
        res=certify(args.delta,args.witness,args.R,args.grid,F(args.hmin),args.budget,args.progress,args.skip_tail,args.skip_interior,args.leaf_out,args.state_file,args.chunk_seconds,args.delta_a,args.delta_b,args.delta_c,args.tail_state_file,args.tail_chunk_seconds,args.dps)
    except (ValueError,KeyError,TypeError,json.JSONDecodeError) as exc:
        res={'verdict':'FAIL','failure':f'invalid checkpoint or input: {exc}'}
    print(json.dumps(res,indent=2))
    if args.json_out: Path(args.json_out).write_text(json.dumps(_strip_timing(res),indent=2)+'\n', newline="\n")
    return 0 if res['verdict']=='PASS' else 2 if res['verdict']=='PAUSED' else 1

if __name__=='__main__': raise SystemExit(main())
