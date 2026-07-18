#!/usr/bin/env python3
"""Reusable open-family certificate driver for planar one-negative-minor networks.

Input is a JSON configuration naming a witness, a compact radius R, an angular
grid, and either a common or anisotropic parameter box.  The driver runs:
  1. exact asymptotic-polygon geometry checks;
  2. primary compact parameter interval certification;
  3. analytic tail transfer;
  4. optional independent Decimal replay of the exported compact leaves;
  5. optional independent Decimal replay of the analytic tail loss.

The current tail lemma applies to planar networks whose Cauchy-Binet expansion
has exactly one negative pair and for which every direction is dominated by at
least one positive pair.  The witness is the first packaged application; the
configuration interface is not hard-coded to its numerical entries.
"""
from __future__ import annotations
import argparse,json,subprocess,sys
from pathlib import Path
from report_validation import ReportValidationError, build_combined

def run(cmd,stdout=None,allowed=(0,)):
    print('+',' '.join(map(str,cmd)),flush=True)
    result=subprocess.run(cmd,check=False,text=True,stdout=stdout)
    if result.returncode not in allowed:
        raise subprocess.CalledProcessError(result.returncode,cmd)
    return result

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('config'); ap.add_argument('--skip-independent',action='store_true'); a=ap.parse_args()
    config_path=Path(a.config).resolve(); cfg=json.load(open(config_path)); root=Path(__file__).resolve().parents[1]; py=sys.executable
    witness=str((config_path.parent/Path(cfg['witness'])).resolve()) if not Path(cfg['witness']).is_absolute() else cfg['witness']
    out=root/'results'/cfg.get('name','open_family'); out.mkdir(parents=True,exist_ok=True)
    da=str(cfg.get('delta_A',cfg.get('delta','0'))); db=str(cfg.get('delta_B',cfg.get('delta','0'))); dc=str(cfg.get('delta_c',cfg.get('delta','0')))
    common=['--delta','0','--delta-a',da,'--delta-b',db,'--delta-c',dc,'--witness',witness]
    # Geometry is exact for the packaged witness verifier. General configurations
    # may supply their own exact geometry command.
    if cfg.get('geometry_command'):
        run(cfg['geometry_command'])
    elif cfg.get('use_witness_geometry',True):
        run([py,str(root/'code'/'verify_witness_polygon.py'),'--witness',witness,'--json-out',str(out/'geometry.json')])
    leaves=out/'compact_leaf_paths.json'; state=out/'compact_state.json'; compact=out/'compact.json'
    # A normal local run can omit chunking. chunk_seconds is useful in constrained environments.
    cmd=[py,str(root/'code'/'certify_parameter_radius.py'),*common,'--R',str(cfg.get('R',582)),'--grid',str(cfg.get('grid',200000)),'--hmin',str(cfg.get('hmin','0.008')),'--skip-tail','--leaf-out',str(leaves),'--state-file',str(state),'--chunk-seconds',str(cfg.get('chunk_seconds',0)),'--json-out',str(compact)]
    run(cmd,allowed=(0,2))
    cres=json.load(open(compact))
    if cres['verdict']=='PAUSED':
        print('PAUSED: rerun the same command to resume'); return 2
    run([py,str(root/'code'/'certify_tail_parameter_transfer.py'),*common,'--R',str(cfg.get('R',582)),'--grid',str(cfg.get('grid',200000)),'--out',str(out/'tail.json')])
    if not a.skip_independent:
        run([py,str(root/'code'/'replay_leaf_paths_decimal.py'),*common,'--paths',str(leaves),'--workers',str(cfg.get('workers',8)),'--json-out',str(out/'compact_independent.json')])
        run([py,str(root/'code'/'independent_tail_transfer.py'),*common,'--R',str(cfg.get('R',582)),'--grid',str(cfg.get('grid',200000)),'--json-out',str(out/'tail_independent.json')])
        run([py,str(root/'code'/'replay_leaf_paths_arb.py'),str(config_path),'--workers',str(cfg.get('workers',8)),'--json-out',str(out/'compact_arb.json')])
    else:
        print('Primary reports generated; independent reports and aggregate were intentionally skipped')
        return 0
    try:
        combined=build_combined(out,config_path)
    except ReportValidationError as exc:
        raise SystemExit(f'semantic report validation failed: {exc}') from exc
    (out/'combined.json').write_text(json.dumps(combined,indent=2)+'\n'); print(json.dumps(combined,indent=2)); return 0
if __name__=='__main__': raise SystemExit(main())
