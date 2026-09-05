#!/usr/bin/env python3
"""Prepare a SEPARATE flat directory for optional, untrusted discovery programs.
Does not overwrite proof inputs and does not mark a generated proposal verified.
"""
from pathlib import Path
import argparse,shutil

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('destination',type=Path,help='New or empty directory outside the proof package')
    args=p.parse_args();root=Path(__file__).resolve().parents[1];dst=args.destination.resolve()
    if root==dst or root in dst.parents:
        p.error('Use a directory outside the proof package')
    if dst.exists() and any(dst.iterdir()):p.error('Destination must be new or empty')
    dst.mkdir(parents=True,exist_ok=True)
    for name in ['geometry19.py','anchor_isolation19.py','interval19.py','exact_force19.py',
                 'exact_arcs.py','candidate_global19.py','anchor_isolation_certificate.json']:
        shutil.copyfile(root/'core'/name,dst/name)
    for path in Path(__file__).resolve().parent.glob('*.py'):
        if path.name!='prepare_workspace.py':shutil.copyfile(path,dst/path.name)
    shutil.copyfile(root/'enumeration'/'survivors19.json',dst/'metric19_residuals.json')
    shutil.copyfile(root/'enumeration'/'candidate19_map.json',dst/'candidate19_map.json')
    print('Untrusted proposal workspace:',dst)
    print('No proof input has been modified. Generated proposals require separate exact replay.')

if __name__=='__main__':main()
