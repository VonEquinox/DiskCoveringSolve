"""Create a separate flat workspace for the optional discovery scripts.

This does not modify the proof inputs. Numerical discovery is not required for
verification. Requires only Python's standard library to prepare the workspace.
"""
from pathlib import Path
import argparse, shutil

def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path,help='New or empty research directory')
    args=parser.parse_args();out=args.directory.resolve()
    package=Path(__file__).resolve().parents[1]
    if out==package or package in out.parents:
        raise ValueError('Choose a research directory outside the proof package')
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('Refusing to overwrite a nonempty research directory')
    out.mkdir(parents=True,exist_ok=True)
    for source in (package/'core',package/'enumeration',package/'discovery'):
        for p in source.iterdir():
            if p.is_file() and p.name not in {'prepare_workspace.py','README.md'} and p.suffix in {'.py','.cpp','.hpp','.json','.gz'}:
                shutil.copy2(p,out/p.name)
    print(f'Discovery workspace: {out}')
    print('No certificate has been accepted by this setup operation.')

if __name__=='__main__':
    main()
