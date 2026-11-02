"""Build a standalone, unadapted IQ-TREE from the pinned vendor source.

The native bridge/memory adapters are deliberately not included. OpenMP is
required: this pinned revision's standalone AliSim crashes in a sequential build
because its per-thread RNG arrays are initialized only when OpenMP is enabled.
"""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('target/alisim-reference'))
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--libomp-prefix', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    source = output / 'source'
    build = output / 'build'
    source.mkdir(parents=True, exist_ok=True)
    prefix = args.libomp_prefix
    if prefix is None and sys.platform == 'darwin':
        prefix = next((p for p in (Path('/opt/homebrew/opt/libomp'), Path('/usr/local/opt/libomp'))
                       if (p / 'include/omp.h').is_file()), None)
    omp = '' if prefix is None else (f'include_directories([=[{prefix.resolve() / "include"}]=])\n'
                                    f'link_directories([=[{prefix.resolve() / "lib"}]=])\n')
    cmake = f'''cmake_minimum_required(VERSION 3.20)
project(pf_alisim_reference LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
set(BUILD_LIB OFF CACHE BOOL "Standalone reference" FORCE)
set(IQTREE_FLAGS "omp" CACHE STRING "Unmodified upstream OpenMP initialization" FORCE)
set(USE_CMAPLE OFF CACHE BOOL "" FORCE)
set(USE_CMAPLE_AA OFF CACHE BOOL "" FORCE)
set(USE_LSD2 OFF CACHE BOOL "" FORCE)
set(USE_TERRAPHAST OFF CACHE BOOL "" FORCE)
set(USE_BOOSTER OFF CACHE BOOL "" FORCE)
{omp}add_subdirectory([=[{root / "dependencies/iqtree2"}]=] iqtree)
'''
    (source / 'CMakeLists.txt').write_text(cmake)
    subprocess.run(['cmake', '-S', str(source), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release'], check=True)
    subprocess.run(['cmake', '--build', str(build), '--target', 'iqtree2', '--parallel', str(args.jobs)], check=True)
    print(build / 'iqtree/iqtree2')


if __name__ == '__main__':
    main()
