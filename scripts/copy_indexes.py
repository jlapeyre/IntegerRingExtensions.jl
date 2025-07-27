#!/usr/bin/env python3
from pathlib import Path
import shutil
import re
import sys

###
### Copy doc index pages from packages to the toplevel docs.
###
### Each package in ../lib has its own ./docs directory so that standalone
### docs, such as they are, can be built. This script copies the doc source
### from the packages to the toplevel ./docs directory. This allows building
### docs at the toplevel that include all packages.
###

PACKAGES = [
    "Angles2",
    "CyclotomicRings",
    "Dyadics",
    "QuadraticRings",
    "RingExtensionsCommon",
    "RingExtensionsUtils",
    "RootOnes",
    "SingletonNumbers",
    "SmallStaticMatrices",
]

def camel_to_snake(name: str) -> str:
    s1 = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    s2 = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s1)
    return s2.lower()

def main() -> int:
    root = Path(__file__).resolve().parent.parent
    print(f"root is {root}")
    dest_dir = root / "docs" / "src" / "packages"
    dest_dir.mkdir(parents=True, exist_ok=True)

    for pkg in PACKAGES:
        src = root / "lib" / pkg / "docs" / "src" / "index.md"
        if not src.is_file():
            print(f"warn: missing {src}", file=sys.stderr)
            continue
        dest = dest_dir / f"{camel_to_snake(pkg)}.md"
        shutil.copyfile(src, dest)
        print(f"{src} -> {dest}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
