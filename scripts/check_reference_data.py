#!/usr/bin/env python
"""Check an unpacked copy of the Boileau et al. 2015 supplementary data (doi:10.1002/cnm.2732).

The data cannot be redistributed here and the publisher does not allow scripted download, so
download the supplementary archive from the article page by hand, unpack it, then run

    python scripts/check_reference_data.py /path/to/SuppMaterial_1DBenchmark_Boileau_etal

It checks that the six benchmark folders are present and that a few key files match the copy
used to validate this package, and prints the export line for ARTERIAL1D_BOILEAU_DIR.
"""
import hashlib
import sys
from pathlib import Path

SHA256 = {
    "Benchmark1_SinglePulseModel/NumData/FVM_inviscid.txt": "b9bde8e113497909ef6dacf3c4d77e5dc8aba1574be2ca701cfe351fae486cfb",
    "Benchmark2_CommonCarotidArtery/NumData/FVM.txt": "8b84bcdbc91051950f2c88fd9bab665e61acbfc495d60845f3a9828e8b6e86fb",
    "Benchmark3_UpperThoracicAorta/NumData/FVM.txt": "28eaafb04b7dcb5e208e8981a1735fcd64ef27881f9fea95a00adbbfe97c3eea",
    "Benchmark4_AorticBifurcation/NumData/FVM_AortaMid.txt": "267cbb8f6c80ebd244a2112183c11877d26ac569f49dedf4d601bd674e200cfc",
    "Benchmark5_37ArteryNetwork/NumData/FVM_P.txt": "70e26e1502dd8e6af22d0ade462f59c3b2a79532def4aa1a3422bda061a2fda9",
    "Benchmark6_ADAN56/NumData/FVM_P.txt": "de1b76ad76abf0a01f9b9c550da824c9fdcd71134f68cf40c5b9497886571a1b",
    "Benchmark6_ADAN56/NumData/FVM_Q.txt": "4e50ed40946504cfa04f1f38dc761a7eb49b7e242506c6d8fdabe260ed1cf007",
    "Benchmark6_ADAN56/inflow.txt": "cfa8cf2850b51039128a7e927c921f192453e35b223f61d550d82b4308deab49",
}


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    root = Path(sys.argv[1]).expanduser().resolve()
    hits = list(root.rglob("Benchmark6_ADAN56"))
    if not (root / "Benchmark6_ADAN56").is_dir() and hits:
        root = hits[0].parent                         # accept the folder above the unpacked one
    bad = 0
    for rel, h in SHA256.items():
        p = root / rel
        if not p.exists():
            print(f"missing   {rel}")
            bad += 1
        elif hashlib.sha256(p.read_bytes()).hexdigest() != h:
            print(f"differs   {rel}  (a different release of the data?)")
            bad += 1
        else:
            print(f"ok        {rel}")
    if bad:
        sys.exit(f"\n{bad} problem(s) in {root}")
    print(f"\nAll present. Use:\n  export ARTERIAL1D_BOILEAU_DIR={root}")


if __name__ == "__main__":
    main()
