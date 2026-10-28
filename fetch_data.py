"""Download the 35 NHANES public-use files analysed (5 cycles x 7 files) into data/ and check each MD5.

The checksums are those of the files analysed (downloaded 2026-10-23 to 2026-10-27). NCHS occasionally revises a
public-use file; a mismatch is reported, not overwritten silently.

  python3 fetch_data.py
"""
import hashlib
import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.abspath(__file__))
YEAR = {'D': '2005', 'E': '2007', 'F': '2009', 'H': '2013', 'J': '2017'}
URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{name}"
MD5 = {
    "DEMO_D.xpt": "7d61923d72db90338b2d75863958f9a0",
    "DXXFEM_D.xpt": "38be7519ca85d0e181f0930c75a97d9b",
    "DXXSPN_D.xpt": "c1ebfd3b8c4a977c20f649f6a75f88c6",
    "RHQ_D.xpt": "a61ad7b2ec9b31c4491ad511966e3b16",
    "OSQ_D.xpt": "dcce2e09fac6d93b7bb49af6bfe5d321",
    "BMX_D.xpt": "336e8aa5e9397aff208d3d69def13153",
    "SMQ_D.xpt": "f5ccab52b9119f89fc85409e1b70eda0",
    "DEMO_E.xpt": "d505cb9087e3597cda5bef9dd4863d71",
    "DXXFEM_E.xpt": "3ab5d899ec55f174441c5dfbaa0ce5b5",
    "DXXSPN_E.xpt": "c1fbc8fcc49fa9fa42e8c44385d00b37",
    "RHQ_E.xpt": "bbc8dfad66e4a1642a80b344e060e334",
    "OSQ_E.xpt": "295fceed3bdd2b4cddde78abaaa4bf84",
    "BMX_E.xpt": "de4d6ae686b3c48f1d90fde505103303",
    "SMQ_E.xpt": "23cd23a1f5f11d538ac023e95dcccf07",
    "DEMO_F.xpt": "aa0311654076bcf7dea333389f3c0635",
    "DXXFEM_F.xpt": "81164955195c6b4e75ad79f1429d40c4",
    "DXXSPN_F.xpt": "c4a0dc49f4eb6d614f65d5b1648cef7e",
    "RHQ_F.xpt": "1990751c726126d90ac5c7b98221f3a4",
    "OSQ_F.xpt": "b86f7f8d0b6e0e01cb402a4c47854844",
    "BMX_F.xpt": "016dc70417a802e659377c604b8dc1c5",
    "SMQ_F.xpt": "6dfbd9a8307c26121a769dd3485f301f",
    "DEMO_H.xpt": "fe59a58ba1a5cc4814d7e6db4ebe820e",
    "DXXFEM_H.xpt": "56cbb66dc10453a855165737fce6cfd0",
    "DXXSPN_H.xpt": "256c47581f36942a0d04d5ca88f3dedb",
    "RHQ_H.xpt": "300068389fcee85c0f5ecaac237f8565",
    "OSQ_H.xpt": "de6af558b270d0cb1b0d4948bd2b9725",
    "BMX_H.xpt": "c5bc2b552c8cd33f10557fd16b747920",
    "SMQ_H.xpt": "4e4ea631cc9d9672a6bc9c30d0986082",
    "DEMO_J.xpt": "3d5e48e4255b2ec80578a0eb45480d4c",
    "DXXFEM_J.xpt": "938b6e28e6f3a4e753a16fb8f5829f0b",
    "DXXSPN_J.xpt": "59336e446b57a5e3a62eadc146ed0bda",
    "RHQ_J.xpt": "f2b0ac697a451a29546876bac450434f",
    "OSQ_J.xpt": "32469ffd2d2a4b84c1bdf65c74f3fa78",
    "BMX_J.xpt": "4c51e08c3b105005a326da3a4fc6b223",
    "SMQ_J.xpt": "2a4fb99269a678cfc9ea7b0bd77a6eae",
}


def main():
    os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
    bad = []
    for name, md5 in MD5.items():
        path = os.path.join(ROOT, "data", name)
        if not os.path.exists(path):
            year = YEAR[name.split("_")[1][0]]
            req = urllib.request.Request(URL.format(year=year, name=name), headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req) as r, open(path, "wb") as f:
                f.write(r.read())
        got = hashlib.md5(open(path, "rb").read()).hexdigest()
        print(("ok      " if got == md5 else "CHANGED ") + name)
        if got != md5:
            bad.append(name)
    if bad:
        sys.exit(f"{len(bad)} file(s) differ from the versions analysed: {', '.join(bad)}")


if __name__ == "__main__":
    main()
