#!/usr/bin/env python3
"""Stage 5 of the CAD->TGeo Geant integration test: compare the geantino material-budget
maps of the exact and tessellated representations.

Usage: compare_matbudget.py <exact/o2sim_matbudget.root> <mesh/o2sim_matbudget.root>

Two representations of the same CAD solid must present the same X0 integral along the
same ray -- this is the physics-level analogue of the X-ray benchmark, and unlike step
counts it has a right answer. Every TH1/TH2 present in both files is compared per bin;
report: mean/max absolute and relative difference, and the count of bins whose relative
difference exceeds 1% and 10%. Needs PyROOT (run inside the sim environment, NOT the
converter's pythonOCC environment).
"""
import sys
import math


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    import ROOT
    ROOT.gROOT.SetBatch(True)
    fa = ROOT.TFile.Open(sys.argv[1])
    fb = ROOT.TFile.Open(sys.argv[2])
    if not fa or fa.IsZombie() or not fb or fb.IsZombie():
        print("cannot open input files")
        return 1
    names = [k.GetName() for k in fa.GetListOfKeys()]
    print(f"{'hist':22s} {'entriesA':>10s} {'meanA':>12s} {'meanB':>12s} "
          f"{'max|d|':>10s} {'mean|d|':>10s} {'maxrel':>9s} {'>1%':>7s} {'>10%':>7s} {'bins':>8s}")
    status = 0
    for name in names:
        ha, hb = fa.Get(name), fb.Get(name)
        if not ha or not hb or not ha.InheritsFrom("TH1"):
            continue
        nx, ny = ha.GetNbinsX(), ha.GetNbinsY()
        if (nx, ny) != (hb.GetNbinsX(), hb.GetNbinsY()):
            print(f"{name:22s} BINNING MISMATCH")
            status = 1
            continue
        sum_a = sum_b = sum_ad = max_ad = max_rel = 0.0
        n_gt1 = n_gt10 = nbins = 0
        for ix in range(1, nx + 1):
            for iy in range(1, ny + 1):
                a = ha.GetBinContent(ix, iy)
                b = hb.GetBinContent(ix, iy)
                d = abs(a - b)
                ref = max(abs(a), abs(b))
                nbins += 1
                sum_a += a
                sum_b += b
                sum_ad += d
                max_ad = max(max_ad, d)
                if ref > 0:
                    rel = d / ref
                    max_rel = max(max_rel, rel)
                    if rel > 0.01:
                        n_gt1 += 1
                    if rel > 0.10:
                        n_gt10 += 1
        print(f"{name:22s} {ha.GetEntries():10.0f} {sum_a / nbins:12.5g} {sum_b / nbins:12.5g} "
              f"{max_ad:10.3g} {sum_ad / nbins:10.3g} {max_rel:9.3g} {n_gt1:7d} {n_gt10:7d} {nbins:8d}")
    return status


if __name__ == "__main__":
    sys.exit(main())
