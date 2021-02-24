import ROOT
import math
from __future__ import division
tf = ROOT.TFile.Open("root://cmsxrootd.fnal.gov//store/user/tmengke/hgcdpg/samples/2photons_E100_eta2p0_nopu/31288926/ntuple/ntuple_0.root”)
tree = tf.Get("hgcAnalyzer/tree”)
h1f = ROOT.TH1F('h1f', 'test', 100, 0, 200)

#"lcHits" store the indices of recHits (link between layerCluster and recHits)
#"trackster*_vertices" store the indices of layerCluster (link between trackster and layerClusters)

for it,t in enumerate(tree):
    for x in t.lcEnergy:
        h1f.Fill(x)
c1 = ROOT.TCanvas( 'c1', 'test', 200, 10, 900, 700 )
h1f.Draw()
c1.SetLogy()
c1.Draw()
c1.SaveAs("test.pdf");
