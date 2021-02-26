from __future__ import division
import ROOT
import math

tf = ROOT.TFile.Open("root://cmsxrootd.fnal.gov//store/user/tmengke/hgcdpg/samples/2photons_E100_eta2p0_nopu/31288926/ntuple/ntuple_0.root")
tree = tf.Get("hgcAnalyzer/tree")

def dist(v1,v2):
    return math.sqrt((v1.X()-v2.X())**2+(v1.Y()-v2.Y())**2+(v1.Z()-v2.Z())**2)

#h1f = ROOT.TH1F('h1f', 'test', 100, 0, 200)
h2f = ROOT.TH2F('h2f', 'test', 100, 0, 12, 100, 0, 230)

#"lcHits" store the indices of recHits (link between layerCluster and recHits)
#"trackster*_vertices" store the indices of layerCluster (link between trackster and layerClusters)

for it,t in enumerate(tree):
    #distance between two genParticles from close-by particle gun
    di = dist(t.genParticlePosition[0],t.genParticlePosition[1])
    #for x in t.lcEnergy:
    #    h1f.Fill(x)
    E_total=0
    #total trackster energy
    for ix, x in enumerate(t.tracksterMerge_raw_energy):
    	E_total+=x
    h2f.Fill(di,E_total)  

c1 = ROOT.TCanvas( 'c1', 'test', 200, 10, 900, 700 )
h2f.Draw('COLZ')
h2f.SetTitle('trackster merge raw energy vs di-pions distance (100GeV)')
h2f.GetXaxis().SetTitle("distance[cm]");
h2f.GetYaxis().SetTitle("energy[GeV]");
#c1.SetLogy()
c1.Draw()
c1.SaveAs("test.pdf");
