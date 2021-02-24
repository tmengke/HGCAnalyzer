#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

//g++ -o ntupleRead ntupleRead.cpp `root-config --cflags --glibs`

int main()
{
	//auto h1f = new TH1F("h1","lc_energy",100,0,200);
	auto h2f = new TH2F("h2","z vs lc_energy", 100, 200, 500, 100, 0, 200);

	auto tf = TFile::Open("root://cmsxrootd.fnal.gov//store/user/tmengke/hgcdpg/samples/2photons_E100_eta2p0_nopu/31288926/ntuple/ntuple_0.root");
   	if (!tf || tf->IsZombie()) {
     		 return 0;
   	}
	
	TTreeReader mytree("hgcAnalyzer/tree", tf);
	TTreeReaderValue<std::vector<double>> lcEnergy(mytree, "lcEnergy");
	TTreeReaderValue<std::vector<ROOT::Math::XYZPointF>> lcPosition(mytree, "lcPosition");
	TTreeReaderValue<std::vector<ROOT::Math::PxPyPzEVector>> genParticles(mytree, "genParticles");


	//"lcHits" store the indices of recHits (link between layerCluster and recHits)
        //"trackster*_vertices" store the indices of layerCluster (link between trackster and layerClusters)
	while (mytree.Next()) {
		//loop layerClusters
		int lcidx=0;
		for (float x: *lcEnergy){
			//h1f->Fill(x);
			h2f->Fill(lcPosition->at(lcidx).Z(),x);
			lcidx++;
		}
 		//loop genParticles
		//for (auto gp: *genParticles){
		//	std::cout<<"genParticles Eta "<<gp.Eta()<<std::endl;
		//}

  	}
	auto c1 = new TCanvas("c1", "test", 200, 10, 900, 700 );
	//h1f->Draw();
	h2f->Draw("colz");
	//c1->SetLogy();
	c1->SaveAs("test.pdf");
	
	return 0;
}
