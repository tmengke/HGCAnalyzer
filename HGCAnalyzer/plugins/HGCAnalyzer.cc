// -*- C++ -*-
//
// Package:    analysis/HGCAnalyzer
// Class:      HGCAnalyzer
//
/**\class HGCAnalyzer HGCAnalyzer.cc analysis/HGCAnalyzer/plugins/HGCAnalyzer.cc

 Description: create ntuple for HGCAL clusters/tracksters study

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tielige Mengke
//         Created:  Thu, 22 Oct 2020 16:50:33 GMT
//
//

// system include files
#include <memory>
#include <numeric>
#include <iomanip>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/Point3D.h"
#include "Math/Vector4Dfwd.h"

using namespace std;

//
// class declaration
//

class HGCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HGCAnalyzer(const edm::ParameterSet&);
  ~HGCAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  //void setRecHitTools(const hgcal::RecHitTools* recHitTools) { rechittools_ = recHitTools; }

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_;
  //edm::EDGetTokenT<std::vector<CaloParticle>> caloParticles_;
  //edm::EDGetTokenT<std::vector<Trackster>> tracksterstrkem_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersem_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterstrk_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstershad_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersmrg_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;
  edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> timelayercluster_;  
   

  //const hgcal::RecHitTools* rechittools_;
  hgcal::RecHitTools rhtools_;
  std::string detector_;
  // ----------member data ---------------------------
  edm::Service<TFileService> fs_;
  TTree* tree_;
  
  std::string s_tracksters[4] = {"tracksterEM", "tracksterHAD", "tracksterMerge", "tracksterTrk"};

  edm::EDGetTokenT<HGCRecHitCollection> recHitsEE_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsFH_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsBH_;

  std::vector<int> PdgId_;
  
  std::vector<float> genPx_;
  std::vector<float> genPy_;
  std::vector<float> genPz_;
  std::vector<float> genEnergy_;
  std::vector<float> genEta_;
  std::vector<float> genPhi_;
  std::vector<float> genPt_;  
  
  std::vector<ROOT::Math::PxPyPzEVector> Particles_;

  std::vector<ROOT::Math::XYZPointF> lcPosition_;
  std::vector<float> lcEnergy_;
  std::vector<uint32_t> lcLayerN_;
  std::vector<float> lcTime_;
  std::vector<float> lcTimeError_;
  std::vector<int> lcId_;
  std::vector<std::vector<int>> lcHits_;   

  std::vector<ROOT::Math::XYZPointF> rhPosition_;
  std::vector<float> rhEnergy_;
  std::vector<float> rhEnergyFrac_;
  std::vector<float> rhTime_;
  std::vector<float> rhTimeError_;
  std::vector<int> rhCluster_; //in which layercluster 
  
  std::vector<std::vector<int>> vertices_[4];
  std::vector<float> time_[4];
  std::vector<float> timeError_[4];
  std::vector<float> regressed_energy_[4];
  std::vector<float> raw_energy_[4];
  std::vector<float> raw_em_energy_[4];
  std::vector<float> raw_pt_[4];
  std::vector<float> raw_em_pt_[4];
  std::vector<int>   temp_;
   
  std::string t_name[8]={"vertices","time","timeError","regressed_energy","raw_energy","raw_em_energy","raw_pt","raw_em_pt"};
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HGCAnalyzer::HGCAnalyzer(const edm::ParameterSet& iConfig)
   : detector_(iConfig.getParameter<std::string>("detector")) 
   {
   clusters_ = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("layer_clusters"));
   timelayercluster_ = consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("time_layerclusters"));
   //tracksterstrkem_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("tracksterstrkem"))
   trackstersem_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersem"));
   trackstershad_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstershad"));
   tracksterstrk_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("tracksterstrk"));
   trackstersmrg_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersmrg"));
   //auto caloParticles = iConfig.getParameter<edm::InputTag>("caloParticles");
   //caloParticles_ = consumes<std::vector<CaloParticle>>(caloParticles);
   auto recHitsEE = iConfig.getParameter<edm::InputTag>("recHitsEE");
   auto recHitsFH = iConfig.getParameter<edm::InputTag>("recHitsFH");
   auto recHitsBH = iConfig.getParameter<edm::InputTag>("recHitsBH");

   recHitsEE_ = consumes<HGCRecHitCollection>(recHitsEE);
   recHitsFH_ = consumes<HGCRecHitCollection>(recHitsFH);
   recHitsBH_ = consumes<HGCRecHitCollection>(recHitsBH);

   genParticles_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen_particles"));
}

HGCAnalyzer::~HGCAnalyzer() {
}

//
// member functions
// 

// ------------ method called for each event  ------------
void HGCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<std::vector<reco::CaloCluster>> clusterHandle;
  iEvent.getByToken(clusters_, clusterHandle); 
  auto const& clusters = *clusterHandle;
  //edm::Handle<std::vector<CaloParticle>> caloParticleHandle; 
  //iEvent.getByToken(caloParticles_, caloParticleHandle);
  //const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;

  edm::Handle<edm::ValueMap<std::pair<float, float>>> timelayerclusterHandle;
  iEvent.getByToken(timelayercluster_,timelayerclusterHandle);
  auto const& clustersTime = *timelayerclusterHandle;  

  //edm::Handle<std::vector<ticl::Trackster>> trackstertrkemHandle;
  //iEvent.getByToken(tracksterstrkem_, trackstertrkemHandle);
  //auto const& tracksterstrkem = *trackstertrkemHandle;

  edm::Handle<std::vector<ticl::Trackster>> tracksteremHandle;
  iEvent.getByToken(trackstersem_, tracksteremHandle);

  edm::Handle<std::vector<ticl::Trackster>> tracksterhadHandle;
  iEvent.getByToken(trackstershad_, tracksterhadHandle);

  edm::Handle<std::vector<ticl::Trackster>> trackstertrkHandle;
  iEvent.getByToken(tracksterstrk_, trackstertrkHandle);

  edm::Handle<std::vector<ticl::Trackster>> trackstermrgHandle;
  iEvent.getByToken(trackstersmrg_, trackstermrgHandle);

  edm::Handle<std::vector<ticl::Trackster>> tracksterHandle[4] = {tracksteremHandle, tracksterhadHandle, trackstertrkHandle, trackstermrgHandle};
  edm::Handle<std::vector<reco::GenParticle>> genParticleHandle;

  edm::Handle<HGCRecHitCollection> recHitHandleEE;
  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  edm::Handle<HGCRecHitCollection> recHitHandleBH;

  edm::ESHandle<CaloGeometry> geom;
  iSetup.get<CaloGeometryRecord>().get(geom);
  rhtools_.setGeometry(*geom);
  int layers_ = rhtools_.lastLayerBH();

  iEvent.getByToken(genParticles_, genParticleHandle);

  
  genPx_.clear();
  genPy_.clear();
  genPz_.clear();
  genEnergy_.clear();
  genEta_.clear();
  genPhi_.clear();
  genPt_.clear();
 
  PdgId_.clear();
  Particles_.clear();  

  auto const& genParticles = *genParticleHandle;
  for (auto const &gp : genParticles){
	
	genPx_.push_back(gp.px());
	genPy_.push_back(gp.py());
	genPz_.push_back(gp.pz());
	genPt_.push_back(gp.pt());
	genEta_.push_back(gp.eta());
        genPhi_.push_back(gp.phi());
        genEnergy_.push_back(gp.energy());
	
	ROOT::Math::PxPyPzEVector v;
	v.SetPxPyPzE(gp.px(),gp.py(),gp.pz(),gp.energy());
	Particles_.push_back(v); 
	PdgId_.push_back(gp.pdgId());
  }

  std::map<DetId, const HGCRecHit*> hitmap;

  iEvent.getByToken(recHitsEE_, recHitHandleEE);
  iEvent.getByToken(recHitsFH_, recHitHandleFH);
  iEvent.getByToken(recHitsBH_, recHitHandleBH);

  const auto& rechitsEE = *recHitHandleEE;
  const auto& rechitsFH = *recHitHandleFH;
  const auto& rechitsBH = *recHitHandleBH;
 
  for (unsigned int i = 0; i < rechitsEE.size(); ++i) {
         hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
        }
        for (unsigned int i = 0; i < rechitsFH.size(); ++i) {
         hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
        }
        for (unsigned int i = 0; i < rechitsBH.size(); ++i) {
         hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
        }
  
  int lc_idx=0;
  int rh_idx=0;
  lcEnergy_.clear();
  lcLayerN_.clear();
  lcTime_.clear();
  lcTimeError_.clear();
  lcId_.clear(); 
  lcPosition_.clear();
  lcHits_.clear();
  rhEnergy_.clear(); 
  rhEnergyFrac_.clear();
  rhTime_.clear();
  rhTimeError_.clear();
  rhPosition_.clear();
  rhCluster_.clear();

  for (auto const &lc : clusters){
        const auto firstHitDetId = lc.hitsAndFractions()[0].first;
	int layerId = rhtools_.getLayerWithOffset(firstHitDetId) + layers_ * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
	
        const auto& hits_and_fractions = lc.hitsAndFractions();

	for (const auto& it_h : hits_and_fractions) {
		const auto rh_detid = it_h.first;
	
		if (hitmap.find(rh_detid) == hitmap.end()) {
			return;
		}
		if ((rh_detid.det() != DetId::Forward) && (rh_detid.det() != DetId::HGCalEE) && (rh_detid.det() != DetId::HGCalHSi) &&
   		   (rh_detid.det() != DetId::HGCalHSc)) {
        		std::cout<<"Not HGCAL detector" << std::endl;
  			break;
		}

		auto global = rhtools_.getPosition(it_h.first);
		ROOT::Math::XYZPointF rh_p;
		rh_p.SetXYZ(global.x(),global.y(),global.z());
		rhPosition_.push_back(rh_p);
		rhTime_.push_back(hitmap[rh_detid]->time());
		rhTimeError_.push_back(hitmap[rh_detid]->timeError());
		rhEnergy_.push_back(hitmap[rh_detid]->energy());
		rhEnergyFrac_.push_back(it_h.second);
		rhCluster_.push_back(lc_idx);
		temp_.push_back(rh_idx);
			
		rh_idx++;
	}
        lcHits_.push_back(temp_);
	temp_.clear();//reuse container
	lcEnergy_.push_back(lc.energy());
        lcLayerN_.push_back(layerId);
 	lcId_.push_back(lc_idx);	
	ROOT::Math::XYZPointF lc_p;
	lc_p.SetXYZ(lc.x(),lc.y(),lc.z());	
        lcPosition_.push_back(lc_p);
	lcTime_.push_back(clustersTime.get(lc_idx).first);
	lcTimeError_.push_back(clustersTime.get(lc_idx).second);

	lc_idx++;
	}

   //tracksters
   for (int i=0;i<4;i++){
	vertices_[i].clear();
	time_[i].clear();
	timeError_[i].clear();
	regressed_energy_[i].clear();
	raw_energy_[i].clear();
	raw_em_energy_[i].clear();
	raw_pt_[i].clear();
	raw_em_pt_[i].clear();
	auto const& tracksters = *tracksterHandle[i];
   	for (auto const &tkr : tracksters){
   		if (!tkr.vertices().empty()) {
			for (auto i : tkr.vertices()){temp_.push_back(i);}
			vertices_[i].push_back(temp_);
			time_[i].push_back(tkr.time());
			timeError_[i].push_back(tkr.timeError());
			regressed_energy_[i].push_back(tkr.regressed_energy());
			raw_energy_[i].push_back(tkr.raw_energy());
			raw_em_energy_[i].push_back(tkr.raw_em_energy());
			raw_pt_[i].push_back(tkr.raw_pt());
			raw_em_pt_[i].push_back(tkr.raw_em_pt());
			temp_.clear();
		}
   
   	}
   }

   tree_->Fill();
   //std::cout<<"finished one event"<<std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void HGCAnalyzer::beginJob() {
 	tree_ = fs_->make<TTree>("tree", "TICL objects");
 	tree_->Branch("genParticles",&Particles_);
		
	tree_->Branch("genPdgId", &PdgId_);
	tree_->Branch("lcEnergy", &lcEnergy_);
  	tree_->Branch("lcLayerN", &lcLayerN_);
  	tree_->Branch("lcId", &lcId_);
        tree_->Branch("lcPosition", &lcPosition_);
	tree_->Branch("lcTime", &lcTime_);
	tree_->Branch("lcTimeError", &lcTimeError_);
	tree_->Branch("lcHits", &lcHits_);
	tree_->Branch("rhTime", &rhTime_);
	tree_->Branch("rhTimeError", &rhTimeError_);
	tree_->Branch("rhEnergy", &rhEnergy_);
        tree_->Branch("rhEnergyFrac", &rhEnergyFrac_);
	tree_->Branch("rhPosition", &rhPosition_);	
 	tree_->Branch("rhCluster", &rhCluster_); 
	
	for(int i=0;i<4;i++){
	std::string s_temp[8]{};
	for(int j=0;j<8;j++){
		s_temp[j]=s_tracksters[i]+"_"+t_name[j];
		}
	
	tree_->Branch(s_temp[0].c_str(), &vertices_[i]);
	tree_->Branch(s_temp[1].c_str(), &time_[i]);
        tree_->Branch(s_temp[2].c_str(), &timeError_[i]); 
        tree_->Branch(s_temp[3].c_str(), &regressed_energy_[i]);
        tree_->Branch(s_temp[4].c_str(), &raw_energy_[i]);
        tree_->Branch(s_temp[5].c_str(), &raw_em_energy_[i]);
        tree_->Branch(s_temp[6].c_str(), &raw_pt_[i]);
        tree_->Branch(s_temp[7].c_str(), &raw_em_pt_[i]);
	
	}

}

// ------------ method called once each job just after ending the event loop  ------------
void HGCAnalyzer::endJob() {
	std::cout<<"finished"<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("time_layerclusters", edm::InputTag("hgcalLayerClusters","timeLayerCluster"));
  desc.add<edm::InputTag>("trackstersem", edm::InputTag("ticlTrackstersEM"));
  desc.add<edm::InputTag>("trackstershad", edm::InputTag("ticlTrackstersHAD"));
  desc.add<edm::InputTag>("trackstersmrg", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("tracksterstrk", edm::InputTag("ticlTrackstersTrk"));
  desc.add<edm::InputTag>("gen_particles", edm::InputTag("genParticles"));
  //desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("recHitsEE", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  desc.add<edm::InputTag>("recHitsFH", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  desc.add<edm::InputTag>("recHitsBH", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
  desc.add<std::string>("detector", "HGCAL");
  descriptions.addDefault(desc);
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCAnalyzer);
