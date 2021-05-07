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
#include <array>

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

#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
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

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticles_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterstrkem_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersem_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterstrk_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstershad_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersmrg_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;
  edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> timelayercluster_;  
  edm::EDGetTokenT<hgcal::LayerClusterToCaloParticleAssociator> LCAssocByEnergyScoreProducer_;
  hgcal::RecHitTools rhtools_;
  // ----------member data ---------------------------
  edm::Service<TFileService> fs_;
  TTree* tree_;
  
  std::string s_tracksters[5] = {"tracksterEM", "tracksterHAD", "tracksterMerge", "tracksterTrk", "tracksterTrkEM"};

  edm::EDGetTokenT<std::unordered_map<DetId, const HGCRecHit*>> hitMap_;

  std::vector<int> gpdgId_;
  std::vector<ROOT::Math::PxPyPzEVector> gp_;
  std::vector<ROOT::Math::XYZPointF> gpPosition_;

  std::vector<ROOT::Math::XYZPointF> lcPosition_;
  std::vector<float> lcEnergy_;
  std::vector<uint32_t> lcLayer_;
  std::vector<float> lcTime_;
  std::vector<float> lcTimeError_;
  std::vector<int> lcId_;
  std::vector<std::vector<int>> lcHits_;   

  std::vector<ROOT::Math::XYZPointF> rhPosition_;
  std::vector<float> rhEnergy_;
  std::vector<float> rhTime_;
  std::vector<float> rhTimeError_;
  std::vector<uint32_t> rhLayer_;  

  std::vector<std::vector<unsigned int>> vertices_[5];
  std::vector<std::vector<unsigned int>> vertex_multiplicity_[5];
  std::vector<int> seedIndex_[5];
  std::vector<float> time_[5];
  std::vector<float> timeError_[5];
  std::vector<float> regressed_energy_[5];
  std::vector<float> raw_energy_[5];
  std::vector<float> raw_em_energy_[5];
  std::vector<float> raw_pt_[5];
  std::vector<float> raw_em_pt_[5];
  std::vector<ROOT::Math::XYZVector> barycenter_[5];
  std::vector<std::vector<float>> sigmas_[5];
  std::vector<std::vector<float>> sigmasPCA_[5];
  std::vector<std::vector<float>> id_probabilities_[5];
  std::vector<float> sig_tmp,sigPCA_tmp,iP_tmp;

  std::vector<int>   temp_;
  std::vector<int>   temp2_; 
  std::vector<float> temp3_;
  std::vector<float> temp4_;
 
  std::string t_name[14]={"vertices","vertex_multiplicity","seedIndex","time","timeError","regressed_energy","raw_energy","raw_em_energy","raw_pt","raw_em_pt","barycenter","sigmas","sigmasPCA","id_probabilities"};

  std::vector<int> cpdgId_;
  std::vector<ROOT::Math::PxPyPzEVector> cp_;
  std::vector<std::vector<int>> cpSC_;
  std::vector<int> cpG4T0evt_;
  std::vector<int> cpG4T0bx_;
  std::vector<ROOT::Math::XYZPointF> cpOrigin_;  

  std::vector<float> scEnergy_;
  std::vector<float> scSimEnergy_;
  std::vector<std::vector<int>> scHits_;
  std::vector<std::vector<float>> scHitsEnergyFrac_;

  std::vector<DetId> rhcontainer_;  
  std::vector<DetId>::iterator it;
  
  
  std::vector<std::vector<float>> lc2cpScore_;
  std::vector<std::vector<int>> lc2cpId_;
  std::vector<std::vector<float>> lc2cpEnergy_; //for testing idx

  std::vector<std::vector<float>> cp2lcScore_;
  std::vector<std::vector<int>> cp2lcId_;
  std::vector<std::vector<float>> cp2lcEnergy_; //for testing idx
  
  int setZside_;
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
HGCAnalyzer::HGCAnalyzer(const edm::ParameterSet& iConfig):
   setZside_(iConfig.getUntrackedParameter<int>("setZside",1)) //set default, 1->zplus, -1->zminus, 0->both
   {
   clusters_ = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("layer_clusters"));
   timelayercluster_ = consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("time_layerclusters"));
   tracksterstrkem_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("tracksterstrkem"));
   trackstersem_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersem"));
   trackstershad_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstershad"));
   tracksterstrk_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("tracksterstrk"));
   trackstersmrg_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersmrg"));
   caloParticles_ = consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("caloParticles"));
   hitMap_ = consumes<std::unordered_map<DetId, const HGCRecHit*>>(iConfig.getParameter<edm::InputTag>("hitMapTag"));

   genParticles_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen_particles"));

   LCAssocByEnergyScoreProducer_ = consumes<hgcal::LayerClusterToCaloParticleAssociator>(iConfig.getParameter<edm::InputTag>("lcAssocByEnergyScoreProducer"));   
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
 
  edm::Handle<std::vector<CaloParticle>> caloParticleHandle; 
  iEvent.getByToken(caloParticles_, caloParticleHandle);
  auto const& caloParticles = *caloParticleHandle;

  edm::Handle<edm::ValueMap<std::pair<float, float>>> timelayerclusterHandle;
  iEvent.getByToken(timelayercluster_,timelayerclusterHandle);
  auto const& clustersTime = *timelayerclusterHandle;  

  edm::Handle<std::vector<ticl::Trackster>> trackstertrkemHandle;
  iEvent.getByToken(tracksterstrkem_, trackstertrkemHandle);

  edm::Handle<std::vector<ticl::Trackster>> tracksteremHandle;
  iEvent.getByToken(trackstersem_, tracksteremHandle);

  edm::Handle<std::vector<ticl::Trackster>> tracksterhadHandle;
  iEvent.getByToken(trackstershad_, tracksterhadHandle);

  edm::Handle<std::vector<ticl::Trackster>> trackstertrkHandle;
  iEvent.getByToken(tracksterstrk_, trackstertrkHandle);

  edm::Handle<std::vector<ticl::Trackster>> trackstermrgHandle;
  iEvent.getByToken(trackstersmrg_, trackstermrgHandle);

  edm::Handle<std::vector<ticl::Trackster>> tracksterHandle[5] = {tracksteremHandle, tracksterhadHandle, trackstermrgHandle, trackstertrkHandle, trackstertrkemHandle};
  edm::Handle<std::vector<reco::GenParticle>> genParticleHandle;

  edm::Handle<std::unordered_map<DetId, const HGCRecHit*>> hitMapHandle;
  iEvent.getByToken(hitMap_, hitMapHandle);
  const auto hitmap = *hitMapHandle;
  edm::Handle<hgcal::LayerClusterToCaloParticleAssociator> LCAssocByEnergyScoreHandle;
  iEvent.getByToken(LCAssocByEnergyScoreProducer_, LCAssocByEnergyScoreHandle);

  hgcal::RecoToSimCollection cpsInLayerClusterMap = LCAssocByEnergyScoreHandle->associateRecoToSim(clusterHandle, caloParticleHandle);
  hgcal::SimToRecoCollection cPOnLayerMap = LCAssocByEnergyScoreHandle->associateSimToReco(clusterHandle, caloParticleHandle);

  edm::ESHandle<CaloGeometry> geom;
  iSetup.get<CaloGeometryRecord>().get(geom);
  rhtools_.setGeometry(*geom);
  int layers_ = rhtools_.lastLayerBH();

  iEvent.getByToken(genParticles_, genParticleHandle);

  //genParticles
  gpdgId_.clear();
  gp_.clear();
  gpPosition_.clear();

  auto const& genParticles = *genParticleHandle;
  for (auto const &gp : genParticles){
	ROOT::Math::PxPyPzEVector p4;
	p4.SetPxPyPzE(gp.px(),gp.py(),gp.pz(),gp.energy());
	gp_.push_back(p4); 
	gpdgId_.push_back(gp.pdgId());
	ROOT::Math::XYZPointF v;
        v.SetXYZ(gp.vx(),gp.vy(),gp.vz());	
  	gpPosition_.push_back(v);
  }
 
  int rh_idx=0;
  lcEnergy_.clear();
  lcLayer_.clear();
  lcTime_.clear();
  lcTimeError_.clear();
  lcId_.clear(); 
  lcPosition_.clear();
  lcHits_.clear();
  rhEnergy_.clear(); 
  rhTime_.clear();
  rhTimeError_.clear();
  rhPosition_.clear();
  rhcontainer_.clear();
  rhLayer_.clear();
  lc2cpScore_.clear();
  lc2cpId_.clear();
  lc2cpEnergy_.clear();
  for (auto const &lc : clusters){
        int lc_idx=&lc-&clusters[0];
	if (setZside_*lc.z()<0) continue;
	const auto firstHitDetId = lc.hitsAndFractions()[0].first;
	int layerId = rhtools_.getLayerWithOffset(firstHitDetId) + layers_ * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
        const auto& hits_and_fractions = lc.hitsAndFractions();
	
	for (const auto& it_h : hits_and_fractions) {
		const auto rh_detid = it_h.first;
		if (hitmap.find(rh_detid) == hitmap.end()) {
			continue;
		}
		if ((rh_detid.det() != DetId::Forward) && (rh_detid.det() != DetId::HGCalEE) && (rh_detid.det() != DetId::HGCalHSi) &&
   		   (rh_detid.det() != DetId::HGCalHSc)) {
        		std::cout<<"Not HGCAL detector" << std::endl;
  			continue;
		}
		auto global = rhtools_.getPosition(it_h.first);
		ROOT::Math::XYZPointF rh_p;
		rh_p.SetXYZ(global.x(),global.y(),global.z());
		rhPosition_.push_back(rh_p);
		rhTime_.push_back(hitmap.at(rh_detid)->time());
		rhTimeError_.push_back(hitmap.at(rh_detid)->timeError());
		rhEnergy_.push_back(hitmap.at(rh_detid)->energy());
		int rhl=rhtools_.getLayerWithOffset(rh_detid) + layers_ * ((rhtools_.zside(rh_detid) + 1) >> 1) - 1;
		rhLayer_.push_back(rhl);		
		temp_.push_back(rh_idx);
		rhcontainer_.push_back(rh_detid);
		rh_idx++;
	}
        lcHits_.push_back(temp_);
	temp_.clear();//reuse container
	lcEnergy_.push_back(lc.energy());
        lcLayer_.push_back(layerId);
 	lcId_.push_back(lc_idx);	
	ROOT::Math::XYZPointF lc_p;
	lc_p.SetXYZ(lc.x(),lc.y(),lc.z());	
        lcPosition_.push_back(lc_p);
	lcTime_.push_back(clustersTime.get(lc_idx).first);
	lcTimeError_.push_back(clustersTime.get(lc_idx).second);

	
	//metrics
	//std::cout<<"lc_idx "<<lc_idx<<std::endl;
	const edm::Ref<std::vector<reco::CaloCluster>> lcRef(clusterHandle, lc_idx);
        const auto& cpsIt = cpsInLayerClusterMap.find(lcRef);
        if (cpsIt == cpsInLayerClusterMap.end()){
                lc2cpScore_.push_back(temp3_);
		lc2cpId_.push_back(temp2_);
		lc2cpEnergy_.push_back(temp4_);
		}
        else{
		const auto& cps = cpsIt->val;
		for (const auto& cpPair : cps) {
				if (cpPair.first->g4Tracks()[0].eventId().event()!=0 or cpPair.first->g4Tracks()[0].eventId().bunchCrossing()!=0) continue;
	
				auto const& cp_linked = 
				std::find_if(std::begin(cPOnLayerMap[cpPair.first]),std::end(cPOnLayerMap[cpPair.first]),
                        	[&lcRef](const std::pair<edm::Ref<reco::CaloClusterCollection>, std::pair<float, float>>& p) {
                        	 return p.first == lcRef;
                       		});
				
      				if (cp_linked == cPOnLayerMap[cpPair.first].end()) continue;
				
				temp3_.push_back(cpPair.second);
                                temp2_.push_back(cpPair.first.index());
				temp4_.push_back(cp_linked->second.first);
				
				}
		lc2cpScore_.push_back(temp3_);
		lc2cpId_.push_back(temp2_);	
		lc2cpEnergy_.push_back(temp4_);
		temp3_.clear();
		temp2_.clear();
		temp4_.clear();
		}
   }
   //caloparticles
   int sc_idx=0;
   cp_.clear();
   cpdgId_.clear();
   cpSC_.clear();
   cpOrigin_.clear();
   cpG4T0evt_.clear();
   cpG4T0bx_.clear();
   scEnergy_.clear();
   scHits_.clear();
   scHitsEnergyFrac_.clear();
   cp2lcScore_.clear();
   cp2lcId_.clear();
   cp2lcEnergy_.clear();
   for (const auto& cp : caloParticles) {
	if (cp.g4Tracks()[0].eventId().event()!=0 or cp.g4Tracks()[0].eventId().bunchCrossing()!=0) continue;  //saving caloparticles from hard scattering only 	
	if (setZside_*cp.g4Tracks()[0].trackerSurfacePosition().Z()<0) continue;
	int cp_idx=&cp-&caloParticles[0];
        cpG4T0evt_.push_back(cp.g4Tracks()[0].eventId().event());
	cpG4T0bx_.push_back(cp.g4Tracks()[0].eventId().bunchCrossing());
	ROOT::Math::PxPyPzEVector p4;
        p4.SetPxPyPzE(cp.px(),cp.py(),cp.pz(),cp.energy());
	cp_.push_back(p4);
        cpdgId_.push_back(cp.pdgId());
	//initial g4track vertex
	ROOT::Math::XYZPointF cp_o;
	cp_o.SetXYZ(cp.g4Tracks()[0].trackerSurfacePosition().X(),cp.g4Tracks()[0].trackerSurfacePosition().Y(),cp.g4Tracks()[0].trackerSurfacePosition().Z());
	cpOrigin_.push_back(cp_o);
	const SimClusterRefVector& simClusterRefVector = cp.simClusters();
	for (const auto& sc : simClusterRefVector) {
		const SimCluster& simCluster = (*(sc));
		const auto& hits_and_fractions = simCluster.hits_and_fractions();	
		for (const auto& it_h : hits_and_fractions) {
			DetId rh_detid = (it_h.first);
			if (!hitmap.count(rh_detid)) continue;
                	if ((rh_detid.det() != DetId::Forward) && (rh_detid.det() != DetId::HGCalEE) && (rh_detid.det() != DetId::HGCalHSi) &&
                   	(rh_detid.det() != DetId::HGCalHSc)) {
                        	std::cout<<"Not HGCAL detector" << std::endl;
                        	continue;
                	}

			it = std::find (rhcontainer_.begin(), rhcontainer_.end(), rh_detid);
			if (it != rhcontainer_.end()) {
				temp_.push_back(it - rhcontainer_.begin());
				temp3_.push_back(it_h.second);
			}
			else {
				auto global = rhtools_.getPosition(it_h.first);
                		ROOT::Math::XYZPointF rh_p;
                		rh_p.SetXYZ(global.x(),global.y(),global.z());
                		rhPosition_.push_back(rh_p);
                		rhTime_.push_back(hitmap.at(rh_detid)->time());
                		rhTimeError_.push_back(hitmap.at(rh_detid)->timeError());
                		rhEnergy_.push_back(hitmap.at(rh_detid)->energy());
                		int rhl=rhtools_.getLayerWithOffset(rh_detid) + layers_ * ((rhtools_.zside(rh_detid) + 1) >> 1) - 1;
                		rhLayer_.push_back(rhl);
				temp_.push_back(rh_idx);
              			temp3_.push_back(it_h.second);
				rh_idx++;
			}
		}
		scEnergy_.push_back(simCluster.energy());
		scSimEnergy_.push_back(simCluster.simEnergy());
		scHits_.push_back(temp_);
		scHitsEnergyFrac_.push_back(temp3_);
		temp_.clear();
		temp3_.clear();
		temp2_.push_back(sc_idx);
	}
	cpSC_.push_back(temp2_);
	temp2_.clear();
   
        
	const edm::Ref<CaloParticleCollection> cpRef(caloParticleHandle, cp_idx);
    	const auto& lcsIt = cPOnLayerMap.find(cpRef);
	if (lcsIt == cPOnLayerMap.end()){
                cp2lcScore_.push_back(temp4_);
                cp2lcId_.push_back(temp2_);
		cp2lcEnergy_.push_back(temp3_);
		}
	else {
        	const auto& lcs = lcsIt->val;
        	for (const auto& lcPair : lcs) {
                	temp4_.push_back(lcPair.second.second);
                        temp2_.push_back(lcPair.first.index());
			temp3_.push_back(lcPair.second.first);	
                	}
        	cp2lcScore_.push_back(temp4_);
        	cp2lcId_.push_back(temp2_);
		cp2lcEnergy_.push_back(temp3_);
        	temp4_.clear();
        	temp2_.clear();
		temp3_.clear();
	}
   }	

   //tracksters
   for (int i=0;i<5;i++){
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
			if (tkr.barycenter().z()*setZside_<0)continue;
			vertices_[i].push_back(tkr.vertices());
			std::vector<unsigned int> temp_vec;
			for (auto& vm:tkr.vertex_multiplicity()) {
				temp_vec.push_back(unsigned(vm));
		        }
			vertex_multiplicity_[i].push_back(temp_vec);
			seedIndex_[i].push_back(tkr.seedIndex());
			time_[i].push_back(tkr.time());
			timeError_[i].push_back(tkr.timeError());
			regressed_energy_[i].push_back(tkr.regressed_energy());
			raw_energy_[i].push_back(tkr.raw_energy());
			raw_em_energy_[i].push_back(tkr.raw_em_energy());
			raw_pt_[i].push_back(tkr.raw_pt());
			raw_em_pt_[i].push_back(tkr.raw_em_pt());
			barycenter_[i].push_back(tkr.barycenter());
			for (int j=0;j<3;j++){
				sig_tmp.push_back(tkr.sigmas()[j]);
				sigPCA_tmp.push_back(tkr.sigmasPCA()[j]); 				
			}
			for (int j=0;j<8;j++)iP_tmp.push_back(tkr.id_probabilities()[j]);
			sigmas_[i].push_back(sig_tmp);
			sigmasPCA_[i].push_back(sigPCA_tmp);
			id_probabilities_[i].push_back(iP_tmp);
			sig_tmp.clear();
			sigPCA_tmp.clear();
			iP_tmp.clear();
		}
   	}
   }

   tree_->Fill();
   //std::cout<<"finished one event"<<std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void HGCAnalyzer::beginJob() {
 	tree_ = fs_->make<TTree>("tree", "TICL objects");
 	//tree_->SetAutoSave(0);
	
	tree_->Branch("genParticle", &gp_);
        tree_->Branch("genParticlePosition", &gpPosition_);		
	tree_->Branch("genPdgId", &gpdgId_);

	tree_->Branch("lcEnergy", &lcEnergy_);
  	tree_->Branch("lcLayer", &lcLayer_);
  	tree_->Branch("lcId", &lcId_);
        tree_->Branch("lcPosition", &lcPosition_);
	tree_->Branch("lcTime", &lcTime_);
	tree_->Branch("lcTimeError", &lcTimeError_);
	tree_->Branch("lcHits", &lcHits_);
	tree_->Branch("lc2cpScore", &lc2cpScore_);
        tree_->Branch("lc2cpId", &lc2cpId_);
        tree_->Branch("lc2cpEnergy", &lc2cpEnergy_);

	tree_->Branch("rhTime", &rhTime_);
	tree_->Branch("rhTimeError", &rhTimeError_);
	tree_->Branch("rhEnergy", &rhEnergy_);
	tree_->Branch("rhPosition", &rhPosition_);	
 	tree_->Branch("rhLayer", &rhLayer_); 
	
	for(int i=0;i<5;i++){
	std::string s_temp[14]{};
	for(int j=0;j<14;j++){
		s_temp[j]=s_tracksters[i]+"_"+t_name[j];
		}
		
	tree_->Branch(s_temp[0].c_str(), &vertices_[i]);
	tree_->Branch(s_temp[1].c_str(), &vertex_multiplicity_[i]);
        tree_->Branch(s_temp[2].c_str(), &seedIndex_[i]);
	tree_->Branch(s_temp[3].c_str(), &time_[i]);
        tree_->Branch(s_temp[4].c_str(), &timeError_[i]); 
        tree_->Branch(s_temp[5].c_str(), &regressed_energy_[i]);
        tree_->Branch(s_temp[6].c_str(), &raw_energy_[i]);
        tree_->Branch(s_temp[7].c_str(), &raw_em_energy_[i]);
        tree_->Branch(s_temp[8].c_str(), &raw_pt_[i]);
        tree_->Branch(s_temp[9].c_str(), &raw_em_pt_[i]);
	tree_->Branch(s_temp[10].c_str(), &barycenter_[i]);
	tree_->Branch(s_temp[11].c_str(), &sigmas_[i]);
	tree_->Branch(s_temp[12].c_str(), &sigmasPCA_[i]);
	tree_->Branch(s_temp[13].c_str(), &id_probabilities_[i]);
        }
	tree_->Branch("cpdgId", &cpdgId_);
	tree_->Branch("cp", &cp_);
	tree_->Branch("cpSC", &cpSC_);
	tree_->Branch("cpOrigin", &cpOrigin_);
	tree_->Branch("cpG4T0evt", &cpG4T0evt_);
	tree_->Branch("cpG4T0bx", &cpG4T0bx_);
	tree_->Branch("cp2lcScore", &cp2lcScore_);
	tree_->Branch("cp2lcId", &cp2lcId_);
	tree_->Branch("cp2lcEnergy", &cp2lcEnergy_);
	
	tree_->Branch("scEnergy", &scEnergy_);
	tree_->Branch("scSimEnergy", &scSimEnergy_);
	tree_->Branch("scHits", &scHits_);
	tree_->Branch("scHitsEnergyFrac", &scHitsEnergyFrac_);	


}

// ------------ method called once each job just after ending the event loop  ------------
void HGCAnalyzer::endJob() {
	std::cout<<"finished"<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  desc.setAllowAnything();
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("time_layerclusters", edm::InputTag("hgcalLayerClusters","timeLayerCluster"));
  desc.add<edm::InputTag>("trackstersem", edm::InputTag("ticlTrackstersEM"));
  desc.add<edm::InputTag>("trackstershad", edm::InputTag("ticlTrackstersHAD"));
  desc.add<edm::InputTag>("trackstersmrg", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("tracksterstrk", edm::InputTag("ticlTrackstersTrk"));
  desc.add<edm::InputTag>("tracksterstrkem", edm::InputTag("ticlTrackstersTrkEM"));
  desc.add<edm::InputTag>("gen_particles", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("hitMapTag", edm::InputTag("hgcalRecHitMapProducer"));
  desc.add<edm::InputTag>("lcAssocByEnergyScoreProducer", edm::InputTag("layerClusterAssociatorByEnergyScore"));
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCAnalyzer);
