// JetValidationMetrics.cc

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"

#include <cmath>
#include <string>

class JetValidationMetrics : public edm::EDAnalyzer {
public:
  explicit JetValidationMetrics(const edm::ParameterSet& iConfig);
  ~JetValidationMetrics() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;

private:
  // inputs
  const edm::InputTag jetsTag_;
  const edm::InputTag genJetsTag_;

  const double recoJetPtThreshold_;
  const double matchGenPtThreshold_;
  const double rThreshold_;
  const double absEtaMax_;

  edm::EDGetTokenT<edm::View<reco::Jet>> jetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;

  // output
  TTree* outTree_ = nullptr;

  // accumulators (over whole job)
  uint64_t nEvents_ = 0;
  uint64_t nGen_ = 0;
  uint64_t nMatched_ = 0;

  double sumResp_ = 0.0;
  double sumResp2_ = 0.0;

  double sumDR_ = 0.0;
  double sumDR2_ = 0.0;

  // tree branches (written once in endJob)
  uint64_t b_nEvents_ = 0;
  uint64_t b_nGen_ = 0;
  uint64_t b_nMatched_ = 0;

  double b_matchEff_ = 0.0;

  double b_meanResp_ = 0.0;
  double b_rmsResp_ = 0.0;

  double b_meanDR_ = 0.0;
  double b_rmsDR_ = 0.0;
};

JetValidationMetrics::JetValidationMetrics(const edm::ParameterSet& iConfig)
    : jetsTag_(iConfig.getParameter<edm::InputTag>("jets")),
      genJetsTag_(iConfig.getParameter<edm::InputTag>("genjets")),
      recoJetPtThreshold_(iConfig.getParameter<double>("recoJetPtThreshold")),
      matchGenPtThreshold_(iConfig.getParameter<double>("matchGenPtThreshold")),
      rThreshold_(iConfig.getParameter<double>("RThreshold")),
      absEtaMax_(iConfig.getParameter<double>("absEtaMax")) {
  jetsToken_ = consumes<edm::View<reco::Jet>>(jetsTag_);
  genJetsToken_ = consumes<reco::GenJetCollection>(genJetsTag_);

  // --- SimpleValidation-style ROOT layout:
  // create a subdirectory named after this module label, then put TTree "output" inside it.
  edm::Service<TFileService> fs;
  const std::string label = moduleDescription().moduleLabel();
  TFileDirectory dir = fs->mkdir(label);

  outTree_ = dir.make<TTree>("output", "output");

  outTree_->Branch("nEvents", &b_nEvents_, "nEvents/l");
  outTree_->Branch("nGen", &b_nGen_, "nGen/l");
  outTree_->Branch("nMatched", &b_nMatched_, "nMatched/l");

  outTree_->Branch("matchEff", &b_matchEff_, "matchEff/D");

  outTree_->Branch("meanResp", &b_meanResp_, "meanResp/D");
  outTree_->Branch("rmsResp", &b_rmsResp_, "rmsResp/D");

  outTree_->Branch("meanDR", &b_meanDR_, "meanDR/D");
  outTree_->Branch("rmsDR", &b_rmsDR_, "rmsDR/D");
}

void JetValidationMetrics::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  ++nEvents_;

  edm::Handle<edm::View<reco::Jet>> jetsH;
  edm::Handle<reco::GenJetCollection> genJetsH;

  iEvent.getByToken(jetsToken_, jetsH);
  iEvent.getByToken(genJetsToken_, genJetsH);

  if (!jetsH.isValid() || !genJetsH.isValid()) return;
  if (jetsH->empty() || genJetsH->empty()) return;

  // loop over gen jets, match to closest reco jet
  for (const auto& gjet : *genJetsH) {
    if (gjet.pt() < matchGenPtThreshold_) continue;
    if (std::abs(gjet.eta()) > absEtaMax_) continue;

    ++nGen_;

    int bestIdx = -1;
    double bestDR2 = 1e99;

    for (size_t i = 0; i < jetsH->size(); ++i) {
      const auto& rjet = (*jetsH)[i];
      if (rjet.pt() < recoJetPtThreshold_) continue;

      const double dR2 = reco::deltaR2(gjet.eta(), gjet.phi(), rjet.eta(), rjet.phi());
      if (dR2 < bestDR2) {
        bestDR2 = dR2;
        bestIdx = static_cast<int>(i);
      }
    }

    if (bestIdx < 0) continue;
    if (bestDR2 >= rThreshold_ * rThreshold_) continue;
    if (gjet.pt() == 0.0) continue;

    ++nMatched_;

    const auto& rjet = (*jetsH)[bestIdx];
    const double resp = rjet.pt() / gjet.pt();
    const double dR = std::sqrt(bestDR2);

    sumResp_ += resp;
    sumResp2_ += resp * resp;

    sumDR_ += dR;
    sumDR2_ += dR * dR;
  }
}

void JetValidationMetrics::endJob() {
  // finalize scalars
  b_nEvents_ = nEvents_;
  b_nGen_ = nGen_;
  b_nMatched_ = nMatched_;

  b_matchEff_ = (nGen_ > 0) ? static_cast<double>(nMatched_) / static_cast<double>(nGen_) : 0.0;

  if (nMatched_ > 0) {
    const double invN = 1.0 / static_cast<double>(nMatched_);

    const double meanResp = sumResp_ * invN;
    const double varResp = (sumResp2_ * invN) - (meanResp * meanResp);
    b_meanResp_ = meanResp;
    b_rmsResp_ = (varResp > 0.0) ? std::sqrt(varResp) : 0.0;

    const double meanDR = sumDR_ * invN;
    const double varDR = (sumDR2_ * invN) - (meanDR * meanDR);
    b_meanDR_ = meanDR;
    b_rmsDR_ = (varDR > 0.0) ? std::sqrt(varDR) : 0.0;
  } else {
    b_meanResp_ = 0.0;
    b_rmsResp_ = 0.0;
    b_meanDR_ = 0.0;
    b_rmsDR_ = 0.0;
  }

  // write exactly one row (SimpleValidation-style: scalar summary for this config)
  outTree_->Fill();
}

void JetValidationMetrics::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets", edm::InputTag("hltAK4PFPuppiJets"));
  desc.add<edm::InputTag>("genjets", edm::InputTag("ak4GenJetsNoNu"));

  desc.add<double>("recoJetPtThreshold", 30.0);
  desc.add<double>("matchGenPtThreshold", 20.0);
  desc.add<double>("RThreshold", 0.4);
  desc.add<double>("absEtaMax", 6.0);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(JetValidationMetrics);
