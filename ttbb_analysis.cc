// -*- C++ -*-
//default rivet-mkanalysis
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
//additional
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/JetAlg.hh"

using namespace std;

namespace Rivet {

  
  /// @brief Add a short analysis description here
  class ttbb_analysis : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ttbb_analysis);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection: 
      // all final-state particles within 
      // the given eta acceptance
      // const FinalState fs(Cuts::abseta < 5);

      Cut eta_full = (Cuts::abseta < 5.0);
      // Lepton cuts
      Cut lep_cuts27 = (Cuts::abseta < 2.5) && (Cuts::pT >= 27*GeV);
      // All final state particles
      FinalState fs(eta_full);

      // the final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      //FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      //declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      //PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      //PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      // Get photons to dress leptons
      PromptFinalState photons(eta_full && Cuts::abspid == PID::PHOTON, true);
      // Projection to find the electrons
      PromptFinalState electrons(eta_full && Cuts::abspid == PID::ELECTRON, true);
      // Projection to find the muons
      PromptFinalState muons(eta_full && Cuts::abspid == PID::MUON, true);

      // dress the prompt bare leptons with prompt photons within dR < 0.1
      // apply some fiducial cuts on the dressed leptons
      //Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      //DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      //declare(dressed_leps, "leptons");
      DressedLeptons dressedelectrons27(photons, electrons, 0.1, lep_cuts27, true);
      DressedLeptons dressedmuons27(photons, muons, 0.1, lep_cuts27, true);
      declare(dressedelectrons27, "elecs");
      declare(dressedmuons27, "muons");

      // From here on we are just setting up the jet clustering
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);

      PromptFinalState jet_photons(eta_full && Cuts::abspid == PID::PHOTON, false);
      DressedLeptons all_dressed_electrons(jet_photons, electrons, 0.1, eta_full, true);
      DressedLeptons all_dressed_muons(jet_photons, muons, 0.1, eta_full, true);

      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(all_dressed_electrons);
      vfs.addVetoOnThisFinalState(all_dressed_muons);
      vfs.addVetoOnThisFinalState(neutrinos);

      // FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::DECAY, JetAlg::Invisibles::DECAY);
      declare(jets, "jets");

      // missing momentum
      // declare(MissingMomentum(fs), "MET");

      vector<double> multiplicity_bins    = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};
      vector<double> pt_bins              = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 500};
      vector<double> ht_bins              = {0, 50, 100, 135, 170, 205, 240, 275, 310, 345, 380, 415, 450, 485, 520, 555, 590, 625, 660, 695, 730, 765, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400, 1500};
      vector<double> m_leading_bins       = {0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 360, 375, 390, 405, 420, 435, 450, 480, 510, 540, 570, 600};
      vector<double> m_closest_bins       = {0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 270, 300, 330, 360, 390, 450, 510};
      vector<double> dr_bins              = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0};

      // Book histograms
      // specify custom binning
      // 3b geq4j category 
      book(_h["N_Jets_3b_geq4j_ljets"],              "N_Jets_3b_geq4j_ljets"              , multiplicity_bins);
      book(_h["N_b_Jets_3b_geq4j_ljets"],            "N_b_Jets_3b_geq4j_ljets"            , multiplicity_bins);
      book(_h["all_bjets_pt_3b_geq4j_ljets"],        "all_bjets_pt_3b_geq4j_ljets"        , pt_bins);

      book(_h["ht_bjets_3b_geq4j_ljets"],            "ht_bjets_3b_geq4j_ljets"            , ht_bins); 
      book(_h["ht_3b_geq4j_ljets"],                  "ht_3b_geq4j_ljets"                  , ht_bins); 
      book(_h["ht_had_3b_geq4j_ljets"],              "ht_had_3b_geq4j_ljets"              , ht_bins);

      book(_h["lead_bjet_pt_3b_geq4j_ljets"],        "lead_bjet_pt_3b_geq4j_ljets"        , pt_bins);
      book(_h["sublead_bjet_pt_3b_geq4j_ljets"],     "sublead_bjet_pt_3b_geq4j_ljets"     , pt_bins);
      book(_h["third_bjet_pt_3b_geq4j_ljets"],       "third_bjet_pt_3b_geq4j_ljets"       , pt_bins);
      book(_h["fourth_bjet_pt_3b_geq4j_ljets"],      "fourth_bjet_pt_3b_geq4j_ljets"      , pt_bins);

      book(_h["m_bb_leading_3b_geq4j_ljets"],        "m_bb_leading_3b_geq4j_ljets"        , m_leading_bins);
      book(_h["pt_bb_leading_3b_geq4j_ljets"],       "pt_bb_leading_3b_geq4j_ljets"       , pt_bins);
      book(_h["dR_bb_leading_3b_geq4j_ljets"],       "dR_bb_leading_3b_geq4j_ljets"       , dr_bins);

      book(_h["m_bb_closest_3b_geq4j_ljets"],        "m_bb_closest_3b_geq4j_ljets"        , m_closest_bins);
      book(_h["pt_bb_closest_3b_geq4j_ljets"],       "pt_bb_closest_3b_geq4j_ljets"       , pt_bins);
      book(_h["dR_bb_closest_3b_geq4j_ljets"],       "dR_bb_closest_3b_geq4j_ljets"       , dr_bins);

      book(_h["dR_bb_average_3b_geq4j_ljets"],       "dR_bb_average_3b_geq4j_ljets"       , dr_bins);

      book(_h["m_bb_leadingVec_3b_geq4j_ljets"],     "m_bb_leadingVec_3b_geq4j_ljets"     , m_leading_bins);
      book(_h["pt_bb_leadingVec_3b_geq4j_ljets"],    "pt_bb_leadingVec_3b_geq4j_ljets"    , pt_bins);
      book(_h["dR_bb_leadingVec_3b_geq4j_ljets"],    "dR_bb_leadingVec_3b_geq4j_ljets"    , dr_bins);


      // geq4b geq4j category 
      book(_h["N_Jets_geq4b_geq4j_ljets"],           "N_Jets_geq4b_geq4j_ljets"           , multiplicity_bins);
      book(_h["N_b_Jets_geq4b_geq4j_ljets"],         "N_b_Jets_geq4b_geq4j_ljets"         , multiplicity_bins);
      book(_h["all_bjets_pt_geq4b_geq4j_ljets"],     "all_bjets_pt_geq4b_geq4j_ljets"     , pt_bins);

      book(_h["ht_bjets_geq4b_geq4j_ljets"],         "ht_bjets_geq4b_geq4j_ljets"         , ht_bins); 
      book(_h["ht_geq4b_geq4j_ljets"],               "ht_geq4b_geq4j_ljets"               , ht_bins); 
      book(_h["ht_had_geq4b_geq4j_ljets"],           "ht_had_geq4b_geq4j_ljets"           , ht_bins);

      book(_h["lead_bjet_pt_geq4b_geq4j_ljets"],     "lead_bjet_pt_geq4b_geq4j_ljets"     , pt_bins);
      book(_h["sublead_bjet_pt_geq4b_geq4j_ljets"],  "sublead_bjet_pt_geq4b_geq4j_ljets"  , pt_bins);
      book(_h["third_bjet_pt_geq4b_geq4j_ljets"],    "third_bjet_pt_geq4b_geq4j_ljets"    , pt_bins);
      book(_h["fourth_bjet_pt_geq4b_geq4j_ljets"],   "fourth_bjet_pt_geq4b_geq4j_ljets"   , pt_bins);

      book(_h["m_bb_leading_geq4b_geq4j_ljets"],     "m_bb_leading_geq4b_geq4j_ljets"     , m_leading_bins);
      book(_h["pt_bb_leading_geq4b_geq4j_ljets"],    "pt_bb_leading_geq4b_geq4j_ljets"    , pt_bins);
      book(_h["dR_bb_leading_geq4b_geq4j_ljets"],    "dR_bb_leading_geq4b_geq4j_ljets"    , dr_bins);

      book(_h["m_bb_closest_geq4b_geq4j_ljets"],     "m_bb_closest_geq4b_geq4j_ljets"     , m_closest_bins);
      book(_h["pt_bb_closest_geq4b_geq4j_ljets"],    "pt_bb_closest_geq4b_geq4j_ljets"    , pt_bins);
      book(_h["dR_bb_closest_geq4b_geq4j_ljets"],    "dR_bb_closest_geq4b_geq4j_ljets"    , dr_bins);

      book(_h["dR_bb_average_geq4b_geq4j_ljets"],    "dR_bb_average_geq4b_geq4j_ljets"    , dr_bins);

      book(_h["m_bb_leadingVec_geq4b_geq4j_ljets"],  "m_bb_leadingVec_geq4b_geq4j_ljets"     , m_leading_bins);
      book(_h["pt_bb_leadingVec_geq4b_geq4j_ljets"], "pt_bb_leadingVec_geq4b_geq4j_ljets"    , pt_bins);
      book(_h["dR_bb_leadingVec_geq4b_geq4j_ljets"], "dR_bb_leadingVec_geq4b_geq4j_ljets"    , dr_bins);


      // dilepton
      book(_h["N_Jets_3b_geq4j_dil"],              "N_Jets_3b_geq4j_dil"              , multiplicity_bins);
      book(_h["N_b_Jets_3b_geq4j_dil"],            "N_b_Jets_3b_geq4j_dil"            , multiplicity_bins);
      book(_h["all_bjets_pt_3b_geq4j_dil"],        "all_bjets_pt_3b_geq4j_dil"        , pt_bins);

      book(_h["ht_bjets_3b_geq4j_dil"],            "ht_bjets_3b_geq4j_dil"            , ht_bins);
      book(_h["ht_3b_geq4j_dil"],                  "ht_3b_geq4j_dil"                  , ht_bins);
      book(_h["ht_had_3b_geq4j_dil"],              "ht_had_3b_geq4j_dil"              , ht_bins);

      book(_h["lead_bjet_pt_3b_geq4j_dil"],        "lead_bjet_pt_3b_geq4j_dil"        , pt_bins);
      book(_h["sublead_bjet_pt_3b_geq4j_dil"],     "sublead_bjet_pt_3b_geq4j_dil"     , pt_bins);
      book(_h["third_bjet_pt_3b_geq4j_dil"],       "third_bjet_pt_3b_geq4j_dil"       , pt_bins);
      book(_h["fourth_bjet_pt_3b_geq4j_dil"],      "fourth_bjet_pt_3b_geq4j_dil"      , pt_bins);

      book(_h["m_bb_leading_3b_geq4j_dil"],        "m_bb_leading_3b_geq4j_dil"        , m_leading_bins);
      book(_h["pt_bb_leading_3b_geq4j_dil"],       "pt_bb_leading_3b_geq4j_dil"       , pt_bins);
      book(_h["dR_bb_leading_3b_geq4j_dil"],       "dR_bb_leading_3b_geq4j_dil"       , dr_bins);

      book(_h["m_bb_closest_3b_geq4j_dil"],        "m_bb_closest_3b_geq4j_dil"        , m_closest_bins);
      book(_h["pt_bb_closest_3b_geq4j_dil"],       "pt_bb_closest_3b_geq4j_dil"       , pt_bins);
      book(_h["dR_bb_closest_3b_geq4j_dil"],       "dR_bb_closest_3b_geq4j_dil"       , dr_bins);

      book(_h["dR_bb_average_3b_geq4j_dil"],       "dR_bb_average_3b_geq4j_dil"       , dr_bins);

      book(_h["m_bb_leadingVec_3b_geq4j_dil"],     "m_bb_leadingVec_3b_geq4j_dil"     , m_leading_bins);
      book(_h["pt_bb_leadingVec_3b_geq4j_dil"],    "pt_bb_leadingVec_3b_geq4j_dil"    , pt_bins);
      book(_h["dR_bb_leadingVec_3b_geq4j_dil"],    "dR_bb_leadingVec_3b_geq4j_dil"    , dr_bins);


      // geq4b geq4j category                                                                                                                                                                                                                                                                                   
      book(_h["N_Jets_geq4b_geq4j_dil"],           "N_Jets_geq4b_geq4j_dil"           , multiplicity_bins);
      book(_h["N_b_Jets_geq4b_geq4j_dil"],         "N_b_Jets_geq4b_geq4j_dil"         , multiplicity_bins);
      book(_h["all_bjets_pt_geq4b_geq4j_dil"],     "all_bjets_pt_geq4b_geq4j_dil"     , pt_bins);

      book(_h["ht_bjets_geq4b_geq4j_dil"],         "ht_bjets_geq4b_geq4j_dil"         , ht_bins);
      book(_h["ht_geq4b_geq4j_dil"],               "ht_geq4b_geq4j_dil"               , ht_bins);
      book(_h["ht_had_geq4b_geq4j_dil"],           "ht_had_geq4b_geq4j_dil"           , ht_bins);

      book(_h["lead_bjet_pt_geq4b_geq4j_dil"],     "lead_bjet_pt_geq4b_geq4j_dil"     , pt_bins);
      book(_h["sublead_bjet_pt_geq4b_geq4j_dil"],  "sublead_bjet_pt_geq4b_geq4j_dil"  , pt_bins);
      book(_h["third_bjet_pt_geq4b_geq4j_dil"],    "third_bjet_pt_geq4b_geq4j_dil"    , pt_bins);
      book(_h["fourth_bjet_pt_geq4b_geq4j_dil"],   "fourth_bjet_pt_geq4b_geq4j_dil"   , pt_bins);

      book(_h["m_bb_leading_geq4b_geq4j_dil"],     "m_bb_leading_geq4b_geq4j_dil"     , m_leading_bins);
      book(_h["pt_bb_leading_geq4b_geq4j_dil"],    "pt_bb_leading_geq4b_geq4j_dil"    , pt_bins);
      book(_h["dR_bb_leading_geq4b_geq4j_dil"],    "dR_bb_leading_geq4b_geq4j_dil"    , dr_bins);

      book(_h["m_bb_closest_geq4b_geq4j_dil"],     "m_bb_closest_geq4b_geq4j_dil"     , m_closest_bins);
      book(_h["pt_bb_closest_geq4b_geq4j_dil"],    "pt_bb_closest_geq4b_geq4j_dil"    , pt_bins);
      book(_h["dR_bb_closest_geq4b_geq4j_dil"],    "dR_bb_closest_geq4b_geq4j_dil"    , dr_bins);

      book(_h["dR_bb_average_geq4b_geq4j_dil"],    "dR_bb_average_geq4b_geq4j_dil"    , dr_bins);

      book(_h["m_bb_leadingVec_geq4b_geq4j_dil"],  "m_bb_leadingVec_geq4b_geq4j_dil"     , m_leading_bins);
      book(_h["pt_bb_leadingVec_geq4b_geq4j_dil"], "pt_bb_leadingVec_geq4b_geq4j_dil"    , pt_bins);
      book(_h["dR_bb_leadingVec_geq4b_geq4j_dil"], "dR_bb_leadingVec_geq4b_geq4j_dil"    , dr_bins);

      book(_h["abs_weight_1000_3b_geq4j_ljets"],    "abs_weight_1000_3b_geq4j_ljets",    1400, 0.0, 1400.0);
      book(_h["abs_weight_1000_geq4b_geq4j_ljets"], "abs_weight_1000_geq4b_geq4j_ljets", 1400, 0.0, 1400.0);
      book(_h["abs_weight_1000_3b_geq4j_dil"],      "abs_weight_1000_3b_geq4j_dil",      1400, 0.0, 1400.0);
      book(_h["abs_weight_1000_geq4b_geq4j_dil"],   "abs_weight_1000_geq4b_geq4j_dil",   1400, 0.0, 1400.0);

      book(_h["abs_weight_100_3b_geq4j_ljets"],     "abs_weight_100_3b_geq4j_ljets",      100, 0.0, 100.0);
      book(_h["abs_weight_100_geq4b_geq4j_ljets"],  "abs_weight_100_geq4b_geq4j_ljets",   100, 0.0, 100.0);
      book(_h["abs_weight_100_3b_geq4j_dil"],       "abs_weight_100_3b_geq4j_dil",        100, 0.0, 100.0);
      book(_h["abs_weight_100_geq4b_geq4j_dil"],    "abs_weight_100_geq4b_geq4j_dil",     100, 0.0, 100.0);

      book(_h["abs_weight_5_3b_geq4j_ljets"],       "abs_weight_5_3b_geq4j_ljets",          5, 0.0, 5.0);
      book(_h["abs_weight_5_geq4b_geq4j_ljets"],    "abs_weight_5_geq4b_geq4j_ljets",       5, 0.0, 5.0);
      book(_h["abs_weight_5_3b_geq4j_dil"],         "abs_weight_5_3b_geq4j_dil",            5, 0.0, 5.0);
      book(_h["abs_weight_5_geq4b_geq4j_dil"],      "abs_weight_5_geq4b_geq4j_dil",         5, 0.0, 5.0);

      book(_h["weight_sign_3b_geq4j_ljets"],        "weight_sign_3b_geq4j_ljets",           2, -1.0, 1.0);
      book(_h["weight_sign_geq4b_geq4j_ljets"],     "weight_sign_geq4b_geq4j_ljets",        2, -1.0, 1.0);
      book(_h["weight_sign_3b_geq4j_dil"],          "weight_sign_3b_geq4j_dil",             2, -1.0, 1.0);
      book(_h["weight_sign_geq4b_geq4j_dil"],       "weight_sign_geq4b_geq4j_dil",          2, -1.0, 1.0);


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      /// @todo Do the event by event analysis here

      // retrieve dressed leptons, sorted by pT
        // vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();
      vector<DressedLepton> leptons;
      for (auto &lep : apply<DressedLeptons>(event, "muons").dressedLeptons()) { leptons.push_back(lep); }
      for (auto &lep : apply<DressedLeptons>(event, "elecs").dressedLeptons()) { leptons.push_back(lep); }
        // retrieve dressed leptons, sorted by pT

        // retrieve clustered jets, sorted by pT, with a minimum pT cut
      // Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);
      // remove all jets within dR < 0.2 of a dressed lepton
      // idiscardIfAnyDeltaRLess(jets, leptons, 0.2);
      const Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      for (const auto& jet : jets) {
        ifilter_discard(leptons, [&](const DressedLepton& lep) { return deltaR(jet, lep) < 0.4; });
      }
      // select jets ghost-associated to B-hadrons with a certain fiducial selection
      // Jets bjets = filter_select(jets, [](const Jet& jet) {
      //   return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      // });
      Jets bjets;
      for (const Jet& jet : jets) {
        if (jet.bTagged(Cuts::pT >= 5*GeV))  bjets += jet;
      }

      size_t njets = jets.size();
      size_t nbjets = bjets.size();
      // veto event if there are no b-jets
      // if (bjets.empty())  vetoEvent;
      // apply a missing-momentum cut
      // if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;
      bool pass_ljets = (leptons.size() == 1 && leptons[0].pT() > 27*GeV);
      bool pass_dil   = (leptons.size() == 2 && leptons[0].pT() > 27*GeV && leptons[1].pT() > 27*GeV);

      if (!(pass_ljets || pass_dil)) vetoEvent;
      if (nbjets < 3 || njets < 4)  vetoEvent;
     
      // fill histogram with leading b-jet pT
      // _h["XXXX"]->fill(bjets[0].pT()/GeV);

      double hthad = sum(jets, pT, 0.0);
      double ht = sum(leptons, pT, hthad);
      FourMomentum jsum = bjets[0].momentum() + bjets[1].momentum();
      double dr_leading = deltaR(bjets[0], bjets[1]);
      size_t ind1, ind2; double mindr = 999.;
      size_t ind1_vec, ind2_vec;
      double vec_pt = 0.0;
      double sum_dr = 0.0;
      size_t sum_n_dr = 0; 
      for (size_t i = 0; i < bjets.size(); ++i) {
        for (size_t j = 0; j < bjets.size(); ++j) {
          if (i == j)  continue;
          double dr = deltaR(bjets[i], bjets[j]);

	  double pt = (bjets[i].momentum()+bjets[j].momentum()).pT();

          if (dr < mindr) {
            ind1 = i;
            ind2 = j;
            mindr = dr;
          }
	  if (pt > vec_pt){

	    ind1_vec = i;
	    ind2_vec = j;

	    vec_pt = pt;

	  }
	  
          sum_dr += dr;
          sum_n_dr += 1; 
        }
      }

      double ht_bjets = 0;
      for (size_t i = 0; i < bjets.size(); ++i) {
          ht_bjets = ht_bjets + bjets[i].pT();
      }

      FourMomentum bb_closest = bjets[ind1].momentum() + bjets[ind2].momentum();
      double dr_closest = deltaR(bjets[ind1], bjets[ind2]);
     
      double dr_leading_vec = deltaR(bjets[ind1_vec], bjets[ind2_vec]);
      double pt_leading_vec = (bjets[ind1_vec].momentum()+bjets[ind2_vec].momentum()).pT();
      double m_leading_vec  = (bjets[ind1_vec].momentum()+bjets[ind2_vec].momentum()).mass();

      // lets do the 3bjets geq4njets category first
      if (pass_ljets && (nbjets == 3 && njets >= 4)) {
          _h["N_Jets_3b_geq4j_ljets"]               -> fill(njets);
          _h["N_b_Jets_3b_geq4j_ljets"]             -> fill(nbjets);
          for (size_t i = 0; i < bjets.size(); ++i) {
              _h["all_bjets_pt_3b_geq4j_ljets"]     -> fill(bjets[i].pT()/GeV);
          }
          _h["ht_bjets_3b_geq4j_ljets"]             -> fill(ht_bjets/GeV);
          // b-jet pTs
          _h["lead_bjet_pt_3b_geq4j_ljets"]         -> fill(bjets[0].pT()/GeV);
          _h["sublead_bjet_pt_3b_geq4j_ljets"]      -> fill(bjets[1].pT()/GeV);
          _h["third_bjet_pt_3b_geq4j_ljets"]        -> fill(bjets[2].pT()/GeV);

          // HT
          _h["ht_3b_geq4j_ljets"]                   -> fill(ht/GeV);
          _h["ht_had_3b_geq4j_ljets"]               -> fill(hthad/GeV);

          // leading bb pair
          _h["m_bb_leading_3b_geq4j_ljets"]         -> fill(jsum.mass()/GeV);
          _h["pt_bb_leading_3b_geq4j_ljets"]        -> fill(jsum.pT()/GeV);
          _h["dR_bb_leading_3b_geq4j_ljets"]        -> fill(dr_leading);

          // closest bb pair
          _h["m_bb_closest_3b_geq4j_ljets"]         -> fill(bb_closest.mass()/GeV);
          _h["pt_bb_closest_3b_geq4j_ljets"]        -> fill(bb_closest.pT()/GeV);
          _h["dR_bb_closest_3b_geq4j_ljets"]        -> fill(dr_closest);

	  // bb pair with highest vectorial sum pt                                            
          _h["m_bb_leadingVec_3b_geq4j_ljets"]         -> fill(m_leading_vec/GeV);
          _h["pt_bb_leadingVec_3b_geq4j_ljets"]        -> fill(pt_leading_vec/GeV);
          _h["dR_bb_leadingVec_3b_geq4j_ljets"]        -> fill(dr_leading_vec);

          // average dR
          _h["dR_bb_average_3b_geq4j_ljets"]        -> fill(sum_dr/sum_n_dr);


	  _h["abs_weight_1000_3b_geq4j_ljets"] -> fill(event.weights()[0]);
	  _h["abs_weight_100_3b_geq4j_ljets"]  -> fill(event.weights()[0]);
	  _h["abs_weight_5_3b_geq4j_ljets"]    -> fill(event.weights()[0]);
	  _h["weight_sign_3b_geq4j_ljets"]     -> fill(event.weights()[0]);


      }
      // lets do the geq4bjets geq4njets category now
      if (pass_ljets && (nbjets >= 4 && njets >= 4)) {
          _h["N_Jets_geq4b_geq4j_ljets"]            -> fill(njets);
          _h["N_b_Jets_geq4b_geq4j_ljets"]          -> fill(nbjets);
          for (size_t i = 0; i < bjets.size(); ++i) {
              _h["all_bjets_pt_geq4b_geq4j_ljets"]  -> fill(bjets[i].pT()/GeV);
          }
          _h["ht_bjets_geq4b_geq4j_ljets"]          -> fill(ht_bjets/GeV);
          // b-jet pTs
          _h["lead_bjet_pt_geq4b_geq4j_ljets"]      -> fill(bjets[0].pT()/GeV);
          _h["sublead_bjet_pt_geq4b_geq4j_ljets"]   -> fill(bjets[1].pT()/GeV);
          _h["third_bjet_pt_geq4b_geq4j_ljets"]     -> fill(bjets[2].pT()/GeV);              
          _h["fourth_bjet_pt_geq4b_geq4j_ljets"]    -> fill(bjets[3].pT()/GeV);

          // HT
          _h["ht_geq4b_geq4j_ljets"]                -> fill(ht/GeV);
          _h["ht_had_geq4b_geq4j_ljets"]            -> fill(hthad/GeV);

          // leading bb pair
          _h["m_bb_leading_geq4b_geq4j_ljets"]      -> fill(jsum.mass()/GeV);
          _h["pt_bb_leading_geq4b_geq4j_ljets"]     -> fill(jsum.pT()/GeV);
          _h["dR_bb_leading_geq4b_geq4j_ljets"]     -> fill(dr_leading);

          // closest bb pair
          _h["m_bb_closest_geq4b_geq4j_ljets"]      -> fill(bb_closest.mass()/GeV);
          _h["pt_bb_closest_geq4b_geq4j_ljets"]     -> fill(bb_closest.pT()/GeV);
          _h["dR_bb_closest_geq4b_geq4j_ljets"]     -> fill(dr_closest);

          // bb pair with highest vectorial sum pt                                                                                                                                                        
          _h["m_bb_leadingVec_geq4b_geq4j_ljets"]   -> fill(m_leading_vec/GeV);
          _h["pt_bb_leadingVec_geq4b_geq4j_ljets"]  -> fill(pt_leading_vec/GeV);
          _h["dR_bb_leadingVec_geq4b_geq4j_ljets"]  -> fill(dr_leading_vec);


          // average dR
          _h["dR_bb_average_geq4b_geq4j_ljets"]     -> fill(sum_dr/sum_n_dr);

          _h["abs_weight_1000_geq4b_geq4j_ljets"] -> fill(event.weights()[0]);
          _h["abs_weight_100_geq4b_geq4j_ljets"]  -> fill(event.weights()[0]);
          _h["abs_weight_5_geq4b_geq4j_ljets"]    -> fill(event.weights()[0]);
          _h["weight_sign_geq4b_geq4j_ljets"]     -> fill(event.weights()[0]);


      }
      
      // dil
      // lets do the 3bjets geq4njets category first                                                                                                                                             
      if (pass_dil && (nbjets == 3 && njets >= 4)) {
	_h["N_Jets_3b_geq4j_dil"]               -> fill(njets);
	_h["N_b_Jets_3b_geq4j_dil"]             -> fill(nbjets);
	for (size_t i = 0; i < bjets.size(); ++i) {
	  _h["all_bjets_pt_3b_geq4j_dil"]     -> fill(bjets[i].pT()/GeV);
	}
	_h["ht_bjets_3b_geq4j_dil"]             -> fill(ht_bjets/GeV);
	// b-jet pTs                                                                                                                                                                                     
	_h["lead_bjet_pt_3b_geq4j_dil"]         -> fill(bjets[0].pT()/GeV);
	_h["sublead_bjet_pt_3b_geq4j_dil"]      -> fill(bjets[1].pT()/GeV);
	_h["third_bjet_pt_3b_geq4j_dil"]        -> fill(bjets[2].pT()/GeV);

	// HT                                                                                                                                                                                             
	_h["ht_3b_geq4j_dil"]                   -> fill(ht/GeV);
	_h["ht_had_3b_geq4j_dil"]               -> fill(hthad/GeV);

	// leading bb pair                                                                                                                                                                               
	_h["m_bb_leading_3b_geq4j_dil"]         -> fill(jsum.mass()/GeV);
	_h["pt_bb_leading_3b_geq4j_dil"]        -> fill(jsum.pT()/GeV);
	_h["dR_bb_leading_3b_geq4j_dil"]        -> fill(dr_leading);

	// closest bb pair                                                                                                                                                                                
	_h["m_bb_closest_3b_geq4j_dil"]         -> fill(bb_closest.mass()/GeV);
	_h["pt_bb_closest_3b_geq4j_dil"]        -> fill(bb_closest.pT()/GeV);
	_h["dR_bb_closest_3b_geq4j_dil"]        -> fill(dr_closest);

	// bb pair with highest vectorial sum pt                                                                                                                                                          
	_h["m_bb_leadingVec_3b_geq4j_dil"]         -> fill(m_leading_vec/GeV);
	_h["pt_bb_leadingVec_3b_geq4j_dil"]        -> fill(pt_leading_vec/GeV);
	_h["dR_bb_leadingVec_3b_geq4j_dil"]        -> fill(dr_leading_vec);

	// average dR                                                                                                                                                                                     
	_h["dR_bb_average_3b_geq4j_dil"]        -> fill(sum_dr/sum_n_dr);

	_h["abs_weight_1000_3b_geq4j_dil"] -> fill(event.weights()[0]);
	_h["abs_weight_100_3b_geq4j_dil"]  -> fill(event.weights()[0]);
	_h["abs_weight_5_3b_geq4j_dil"]    -> fill(event.weights()[0]);
	_h["weight_sign_3b_geq4j_dil"]     -> fill(event.weights()[0]);


      }
      // lets do the geq4bjets geq4njets category now                                                                                                                                                     
      if (pass_dil && (nbjets >= 4 && njets >= 4)) {
	_h["N_Jets_geq4b_geq4j_dil"]            -> fill(njets);
	_h["N_b_Jets_geq4b_geq4j_dil"]          -> fill(nbjets);
	for (size_t i = 0; i < bjets.size(); ++i) {
	  _h["all_bjets_pt_geq4b_geq4j_dil"]  -> fill(bjets[i].pT()/GeV);
	}
	_h["ht_bjets_geq4b_geq4j_dil"]          -> fill(ht_bjets/GeV);
	// b-jet pTs                                                                                                                                                                                      
	_h["lead_bjet_pt_geq4b_geq4j_dil"]      -> fill(bjets[0].pT()/GeV);
	_h["sublead_bjet_pt_geq4b_geq4j_dil"]   -> fill(bjets[1].pT()/GeV);
	_h["third_bjet_pt_geq4b_geq4j_dil"]     -> fill(bjets[2].pT()/GeV);
	_h["fourth_bjet_pt_geq4b_geq4j_dil"]    -> fill(bjets[3].pT()/GeV);

	// HT                                                                                                                                                                                             
	_h["ht_geq4b_geq4j_dil"]                -> fill(ht/GeV);
	_h["ht_had_geq4b_geq4j_dil"]            -> fill(hthad/GeV);

	// leading bb pair                                                                                                                                                                                
	_h["m_bb_leading_geq4b_geq4j_dil"]      -> fill(jsum.mass()/GeV);
	_h["pt_bb_leading_geq4b_geq4j_dil"]     -> fill(jsum.pT()/GeV);
	_h["dR_bb_leading_geq4b_geq4j_dil"]     -> fill(dr_leading);

	// closest bb pair                                                                                                                                                                                
	_h["m_bb_closest_geq4b_geq4j_dil"]      -> fill(bb_closest.mass()/GeV);
	_h["pt_bb_closest_geq4b_geq4j_dil"]     -> fill(bb_closest.pT()/GeV);
	_h["dR_bb_closest_geq4b_geq4j_dil"]     -> fill(dr_closest);

	// bb pair with highest vectorial sum pt                                                                                                                                                          
	_h["m_bb_leadingVec_geq4b_geq4j_dil"]   -> fill(m_leading_vec/GeV);
	_h["pt_bb_leadingVec_geq4b_geq4j_dil"]  -> fill(pt_leading_vec/GeV);
	_h["dR_bb_leadingVec_geq4b_geq4j_dil"]  -> fill(dr_leading_vec);


	// average dR                                                                                                                                                                                     
	_h["dR_bb_average_geq4b_geq4j_dil"]     -> fill(sum_dr/sum_n_dr);

	_h["abs_weight_1000_geq4b_geq4j_dil"] -> fill(event.weights()[0]);
	_h["abs_weight_100_geq4b_geq4j_dil"]  -> fill(event.weights()[0]);
	_h["abs_weight_5_geq4b_geq4j_dil"]    -> fill(event.weights()[0]);
	_h["weight_sign_geq4b_geq4j_dil"]     -> fill(event.weights()[0]);
	

      }



    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize(_h["YYYY"]); // normalize to unity
      // scale(_h["ZZZZ"], crossSection()/picobarn/sumOfWeights()); // norm to cross section
      const double sf = crossSection() / picobarn / sumOfWeights();
      for (auto const& h : _h) {
          scale(h.second, sf);
          for (size_t i =0; i < h.second -> numBins(); i++) {
            h.second -> bin(i).scaleW(1/h.second -> bin(i).width());
          }
          // normalize(h.second, 1.0);
      }
    }

    //@}


    /// @name Histograms
    //@{

    map<string, Histo1DPtr> _h;
    // map<string, Profile1DPtr> _p;
    // map<string, CounterPtr> _c;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ttbb_analysis);


}
