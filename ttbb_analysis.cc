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

      vector<double> genNjets                     = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};
      vector<double> genNbjets                    = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};
      vector<double> genJetPt                     = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500};
      vector<double> genbjetspt                   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500};
      vector<double> genHTjets                    = {0.00000, 34.09091, 68.18182, 102.27275, 136.36365, 170.45455, 204.54545, 238.63635, 272.72725, 306.81815, 340.90905, 375.00005, 409.09095, 443.18185, 477.27275, 511.36365, 545.45455, 579.54545, 613.63635, 647.72725, 681.81815, 715.90905, 750.00005, 784.09095, 818.18185, 852.27275, 886.36365, 920.45455, 954.54545, 988.63655, 1022.72755, 1056.81855, 1090.90955, 1124.99955, 1159.09055, 1193.18155, 1227.27255, 1261.36355, 1295.45455, 1329.54555, 1363.63655, 1397.72755, 1431.81855, 1465.90955, 1500.00045};
      vector<double> genEvt_M_MinDeltaRGenBJets   = {0.000, 15.625, 31.250, 46.875, 62.500, 78.125, 93.750, 109.375, 125.000, 140.625, 156.250, 171.875, 187.500, 203.125, 218.750, 234.375, 250.000, 265.625, 281.250, 296.875, 312.500, 328.125, 343.750, 359.375, 375.000, 390.625, 406.250, 421.875, 437.500, 453.125, 468.750, 484.375, 500.000};
      vector<double> genEvt_M_HardestGenBJets     = {0.000, 15.625, 31.250, 46.875, 62.500, 78.125, 93.750, 109.375, 125.000, 140.625, 156.250, 171.875, 187.500, 203.125, 218.750, 234.375, 250.000, 265.625, 281.250, 296.875, 312.500, 328.125, 343.750, 359.375, 375.000, 390.625, 406.250, 421.875, 437.500, 453.125, 468.750, 484.375, 500.000};
      vector<double> genEvt_Dr_MinDeltaRGenBJets  = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0};
      vector<double> genEvt_Dr_HardestGenBJets    = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0};

      // Book histograms
      // specify custom binning
      // 3b geq4j category 
      book(_h["N_Jets_3b_geq4j_ljets"],              "N_Jets_3b_geq4j_ljets"              , genNjets);
      book(_h["N_b_Jets_3b_geq4j_ljets"],            "N_b_Jets_3b_geq4j_ljets"            , genNbjets);
      book(_h["all_bjets_pt_3b_geq4j_ljets"],        "all_bjets_pt_3b_geq4j_ljets"        , genbjetspt);

      book(_h["ht_bjets_3b_geq4j_ljets"],            "ht_bjets_3b_geq4j_ljets"            , genHTjets); 
      book(_h["ht_3b_geq4j_ljets"],                  "ht_3b_geq4j_ljets"                  , genHTjets); 
      book(_h["ht_had_3b_geq4j_ljets"],              "ht_had_3b_geq4j_ljets"              , genHTjets);

      book(_h["lead_bjet_pt_3b_geq4j_ljets"],        "lead_bjet_pt_3b_geq4j_ljets"        , genJetPt);
      book(_h["sublead_bjet_pt_3b_geq4j_ljets"],     "sublead_bjet_pt_3b_geq4j_ljets"     , genJetPt);
      book(_h["third_bjet_pt_3b_geq4j_ljets"],       "third_bjet_pt_3b_geq4j_ljets"       , genJetPt);
      book(_h["fourth_bjet_pt_3b_geq4j_ljets"],      "fourth_bjet_pt_3b_geq4j_ljets"      , genJetPt);

      book(_h["m_bb_leading_3b_geq4j_ljets"],        "m_bb_leading_3b_geq4j_ljets"        , genEvt_M_HardestGenBJets);
      book(_h["pt_bb_leading_3b_geq4j_ljets"],       "pt_bb_leading_3b_geq4j_ljets"       , genJetPt);
      book(_h["dR_bb_leading_3b_geq4j_ljets"],       "dR_bb_leading_3b_geq4j_ljets"       , genEvt_Dr_HardestGenBJets);

      book(_h["m_bb_closest_3b_geq4j_ljets"],        "m_bb_closest_3b_geq4j_ljets"        , genEvt_M_MinDeltaRGenBJets);
      book(_h["pt_bb_closest_3b_geq4j_ljets"],       "pt_bb_closest_3b_geq4j_ljets"       , genJetPt);
      book(_h["dR_bb_closest_3b_geq4j_ljets"],       "dR_bb_closest_3b_geq4j_ljets"       , genEvt_Dr_MinDeltaRGenBJets);

      // geq4b geq4j category 
      book(_h["N_Jets_geq4b_geq4j_ljets"],           "N_Jets_geq4b_geq4j_ljets"           , genNjets);
      book(_h["N_b_Jets_geq4b_geq4j_ljets"],         "N_b_Jets_geq4b_geq4j_ljets"         , genNbjets);
      book(_h["all_bjets_pt_geq4b_geq4j_ljets"],     "all_bjets_pt_geq4b_geq4j_ljets"     , genbjetspt);

      book(_h["ht_bjets_geq4b_geq4j_ljets"],         "ht_bjets_geq4b_geq4j_ljets"         , genHTjets); 
      book(_h["ht_geq4b_geq4j_ljets"],               "ht_geq4b_geq4j_ljets"               , genHTjets); 
      book(_h["ht_had_geq4b_geq4j_ljets"],           "ht_had_geq4b_geq4j_ljets"           , genHTjets);

      book(_h["lead_bjet_pt_geq4b_geq4j_ljets"],     "lead_bjet_pt_geq4b_geq4j_ljets"     , genJetPt);
      book(_h["sublead_bjet_pt_geq4b_geq4j_ljets"],  "sublead_bjet_pt_geq4b_geq4j_ljets"  , genJetPt);
      book(_h["third_bjet_pt_geq4b_geq4j_ljets"],    "third_bjet_pt_geq4b_geq4j_ljets"    , genJetPt);
      book(_h["fourth_bjet_pt_geq4b_geq4j_ljets"],   "fourth_bjet_pt_geq4b_geq4j_ljets"   , genJetPt);

      book(_h["m_bb_leading_geq4b_geq4j_ljets"],     "m_bb_leading_geq4b_geq4j_ljets"     , genEvt_M_HardestGenBJets);
      book(_h["pt_bb_leading_geq4b_geq4j_ljets"],    "pt_bb_leading_geq4b_geq4j_ljets"    , genJetPt);
      book(_h["dR_bb_leading_geq4b_geq4j_ljets"],    "dR_bb_leading_geq4b_geq4j_ljets"    , genEvt_Dr_HardestGenBJets);

      book(_h["m_bb_closest_geq4b_geq4j_ljets"],     "m_bb_closest_geq4b_geq4j_ljets"     , genEvt_M_MinDeltaRGenBJets);
      book(_h["pt_bb_closest_geq4b_geq4j_ljets"],    "pt_bb_closest_geq4b_geq4j_ljets"    , genJetPt);
      book(_h["dR_bb_closest_geq4b_geq4j_ljets"],    "dR_bb_closest_geq4b_geq4j_ljets"    , genEvt_Dr_MinDeltaRGenBJets);


      // book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      // book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      // book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      // book(_h["AAAA"], 1, 1, 1);
      // book(_p["BBBB"], 2, 1, 1);
      // book(_c["CCCC"], 3, 1, 1);

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
      if (!pass_ljets) vetoEvent;
      if (pass_ljets && (nbjets < 3 || njets < 4))  vetoEvent;

      // fill histogram with leading b-jet pT
      // _h["XXXX"]->fill(bjets[0].pT()/GeV);

      double hthad = sum(jets, pT, 0.0);
      double ht = sum(leptons, pT, hthad);
      FourMomentum jsum = bjets[0].momentum() + bjets[1].momentum();
      double dr_leading = deltaR(bjets[0], bjets[1]);
      size_t ind1, ind2; double mindr = 999.;
      for (size_t i = 0; i < bjets.size(); ++i) {
        for (size_t j = 0; j < bjets.size(); ++j) {
          if (i == j)  continue;
          double dr = deltaR(bjets[i], bjets[j]);
          if (dr < mindr) {
            ind1 = i;
            ind2 = j;
            mindr = dr;
          }
        }
      }
      double ht_bjets = 0;
      for (size_t i = 0; i < bjets.size(); ++i) {
          ht_bjets = ht_bjets + bjets[i].pT();
      }

      FourMomentum bb_closest = bjets[ind1].momentum() + bjets[ind2].momentum();
      double dr_closest = deltaR(bjets[ind1], bjets[ind2]);
     
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
          if (nbjets >= 4) { 
              _h["fourth_bjet_pt_3b_geq4j_ljets"]   -> fill(bjets[3].pT()/GeV);
          }
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
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize(_h["YYYY"]); // normalize to unity
      // scale(_h["ZZZZ"], crossSection()/picobarn/sumOfWeights()); // norm to cross section
      const double sf = crossSection() / picobarn / sumOfWeights();
      for (auto const& h : _h) {
          scale(h.second, sf);
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
