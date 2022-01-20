#include "RecoTauTag/RecoTau/interface/MVADMDiscriminationProducerBase.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include <Math/VectorUtil.h>
#include "TMVA/Reader.h"

/* class PATTauDiscriminationMVADM
 *
 *  Seperate tau decay types using MVA
 *  Returns MVA score for each class
 *
 */

namespace {

template <class T, class U>
bool sortStrips (std::pair<T,U> i, std::pair<T,U> j) {
  return (i.first.pt() > j.first.pt());
}

template <class T>
bool sortByPT (T i, T j) {
  return (i->pt() > j->pt());
}


class PATTauDiscriminationMVADM final : public PATTauMVADMDiscriminationProducerBase  {
  public:
    explicit PATTauDiscriminationMVADM(const edm::ParameterSet& iConfig)
        :PATTauMVADMDiscriminationProducerBase(iConfig){
          version_ = iConfig.getParameter<std::string>("version");
          TString input_name_dm_10_applytoeven;
          TString input_name_dm_10_applytoodd;
          TString input_name_dm_0_1_applytoeven;
          TString input_name_dm_0_1_applytoodd;

          if(version_ == "MVADM_2016_v1") {
            input_name_dm_10_applytoeven = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_applytoeven_2016v1_dm10.xml";
            input_name_dm_10_applytoodd = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_applytoodd_2016v1_dm10.xml";
            input_name_dm_0_1_applytoeven = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_applytoeven_2016v1_dm0_dm1.xml";
            input_name_dm_0_1_applytoodd = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_applytoodd_2016v1_dm0_dm1.xml";
          }
          if(version_ == "MVADM_2017_v1") {
            input_name_dm_10_applytoeven = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_applytoeven_2017v1_dm10.xml";
            input_name_dm_10_applytoodd = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_applytoodd_2017v1_dm10.xml";
            input_name_dm_0_1_applytoeven = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_applytoeven_2017v1_dm0_dm1.xml";
            input_name_dm_0_1_applytoodd = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_applytoodd_2017v1_dm0_dm1.xml";
          }
          else {
            cms::Exception("MVA DM version not found") << "Requested version of ID does not exist."; 
          }

          reader_even_ = new TMVA::Reader();
          reader_odd_ = new TMVA::Reader();
          reader_dm10_even_ = new TMVA::Reader();
          reader_dm10_odd_ = new TMVA::Reader();
   
          for(unsigned i=0; i<(unsigned)var_names_.size(); ++i){
            reader_even_->AddVariable( var_names_[i], &(vars_[i]) );
            reader_odd_->AddVariable( var_names_[i], &(vars_[i]) );
          }
          for(unsigned i=0; i<(unsigned)var_names_dm10_.size(); ++i){
            reader_dm10_even_->AddVariable( var_names_dm10_[i], &(vars_dm10_[i]) );
            reader_dm10_odd_->AddVariable( var_names_dm10_[i], &(vars_dm10_[i]) );
          }

          reader_even_->BookMVA( "BDT method", input_name_dm_0_1_applytoeven );
          reader_odd_->BookMVA( "BDT method", input_name_dm_0_1_applytoodd );
          reader_dm10_even_->BookMVA( "BDT method", input_name_dm_10_applytoeven );
          reader_dm10_odd_->BookMVA( "BDT method", input_name_dm_10_applytoodd );

        }
    ~PATTauDiscriminationMVADM() override{}
    std::vector<double> read_mva_score(std::vector<float> vars, int decay_mode);

    void beginEvent(const edm::Event&, const edm::EventSetup&) override;
    std::vector<float> discriminate(const TauRef& tau) const override;


  private:


    TMVA::Reader *reader_even_;
    TMVA::Reader *reader_odd_;
    TMVA::Reader *reader_dm10_even_;
    TMVA::Reader *reader_dm10_odd_;

    std::string version_ = "MVADM_2016_v1";

    mutable std::vector<float> vars_ = std::vector<float>(24);
    mutable std::vector<float> vars_dm10_ = std::vector<float>(40);

    std::vector<TString> var_names_      = { "Egamma1_tau", "Egamma2_tau", "Epi_tau", "rho_dEta_tau", "rho_dphi_tau",
                                             "gammas_dEta_tau", "gammas_dR_tau", "DeltaR2WRTtau_tau", "tau_decay_mode",
                                             "eta", "pt", "Epi0", "Epi", "rho_dEta", "rho_dphi", "gammas_dEta", "Mrho", 
                                             "Mpi0", "DeltaR2WRTtau", "Mpi0_TwoHighGammas", "Mrho_OneHighGammas",
                                             "Mrho_TwoHighGammas", "Mrho_subleadingGamma", "strip_pt" };
    std::vector<TString> var_names_dm10_ = { "E1_overEa1", "E2_overEa1", "E1_overEtau", "E2_overEtau", "E3_overEtau",
                                             "a1_pi0_dEta_timesEtau", "a1_pi0_dphi_timesEtau", "h1_h2_dphi_timesE12",
                                             "h1_h2_dEta_timesE12", "h1_h3_dphi_timesE13", "h1_h3_dEta_timesE13",
                                             "h2_h3_dphi_timesE23", "h2_h3_dEta_timesE23", "gammas_dEta_timesEtau",
                                             "gammas_dR_timesEtau", "tau_decay_mode", "mass0", "mass1", "mass2", "E1", "E2",
                                             "E3", "strip_E", "a1_pi0_dEta", "a1_pi0_dphi", "strip_pt", "pt", "eta", "E",
                                             "h1_h2_dphi", "h1_h3_dphi", "h2_h3_dphi", "h1_h2_dEta", "h1_h3_dEta",
                                             "h2_h3_dEta", "Egamma1", "Egamma2", "gammas_dEta", "Mpi0",
                                             "Mpi0_TwoHighGammas"};

    bool isEven_ = false;
    unsigned long long event_;
    mutable std::vector<reco::CandidatePtr> gammas_;

    typedef ROOT::Math::PtEtaPhiEVector Vector;

    edm::Handle<TauCollection> taus_;

    std::vector<float> read_mva_score(int decay_mode) const {
      std::vector<float> mva_scores = {};
      if(decay_mode==0 || decay_mode==1) {
        if(isEven_) mva_scores = reader_even_->EvaluateMulticlass("BDT method");
        else       mva_scores = reader_odd_->EvaluateMulticlass("BDT method");
      } else if(decay_mode==10 || decay_mode==11) {
        if(isEven_) mva_scores = reader_dm10_even_->EvaluateMulticlass("BDT method");
        else       mva_scores = reader_dm10_odd_->EvaluateMulticlass("BDT method");
      }
      return mva_scores;
    }


    Vector GetPi0 (std::vector<reco::CandidatePtr> gammas, bool leadEtaPhi) const {
      Vector pi0;
      if(gammas.size()>0) {
        double E = 0.;
        double phi = 0.;
        double eta = 0.;
        for(auto g: gammas) {
          E+=g->energy();
          phi+=g->energy()*g->phi();
          eta+=g->energy()*g->eta();
        }
        eta/=E;
        phi/=E;
  
        if(leadEtaPhi){
          // if true sets the eta and phi of the pi0 to that of the leading gamma rather than using the weighted average
          eta = gammas[0]->eta();
          phi = gammas[0]->phi();
        }
  
        double mass = 0.1349;
        double p = sqrt(E*E-mass*mass);
        double theta = atan(exp(-eta))*2;
        double pt = p*sin(theta);
        pi0 = Vector(pt,eta,phi,E);
      }
      return pi0;
    }


    std::vector<std::pair<Vector,std::vector<reco::CandidatePtr>>> HPSGammas (std::vector<reco::CandidatePtr> cands) const {
      std::vector<std::pair<Vector,std::vector<reco::CandidatePtr>>> strips;   
      while(!cands.empty()) {
  
        std::vector<reco::CandidatePtr> Associated = {};
        std::vector<reco::CandidatePtr> notAssociated = {};
  
        Vector stripVector(0,0,0,0);
        stripVector=cands[0]->p4();
        Associated.push_back(cands[0]);
 
        bool repeat = true;
        unsigned int mini=1;
        while (repeat) {
          repeat = false;
          for(unsigned int i=mini;i<cands.size();++i) {
            double etaAssociationDistance = 0.20*pow(cands[i]->pt(),-0.66) + 0.20*pow(stripVector.Pt(),-0.66);
            double phiAssociationDistance = 0.35*pow(cands[i]->pt(),-0.71) + 0.35*pow(stripVector.Pt(),-0.71);
            etaAssociationDistance = std::min(etaAssociationDistance, 0.15);
            etaAssociationDistance = std::max(etaAssociationDistance, 0.05);
            phiAssociationDistance = std::min(phiAssociationDistance, 0.30);
            phiAssociationDistance = std::max(phiAssociationDistance, 0.05);
  
            if(fabs(cands[i]->eta()-stripVector.eta())<etaAssociationDistance &&
              fabs(ROOT::Math::VectorUtil::DeltaPhi(cands[i]->p4(),stripVector))<phiAssociationDistance) {
              stripVector+=cands[i]->p4();
              Associated.push_back(cands[i]);
              repeat = true;
            }
            else {
              notAssociated.push_back(cands[i]);
            }
          }
          cands.swap(notAssociated);
          notAssociated.clear(); 
          mini=0;
        }

        Vector strip = GetPi0(Associated, false);
        strips.push_back(std::make_pair(strip, Associated));
 
      }
      std::sort(strips.begin(), strips.end(), sortStrips<Vector, std::vector<reco::CandidatePtr>>);
  
      return strips;
    }


    std::pair<Vector,Vector> GetRho (const TauRef& tau, double gammas_pt_cut) const {
      Vector pi;
      Vector pi0;
      gammas_.clear();

      std::vector<reco::CandidatePtr> gammas;
      for (auto g: tau->signalGammaCands()) if(g->pt()>gammas_pt_cut) gammas.push_back(g);
      reco::CandidatePtrVector hads = tau->signalChargedHadrCands();

      if(hads.size()>0) pi = hads[0]->p4();

      double cone_size = std::max(std::min(0.1, 3./tau->pt()),0.05);
      std::vector<std::pair<Vector, std::vector<reco::CandidatePtr>>> strip_pairs = HPSGammas(gammas);
      std::vector<std::pair<Vector, std::vector<reco::CandidatePtr>>> strips_incone;
      for(auto s : strip_pairs) {
        if(std::fabs(ROOT::Math::VectorUtil::DeltaR(s.first,tau->p4()))<cone_size) strips_incone.push_back(s);
      }
      if(tau->decayMode()==0) {
        if(strips_incone.size()>0) {
          gammas = strips_incone[0].second;
        } else if(strip_pairs.size()>0) {
          gammas = strip_pairs[0].second;
        }
      }
      if(tau->decayMode()==1 && strip_pairs.size()>0) pi0 = GetPi0(strip_pairs[0].second, true);
      else {
        pi0 = GetPi0(gammas, true);
      }
      std::sort(gammas.begin(), gammas.end(), sortByPT<reco::CandidatePtr>); 
      gammas_ = gammas;
      return std::make_pair(pi,pi0);
    }

  std::pair<std::vector<Vector>, Vector> GetA1 (const TauRef& tau, float gammas_pt_cut) const {
    std::vector<Vector> prongs;
    Vector pi0;
    std::vector<reco::CandidatePtr> hads;
    for (auto h : tau->signalChargedHadrCands()) hads.push_back(h);
    if(hads.size()==3) {
      // arrange hadrons so the oppositly charged hadron is contained in the first element
      if(hads[1]->charge()!=hads[0]->charge()&&hads[1]->charge()!=hads[2]->charge()){
        auto temp = hads[1];
        hads[1] = hads[0];
        hads[0] = temp;
      }
      else if(hads[2]->charge()!=hads[0]->charge()&&hads[2]->charge()!=hads[1]->charge()){
        auto temp = hads[2];
        hads[2] = hads[0];
        hads[0] = temp;
      } 
      // from the two same sign hadrons place the one that gives the mass most similar to the rho meson as the second element
      double rho_mass = 0.7755;
      double dM1 = std::fabs((hads[0]->p4()+hads[1]->p4()).M()-rho_mass);
      double dM2 = std::fabs((hads[0]->p4()+hads[2]->p4()).M()-rho_mass);
      if(dM2<dM1){
        auto temp = hads[2];
        hads[2] = hads[1];
        hads[1] = temp;
      }
    }

    std::vector<reco::CandidatePtr> gammas_merge;
    for (auto g: tau->signalGammaCands()) if(g->pt()>gammas_pt_cut) gammas_merge.push_back(g);
    if(tau->decayMode()!=11){
      for (auto g: tau->isolationGammaCands()) if(g->pt()>gammas_pt_cut) gammas_merge.push_back(g); 
    } 
    std::sort(gammas_merge.begin(), gammas_merge.end(), sortByPT<reco::CandidatePtr>);
    std::vector<reco::CandidatePtr> gammas = {};
    for(auto g: gammas_merge) gammas.push_back(g);
    double cone_size = std::max(std::min(0.1, 3./tau->pt()),0.05);
    std::vector<std::pair<Vector, std::vector<reco::CandidatePtr>>> strip_pairs = HPSGammas(gammas);
    std::vector<std::pair<Vector, std::vector<reco::CandidatePtr>>> strips_incone; 
    for(auto s : strip_pairs) if(std::fabs(ROOT::Math::VectorUtil::DeltaR(s.first,tau->p4()))<cone_size) strips_incone.push_back(s);
     
    std::vector<reco::CandidatePtr> signal_gammas = {};

    if(strips_incone.size()>0) {
      signal_gammas = strips_incone[0].second;
    } else if(strip_pairs.size()>0) {
      signal_gammas = strip_pairs[0].second;
    }
    pi0 = GetPi0(signal_gammas, true);
    std::sort(signal_gammas.begin(), signal_gammas.end(), sortByPT<reco::CandidatePtr>);
    gammas_ = signal_gammas;

    for (auto h : hads) prongs.push_back((Vector)h->p4());

    return std::make_pair(prongs, pi0);
  }

 
};

void PATTauDiscriminationMVADM::beginEvent(const edm::Event& evt, const edm::EventSetup& es) {
  isEven_ = evt.id().event() % 2 == 0;
  event_ = evt.id().event();
  evt.getByToken(Tau_token, taus_);
}


std::vector<float> PATTauDiscriminationMVADM::discriminate(const TauRef& tau) const {
  double gammas_pt_cut;
  if (version_=="MVADM_2016_v1") gammas_pt_cut=0.5;
  if (version_=="MVADM_2017_v1") gammas_pt_cut=1.0;

  std::vector<float> scores= {};
  gammas_.clear();
  // define all variables used by MVA
  float tau_decay_mode = tau->decayMode();

  if (tau_decay_mode>11 || (tau_decay_mode>1&&tau_decay_mode<10)) return scores;

  Vector pi0;
  Vector pi;
  std::pair<Vector,Vector> rho;
  std::vector<Vector> a1_daughters = {};

  if(tau_decay_mode>=10) {
    std::pair<std::vector<Vector>, Vector>  a1 = GetA1(tau, gammas_pt_cut);
    a1_daughters  = a1.first;
    pi0 = a1.second;
  } else {
    for (auto g: tau->signalGammaCands()) if(g->pt()>gammas_pt_cut) gammas_.push_back(g);
    rho = GetRho (tau, gammas_pt_cut);
    pi0 = rho.second;
  }
 
  float strip_pt = pi0.pt();
  float E = tau->energy();

  float E1=-1;
  float E2=-1;
  float E3=-1;
  float a1_pi0_dEta=-1;
  float a1_pi0_dphi=-1;
  float a1_pi0_dEta_timesEtau=-1;
  float a1_pi0_dphi_timesEtau=-1;
  float h1_h2_dEta=-1;
  float h1_h2_dphi=-1;
  float h1_h3_dEta=-1;
  float h1_h3_dphi=-1;
  float h2_h3_dEta=-1;
  float h2_h3_dphi=-1;
  float h1_h2_dphi_timesE12=-1;
  float h1_h3_dphi_timesE13=-1;
  float h2_h3_dphi_timesE23=-1;
  float h1_h2_dEta_timesE12=-1;
  float h1_h3_dEta_timesE13=-1;
  float h2_h3_dEta_timesE23=-1;
  float mass0=-1;
  float mass1=-1;
  float mass2=-1;
  float strip_E=-1;
  float E1_overEa1=-1;
  float E2_overEa1=-1;
  float E1_overEtau=-1;
  float E2_overEtau=-1;
  float E3_overEtau=-1;

  if(tau_decay_mode>9 && a1_daughters.size()>2) {
    strip_E = pi0.energy();
    mass0 = (a1_daughters[0] + a1_daughters[1] + a1_daughters[2]).M();
    mass1 = (a1_daughters[0] + a1_daughters[1]).M();
    mass2 = (a1_daughters[0] + a1_daughters[2]).M();
    E1 = a1_daughters[0].energy();
    E2 = a1_daughters[1].energy();
    E3 = a1_daughters[2].energy();

    if(strip_pt>0) {
      a1_pi0_dEta = std::fabs(pi0.eta()-tau->eta());
      a1_pi0_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(pi0,tau->p4()));
    }

    h1_h2_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(a1_daughters[0],a1_daughters[1]));
    h1_h3_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(a1_daughters[0],a1_daughters[2]));
    h2_h3_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(a1_daughters[1],a1_daughters[2]));
    h1_h2_dEta = std::fabs(a1_daughters[0].eta()-a1_daughters[1].eta());
    h1_h3_dEta = std::fabs(a1_daughters[0].eta()-a1_daughters[2].eta());
    h2_h3_dEta = std::fabs(a1_daughters[1].eta()-a1_daughters[2].eta());
  }

  float Ea1 = E1+E2+E3;
  E1_overEa1 = E1/Ea1;
  E2_overEa1 = E2/Ea1;
  float Etau = Ea1+strip_E;
  E1_overEtau = E1/Etau;
  E2_overEtau = E2/Etau;
  E3_overEtau = E3/Etau;
  a1_pi0_dEta_timesEtau=a1_pi0_dEta*Etau;
  a1_pi0_dphi_timesEtau=a1_pi0_dphi*Etau;
  h1_h2_dphi_timesE12=h1_h2_dphi*(E1+E2);
  h1_h3_dphi_timesE13=h1_h3_dphi*(E1+E3);
  h2_h3_dphi_timesE23=h2_h3_dphi*(E2+E3);
  h1_h2_dEta_timesE12=h1_h2_dEta*(E1+E2);
  h1_h3_dEta_timesE13=h1_h3_dEta*(E1+E3);
  h2_h3_dEta_timesE23=h2_h3_dEta*(E2+E3);


  if (tau_decay_mode<12) {
    pi = rho.first; 
  }

  float Egamma1=-1, Egamma2=-1;
  float Epi = pi.energy();
  float Epi0 = pi0.energy();

  if(gammas_.size()>=1) Egamma1 = gammas_[0]->energy();
  if(gammas_.size()>=2) Egamma2 = gammas_[1]->energy();

  float Egamma1_tau = Egamma1/E;
  float Egamma2_tau = Egamma2/E;

  float Epi_tau = Epi/E;

  float pt = tau->pt();
  float eta = tau->eta();

  float rho_dEta=-1, rho_dphi=-1, gammas_dEta = -1., gammas_dphi = -1.;

  if(Epi0>0) {
    rho_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(pi,pi0));
    rho_dEta = std::fabs(pi.eta()-pi0.eta());
  }
  float rho_dEta_tau = rho_dEta*E;
  float rho_dphi_tau = rho_dphi*E;

  if(gammas_.size()>1) {
    gammas_dphi =  std::fabs(ROOT::Math::VectorUtil::DeltaPhi(gammas_[0]->p4(),gammas_[1]->p4()));
    gammas_dEta =  std::fabs(gammas_[0]->eta()-gammas_[1]->eta());
  }
  float gammas_dEta_tau = gammas_dEta* E;
  float gammas_dR_tau =  sqrt(gammas_dEta*gammas_dEta + gammas_dphi*gammas_dphi)*E;

  float Mpi0=-1, Mpi0_TwoHighGammas=-1;
  Vector gammas_vector;
  for (auto g : gammas_) gammas_vector+=g->p4();
  Mpi0 = gammas_vector.M();
  if(gammas_.size()>=2) Mpi0_TwoHighGammas = (gammas_[0]->p4() + gammas_[1]->p4()).M();

  float Mrho=-1, Mrho_OneHighGammas=-1, Mrho_TwoHighGammas=-1, Mrho_subleadingGamma=-1;
  Mrho = (pi + pi0).M();
  if(gammas_.size()>=1) Mrho_OneHighGammas=(pi + gammas_[0]->p4() ).M();
  if(gammas_.size()>=2) Mrho_TwoHighGammas=(pi + gammas_[0]->p4() + gammas_[1]->p4()).M();
  if(gammas_.size()>=2) Mrho_subleadingGamma=(pi + gammas_[1]->p4()).M();

  float DeltaR2WRTtau=-999;
  if(gammas_.size()>=1){
    DeltaR2WRTtau=0;
    double SumPt=0;
    DeltaR2WRTtau=std::pow(ROOT::Math::VectorUtil::DeltaR(pi,tau->p4()),2)*std::pow(pi.pt(),2);
    SumPt=std::pow(pi.pt(),2);
    for(auto g : gammas_){
      DeltaR2WRTtau+=std::pow(ROOT::Math::VectorUtil::DeltaR(g->p4(),tau->p4()),2)*std::pow(g->pt(),2);
      SumPt+=std::pow(g->pt(),2);
    }
    DeltaR2WRTtau/=SumPt;
  }
  float DeltaR2WRTtau_tau = DeltaR2WRTtau*E*E;

  // once the variables are computed they need to be stored in the order expected by TMVA reader
  std::vector<float> inputs;
  if(tau_decay_mode<2) {

    inputs.resize(24);

    vars_[0] = Egamma1_tau;
    vars_[1] = Egamma2_tau;
    vars_[2] = Epi_tau;
    vars_[3] = rho_dEta_tau;
    vars_[4] = rho_dphi_tau;
    vars_[5] = gammas_dEta_tau;
    vars_[6] = gammas_dR_tau;
    vars_[7] = DeltaR2WRTtau_tau;
    vars_[8] = tau_decay_mode;
    vars_[9] = eta;
    vars_[10] = pt;
    vars_[11] = Epi0;
    vars_[12] = Epi;
    vars_[13] = rho_dEta;
    vars_[14] = rho_dphi;
    vars_[15] = gammas_dEta;
    vars_[16] = Mrho;
    vars_[17] = Mpi0;
    vars_[18] = DeltaR2WRTtau;
    vars_[19] = Mpi0_TwoHighGammas;
    vars_[20] = Mrho_OneHighGammas;
    vars_[21] = Mrho_TwoHighGammas;
    vars_[22] = Mrho_subleadingGamma;
    vars_[23] = strip_pt;
  }
  if(tau_decay_mode>9) {

    inputs.resize(40);

    vars_dm10_[0] = E1_overEa1;
    vars_dm10_[1] = E2_overEa1;
    vars_dm10_[2] = E1_overEtau;
    vars_dm10_[3] = E2_overEtau;
    vars_dm10_[4] = E3_overEtau;
    vars_dm10_[5] = a1_pi0_dEta_timesEtau;
    vars_dm10_[6] = a1_pi0_dphi_timesEtau;
    vars_dm10_[7] = h1_h2_dphi_timesE12;
    vars_dm10_[8] = h1_h2_dEta_timesE12;
    vars_dm10_[9] = h1_h3_dphi_timesE13;
    vars_dm10_[10] = h1_h3_dEta_timesE13;
    vars_dm10_[11] = h2_h3_dphi_timesE23;
    vars_dm10_[12] = h2_h3_dEta_timesE23;
    vars_dm10_[13] = gammas_dEta_tau;
    vars_dm10_[14] = gammas_dR_tau;
    vars_dm10_[15] = tau_decay_mode;
    vars_dm10_[16] = mass0;
    vars_dm10_[17] = mass1;
    vars_dm10_[18] = mass2;
    vars_dm10_[19] = E1;
    vars_dm10_[20] = E2;
    vars_dm10_[21] = E3;
    vars_dm10_[22] = strip_E;
    vars_dm10_[23] = a1_pi0_dEta;
    vars_dm10_[24] = a1_pi0_dphi;
    vars_dm10_[25] = strip_pt;
    vars_dm10_[26] = pt;
    vars_dm10_[27] = eta;
    vars_dm10_[28] = E;
    vars_dm10_[29] = h1_h2_dphi;
    vars_dm10_[30] = h1_h3_dphi;
    vars_dm10_[31] = h2_h3_dphi;
    vars_dm10_[32] = h1_h2_dEta;
    vars_dm10_[33] = h1_h3_dEta;
    vars_dm10_[34] = h2_h3_dEta;
    vars_dm10_[35] = Egamma1;
    vars_dm10_[36] = Egamma2;
    vars_dm10_[37] = gammas_dEta;
    vars_dm10_[38] = Mpi0;
    vars_dm10_[39] = Mpi0_TwoHighGammas;
  }

  scores = read_mva_score((int)tau_decay_mode);

  return scores;
}


DEFINE_FWK_MODULE(PATTauDiscriminationMVADM);
}
