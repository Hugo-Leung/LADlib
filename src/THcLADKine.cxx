
#include "THcLADKine.h"
#include "THcGlobals.h"
#include "THcHallCSpectrometer.h"
#include "THcLADGEM.h"
#include "THcLADHodoscope.h"
#include "THcParmList.h"
#include "THcPrimaryKine.h"
#include "THcReactionPoint.h"
#include "TVector3.h"
#include "VarDef.h"
#include "VarType.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


ClassImp(THcLADKine)
    //_____________________________________________________________________________
    THcLADKine::THcLADKine(const char *name, const char *description, const char *spectro, const char *primary_kine,
                           const char *vertex_module)
    : THcPrimaryKine(name, description, spectro, primary_kine, 0.938), fSpecName(spectro), fGEM(nullptr),
      fHodoscope(nullptr), fVertexModuleName(vertex_module), fVertexModule(nullptr), fTrack(nullptr) {
  // Constructor
  fGoodLADHits    = new TClonesArray("THcGoodLADHit", MAXGOODHITS);
  goodhit_n       = 0;
  fTVertex        = kBig;
  fRFTime         = kBig;
  fTVertex_RFcorr = kBig;
  fFixed_z        = nullptr;
  rf_offset       = nullptr;
  n_rf_offsets    = 0;
  rf_period       = 4.00801; // Default RF period in ns

  fZCellMin  = -15.0;   // cm
  fZCellMax  = +15.0;   // cm
  fThetaMin  =  60.0*TMath::DegToRad();   // rad
  fThetaMax  = 170.0*TMath::DegToRad();   // rad
  fPhiMin    = -50.0*TMath::DegToRad();  // rad
  fPhiMax    = +50.0*TMath::DegToRad();  // rad
  fchisq_cut[0] = 20.0; // default chi2 cut for tracks with 2 hodo hits
  fchisq_cut[1] = 15.0; // default chi2 cut for tracks with 1 hodo hit
  fchisq_cut[2] = 10.0; // default chi2 cut for tracks with no hodo hits
  fSigma_GEM = 0.1; // default GEM resolution in cm
  fSigma_Hodo = 10; // default Hodoscope resolution in cm


}
//_____________________________________________________________________________
THcLADKine::~THcLADKine() {
  // Destructor
  DefineVariables(kDelete);
  fGoodLADHits->Delete();
  delete fGoodLADHits;
  delete[] fFixed_z;
  delete[] rf_offset;
}
//_____________________________________________________________________________
void THcLADKine::Clear(Option_t *opt) {
  // Clear the object
  THcPrimaryKine::Clear(opt);
  fGoodLADHits->Delete();
  goodhit_n = 0;
  for (Int_t i = 0; i < n_rf_offsets; i++) {
    rf_offset[i] = 0.0;
  }
  fTVertex        = kBig;
  fRFTime         = kBig;
  fTVertex_RFcorr = kBig;
}
//_____________________________________________________________________________
THaAnalysisObject::EStatus THcLADKine::Init(const TDatime &run_time) {

  // Initialize the object
  if (fSpectro == nullptr) {
    fSpectro = dynamic_cast<THcHallCSpectrometer *>(FindModule(fSpecName.Data(), "THcHallCSpectrometer"));
    if (fSpectro == nullptr) {
      Error("Init", "Cannot find spectrometer %s", fSpecName.Data());
      return kInitError;
    }
  }

  if (fHodoscope == nullptr) {
    if (fSpectro) {
      fHodoscope = dynamic_cast<THcLADHodoscope *>(fSpectro->GetDetector("ladhod"));
    }
    if (fHodoscope == nullptr) {
      Error("Init", "No Hodoscope module found in spectrometer");
      return kInitError;
    }
  }
  if (fGEM == nullptr) {
    if (fSpectro) {
      fGEM = dynamic_cast<THcLADGEM *>(fSpectro->GetDetector("gem"));
    }
    if (fGEM == nullptr) {
      Error("Init", "No GEM module found in spectrometer");
      return kInitError;
    }
  }

  if (fVertexModuleName != "") {
    fVertexModule = static_cast<THcReactionPoint *>(FindModule(fVertexModuleName, "THcReactionPoint"));
    if (fVertexModule == nullptr) {
      Error("Init", "No Vertex module found");
      return kInitError;
    }
  }

  TString fTrigDetName;
  // Get the apparatus prefix
  char aparatus_prefix[2];
  aparatus_prefix[0] = tolower(fSpecName[0]);
  aparatus_prefix[1] = '\0';
  if (aparatus_prefix[0] == 'h') {
    fTrigDetName = "T.hms";
  } else if (aparatus_prefix[0] == 'p') {
    fTrigDetName = "T.shms";
  } else {
    Error("Init", "Invalid apparatus prefix: %s", aparatus_prefix);
    return kInitError;
  }

  fTrigDet = dynamic_cast<THcTrigDet *>(FindModule(fTrigDetName.Data(), "THcTrigDet"));
  if (!fTrigDet) {
    cout << "THcCoinTime module  Cannnot find TrigDet = " << fTrigDetName.Data() << endl;
    fStatus = kInitError;
    return fStatus;
  }

  return THaPhysicsModule::Init(run_time);
}
//_____________________________________________________________________________
Int_t THcLADKine::ReadDatabase(const TDatime &date) {
  // Read the database
  Int_t err = THcPrimaryKine::ReadDatabase(date);
  if (err)
    return err;

  char prefix[2];
  prefix[0] = 'l';
  prefix[1] = '\0';

  fD0Cut_wVertex        = 0.0;
  fD0Cut_noVertex       = 0.0;
  fMax_dTrk_horiz_match = 0.0;
  fMax_dTrk_vert_match  = 0.0;
  fNfixed_z             = 0;
  fglobal_time_offset   = 0.0;
  fTrk_dtCut            = 10.0;

  cout << "Reading LAD Kinematics parameters from database..." << endl;

  DBRequest list[] = {{"d0_cut_wVertex", &fD0Cut_wVertex, kDouble, 0, 1},
                      {"d0_cut_noVertex", &fD0Cut_noVertex, kDouble, 0, 1},
                      {"max_dTrk_horiz_hitMatch", &fMax_dTrk_horiz_match, kDouble, 0, 1},
                      {"max_dTrk_vert_hitMatch", &fMax_dTrk_vert_match, kDouble, 0, 1},
                      {"nfixed_z", &fNfixed_z, kInt, 0, 1},
                      {"global_time_offset", &fglobal_time_offset, kDouble, 0, 1},
                      {"trk_dt_cut", &fTrk_dtCut, kDouble, 0, 1},
                      {"_rf_period", &rf_period, kDouble, 0, 1},
                      {0}};
  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  DBRequest list_rf[] = {{"_nr_rf_offsets", &n_rf_offsets, kInt, 0, 1}, {0}};
  gHcParms->LoadParmValues((DBRequest *)&list_rf, prefix);
  rf_offset = new Double_t[n_rf_offsets * 3];
  for (Int_t i = 0; i < n_rf_offsets * 3; i++) {
    rf_offset[i] = 0.0;
  }
  DBRequest list_rf_offsets[] = {{"_rf_offset", &rf_offset[0], kDouble, (UInt_t)(n_rf_offsets * 3)}, {0}};
  gHcParms->LoadParmValues((DBRequest *)&list_rf_offsets, prefix);

  delete[] fFixed_z;
  if (fNfixed_z > 0) {
    fFixed_z          = new Double_t[fNfixed_z];
    DBRequest list2[] = {{"fixed_z_pos", fFixed_z, kDouble, fNfixed_z}, {0}};
    gHcParms->LoadParmValues((DBRequest *)&list2, prefix);
  } else {
    fFixed_z = nullptr;
  }

  DBRequest track_constraints[] = {{ "ladkin.z_cell_min", &fZCellMin, kDouble, 0, 1 },
                                   { "ladkin.z_cell_max", &fZCellMax, kDouble, 0, 1 },
                                   { "ladkin.theta_min", &fThetaMin, kDouble, 0, 1 },
                                   { "ladkin.theta_max", &fThetaMax, kDouble, 0, 1 },
                                   { "ladkin.phi_min", &fPhiMin, kDouble, 0, 1 },
                                   { "ladkin.phi_max", &fPhiMax, kDouble, 0, 1 },
                                   { "ladkin.chisq_cut", fchisq_cut, kDouble, 3, 1 },
                                   { 0 } };
  gHcParms->LoadParmValues((DBRequest *)&track_constraints, prefix);

  if (fZCellMin > fZCellMax) {
    Double_t tmp = fZCellMin; fZCellMin = fZCellMax; fZCellMax = tmp;
  }
  if (fThetaMin > fThetaMax) {
    Double_t tmp = fThetaMin; fThetaMin = fThetaMax; fThetaMax = tmp;
  }
  if (fPhiMin > fPhiMax) {
    Double_t tmp = fPhiMin; fPhiMin = fPhiMax; fPhiMax = tmp;
  }

  fThetaMin = TMath::Max(0.0,   TMath::Min(TMath::Pi(), fThetaMin));
  fThetaMax = TMath::Max(0.0,   TMath::Min(TMath::Pi(), fThetaMax));
  fPhiMin   = TMath::Max(-TMath::Pi(),TMath::Min(TMath::Pi(), fPhiMin));
  fPhiMax   = TMath::Max(-TMath::Pi(),TMath::Min(TMath::Pi(), fPhiMax));


  return err;
}
//_____________________________________________________________________________
Int_t THcLADKine::Process(const THaEvData &evdata) {
  // Get Golden Track (for tof calculation later)
  fTrack = fSpectro->GetGoldenTrack();
  CalculateTVertex();

  //Get LADGoodHits for track matching 
  TClonesArray *LADHits_unfiltered = fHodoscope->GetLADGoodHits();
  Int_t nHits                      = LADHits_unfiltered->GetLast() + 1;
  goodhit_n = 0;
  //////////////////////////////////////////////////////////////////////////////
  // Check track projection to vertex
  fGEMTracks    = fGEM->GetTracks();
  Int_t ntracks = fGEMTracks->GetLast() + 1;
  std::vector<bool> isGoodTrack(ntracks, false);
  for (Int_t iTrack = 0; iTrack < ntracks; iTrack++) {
    THcLADGEMTrack *track = static_cast<THcLADGEMTrack *>(fGEMTracks->At(iTrack));
    if (track == nullptr)
      continue;

    // Get the hit positions
    TVector3 v_hit1(track->GetX1(), track->GetY1(), track->GetZ1());
    TVector3 v_hit2(track->GetX2(), track->GetY2(), track->GetZ2());

    TVector3 vertex;
    if (fVertexModule->HasVertex()) {
      vertex = fVertexModule->GetVertex();
      track->SetGoodD0(kTRUE);
    } else {
      vertex.SetXYZ(0, 0, 0);
    }
    
    // Fix track vertex (for improved resolution on multifoils)
    if (fNfixed_z > 0 && track->GetGoodD0()) {
      std::vector<double> distances;
      for (Int_t j = 0; j < fNfixed_z; j++) {
        distances.push_back(abs(vertex.Z() - fFixed_z[j]));
      }
      auto min_it     = std::min_element(distances.begin(), distances.end());
      Int_t min_index = std::distance(distances.begin(), min_it);
      vertex.SetZ(fFixed_z[min_index]);
    }

    Double_t gemdir[3];
    gemdir[0] = TMath::ACos((v_hit2.Z() - v_hit1.Z()) / (v_hit2 - v_hit1).Mag()) ;
    gemdir[1] = TMath::ATan2((v_hit2.Y() - v_hit1.Y()) , (v_hit2.X() - v_hit1.X()));
    gemdir[2] = vertex.Z(); //Hopefully the GEM track will be the same as elecctron vertex Z

    Double_t bestchisq[2]={kBig,kBig};// 2 hodo hits, 1 hodo hit 
    Double_t  best_gemdir[2][3];
    int bestHodoHitIndex[2]={-1,-1};
    //loop over hodoscope hits to find the best match to the track projection
    for (int iHodoHit=0; iHodoHit<nHits; iHodoHit++) {
      THcGoodLADHit *hodo_hit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(iHodoHit));
      if (hodo_hit == nullptr)
        continue;
      // There might be one or two hits with a signle goodhit, need to check both
      std::vector<TVector3> v_hits;
      std::vector<Double_t> v_resolutions;
      v_hits.push_back(v_hit1);
      v_hits.push_back(v_hit2);
      v_resolutions.push_back(fSigma_GEM);
      v_resolutions.push_back(fSigma_GEM);
      int num_hodo_hits = 0;
      TVector3 v_hodo_hit;
      if (hodo_hit->GetPlaneHit0() >= 0) {
        v_hodo_hit = fHodoscope->GetHitPositionLab(hodo_hit->GetPlaneHit0(), hodo_hit->GetPaddleHit0(), hodo_hit->GetHitYPosHit0());
        v_hits.push_back(v_hodo_hit);
        v_resolutions.push_back(fSigma_Hodo);
        num_hodo_hits++;
      } 
      if (hodo_hit->GetPlaneHit1() >= 0) {
        v_hodo_hit = fHodoscope->GetHitPositionLab(hodo_hit->GetPlaneHit1(), hodo_hit->GetPaddleHit1(), hodo_hit->GetHitYPosHit1());
        v_hits.push_back(v_hodo_hit);
        v_resolutions.push_back(fSigma_Hodo);
        num_hodo_hits++;
      }
      Double_t dir[3];
      dir[0]=gemdir[0];
      dir[1]=gemdir[1];
      dir[2]=gemdir[2];
      Double_t chisq = FitTrack(vertex, v_hits, v_resolutions, dir);
      // Check if this is the best match so far, we only want to keep the best hodoscope hit for each track to avoid double counting hits
      if (chisq>=0&&chisq<bestchisq[2-num_hodo_hits]) {
        bestchisq[2-num_hodo_hits]=chisq;
        best_gemdir[2-num_hodo_hits][0]=dir[0];
        best_gemdir[2-num_hodo_hits][1]=dir[1];
        best_gemdir[2-num_hodo_hits][2]=dir[2];
        bestHodoHitIndex[2-num_hodo_hits]=iHodoHit;
      }
      v_hits.clear();
      v_resolutions.clear();
    }
    if(bestchisq[0]>=0&&bestchisq[0]<fchisq_cut[0]) {
      isGoodTrack[iTrack] = true;
      track->SetAngles(best_gemdir[0][0], best_gemdir[0][1]);
      track->SetProjVertex(vertex.X(), vertex.Y(), best_gemdir[0][2]);
      track->SetChisq(bestchisq[0]);
      track->SetGoodD0(kTRUE);
      track->SetD0(TMath::Abs(best_gemdir[0][2]-vertex.Z())*TMath::Sin(best_gemdir[0][0]));
      THcGoodLADHit *bestHodoHit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(bestHodoHitIndex[0]));
      track->SetBestHodoHit(bestHodoHit);
      track->SetHasHodoHit(2);
      THcGoodLADHit *goodhit = new ((*fGoodLADHits)[goodhit_n]) THcGoodLADHit();
      goodhit->CopyHit(0,0,bestHodoHit);
      goodhit->CopyHit(1,1,bestHodoHit);
      goodhit->SetTrackID(0, iTrack);
      goodhit->SetTrackID(1, iTrack);
      goodhit_n++;
    }else if(bestchisq[1]>=0&&bestchisq[1]<fchisq_cut[1]) {
      isGoodTrack[iTrack] = true;
      track->SetAngles(best_gemdir[1][0], best_gemdir[1][1]);
      track->SetProjVertex(vertex.X(), vertex.Y(), best_gemdir[1][2]);
      track->SetChisq(bestchisq[1]);
      track->SetGoodD0(kTRUE);
      track->SetD0(TMath::Abs(best_gemdir[1][2]-vertex.Z())*TMath::Sin(best_gemdir[1][0]));
      THcGoodLADHit *bestHodoHit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(bestHodoHitIndex[1]));
      track->SetBestHodoHit(bestHodoHit);
      track->SetHasHodoHit(1);
      THcGoodLADHit *goodhit = new ((*fGoodLADHits)[goodhit_n]) THcGoodLADHit();
      goodhit->CopyHit(0,0,bestHodoHit);
      goodhit->CopyHit(1,1,bestHodoHit);
      goodhit->SetTrackID(0, iTrack);
      goodhit->SetTrackID(1, iTrack);
      goodhit_n++;
    }else {
      Double_t dir[3];
      dir[0]=gemdir[0];
      dir[1]=gemdir[1];
      dir[2]=gemdir[2];// if the track doesn't have a good hodoscope match, just fit to the GEM hits and vertex
      Double_t chisq= FitTrack(vertex, {v_hit1, v_hit2}, {fSigma_GEM, fSigma_GEM}, dir);
      track->SetChisq(chisq);// set chisq even if the fit fails so we can see why
      if (chisq>=0&&chisq<fchisq_cut[2]) {
        track->SetGoodD0(kTRUE);
        isGoodTrack[iTrack] = true;
        track->SetAngles(dir[0], dir[1]);
        track->SetProjVertex(vertex.X(), vertex.Y(), dir[2]);
        track->SetD0(TMath::Abs(dir[2]-vertex.Z())*TMath::Sin(dir[0]));
        track->SetBestHodoHit(nullptr);
        track->SetHasHodoHit(0);
      }else {
        //Still bad track, mark as bad and set parameters to default values
        isGoodTrack[iTrack] = false;
        //cout<<"THcLADKine: Track fit failed for track "<<iTrack<<" with error code "<<chisq<<endl;
        track->SetAngles(-kBig, -kBig);
        track->SetProjVertex(-kBig, -kBig, -kBig);
        track->SetD0(-kBig);
        track->SetGoodD0(kFALSE);
        track->SetHasHodoHit(0);
        track->SetBestHodoHit(nullptr);
        track->SetIsGoodTrack(false);
        continue;
      }
    }

    if (track->GetGoodD0()) {
      if (track->GetD0()>0 && track->GetD0() < fD0Cut_wVertex) {
        isGoodTrack[iTrack] = isGoodTrack[iTrack] && true;
      }
    } else {
      if (track->GetD0() > 0 && track->GetD0() < fD0Cut_noVertex) {
        isGoodTrack[iTrack] = isGoodTrack[iTrack] && true;
      }
    }
    if(track->GetChisq()<0) {
      isGoodTrack[iTrack] = false;
    }
    if (track->GetdT() > fTrk_dtCut) {
      isGoodTrack[iTrack] = false;
    }
    track->SetIsGoodTrack(isGoodTrack[iTrack]);
  }
  // Calculate beta, alpha, tof, etc. for the hits
  for (Int_t i = 0; i < goodhit_n; i++) {
    THcGoodLADHit *goodhit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(i));
    if (goodhit == nullptr)
      continue;

    Double_t tof, tof_rfcorr, beta, alpha;
    if (goodhit->GetPlaneHit0() >= 0 && goodhit->GetPlaneHit0() < 999) {
      tof = CalculateToF(goodhit->GetHitTimeHit0());
      goodhit->SetHitTOF(0, tof);
      tof_rfcorr = CalculateTOFRFcorr(goodhit->GetHitTimeHit0());
      goodhit->SetHitTOFRFcorr(0, tof_rfcorr);
      // Calculate beta, alpha, etc. for the first plane
    }
    if (goodhit->GetPlaneHit1() >= 0 && goodhit->GetPlaneHit1() < 999) {
      tof = CalculateToF(goodhit->GetHitTimeHit1());
      goodhit->SetHitTOF(1, tof);
      tof_rfcorr = CalculateTOFRFcorr(goodhit->GetHitTimeHit1());
      goodhit->SetHitTOFRFcorr(1, tof_rfcorr);
      // Calculate beta, alpha, etc. for the second plane
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Calculate tof for raw planes
  // Get the number of planes
  Int_t nPlanes = fHodoscope->GetNPlanes();

  // Loop through the planes and get LAD hits
  for (Int_t plane = 0; plane < nPlanes; plane++) {
    TClonesArray *planeHits = fHodoscope->GetLADHits(plane);
    Int_t nPlaneHits        = planeHits->GetLast() + 1;

    for (Int_t i = 0; i < nPlaneHits; i++) {
      THcLADHodoHit *hit = static_cast<THcLADHodoHit *>(planeHits->At(i));
      if (hit == nullptr)
        continue;

      // Calculate time of flight for each hit
      Double_t tof_avg = CalculateToF(hit->GetScinCorrectedTime());
      Double_t tof_top = CalculateToF(hit->GetTopCorrectedTime());
      Double_t tof_btm = CalculateToF(hit->GetBtmCorrectedTime());
      hit->SetTOFCorrectedTimes(tof_top, tof_btm, tof_avg);

      tof_avg = CalculateTOFRFcorr(hit->GetScinCorrectedTime());
      tof_top = CalculateTOFRFcorr(hit->GetTopCorrectedTime());
      tof_btm = CalculateTOFRFcorr(hit->GetBtmCorrectedTime());
      hit->SetTOF_RF_CorrectedTimes(tof_top, tof_btm, tof_avg);
    }
  }

  // LADHits_unfiltered->Clear(); // Clear the unfiltered hits
  // fGEMTracks->Clear();      // Clear the tracks
  return kOK;
}
//_____________________________________________________________________________
void THcLADKine::CalculateTVertex() {
  // Calculate the time of flight
  if (fTrack == nullptr) {
    // Error("CalculateToF", "No track found");
    return;
  }

  Double_t HMScentralPathLen  = 22.0 * 100.;
  Double_t SHMScentralPathLen = 18.1 * 100.;
  Double_t elecMass           = 0.510998 / 1000.0; // electron mass in GeV/c^2
  Double_t lightSpeed         = 29.9792;           // in cm/ns

  Double_t fPtime = fTrack->GetFPTime();
  Double_t xptar  = fTrack->GetTTheta();
  Double_t dP     = fTrack->GetDp();
  Double_t elec_P = fTrack->GetP();
  if (fPtime == -2000 || fPtime == -1000) {
    return;
  }

  // Leave these here for now, but don't think I need them
  // Double_t pSHMS_TdcTime_ROC1 = fTrigDet->Get_CT_Trigtime(0); // SHMS
  // Double_t pHMS_TdcTime_ROC1  = fTrigDet->Get_CT_Trigtime(1); // HMS
  // Double_t pSHMS_TdcTime_ROC2 = fTrigDet->Get_CT_Trigtime(2); // SHMS pTrig1
  // Double_t pHMS_TdcTime_ROC2  = fTrigDet->Get_CT_Trigtime(3); // HMS pTrig3

  char aparatus_prefix[2];
  aparatus_prefix[0] = tolower(fSpecName[0]);
  aparatus_prefix[1] = '\0';

  Double_t DeltaPathLength;
  if (aparatus_prefix[0] == 'h') {
    DeltaPathLength = (.12 * xptar * 1000 + 0.17 * dP / 100.);
    DeltaPathLength += HMScentralPathLen;
  } else if (aparatus_prefix[0] == 'p') {
    DeltaPathLength = (.11 * xptar * 1000 + 0.057 * dP / 100.);
    DeltaPathLength += SHMScentralPathLen;
  } else {
    Error("CalculateToF", "Invalid apparatus prefix: %s", aparatus_prefix);
    return;
  }

  Double_t beta_calc      = elec_P / sqrt(elec_P * elec_P + elecMass * elecMass);
  Double_t PathLengthCorr = DeltaPathLength / (lightSpeed * beta_calc);

  fTVertex = fPtime - PathLengthCorr;

  THaVar *varptr;
  varptr = gHcParms->Find("gen_run_number");
  Int_t *runnum;
  if (varptr) {
    runnum = (Int_t *)varptr->GetValuePointer(); // Assume correct type
  }
  int fRFTimeIndex = -1;
  int run_idx      = 0;
  while (fRFTimeIndex < 0) {
    if (run_idx >= 3 * n_rf_offsets || rf_offset[run_idx] > *runnum) {
      fRFTimeIndex = run_idx - 3;
    }
    run_idx += 3;
  }
  fRFTimeIndex += (aparatus_prefix[0] == 'p' ? 1 : 2);
  // (aparatus_prefix[0] == 'p') ? 0 : 1;

  fRFTime    = fTrigDet->Get_RF_TrigTime(fRFTimeIndex) + rf_offset[fRFTimeIndex];
  double tmp = fmod(fTVertex - fRFTime, rf_period);
  if (tmp > 2)
    tmp -= rf_period;
  fTVertex_RFcorr = fTVertex - tmp;
  cout << "Run: " << *runnum << ", RF time: " << fRFTime << ", Vertex time (uncorrected): " << fTVertex
       << ", Vertex time (RF corrected): " << fTVertex_RFcorr << endl;

  return;
}
//_____________________________________________________________________________
Double_t THcLADKine::CalculateToF(Double_t t_raw) {

  Double_t tof = t_raw - fTVertex + fglobal_time_offset;

  return tof;
}
Double_t THcLADKine::CalculateTOFRFcorr(Double_t t_raw) {

  Double_t tof = t_raw - fTVertex_RFcorr + fglobal_time_offset;

  return tof;
}
//_____________________________________________________________________________
Int_t THcLADKine::DefineVariables(EMode mode) {
  if (mode == kDefine && fIsSetup)
    return kOK;
  fIsSetup = (mode == kDefine);

  // Define variables for the analysis tree
  if (mode == kDefine) {
    RVarDef vars[] = {
        {"t_vertex", "Calculated vertex time", "fTVertex"},
        {"t_vertex_RFcorr", "Calculated vertex time (RF corrected)", "fTVertex_RFcorr"},
        {"n_goodhits", "Number of good hits", "goodhit_n"},
        {"goodhit_plane_0", "Good hit plane", "fGoodLADHits.THcGoodLADHit.GetPlaneHit0()"},
        {"goodhit_paddle_0", "Good hit paddle", "fGoodLADHits.THcGoodLADHit.GetPaddleHit0()"},
        {"goodhit_trackid_0", "Good hit track ID", "fGoodLADHits.THcGoodLADHit.GetTrackIDHit0()"},
        {"goodhit_beta_0", "Good hit beta", "fGoodLADHits.THcGoodLADHit.GetBetaHit0()"},
        {"goodhit_dTrkHoriz_0", "Good hit horizontal trk proj - hit position",
         "fGoodLADHits.THcGoodLADHit.GetdTrkHorizHit0()"},
        {"goodhit_dTrkVert_0", "Good hit vertical trk proj - hit position",
         "fGoodLADHits.THcGoodLADHit.GetdTrkVertHit0()"},
        {"goodhit_hittime_0", "Good hit time", "fGoodLADHits.THcGoodLADHit.GetHitTimeHit0()"},
        {"goodhit_hittheta_0", "Good hit theta", "fGoodLADHits.THcGoodLADHit.GetHitThetaHit0()"},
        {"goodhit_hitphi_0", "Good hit phi", "fGoodLADHits.THcGoodLADHit.GetHitPhiHit0()"},
        {"goodhit_hitedep_0", "Good hit energy deposition", "fGoodLADHits.THcGoodLADHit.GetHitEdepHit0()"},
        {"goodhit_hitedep_amp_0", "Good hit energy deposition (amplitude)",
         "fGoodLADHits.THcGoodLADHit.GetHitEdepAmpHit0()"},
        {"goodhit_plane_1", "Good hit plane (second plane)", "fGoodLADHits.THcGoodLADHit.GetPlaneHit1()"},
        {"goodhit_paddle_1", "Good hit paddle (second plane)", "fGoodLADHits.THcGoodLADHit.GetPaddleHit1()"},
        {"goodhit_trackid_1", "Good hit track ID (second plane)", "fGoodLADHits.THcGoodLADHit.GetTrackIDHit1()"},
        {"goodhit_beta_1", "Good hit beta (second plane)", "fGoodLADHits.THcGoodLADHit.GetBetaHit1()"},
        {"goodhit_dTrkHoriz_1", "Good hit horizontal trk proj - hit position (second plane)",
         "fGoodLADHits.THcGoodLADHit.GetdTrkHorizHit1()"},
        {"goodhit_dTrkVert_1", "Good hit vertical trk proj - hit position (second plane)",
         "fGoodLADHits.THcGoodLADHit.GetdTrkVertHit1()"},
        {"goodhit_hittime_1", "Good hit time (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitTimeHit1()"},
        {"goodhit_hittheta_1", "Good hit theta (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitThetaHit1()"},
        {"goodhit_hitphi_1", "Good hit phi (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitPhiHit1()"},
        {"goodhit_hitedep_1", "Good hit energy deposition (second plane)",
         "fGoodLADHits.THcGoodLADHit.GetHitEdepHit1()"},
        {"goodhit_hitedep_amp_1", "Good hit energy deposition (amplitude) (second plane)",
         "fGoodLADHits.THcGoodLADHit.GetHitEdepAmpHit1()"},
        {"goodhit_ypos_0", "Good hit y position", "fGoodLADHits.THcGoodLADHit.GetHitYPosHit0()"},
        {"goodhit_ypos_1", "Good hit y position (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitYPosHit1()"},
        {"goodhit_tof_0", "Good hit time of flight", "fGoodLADHits.THcGoodLADHit.GetHitTOFHit0()"},
        {"goodhit_tof_1", "Good hit time of flight (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitTOFHit1()"},
        {"goodhit_alpha_0", "Good hit alpha", "fGoodLADHits.THcGoodLADHit.GetHitAlphaHit0()"},
        {"goodhit_alpha_1", "Good hit alpha (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitAlphaHit1()"},
        {0}};
    return DefineVarsFromList(vars, mode);
  }
  return kOK;
}
//_____________________________________________________________________________

Double_t THcLADKine::FitTrack(TVector3 vertex, std::vector<TVector3> sp_positions, std::vector<double> sp_resolutions, double dir[3]){
  if (dir == nullptr || sp_positions.empty() || sp_resolutions.empty()) {
    return -1; // Invalid input
  }
  int nPoints = sp_positions.size();
  if (nPoints < 2 || nPoints > 4) {
    return -1; // Not enough points to fit a track
  }
  if (nPoints != sp_resolutions.size()) {
    return -2; // Mismatch in number of points and resolutions
  }
  for(int i = 0; i < nPoints; i++) {
    if (sp_resolutions[i] <= 0) {
      return -3; // Invalid resolution value
    }
  }
  if (dir[0] < 0 || dir[0] > TMath::Pi() || dir[1] < -TMath::Pi() || dir[1] > TMath::Pi()) {
    return -4; // Invalid initial direction values
  }
  //check if the gem track is pointing to the hodoscope hits if there are any, return negative chisq if not
  if (nPoints >2){
    double sx = TMath::Sin(dir[0]) * TMath::Cos(dir[1]);
    double sy = TMath::Sin(dir[0]) * TMath::Sin(dir[1]);
    double sz = TMath::Cos(dir[0]);
    for (int i = 2; i < nPoints; i++) {
      //closest approach of the track to the point
      double t = ((sp_positions[i].X() - sp_positions[0].X()) * sx +
                  (sp_positions[i].Y() - sp_positions[0].Y()) * sy +
                  (sp_positions[i].Z() - sp_positions[0].Z()) * sz);
      double x_closest = sp_positions[0].X() + t * sx;
      double y_closest = sp_positions[0].Y() + t * sy;
      double z_closest = sp_positions[0].Z() + t * sz;
      double dx = sp_positions[i].X() - x_closest;
      double dy = sp_positions[i].Y() - y_closest;
      double dz = sp_positions[i].Z() - z_closest;
      double dist2 = dx*dx +  dz*dz;
      if (dist2 > (22*22)) { // if the track is more than 22cm away from the hodoscope hit, return negative chisq
        return -5; // Track does not point to the hodoscope hit
      }
    }
  }


  //requireing the track to originate from vertex (x,y), and only fitting for the track direction (theta, phi) and z vertex position. 
  double chi2 = -kBig;
  ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimizer->SetMaxFunctionCalls(1000);
  minimizer->SetTolerance(1e-6);
  ROOT::Math::Functor f([=](const double* params) { 
    double theta = params[0]; //  in radians
    double phi   = params[1]; //  in radians
    double z     = params[2];
    double chi2_local = 0;
    double sx = TMath::Sin(theta) * TMath::Cos(phi);
    double sy = TMath::Sin(theta) * TMath::Sin(phi);
    double sz = TMath::Cos(theta);
    for (int i = 0; i < nPoints; i++) {
      //closest approach of the track to the point
      double t = ((sp_positions[i].X() - vertex.X()) * sx +
                  (sp_positions[i].Y() - vertex.Y()) * sy +
                  (sp_positions[i].Z() - z) * sz);
      double x_closest = vertex.X() + t * sx;
      double y_closest = vertex.Y() + t * sy;
      double z_closest = z + t * sz;
      double dx = sp_positions[i].X() - x_closest;
      double dy = sp_positions[i].Y() - y_closest;
      double dz = sp_positions[i].Z() - z_closest;
      double dist2 = dx*dx +dz*dz;
      if (i>2){
        // Assume these are hodoscope hit, so no chisq penalty if the hit is within the width of the paddle
        double paddle_width = 22; // in cm, TODO: get actual paddle width from database
        if (dist2 < (paddle_width/2.0)*(paddle_width/2.0)) {
          dist2 = 0;
        }
      }
      chi2_local += (dist2+ dy*dy) / (sp_resolutions[i]*sp_resolutions[i]);//do we want to weight the vertical and horizontal residuals differently based on detector performance?
    }
    return chi2_local;
  },3);
  minimizer->SetFunction(f);

  double initial_params[3] = {dir[0], dir[1], vertex.Z()};
  minimizer->SetLimitedVariable(0, "theta", initial_params[0], 0.1, fThetaMin, fThetaMax);// Limit theta  to avoid unphysical solutions (tracks going backwards)
  minimizer->SetLimitedVariable(1, "phi", initial_params[1], 0.1, fPhiMin, fPhiMax);//Limit phi as the GEMs are only on the left side of the target in beam direction
  minimizer->SetLimitedVariable(2, "z_vertex", initial_params[2], 0.1, fZCellMin, fZCellMax); //the target length is 20cm
  minimizer->Minimize();
  if (minimizer->Status() != 0) {
    return -6; // Fit did not converge
  }
  const double *best_params = minimizer->X();
  dir[0] = best_params[0]; // theta in radians
  dir[1] = best_params[1]; // phi in radians
  dir[2] = best_params[2]; // z vertex position
  chi2 = minimizer->MinValue();
  //clean up
  delete minimizer;

  //return the chi2 of the fit
  return chi2;
}

