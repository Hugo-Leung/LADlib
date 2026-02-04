
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
  rf_offset[0]    = 0.0;
  rf_offset[1]    = 0.0;
}
//_____________________________________________________________________________
THcLADKine::~THcLADKine() {
  // Destructor
  DefineVariables(kDelete);
  fGoodLADHits->Delete();
  delete fGoodLADHits;
  delete[] fFixed_z;
}
//_____________________________________________________________________________
void THcLADKine::Clear(Option_t *opt) {
  // Clear the object
  THcPrimaryKine::Clear(opt);
  fGoodLADHits->Delete();
  goodhit_n = 0;
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

  DBRequest list[] = {{"d0_cut_wVertex", &fD0Cut_wVertex, kDouble, 0, 1},
                      {"d0_cut_noVertex", &fD0Cut_noVertex, kDouble, 0, 1},
                      {"max_dTrk_horiz_hitMatch", &fMax_dTrk_horiz_match, kDouble, 0, 1},
                      {"max_dTrk_vert_hitMatch", &fMax_dTrk_vert_match, kDouble, 0, 1},
                      {"nfixed_z", &fNfixed_z, kInt, 0, 1},
                      {"global_time_offset", &fglobal_time_offset, kDouble, 0, 1},
                      {"trk_dt_cut", &fTrk_dtCut, kDouble, 0, 1},
                      {"_rf_offset", rf_offset, kDouble, 2},
                      {0}};
  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  delete[] fFixed_z;
  if (fNfixed_z > 0) {
    fFixed_z          = new Double_t[fNfixed_z];
    DBRequest list2[] = {{"fixed_z_pos", fFixed_z, kDouble, fNfixed_z}, {0}};
    gHcParms->LoadParmValues((DBRequest *)&list2, prefix);
  } else {
    fFixed_z = nullptr;
  }

  return err;
}
//_____________________________________________________________________________
Int_t THcLADKine::Process(const THaEvData &evdata) {
  // Get Golden Track (for tof calculation later)
  fTrack = fSpectro->GetGoldenTrack();
  CalculateTVertex();

  //////////////////////////////////////////////////////////////////////////////
  // Check track projection to vertex
  fGEMTracks    = fGEM->GetTracks();
  Int_t ntracks = fGEMTracks->GetLast() + 1;
  std::vector<bool> isGoodTrack(ntracks, false);
  for (Int_t i = 0; i < ntracks; i++) {
    THcLADGEMTrack *track = static_cast<THcLADGEMTrack *>(fGEMTracks->At(i));
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

    // Re-calculate d0
    double numer = ((vertex - v_hit1).Cross((vertex - v_hit2))).Mag();
    double denom = (v_hit2 - v_hit1).Mag();
    // here we can put a range/fiducial cut on d0 taking into account the target size
    double d0 = numer / denom;
    track->SetD0(d0);

    // Re-calculate vpz
    double vpz = v_hit1.Z() + (vertex.X() - v_hit1.X()) * (v_hit2.Z() - v_hit1.Z()) / (v_hit2.X() - v_hit1.X());
    track->SetZVertex(vpz);

    // Re-calculate vpy
    double vpy = v_hit1.Y() + (vertex.X() - v_hit1.X()) * (v_hit2.Y() - v_hit1.Y()) / (v_hit2.X() - v_hit1.X());
    track->SetYVertex(vpy);

    if (track->GetGoodD0()) {
      if (d0 < fD0Cut_wVertex) {
        isGoodTrack[i] = true;
      }
    } else {
      if (d0 < fD0Cut_noVertex) {
        isGoodTrack[i] = true;
      }
    }

    if (track->GetdT() > fTrk_dtCut) {
      isGoodTrack[i] = false;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Remove hits with bad tracks and duplicate hits
  TClonesArray *LADHits_unfiltered = fHodoscope->GetLADGoodHits();
  Int_t nHits                      = LADHits_unfiltered->GetLast() + 1;
  for (Int_t i = 0; i < nHits; i++) {
    THcGoodLADHit *hit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(i));
    if (hit == nullptr)
      continue;

    // Check if the hit has a track that points back to the origin.
    Int_t track_id = hit->GetTrackIDHit0();
    if (track_id < 0 || track_id >= ntracks || !isGoodTrack[track_id])
      continue;

    // Get the track parameters
    Int_t plane  = hit->GetPlaneHit0();
    Int_t paddle = hit->GetPaddleHit0();

    Int_t plane_loc; // Front vs back plane
    if (plane == 0 || plane == 2 || plane == 4)
      plane_loc = 0;
    else
      plane_loc = 1;

    Int_t matching_hit_index = -1;
    Int_t partner_hit_index  = -1;
    for (Int_t j = 0; j < goodhit_n; j++) {
      THcGoodLADHit *goodhit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(j));
      if (goodhit == nullptr)
        continue;

      // Check if the hit is already in the list
      Int_t good_hit_plane  = plane_loc ? goodhit->GetPlaneHit1() : goodhit->GetPlaneHit0();
      Int_t good_hit_paddle = plane_loc ? goodhit->GetPaddleHit1() : goodhit->GetPaddleHit0();
      if (good_hit_plane == plane && good_hit_paddle == paddle) {
        matching_hit_index = j;
        break;
      }
      // Check if the hit has a partner hit
      Int_t partner_hit_track_id = plane_loc ? goodhit->GetTrackIDHit1() : goodhit->GetTrackIDHit0();
      if (partner_hit_track_id == track_id) {
        partner_hit_index = j;
      }
    }
    if (matching_hit_index == -1 && partner_hit_index == -1) {
      THcGoodLADHit *goodhit = new ((*fGoodLADHits)[goodhit_n]) THcGoodLADHit();
      goodhit_n++;
      goodhit->CopyHit(plane_loc, 0, hit); // Copy the new hit into the good hit
      goodhit->SetTrackID(plane_loc, track_id);
    }
    if (matching_hit_index != -1) {
      THcGoodLADHit *goodhit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(matching_hit_index));
      if (abs(hit->GetdTrkHorizHit0()) <
          (plane_loc ? abs(goodhit->GetdTrkHorizHit1()) : abs(goodhit->GetdTrkHorizHit0()))) {

        goodhit->CopyHit(plane_loc, 0, hit); // Copy the new hit into the good hit
        goodhit->SetTrackID(plane_loc, track_id);
      }
    }
    // Todo: Calculate beta, alpha, etc. for the hit
    // FIXME. One track can still have multiple (distinct) hodo hits that are counted as good
  }

  //////////////////////////////////////////////////////////////////////////////
  // Find matching front + back plane hits
  // Create a vector of pairs to store matching hits for each track
  std::vector<std::pair<Int_t, Int_t>> matchingHits(ntracks, {-1, -1});
  std::vector<double> dTrk_horiz_values;

  // TODO: This can probably be incorporated into the loop above
  for (Int_t i = 0; i < goodhit_n; i++) {
    THcGoodLADHit *goodhit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(i));
    if (goodhit == nullptr)
      continue;

    // Add dTrk_horiz values for both planes (if valid) to the vector
    if (goodhit->GetdTrkHorizHit0() != 0) {
      dTrk_horiz_values.push_back(goodhit->GetdTrkHorizHit0());
    }
    if (goodhit->GetdTrkHorizHit1() != 0) {
      dTrk_horiz_values.push_back(goodhit->GetdTrkHorizHit1());
    }

    Int_t trackID0 = goodhit->GetTrackIDHit0();
    Int_t trackID1 = goodhit->GetTrackIDHit1();

    if (trackID0 >= 0 && trackID0 < ntracks) {
      // Check if the track ID is already in the vector
      if (matchingHits[trackID0].first == -1) {
        // If it is, set the second hit index to the current index
        matchingHits[trackID0].first = i;
      } else if (dTrk_horiz_values[matchingHits[trackID0].first] > goodhit->GetdTrkHorizHit0()) {
        matchingHits[trackID0].first = i;
      }
    }
    if (trackID1 >= 0 && trackID1 < ntracks) {
      // Check if the track ID is already in the vector
      if (matchingHits[trackID1].second == -1) {
        // If it is, set the second hit index to the current index
        matchingHits[trackID1].second = i;
      } else if (dTrk_horiz_values[matchingHits[trackID1].second] > goodhit->GetdTrkHorizHit1()) {
        matchingHits[trackID1].second = i;
      }
    }
  }

  // Loop over the matching hits and set the partner hit index
  for (auto &match : matchingHits) {
    if (match.first != -1 && match.second != -1) {
      THcGoodLADHit *firstHit  = static_cast<THcGoodLADHit *>(fGoodLADHits->At(match.first));
      THcGoodLADHit *secondHit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(match.second));
      if (firstHit && secondHit) {
        firstHit->CopyHit(1, 0, secondHit); // Copy the back (1) plane of the second hit into the first hit (0)

        // Remove the second hit from the good hits array
        fGoodLADHits->RemoveAt(match.second);
        goodhit_n--;
      }
    }
  }
  fGoodLADHits->Compress(); // Compress the array to remove null entries

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

  int fRFTimeIndex = (aparatus_prefix[0] == 'p') ? 0 : 1;
  fRFTime          = fTrigDet->Get_RF_TrigTime(fRFTimeIndex) + rf_offset[fRFTimeIndex];
  double tmp       = fmod(fTVertex - fRFTime, 4);
  if (tmp > 2)
    tmp -= 4;
  fTVertex_RFcorr = fTVertex - tmp;

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
