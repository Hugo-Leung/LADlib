#include "THcLADHodoscope.h"
#include "THaCutList.h"
#include "THaEvData.h"
#include "THaGlobals.h"
#include "THaSpectrometer.h"
#include "THaTrack.h"
#include "THcDetectorMap.h"
#include "THcGlobals.h"
#include "THcHitList.h"
#include "THcLADSpectrometer.h"
#include "THcParmList.h"
#include "THcSignalHit.h"
#include "VarDef.h"
#include "VarType.h"

//_________________________________________________________________
THcLADHodoscope::THcLADHodoscope(const char *name, const char *description, THaApparatus *apparatus)
    : THaNonTrackingDetector(name, description, apparatus) {
  // Constructor
  fNPlanes     = 1;
  fGoodLADHits = new TClonesArray("THcGoodLADHit", MAXGOODHITS);
}

//_________________________________________________________________
THcLADHodoscope::~THcLADHodoscope() {
  // Destructor

  for (int ip = 0; ip < fNPlanes; ip++)
    delete fPlanes[ip];
  delete[] fPlanes;

  delete[] fNPaddle;
  fNPaddle = NULL;
  delete[] fTdcOffset;
  fTdcOffset = NULL;
  delete[] fAdcTdcOffset;
  fAdcTdcOffset = NULL;
  delete[] fHodoSlop;
  fHodoSlop = NULL;

  delete[] fHodoBtmAdcTimeWindowMin;
  fHodoBtmAdcTimeWindowMin = NULL;
  delete[] fHodoBtmAdcTimeWindowMax;
  fHodoBtmAdcTimeWindowMax = NULL;
  delete[] fHodoTopAdcTimeWindowMin;
  fHodoTopAdcTimeWindowMin = NULL;
  delete[] fHodoTopAdcTimeWindowMax;
  fHodoTopAdcTimeWindowMax = NULL;

  delete[] fPlaneNames;
  fPlaneNames = NULL;

  delete[] fHodoVelLight;
  fHodoVelLight = NULL;

  // Time walk
  delete[] fHodoVelFit;
  fHodoVelFit = NULL;
  delete[] fHodoCableFit;
  fHodoCableFit = NULL;
  delete[] fHodo_LCoeff;
  fHodo_LCoeff = NULL;
  delete[] fHodoVelFit_FADC;
  fHodoVelFit_FADC = NULL;
  delete[] fHodoCableFit_FADC;
  fHodoCableFit_FADC = NULL;
  delete[] fHodo_LCoeff_FADC;
  fHodo_LCoeff_FADC = NULL;
  delete[] fHodoTop_c1;
  fHodoTop_c1 = NULL;
  delete[] fHodoBtm_c1;
  fHodoBtm_c1 = NULL;
  delete[] fHodoTop_c2;
  fHodoTop_c2 = NULL;
  delete[] fHodoBtm_c2;
  fHodoBtm_c2 = NULL;

  fGoodLADHits->Delete();
  delete fGoodLADHits;
}

//_________________________________________________________________
void THcLADHodoscope::Clear(Option_t *opt) {

  /*! \brief Clears variables
   *
   *  Called by  THcHodoscope::Decode
   *
   */

  if (*opt != 'I') {
    for (Int_t ip = 0; ip < fNPlanes; ip++) {
      fPlanes[ip]->Clear();
    }
  }

  // Hodo good hit variables
  goodhit_n            = 0;
  num_unique_good_hits = 0;
  num_unique_hits      = 0;
}
//_________________________________________________________________

void THcLADHodoscope::Setup(const char *name, const char *description) {
  /**
     Create the scintillator plane objects for the hodoscope.

     Uses the Xhodo_num_planes and Xhodo_plane_names to get the number of
     planes and their names.

  */
  if (IsZombie())
    return;

  // fDebug = 1;  // Keep this at one while we're working on the code

  char prefix[2];

  prefix[0] = tolower(GetApparatus()->GetName()[0]);
  prefix[1] = '\0';
  // Fix to prevent different param files for SHMS & HMS LAD hodoscopes
  char prefix_lad[2];
  prefix_lad[0] = 'l';
  prefix_lad[1] = '\0';

  string planenamelist;
  DBRequest listextra[] = {{"ladhodo_num_planes", &fNPlanes, kInt},
                           {"ladhodo_plane_names", &planenamelist, kString},
                           {"ladhodo_tdcrefcut", &fTDC_RefTimeCut, kInt, 0, 1},
                           {"ladhodo_adcrefcut", &fADC_RefTimeCut, kInt, 0, 1},
                           {0}};
  // fNPlanes = 4; 		// Default if not defined
  fTDC_RefTimeCut = 0; // Minimum allowed reference times
  fADC_RefTimeCut = 0;
  gHcParms->LoadParmValues((DBRequest *)&listextra, prefix_lad);

  cout << "Plane Name List : " << planenamelist << endl;

  vector<string> plane_names = Podd::vsplit(planenamelist);
  // Plane names
  if (plane_names.size() != (UInt_t)fNPlanes) {
    cout << "ERROR: Number of planes " << fNPlanes << " doesn't agree with number of plane names " << plane_names.size()
         << endl;
    // Should quit.  Is there an official way to quit?
  }

  fPlaneNames = new char *[fNPlanes];
  for (Int_t i = 0; i < fNPlanes; i++) {
    fPlaneNames[i] = new char[plane_names[i].length() + 1];
    strcpy(fPlaneNames[i], plane_names[i].c_str());
  }

  // Probably shouldn't assume that description is defined
  char *desc = new char[strlen(description) + 100];
  fPlanes    = new THcLADHodoPlane *[fNPlanes];
  for (Int_t i = 0; i < fNPlanes; i++) {
    strcpy(desc, description);
    strcat(desc, " Plane ");
    strcat(desc, fPlaneNames[i]);
    fPlanes[i] = new THcLADHodoPlane(fPlaneNames[i], desc, i + 1, this); // Number planes starting from zero!!
    cout << "Created Scintillator Plane " << fPlaneNames[i] << ", " << desc << endl;
  }

  delete[] desc;
}

//_________________________________________________________________
THaAnalysisObject::EStatus THcLADHodoscope::Init(const TDatime &date) {
  Setup(GetName(), GetTitle());

  char EngineDID[] = "xLADSCIN"; // LADSCIN avoids confusion with spectrometer scintillators
  EngineDID[0]     = toupper(GetApparatus()->GetName()[0]);
  if (gHcDetectorMap->FillMap(fDetMap, EngineDID) < 0) {
    static const char *const here = "Init()";
    Error(Here(here), "Error filling detectormap for %s.", EngineDID);
    return kInitError;
  }

  // Should probably put this in ReadDatabase as we will know the
  // maximum number of hits after setting up the detector map
  // But it needs to happen before the sub detectors are initialized
  // so that they can get the pointer to the hitlist.
  cout << " Hodo tdc ref time cut = " << fTDC_RefTimeCut << " " << fADC_RefTimeCut << endl;

  InitHitList(fDetMap, "THcRawHodoHit", fDetMap->GetTotNumChan() + 1, fTDC_RefTimeCut, fADC_RefTimeCut);

  EStatus status;
  if ((status = THaNonTrackingDetector::Init(date)))
    return fStatus = status;

  for (Int_t ip = 0; ip < fNPlanes; ip++) {
    if ((status = fPlanes[ip]->Init(date))) {
      return fStatus = status;
    }
  }

  fPresentP        = 0;
  THaVar *vpresent = gHaVars->Find(Form("%s.present", GetApparatus()->GetName()));
  if (vpresent) {
    fPresentP = (Bool_t *)vpresent->GetValuePointer();
  }

  return kOK;
}

//_________________________________________________________________
Int_t THcLADHodoscope::End(THaRunBase *run) {
  // Do we really need this function?

  return 0;
}

//_________________________________________________________________
Int_t THcLADHodoscope::DefineVariables(EMode mode) {
  // Initialize variables to read out
  RVarDef vars[] = {
      {"goodhit_n", "Number of good hits", "goodhit_n"},
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
      {"goodhit_hitedep_1", "Good hit energy deposition (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitEdepHit1()"},
      {"good_hit_n_unique", "Number of unique good hits", "num_unique_good_hits"},
      {"all_hits_n_unique", "Number of all hits, not just with tracks", "num_unique_hits"},
      {0}};

  return DefineVarsFromList(vars, mode);
}

//_________________________________________________________________
Int_t THcLADHodoscope::ReadDatabase(const TDatime &date) {

  cout << "THcLADHodoscope::ReadDatabase()" << endl;
  char prefix[2];
  prefix[0] = tolower(GetApparatus()->GetName()[0]); // "lad"
  prefix[1] = '\0';
  // Fix to prevent different param files for SHMS & HMS LAD hodoscopes
  char prefix_lad[2];
  prefix_lad[0] = 'l';
  prefix_lad[1] = '\0';

  // since we define each hodoscope as a separate detector
  // we will need to use the detector name to load parameters
  // for each detector -- to be updated
  fNPaddle = new Int_t[fNPlanes];
  for (int ip = 0; ip < fNPlanes; ip++) {
    string parname    = "ladhodo_" + string(fPlanes[ip]->GetName()) + "_nr";
    DBRequest list2[] = {{parname.c_str(), &fNPaddle[ip], kInt}, {0}};

    gHcParms->LoadParmValues((DBRequest *)&list2, prefix_lad);
  }

  // for all planes
  fTdcOffset    = new Int_t[fNPlanes];
  fAdcTdcOffset = new Double_t[fNPlanes];
  fHodoSlop     = new Double_t[fNPlanes];
  for (int ip = 0; ip < fNPlanes; ip++) {
    fTdcOffset[ip]    = 0;
    fAdcTdcOffset[ip] = 0;
    fHodoSlop[ip]     = 0;
  }

  // for all elements
  fMaxHodoScin             = fNPaddle[0] * fNPlanes;
  fHodoTopAdcTimeWindowMin = new Double_t[fMaxHodoScin];
  fHodoTopAdcTimeWindowMax = new Double_t[fMaxHodoScin];
  fHodoBtmAdcTimeWindowMin = new Double_t[fMaxHodoScin];
  fHodoBtmAdcTimeWindowMax = new Double_t[fMaxHodoScin];

  fHodoVelLight = new Double_t[fMaxHodoScin];

  // Time walk
  fHodoVelFit        = new Double_t[fMaxHodoScin];
  fHodoCableFit      = new Double_t[fMaxHodoScin];
  fHodo_LCoeff       = new Double_t[fMaxHodoScin];
  fHodoVelFit_FADC   = new Double_t[fMaxHodoScin];
  fHodoCableFit_FADC = new Double_t[fMaxHodoScin];
  fHodo_LCoeff_FADC  = new Double_t[fMaxHodoScin];
  fHodoTop_c1        = new Double_t[fMaxHodoScin];
  fHodoBtm_c1        = new Double_t[fMaxHodoScin];
  fHodoTop_c2        = new Double_t[fMaxHodoScin];
  fHodoBtm_c2        = new Double_t[fMaxHodoScin];
  fEdep2MeV_int      = new Double_t[fMaxHodoScin];
  fEdep2MeV_amp      = new Double_t[fMaxHodoScin];

  for (int ii = 0; ii < fMaxHodoScin; ii++) {
    fHodoTopAdcTimeWindowMin[ii] = -1000.;
    fHodoTopAdcTimeWindowMax[ii] = 1000.;
    fHodoBtmAdcTimeWindowMin[ii] = -1000.;
    fHodoBtmAdcTimeWindowMax[ii] = 1000.;
    fEdep2MeV_int[ii]            = 1.0; // Default value
    fEdep2MeV_amp[ii]            = 1.0; // Default value
  }

  DBRequest list3[]  = {{"ladcosmicflag", &fCosmicFlag, kInt, 0, 1},
                        {"ladhodo_tdc_min", &fScinTdcMin, kDouble},
                        {"ladhodo_tdc_max", &fScinTdcMax, kDouble},
                        {"ladhodo_tdc_to_time", &fScinTdcToTime, kDouble},
                        {"ladhodo_TopAdcTimeWindowMin", fHodoTopAdcTimeWindowMin, kDouble, (UInt_t)fMaxHodoScin, 1},
                        {"ladhodo_TopAdcTimeWindowMax", fHodoTopAdcTimeWindowMax, kDouble, (UInt_t)fMaxHodoScin, 1},
                        {"ladhodo_BtmAdcTimeWindowMin", fHodoBtmAdcTimeWindowMin, kDouble, (UInt_t)fMaxHodoScin, 1},
                        {"ladhodo_BtmAdcTimeWindowMax", fHodoBtmAdcTimeWindowMax, kDouble, (UInt_t)fMaxHodoScin, 1},
                        {"ladhodo_matching_paddle_tol", &fMatchingPaddleTol, kInt, 0, 1},
                        {"ladhodo_matching_time_tol", &fMatchingTimeTol, kDouble, 0, 1},
                        {"ladhodo_matching_dy_tol", &fMatchingDyTol, kDouble, 0, 1},
                        {0}};
  fMatchingPaddleTol = 0;
  fMatchingTimeTol   = 10.0;
  fMatchingDyTol     = 20.0;
  fCosmicFlag        = 0;
  fScinTdcMin        = 0;
  fScinTdcMax        = 0;
  fScinTdcToTime     = 0;

  gHcParms->LoadParmValues((DBRequest *)&list3, prefix_lad);

  DBRequest list6[] = {{"ladhodo_adc_tdc_offset", fAdcTdcOffset, kDouble, (UInt_t)fNPlanes, 1},
                       {"ladhodo_tdc_offset", fTdcOffset, kInt, (UInt_t)fNPlanes, 1},
                       {0}};
  gHcParms->LoadParmValues((DBRequest *)&list6, prefix);

  DBRequest list5[] = {{"is_mc", &fIsMC, kInt, 0, 1}, {0}};
  fIsMC             = 0;
  gHcParms->LoadParmValues((DBRequest *)&list5, "");

  DBRequest list[] = {{"ladhodo_vel_light", &fHodoVelLight[0], kDouble, (UInt_t)fMaxHodoScin, 1}, {0}};
  gHcParms->LoadParmValues((DBRequest *)&list, prefix_lad);

  DBRequest list4[] = {{"ladhodo_velFit", &fHodoVelFit[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_cableFit", &fHodoCableFit[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_LCoeff", &fHodo_LCoeff[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_velFit_FADC", &fHodoVelFit_FADC[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_cableFit_FADC", &fHodoCableFit_FADC[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_LCoeff_FADC", &fHodo_LCoeff_FADC[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_c1_Top", &fHodoTop_c1[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_c1_Btm", &fHodoBtm_c1[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_c2_Top", &fHodoTop_c2[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_c2_Btm", &fHodoBtm_c2[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_adcAmp2MeV", &fEdep2MeV_int[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_adcInt2MeV", &fEdep2MeV_amp[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"TDC_threshold", &fTdc_Thrs, kDouble, 0, 1},
                       {0}};

  fTdc_Thrs = 1.0;
  // Set Default Values if NOT defined in param file
  for (int i = 0; i < fMaxHodoScin; i++) {

    // Turn OFF Time-Walk Correction if param file NOT found
    fHodoTop_c1[i] = 0.0;
    fHodoTop_c2[i] = 0.0;
    fHodoBtm_c1[i] = 0.0;
    fHodoBtm_c2[i] = 0.0;
  }
  for (int i = 0; i < fMaxHodoScin; i++) {
    // Set scin Velocity/Cable to default
    fHodoCableFit[i]    = 0.0;
    fHodoVelFit[i]      = 15.0;
    fHodoVelFit_FADC[i] = 15.0;
    // set time coeff between paddles to default
    fHodo_LCoeff[i] = 0.0;
  }

  gHcParms->LoadParmValues((DBRequest *)&list4, prefix_lad);

  return kOK;
}

//_________________________________________________________________
Int_t THcLADHodoscope::Decode(const THaEvData &evdata) {

  // Decode raw data and pass it to hitlist
  // Read raw data -- THcHitList::DecodeToHitList
  // Processing hitlist for each plane -- THcLADHodoPlane::ProcessHits

  // If one spectrometer triggers, don't process anything related to LAD in the other spectrometer. Regular spectrometer
  // values will still be processed.
  // Get the spectrometer prefix
  TString prefix = GetApparatus()->GetName();
  prefix.ToLower();
  if ((gHaCuts->Result("SHMS_event") && prefix == "h") || (gHaCuts->Result("HMS_event") && prefix == "p")) {
    return 0;
  }

  Bool_t present = kTRUE;
  if (fPresentP) {
    present = *fPresentP;
  }

  if (fIsMC) {
    fNSA          = 1;
    fNSB          = 0;
    fNPED         = 1;
    fHaveFADCInfo = true;
  }

  fNHits = DecodeToHitList(evdata, !present);

  // To analyze pedestal events -- Must define "Pedestal_event" cut in the cuts .def file
  // do we want to do this? or calculate pedestal for each event (using the first # of samples, e.g)
  // keeping it for now
  if (gHaCuts->Result("Pedestal_event")) {
    Int_t nexthit = 0;
    for (Int_t ip = 0; ip < fNPlanes; ip++) {
      nexthit = fPlanes[ip]->AccumulatePedestals(fRawHitList, nexthit);
    }
    fAnalyzePedestals = 1; // Analyze pedestals first normal events
    return (0);
  }
  if (fAnalyzePedestals) {
    for (Int_t ip = 0; ip < fNPlanes; ip++) {
      fPlanes[ip]->CalculatePedestals();
    }
    fAnalyzePedestals = 0; // Don't analyze pedestals next event
  }

  Int_t nexthit = 0;
  for (Int_t iplane = 0; iplane < fNPlanes; iplane++) {

    nexthit = fPlanes[iplane]->ProcessHits(fRawHitList, nexthit);
  }

  return fNHits;
}

//_________________________________________________________________
Int_t THcLADHodoscope::CoarseProcess(TClonesArray &tracks) { return 0; }

//_________________________________________________________________
Int_t THcLADHodoscope::FineProcess(TClonesArray &tracks) {

  fGoodLADHits->Delete();
  // Hodo tracking coes in FineProcess, to ensure it comes after GEM Coarse Process tracking, no matter what.

  fSpectro = dynamic_cast<THaApparatus *>(GetApparatus());
  fGEM     = dynamic_cast<THcLADGEM *>(fSpectro->GetDetector("gem"));
  if (!fGEM) {
    return 0;
  }
  TClonesArray *fGEMTracks = fGEM->GetTracks();

  if (!fGEMTracks) {
    return 0;
  }

  // TClonesArray* gemTracks = spectrometer->GetGEM()->GetTracks();

  Int_t ntracks = fGEMTracks->GetLast() + 1;

  vector<Double_t> nPmtHit(ntracks);
  vector<Double_t> timeAtFP(ntracks);
  // Loop over all tracks and get corrected time, tof, beta...
  num_unique_hits = 0;
  for (Int_t ip = 0; ip < fNPlanes; ip++) {
    if (strcmp(fPlanes[ip]->GetName(), "REFBAR") == 0) {
      continue;
    }
    num_unique_hits += fPlanes[ip]->GetNScinHits();
  }
  num_unique_good_hits = 0;
  std::map<std::pair<int, int>, bool> hitUsedMap;

  // Initialize the map for all hits
  for (Int_t ip = 0; ip < fNPlanes; ip++) {
    for (Int_t iphit = 0; iphit < fPlanes[ip]->GetNScinHits(); iphit++) {
      THcLADHodoHit *hit                     = (THcLADHodoHit *)fPlanes[ip]->GetHits()->At(iphit);
      int paddle                             = hit->GetPaddleNumber() - 1;
      hitUsedMap[std::make_pair(ip, paddle)] = false;
    }
  }

  // Pair front and back hits to make good hits
  for (Int_t ip = 0; ip < fNPlanes; ip += 2) {
    for (Int_t iphit = 0; iphit < fPlanes[ip]->GetNScinHits(); iphit++) {
      THcLADHodoHit *hit = (THcLADHodoHit *)fPlanes[ip]->GetHits()->At(iphit);
      if (hitUsedMap[std::make_pair(ip, hit->GetPaddleNumber() - 1)]) {
        continue;
      }
      if (ip + 1 < fNPlanes) {
        for (Int_t iphit2 = 0; iphit2 < fPlanes[ip + 1]->GetNScinHits(); iphit2++) {
          THcLADHodoHit *hit2 = (THcLADHodoHit *)fPlanes[ip + 1]->GetHits()->At(iphit2);
          if (hitUsedMap[std::make_pair(ip + 1, hit2->GetPaddleNumber() - 1)]) {
            continue;
          }

          bool same_paddle = fabs(hit2->GetPaddleNumber() - hit->GetPaddleNumber()) <= fMatchingPaddleTol;
          bool time_match  = fabs(hit2->GetScinCorrectedTime() - hit->GetScinCorrectedTime()) < fMatchingTimeTol;
          bool dy_match    = fabs(hit2->GetCalcPosition() - hit->GetCalcPosition()) < fMatchingDyTol;
          if (same_paddle && time_match && dy_match) {
            hitUsedMap[std::make_pair(ip, hit->GetPaddleNumber() - 1)]      = true;
            hitUsedMap[std::make_pair(ip + 1, hit2->GetPaddleNumber() - 1)] = true;
            THcGoodLADHit *goodhit = new ((*fGoodLADHits)[goodhit_n]) THcGoodLADHit();
            goodhit_n++;
            goodhit->SetPlane(0, ip);
            goodhit->SetPaddle(0, hit->GetPaddleNumber() - 1);
            // goodhit->SetTrackID(0, -1); // Set to -1 for now, will be updated in THcLADKine
            // goodhit->SetdTrkHoriz(0,
            //                       hit->GetCalcPosition()); // Set to hit position for now, will be updated in
            //                       THcLADKine
            // goodhit->SetdTrkVert(0,
            //                      hit->GetCalcPosition()); // Set to hit position for now, will be updated in
            //                      THcLADKine
            // goodhit->SetHitTheta(0, 0);                   // Set to 0 for now, will be updated in THcLADKine
            // goodhit->SetHitPhi(0, 0);                     // Set to 0 for now, will be updated in THcLADKine
            goodhit->SetHitTime(0, hit->GetScinCorrectedTime());
            goodhit->SetHitEdep(0, hit->GetPaddleADC());
            goodhit->SetHitEdepAmp(0, hit->GetPaddleADCpeak());
            goodhit->SetHitYPos(0, hit->GetCalcPosition());
            goodhit->SetPlane(1, ip + 1);
            goodhit->SetPaddle(1, hit2->GetPaddleNumber() - 1);
            // goodhit->SetTrackID(1, -1); // Set to -1 for now, will be updated in THcLADKine
            // goodhit->SetdTrkHoriz(1, hit2->GetCalcPosition()); // Set to hit position for now, will be updated in
            // THcLADKine goodhit->SetdTrkVert(1, hit2->GetCalcPosition());  // Set to hit position for now, will be
            // updated in THcLADKine goodhit->SetHitTheta(1, 0); // Set to 0 for now, will be updated in THcLADKine
            // goodhit->SetHitPhi(1, 0);   // Set to 0 for now, will be updated in THcLADKine
            goodhit->SetHitTime(1, hit2->GetScinCorrectedTime());
            goodhit->SetHitEdep(1, hit2->GetPaddleADC());
            goodhit->SetHitEdepAmp(1, hit2->GetPaddleADCpeak());
            goodhit->SetHitYPos(1, hit2->GetCalcPosition());
          }
        }
      }
    }
  }

  for (Int_t ip = 0; ip < fNPlanes; ip++) {
    // skip reference bar plane
    if (strcmp(fPlanes[ip]->GetName(), "REFBAR") == 0)
      continue;
    // only include planes 0,2,4
    if (ip != 0 && ip != 2 && ip != 4)
      continue;

    for (Int_t iphit = 0; iphit < fPlanes[ip]->GetNScinHits(); iphit++) {
      THcLADHodoHit *hit = (THcLADHodoHit *)fPlanes[ip]->GetHits()->At(iphit);
      int paddle         = hit->GetPaddleNumber() - 1;
      auto key           = std::make_pair(ip, paddle);
      if (hitUsedMap[key])
        continue;

      // mark as used and create a single-plane good hit
      hitUsedMap[key]        = true;
      THcGoodLADHit *goodhit = new ((*fGoodLADHits)[goodhit_n]) THcGoodLADHit();
      goodhit_n++;
      goodhit->SetPlane(0, ip);
      goodhit->SetPaddle(0, paddle);
      goodhit->SetHitTime(0, hit->GetScinCorrectedTime());
      goodhit->SetHitEdep(0, hit->GetPaddleADC());
      goodhit->SetHitEdepAmp(0, hit->GetPaddleADCpeak());
      goodhit->SetHitYPos(0, hit->GetCalcPosition());
    }
  }

  num_unique_good_hits = 0;
  for (const auto &hit : hitUsedMap) {
    if (hit.second) {
      num_unique_good_hits++;
    }
  }

  return 0;
}

TVector3 THcLADHodoscope::GetHitPositionLab(int plane, int paddle, double ypos) {
  // Get the position of a hit in the hodoscope, given the plane and paddle number
  // This is a placeholder function and should be implemented based on the actual geometry of the hodoscope
  if (plane < 0 || plane >= fNPlanes) {
    cout << "[THcLADHodoscope] Error: Invalid plane number" << endl;
    return TVector3(0, 0, 0);
  }
  if (paddle < 0 || paddle >= fNPaddle[plane]) {
    cout << "[THcLADHodoscope] Error: Invalid paddle number" << endl;
    return TVector3(0, 0, 0);
  }
  Double_t offset  = fPlanes[plane]->GetPosCenter(paddle);
  TVector3 hit = TVector3(offset, ypos , fPlanes[plane]->GetZpos());
  hit.RotateY(fPlanes[plane]->GetTheta());
  return hit;
}


ClassImp(THcLADHodoscope)
