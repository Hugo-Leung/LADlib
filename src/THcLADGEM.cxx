#include "THcLADGEM.h"
#include "THaApparatus.h"
#include "THaCutList.h"
#include "THaEvData.h"
#include "THaGlobals.h"
#include "THaTrack.h"
#include "THcGlobals.h"
#include "THcHallCSpectrometer.h"
#include "THcParmList.h"
#include "VarDef.h"
#include "VarType.h"
#include "THaRunBase.h"

using namespace std;

//____________________________________________________________
THcLADGEM::THcLADGEM(const char *name, const char *description, THaApparatus *apparatus)
    : THaNonTrackingDetector(name, description, apparatus) {
  // constructor

  fModules.clear();
  fModulesInitialized = false;

  fNModules = 0;
  fNLayers  = 0;
  fNhits    = 0;

  fGEMTracks    = new TClonesArray("THcLADGEMTrack", MAXTRACKS);
  fVertexModule = nullptr;
  fIgnoreVertex = false;
}

//____________________________________________________________
THcLADGEM::~THcLADGEM() {
  // Destructor
  for (auto module : fModules) {
    delete module;
  }
  fModules.clear();

  delete fGEMTracks;
}

//____________________________________________________________
void THcLADGEM::Clear(Option_t *opt) {

  //  cout << "THcLADGEM::Clear" << endl;
  fClusOutData.clear();
  fNClusters = 0;
  fNTracks   = 0;

  for (auto &module : fModules)
    module->Clear();

  for (int i = 0; i < fNLayers; i++)
    f2DHits[i].clear();

  fNhits = 0;
  fPosX.clear();
  fPosY.clear();
  fPosZ.clear();
  fTimeMean.clear();
  fADCMean.clear();
  fADCAsym.clear();
  fTimeDiff.clear();
  fTimeCorr.clear();
  fIsGoodHit.clear();
  fClusID0.clear();
  fClusID1.clear();
  fSPID.clear();
  fLayer.clear();
}

//____________________________________________________________
THaAnalysisObject::EStatus THcLADGEM::Init(const TDatime &date) {

  //  cout << "THcLADGEM::Init" << endl;

  EStatus status;
  if ((status = THaNonTrackingDetector::Init(date)))
    return fStatus = status;

  fPresentP        = nullptr;
  THaVar *vpresent = gHaVars->Find(Form("%s.present", GetApparatus()->GetName()));
  if (vpresent) {
    fPresentP = (Bool_t *)vpresent->GetValuePointer();
  }

  // Call Hall C style DetectorMap
  // GEM channel maps to be defined in MAPS/LAD/detector.map
  // e.g. gHcDetectorMap->FillMap(fDetMap, EngineID)
  // InitHitList(fDetMap, "THcLADGEMHit", fDetMap->GetTotNumChan()+1,0, RefTimeCut)

  // Init Subdetectors
  for (auto &module : fModules) {
    status = module->Init(date);
    if (status != kOK)
      return fStatus = status;
  }

  // Load pedestal and CM if files are provided
  if (!fPedFilename.empty())
    LoadPedestals();
  if (!fCMFilename.empty())
    LoadCM();

  return fStatus = kOK;
}

//____________________________________________________________________
Int_t THcLADGEM::Begin( THaRunBase* run ) {
  for( auto& module: fModules ) {
    module->Begin(run);
  }
  return 0;
}

//____________________________________________________________
Int_t THcLADGEM::DefineVariables(EMode mode) {
  if (mode == kDefine && fIsSetup)
    return kOK;
  fIsSetup = (mode == kDefine);

  // Cluster variables
  RVarDef vars_clus[] = {{"clust.layer", "GEM Layer", "fClusOutData.layer"},
                         {"clust.module", "GEM Module ID", "fClusOutData.imodule"},
                         {"clust.axis", "U/V axis", "fClusOutData.axis"},
                         {"clust.mpd", "MPD ID", "fClusOutData.mpdid"},
                         {"clust.nstrip", "Number of strips in cluster", "fClusOutData.nstrip"},
                         {"clust.maxstrip", "Max strip of the given cluster", "fClusOutData.maxstrip"},
                         {"clust.index", "Cluster index", "fClusOutData.clindex"},
                         {"clust.adc", "Cluster ADC sum", "fClusOutData.adc"},
                         {"clust.time", "Cluster Time mean", "fClusOutData.time"},
                         {"clust.pos", "Weighted cluster position", "fClusOutData.pos"},
                         {"clust.maxpos", "Max strip pos of the cluster", "fClusOutData.mpos"},
                         {"clust.maxsamp", "Time sample with max ADC", "fClusOutData.maxsamp"},
                         {"clust.maxadc", "Max strip ADC", "fClusOutData.maxadc"},
                         {0}};

  DefineVarsFromList(vars_clus, mode);

  // Hit variables
  RVarDef vars_hit[] = {{"hit.nlayer", "number of layers with any strip fired", "fNlayers_hit"},
                        {"hit.nlayeru", "number of layers with U strip fired", "fNlayers_hitU"},
                        {"hit.nlayerv", "number of layers with V strip fired", "fNlayers_hitV"},
                        {"hit.nlayeruv", "number of layers with 2d hits", "fNlayers_hitUV"},
                        {"hit.nstripsu_layer", "total number of U strips fired by layer", "fNstripsU_layer"},
                        {"hit.nstripsv_layer", "total number of V strips fired by layer", "fNstripsV_layer"},
                        {"hit.nclustu_layer", "total number of U clusters by layer", "fNclustU_layer"},
                        {"hit.nclustv_layer", "total number of V clusters by layer", "fNclustV_layer"},
                        {0}};
  DefineVarsFromList(vars_hit, mode);

  RVarDef vars_sp[] = {{"sp.nhits", "Number of hits in GEM layer 0", "nhits"},
                       {"sp.posX", "X position of GEM hit in layer 0", "fPosX"},
                       {"sp.posY", "Y position of GEM hit in layer 0", "fPosY"},
                       {"sp.posZ", "Z position of GEM hit in layer 0", "fPosZ"},
                       {"sp.time", "Time mean of GEM hit in layer 0", "fTimeMean"},
                       {"sp.adc", "ADC mean of GEM hit in layer 0", "fADCMean"},
                       {"sp.asym", "ADC asym of GEM hit in layer 0", "fADCAsym"},
                       {"sp.dt", "Time difference of GEM hit in layer 0", "fTimeDiff"},
                       {"sp.ct", "Corrected time of GEM hit in layer 0", "fTimeCorr"},
                       {"sp.isgoodhit", "Good hit flag of GEM hit in layer 0", "fIsGoodHit"},
                       {"sp.clusID1", "Cluster ID1 of GEM hit in layer 0", "fClusID0"},
                       {"sp.clusID2", "Cluster ID2 of GEM hit in layer 0", "fClusID1"},
                       {"sp.spID", "Track ID of GEM hit in layer 0", "fSPID"},
                       {"sp.layer", "2D hit layer", "fLayer"},
                       {0}};
  DefineVarsFromList(vars_sp, mode);

  //  Place holder for space point output
  /*
    RVarDef vars_sp[] = {
      {"sp.layer",      "2D hit layer",                             ""},
      {"sp.trackid",    "Track ID associated with the space point", ""},
      {"sp.clusid",     "Cluster IDs associated with the sp",       ""},
      {"sp.adc",        "2D hit ADC mean",                          ""},
      {"sp.asy",        "2D hit ADC asym",                          ""},
      {"sp.x",          "2D hit x pos",                             ""},
      {"sp.y",          "2D hit y pos",                             ""},
      {"sp.z",          "2D hit z pos",                             ""},
      {"sp.t",          "2D hit avg time",                          ""},
      {"sp.dt",         "2D hit time difference",                   ""},
      {"sp.ct",         "2D hit corrected time",                    ""},
      {"sp.isgoodhit",  "2D hit good hit flag",                     ""},
      { 0 }
    };
    DefineVarsFromList( vars_sp, mode );
  */

  // track variables only available when there are at least two layers
  if (fNLayers > 1) {
    // Track/Space point variables
    RVarDef vars_trk[] = {
        {"trk.ntracks", "Number of GEM track candidates", "fNTracks"},
        {"trk.id", "GEM Track ID", "fGEMTracks.THcLADGEMTrack.GetTrackID()"},
        {"trk.spID_0u", "Space Point ID for Layer 0 U", "fGEMTracks.THcLADGEMTrack.GetSpacePointID_0U()"},
        {"trk.spID_0v", "Space Point ID for Layer 0 V", "fGEMTracks.THcLADGEMTrack.GetSpacePointID_0V()"},
        {"trk.spID_1u", "Space Point ID for Layer 1 U", "fGEMTracks.THcLADGEMTrack.GetSpacePointID_1U()"},
        {"trk.spID_1v", "Space Point ID for Layer 1 V", "fGEMTracks.THcLADGEMTrack.GetSpacePointID_1V()"},
        {"trk.adc1", "2D hit ADC mean for 1st layer", "fGEMTracks.THcLADGEMTrack.GetADCMean_Sp1()"},
        {"trk.adc2", "2D hit ADC mean for 2nd layer", "fGEMTracks.THcLADGEMTrack.GetADCMean_Sp2()"},
        {"trk.asy1", "2D hit ADC asym for 1st layer", "fGEMTracks.THcLADGEMTrack.GetADCasym_Sp1()"},
        {"trk.asy2", "2D hit ADC asym for 2nd layer", "fGEMTracks.THcLADGEMTrack.GetADCasym_Sp2()"},
        {"trk.x1", "Space point1 X", "fGEMTracks.THcLADGEMTrack.GetX1()"},
        {"trk.y1", "Space point1 Y", "fGEMTracks.THcLADGEMTrack.GetY1()"},
        {"trk.z1", "Space point1 Z", "fGEMTracks.THcLADGEMTrack.GetZ1()"},
        {"trk.x2", "Space point2 X", "fGEMTracks.THcLADGEMTrack.GetX2()"},
        {"trk.y2", "Space point2 Y", "fGEMTracks.THcLADGEMTrack.GetY2()"},
        {"trk.z2", "Space point2 Z", "fGEMTracks.THcLADGEMTrack.GetZ2()"},
        {"trk.x1_local", "Space point1 X Local", "fGEMTracks.THcLADGEMTrack.GetX1_local()"},
        {"trk.y1_local", "Space point1 Y Local", "fGEMTracks.THcLADGEMTrack.GetY1_local()"},
        {"trk.x2_local", "Space point2 X Local", "fGEMTracks.THcLADGEMTrack.GetX2_local()"},
        {"trk.y2_local", "Space point2 Y Local", "fGEMTracks.THcLADGEMTrack.GetY2_local()"},
        {"trk.t", "Avg time", "fGEMTracks.THcLADGEMTrack.GetT()"},
        {"trk.dt", "Time difference between two sp", "fGEMTracks.THcLADGEMTrack.GetdT()"},
        {"trk.d0", "Track dist from vertex", "fGEMTracks.THcLADGEMTrack.GetD0()"},
        {"trk.d0_good", "Track dist from true vertex", "fGEMTracks.THcLADGEMTrack.GetGoodD0()"},
        {"trk.projz", "Projected z-vertex", "fGEMTracks.THcLADGEMTrack.GetProjVz()"},
        {"trk.projy", "Projected y-vertex", "fGEMTracks.THcLADGEMTrack.GetProjVy()"},
        {"trk.has_hodo_hit", "Track has hodoscope hit", "fGEMTracks.THcLADGEMTrack.GetHasHodoHit()"},
        // {"trk.theta", "Track theta", "fGEMTracks.THcLADGEMTrack.GetTheta()"},
        // {"trk.phi", "Track phi", "fGEMTracks.THcLADGEMTrack.GetPhi()"},
        {0}};
    DefineVarsFromList(vars_trk, mode);
  }

  return kOK;
}

//____________________________________________________________
Int_t THcLADGEM::ReadDatabase(const TDatime &date) {
  //  cout << "THcLADGEM::ReadDatabase" << endl;
  // Called by THaDetectorBase::Init()
  // Read parameters from THcParmList

  char prefix[2];

  prefix[0] = std::tolower(GetApparatus()->GetName()[0]);
  // Fix to prevent different param files for SHMS & HMS LAD hodoscopes
  prefix[0] = 'l';
  prefix[1] = '\0';

  // initial values
  fD0Cut       = 100.0; // DCA cut in cm
  fPedFilename = "";
  fCMFilename  = "";
  fPedestalMode = 0;
  DBRequest list[] = {{"gem_num_modules", &fNModules, kInt}, // should be defined in DB file
                      {"gem_num_layers", &fNLayers, kInt},
                      {"gem_pedfile", &fPedFilename, kString, 0, 1},
                      {"gem_cmfile", &fCMFilename, kString, 0, 1},
                      {"gem_d0_cut", &fD0Cut, kDouble, 0, 1},
                      {0}

  };
  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  DBRequest list2[] = {{"lgem_pedestal_mode", &fPedestalMode, kInt, 0, 1}};
  gHcParms->LoadParmValues((DBRequest *)&list2,"");
  std::cout<< "THcLADGEM: Pedestal mode set to " << fPedestalMode << std::endl;
  // Define GEM Modules
  for (int imod = 0; imod < fNModules; imod++) {
    THcLADGEMModule *new_module = new THcLADGEMModule(Form("m%d", imod), Form("m%d", imod), imod, this);
    fModules.push_back(new_module);
  }

  f2DHits.resize(fNLayers);
  fNstripsU_layer.resize(fNLayers);
  fNstripsV_layer.resize(fNLayers);
  fNclustU_layer.resize(fNLayers);
  fNclustV_layer.resize(fNLayers);
  return kOK;
}

//____________________________________________________________
Int_t THcLADGEM::Decode(const THaEvData &evdata) {
  // If one spectrometer triggers, don't process anything related to LAD in the other spectrometer. Regular spectrometer
  // values will still be processed.
  TString prefix = GetApparatus()->GetName();
  prefix.ToLower();
  if ((gHaCuts->Result("SHMS_event") && prefix == "h") || (gHaCuts->Result("HMS_event") && prefix == "p")) {
    return 0;
  }

  // Decode MPD data
  for (auto &module : fModules) {
    module->Decode(evdata);
  }

  return 0;
}

//____________________________________________________________
Int_t THcLADGEM::CoarseProcess(TClonesArray &tracks) {
  //  cout << "THcLADGEM::CoarseProcess" << endl;

  fNTracks = 0;
  fGEMTracks->Delete();
  // delete fGEMTracks;
  // fGEMTracks = new TClonesArray("THcLADGEMTrack", MAXTRACKS);

  fNlayers_hit   = 0;
  fNlayers_hitU  = 0;
  fNlayers_hitV  = 0;
  fNlayers_hitUV = 0;

  for (int i = 0; i < fNLayers; i++) {
    fNstripsU_layer[i] = 0;
    fNstripsV_layer[i] = 0;
    fNclustU_layer[i]  = 0;
    fNclustV_layer[i]  = 0;
  }

  for (auto module : fModules) {
    module->CoarseProcess(tracks); // X/Y clustering, form 2D hits

    // Counters, number of layers with strip fired (each layer has one module for LAD)
    fNstripsU_layer[module->GetLayerNum()] += module->GetNStripsHitU();
    fNstripsV_layer[module->GetLayerNum()] += module->GetNStripsHitV();

    if (module->GetNStripsHitU() > 0)
      fNlayers_hitU++;
    if (module->GetNStripsHitV() > 0)
      fNlayers_hitV++;
    if (module->GetNStripsHitU() + module->GetNStripsHitV() > 0)
      fNlayers_hit++;
    if (module->GetN2DHits() > 0)
      fNlayers_hitUV++;

    fNclustU_layer[module->GetLayerNum()] += module->GetNClusters(0);
    fNclustV_layer[module->GetLayerNum()] += module->GetNClusters(1);

    // Cluster output handling
    for (int i = 0; i < 2; i++) {
      // 0:U(X) cluster 1:V(Y) cluster
      for (auto &cluster : module->GetClusters(i)) {
        fClusOutData.layer.push_back(cluster.GetLayer());
        fClusOutData.imodule.push_back(module->GetModuleID());
        fClusOutData.mpdid.push_back(cluster.GetMPD());
        fClusOutData.axis.push_back(cluster.GetAxis());
        fClusOutData.nstrip.push_back(cluster.GetNStrips());
        fClusOutData.maxstrip.push_back(cluster.GetStripMax());
        fClusOutData.clindex.push_back(cluster.GetCLIndex());
        fClusOutData.adc.push_back(cluster.GetADCsum());
        fClusOutData.time.push_back(cluster.GetTime());
        fClusOutData.pos.push_back(cluster.GetPos());
        fClusOutData.mpos.push_back(cluster.GetPosMax());
        fClusOutData.maxsamp.push_back(cluster.GetSampMax());
        fClusOutData.maxadc.push_back(cluster.GetADCMax());
        fNClusters++;
      }
    }
  }

  // Loop over all 2D hits and find track candidates
  // Using only two layers to define a track candidate
  // LAD has only two layers...If more than two layers, use the outer two

  nhits = 0;
  for (int layer = 0; layer < fNLayers; ++layer) {
    nhits += f2DHits[layer].size();
  }
  // Loop through f2DHits and fill the vectors for all layers

  fPosX.reserve(nhits);
  fPosY.reserve(nhits);
  fPosZ.reserve(nhits);
  fTimeMean.reserve(nhits);
  fADCMean.reserve(nhits);
  fADCAsym.reserve(nhits);
  fTimeDiff.reserve(nhits);
  fTimeCorr.reserve(nhits);
  fIsGoodHit.reserve(nhits);
  fClusID0.reserve(nhits);
  fClusID1.reserve(nhits);
  fSPID.reserve(nhits);
  fLayer.reserve(nhits);

  for (int layer = 0; layer < fNLayers; ++layer) {
    for (const auto &hit : f2DHits[layer]) {
      fPosX.push_back(hit.posX);
      fPosY.push_back(hit.posY);
      fPosZ.push_back(hit.posZ);
      fTimeMean.push_back(hit.TimeMean);
      fADCMean.push_back(hit.ADCMean);
      fADCAsym.push_back(hit.ADCasym);
      fTimeDiff.push_back(hit.TimeDiff);
      fTimeCorr.push_back(hit.TimeCorr);
      fIsGoodHit.push_back(hit.IsGoodHit);
      fClusID0.push_back(hit.clusID[0]);
      fClusID1.push_back(hit.clusID[1]);
      fSPID.push_back(hit.spID);
      fLayer.push_back(layer);
    }
  }

  // if we have less than two layers, no tracking can be done
  if (fNLayers < 2)
    return 0;

  for (auto gemhit1 : f2DHits[fNLayers - 2]) {
    int gemhit1_id = 0;
    // LHE. The -X is a hard-coded fix to get the right coordinate system. This should be really easy to fix in the
    // param file instead, but we're currently trying to debug low-level gem issues, and doing this is one less moving
    // part.
    TVector3 v_hit1(gemhit1.posX, gemhit1.posY, gemhit1.posZ);
    for (auto gemhit2 : f2DHits[fNLayers - 1]) {
      int gemhit2_id = 0;

      if (!gemhit1.IsGoodHit || !gemhit2.IsGoodHit)
        continue;

      double tdiff = gemhit1.TimeMean - gemhit2.TimeMean;         // time difference (TimeMean1 - TimeMean2)
      double tmean = (gemhit1.TimeMean + gemhit2.TimeMean) * 0.5; // average time

      TVector3 v_hit2(gemhit2.posX, gemhit2.posY, gemhit2.posZ);

      // THaTrack* this_track = nullptr;
      // this_track = AddTrack(tracks, 0.0, 0.0, 0.0, 0.0); // AddTrack is func of THaTrackingDetector
      // FIXME: theta, phi might be defined differently in TVector3 and THaTrack
      // this_track->SetD(v_hit1.X(), v_hit1.Y(), v_hit1.Theta(), v_hit1.Phi() ); // DCS x, y , theta, phi

      // d0: DCAr from the primary vertex, assume (0,0,0) for now
      // we want to get the prima(0., 0., 0.);

      // LHE: Uncomment the lines below when runtime starts. Curently ok (doesn't cause crash, but throws many errors)
      TVector3 v_prim;
      TString fVertexModuleName = TString(GetApparatus()->GetName()) + ".react"; // Name is currently hard-coded to
      // be "react". Probably not worth changing

      if (fIgnoreVertex) {
        v_prim.SetXYZ(0., 0., 0.);
      } else {
        fVertexModule = dynamic_cast<THcReactionPoint *>(FindModule(fVertexModuleName.Data(), "THcReactionPoint"));

        if (fVertexModule && fVertexModule->HasVertex()) {
          v_prim = fVertexModule->GetVertex();
        } else {
          // Need to be carful that 0,0,0 doesn't get called during the run (or doesn't make it into the data)
          v_prim.SetXYZ(0., 0., 0.);
        }
      }

      double numer = ((v_prim - v_hit1).Cross((v_prim - v_hit2))).Mag();
      double denom = (v_hit2 - v_hit1).Mag();
      // here we can put a range/fiducial cut on d0 taking into account the target size
      double d0 = numer / denom;

      // Calculate d0 in the x-z plane (y=0)
      // double numer_xz = std::abs((v_prim.X() - v_hit1.X()) * (v_prim.Z() - v_hit2.Z()) -
      //            (v_prim.Z() - v_hit1.Z()) * (v_prim.X() - v_hit2.X()));
      // double denom_xz = std::sqrt(std::pow(v_hit2.X() - v_hit1.X(), 2) +
      //             std::pow(v_hit2.Z() - v_hit1.Z(), 2));
      // double d0 = numer_xz / denom_xz;

      if (d0 > fD0Cut) {
        // cout << "d0 too large: " << d0 << endl;
        continue;
      }
      // DCAz, projected z-vertex
      // First check if it intercepts with z-axis (within 1 cm)?
      double vpz;
      double vpy;
      double t1 = -v_hit1[0] / (v_hit2[0] - v_hit1[0]);
      double t2 = -v_hit1[1] / (v_hit2[1] - v_hit1[1]);
      if (abs(t1 - t2) > 200000) {
        vpz = -999999.;
        vpy = -999999.;
      } else {
        vpz = -v_hit1[0] * (v_hit2[2] - v_hit1[2]) / (v_hit2[0] - v_hit1[0]) + v_hit1[2];
        vpy = -v_hit1[0] * (v_hit2[1] - v_hit1[1]) / (v_hit2[0] - v_hit1[0]) + v_hit1[1];
      }

      gemhit1.spID = fNTracks;
      gemhit2.spID = fNTracks;

      if (fNTracks < MAXTRACKS) {
        // Add track object
        THcLADGEMTrack *theGEMTrack = new ((*fGEMTracks)[fNTracks]) THcLADGEMTrack(fNLayers);
        theGEMTrack->SetTrackID(fNTracks);
        theGEMTrack->AddSpacePoint(gemhit1);
        theGEMTrack->AddSpacePoint(gemhit2);
        theGEMTrack->SetTime(tmean, tdiff);
        theGEMTrack->SetD0(d0);
        theGEMTrack->SetZVertex(vpz);
        theGEMTrack->SetYVertex(vpy);
      }
      fNTracks++;
      gemhit2_id++;
    }
    gemhit1_id++;
  }

  return 0;
}
//____________________________________________________________
void THcLADGEM::RotateToLab(Double_t angle, TVector3 &vect) {

  // Initially defined in the local detector coordinate system (this is not transport coord)
  // Rotate along y-axis
  Double_t x_lab = vect.X() * cos(angle) + vect.Z() * sin(angle);
  Double_t y_lab = vect.Y();
  Double_t z_lab = -vect.X() * sin(angle) + vect.Z() * cos(angle);

  // Redefine the vector
  vect.SetXYZ(x_lab, y_lab, z_lab);
}

//____________________________________________________________
GEM2DHits* THcLADGEM::Add2DHits(Int_t ilayer, Double_t x, Double_t y, Double_t z, Double_t t, Double_t dt, Double_t tc,
                          Bool_t goodhit, Double_t adc, Double_t adcasy, Int_t clust_id1, Int_t clust_id2,
                          Int_t sp_index) {
  // FIXME:Add flag for filtering good hits?

  GEM2DHits gemhit;
  gemhit.Set(ilayer, x, y, z, t, dt, tc, goodhit, adc, adcasy);
  gemhit.SetClusterIDs(clust_id1, clust_id2);
  // gemhit.SetSpacePointIndex(sp_index);
  gemhit.SetSPID(sp_index);
  f2DHits[ilayer].push_back(gemhit);
  // Set initial spID = -1
  //  f2DHits[ilayer].push_back( {ilayer, x, y, z, t, dt, tc, goodhit, adc, adcasy, -1} );
  return &(f2DHits[ilayer].back());
}

//____________________________________________________________
Int_t THcLADGEM::FineProcess(TClonesArray &tracks) {
  //  cout << "THcLADGEM::FineProcess" << endl;

  // for (Int_t i = 0; i < std::min(fNTracks, MAXTRACKS); i++) {
  //   delete fGEMTracks->At(i);
  // }

  return 0;
}

//____________________________________________________________
void THcLADGEM::LoadPedestals() {

  std::ifstream pedfile(fPedFilename);
  if (!pedfile.good()) {
    std::cout << "THcLADGEM Warning: could not find ped file " << fPedFilename.c_str() << std::endl;
    return;
  } else {
    std::cout << "THcLADGEM: Read pedestal file " << fPedFilename.c_str() << std::endl;
  }

  // omg..
  std::map<int, std::map<int, std::map<int, std::vector<double>>>> PedMean;
  std::map<int, std::map<int, std::map<int, std::vector<double>>>> PedRMS;
  std::map<int, std::map<int, std::map<int, std::vector<int>>>> APVChan;

  int crate, slot, mpd, adc_ch;
  std::string thisline;
  while (std::getline(pedfile, thisline)) {
    if (pedfile.eof())
      break;
    if (thisline[0] != '#') {
      std::istringstream is(thisline);
      string dummy;

      if (thisline.find("APV") == 0) {
        is >> dummy >> crate >> slot >> mpd >> adc_ch;
      } else {
        int index = adc_ch + 16 * mpd;
        int apvchan;
        double mean, rms;

        is >> apvchan >> mean >> rms;

        PedMean[crate][slot][index].push_back(mean);
        PedRMS[crate][slot][index].push_back(rms);
        APVChan[crate][slot][index].push_back(apvchan);
      }
    }
  }

  for (int module = 0; module < fNModules; module++) {
    for (auto it = fModules[module]->fMPDmap.begin(); it != fModules[module]->fMPDmap.end(); ++it) {

      int this_crate = it->crate;
      int this_index = it->adc_id + 16 * it->mpd_id;
      int this_slot  = it->slot;

      if (PedMean.find(this_crate) != PedMean.end()) {
        if (PedMean[this_crate].find(this_slot) != PedMean[this_crate].end()) {
          if (PedMean[this_crate][this_slot].find(this_index) != PedMean[this_crate][this_slot].end()) {
            for (int i = 0; i < 128; i++) {
              int this_apvchan = APVChan[this_crate][this_slot][this_index][i];
              double this_mean = PedMean[this_crate][this_slot][this_index][i];
              double this_rms  = PedRMS[this_crate][this_slot][this_index][i];
              int this_strip   = fModules[module]->GetStripNumber(this_apvchan, it->axis, it->pos, it->invert);

              if (it->axis == LADGEM::kUaxis) {
                fModules[module]->fPedestalU[this_strip] = this_mean;
                fModules[module]->fPedRMSU[this_strip]   = this_rms;
              } else {
                fModules[module]->fPedestalV[this_strip] = this_mean;
                fModules[module]->fPedRMSV[this_strip]   = this_rms;
              }
            }
          }
        }
      }
    }
  } // loop over modules
}

//____________________________________________________________
void THcLADGEM::LoadCM() {

  std::ifstream cmfile(fCMFilename);
  if (!cmfile.good()) {
    std::cout << "THcLADGEM Warning: could not find cm file " << fCMFilename.c_str() << std::endl;
    return;
  } else {
    std::cout << "THcLADGEM: Read CM file " << fCMFilename.c_str() << std::endl;
  }

  std::map<int, std::map<int, std::map<int, double>>> CMMean;
  std::map<int, std::map<int, std::map<int, double>>> CMRMS;

  int crate, slot, mpd, adc_ch;
  double mean, rms;

  std::string thisline;
  while (std::getline(cmfile, thisline)) {
    if (cmfile.eof())
      break;

    std::istringstream is(thisline);

    is >> crate >> slot >> mpd >> adc_ch >> mean >> rms;

    int index = adc_ch + 16 * mpd;

    CMMean[crate][slot][index] = mean;
    CMRMS[crate][slot][index]  = rms;
  }

  for (int module = 0; module < fNModules; module++) {
    for (auto it = fModules[module]->fMPDmap.begin(); it != fModules[module]->fMPDmap.end(); ++it) {

      int this_crate = it->crate;
      int this_index = it->adc_id + 16 * it->mpd_id;
      int this_slot  = it->slot;
      int this_apv   = it->pos;

      if (CMMean.find(this_crate) != CMMean.end()) {
        if (CMMean[this_crate].find(this_slot) != CMMean[this_crate].end()) {
          if (CMMean[this_crate][this_slot].find(this_index) != CMMean[this_crate][this_slot].end()) {

            double this_mean = CMMean[this_crate][this_slot][this_index];
            double this_rms  = CMRMS[this_crate][this_slot][this_index];

            if (it->axis == LADGEM::kUaxis) {
              fModules[module]->fCommonModeMeanU[this_apv] = this_mean;
              fModules[module]->fCommonModeRMSU[this_apv]  = this_rms;
            } else {
              fModules[module]->fCommonModeMeanV[this_apv] = this_mean;
              fModules[module]->fCommonModeRMSV[this_apv]  = this_rms;
            }
          }
        }
      }
    }
  } // module loop
}

Int_t THcLADGEM::End(THaRunBase* r) {
  UInt_t runnum = r->GetNumber();
  std::cout << "THcLADGEM::End() fPedestalMode: " << fPedestalMode << std::endl; 
  if( fPedestalMode ){
    TString fname_dbped, fname_daqped, fname_dbcm, fname_daqcm;
    
    
    TString specname = GetApparatus()->GetName();
    TString detname = GetName();
    
    fname_dbped.Form( "db_ped_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_dbcm.Form( "db_cmr_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_daqped.Form( "daq_ped_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    fname_daqcm.Form( "daq_cmr_%s_%s_run%d.dat", specname.Data(), detname.Data(), runnum );
    std::cout<<"\n\n\n"<<specname<<" "<<detname<<"\n\n\n";
    
    fCMfile_dbase.open( fname_dbcm.Data() );
    fpedfile_daq.open( fname_daqped.Data() );
    fCMfile_daq.open( fname_daqcm.Data() );
    
    
    TString sdate = r->GetDate().AsString();
    sdate.Prepend( "#" );
    
    TString message;

    message.Form( "# Copy file into sbs-onl@sbsvtp#:~/cfg/pedestals for online pedestal subtraction" );
    fCMfile_daq << sdate << std::endl;
    fCMfile_daq << message << std::endl;
    fCMfile_daq << "# format = crate, slot, mpd, adc_ch, CM min, CM max"
    		  << std::endl;

    message.Form( "# Copy this file into $DB_DIR/gemped to use these pedestals for analysis");
    fCMfile_dbase << sdate << std::endl;
    fCMfile_dbase << message << std::endl;
    fCMfile_dbase << "# format = crate, slot, mpd, adc_ch, CM mean, CM RMS"
		  << std::endl;
    
    message.Form( "# Copy file into sbs-onl@sbsvtp#:~/cfg/pedestals for online pedestal subtraction" );
    
    fpedfile_daq << sdate << std::endl;
    fpedfile_daq <<  message << std::endl;
    fpedfile_daq << "# format = APV        crate       slot       mpd_id       adc_ch followed by " << std::endl
     		 << "# APV channel number      pedestal mean      pedestal rms " << std::endl;

  }

  for (auto &module : fModules) {
    if( fPedestalMode ) { module->PrintPedestals(fCMfile_dbase, fpedfile_daq, fCMfile_daq); }
    module->End(r);
  }
  return 0;
}

//____________________________________________________________

ClassImp(THcLADGEM)
