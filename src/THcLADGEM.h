#ifndef THcLADGEM_h
#define THcLADGEM_h

#include "THaNonTrackingDetector.h"
#include "THcHitList.h"
#include "THcLADGEMModule.h"
#include "THcLADGEMTrack.h"
#include "THcReactionPoint.h"

class GEM2DHits;

struct ClusterOutputData {
  // handle 1d cluster output variables
  std::vector<Int_t> layer;
  std::vector<Int_t> imodule;
  std::vector<Int_t> mpdid;
  std::vector<Int_t> axis;
  std::vector<Int_t> nstrip;
  std::vector<Int_t> maxstrip;
  std::vector<Int_t> clindex;
  std::vector<Double_t> adc; // adc sum
  std::vector<Double_t> time;
  std::vector<Double_t> pos;
  std::vector<Double_t> mpos;    // max strip pos
  std::vector<Double_t> maxsamp; // max time sample
  std::vector<Double_t> maxadc;  // max adc

  void clear() {
    layer.clear();
    imodule.clear();
    mpdid.clear();
    axis.clear();
    nstrip.clear();
    maxstrip.clear();
    clindex.clear();
    adc.clear();
    time.clear();
    pos.clear();
    mpos.clear();
    maxsamp.clear();
    maxadc.clear();
  }
};

class THcLADGEM : public THaNonTrackingDetector, public THcHitList {
public:
  THcLADGEM(const char *name, const char *description = "", THaApparatus *apparatus = nullptr);
  virtual ~THcLADGEM();

  virtual Int_t Decode(const THaEvData &);
  virtual EStatus Init(const TDatime &date);
  Int_t Begin( THaRunBase* run ) ;
  virtual Int_t CoarseProcess(TClonesArray &tracks);
  virtual Int_t FineProcess(TClonesArray &tracks);
  virtual void Clear(Option_t *opt = "");
  virtual Int_t End(THaRunBase *r = nullptr);
  Int_t GetNTracks() const { return fNTracks; }
  TClonesArray *GetTracks() const { return fGEMTracks; }

protected:
  virtual Int_t DefineVariables(EMode mode = kDefine);
  virtual Int_t ReadDatabase(const TDatime &date);

  void LoadPedestals();
  void LoadCM();
  void RotateToLab(Double_t angle, TVector3 &vect);

  static const Int_t MAXTRACKS = 500; // Was 100

  Int_t fNModules; // total number of modules
  Int_t fNLayers;  // total number of GEM layers
  Int_t fNhits;
  Int_t fNClusters;
  Int_t fNTracks;

  std::string fPedFilename;
  std::string fCMFilename;

  Double_t fD0Cut;

  TClonesArray *fGEMTracks;

  THcReactionPoint *fVertexModule; // Vertex module
  // pointer to global var indicatiing whether this spectrometer is triggered
  // for this event
  Bool_t *fPresentP;

  Bool_t fPedestalMode;

  vector<THcLADGEMModule *> fModules;

  bool fModulesInitialized;

  // Cluster data output handling
  ClusterOutputData fClusOutData;

  std::vector<vector<GEM2DHits>> f2DHits;

  // Hit output variables
  Int_t fNlayers_hit;
  Int_t fNlayers_hitU;
  Int_t fNlayers_hitV;
  Int_t fNlayers_hitUV;

  std::vector<Int_t> fNstripsU_layer;
  std::vector<Int_t> fNstripsV_layer;
  std::vector<Int_t> fNclustU_layer;
  std::vector<Int_t> fNclustV_layer;

  // spacepoint data vectors (for writing out only)
  // LHE: It's dumb to have separate vector for both layers, but def file cut based on hit wasn't working
  // Layer 0 spacepoint data
  int nhits;
  std::vector<Double_t> fPosX;
  std::vector<Double_t> fPosY;
  std::vector<Double_t> fPosZ;
  std::vector<Double_t> fTimeMean;
  std::vector<Double_t> fADCMean;
  std::vector<Double_t> fADCAsym;
  std::vector<Double_t> fTimeDiff;
  std::vector<Double_t> fTimeCorr;
  std::vector<Int_t> fIsGoodHit;
  std::vector<Int_t> fClusID0;
  std::vector<Int_t> fClusID1;
  std::vector<Int_t> fSPID;
  std::vector<Int_t> fLayer;

  bool fIgnoreVertex;
  std::ofstream fpedfile_dbase, fCMfile_dbase, fpedfile_daq, fCMfile_daq, fCMbiasfile_dbase;

public:
  GEM2DHits* Add2DHits(Int_t ilayer, Double_t x, Double_t y, Double_t z, Double_t t, Double_t dt, Double_t tc, Bool_t goodhit,
                 Double_t adc, Double_t adcasy, Int_t clust_id1, Int_t clust_id2, Int_t sp_index);
  std::vector<GEM2DHits> Get2DHits(int layer) { return f2DHits[layer]; }

  void setIgnoreVertex(bool val) { fIgnoreVertex = val; }

  ClassDef(THcLADGEM, 0)
};

#endif /* THcLADGEM_h */
