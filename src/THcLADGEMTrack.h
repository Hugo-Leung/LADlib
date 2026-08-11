#ifndef THcLADGEMTrack_h
#define THcLADGEMTrack_h

// LAD GEM Track Object

#include "DataType.h"
#include "THcGoodLADHit.h"
#include "TObject.h"
#include "TVector3.h"

class GEM2DHits {
public:
  GEM2DHits()
      : layer(-1), posX(kBig), posY(kBig), posZ(kBig), posX_local(kBig), posY_local(kBig), TimeMean(kBig),
        TimeDiff(kBig), TimeCorr(kBig), IsGoodHit(kFALSE), ADCMean(kBig), ADCasym(kBig), spID(-1), clusID{-1, -1},
        corrcoeff(0.0), corrcoeff_deconv(0.0), corrcoeff_strip(0.0), corrcoeff_strip_deconv(0.0), trackID(-1) {}
  GEM2DHits(const GEM2DHits &other)
      : layer(other.layer), posX(other.posX), posY(other.posY), posZ(other.posZ), posX_local(other.posX_local),
        posY_local(other.posY_local), TimeMean(other.TimeMean), TimeDiff(other.TimeDiff), TimeCorr(other.TimeCorr),
        IsGoodHit(other.IsGoodHit), ADCMean(other.ADCMean), ADCasym(other.ADCasym), spID(other.spID),
        corrcoeff(other.corrcoeff), corrcoeff_deconv(other.corrcoeff_deconv), corrcoeff_strip(other.corrcoeff_strip),
        corrcoeff_strip_deconv(other.corrcoeff_strip_deconv), trackID(other.trackID) {
    clusID[0] = other.clusID[0];
    clusID[1] = other.clusID[1];
  }
  virtual ~GEM2DHits() = default;
  void Set(Int_t _layer, Double_t _x, Double_t _y, Double_t _z, Double_t _t, Double_t _dt, Double_t _ct, Bool_t _flag,
           Double_t _adc, Double_t _asym) {
    layer      = _layer;
    posX       = _x;
    posY       = _y;
    posZ       = _z;
    posX_local = _x;
    posY_local = _y;
    TimeMean   = _t;
    TimeDiff   = _dt;
    TimeCorr   = _ct;
    IsGoodHit  = _flag;
    ADCMean    = _adc;
    ADCasym    = _asym;
  }
  void SetClusterIDs(Int_t id1, Int_t id2) {
    clusID[0] = id1;
    clusID[1] = id2;
  }
  void SetSPID(Int_t id) { spID = id; }

  Int_t layer;
  Double_t posX;
  Double_t posY;
  Double_t posZ;
  Double_t posX_local;
  Double_t posY_local;
  Double_t TimeMean; // average time
  Double_t TimeDiff;
  Double_t TimeCorr;
  Bool_t IsGoodHit;
  Double_t ADCMean; // average adc sum
  Double_t ADCasym;
  Int_t spID; // this Space point's ID
  Int_t trackID;
  Int_t clusID[2]; // associated cluster IDs
  double corrcoeff;
  double corrcoeff_deconv;
  double corrcoeff_strip;
  double corrcoeff_strip_deconv;
};

class THcLADGEMTrack : public TObject {
public:
  THcLADGEMTrack(Int_t nlayers);
  virtual ~THcLADGEMTrack();
  virtual void Clear(Option_t *opt = "");

  // THcLADGEMHit* GetHit(Int_t ihit) {};

  Int_t GetTrackID() const { return fTrackID; }
  // Space points
  Int_t GetNSpacePoints() const { return fNSp; }
  Int_t *GetClusID_Sp1() const { return fSp[0].clusID; }
  Int_t *GetClusID_Sp2() const { return fSp[1].clusID; }
  Double_t GetADCMean_Sp1() const { return fSp[0].ADCMean; }
  Double_t GetADCMean_Sp2() const { return fSp[1].ADCMean; }
  Double_t GetADCasym_Sp1() const { return fSp[0].ADCasym; }
  Double_t GetADCasym_Sp2() const { return fSp[1].ADCasym; }
  Double_t GetX1() const { return fSp[0].posX; }
  Double_t GetY1() const { return fSp[0].posY; }
  Double_t GetZ1() const { return fSp[0].posZ; }
  Double_t GetX2() const { return fSp[1].posX; }
  Double_t GetY2() const { return fSp[1].posY; }
  Double_t GetZ2() const { return fSp[1].posZ; }
  Double_t GetX1_local() const { return fSp[0].posX_local; }
  Double_t GetY1_local() const { return fSp[0].posY_local; }
  Double_t GetX2_local() const { return fSp[1].posX_local; }
  Double_t GetY2_local() const { return fSp[1].posY_local; }
  // Track quantities
  Double_t GetProjVz() const { return fProjVz; } // projected z-vertex
  Double_t GetProjVy() const { return fProjVy; } // projected y-vertex
  Double_t GetProjVx() const { return fProjVx; } // projected x-vertex
  Double_t GetD0() const { return fD0; }
  Bool_t GetGoodD0() const { return fhasGoodD0; }
  Double_t GetT() const { return fT; }
  Double_t GetdT() const { return fTdiff; }
  short GetHasHodoHit() const { return fHasHodoHit; }
  THcGoodLADHit *GetBestHodoHit() const { return fBestHodoHit; }
  void SetBestHodoHit(THcGoodLADHit *hit) { fBestHodoHit = hit; }

  void SetHasHodoHit(short hasHodoHit) { fHasHodoHit = hasHodoHit; }
  void SetTrackID(int itrk) { fTrackID = itrk; }
  void SetD0(Double_t d0) { fD0 = d0; }
  void SetGoodD0(Bool_t good) { fhasGoodD0 = good; }
  void SetZVertex(Double_t vz) { fProjVz = vz; }
  void SetYVertex(Double_t vy) { fProjVy = vy; }
  void SetXVertex(Double_t vx) { fProjVx = vx; }
  void SetProjVertex(Double_t vx, Double_t vy, Double_t vz) {
    fProjVx = vx;
    fProjVy = vy;
    fProjVz = vz;
  }
  void SetAngles(Double_t theta, Double_t phi) {
    ftheta = theta;
    fphi   = phi;
  }
  void GetAngles(Double_t &theta, Double_t &phi) const {
    theta = ftheta;
    phi   = fphi;
  }
  Double_t GetTheta() const { return ftheta; }
  Double_t GetPhi() const { return fphi; }
  void SetChisq(Double_t c) { chisq = c; }
  Double_t GetChisq() const { return chisq; }
  void SetTime(Double_t t, Double_t dt) {
    fT     = t;
    fTdiff = dt;
  }

  bool IsGoodTrack() const { return fIsGoodTrack; }
  void SetIsGoodTrack(bool isGood) { fIsGoodTrack = isGood; }

  // ---- no-vertex tracking (fit does NOT use the target vertex) ----
  void SetChisq_noTrackVertex(Double_t c) { chisq_noTrackVertex = c; }
  Double_t GetChisq_noTrackVertex() const { return chisq_noTrackVertex; }
  void SetAngles_noTrackVertex(Double_t theta, Double_t phi) {
    ftheta_noTrackVertex = theta;
    fphi_noTrackVertex   = phi;
  }
  void GetAngles_noTrackVertex(Double_t &theta, Double_t &phi) const {
    theta = ftheta_noTrackVertex;
    phi   = fphi_noTrackVertex;
  }
  Double_t GetTheta_noTrackVertex() const { return ftheta_noTrackVertex; }
  Double_t GetPhi_noTrackVertex() const { return fphi_noTrackVertex; }
  void SetProjVertex_noTrackVertex(Double_t vx, Double_t vy, Double_t vz) {
    fProjVx_noTrackVertex = vx;
    fProjVy_noTrackVertex = vy;
    fProjVz_noTrackVertex = vz;
  }
  Double_t GetProjVx_noTrackVertex() const { return fProjVx_noTrackVertex; }
  Double_t GetProjVy_noTrackVertex() const { return fProjVy_noTrackVertex; }
  Double_t GetProjVz_noTrackVertex() const { return fProjVz_noTrackVertex; }
  void SetD0_noTrackVertex(Double_t d0) { fD0_noTrackVertex = d0; }
  Double_t GetD0_noTrackVertex() const { return fD0_noTrackVertex; }
  void SetHasHodoHit_noTrackVertex(short hasHodoHit) { fHasHodoHit_noTrackVertex = hasHodoHit; }
  short GetHasHodoHit_noTrackVertex() const { return fHasHodoHit_noTrackVertex; }
  void SetBestHodoHit_noTrackVertex(THcGoodLADHit *hit) { fBestHodoHit_noTrackVertex = hit; }
  THcGoodLADHit *GetBestHodoHit_noTrackVertex() const { return fBestHodoHit_noTrackVertex; }
  void SetIsGoodTrack_noTrackVertex(bool isGood) { fIsGoodTrack_noTrackVertex = isGood; }
  bool IsGoodTrack_noTrackVertex() const { return fIsGoodTrack_noTrackVertex; }

  // ---- x-z tracking (vertex-constrained fit; chi-square ignores the y position) ----
  void SetChisq_xz(Double_t c) { chisq_xz = c; }
  Double_t GetChisq_xz() const { return chisq_xz; }
  void SetAngles_xz(Double_t theta, Double_t phi) {
    ftheta_xz = theta;
    fphi_xz   = phi;
  }
  void GetAngles_xz(Double_t &theta, Double_t &phi) const {
    theta = ftheta_xz;
    phi   = fphi_xz;
  }
  Double_t GetTheta_xz() const { return ftheta_xz; }
  Double_t GetPhi_xz() const { return fphi_xz; }
  void SetProjVertex_xz(Double_t vx, Double_t vy, Double_t vz) {
    fProjVx_xz = vx;
    fProjVy_xz = vy;
    fProjVz_xz = vz;
  }
  Double_t GetProjVx_xz() const { return fProjVx_xz; }
  Double_t GetProjVy_xz() const { return fProjVy_xz; }
  Double_t GetProjVz_xz() const { return fProjVz_xz; }
  void SetD0_xz(Double_t d0) { fD0_xz = d0; }
  Double_t GetD0_xz() const { return fD0_xz; }
  void SetHasHodoHit_xz(short hasHodoHit) { fHasHodoHit_xz = hasHodoHit; }
  short GetHasHodoHit_xz() const { return fHasHodoHit_xz; }
  void SetBestHodoHit_xz(THcGoodLADHit *hit) { fBestHodoHit_xz = hit; }
  THcGoodLADHit *GetBestHodoHit_xz() const { return fBestHodoHit_xz; }
  void SetIsGoodTrack_xz(bool isGood) { fIsGoodTrack_xz = isGood; }
  bool IsGoodTrack_xz() const { return fIsGoodTrack_xz; }

  // ---- x-z tracking, no vertex (free-line fit; chi-square ignores the y position) ----
  void SetChisq_noTrackVertex_xz(Double_t c) { chisq_noTrackVertex_xz = c; }
  Double_t GetChisq_noTrackVertex_xz() const { return chisq_noTrackVertex_xz; }
  void SetAngles_noTrackVertex_xz(Double_t theta, Double_t phi) {
    ftheta_noTrackVertex_xz = theta;
    fphi_noTrackVertex_xz   = phi;
  }
  void GetAngles_noTrackVertex_xz(Double_t &theta, Double_t &phi) const {
    theta = ftheta_noTrackVertex_xz;
    phi   = fphi_noTrackVertex_xz;
  }
  Double_t GetTheta_noTrackVertex_xz() const { return ftheta_noTrackVertex_xz; }
  Double_t GetPhi_noTrackVertex_xz() const { return fphi_noTrackVertex_xz; }
  void SetProjVertex_noTrackVertex_xz(Double_t vx, Double_t vy, Double_t vz) {
    fProjVx_noTrackVertex_xz = vx;
    fProjVy_noTrackVertex_xz = vy;
    fProjVz_noTrackVertex_xz = vz;
  }
  Double_t GetProjVx_noTrackVertex_xz() const { return fProjVx_noTrackVertex_xz; }
  Double_t GetProjVy_noTrackVertex_xz() const { return fProjVy_noTrackVertex_xz; }
  Double_t GetProjVz_noTrackVertex_xz() const { return fProjVz_noTrackVertex_xz; }
  void SetD0_noTrackVertex_xz(Double_t d0) { fD0_noTrackVertex_xz = d0; }
  Double_t GetD0_noTrackVertex_xz() const { return fD0_noTrackVertex_xz; }
  void SetHasHodoHit_noTrackVertex_xz(short hasHodoHit) { fHasHodoHit_noTrackVertex_xz = hasHodoHit; }
  short GetHasHodoHit_noTrackVertex_xz() const { return fHasHodoHit_noTrackVertex_xz; }
  void SetBestHodoHit_noTrackVertex_xz(THcGoodLADHit *hit) { fBestHodoHit_noTrackVertex_xz = hit; }
  THcGoodLADHit *GetBestHodoHit_noTrackVertex_xz() const { return fBestHodoHit_noTrackVertex_xz; }
  void SetIsGoodTrack_noTrackVertex_xz(bool isGood) { fIsGoodTrack_noTrackVertex_xz = isGood; }
  bool IsGoodTrack_noTrackVertex_xz() const { return fIsGoodTrack_noTrackVertex_xz; }

  GEM2DHits GetSpacePoint(int isp) { return fSp[isp]; }
  virtual void AddSpacePoint(GEM2DHits sp);
  Int_t GetSpID_0() const { return fSp[0].spID; };
  Int_t GetSpID_1() const { return fSp[1].spID; };
  int GetSpacePointID_0U() { return fSp[0].clusID[0]; }
  int GetSpacePointID_0V() { return fSp[0].clusID[1]; }
  int GetSpacePointID_1U() { return fSp[1].clusID[0]; }
  int GetSpacePointID_1V() { return fSp[1].clusID[1]; }

protected:
  Int_t fNSp;
  Double_t fProjVz;
  Double_t fProjVy;
  Double_t fProjVx;
  Double_t fD0;
  Bool_t fhasGoodD0;
  Int_t fTrackID;
  Int_t fClustID0;
  Int_t fClustID1;
  Double_t fT;
  Double_t fTdiff;
  short fHasHodoHit;
  Double_t chisq;
  Double_t ftheta; // track angle between the track and z-azis in radians
  Double_t fphi;   // track angle between the track and x-azis in radians
  bool fIsGoodTrack;
  GEM2DHits *fSp;
  THcGoodLADHit *fBestHodoHit; // pointer to the best hodoscope hit associated with this track, if any

  // no-vertex tracking quantities (parallel to the vertex-constrained ones above)
  Double_t chisq_noTrackVertex;
  Double_t ftheta_noTrackVertex;
  Double_t fphi_noTrackVertex;
  Double_t fProjVx_noTrackVertex;
  Double_t fProjVy_noTrackVertex;
  Double_t fProjVz_noTrackVertex;
  Double_t fD0_noTrackVertex;
  short fHasHodoHit_noTrackVertex;
  bool fIsGoodTrack_noTrackVertex;
  THcGoodLADHit *fBestHodoHit_noTrackVertex;

  // x-z tracking quantities (vertex-constrained fit; chi-square ignores y)
  Double_t chisq_xz;
  Double_t ftheta_xz;
  Double_t fphi_xz;
  Double_t fProjVx_xz;
  Double_t fProjVy_xz;
  Double_t fProjVz_xz;
  Double_t fD0_xz;
  short fHasHodoHit_xz;
  bool fIsGoodTrack_xz;
  THcGoodLADHit *fBestHodoHit_xz;

  // x-z tracking quantities, no vertex (free-line fit; chi-square ignores y)
  Double_t chisq_noTrackVertex_xz;
  Double_t ftheta_noTrackVertex_xz;
  Double_t fphi_noTrackVertex_xz;
  Double_t fProjVx_noTrackVertex_xz;
  Double_t fProjVy_noTrackVertex_xz;
  Double_t fProjVz_noTrackVertex_xz;
  Double_t fD0_noTrackVertex_xz;
  short fHasHodoHit_noTrackVertex_xz;
  bool fIsGoodTrack_noTrackVertex_xz;
  THcGoodLADHit *fBestHodoHit_noTrackVertex_xz;

private:
  THcLADGEMTrack(const THcLADGEMTrack &);
  THcLADGEMTrack &operator=(const THcLADGEMTrack &);

  ClassDef(THcLADGEMTrack, 0)
};

#endif
