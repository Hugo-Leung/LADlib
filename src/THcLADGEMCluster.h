#ifndef THcLADGEMCluster_h
#define THcLADGEMCluster_h

#include "TObject.h"

class THcLADGEMCluster : public TObject {
 public:
  THcLADGEMCluster();
  virtual ~THcLADGEMCluster() {}
  virtual void Clear( Option_t* opt );
  virtual void Print( Option_t* opt="" ) const;
  /*
  void AddHit(THcLADGEMHit* hit) {
    fHits.push_back(hit);
    fNHits++;
  }

  THcLADGEMHits* GetHit(Int_t ihit) { return fHits[ihit]; }
  */

  Int_t    GetNHits() { return fNHits; }
  Int_t    GetNStrips() { return fNStrips; }
  Int_t    GetStripLow() { return fStripLow;}
  Int_t    GetStripHigh() { return fStripHigh; }
  Int_t    GetStripMax() { return fStripMax; }

  Double_t GetADCsum()   { return fADCsum; }
  Double_t GetADCMax()   { return fADCMax; }
  Double_t GetPos()      { return fPos; }
  Double_t GetPosMax()   { return fPosMax; }
  Double_t GetPosSigma() { return fPosSigma; }
  Double_t GetMoments()  { return fMoments; }
  Double_t GetPosDiff()  { return fPosDiff; }
  Double_t GetE()        { return fE; }
  Double_t GetTime()     { return fTime; }
  Double_t GetTimeFit()  { return fTimeFit; }
  Double_t GetTDeconv()  { return fTDeconv; }
  Double_t GetADCsumDeconv() { return fADCsumDeconv; }
  Double_t GetADCsumDeconvMaxCombo() { return fADCsumDeconvMaxCombo; }
  const std::vector<Double_t>& GetStripADCsum() { return stripADCsum; }
  const std::vector<Double_t>& GetDeconvADCsum() { return DeconvADCsum; }
  const std::vector<Double_t>& GetADCsamples() { return ADCsamples; }
  const std::vector<Double_t>& GetDeconvADCsamples() { return DeconvADCsamples; }
  Int_t    GetSampMaxDeconv()  { return fSampMaxDeconv; }
  Int_t    GetComboMaxDeconv()  { return fComboMaxDeconv; }
  Int_t    GetSampMax()  { return fSampMax; }

  Double_t GetMPD() { return fMPD; }
  Double_t GetAPV() { return fAPV; }
  Double_t GetLayer() { return fLayer; }
  Double_t GetAxis() { return fAxis; }
  Int_t    GetCLIndex() { return fCLIndex; }
  std::vector<UInt_t> GetHitIndex() { return hitindex; }
  Int_t    GetRawStrip() { return rawstrip; }

  void SetMode(int this_value) { fClusteringFlag = this_value; }
  void SetNStrips(int this_value) { fNStrips = this_value; }
  void SetStripLow(int this_value) { fStripLow = this_value; }
  void SetStripHigh(int this_value) { fStripHigh = this_value; }
  void SetStripMax(int this_value) { fStripMax = this_value; }  
  void SetStrips(int nstrip, int striplo, int striphi, int stripmax)
  { fNStrips = nstrip, fStripLow = striplo, fStripHigh = striphi, fStripMax = stripmax; }

  void SetPosition(double this_value) { fPos = this_value; }
  void SetPosMaxStrip(double this_value) { fPosMax = this_value; }
  void SetPosSigma(double this_value) { fPosSigma = this_value; }
  void SetMoments(double this_value) { fMoments = this_value; }

  void SetTime(double this_value) { fTime = this_value; }
  void SetTimeFit(double this_value) { fTimeFit = this_value; }
  void SetTime(double value1, double value2) { fTime = value1; fTSigma = value2; }
  void SetTimeDeconv(double this_value) { fTDeconv = this_value; }
  void SetADCsum(double this_value) { fADCsum = this_value; }
  void SetADCMax(double this_value) { fADCMax = this_value; }
  void SetSampMax(int this_value) { fSampMax = this_value; }
  void SetSampMaxDeconv(int this_value) { fSampMaxDeconv = this_value; }
  void SetComboMaxDeconv(int this_value) { fComboMaxDeconv = this_value; }
  void SetADCsamples( const std::vector<Double_t> &samples ) { ADCsamples = samples; }
  void SetDeconvADCsamples( const std::vector<Double_t> &samples ) { DeconvADCsamples = samples; }
  void SetADCsumDeconv( double this_value ) { fADCsumDeconv = this_value; }
  void SetADCsumDeconvMaxCombo( double this_value ) { fADCsumDeconvMaxCombo = this_value; }

  void SetLayer(int this_value ) { fLayer = this_value; }
  void SetMPD(int this_value ) { fMPD = this_value; }
  void SetAPV(int this_value ) { fAPV = this_value; }
  void SetAxis(int this_value) { fAxis = this_value; }
  void SetCLIndex(int this_value) { fCLIndex = this_value; }
  void SetHitIndex( const std::vector<UInt_t> &indices ) { hitindex = indices; }
  void SetRawStrip(int this_value) { rawstrip = this_value; }
  void SetStripADCsum( const std::vector<Double_t> &adc_sums ) { stripADCsum = adc_sums; }
  void AddStripADCsum( Double_t adc_sum ) { stripADCsum.push_back( adc_sum ); }
  void SetDeconvADCsum( const std::vector<Double_t> &adc_sums ) { DeconvADCsum = adc_sums; }  
  void AddDeconvADCsum( Double_t adc_sum ) { DeconvADCsum.push_back( adc_sum ); }
 protected:

  Int_t    fClusteringFlag;

  Int_t    fNStrips;
  Int_t    fStripLow;
  Int_t    fStripHigh;
  Int_t    fStripMax;

  Int_t    fLayer;
  Int_t    fMPD;
  Int_t    fAPV; // APV adc id
  Int_t    fAxis; // U/V, X/Y 
  Int_t    fCLIndex; // cluster index, associated with a particular 2D hit

  Double_t fPos; // ADC weighted mean coordinate along the direction measured by the strip
  Double_t fPosMax;  // Max strip position
  Double_t fPosDiff;
  Double_t fPosSigma; // RMS coordinate deviation from the mean 
  Double_t fMoments;      // clustermoments
  Int_t    fSampMax; // time sample in which peak of cluster summed ADC values occurs
  Int_t    fSampMaxDeconv; //time sample in which the deconvoluted cluster-summed ADC samples peaks.
  Int_t    fComboMaxDeconv; //2nd time sample of max two-sample combo for cluster-summed deconvoluted ADC samples
  std::vector<Double_t> ADCsamples; //cluster-summed ADC samples (accounting for split fraction)
  std::vector<Double_t> DeconvADCsamples; //cluster-summed deconvoluted ADC samples (accounting for split fraction)

  Double_t fADCsum;
  Double_t fADCMax;
  Double_t fADCsumDeconv; //sum of all deconvoluted ADC samples over all strips in the cluster
  Double_t fADCsumDeconvMaxCombo; //sum over all strips in the cluster of max two-sample combo
  std::vector<Double_t> stripADCsum; //Sum of individual strip ADCs over all samples on all strips; accounting for split fraction
  std::vector<Double_t> DeconvADCsum; //Sum of individual deconvoluted ADC samples over all samples on all strips; accounting for split fraction

  Int_t    fNHits;
  Double_t fE;
  Double_t fTime;   // mean time
  Double_t fTSigma; 
  Double_t fTimeFit;
  Double_t fTDeconv; //cluster-summed mean deconvoluted hit time.

  //  std::vector<THcLADGEMHit*> fHits;
  std::vector<UInt_t> hitindex; //position in decoded hit array of each strip in the cluster:
  Int_t rawstrip; //Raw APV strip number before decoding 
  Int_t rawMPD; //Raw MPD number before decoding 
  Int_t rawAPV; //Raw APV number before decoding 


  ClassDef(THcLADGEMCluster,0)
};

#endif
