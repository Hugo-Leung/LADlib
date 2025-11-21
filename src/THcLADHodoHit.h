#ifndef THcLADHodoHit_h
#define THcLADHodoHit_h

#include "THcLADHodoPlane.h"
#include "THcRawHodoHit.h"
#include "TObject.h"

class THcLADHodoHit : public TObject {
public:
  THcLADHodoHit(Int_t toptdc, Int_t btmtdc, Double_t topadc, Double_t btmadc, Int_t ipad, THcLADHodoPlane *hp);
  virtual ~THcLADHodoHit();

  Double_t GetTopADC() const { return fTopADC_Ped; }
  Double_t GetBtmADC() const { return fBtmADC_Ped; }
  Double_t GetTopADCpeak() const { return fTopADC_Peak; }
  Double_t GetBtmADCpeak() const { return fBtmADC_Peak; }
  Double_t GetTopADCtime() const { return fTopADC_Time; }
  Double_t GetBtmADCtime() const { return fBtmADC_Time; }
  Double_t GetTopADCCorrtime() const { return fTopADC_CorrTime; }
  Double_t GetBtmADCCorrtime() const { return fBtmADC_CorrTime; }
  Double_t GetCalcPosition() const { return fCalcPosition; }
  Double_t GetCalcPosition_FADC() const { return fCalcPosition_FADC; }
  Int_t GetTopTDC() const { return fTopTDC; }
  Int_t GetBtmTDC() const { return fBtmTDC; }
  Double_t GetTopCorrectedTime() const { return fTopCorrectedTime; }
  Double_t GetBtmCorrectedTime() const { return fBtmCorrectedTime; }
  Double_t GetScinCorrectedTime() const { return fScinCorrectedTime; }
  Double_t GetScinCorrectedTime_FADC() const { return fScinCorrectedTime_FADC; }
  Bool_t GetTwoGoodTimes() const { return fTwoGoodTimes; }
  Bool_t GetHasCorrectedTimes() const { return fHasCorrectedTimes; }
  Int_t GetPaddleNumber() const { return fPaddleNumber; }
  Int_t GetPaddleCenter() const { return fPaddleCenter; }
  Double_t GetPaddleADC() const { return sqrt(fTopADC_Ped * fBtmADC_Ped); }
  Double_t GetPaddleADCpeak() const { return sqrt(fTopADC_Peak * fBtmADC_Peak); }
  Double_t GetPaddleADC_MeV() const { return fPlane->GetEdep2MeV_int(fPaddleNumber) * GetPaddleADC(); }
  Double_t GetPaddleADCpeak_MeV() const { return fPlane->GetEdep2MeV_amp(fPaddleNumber) * GetPaddleADCpeak(); }
  Double_t GetTopTOFCorrectedTime() const { return fTopTOFCorrectedTime; }
  Double_t GetBtmTOFCorrectedTime() const { return fBtmTOFCorrectedTime; }
  Double_t GetScinTOFCorrectedTime() const { return fScinTOFCorrectedTime; }

  void SetPaddleCenter(Double_t padcenter) { fPaddleCenter = padcenter; }
  void SetTopADCpeak(Double_t adc) { fTopADC_Peak = adc; }
  void SetBtmADCpeak(Double_t adc) { fBtmADC_Peak = adc; }
  void SetTopADCtime(Double_t ptime) { fTopADC_Time = ptime; }
  void SetTopADCCorrtime(Double_t ptime) { fTopADC_CorrTime = ptime; }
  void SetBtmADCCorrtime(Double_t ptime) { fBtmADC_CorrTime = ptime; }
  void SetCalcPosition(Double_t calcpos) { fCalcPosition = calcpos; }
  void SetCalcPosition_FADC(Double_t calcpos) { fCalcPosition_FADC = calcpos; }
  void SetBtmADCtime(Double_t ptime) { fBtmADC_Time = ptime; }
  void SetCorrectedTimes(Double_t top, Double_t btm);
  void SetCorrectedTimes(Double_t top, Double_t btm, Double_t timeave);
  void SetCorrectedTimes_FADC(Double_t top, Double_t btm);
  void SetCorrectedTimes_FADC(Double_t top, Double_t btm, Double_t timeave);
  void SetTOFCorrectedTimes(Double_t top, Double_t btm, Double_t timeave);

protected:
  Int_t fTopTDC;
  Int_t fBtmTDC;
  Double_t fTopADC_Ped; // Pedestal subtracted ADC
  Double_t fBtmADC_Ped;
  Double_t fTopADC_Peak; // ADC peak amplitude
  Double_t fBtmADC_Peak; // ADC peak amplitude
  Double_t fTopADC_Time;
  Double_t fBtmADC_Time;
  Double_t fTopADC_CorrTime; // ADC time
  Double_t fBtmADC_CorrTime; // ADC time
  Double_t fCalcPosition;    // Position along paddle calculated by time diff
  Double_t fCalcPosition_FADC; // Position along paddle calculated by time diff using FADC times
  Int_t fPaddleNumber;

  Double_t fTopCorrectedTime;  // Pulse height corrected time
  Double_t fBtmCorrectedTime;  // Pulse height corrected time
  Double_t fScinCorrectedTime; // Time average corrected for position
  // based on ADCs.
  Double_t fScinCorrectedTime_FADC; // Time average corrected for position
                                   // based on ADCs.

  Double_t fTopTOFCorrectedTime;  // Times corrected for z position
  Double_t fBtmTOFCorrectedTime;  // using nominal beta
  Double_t fScinTOFCorrectedTime; // Times corrected for z position
                                  // using nominal beta

  Bool_t fHasCorrectedTimes;
  Bool_t fTwoGoodTimes;
  Double_t fPaddleCenter;

  THcLADHodoPlane *fPlane;

private:
  THcLADHodoHit(const THcLADHodoHit &);
  THcLADHodoHit &operator=(const THcLADHodoHit &);

  ClassDef(THcLADHodoHit, 0)
};

#endif
