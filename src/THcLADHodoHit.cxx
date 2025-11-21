#include "THcLADHodoHit.h"

#include <iostream>

//___________________________________________________________________
THcLADHodoHit::THcLADHodoHit(Int_t toptdc, Int_t btmtdc, Double_t topadc, Double_t btmadc, Int_t ipad,
                             THcLADHodoPlane *hp)
    : fTopTDC(toptdc), fBtmTDC(btmtdc), fTopADC_Ped(topadc), fBtmADC_Ped(btmadc), fPaddleNumber(ipad),fTwoGoodTimes(kFALSE), fPlane(hp) {}

//___________________________________________________________________
THcLADHodoHit::~THcLADHodoHit() {}

void THcLADHodoHit::SetCorrectedTimes(Double_t top, Double_t btm, Double_t timeave) {
  fTopCorrectedTime = top;
  fBtmCorrectedTime = btm;
  fScinCorrectedTime = timeave;
  fHasCorrectedTimes = kTRUE;
}

void THcLADHodoHit::SetCorrectedTimes(Double_t top, Double_t btm) {
  fTopCorrectedTime = top;
  fBtmCorrectedTime = btm;
  fHasCorrectedTimes = kFALSE;
}

void THcLADHodoHit::SetCorrectedTimes_FADC(Double_t top, Double_t btm, Double_t timeave) {
  fTopADC_Time = top;
  fBtmADC_Time = btm;
  fScinCorrectedTime_FADC = timeave;
  fHasCorrectedTimes = kTRUE;
}
void THcLADHodoHit::SetCorrectedTimes_FADC(Double_t top, Double_t btm) {
  fTopADC_Time = top;
  fBtmADC_Time = btm;
  fHasCorrectedTimes = kFALSE;
}

void THcLADHodoHit::SetTOFCorrectedTimes(Double_t top, Double_t btm, Double_t timeave) {
  fTopTOFCorrectedTime = top;
  fBtmTOFCorrectedTime = btm;
  fScinTOFCorrectedTime = timeave;
  fHasCorrectedTimes = kTRUE;
}

ClassImp(THcLADHodoHit)
