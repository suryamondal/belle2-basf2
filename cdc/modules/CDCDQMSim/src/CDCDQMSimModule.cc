#include <cdc/modules/CDCDQMSim/CDCDQMSimModule.h>

#include <cdc/geometry/CDCGeometryPar.h>
#include <cdc/dataobjects/WireID.h>
#include <framework/core/ModuleParam.templateDetails.h>

#include <TMath.h>
#include <TDirectory.h>

#include <cmath>
#include <set>

using namespace Belle2;
using namespace Belle2::CDC;

REG_MODULE(CDCDQMSim);

CDCDQMSimModule::CDCDQMSimModule() : HistoModule()
{
  setDescription("CDC DQM for simulation — same as CDCDQM but without trigger/HLT requirements.");
  addParam("MinHits",         m_minHits,         "Minimum CDC hits per event",  0);
  addParam("MinNDF",          m_minNdf,           "Minimum track NDF",          20);
  addParam("MinPt",           m_minPt,            "Minimum track pT [GeV/c] for wire efficiency", 0.0);
  addParam("AdjustWireShift", m_adjustWireShift,  "Correct phi for stereo layers", true);
  setPropertyFlags(c_ParallelProcessingCertified);
}

void CDCDQMSimModule::defineHisto()
{
  TDirectory* oldDir = gDirectory;
  oldDir->mkdir("CDC")->cd();

  m_hNEvents         = new TH1F("hNEvents",         "hNEvents",                    10,    0,    10);
  m_hNEvents->GetXaxis()->SetBinLabel(1, "number of events");
  m_hOcc             = new TH1F("hOcc",             "hOccupancy",                 150,    0,   1.5);
  m_hHit             = new TH2F("hHit",             "CDC-hits;layer;wire",         56,    0,    56,  400, 0, 400);
  m_hPhi             = new TH1F("hPhi",             "CDC track #phi (IP);#phi [deg];entries", 360, -180, 180);
  m_hPhiEff          = new TH2F("hPhiEff",          "CDC track #phi;#phi [deg];nCDChits",     360, -180, 180, 100, 0, 100);
  m_hPhiHit          = new TH2F("h2HitPhi",         "CDC hits map;#phi [deg];layer",           90, -180, 180,  56, 0,  56);
  m_hPhiNCDC         = new TH2F("hPhiNCDC",         "nCDCHits vs #phi;#phi [deg];nCDCHits",    45, -180, 180,  61, -0.5, 60.5);
  m_hADCBoard        = new TH2F("hADCBoard",        "ADC vs board;board index;ADC",           300,    0,   300, 200,    0, 1000);
  m_hADCLayer        = new TH2F("hADCLayer",        "ADC vs layer;layer index;ADC",            56,    0,    56, 200,    0, 1000);
  m_hTDC             = new TH2F("hTDC",             "TDC vs board;board index;TDC",           300,    0,   300, 1000, 4200, 5200);
  m_hTrackingWireEff = new TH2F("hTrackingWireEff",
                                "Attached vs Expected wires (backplate view);wire bin;layer",
                                400, 0.5, 400.5, 56 * 2, -0.5, 56 * 2 - 0.5);

  oldDir->cd();
}

void CDCDQMSimModule::initialize()
{
  REG_HISTOGRAM
  m_cdcHits.isOptional();
  m_Tracks.isOptional();
  m_RecoTracks.isOptional();
}

void CDCDQMSimModule::beginRun()
{
  m_nEvents = 0;
  if (m_hNEvents)         m_hNEvents->Reset();
  if (m_hOcc)             m_hOcc->Reset();
  if (m_hHit)             m_hHit->Reset();
  if (m_hPhi)             m_hPhi->Reset();
  if (m_hPhiEff)          m_hPhiEff->Reset();
  if (m_hPhiHit)          m_hPhiHit->Reset();
  if (m_hPhiNCDC)         m_hPhiNCDC->Reset();
  if (m_hADCBoard)        m_hADCBoard->Reset();
  if (m_hADCLayer)        m_hADCLayer->Reset();
  if (m_hTDC)             m_hTDC->Reset();
  if (m_hTrackingWireEff) m_hTrackingWireEff->Reset();
}

void CDCDQMSimModule::event()
{
  static CDCGeometryPar& geo = CDCGeometryPar::Instance();
  const int nWires = 14336;

  if (m_cdcHits.getEntries() < m_minHits) return;

  m_nEvents++;
  m_hOcc->Fill(static_cast<float>(m_cdcHits.getEntries()) / nWires);

  for (const auto& hit : m_cdcHits)
    m_hHit->Fill(hit.getICLayer(), hit.getIWire());

  for (const auto& b2track : m_Tracks) {
    const TrackFitResult* fit = b2track.getTrackFitResultWithClosestMass(Const::pion);
    if (!fit) continue;

    RecoTrack* rtrack = b2track.getRelatedTo<RecoTrack>();
    if (!rtrack) continue;

    const genfit::FitStatus* fs = rtrack->getTrackFitStatus();
    if (!fs) continue;

    if (std::fabs(fit->getD0()) < 1.0 && std::fabs(fit->getZ0()) < 1.0 &&
        fit->getMomentum().Rho() > m_minPt) {
      std::set<int> hitLayers;
      for (const RecoHitInformation::UsedCDCHit* hit : rtrack->getCDCHitList()) {
        if (!rtrack->getRecoHitInformation(hit)) continue;
        hitLayers.insert(hit->getICLayer());
      }
      if (hitLayers.empty()) continue;
      auto helix = fit->getHelix();
      int nSL = geo.getNumberOfSenseLayers();
      for (int lay = 0; lay < nSL; lay++) {
        double arcL = helix.getArcLength2DAtCylindricalR(geo.senseWireR(lay));
        if (std::isnan(arcL)) continue;
        auto pos = helix.getPositionAtArcLength2D(arcL);
        if (pos.Z() > geo.senseWireFZ(lay) || pos.Z() < geo.senseWireBZ(lay)) continue;
        int bin = findPhiBin(getShiftedPhi(pos, lay), lay);
        m_hTrackingWireEff->Fill(bin, lay);
        if (hitLayers.count(lay))
          m_hTrackingWireEff->Fill(bin, lay + nSL);
      }
    }

    double phi = fit->getPhi() / Unit::deg;
    m_hPhiNCDC->Fill(phi, std::min(static_cast<int>(rtrack->getNumberOfCDCHits()), 60));

    if (fs->getNdf() < m_minNdf) continue;

    if (std::fabs(fit->getD0()) < 1.0 && std::fabs(fit->getZ0()) < 1.0) {
      m_hPhi->Fill(phi);
      double nsvd = fit->getHitPatternVXD().getNSVDHits();
      double ncdc = fit->getHitPatternCDC().getNHits();
      if (nsvd > 6)
        m_hPhiEff->Fill(phi, std::min(ncdc, 99.5));
    }

    for (const RecoHitInformation::UsedCDCHit* hit : rtrack->getCDCHitList()) {
      if (!rtrack->getRecoHitInformation(hit)) continue;
      UChar_t  lay = hit->getICLayer();
      UShort_t adc = hit->getADCCount();
      UShort_t tdc = hit->getTDCCount();
      unsigned short bid = geo.getBoardID(WireID(lay, hit->getIWire()));
      m_hADCBoard->Fill(bid, adc);
      m_hADCLayer->Fill(lay, adc);
      m_hTDC->Fill(bid, tdc);
      m_hPhiHit->Fill(phi, lay);
    }
  }
}

void CDCDQMSimModule::endRun()
{
  m_hNEvents->SetBinContent(1, m_nEvents);
}

double CDCDQMSimModule::getShiftedPhi(const ROOT::Math::XYZVector& pos, int lay)
{
  static CDCGeometryPar& geo = CDCGeometryPar::Instance();
  double phi = TMath::ATan2(pos.Y(), pos.X());
  if (m_adjustWireShift) {
    int nShifts = geo.nShifts(lay);
    if (nShifts) {
      int    nW   = geo.nWiresInLayer(lay);
      double fZ   = geo.senseWireFZ(lay);
      double bZ   = geo.senseWireBZ(lay);
      double R    = geo.senseWireR(lay);
      double phiF = (2 * TMath::Pi() / nW) * 0.5 * nShifts;
      B2Vector3D f(R * TMath::Cos(phiF), R * TMath::Sin(phiF), fZ);
      B2Vector3D b(R, 0, bZ);
      B2Vector3D u = (f - b).Unit();
      double beta = (pos.Z() - b.Z()) / u.Z();
      B2Vector3D p = b + beta * u;
      phi -= TMath::ATan2(p.Y(), p.X());
    }
  }
  while (phi <  0)             phi += 2 * TMath::Pi();
  while (phi >= 2 * TMath::Pi()) phi -= 2 * TMath::Pi();
  return phi;
}

int CDCDQMSimModule::findPhiBin(double phi, int lay)
{
  static CDCGeometryPar& geo = CDCGeometryPar::Instance();
  int    nW     = geo.nWiresInLayer(lay);
  double shift  = (geo.offset(lay) - 0.5) * 2 * TMath::Pi() / nW;
  double lo     = shift;
  double hi     = 2 * TMath::Pi() + shift;
  if (phi <  lo) phi += 2 * TMath::Pi();
  else if (phi >= hi) phi -= 2 * TMath::Pi();
  return static_cast<int>((phi - lo) / ((hi - lo) / nW)) + 1;
}
