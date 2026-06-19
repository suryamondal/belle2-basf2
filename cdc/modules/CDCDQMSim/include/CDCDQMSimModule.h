#pragma once

#include <framework/core/HistoModule.h>
#include <framework/datastore/StoreArray.h>
#include <framework/datastore/StoreObjPtr.h>

#include <cdc/dataobjects/CDCHit.h>
#include <mdst/dataobjects/Track.h>
#include <mdst/dataobjects/TrackFitResult.h>
#include <mdst/dataobjects/HitPatternVXD.h>
#include <mdst/dataobjects/HitPatternCDC.h>
#include <tracking/dataobjects/RecoTrack.h>

#include <TH1F.h>
#include <TH2F.h>

namespace Belle2 {

  /**
   * CDC DQM module for simulation — identical to CDCDQMModule but with all
   * TRGSummary and SoftwareTriggerResult checks removed so it works without
   * trigger or HLT simulation.
   */
  class CDCDQMSimModule : public HistoModule {

  public:
    CDCDQMSimModule();
    ~CDCDQMSimModule() = default;

    void defineHisto() override; /**< book DQM histograms */
    void initialize() override;  /**< register StoreArrays and set up cuts */
    void beginRun() override;    /**< reset histograms at run start */
    void event() override;       /**< fill histograms per event */
    void endRun() override;      /**< normalise histograms at run end */
    void terminate() override {} /**< no-op */

  private:
    /** Compute azimuthal angle of the wire at the given XY position in a given layer. */
    double getShiftedPhi(const ROOT::Math::XYZVector& position, int lay);
    /** Return the phi histogram bin index for the given phi and layer. */
    int    findPhiBin(double phi, int lay);

    StoreArray<CDCHit>    m_cdcHits;    /**< CDC hits */
    StoreArray<Track>     m_Tracks;     /**< reconstructed tracks */
    StoreArray<RecoTrack> m_RecoTracks; /**< reconstructed track objects */

    int    m_minHits{0};           /**< minimum CDC hits per track */
    int    m_minNdf{20};           /**< minimum NDF for track fit quality */
    double m_minPt{0.0};           /**< minimum transverse momentum [GeV/c] */
    bool   m_adjustWireShift{true}; /**< apply per-layer wire offset correction */
    Long64_t m_nEvents{0};          /**< processed event counter */

    TH1F* m_hNEvents        {nullptr}; /**< number of events */
    TH1F* m_hOcc            {nullptr}; /**< CDC occupancy per event */
    TH2F* m_hHit            {nullptr}; /**< hit map (layer vs wire) */
    TH1F* m_hPhi            {nullptr}; /**< track phi distribution */
    TH2F* m_hPhiEff         {nullptr}; /**< phi efficiency (layer vs phi bin) */
    TH2F* m_hPhiHit         {nullptr}; /**< phi hit count (layer vs phi bin) */
    TH2F* m_hPhiNCDC        {nullptr}; /**< number of CDC hits per track vs phi */
    TH2F* m_hADCBoard       {nullptr}; /**< ADC count vs board number */
    TH2F* m_hADCLayer       {nullptr}; /**< ADC count vs layer */
    TH2F* m_hTDC            {nullptr}; /**< TDC count vs layer */
    TH2F* m_hTrackingWireEff{nullptr}; /**< wire-by-wire tracking efficiency (expected/observed) */
  };

} // namespace Belle2
