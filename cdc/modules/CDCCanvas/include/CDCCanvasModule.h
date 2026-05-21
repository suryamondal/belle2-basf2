#pragma once

#include <framework/core/Module.h>
#include <string>

class TH2Poly;
class TH2F;

namespace Belle2 {

  /**
   * Produces the 4-panel CDC wire-efficiency canvas from a CDCDQMSim output file.
   *
   * Panel layout:
   *   (a) Observed hits        — TH2Poly backplate view
   *   (b) Expected crossings   — TH2Poly backplate view
   *   (c) Wire efficiency      — TH2Poly backplate view
   *   (d) Efficiency 1D dist.  — TH1F
   *
   * TH2Poly bins use the tangent formula from DQMHistAnalysisCDCMonObj so that
   * adjacent bins share their inner boundary exactly, with no triangular gaps.
   *
   * Usage (in 03_canvas.py):
   *   main.add_module('CDCCanvas',
   *                   InputFile='../output/dqm_merged.root',
   *                   Output='../output/wire_eff')
   */
  class CDCCanvasModule : public Module {

  public:
    CDCCanvasModule();
    ~CDCCanvasModule() = default;

    void initialize() override; /**< open input file and read hTrackingWireEff */
    void event() override;      /**< build TH2Poly histograms and write the canvas */

  private:
    /** Build a blank TH2Poly with one bin per CDC wire (tangent inner corners). */
    TH2Poly* makePoly(const std::string& name, const std::string& title) const;

    std::string m_inputFile; /**< Path to CDCDQMSim output ROOT file. */
    std::string m_output;    /**< Output path prefix (.pdf and .root appended). */
    std::string m_label;     /**< Main label line, e.g. "Belle II Simulation". */
    std::string m_sublabel;  /**< Second label line, e.g. decay and conditions. */
    std::string m_release;   /**< Third label line, e.g. software release tag. */

    TH2F* m_hSrc{nullptr};  /**< hTrackingWireEff read from input file. */
  };

} // namespace Belle2
