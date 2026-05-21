#include <cdc/modules/CDCWireEfficiencyCanvas/CDCWireEfficiencyCanvasModule.h>
#include <cdc/geometry/CDCGeometryPar.h>
#include <framework/logging/Logger.h>

#include <TFile.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>

#include <cmath>

using namespace Belle2;

REG_MODULE(CDCWireEfficiencyCanvas);

CDCWireEfficiencyCanvasModule::CDCWireEfficiencyCanvasModule() : Module()
{
  setDescription("Produces a 4-panel CDC wire-efficiency canvas from CDCDQMSim output.");
  addParam("InputFile", m_inputFile, "CDCDQMSim ROOT file containing CDC/hTrackingWireEff",
           std::string("../output/dqm_merged.root"));
  addParam("Output", m_output, "Output path prefix — .pdf and .root are appended",
           std::string("../output/wire_eff"));
  addParam("Label",    m_label,    "Main label drawn top-right of each panel",
           std::string("Belle II Simulation"));
  addParam("Sublabel", m_sublabel, "Second label line (decay mode / conditions)",
           std::string(""));
  addParam("Release",  m_release,  "Third label line (software release tag)",
           std::string(""));
}

void CDCWireEfficiencyCanvasModule::initialize()
{
  TFile* fin = TFile::Open(m_inputFile.c_str(), "READ");
  if (!fin || fin->IsZombie())
    B2FATAL("Cannot open input file: " << m_inputFile);

  m_hSrc = dynamic_cast<TH2F*>(fin->Get("CDC/hTrackingWireEff"));
  if (!m_hSrc)
    B2FATAL("CDC/hTrackingWireEff not found in " << m_inputFile);
  m_hSrc->SetDirectory(0);
  fin->Close();
}

// ── simulation label ──────────────────────────────────────────────────────────
static void drawLabel(TPad* p,
                      const std::string& label,
                      const std::string& sublabel,
                      const std::string& release,
                      double xshift = 0.0,
                      double yshift = 0.0)
{
  p->cd();
  TLatex tex;
  tex.SetNDC(kTRUE);
  tex.SetTextAlign(33);  // top-right
  tex.SetTextSize(0.036);
  const double dy = 0.050;
  double x = 0.84 + xshift, y = 0.92 + yshift;
  if (!label.empty()) {
    tex.SetTextFont(72);   // italic bold (Belle II house style)
    tex.DrawLatex(x, y, label.c_str());
    y -= dy;
  }
  tex.SetTextFont(42);
  if (!sublabel.empty()) {
    tex.DrawLatex(x, y, sublabel.c_str());
    y -= dy;
  }
  if (!release.empty()) {
    tex.DrawLatex(x, y, release.c_str());
  }
}

// ── polygon builder ────────────────────────────────────────────────────────────
TH2Poly* CDCWireEfficiencyCanvasModule::makePoly(const std::string& name,
                                                 const std::string& title) const
{
  const CDC::CDCGeometryPar& geo = CDC::CDCGeometryPar::Instance();
  const int nSL = geo.getNumberOfSenseLayers();
  const double maxR = geo.senseWireR(nSL - 1);

  auto* h = new TH2Poly(name.c_str(), title.c_str(),
                        -maxR * 1.05, maxR * 1.05,
                        -maxR * 1.05, maxR * 1.05);

  for (int lay = 0; lay < nSL; ++lay) {
    const int    nW  = geo.nWiresInLayer(lay);
    const double off = geo.offset(lay);
    const double R   = geo.senseWireR(lay);

    // inner/outer boundary radii — midpoints between adjacent sense layers
    double r1, r2;
    if (lay == 0) {
      const double dr = geo.senseWireR(1) - R;
      r1 = R - dr / 2.0;  r2 = R + dr / 2.0;
    } else if (lay == nSL - 1) {
      const double dr = R - geo.senseWireR(lay - 1);
      r1 = R - dr / 2.0;  r2 = R + dr / 2.0;
    } else {
      r1 = R - (R - geo.senseWireR(lay - 1)) / 2.0;
      r2 = R + (geo.senseWireR(lay + 1) - R)  / 2.0;
    }

    const double dPhi = 2.0 * M_PI / nW;

    for (int wire = 0; wire < nW; ++wire) {
      const double phi  = dPhi * (wire + off);
      const double phi1 = phi - dPhi * 0.5;   // left boundary
      const double phi2 = phi + dPhi * 0.5;   // right boundary

      // Simple formula (DQMHistAnalysisCDCEpics::createEffiTH2Poly):
      // corners lie exactly on the inner/outer circles at phi1 and phi2.
      // Adjacent bins share corner coordinates exactly → no gap at any phi
      // boundary regardless of layer offset (0 or 0.5).
      Double_t xx[] = {
        r1* std::cos(phi1),    // inner-left
        r2* std::cos(phi1),    // outer-left
        r2* std::cos(phi2),    // outer-right
        r1* std::cos(phi2)     // inner-right
      };
      Double_t yy[] = {
        r1* std::sin(phi1),    // inner-left
        r2* std::sin(phi1),    // outer-left
        r2* std::sin(phi2),    // outer-right
        r1* std::sin(phi2)     // inner-right
      };
      h->AddBin(4, xx, yy);
    }
  }
  return h;
}

// ── main event ────────────────────────────────────────────────────────────────
void CDCWireEfficiencyCanvasModule::event()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  const CDC::CDCGeometryPar& geo = CDC::CDCGeometryPar::Instance();
  const int nSL = geo.getNumberOfSenseLayers();

  // ── build empty polys ────────────────────────────────────────────────────
  B2INFO("CDCCanvas: building TH2Poly geometry ...");
  TH2Poly* hObs = makePoly("hist_attachedWires", ";X [cm];Y [cm]");
  TH2Poly* hExp = makePoly("hist_expectedWires",  ";X [cm];Y [cm]");
  TH2Poly* hEff = makePoly("hist_wireAttachEff",  ";X [cm];Y [cm]");
  auto*    h1d  = new TH1F("hist_wire_attach_eff_1d",
                           ";Efficiency;Wires / bin", 104, -0.02, 1.02);

  // ── fill from hTrackingWireEff ───────────────────────────────────────────
  B2INFO("CDCCanvas: filling histograms ...");
  for (int lay = 0; lay < nSL; ++lay) {
    const int    nW      = geo.nWiresInLayer(lay);
    const double off     = geo.offset(lay);
    const double R       = geo.senseWireR(lay);
    const double dPhi    = 2.0 * M_PI / nW;
    const int    expYBin = lay + 1;
    const int    obsYBin = lay + nSL + 1;

    for (int wire = 0; wire < nW; ++wire) {
      const double phi  = dPhi * (wire + off);
      const double x    = R * std::cos(phi);
      const double y    = R * std::sin(phi);
      const double expV = m_hSrc->GetBinContent(wire + 1, expYBin);
      const double obsV = m_hSrc->GetBinContent(wire + 1, obsYBin);

      hObs->Fill(x, y, obsV);
      hExp->Fill(x, y, expV);
      if (expV > 0) {
        const double eff = obsV / expV;
        hEff->Fill(x, y, eff);
        h1d->Fill(eff);
      }
    }
  }

  // ── canvas layout ────────────────────────────────────────────────────────
  auto* c = new TCanvas("canvas_wire_eff", "CDC Wire Tracking Efficiency", 1600, 1400);

  const double gap = 0.005;
  auto makePad = [&](const char* name, double x1, double y1, double x2, double y2) {
    auto* p = new TPad(name, name, x1, y1, x2, y2);
    p->Draw();
    return p;
  };

  TPad* p1 = makePad("p1", gap,       0.5 + gap, 0.5 - gap, 1.0 - gap);
  TPad* p2 = makePad("p2", 0.5 + gap, 0.5 + gap, 1.0 - gap, 1.0 - gap);
  TPad* p3 = makePad("p3", gap,       gap,        0.5 - gap, 0.5 - gap);
  TPad* p4 = makePad("p4", 0.5 + gap, gap,        1.0 - gap, 0.5 - gap);

  auto setupCirclePad = [](TPad * p) {
    // left+right = top+bottom → equal plot-area fraction → 1:1 axis scale
    p->SetLeftMargin(0.13);
    p->SetRightMargin(0.15);
    p->SetTopMargin(0.11);
    p->SetBottomMargin(0.07);
  };

  TLatex lbl;
  lbl.SetTextSize(0.060);
  lbl.SetTextFont(62);
  lbl.SetNDC(kTRUE);

  // (a) observed
  p1->cd();
  setupCirclePad(p1);
  hObs->SetStats(0);
  hObs->GetYaxis()->SetTitleOffset(1.2);
  hObs->GetZaxis()->SetRangeUser(0, hObs->GetMaximum());
  hObs->Draw("COLZ");
  lbl.DrawLatex(0.15, 0.82, "(a)");
  drawLabel(p1, m_label, m_sublabel, m_release);

  // (b) expected
  p2->cd();
  setupCirclePad(p2);
  hExp->SetStats(0);
  hExp->GetYaxis()->SetTitleOffset(1.2);
  hExp->GetZaxis()->SetRangeUser(0, hExp->GetMaximum());
  hExp->Draw("COLZ");
  lbl.DrawLatex(0.15, 0.82, "(b)");
  drawLabel(p2, m_label, m_sublabel, m_release);

  // (c) efficiency
  p3->cd();
  setupCirclePad(p3);
  p3->SetTopMargin(0.03);
  p3->SetBottomMargin(0.15);
  hEff->SetStats(0);
  hEff->GetYaxis()->SetTitleOffset(1.2);
  hEff->GetZaxis()->SetTitle("Efficiency");
  hEff->GetZaxis()->SetTitleOffset(1.3);
  hEff->GetZaxis()->SetRangeUser(0, 1);
  hEff->Draw("COLZ");
  lbl.DrawLatex(0.15, 0.87, "(c)");
  drawLabel(p3, m_label, m_sublabel, m_release, 0.0, 0.08);

  // ── compute efficiency category fractions from h1d ──────────────────────
  const double total  = h1d->Integral(1, h1d->GetNbinsX());
  const int    b08    = h1d->FindBin(0.08);
  const int    b72    = h1d->FindBin(0.72);
  const double fDead  = (total > 0) ? h1d->Integral(1,    b08 - 1)              / total * 100.0 : 0.0;
  const double fMid   = (total > 0) ? h1d->Integral(b08,  b72 - 1)              / total * 100.0 : 0.0;
  const double fGood  = (total > 0) ? h1d->Integral(b72,  h1d->GetNbinsX())     / total * 100.0 : 0.0;
  const double meanEff = h1d->GetMean();

  char sDead[64], sMid[64], sGood[64], sMean[64];
  std::snprintf(sDead, sizeof(sDead), "%.2f%% : #varepsilon < 0.08",           fDead);
  std::snprintf(sMid,  sizeof(sMid),  "%.2f%% : 0.08 #leq #varepsilon < 0.72", fMid);
  std::snprintf(sGood, sizeof(sGood), "%.2f%% : #varepsilon #geq 0.72",        fGood);
  std::snprintf(sMean, sizeof(sMean), "Mean #varepsilon = %.3f",               meanEff);

  // (d) 1D distribution
  p4->cd();
  p4->SetLeftMargin(0.14);
  p4->SetRightMargin(0.06);
  p4->SetTopMargin(0.03);
  p4->SetBottomMargin(0.15);
  h1d->SetStats(0);
  h1d->SetLineWidth(2);
  h1d->Draw();
  lbl.DrawLatex(0.16, 0.87, "(d)");
  drawLabel(p4, m_label, m_sublabel, m_release, -0.08);

  // draw stats inside panel (d)
  // percentages right-aligned at sx_r; descriptions left-aligned at sx_l → colons align
  char sPct1[32], sPct2[32], sPct3[32];
  std::snprintf(sPct1, sizeof(sPct1), "%.2f%%", fGood);
  std::snprintf(sPct2, sizeof(sPct2), "%.2f%%", fMid);
  std::snprintf(sPct3, sizeof(sPct3), "%.2f%%", fDead);

  TLatex stats;
  stats.SetNDC(kTRUE);
  stats.SetTextSize(0.036);
  stats.SetTextFont(42);

  const double sx_def_r = 0.45;   // right edge of definition column
  const double sx_pct_l = 0.465; // left edge of ": XX.XX%" column
  double sy = 0.58;
  const double sdy = 0.045;

  // mean first, left-aligned
  stats.SetTextAlign(11);
  stats.DrawLatex(0.20, sy, sMean);
  sy -= sdy;

  // three category lines: definition right-aligned | ": XX.XX%" left-aligned
  stats.SetTextAlign(31);
  stats.DrawLatex(sx_def_r, sy,          "#varepsilon #geq 0.72");
  stats.DrawLatex(sx_def_r, sy - sdy,    "0.08 #leq #varepsilon < 0.72");
  stats.DrawLatex(sx_def_r, sy - 2 * sdy,  "#varepsilon < 0.08");

  char sColon1[32], sColon2[32], sColon3[32];
  std::snprintf(sColon1, sizeof(sColon1), ": %.2f%%", fGood);
  std::snprintf(sColon2, sizeof(sColon2), ": %.2f%%", fMid);
  std::snprintf(sColon3, sizeof(sColon3), ": %.2f%%", fDead);

  stats.SetTextAlign(11);
  stats.DrawLatex(sx_pct_l, sy,          sColon1);
  stats.DrawLatex(sx_pct_l, sy - sdy,    sColon2);
  stats.DrawLatex(sx_pct_l, sy - 2 * sdy,  sColon3);

  c->cd();
  c->SaveAs((m_output + ".pdf").c_str());
  B2INFO("CDCCanvas: saved " << m_output << ".pdf");

  // ── write ROOT file ──────────────────────────────────────────────────────
  TFile fout((m_output + ".root").c_str(), "RECREATE");
  hObs->Write();
  hExp->Write();
  hEff->Write();
  h1d->Write();
  fout.Close();
  B2INFO("CDCCanvas: saved " << m_output << ".root");

  delete hObs; delete hExp; delete hEff; delete h1d; delete c;
}
