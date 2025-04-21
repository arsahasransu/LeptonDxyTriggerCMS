import ROOT
from ROOT import TH1D, TCanvas

file = ROOT.TFile('out_histos.root')

def plot_seedtime():

    data_prompt = file.Get('eleEB_el0_seedtime')
    data_all = file.Get('eleEB_el0_seedtime')

    data_prompt.Scale(1.0 / data_prompt.Integral())
    data_all.Scale(1.0 / data_all.Integral())

    c = TCanvas('c', 'c', 800, 600)
    data_prompt.SetLineColor(ROOT.kRed)
    data_prompt.SetLineWidth(2)
    data_all.SetLineColor(ROOT.kBlue)
    data_all.SetLineWidth(2)
    data_prompt.Draw('hist')
    # data_all.Draw('hist same')
    c.SaveAs('seedtime.png')


if __name__ == '__main__':
    plot_seedtime()