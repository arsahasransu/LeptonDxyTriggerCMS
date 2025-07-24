import sys

import ROOT


base_variables = {'e': (100, 0, 200),
                  'pt': (100, 0, 100),
                  'eta': (54, -2.7, 2.7),
                  'phi': (66, -3.3, 3.3),
                  'IDloose': (4, -1, 3),
                  'IDtight': (4, -1, 3),
                  'seedtime': (300, -10, 20)}


def add_plots(df, hists, prefix, filter=''):

    for var, (xbins, xlow, xup) in base_variables.items():
        hists.append(df.Histo1D((f'{prefix}{filter}_{var}', f'{var}', xbins, xlow, xup), f'{prefix}_{var}'))


STRCPPFUNC_findelprobe_inZwindow = '''
    int findelprobe_inZwindow(ROOT::VecOps::RVec<double> &pt,
                                                  ROOT::VecOps::RVec<double> &eta,
                                                  ROOT::VecOps::RVec<double> &phi,
                                                  ROOT::VecOps::RVec<double> &e,
                                                  const double &pt0,
                                                  const double &eta0,
                                                  const double &phi0,
                                                  const double &e0,
                                                  const double Zlow,
                                                  const double Zhigh) {

        int elprobeidx = -1;

        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> el0(pt0, eta0, phi0, e0);
        for (int i = 0; i < pt.size(); i++) {
            if(pt[i] == pt0 && eta[i] == eta0 && phi[i] == phi0 && e[i] == e0) continue;
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> el1(pt[i], eta[i], phi[i], e[i]);
            double mass = (el0 + el1).M();
            if (mass > Zlow && mass < Zhigh) {
                elprobeidx = i;
                break;
            }
        }

        return elprobeidx;
    }
'''
ROOT.gInterpreter.Declare(STRCPPFUNC_findelprobe_inZwindow)


STRCPPFUNC_get_inv_mass = '''
    double get_inv_mass(const double e0, const double pt0, const double eta0, const double phi0,
                        const double e1, const double pt1, const double eta1, const double phi1) {
        double mass = -1.0;
        if (e0 > 0 && e1 > 0) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> el0(pt0, eta0, phi0, e0);
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> el1(pt1, eta1, phi1, e1);
            mass = (el0 + el1).M();
        }
        return mass;
    }
'''
ROOT.gInterpreter.Declare(STRCPPFUNC_get_inv_mass)

def analyser(df, outfilename):

    histograms = []

    add_plots(df, histograms, 'ele')

    # find one very good electron
    # assume pt sorted electron collection
    tag_electron_selection = 'ele_pt > 20 && abs(ele_eta) < 2.1 && ele_IDtight == 1 && ele_seedtime > 1.0'
    for var in base_variables.keys():
        df = df.Define(f'eletag_{var}', f'ele_{var}[{tag_electron_selection}]')
    df = df.Filter('eletag_e.size() > 0')
    add_plots(df, histograms, 'eletag')

    for var in base_variables.keys():
        df = df.Define(f'eletag0_{var}', f'eletag_{var}[0]')
    add_plots(df, histograms, 'eletag0')

    # find other electron in Z window
    Zlow = 84
    Zhigh = 96
    df = df.Define('elprobe_inZwindow_idx', f'findelprobe_inZwindow(ele_pt, ele_eta, ele_phi, ele_e, eletag0_pt, eletag0_eta, eletag0_phi, eletag0_e, {Zlow}, {Zhigh})')
    df = df.Filter('elprobe_inZwindow_idx >= 0')
    for var in base_variables.keys():
        df = df.Define(f'elprobe_inZwindow_{var}', f'ele_{var}[elprobe_inZwindow_idx]')
    add_plots(df, histograms, 'elprobe_inZwindow')

    # calculate the invariant mass of the two electrons
    df = df.Define('invmass', f'get_inv_mass(eletag0_e, eletag0_pt, eletag0_eta, eletag0_phi, elprobe_inZwindow_e, elprobe_inZwindow_pt, elprobe_inZwindow_eta, elprobe_inZwindow_phi)')
    histograms.append(df.Histo1D(('invmass', 'invM', 150, 0, 150), 'invmass'))

    # For seeding time trigger efficiency
    dftightid = df.Filter('elprobe_inZwindow_IDtight == 1')
    for var in base_variables.keys():
        dftightid = dftightid.Define(f'tid_elprobe_inZwindow_{var}', f'elprobe_inZwindow_{var}')
    add_plots(dftightid, histograms, 'tid_elprobe_inZwindow')

    dftidhlt = dftightid.Filter('HLT_DiPhoton10Time1ns == 1')
    for var in base_variables.keys():
        dftidhlt = dftidhlt.Define(f'tidhlt_elprobe_inZwindow_{var}', f'tid_elprobe_inZwindow_{var}')
    add_plots(dftidhlt, histograms, 'tidhlt_elprobe_inZwindow')

    # For trigger efficiency of other kinematics
    dftightid1ns = dftightid.Filter('tid_elprobe_inZwindow_seedtime > 1.0')
    for var in base_variables.keys():
        dftightid1ns = dftightid1ns.Define(f'tid1ns_elprobe_inZwindow_{var}', f'tid_elprobe_inZwindow_{var}')
    add_plots(dftightid1ns, histograms, 'tid1ns_elprobe_inZwindow')

    dftid1nshlt = dftightid1ns.Filter('HLT_DiPhoton10Time1ns == 1')
    for var in base_variables.keys():
        dftid1nshlt = dftid1nshlt.Define(f'tid1nshlt_elprobe_inZwindow_{var}', f'tid1ns_elprobe_inZwindow_{var}')
    add_plots(dftid1nshlt, histograms, 'tid1nshlt_elprobe_inZwindow')

    # Separate probe to leading and subleading
    dftid1ns_leadprobe = dftightid1ns.Filter('tid1ns_elprobe_inZwindow_pt > eletag0_pt')
    add_plots(dftid1ns_leadprobe, histograms, 'tid1ns_elprobe_inZwindow', 'lead')

    dftid1nshlt_leadprobe = dftid1ns_leadprobe.Filter('HLT_DiPhoton10Time1ns == 1')
    for var in base_variables.keys():
        dftid1nshlt_leadprobe = dftid1nshlt_leadprobe.Define(f'tid1nshlt_leadelprobe_inZwindow_{var}', f'tid1ns_elprobe_inZwindow_{var}')
    add_plots(dftid1nshlt_leadprobe, histograms, 'tid1nshlt_leadelprobe_inZwindow')

    outfile = ROOT.TFile(outfilename, 'RECREATE')
    for hist in histograms:
        hist.Write()
    outfile.Close()


if __name__ == "__main__":
    
    ROOT.EnableImplicitMT()
    # Get input and output file names
    infilename = sys.argv[1]
    outfilename = sys.argv[2]

    print(f'Starting analysis on {infilename} and saving to {outfilename}')

    df = ROOT.RDataFrame('demo/tree', infilename)
    print('Entries in the tree to process:', df.Count().GetValue())
    analyser(df, outfilename)