import ROOT

import mycpputils

def plot_some_vars_beforeselection(df, histograms):
    histograms.append(df.Histo1D(('ele_seedtime', 'ele_seedtime', 1000, -25, 25), 'ele_seedtime'))
    histograms.append(df.Histo1D(('pho_seedtime', 'pho_seedtime', 1000, -25, 25), 'pho_seedtime'))

    return df


def filter_gen_electrons(df, histograms):

    df = df.Define('genmask', 'abs(genpart_pdg) == 11 && abs(genpart_mompdg) == 9000007')
    df = df.Define('vx_e_genms', 'genpart_vx[genmask]')
    df = df.Define('vy_e_genms', 'genpart_vy[genmask]')
    df = df.Define('vz_e_genms', 'genpart_vz[genmask]')
    df = df.Define('vt_e_genms', 'sqrt(vx_e_genms*vx_e_genms + vy_e_genms*vy_e_genms)')
    df = df.Define('pt_e_genms', 'genpart_pt[genmask]')
    df = df.Define('eta_e_genms', 'genpart_eta[genmask]')
    df = df.Define('phi_e_genms', 'genpart_phi[genmask]')
    df = df.Define('ecaleta_e_genms', 'genpart_ecaleta[genmask]')
    df = df.Define('ecalphi_e_genms', 'genpart_ecalphi[genmask]')
    df = df.Define('n_e_genms', 'pt_e_genms.size()')
    df = df.Filter('n_e_genms > 0')

    h_vt = df.Histo1D(('vt_e', 'vt', 10000, 0, 100), 'vt_e_genms')
    histograms.append(h_vt)
    h_pt = df.Histo1D(('pt_e', 'pt', 100, 0, 100), 'pt_e_genms')
    histograms.append(h_pt)
    h_eta = df.Histo1D(('eta_e', 'eta', 54, -2.7, 2.7), 'eta_e_genms')
    histograms.append(h_eta)
    h_phi = df.Histo1D(('phi_e', 'phi', 66, -3.3, 3.3), 'phi_e_genms')
    histograms.append(h_phi)
    h_n = df.Histo1D(('n_e', 'multiplicity', 10, 0, 10), 'n_e_genms')
    histograms.append(h_n)

    return df


def do_gen_matching(df, histograms):

    STR_getminangdiffs_el = 'getminangs(eta_e_genms, phi_e_genms, ele_eta, ele_phi)'
    df = df.Define('gen_el_deta', f'std::get<0>({STR_getminangdiffs_el})')
    df = df.Define('gen_el_dphi', f'std::get<1>({STR_getminangdiffs_el})')
    df = df.Define('gen_el_dR', f'std::get<2>({STR_getminangdiffs_el})')
    histograms.append(df.Histo1D(('gen_el_deta', 'gen_el_deta', 20000, -1, 1), 'gen_el_deta'))
    histograms.append(df.Histo1D(('gen_el_dphi', 'gen_el_dphi', 20000, -10, 10), 'gen_el_dphi'))
    histograms.append(df.Histo1D(('gen_el_dR', 'gen_el_dR', 10000, 0, 10), 'gen_el_dR'))

    STR_getminangdiffs_ph = 'getminangs(eta_e_genms, phi_e_genms, pho_eta, pho_phi)'
    df = df.Define('gen_ph_deta', f'std::get<0>({STR_getminangdiffs_ph})')
    df = df.Define('gen_ph_dphi', f'std::get<1>({STR_getminangdiffs_ph})')
    df = df.Define('gen_ph_dR', f'std::get<2>({STR_getminangdiffs_ph})')
    histograms.append(df.Histo1D(('gen_ph_deta', 'gen_ph_deta', 20000, -1, 1), 'gen_ph_deta'))
    histograms.append(df.Histo1D(('gen_ph_dphi', 'gen_ph_dphi', 20000, -10, 10), 'gen_ph_dphi'))
    histograms.append(df.Histo1D(('gen_ph_dR', 'gen_ph_dR', 10000, 0, 10), 'gen_ph_dR'))

    STR_getminangdiffs_ecalph = 'getminangs(ecaleta_e_genms, ecalphi_e_genms, pho_eta, pho_phi)'
    df = df.Define('gen_ecalph_deta', f'std::get<0>({STR_getminangdiffs_ecalph})')
    df = df.Define('gen_ecalph_dphi', f'std::get<1>({STR_getminangdiffs_ecalph})')
    df = df.Define('gen_ecalph_dR', f'std::get<2>({STR_getminangdiffs_ecalph})')
    histograms.append(df.Histo1D(('gen_ecalph_deta', 'gen_ecalph_deta', 20000, -1, 1), 'gen_ecalph_deta'))
    histograms.append(df.Histo1D(('gen_ecalph_dphi', 'gen_ecalph_dphi', 20000, -10, 10), 'gen_ecalph_dphi'))
    histograms.append(df.Histo1D(('gen_ecalph_dR', 'gen_ecalph_dR', 10000, 0, 10), 'gen_ecalph_dR'))

    STR_getmatchedidxs_el = 'getmatchedidxs(eta_e_genms, phi_e_genms, ele_eta, ele_phi, 0.05)'
    df = df.Define('genmatched_el_idx', f'std::get<1>({STR_getmatchedidxs_el})')
    df = df.Define('genmatched_el_seedtime', 'ele_seedtime[genmatched_el_idx != -1]')
    histograms.append(df.Histo1D(('genmatched_el_seedtime', 'genmatched_el_seedtime', 1000, -25, 25), 'genmatched_el_seedtime'))

    STR_getmatchedidxs_ph = 'getmatchedidxs(ecaleta_e_genms, ecalphi_e_genms, pho_eta, pho_phi, 0.17)'
    df = df.Define('genmatched_ph_idx', f'std::get<1>({STR_getmatchedidxs_ph})')
    df = df.Define('genmatched_ph_seedtime_eb', 'pho_seedtime[genmatched_ph_idx != -1 and abs(pho_eta) < 1.479]')
    df = df.Define('genmatched_ph_seedtime_ee', 'pho_seedtime[genmatched_ph_idx != -1 and abs(pho_eta) > 1.479]')
    histograms.append(df.Histo1D(('genmatched_ph_seedtime_eb', 'genmatched_ph_seedtime_eb', 1000, -25, 25), 'genmatched_ph_seedtime_eb'))
    histograms.append(df.Histo1D(('genmatched_ph_seedtime_ee', 'genmatched_ph_seedtime_ee', 1000, -25, 25), 'genmatched_ph_seedtime_ee'))

    return df


def analyse_mcfiles(infile, outfile):
    df = ROOT.RDataFrame('demo/tree', infile)
    print(f'Entries in the {infile}: {df.Count().GetValue()}')

    histograms = []
    df = plot_some_vars_beforeselection(df, histograms)
    df = filter_gen_electrons(df, histograms)
    df = do_gen_matching(df, histograms)

    outfile = ROOT.TFile(outfile, 'RECREATE')
    for hist in histograms:
        hist.Write()
    outfile.Close()


if __name__=='__main__':

    ROOT.gInterpreter.Declare(mycpputils.STRCPPFUNC_getminangs)
    ROOT.gInterpreter.Declare(mycpputils.STRCPPFUNC_getmatchedidxs)

    # analyse_mcfiles('./data/NTuples_250418_3cm.root', './hists/hist_3cm.root')
    # analyse_mcfiles('./data/NTuples_250418_30cm.root', './hists/hist_30cm.root')
    # analyse_mcfiles('./data/NTuples_250418_1m.root', './hists/hist_1m.root')
    # analyse_mcfiles('./data/NTuples_250418_3m.root', './hists/hist_3m.root')

    analyse_mcfiles('./data/NTuples_250419_3cm.root', './hists/hist_3cm.root')
    analyse_mcfiles('./data/NTuples_250419_30cm.root', './hists/hist_30cm.root')
    analyse_mcfiles('./data/NTuples_250419_1m.root', './hists/hist_1m.root')
    analyse_mcfiles('./data/NTuples_250419_3m.root', './hists/hist_3m.root')
