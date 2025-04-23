import ROOT

def filter_gen_electrons(df, histograms):

    df = df.Define('genmask', 'abs(genpart_pdg) == 11 && abs(genpart_mompdg) == 9000007')
    df = df.Define('vx_e_genms', 'genpart_vx[genmask]')
    df = df.Define('vy_e_genms', 'genpart_vy[genmask]')
    df = df.Define('vz_e_genms', 'genpart_vz[genmask]')
    df = df.Define('vt_e_genms', 'sqrt(vx_e_genms*vx_e_genms + vy_e_genms*vy_e_genms)')
    df = df.Define('pt_e_genms', 'genpart_pt[genmask]')
    df = df.Define('eta_e_genms', 'genpart_eta[genmask]')
    df = df.Define('phi_e_genms', 'genpart_phi[genmask]')
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

    STR_getminangdiffs = 'getminangs(pt_e_genms, eta_e_genms, phi_e_genms, ele_pt, ele_eta, ele_phi)'
    df = df.Define('gen_el_deta', f'std::get<0>({STR_getminangdiffs})')
    df = df.Define('gen_el_dphi', f'std::get<1>({STR_getminangdiffs})')
    df = df.Define('gen_el_dR', f'std::get<2>({STR_getminangdiffs})')
    histograms.append(df.Histo1D(('gen_el_deta', 'gen_el_deta', 20000, -1, 1), 'gen_el_deta'))
    histograms.append(df.Histo1D(('gen_el_dphi', 'gen_el_dphi', 20000, -10, 10), 'gen_el_dphi'))
    histograms.append(df.Histo1D(('gen_el_dR', 'gen_el_dR', 10000, 0, 10), 'gen_el_dR'))

    return df

def analyse_mcfiles(infile, outfile):
    df = ROOT.RDataFrame('demo/tree', infile)
    print(f'Entries in the {infile}: {df.Count().GetValue()}')

    histograms = []
    df = filter_gen_electrons(df, histograms)
    df = do_gen_matching(df, histograms)

    outfile = ROOT.TFile(outfile, 'RECREATE')
    for hist in histograms:
        hist.Write()
    outfile.Close()


if __name__=='__main__':

    STRCPPFUNC_getminangs = """
    std::tuple< ROOT::VecOps::RVec<double>, ROOT::VecOps::RVec<double>,
        ROOT::VecOps::RVec<double> > getminangs(ROOT::VecOps::RVec<double> &gpt,
                                               ROOT::VecOps::RVec<double> &geta,
                                               ROOT::VecOps::RVec<double> &gphi,
                                               ROOT::VecOps::RVec<double> &pt,
                                               ROOT::VecOps::RVec<double> &eta,
                                               ROOT::VecOps::RVec<double> &phi) {

        ROOT::VecOps::RVec<double> mindeta;
        ROOT::VecOps::RVec<double> mindphi;
        ROOT::VecOps::RVec<double> mindR;

        for (int i = 0; i < gpt.size(); i++) {
            TVector3 gvec;
            gvec.SetPtEtaPhi(gpt[i], geta[i], gphi[i]);
            float min_dR = 99999;
            float min_deta = 99999;
            float min_dphi = 99999;
            for (int j = 0; j < pt.size(); j++) {
                TVector3 vec;
                vec.SetPtEtaPhi(pt[j], eta[j], phi[j]);
                float dR = gvec.DeltaR(vec);
                float deta = gvec.Eta() - vec.Eta();
                float dphi = gvec.DeltaPhi(vec);
                if (dR < min_dR) min_dR = dR;
                if (fabs(deta) < fabs(min_deta)) min_deta = deta;
                if (dphi < min_dphi) min_dphi = dphi;
            }
            mindeta.push_back(min_deta);
            mindphi.push_back(min_dphi);
            mindR.push_back(min_dR);
        }
        return std::make_tuple(mindeta, mindphi, mindR);
    }
    """
    ROOT.gInterpreter.Declare(STRCPPFUNC_getminangs)

    analyse_mcfiles('./data/NTuples_250418_3cm.root', './hists/hist_3cm.root')
    analyse_mcfiles('./data/NTuples_250418_30cm.root', './hists/hist_30cm.root')
    analyse_mcfiles('./data/NTuples_250418_1m.root', './hists/hist_1m.root')
    analyse_mcfiles('./data/NTuples_250418_3m.root', './hists/hist_3m.root')
