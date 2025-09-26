import ROOT


calculate_inv_mass_str = '''
    #include "Math/Vector4D.h"

    double get_inv_mass(ROOT::VecOps::RVec<double> energy,
                        ROOT::VecOps::RVec<double> pt,
                        ROOT::VecOps::RVec<double> eta,
                        ROOT::VecOps::RVec<double> phi) {

        double mass = -1.0;

        if(energy.size()>=2) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> el0(pt[0], eta[0], phi[0], energy[0]);
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> el1(pt[1], eta[1], phi[1], energy[1]);
            mass = (el0+el1).M();
        }
        return mass;
    }
'''
ROOT.gInterpreter.Declare(calculate_inv_mass_str)


def add_new_variables(df, selId):

    df = df.Define(f'ele{selId}_invmass', f'get_inv_mass(ele{selId}_e, ele{selId}_pt, ele{selId}_eta, ele{selId}_phi)')

    return df

def add_plots(df, prevSelId, selId, selStr):

    histograms = []

    base_variables = {'e': (100, 0, 200),
                      'pt': (100, 0, 100),
                      'eta': (54, -2.7, 2.7),
                      'phi': (66, -3.3, 3.3),
                      'IDtight': (4, -1, 3),
                      'seedtime': (300, -10, 20)}
    
    new_variables = {'invmass':(150, 0, 150)}
    
    for var, (xbins, xlow, xup) in base_variables.items():
        df = df.Define(f'ele{selId}_{var}', f'ele{prevSelId}_{var}[{selStr}]')

    df = df.Define(f'ele{selId}_n', f'ele{selId}_e.size()')
    histograms.append(df.Histo1D((f'ele{selId}_n', 'multiplicity', 10, 0, 10), f'ele{selId}_n'))
    df = add_new_variables(df, selId)

    dfn2 = df.Filter(f'ele{selId}_n >= 2')
    for var, (xbins, xlow, xup) in base_variables.items():
        dfn2 = dfn2.Define(f'ele{selId}_el0_{var}', f'ele{selId}_{var}[0]')
        histograms.append(dfn2.Histo1D((f'ele{selId}_el0_{var}', f'{var}', xbins, xlow, xup), f'ele{selId}_el0_{var}'))

        dfn2 = dfn2.Define(f'ele{selId}_el1_{var}', f'ele{selId}_{var}[1]')
        histograms.append(dfn2.Histo1D((f'ele{selId}_el1_{var}', f'{var}', xbins, xlow, xup), f'ele{selId}_el1_{var}'))

    # dfn2 = add_new_variables(dfn2, selId)

    dfn2_invmass_filtered = dfn2.Filter(f'ele{selId}_invmass > 0')
    histograms.append(dfn2_invmass_filtered.Histo1D((f'ele{selId}_invmass', 'invM', new_variables['invmass'][0],
                        new_variables['invmass'][1], new_variables['invmass'][2]), f'ele{selId}_invmass'))

    return (df, histograms)

def analyser(df):

    histograms = []

    (df_eb, hist_list) = add_plots(df, '', 'EB', 'abs(ele_eta)<1.2')
    histograms.extend(hist_list)
    (df_eb, hist_list) = add_plots(df_eb, 'EB', 'EBID', 'abs(eleEB_eta)<1.2 && eleEB_IDtight == 1')
    histograms.extend(hist_list)
    df_eb = df_eb.Filter('eleEBID_invmass > 84 and eleEBID_invmass < 96')
    (df_eb, hist_list) = add_plots(df_eb, 'EBID', 'EBZ', 'abs(eleEBID_eta)<1.2 && eleEBID_IDtight == 1')
    histograms.extend(hist_list)

    (df_ee, hist_list) = add_plots(df, '', 'EE', 'abs(ele_eta)>1.6 && abs(ele_eta)<2.1 ')
    histograms.extend(hist_list)
    (df_ee, hist_list) = add_plots(df_ee, 'EE', 'EEID', 'abs(eleEE_eta)>1.6 && abs(eleEE_eta)<2.1 && eleEE_IDtight == 1')
    histograms.extend(hist_list)
    df_ee = df_ee.Filter('eleEEID_invmass > 84 and eleEEID_invmass < 96')
    (df_ee, hist_list) = add_plots(df_ee, 'EEID', 'EEZ', 'abs(eleEEID_eta)>1.6 && abs(eleEEID_eta)<2.1 && eleEEID_IDtight == 1')
    histograms.extend(hist_list)

    outfile = ROOT.TFile('data_histos_newseedtimeplot.root', 'RECREATE')
    for hist in histograms:
        hist.Write()
    outfile.Close()



if __name__ == "__main__":
    print('Starting analysis...')
    # df = ROOT.RDataFrame('demo/tree', './data/EGamma0_EXOLLPTRG_Nano.root')
    df = ROOT.RDataFrame('demo/tree', './data/DYTo2L_Run3Winter25_EXOLLPTRG_250923Nano.root')
    print('Entries in the tree to process:', df.Count().GetValue())
    analyser(df)
