import ternary
import numpy as np
from scipy.interpolate import interp1d
import math
import seaborn as sns


def color_point(x, y, z, scale):
    w = 255
    x_color = x * w / float(scale)
    y_color = y * w / float(scale)
    z_color = z * w / float(scale)
    r = math.fabs(w - y_color) / w
    g = math.fabs(w - x_color) / w
    b = math.fabs(w - z_color) / w
    return (r, g, b, .3)


def generate_heatmap_data(scale=5):
    from ternary.helpers import simplex_iterator
    d = dict()
    for (i, j, k) in simplex_iterator(scale):
        d[(i, j, k)] = color_point(i, j, k, scale)
    return d


def get_drug_trajectories(dfi, cell_line, drugs):
    drug_trajectories = {}
    drug_scatter = {}
    drug_conc = {}
    dfs = dfi[dfi.cell_line == cell_line].copy()
    dfs = recompute_fractions(dfs)
    for drug in drugs:
        print(drug)
        dfsd = dfs[dfs.agent == drug].copy()
        dfsd = dfsd.sort_values('concentration')
        dfc = dfsd[['g2', 'S', 'g1']].copy()
        points = dfc.values
        drug_trajectories[drug] = smooth_points(points)
        drug_scatter[drug] = points
        drug_conc[drug] = [20 * d for d in
                           dfsd.concentration.tolist()]
    return drug_trajectories, drug_scatter, drug_conc


def plot_trajectories(dfi, cell_line, drugs):
    dfi['agent'] = dfi['agent'].fillna('control')
    colrs = sns.color_palette("hls", len(drugs))
    colrs[np.argmax(drugs == 'control')] = (0, 0, 0)
    colr_dict = {d: c for d, c in zip(drugs, colrs)}
    drug_trajectories, drug_scatter, drug_conc = get_drug_trajectories(
        dfi, cell_line, drugs)
    # data = generate_heatmap_data(scale=10)
    fig, tax = ternary.figure(scale=10)
    # tax.heatmap(data, colormap=False)
    tax.boundary()
    tax.gridlines(multiple=1, color='black')
    for drug in list(drug_trajectories.keys()):
        amin = np.min(drug_conc[drug])
        amax = np.max(drug_conc[drug])
        facecolor = np.zeros((len(drug_conc[drug]), 4))
        facecolor[:, :3] = colr_dict[drug]
        if drug != 'control':
            fc = [1 - (ds-amin)/amax for ds in drug_conc[drug]]
            facecolor[:, 3] = fc 
        else:
            drug_conc[drug] = [1] * len(drug_conc[drug])
            fc = [0.5] * len(drug_conc[drug])
            facecolor[:, 3] = fc 
                                
        #tax.plot(10 * drug_trajectories[drug], linewidth=0.5,
        #          label=drug, linestyle='--', color=colr_dict[drug])
        tax.scatter(10 * drug_scatter[drug], label=drug,
                    facecolor=facecolor, edgecolor=colr_dict[drug],
                    s=drug_conc[drug])
    tax.legend()
    tax.top_corner_label("S", fontsize=20)
    tax.left_corner_label("G1", fontsize=20)
    tax.right_corner_label("G2", fontsize=20)
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    tax.show()
    return colr_dict
    

def recompute_fractions(dfs):
    df = dfs.copy()
    df['g1'] = [g1+0.5 * sd for g1, sd
                in zip(df['G1'], df['S_dropout'])]
    df['g2'] = [g2+0.5 * sd for g2, sd
                in zip(df['G2'], df['S_dropout'])]
    df['live_frac'] = df[['g1', 'g2', 'S']].sum(axis=1)
    df[['g1', 'g2', 'S']] = df[['g1', 'g2', 'S']].div(df['live_frac'],
                                                      axis=0)
    return df


def smooth_points(points):
    tpoints = points.T
    interp_points = []
    for ftr in tpoints:
        x = np.linspace(0, 1, len(points))
        f = interp1d(x, ftr, kind='cubic')
        xnew = np.linspace(0, 1, 100)
        interp_ftr = f(xnew)
        interp_points.append(interp_ftr)
    interp_points = np.array(interp_points).T
    return interp_points
