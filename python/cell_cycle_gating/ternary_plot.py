import ternary
import numpy as np
from scipy.interpolate import interp1d


def get_drug_trajectories(dfi, cell_line, drugs):
    drug_trajectories = {}
    drug_scatter = {}
    dfs = dfi[dfi.cell_line == cell_line].copy()
    dfs = recompute_fractions(dfs)
    for drug in drugs:
        print(drug)
        dfsd = dfs[dfs.agent == drug].copy()
        dfsd = dfsd.sort_values('concentration')
        dfc = dfsd[['g2', 'S', 'g1']].copy()
        points = dfc.values
        drug_trajectories[drug] = smooth_points(points)
        drug_scatter['drug'] = points
    return drug_trajectories, drug_scatter


def plot_trajectories(dfi, cell_line, drugs):
    drug_trajectories, drug_scatter = get_drug_trajectories(
        dfi, cell_line, drugs)
    fig, tax = ternary.figure(scale=1.0)
    tax.boundary()
    tax.gridlines(multiple=0.1, color='black')
    for drug in list(drug_trajectories.keys()):
        tax.plot(drug_trajectories[drug], linewidth=1.0,
                 label=drug, linestyle='-')
    tax.legend()
    tax.top_corner_label("S", fontsize=20)
    tax.left_corner_label("G1", fontsize=20)
    tax.right_corner_label("G2", fontsize=20)
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    tax.show()


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
