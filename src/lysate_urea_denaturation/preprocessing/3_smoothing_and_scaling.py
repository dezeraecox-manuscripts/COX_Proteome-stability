import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from GEN_Utils import FileHandling
from pykalman import KalmanFilter
from scipy.interpolate import UnivariateSpline
from scipy.stats import ttest_1samp

from loguru import logger
logger.info('Import OK')


def one_sample_ttest(compiled, sample_cols, popmean=1):
    df = compiled.copy()
    ttest_results = []
    for sequence, df in df.groupby('Sequence'):
        results = []
        for col in sample_cols:
            test_vals = df[col].values
            if len(test_vals) > 1:
                results.append(tuple(ttest_1samp(test_vals, popmean=popmean, nan_policy='omit')))
            else:
                results.append(tuple([np.nan, np.nan]))
        results = pd.DataFrame(results)
        results.columns = ['t-stat', 'p-val']
        results['conc'] = sample_cols
        results['Sequence'] = sequence
        ttest_results.append(results)
    ttest_results = pd.concat(ttest_results)
    ttest_results[['t-stat', 'p-val']] = ttest_results[['t-stat', 'p-val']].astype(float)

    return ttest_results # 5% of the points detected as significantly different


def coverage_filter(compiled, channel_threshold):
    """Drop peptides with < threshold coverage"""
    df = compiled.copy()
    df['coverage'] = df[sample_cols].count(axis=1)
    return df[df['coverage'] >= channel_threshold] # 458 peptides covering 265 proteins


def replicate_filter(compiled, replicate_threshold):
    """collect only those seen in threshold replicates"""
    df = compiled.copy()
    replicates = df.groupby('Sequence').count()['Proteins']
    rep_sequences = replicates[replicates == replicate_threshold].reset_index()['Sequence']
    return df[df['Sequence'].isin(rep_sequences)] # 92 proteins, covered by 125 peptides in three replicates


def significance_filter(compiled, test_results, threshold=0.05):
    """Collects peptides containing at least one significant sample from long-form ttest results and filters from compiled df"""
    
    ttest_results = one_sample_ttest(compiled, sample_cols)
    sig_peptides = list(set(ttest_results[ttest_results['p-val'] < threshold]['Sequence'])) # 519 peptides
    
    return compiled[compiled['Sequence'].isin(sig_peptides)]


def symmetrical_normalisation(compiled, sample_cols):
    """obtain symmetric values for equal increases and decreases"""
    df = compiled.copy()

    def symmetric_ratio(ratio):
        if ratio > 1:
            return (-1/ratio) + 1
        elif ratio < 1:
            return (ratio) - 1
        elif ratio == 1:
            return 0
        else:
            return np.nan

    df[sample_cols] = df[sample_cols].applymap(symmetric_ratio)
    return df


def min_max_scaler(compiled, sample_cols, max_range=1):

    def scaler(ratios, max_range=1):
        if any(abs(ratio) > max_range for ratio in ratios):
            return ratios / np.max(abs(ratios))
        else:
            return ratios
    scaled = compiled.copy()
    scaled[sample_cols] = scaled[sample_cols].T.apply(scaler, args=(max_range,)).T

    return scaled


def st_dev_normalisation(compiled, sample_cols):
    """Standardize by dividing the peptide ratio for each concentration by the SD for this peptide overall"""
    df = compiled.copy()
    df[sample_cols] = (df[sample_cols].T / df[sample_cols].std(axis=1)).T
    return df


def ewm_smoothing(compiled, max_smooth, min_smooth=2, plot=False):
    smooth = []
    for sequence, df in compiled.groupby('Sequence'):
        for_smoothing = pd.melt(df, id_vars=None, value_vars=sample_cols, var_name='channel', value_name='ratio')
        for_smoothing['channel'] = for_smoothing['channel'].astype(int)
        for_smoothing['Sequence'] = sequence
        for_smoothing = for_smoothing.sort_values('channel')
        # sns.lineplot(data=for_smoothing, x='channel', y='ratio',  ci='sd', label='raw')
        for smooth_factor in np.arange(min_smooth, max_smooth):
            for_smoothing[f'EMA_{smooth_factor}'] = for_smoothing.loc[:,'ratio'].ewm(com=smooth_factor, adjust=True).mean()
            if plot:
                sns.lineplot(data=for_smoothing, x='channel', y=f'EMA_{smooth_factor}',  ci='sd', label=f'{smooth_factor}')
        if plot:
            plt.title(f'{sequence} - com')
            plt.legend(bbox_to_anchor=(1.0, 1.0))
            # plt.ylim(0, 2)
            plt.show()
        smooth.append(for_smoothing)
    smooth = pd.concat(smooth)

    return pd.pivot_table(smooth, values='EMA_3', index='Sequence', columns='channel')


def loess_smoothing(dataframe, degrees=(2, 5), plot=False, v=None):
    # testing custom LOWESS method
    def loc_eval(x, b):
        loc_est = 0
        for i in enumerate(b): loc_est+=i[1]*(x**i[0])
        return(loc_est)


    def loess(xvals, yvals, data, alpha, poly_degree=1, v=None):
        all_data = sorted(zip(data[xvals].tolist(), data[yvals].tolist()), key=lambda x: x[0])
        xvals, yvals = zip(*all_data)
        evalDF = pd.DataFrame(columns=['v','g'])
        n = len(xvals)
        m = n + 1
        q = int(np.floor(n * alpha) if alpha <= 1.0 else n)
        if not v:
            avg_interval = ((max(xvals)-min(xvals))/len(xvals))
            v_lb = min(xvals)-(.5*avg_interval)
            v_ub = (max(xvals)+(.5*avg_interval))
            v = enumerate(np.linspace(start=v_lb, stop=v_ub, num=m), start=1)
        xcols = [np.ones_like(xvals)]
        for j in range(1, (poly_degree + 1)):
            xcols.append([i ** j for i in xvals])
        X = np.vstack(xcols).T
        for i in v:
            iterpos = i[0]
            iterval = i[1]
            iterdists = sorted([(j, np.abs(j-iterval)) for j in xvals], key=lambda x: x[1])
            _, raw_dists = zip(*iterdists)
            scale_fact = raw_dists[q-1]
            scaled_dists = [(j[0],(j[1]/scale_fact)) for j in iterdists]
            weights = [(j[0],((1-np.abs(j[1]**3))**3 if j[1]<=1 else 0)) for j in scaled_dists]
            _, weights      = zip(*sorted(weights,     key=lambda x: x[0]))
            _, raw_dists    = zip(*sorted(iterdists,   key=lambda x: x[0]))
            _, scaled_dists = zip(*sorted(scaled_dists,key=lambda x: x[0]))
            W         = np.diag(weights)
            b         = np.linalg.inv(X.T @ W @ X) @ (X.T @ W @ yvals)
            local_est = loc_eval(iterval, b)
            iterDF2   = pd.DataFrame({
                        'v'  :[iterval],
                        'g'  :[local_est]
                        })
            evalDF = pd.concat([evalDF, iterDF2])
        evalDF = evalDF[['v','g']]
        return(evalDF)


    avg_interval = ((max(sample_cols)-min(sample_cols))/len(sample_cols))
    v_lb = min(sample_cols)-(.5*avg_interval)
    v_ub = (max(sample_cols)+(.5*avg_interval))
    loess_vals = []
    for sequence, df in dataframe.groupby('Sequence'):
        test = pd.melt(df, id_vars=None, value_vars=sample_cols, var_name='channel', value_name='ratio')
        test = test.sort_values('channel')
        test.dropna(inplace=True, axis=0)
        
        if plot:
            fig, ax = plt.subplots()
        for degree in np.arange(*degrees):
            v = enumerate(np.linspace(start=v_lb, stop=v_ub, num=100), start=1)
            evalDF = loess(xvals='channel', yvals='ratio', data=test, alpha=0.7, poly_degree=degree, v=v)
            if plot:
                plt.plot(evalDF['v'], evalDF['g'], linewidth= 3, label=degree)
        evalDF['Sequence'] = sequence
        loess_vals.append(evalDF)
        if plot:
            sns.lineplot(test['channel'], test['ratio'], color='grey', label='mean')
            sns.scatterplot(test['channel'], test['ratio'], color='grey', label='raw')
            plt.title(sequence)
            plt.legend()
            plt.tight_layout()
            # plt.savefig(f'{output_folder}degree_3/{sequence}.png')
            plt.show()
            plt.clf()

    loess_smoothed = pd.concat(loess_vals)
    loess_smoothed = pd.pivot_table(loess_smoothed, values='g', index='v', columns='Sequence')
    # reset intercept to be 0
    loess_smoothed = loess_smoothed - loess_smoothed.iloc[0]

    return loess_smoothed.T


def kalman_smoothing(compiled, sample_cols, plot=False):
    for sequence, df in compiled.groupby('Sequence'):
        test = pd.melt(df, id_vars=None, value_vars=sample_cols, var_name='channel', value_name='ratio')
        test = test.sort_values('channel')
        test_vals = test.groupby('channel').mean()['ratio'].values
        kf = UnscentedKalmanFilter()
        (filtered_state_means, filtered_state_covariances) = kf.filter(test_vals)
        (smoothed_state_means, smoothed_state_covariances) = kf.smooth(test_vals)
        if plot:
            sns.lineplot(x=sample_cols, y=pd.DataFrame(test_vals)[0], label='raw')
            sns.lineplot(x=sample_cols, y=pd.DataFrame(filtered_state_means)[0].values, label='filtered')
            sns.lineplot(x=sample_cols, y=pd.DataFrame(smoothed_state_means)[0].values, label='smooth')
            plt.title(sequence)
            plt.legend()
            plt.ylim(0, 2)
            plt.show()
            plt.clf()
        return smoothed_state_means


def desplined_smoothing(compiled, sample_cols, plot=False):
    desplined = {}
    xs = np.linspace(np.min(sample_cols), np.max(sample_cols), 100)
    test_sequences = list(set(compiled['Sequence']))
    for sequence, df in compiled.groupby('Sequence'):
        test = pd.melt(df, id_vars=None, value_vars=sample_cols, var_name='channel', value_name='ratio')
        test = test.sort_values('channel').dropna()
        y = test['ratio']
        x = test['channel']
        s = UnivariateSpline(x, y, k=4)
        ys = s(xs)

        desplined[sequence] = ys

        if plot:
            sns.lineplot(test['channel'], test['ratio'], ci='sd')
            plt.plot(xs, ys)
            plt.title(sequence)
            plt.show()
            plt.clf()
    desplined = pd.DataFrame(desplined).T
    desplined.columns = xs

    return desplined


def updown_smoothing(compiled, sample_cols, threshold=0.05, center=0):
    ttest_results = one_sample_ttest(compiled, sample_cols, popmean=center)
    p_vals = pd.pivot_table(ttest_results, values='p-val', index='Sequence', columns='conc')
    # if significantly different (i.e. < 0.05) then record point as increase or decrease
    updown = compiled.groupby('Sequence').mean().copy().sort_values('Sequence')
    updown = updown.where(p_vals.sort_values('Sequence') < threshold)
    def direction(value):
        if value > center:
            return center + 1
        elif value < center:
            return center - 1
        else:
            return np.nan

    updown[sample_cols] = updown[sample_cols].applymap(direction).fillna(center)
    return updown


def pval_smoothing(compiled, sample_cols):
    """Record point as mean * 1/p-val"""
    ttest_results = one_sample_ttest(compiled, sample_cols)
    ttest_results['exp_p-val'] = - 20**ttest_results['p-val']
    p_vals = pd.pivot_table(ttest_results, values='exp_p-val', index='Sequence', columns='conc')

    proportional_pval = compiled.groupby('Sequence').mean().copy().sort_values('Sequence')
    proportional_pval[sample_cols] = (1 - proportional_pval[sample_cols]) * (1 / p_vals.sort_values('Sequence'))
    proportional_pval[0.0] = 0

    return proportional_pval


def stdev_smoothing(compiled, sample_cols):
    """Record point as mean * 1/SD of that sample for that peptide"""

    proportional_SD = compiled.groupby('Sequence').mean().copy().sort_values('Sequence')
    proportional_SD[sample_cols] = proportional_SD[sample_cols] * \
        (1 / compiled.groupby('Sequence').std().copy().sort_values('Sequence')[sample_cols])

    return proportional_SD


def sliding_pval_smoothing(compiled, sample_cols, win_width=3):
    """Apply sliding window around central point for smoothing, then calculate SE of change
    - this helps with peptides only found in one replicate for which SE cannot be calculated"""

    ttest_results = one_sample_ttest(compiled, sample_cols)
    ttest_results['exp_p-val'] = - 20**ttest_results['p-val']
    smoothed_window = []
    for sequence in ttest_results['Sequence'].unique():
        df = ttest_results[ttest_results['Sequence'] == sequence]
        first_val = np.mean(df['exp_p-val'][:3])
        final_val = np.mean(df['exp_p-val'][-2:])
        df['smooth'] = df['exp_p-val'].rolling(win_width, win_type ='triang', center=True, min_periods=1).mean()
        smoothed_window.append(df)
    smoothed_window = pd.concat(smoothed_window)
    p_vals = pd.pivot_table(smoothed_window, values='smooth', index='Sequence', columns='conc')

    proportional_pval = compiled.groupby('Sequence').mean().copy().sort_values('Sequence')
    proportional_pval[sample_cols] = (1 - proportional_pval[sample_cols]) * (1 / p_vals.sort_values('Sequence'))
    # remove any peptides which are all NaN (can be either single-replicate peptides, or those that were only quantified in a conc once)
    proportional_pval.dropna(inplace=True, subset=sample_cols, how='all')
    proportional_pval[0.0] = 0

    return proportional_pval


def visualise_transformation(test_dfs, sequences, range_vals=None):
    """Visualise what each transformation does to the data
    test_df: dict mapping str(data_type): (df, cols)"""
    for sequence in sequences:
        fig, ax = plt.subplots()
        for data_type, (df, cols) in test_dfs.items():
            data = df[df['Sequence'] == sequence]
            if data.shape[0] > 0:
                sns.lineplot(x='channel', y='ratio', data=pd.melt(data, id_vars=None, value_vars=cols, var_name='channel', value_name='ratio'), ci='sd', label=data_type)
        if range_vals:
            plt.ylim(*range_vals)
        plt.legend()
        plt.title(sequence)
        ax.axhline(1, c='r', linestyle='--')
        ax.axhline(0, c='grey', linestyle='--')
        plt.tight_layout()
        plt.show()
        plt.clf()


def population_plot(df, cols, title, yrange):
    df[cols].T.plot(legend=None)
    plt.ylim(*yrange)
    plt.title(title)
    plt.savefig(f'{output_folder}{title}.png')
    plt.show()


if __name__ == "__main__":
    pass

    input_path = f'results/lysate_denaturation/normalised/normalised_summary.xlsx'
    output_folder = f'results/lysate_denaturation/smoothing/'

    channel_threshold = 12
    replicate_threshold = 3
    filter_cols = ['12']

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    # read in raw data
    raw_data = pd.read_excel(input_path, sheet_name='control_norm_cys_ratio')
    compiled = raw_data.copy().drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1)
    compiled.drop('cys_rank', inplace=True, axis=1)
    # Drop channel 12 as it is a MESS
    sample_cols = [str(x) for x in np.arange(1, 14) if str(x) not in filter_cols]
    info_cols = ['Sequence', 'Proteins', 'replicate', 'coverage']
    compiled = compiled.set_index([col for col in compiled.columns.tolist() if col in info_cols])[sample_cols]
    # relabel columns with urea conc
    urea_conc = {'1': 0.0, '2': 0.5, '3': 1.0, '4': 1.5, '5': 2.0, '6': 2.5, '7': 3.0, '8': 3.5, '9': 4.0, '10': 4.5, '11': 5.0, '13': 6.0}
    compiled.rename(columns=urea_conc, inplace=True)
    sample_cols = compiled.columns.tolist()
    compiled.reset_index(inplace=True)

    raw = compiled.copy() #1575 peptides covering 750 peptides


    """--------------SMOOTHING-----------------------"""

    raw = compiled.copy()
    # population_plot(raw, sample_cols, 'raw_data', (0, 2))

    # Complete scaling pval methods
    # Gather smoothed data for all peptides
    pval_smooth = sliding_pval_smoothing(raw, sample_cols, win_width=3)
    pval_smooth = coverage_filter(pval_smooth, channel_threshold=6)
    # population_plot(pval_smooth, sample_cols, 'pval_smooth', (-2, 2))

    pval_scaled = min_max_scaler(pval_smooth, sample_cols, max_range=0.1)
    # population_plot(pval_scaled, sample_cols, 'pval_scaled', (-2, 2))

    pval_loess = loess_smoothing(dataframe=pval_scaled, degrees=(2, 3), plot=False)
    pval_loess_cols = pval_loess.columns.tolist()
    pval_loess = pval_loess.reset_index()
    # population_plot(pval_loess, pval_loess_cols, 'pval_loess', (-2, 2))

    # Visualise the transformation on individual sequences as examples
    test_sequences = pval_loess['Sequence'][0:10]
    visualise_transformation(
        test_dfs={
            'raw': (raw, sample_cols),   
            'pval_smooth': (pval_smooth.reset_index(), sample_cols),
            'pval_scaled': (pval_scaled.reset_index(), sample_cols),
            'pval_loess': (pval_loess.reset_index(), pval_loess_cols),         
            },
        sequences=test_sequences,
        range_vals=(-2, 2)
        )


    # save original and processed data to excel
    data_dict = {'raw': raw, 'pval_smooth': pval_smooth, 'pval_scaled': pval_scaled, 'pval_loess': pval_loess}
    FileHandling.df_to_excel(output_path=f'{output_folder}processed_data.xlsx', sheetnames=list(data_dict.keys()), data_frames=list(data_dict.values()))

