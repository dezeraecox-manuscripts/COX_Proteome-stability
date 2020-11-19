import os, re, string
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import bokeh.palettes as bp

# from pandas_profiling import ProfileReport

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

# ## Import and cleanup raw plate reader data for each timepoint
input_folder = 'raw_data/recombinant_denaturation/'
output_path = 'results/recombinant_denaturation/initial_cleanup/'

if not os.path.exists(output_path):
    os.makedirs(output_path)

## Generating standard information - labels and colours
samples = ['TRP', 'TRP', 'TRP', 'TPE', 'TPE', 'TPE', 'TRP_control', 'TPE_control']
sample_pos = [x for x in string.ascii_uppercase[0:len(samples)]]
sample_map = dict(zip(sample_pos, samples))

urea_conc = list(np.arange(0, 6.5, 0.5))
urea_conc.pop(-2)
urea_pos = np.arange(1, len(urea_conc)+1, 1)
urea_map = dict(zip(urea_pos, urea_conc))
colours = list(reversed(bp.magma(12)))

channel_dict = {'350-20/465-30': "TPE", 'F: 295-10/F: 360-20 2': 'TRP'}

# generating list of files, selecting kinetic read for 4 h timepoint
file_list = [filename for filename in os.listdir(input_folder) if 'lac' in filename]

for filename in file_list:
    sample_name = filename.split('_')[0]
    # import raw data, label dataframe with appropriate indices
    raw_data = pd.read_excel(f'{input_folder}{filename}', skiprows=np.arange(0, 13, 1))
    raw_data.columns.tolist()
    # group by channel info, adjust column names etc
    sample_data = {}
    for key, value in channel_dict.items():
        # select only columns containing this channel
        channel_cols = [col for col in raw_data.columns.tolist() if key in col]
        channel_data = raw_data.set_index(['Well\nRow', 'Well\nCol'])[channel_cols].copy()
        time_cols = [time.split(' ')[0] for time in list(channel_data.reset_index().iloc[0].dropna())]
        channel_data.columns = time_cols
        channel_data = channel_data.reset_index().drop(0)
        channel_data[time_cols] = channel_data[time_cols].astype(float)
        channel_data['urea'] = channel_data['Well\nCol'].map(urea_map)
        channel_data['samples'] =channel_data['Well\nRow'].map(sample_map)
        # select only samples for which that channel is relevant
        channel_data = channel_data[channel_data.samples.str.contains(value,case=False)]
        sample_data[value] = channel_data

    # Save cleaned df to excel
    data_frames = [df.sort_values(['samples', 'urea']) for df in sample_data.values()]
    sheetnames = [key for key in sample_data.keys()]
    FileHandling.df_to_excel(data_frames=data_frames, sheetnames=sheetnames, output_path=f'{output_path}{sample_name}_cleaned_kinetic.xlsx')
    