import os
import glob
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lsst.obs.lsst
from focal_plane_plotting import plot_focal_plane


__all__ = ['extract_amp_metric', 'plot_amp_metric']


CAMERA = lsst.obs.lsst.LsstCam.getCamera()

CHANNELS = {amp: channel for amp, channel in
            enumerate('''C10 C11 C12 C13 C14 C15 C16 C17
                         C07 C06 C05 C04 C03 C02 C01 C00'''.split(), 1)}


def extract_amp_metric(stat_data_dir, suffix='.pickle'):
    amp_data = defaultdict(dict)
    for detector, det in enumerate(CAMERA):
        det_name = det.getName()
        stat_file = glob.glob(os.path.join(stat_data_dir,
                                           f'{det_name}_*' + suffix))[0]
        df0 = pd.read_pickle(stat_file).query('tseqnum > 19')
        amps = sorted(list(set(df0['amp'])))
        for amp in amps:
            channel = CHANNELS[amp]
            df = df0.query(f'amp=={amp}')
            amp_data[det_name][channel] = np.std(df['meanclip'])
    return amp_data


def plot_amp_metric(amp_data, stat_name, z_range=(0, 1.3), figsize=(9, 8),
                    outfile=None):
    plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    plot_focal_plane(ax, amp_data, camera=CAMERA, z_range=z_range)
    plt.title(f'{stat_name}, stdev(meanclip)')
    if outfile is None:
        outfile = f'{stat_name}_stdev_meanclip.png'
        outfile = outfile.replace(' ', '_')
    plt.savefig(outfile)
