from collections import defaultdict
import operator
from functools import reduce
import warnings
import matplotlib.pyplot as plt
import pandas as pd
import lsst.afw.math as afw_math


__all__ = ['compute_cpBias_stats', 'compute_eotest_stats', 'plot_amp_stats']


def amp_stats(data, amp, amp_image, md):
    flags = dict(mean=afw_math.MEAN, stdev=afw_math.STDEV,
                 meanclip=afw_math.MEANCLIP, stdevclip=afw_math.STDEVCLIP)
    data['amp'].append(amp)
    data['tseqnum'].append(md.get('TSEQNUM'))
    stats = afw_math.makeStatistics(amp_image,
                                    reduce(operator.ior, flags.values()))
    for flag_name, stat in flags.items():
        data[flag_name].append(stats.getValue(stat))
    return data


def compute_cpBias_stats(butler, dsrefs, det, outfile=None):
    data = defaultdict(list)
    for dsref in dsrefs:
        exposure = butler.get(dsref)
        image = exposure.getImage()
        md = exposure.getMetadata()
        for amp, amp_info in enumerate(det, 1):
            amp_image = image.Factory(image, amp_info.getBBox())
            data = amp_stats(data, amp, amp_image, md)
    df = pd.DataFrame(data=data)
    if outfile is not None:
        df.to_pickle(outfile)
    return df


def compute_eotest_stats(fits_files, outfile=None):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import lsst.eotest.sensor as sensorTest
    data = defaultdict(list)
    for item in fits_files:
        ccd = sensorTest.MaskedCCD(item)
        for amp in ccd:
            image = ccd[amp].getImage()
            amp_image = image.Factory(image, ccd.amp_geom.imaging)
            data = amp_stats(data, amp, amp_image, ccd.md)
    df = pd.DataFrame(data=data)
    if outfile is not None:
        df.to_pickle(outfile)
    return df


def plot_amp_stats(df0, title=None, outfile=None):
    columns = 'mean stdev meanclip stdevclip'.split()
    amps = sorted(list(set(df0['amp'])))
    plt.figure(figsize=(8, 8))
    for i, column in enumerate(columns, 1):
        plt.subplot(2, 2, i)
        for amp in amps:
            df = df0.query(f'amp=={amp}')
            plt.scatter(df['tseqnum'], df[column], s=2, label=str(amp))
        plt.xlabel('TSEQNUM')
        plt.ylabel(column)
        plt.legend(fontsize='x-small', ncol=4)
    plt.suptitle(title)
    plt.tight_layout()
    if outfile is not None:
        plt.savefig(outfile)
