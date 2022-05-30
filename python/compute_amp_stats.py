from collections import defaultdict
import operator
from functools import reduce
import warnings
import matplotlib.pyplot as plt
import pandas as pd
from lsst.afw.cameraGeom import ReadoutCorner
import lsst.afw.math as afw_math
import lsst.geom


__all__ = ['compute_cpBias_stats', 'compute_eotest_stats', 'plot_amp_stats']


def output_corner_bbox(amp, dx=200, dy=200):
    bbox = amp.getBBox()
    if amp.getReadoutCorner() == ReadoutCorner.LL:
        xmin, xmax = bbox.minX, bbox.minX + dx
        ymin, ymax = bbox.minY, bbox.minY + dy
    elif amp.getReadoutCorner() == ReadoutCorner.LR:
        xmin, xmax = bbox.maxX - dx, bbox.maxX
        ymin, ymax = bbox.minY, bbox.minY + dy
    elif amp.getReadoutCorner() == ReadoutCorner.UR:
        xmin, xmax = bbox.maxX - dx, bbox.maxX
        ymin, ymax = bbox.maxY - dy, bbox.maxY
    return lsst.geom.Box2I(lsst.geom.Point2I(xmin, ymin),
                           lsst.geom.Point2I(xmax, ymax))


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


def compute_cpBias_stats(butler, dsrefs, outfile=None, output_corner=False):
    data = defaultdict(list)
    for dsref in dsrefs:
        exposure = butler.get(dsref)
        det = exposure.getDetector()
        image = exposure.getImage()
        md = exposure.getMetadata()
        for amp, amp_info in enumerate(det, 1):
            if output_corner:
                bbox = output_corner_bbox(amp_info)
            else:
                bbox = amp_info.getBBox()
            amp_image = image.Factory(image, bbox)
            data = amp_stats(data, amp, amp_image, md)
    df = pd.DataFrame(data=data)
    if outfile is not None:
        df.to_pickle(outfile)
    return df


def compute_eotest_stats(fits_files, outfile=None, output_corner=False):
    output_corner_bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                         lsst.geom.Point2I(200, 200))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import lsst.eotest.sensor as sensorTest
    data = defaultdict(list)
    for item in fits_files:
        ccd = sensorTest.MaskedCCD(item)
        if output_corner:
            bbox = output_corner_bbox
        else:
            bbox = ccd.amp_geom.imaging
        for amp in ccd:
            image = ccd[amp].getImage()
            amp_image = image.Factory(image, bbox)
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
