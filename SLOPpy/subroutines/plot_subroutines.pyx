from __future__ import print_function, division
from SLOPpy.subroutines.common import *

def make_color_array(lists, input_data):
    """ Creation of the color array, based on the BJD of the observations
    """
    bjd = []
    am = []

    for obs in lists['observations']:
        bjd.append(input_data[obs]['BJD'] - 2450000.0)
        am.append(input_data[obs]['AIRMASS'])

    colors = np.asarray(bjd)
    cmap = plt.cm.viridis
    #cmap = plt.cm.Spectral
    line_colors = cmap(np.linspace(0, 1, len(lists['observations'])))

    return colors, cmap, line_colors

def make_color_array_matplotlib3(lists, input_data):
    """ Creation of the color array, based on the BJD of the observations
    """
    bjd = []
    mbjd = []
    am = []


    for obs in lists['observations']:
        bjd.append(input_data[obs]['BJD'])
        mbjd.append(input_data[obs]['mBJD'])
        am.append(input_data[obs]['AIRMASS'])

    color_cmap = plt.cm.viridis
    color_bjd_norm = plt.Normalize(vmin=bjd[0], vmax=bjd[-1])

    colors_bjd = color_cmap(color_bjd_norm(np.asarray(bjd)))

    color_am_norm = plt.Normalize(vmin=np.amin(am), vmax=np.amax(am))
    colors_am = color_cmap(color_am_norm(np.asarray(am)))

    colors_plot = {
        'BJD' : {},
        'mBJD' : {},
        'AIRMASS' : {}
    }
    colors_scatter = {
        'BJD' : {},
        'mBJD' : {},
        'AIRMASS' : {}
    }
    colors_properties = {
        'norm' : {
            'BJD': plt.Normalize(vmin=bjd[0], vmax=bjd[-1]),
            'mBJD': plt.Normalize(vmin=mbjd[0], vmax=mbjd[-1]),
            'AIRMASS': plt.Normalize(vmin=np.amin(am), vmax=np.amax(am))
        },
        'cmap' : plt.cm.viridis
    }

    for obs in lists['observations']:
        colors_plot['mBJD'][obs] = colors_properties['cmap'](
            colors_properties['norm']['mBJD'](input_data[obs]['mBJD']))
        colors_plot['BJD'][obs] = colors_properties['cmap'](
            colors_properties['norm']['BJD'](input_data[obs]['BJD']))
        colors_plot['AIRMASS'][obs] = colors_properties['cmap'](
            colors_properties['norm']['AIRMASS'](input_data[obs]['AIRMASS']))

        colors_scatter['mBJD'][obs] = [colors_properties['cmap'](
            colors_properties['norm']['mBJD'](input_data[obs]['mBJD']))[:-1]]
        colors_scatter['BJD'][obs] = [colors_properties['cmap'](
            colors_properties['norm']['BJD'](input_data[obs]['BJD']))[:-1]]
        colors_scatter['AIRMASS'][obs] = [colors_properties['cmap'](
            colors_properties['norm']['AIRMASS'](input_data[obs]['AIRMASS']))[:-1]]

    return colors_properties, colors_plot, colors_scatter

def grid_1plot():
    fig = plt.figure(figsize=(12, 6))
    gs = GridSpec(1, 2, width_ratios=[50, 1])

    ax = plt.subplot(gs[0, 0])
    cbax1 = plt.subplot(gs[0, 1])

    fig.subplots_adjust(wspace=0.04, hspace=0.25)

    return fig, gs, cbax1, ax

def grid_2plot(sharex=True, sharey=True):
    fig = plt.figure(figsize=(12, 6))
    gs = GridSpec(2, 2, width_ratios=[50, 1])

    ax1 = plt.subplot(gs[0, 0])

    if sharey and sharex:
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    elif sharex:
        ax2 = plt.subplot(gs[1, 0], sharex=ax1)
    elif sharey:
        ax2 = plt.subplot(gs[1, 0], sharey=ax1)
    else:
        ax2 = plt.subplot(gs[1, 0])

    cbax1 = plt.subplot(gs[:, 1])

    fig.subplots_adjust(wspace=0.04, hspace=0.25)

    return fig, gs, cbax1, ax1, ax2

def grid_3plot_small(sharex=False, sharey=False, partial_share=True):
    fig = plt.figure(figsize=(12, 9))
    gs = GridSpec(3, 2, width_ratios=[50, 1], height_ratios = [3, 1, 3])

    ax1 = plt.subplot(gs[0, 0])

    if sharey and sharex:
        ax2 = plt.subplot(gs[1, 0], sharex=ax1, sharey=ax1)
        ax3 = plt.subplot(gs[2, 0], sharex=ax1, sharey=ax1)
    elif sharex:
        ax2 = plt.subplot(gs[1, 0], sharex=ax1)
        ax3 = plt.subplot(gs[2, 0], sharex=ax1)
    elif sharey:
        ax2 = plt.subplot(gs[1, 0], sharey=ax1)
        ax3 = plt.subplot(gs[2, 0], sharey=ax1)
    elif partial_share:
        ax2 = plt.subplot(gs[1, 0], sharex=ax1)
        ax3 = plt.subplot(gs[2, 0], sharex=ax1, sharey=ax1)
    else:
        ax2 = plt.subplot(gs[1, 0])
        ax3 = plt.subplot(gs[2, 0])

    cbax1 = plt.subplot(gs[:, 1])

    fig.subplots_adjust(wspace=0.04, hspace=0.25)

    return fig, gs, cbax1, ax1, ax2, ax3