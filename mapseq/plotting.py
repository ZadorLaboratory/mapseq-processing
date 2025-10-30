import logging
import os

import datetime as dt
import math

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from natsort import natsorted

logging.getLogger('matplotlib.font_manager').disabled = True

def make_counts_plots(df, 
                      outdir=None, 
                      groupby='label', 
                      type=None, 
                      column='read_count',
                      min_count=1,  
                      cp=None):
    '''
    take standard aggregated, readtable or vbctable DFs and create 
    read_count or umi_count frequency plots for all real targets.  
    
    confirm that groupby column value exists. 
    
    if type is None, include all. 
    
    '''
    project_id = cp.get('project','project_id')
    
    min_count = int(min_count)
    if min_count > 1:
        logging.debug(f'min_count > 1, thresholding...')
        df = df[df[column] >= min_count]
        df.reset_index(inplace=True, drop=True)
    
    before = len(df)
    df = df[df[groupby] != '']
    after = len(df)
    removed = before - after
    logging.debug(f'removed {removed} rows ( {before} - {after}) with no value for {groupby}')

    if type != None:
        before = len(df)
        df = df[df['type'] != type]
        after = len(df)
        removed = before - after
        logging.debug(f'removed {removed} rows ( {before} - {after}) with type != {type}')
    else:
        type = 'all'    
    
    make_freqplot_combined_sns(df, 
                               title=f'{project_id}:{type} {column} frequency',  
                               outfile=os.path.join(outdir, f'{project_id}_{type}_{column}_by{groupby}_freq.c{min_count}.pdf'),
                               groupby=groupby, 
                               column=column,
                               scale='log10' )


def make_freqplot_combined_sns(df, 
                           title='Frequency',  
                           outfile='frequency-plots.pdf',
                           groupby='label', 
                           column='read_count',
                           scale=None ):    
    '''
     makes combined figure with all plots.
      
    scale = log10 | log2 | None  
      
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    datestr = dt.datetime.now().strftime("%Y%m%d%H%M")

    if scale is not None:
        title = f'{title} ({scale})'
    
    groups=natsorted( list(df[groupby].unique()) )
    page_dims = (11.7, 8.27)

    with pdfpages(outfile) as pdfpages:
        #fig_n = math.ceil( math.sqrt(len(filelist)) )
        #fig, axes = plt.subplots(nrows=fig_n, ncols=fig_n, figsize=a4_dims,  layout='constrained')
        plots_per_page = 9
        num_figs = float(len(groups)) / float(plots_per_page)
        if num_figs % 9 == 0:
            num_figs = int(num_figs)
        else:
            num_figs = int(num_figs) + 1
        logging.debug(f'with {plots_per_page} plots/page, need {num_figs} for {len(groups)} groups.')
        
        figlist = []
        axlist = []
        for i in range(0,num_figs):
            fig,axes = plt.subplots(nrows=3, ncols=3, figsize=page_dims,  layout='constrained')
            fig.suptitle(title)
            figlist.append(fig)
            # numpy.flatirator doesn't handle indexing
            for a in axes.flat:
                axlist.append(a)
        logging.debug(f'created {len(groups)} figures to go on {num_figs} pages. ')
        
 
        for i, group in enumerate(groups):
            gdf = df[ df[groupby] == group ]                
            logging.debug(f'handling {group} length={len(gdf)}')
            ax = axlist[i]
            counts_axis_plot_sns(ax, gdf, column=column, title=f'{group} {column}', scale=scale)
                
        for f in figlist:
            pdfpages.savefig(f)
    logging.info(f'saved plot PDF to {outfile}')

def make_freqplot_single_sns(df, 
                           title='Frequency',  
                           outfile='frequency-plot.pdf',
                           column='read_count',
                           scale=None ):    
    '''
    single figure with one plot. 
    scale = log10 | log2 | None  
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    datestr = dt.datetime.now().strftime("%Y%m%d%H%M")

    if scale is not None:
        title = f'{title} ({scale})'
    
    page_dims = (8, 6)
    with pdfpages(outfile) as pdfpages:
        fig, axes = plt.subplots(figsize=page_dims)
        fig.suptitle(title)
        counts_axis_plot_sns(axes, df, column=column, title=f'{column} freqplot', scale=scale)
        pdfpages.savefig(fig)
    logging.info(f'saved plot PDF to {outfile}')

def counts_axis_plot_sns(ax, df, scale=None, column='read_count', title='counts frequency' ) :
    '''
    Creates individual axes for single plot within figure. 
    scale = None | log10  | log2
        # Set the y-axis to a logarithmic scale
        plt.yscale('log')
        
        # Customize y-axis ticks (optional)
        # For major ticks:
        major_ticks = [1, 10, 100, 1000, 10000, 100000, 1000000]
        plt.yticks(major_ticks)
        
        # For minor ticks:
        ax = plt.gca() # Get current axes
        locmin = mticker.LogLocator(base=10.0, subs=np.arange(0.1, 1, 0.1), numticks=10)
        ax.yaxis.set_minor_locator(locmin)
        ax.yaxis.set_minor_formatter(mticker.NullFormatter()) # Hide minor tick labels if desired
        
        plt.title('Seaborn Line Plot with Logarithmic Y-axis')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis (Log Scale)')
        plt.grid(True, which="both", ls="-") # Add grid for both major and minor ticks
        plt.show()    
    

    '''
    logging.debug(f'column={column} scale={scale} title={title}')
    df.sort_values(by=[column], ascending=False, inplace=True)
    df.reset_index(inplace=True, drop=True)
    df.reset_index(inplace=True)
    
    s = df[column].sum()
    n = len(df)
    t = df[column].max()
    r = df['index'].max()
    h = calc_freq_threshold(df, fraction=0.85, column = column)

    logging.debug(f'making lineplot...')
    sns.lineplot(ax=ax, data=df, x=df['index'], y=df[column])
    logging.debug(f'done. calc bounds...')

    lx = df['index'].max()
    ly = df[column].max()

    ax.set_xlabel('Rank')

    if scale == 'log10':
        ax.set_yscale('log')
        #major_ticks = [1, 10, 100, 1000, 10000, 100000, 1000000]
        major_ticks = make_logticks(ly)
        ax.set_yticks(major_ticks)
        ax.set_ylabel(f'log10( {column} )')
        title = f'{title} log10().'
        #ly = math.log10(ly)       
        logging.debug(f'made axis with log scale y-axis.')
    else:
        ax.set_ylabel(f'{column}')
        logging.debug(f'made axis with no scaling.')
    
    ax.text(lx, ly, s=f"n={n}\ntop={t}\nsum={s}\nest_0.85_threshold={h}",
            fontsize=11, 
            horizontalalignment='right',
            verticalalignment='top'            
            )
    ax.set_title(title, fontsize=10)
    
def make_logticks(max_value):
    '''
    #major_ticks = [1, 10, 100, 1000, 10000, 100000, 1000000]
    '''
    ticklist = []
    i = 1
    while i < max_value:
        ticklist.append(i)
        i *= 10
    ticklist.append(i)
    return ticklist

def counts_axis_plot_sns_old(ax, df, scale=None, column='read_count', title='counts frequency' ) :
    '''
    Creates individual axes for single plot within figure. 
    scale = None | log10  | log2

    '''
    df.sort_values(by=[column], ascending=False, inplace=True)
    df.reset_index(inplace=True, drop=True)
    df.reset_index(inplace=True)
    
    s = df[column].sum()
    n = len(df)
    t = df[column].max()
    r = df['index'].max()
    h = calc_freq_threshold(df, fraction=0.9, column = column)
    
    if scale is None:
        sns.lineplot(ax=ax, x=df['index'], y=df[column] )
        lx = df['index'].max()
        ly = df[column].max()
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\nest_90pct_threshold={h}", 
                fontsize=11, 
                horizontalalignment='right',
                verticalalignment='top',) #add text
        logging.debug(f'made axis without scale.') 

    elif scale == 'log10':
        #  Avoid divide by 0 runtime warning...
        #  Switch to non-log-scaled X axis
        df['n_index'] = df['index'] + 1
        df['log10counts'] = np.log10(df[column])
        sns.lineplot(ax=ax, x=df['n_index'], y=df['log10counts'] )
        #lx = 0.05 * np.log10(t) 
        #ly = 0.05 * np.log10(r)        
        #lx = 0.1
        lx = df['index'].max() * 0.8
        ly = df['log10counts'].max() * 0.8
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\nest_90pct_threshold={h}", 
                fontsize=11,
                horizontalalignment='right',
                verticalalignment='top',) #add text
        title = f'{title} log10().'
        logging.debug(f'made axis with log10 scale.')       

    ax.set_title(title, fontsize=10)


def calc_freq_threshold(df, fraction=0.9, column = 'read_count'):
    '''
    sorts column of input column
    calculates index of point at which <fraction> of data points are less 
    returns column value at that point + 1 
    '''
    ser = df[column].copy()
    ser.sort_values(ascending = False, inplace=True)
    ser.reset_index(drop=True, inplace=True)
    idx = int(len(ser) * fraction)
    t = int( ser.iloc[idx] + 1 )    
    return t


def counts_freq(matlabfile, logscale = 'log10', logcol = 'counts' ):
    '''
    '''
    df = pd.read_csv('barcodematrix.tsv',sep='\t', header=None)
    rowsum = df.sum(axis=1)
    rowsort = rowsum.sort_values(ascending=False)
    df = pd.DataFrame(data=rowsort)
    df.columns = ['counts']
    df  = df.reset_index(drop=True)
       
    df[f'log_{logcol}'] = np.log10(df[f'{logcol}'])

    
    #if logscale == 'log2':
    #    df = np.log2(df)
    #elif logscale == 'log10':
    #    df = np.log10(df)
    
    ax = sns.lineplot(data=df, x=df.index, y=df[ f'log_{logcol}' ])
    ax.set_title('HZ120Vamp2 counts frequency')
    ax.set(xlabel='Sequence Rank', ylabel='log10(BC molecule count)')
    
    #   OR 
    # plt.xlabel('x-axis label')
    # plt.ylabel('y-axis label')
    
    plt.savefig('Z120Vamp2.log10.countsfreq.png')
    #plt.show()


def make_simple_freqplot(df, 
                         column='read_count', 
                         title='freqplot', 
                         outfile = 'freqplot.pdf',
                         xscale=None,
                         yscale='log10' ):
    '''
    Plot a single column as frequency plot, with y-axis scale. 
    
    scale = log10 | log2 | None 

    '''
    ax = None
    
    df['n_index'] = pd.Series( df.index ) + 1

    s = df[column].sum()
    n = len(df)
    t = df[column].max()
    r = df['n_index'].max()

    df['log10index'] = np.log10(df['n_index'])
    df['log10counts'] = np.log10( df[column] )

    if yscale is None:
        logging.debug(f'making plot with no y scaling.')
        if xscale is None:
            ax = sns.lineplot(x=df['n_index'], y=df[column] )
            logging.debug(f'made axis without x or y axis scaling.')
            ax.set(xlabel='Sequence Rank')
        elif xscale == 'log10':
            ax = sns.lineplot(x=df['log10index'], y=df[column] )
            ax.set(xlabel='log10( sequence rank )') 
            logging.debug(f'made axis with x-axis only scaling.')
        ax.set(ylabel=column)
        lx = df['index'].max()
        ly = df[column].max()
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\n", 
                fontsize=11, 
                horizontalalignment='right',
                verticalalignment='top',) #add text
 
    elif yscale == 'log10':
        #  Avoid divide by 0 runtime warning...
        #  Switch to non-log-scaled X axis
        logging.debug(f'making plot with y scaling = {yscale} .')
        if xscale is None:
            ax = sns.lineplot(x=df['n_index'], y=df['log10counts'] )
            logging.debug(f'made axis with y-axis only scaling.')     
            ax.set(xlabel='sequence rank')
        elif xscale == 'log10':
            ax = sns.lineplot(x=df['log10index'], y=df['log10counts'] )
            logging.debug(f'made axis with x-axis and y-axis scaling.')
            ax.set(xlabel='log10( sequence rank )')
        ax.set(ylabel=f'log10( {column} )')
        lx = 0.05 * np.log10(t) 
        ly = 0.05 * np.log10(r)        
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\n", fontsize=11) #add text
              
    ax.set_title(title)
    
    of = os.path.abspath( outfile )
    logging.info(f'saving to {of} ')       
    plt.savefig(of)
    plt.clf()
    logging.debug(f'done')
    return ax

   