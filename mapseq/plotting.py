import logging
import os

import datetime as dt

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from natsort import natsorted

def make_counts_plots(df, outdir=None, groupby='label', type=None, column='read_count', cp=None):
    '''
    take standard aggregated, readtable or vbctable DFs and create 
    read_count or umi_count frequency plots for all real targets.  
    
    confirm that groupby column value exists. 
    
    if type is None, include all. 
    
    '''
    project_id = cp.get('project','project_id')
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
                               outfile=os.path.join(outdir, f'{project_id}_{type}_{column}_by{groupby}_frequency.pdf'),
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
    
    groups=natsorted( list(df[groupby].unique()) )

    if scale is not None:
        title = f'{title} ({scale})'

    # do nine per figure...
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



    
def counts_axis_plot_sns(ax, df, scale=None, column='read_count', title='counts frequency' ) :
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
    
    if scale is None:
        h = calc_freq_threshold(df, fraction=0.9, column = column)
        sns.lineplot(ax=ax, x=df['index'], y=df[column] )
        title = title

        lx = df['index'].max()
        ly = df[column].max()
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\nestimated_threshold={h}", 
                fontsize=11, 
                horizontalalignment='right',
                verticalalignment='top',) #add text
        logging.debug(f'made axis without scale.') 

    elif scale == 'log10':
        # avoid divide by 0 runtime warning...
        df['log10index'] = np.log10(df['index'] + 1)
        df['log10counts'] = np.log10(df[column])
        sns.lineplot(ax=ax, x=df['log10index'], y=df['log10counts'] )
        #lx = 0.05 * np.log10(t) 
        #ly = 0.05 * np.log10(r)        
        lx = 0.1
        ly = 0.1
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\n", fontsize=11) #add text
        title = f'{title} log10().'
        logging.debug(f'made axis with log10 scale.')       

    ax.set_title(title, fontsize=10)


   