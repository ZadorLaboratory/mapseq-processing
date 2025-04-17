import logging
import os

import datetime as dt

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from natsort import natsorted


def make_counts_plots(df, outdir=None, groupby='label', column='read_count', cp=None):
    '''
    take standard aggregated, readtable or vbctable DFs and create 
    read_count or umi_count frequency plots for all real targets.  
    
    
    '''
    make_freqplot_combined_sns(df, 
                               title=f'{column} frequency',  
                               outfile=os.path.join(outdir, f'{column}_by{groupby}_frequency.pdf'),
                               groupby=groupby, 
                               column=column,
                               scale='log10' )


def make_clustered_heatmap(df, outprefix, columns=None ):
    '''
    
    Caller should edit columns in order to exclude injection areas from plot. 
    '''
    camp = 'Reds'
    g = sns.clustermap(df, cmap=camp, yticklabels=False, col_cluster=False, standard_scale=0)
    g.fig.subplots_adjust(right=0.7)
    g.ax_cbar.set_position((0.8, .2, .03, .4))
    plt.title(f'{prefix}\nCounts')
    plt.savefig(f'{outprefix}.heatmap.pdf')
    logging.info(f'done making {outprefix}.heatmap.pdf ')
    

def make_read_countplot(config, df, outfile, title=None ): 
    '''
    makes individual read count plot from sequence read_count DF 
    assumes 'sequence' and 'read_count' columns. 
    
    '''   
    plt.set_loglevel (level = 'warning')
    
    logging.debug(f'handling sequence df len={len(df)}')
    outfile = os.path.abspath(outfile)    

    if title is None:
        title='Read count frequence plot.'

    df.sort_values(by='read_count', ascending=False, inplace=True)
    df.reset_index(inplace=True)
    df['index'] = df['index'].astype('int64')

    plt.figure()
    plt.plot(np.log10(df['index']), np.log10(df['read_count']))
    plt.title(title)
    plt.xlabel("log10(index)")
    plt.ylabel("log10(read_count)")
    logging.info(f'Saving count plot to {outfile}')
    plt.savefig(outfile)


def make_read_countplot_sns(cp, df, outfile='count-frequency-plot.pdf' ):
    '''
    makes two figures, one log-normalized, one straight for same counts df. 
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    df.sort_values(by='read_count', ascending=False, inplace=True)
    df.reset_index(inplace=True)
    df['index'] = df.index.astype('int64')

    #page_dims = 
    #  A4 landscape   (11.69,8.27) 
    #  A3 landscape  (16.53,11.69)
    #  A4 portrait  (8.27, 11.69)  
    page_dims = (8.27, 11.69)
    with pdfpages(outfile) as pdfpages:
        axlist = []
        fig,axes = plt.subplots(nrows=2, ncols=1, figsize=page_dims, layout='constrained') 
        fig.suptitle(f'Read counts frequency plots.')
        for a in axes.flat:
            axlist.append(a)
        ax = axlist[0]
        counts_axis_plot_sns(ax, df, scale=None)
        ax = axlist[1]
        counts_axis_plot_sns(ax, df, scale='log10')
        
        pdfpages.savefig(fig)


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
        title = f'{title} (scale)'

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




def make_shoulder_plot_sns(df, 
                           title='frequency plot', 
                           site='target', 
                           outfile='frequency-plot.pdf', 
                           column='read_count' ):
    '''
    makes two figures, one log-normalized, one straight for same counts df. 
    '''
    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    
    cdf = pd.DataFrame( df[df['site'] == site ][column].copy(), columns=[column])
    cdf.sort_values(by=column, ascending=False, inplace=True)
    cdf.reset_index(inplace=True,  drop=True)
    cdf['index'] = cdf.index.astype('int64')

    #page_dims = 
    #  A4 landscape   (11.69,8.27) 
    #  A3 landscape  (16.53,11.69)
    #  A4 portrait  (8.27, 11.69)  
    page_dims = (8.27, 11.69)
    with pdfpages(outfile) as pdfpages:
        axlist = []
        fig,axes = plt.subplots(nrows=2, ncols=1, figsize=page_dims, layout='constrained') 
        fig.suptitle(f'{column} freq: site={site}')
        for a in axes.flat:
            axlist.append(a)
        ax = axlist[0]
        counts_axis_plot_sns(ax, cdf, scale=None, column=column)
        ax = axlist[1]
        counts_axis_plot_sns(ax, cdf, scale='log10', column=column)        
        pdfpages.savefig(fig)


        
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
        lx = 0.2
        ly = 0.2
        ax.text(lx, ly, f"n={n}\ntop={t}\nsum={s}\n", fontsize=11) #add text
        title = f'{title} log10().'
        logging.debug(f'made axis with log10 scale.')       

    ax.set_title(title, fontsize=10)

def make_shoulder_plots(df, outdir=None, cp=None):
    # make shoulder plots. injection, target
    logging.info('making shoulder plots...')
    if outdir is None:
        outdir = os.path.abspath('./')
    logging.getLogger('matplotlib.font_manager').disabled = True
    if len(df[df['site'] == 'injection'] ) > 1:
        make_shoulder_plot_sns(df, site='injection', outfile=f'{outdir}/inj-counts.pdf')
    else:
        logging.info(f'no injection sites, so no plot.')
    if len(df[df['site'] == 'target'] ) > 1:    
        make_shoulder_plot_sns(df, site='target', outfile=f'{outdir}/target-counts.pdf')   
    else:
        logging.info(f'no target sites, so no plot.')
   