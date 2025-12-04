
import logging
import os
import sys

import pandas as pd
import pyarrow as pa
import seaborn as sns
import kneed
import matplotlib
import matplotlib.pyplot as plt

from mapseq.core import *
 

QPCR_COLNAME_MAP = {
     'Well'          :  'well',
     'Well Position' :  'well_pos',
     'Cycle'         :  'cycle',
     'Target Name'   :  'target',              
     'Rn'            :  'rn',        
     'Delta Rn'      :  'rn_delta'
    }

QPCR_SAMPLE_SHEET = 'Sample Setup'
QPCR_AMP_SHEET = 'Amplification Data'

QPCR_OUT_COLUMNS = [ 'well_id','low','high','mean','error']
QPCR_GROUP_COLUMNS = [ 'group_id','well_id','mean']



def load_amp_sheet(infile, cp=None):
    amp_sheet = QPCR_AMP_SHEET
    edf = pd.read_excel(infile, 
                        sheet_name = amp_sheet, 
                        header=None, 
                        dtype=str) 
    edf = edf[~edf.isna().all(axis=1)]    
    df = edf[ ~edf.iloc[:,3].isna()]
    columns = list( df.iloc[0,:] )
    df = df.iloc[1:,:]
    df.columns = columns
    df = df.rename(columns = QPCR_COLNAME_MAP)
    df.reset_index(inplace=True, drop=True)
    df['cycle'] = df['cycle'].astype(int)
    df['rn'] = df['rn'].astype(float)
    df['rn_delta'] = df['rn_delta'].astype(float)
    return df

def make_groups_df(df, cp=None):
    '''
    @arg df    cycle df. columns=[ 'well_id','low','high','mean','error']
    
    '''
    if cp is None:
        cp = get_default_config()
    
    group_span = cp.getint('qpcr','group_span', fallback=3)
    
    n_init = len(df)
    df = df[df['error'] != True ]
    df.reset_index(inplace=True, drop=True)
    n_noerror = len(df)
    n_error = n_init - n_noerror
    if n_error > 0:
        logging.info(f'removed {n_error} entries.')
    
    sdf = df.sort_values(by='mean')
    sdf = sdf[['well_id','mean']]
    slist = sdf.values.tolist()
    # Assign groups. 
    group_id = 1
    floor = None
    
    # LOL for output df. 
    grouped_list = []
    
    for well_id, mean in slist:
        print(f'{well_id} {mean}')
        if floor is None:
            floor = mean
            grouped_list.append( [f'group{group_id}', well_id, mean ])
        elif mean <= (floor + group_span):
            grouped_list.append( [ f'group{group_id}', well_id, mean])
        else:
            floor = mean
            group_id += 1
            grouped_list.append( [f'group{group_id}', well_id, mean ])
    
    gdf = pd.DataFrame(grouped_list, columns=QPCR_GROUP_COLUMNS)    
    
    gmdf = gdf.groupby('group_id').agg( {'mean':'mean'})
    gmdf['cycles_rec'] = gmdf['mean'].apply(np.ceil)
    gmdf['cycles_rec'] = gmdf['cycles_rec'].astype(int)
    gmdf.drop(['mean'], axis=1, inplace=True)
    gdf = pd.merge(gdf, gmdf, on='group_id')
    return gdf


def pad_cycles(df, column='cycle', prepend=3):
    '''
    make additional rows below lowest value of column, but 
    with all other values the same. 
    
    '''
    min_cval = df[column].min()
    logging.debug(f'minimum {column} value = {min_cval}')
    mdf = df[df['cycle'] == min_cval ]
    mrow = mdf.iloc[0].to_dict()
    newrows = []
    
    for i in range(1, prepend):
        logging.debug(f'handling prepending i={i}')
        nrow = mrow.copy()
        nrow[column] = min_cval - i
        newrows.append(nrow)
    ndf = pd.DataFrame(newrows)
    logging.debug(f'made padding DF:\n{ndf}')
    df = pd.concat( [df, ndf], copy=False, ignore_index=True )
    df.sort_values(by=column, inplace=True)
    logging.debug(f'new padded DF:\{df}')
    return df


def qpcr_check_wells(infile, 
                     outfile,
                     sensitivity, 
                     polynomial, 
                     cp = None):
    '''
    parse qpcr report XLS, assess optimal cycles, output report XLSX.
    
    for each well, determine which cycle numbers correspond
    to 0.1 -> 2.0 (or plateau) for rn_delta. 
    
    Rn = Reporter
    
    '''
    logging.info(f'QPCR check wells. infile={infile} outfile={outfile}')

    if cp is None:
        cp = get_default_config()

    from matplotlib.backends.backend_pdf import PdfPages as pdfpages
    from matplotlib.ticker import FixedLocator, MaxNLocator

    # get values to create other files in outdir. 
    outfile = os.path.abspath(outfile)
    filepath = os.path.abspath(outfile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (base, ext) = os.path.splitext(filename)   
    head = base.split('.')[0]
    outdir = dirname    

    df = load_amp_sheet(infile, cp)    
    wells = list( df['well_pos'].unique() )

    # Get params
    project_id = cp.get('project','project_id')
    if sensitivity is None:
        sensitivity = cp.getfloat('qpcr','kneed_sensitivity', fallback=2.0)
    if polynomial is None:
        polynomial = cp.getint('qpcr','kneed_polynomial', fallback=2)
    
    low_offset = cp.getint('qpcr','cycle_low_offset', fallback = 1)
    high_offset = cp.getint('qpcr','cycle_high_offset', fallback = 0)


    plotfile = f'{outdir}/{project_id}.qpcrwells.pdf'
    os.makedirs(outdir, exist_ok=True)

    page_dims = (11.7, 8.27)
    title = f'{project_id}: cycle estimation.\nS={sensitivity} poly={polynomial}'
    
    outdflist = []
        
    with pdfpages(plotfile) as pdfpages:
        
        # set up plot pages
        plots_per_page = 9
        num_figs = float(len(wells)) / float(plots_per_page)
        if num_figs % 9 == 0:
            num_figs = int(num_figs)
        else:
            num_figs = int(num_figs) + 1
        logging.debug(f'with {plots_per_page} plots/page, need {num_figs} for {len(wells)} wells.')
        
        figlist = []
        axlist = []
        for i in range(0, num_figs):
            fig, axes = plt.subplots(nrows=3, ncols=3, figsize=page_dims,  layout='constrained')
            fig.suptitle(title)
            figlist.append(fig)
            # numpy.flatirator doesn't handle indexing
            for a in axes.flat:
                axlist.append(a)
        logging.debug(f'created {len(wells)} figures to go on {num_figs} pages. ')        
        
        
        for i, well_id in enumerate(wells):
            logging.debug(f'handling well_id = {well_id} idx={i}')
            wdf =  df[df['well_pos'] == well_id ]
            wdf = pad_cycles(wdf)
            x = list( wdf['cycle'].astype(int) )
            y = list( wdf['rn'].astype(float) )
            min_x = min(x)
            max_x = max(x)
            min_y = min(y)
            max_y = max(y)
            logging.debug(f'[{well_id}] Calling kneed. S={sensitivity} poly={polynomial} x_range=[{min_x},{max_x}] y_range=[{min_y},{max_y}]')
            kneedle = kneed.KneeLocator(x, y, 
                                  S=sensitivity, 
                                  curve="convex", 
                                  online=True, 
                                  direction="increasing", 
                                  polynomial_degree=polynomial)
            floor = kneedle.knee
            logging.debug(f'[{well_id}]: floor={floor}')
            if floor is None:
                floor = min_x
            else:
                floor = floor + low_offset
               
            kneedle = kneed.KneeLocator(x, y, 
                                  S=sensitivity, 
                                  curve="concave", 
                                  online=True, 
                                  direction="increasing", 
                                  polynomial_degree=polynomial)
    
            ceiling = kneedle.knee
            logging.debug(f'[{well_id}]: floor={ceiling}')
            if ceiling is None:
                ceiling = max_x
            else:
                ceiling = ceiling + high_offset

            meanval = (( ceiling + floor ) / 2 )
            
            #Build plot
            sns.lineplot(ax=axlist[i], x=x, y=y)
            axlist[i].axvline(floor)
            axlist[i].axvline(ceiling)
            axlist[i].set_ylabel(f'Rn')
            axlist[i].set_xlabel(f'Cycle')
            axlist[i].xaxis.set_minor_locator(FixedLocator(range(min_x, max_x)))             
            axlist[i].set_title(f'{well_id}: low={floor} high={ceiling} ')        

            logging.debug(f'checking for errors. well_id = {well_id}')
            error_alert = False
            if floor > ceiling:
                error_alert = True
            if floor == min_x + low_offset:
                error_alert = True
            if ceiling == max_x + high_offset:
                error_alert = True
            # Add other checks here.
            
            
            if error_alert:
                xpos = max_x / 2
                ypos =  min_y + (( max_y - min_y ) / 2)
                logging.error(f'error position: max_x={max_x} min_x={min_x} max_y={max_y} min_y={min_y} -> ({xpos},{ypos})')
                axlist[i].text(xpos,
                               ypos, 
                               s=f"ERROR", 
                               fontsize=16,
                               color="Red",             
            )
            logging.info(f'well_id={well_id} range=[{min_x},{max_x}] floor={floor} ceiling={ceiling} ')

            outdflist.append( [ well_id, floor, ceiling, meanval, error_alert])


        for f in figlist:
            pdfpages.savefig(f)
    logging.info(f'Made cycle plot in {plotfile}')
    
    cycle_df = pd.DataFrame(outdflist, columns=QPCR_OUT_COLUMNS)
    logging.debug(f'cycle_df = \n{cycle_df}')

    group_df = make_groups_df(cycle_df)
    
    mdf = pd.merge( cycle_df, group_df, on='well_id',how='outer')
    mdf.drop(['mean_y'], axis=1, inplace=True)
    mdf.rename({'mean_x':'mean'}, inplace=True, axis=1)
    mdf.sort_values(by=['group_id','mean'], inplace=True)

    logging.info(f'writing out XLSX report: {outfile}')
    with pd.ExcelWriter(outfile) as writer:
        mdf.to_excel(writer, sheet_name='Estimated Cycle Range')
        #group_df.to_excel(writer, sheet_name='Cycle Groups')
    return cycle_df

