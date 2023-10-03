# This is a sample Python script.
# to QC and graph heatmap for Mapseq data

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy
import os
import openpyxl
import sys
sys.setrecursionlimit(100000)

from configparser import ConfigParser

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *



def filter_counts(countsdata, inject_sites_bc, neg_controls_bc, tech_controls_bc, project_sites_bc, min_inj=5, min_proj_sum=0, max_proj_sites=1, max_controls=2):
    # inject sites > min_inj
    brain_data = countsdata[inject_sites_bc + neg_controls_bc + tech_controls_bc + project_sites_bc]
    brain_data = brain_data[brain_data[inject_sites_bc[0]] > min_inj]
    # remove the barcodes have 0 counts in all project sites
    brain_data = brain_data[brain_data[project_sites_bc].sum(1) > min_proj_sum]
    # max(project_sites) > 1 and max(control) < 2
    brain_data = brain_data[brain_data[project_sites_bc].max(1) > max_proj_sites]
    brain_data = brain_data[brain_data[tech_controls_bc + neg_controls_bc].max(1) < max_controls]
    brain_data_proj = brain_data[neg_controls_bc + project_sites_bc]
    return brain_data, brain_data_proj


def spike_norm(spikedata, inject_sites_bc, neg_controls_bc, tech_controls_bc, project_sites_bc):
    # spike counts data
    brain_spikedata = spikedata[inject_sites_bc + neg_controls_bc + tech_controls_bc + project_sites_bc]
    norm_factor = brain_spikedata.sum()
    norm_factor_proj = norm_factor[neg_controls_bc + project_sites_bc]
    return norm_factor, norm_factor_proj


def graph_heatmap(outputdir, counts_data, sampleID='sample', camp='Reds'):
    # plotting heatmap of counts in each BC
    # counts
    # assign 0 to Nan
    counts_data = counts_data.fillna(0)
    g = sns.clustermap(counts_data, cmap=camp, yticklabels=False, col_cluster=False, standard_scale=0)
    g.fig.subplots_adjust(right=0.7)
    g.ax_cbar.set_position((0.8, .2, .03, .4))
    plt.title(sampleID+'\nCounts')
    plt.savefig(outputdir+sampleID+'.count.pdf')


def main():
# read data
    countsdata = pd.read_csv(workingdir+'M230.all.tsv', sep='\t')
    spikedata = pd.read_csv(workingdir+'M230.all.spike.tsv', sep='\t')
    # Load the xlsx file
    excel_data = pd.read_excel('FinalM230sampleinformaitonlist_SZ_HZwRTprimerInfo.xlsx')
    colnames = excel_data.iloc[0]
    excel_data = excel_data.iloc[1:excel_data.shape[0]]
    excel_data.columns = colnames
    # Read the values of the file in the dataframe
    data = pd.DataFrame(excel_data, columns=['Brain', 'RT primers for MAPseq', 'Region name'])
    data['RT primers for MAPseq'] = 'BC' + data['RT primers for MAPseq'].astype(str)
    tech_controls = data['RT primers for MAPseq'][data['Region name'].str.contains('H2O control')]
    tech_controls = tech_controls.values.tolist()
    allnormmatrix = pd.DataFrame()

    for brainId in data['Brain'].dropna().unique():
        # find barcodes
        sampledata = data[data['Brain'] == brainId]
        inject_sites = sampledata['RT primers for MAPseq'][sampledata['Region name']=='Inj']
        inject_sites = inject_sites.values.tolist()
        neg_controls = sampledata['RT primers for MAPseq'][sampledata['Region name']=='NG']
        neg_controls = neg_controls.values.tolist()
        project_sites = sampledata['RT primers for MAPseq']
        project_sites = project_sites.values.tolist()
        project_sites = [ele for ele in project_sites if ele not in inject_sites+neg_controls]
        # filter data
        brain_data, brain_data_proj = filter_counts(countsdata, inject_sites, neg_controls, tech_controls, project_sites, min_inj=5, min_proj_sum=0, max_proj_sites=1, max_controls=2)
        regionnames = []
        for bc in brain_data_proj.columns:
            regionnames.append(sampledata['Region name'][sampledata['RT primers for MAPseq']==bc].to_string(index=False))
        brain_data_proj.columns = regionnames
        # norm factor
        norm_factor, norm_factor_proj = spike_norm(spikedata, inject_sites, neg_controls, tech_controls, project_sites)
        norm_factor_proj.index = regionnames
        brain_data_norm = brain_data_proj/norm_factor_proj
        # merge
        allnormmatrix = pd.concat([allnormmatrix, brain_data_proj/norm_factor_proj ], ignore_index=True, axis=0)
        # graph
        graph_heatmap(workingdir, brain_data_proj, sampleID='Brain'+str(brainId))
        graph_heatmap(workingdir, brain_data_norm, sampleID='Brain' + str(brainId) + '.norm')

    graph_heatmap(workingdir, allnormmatrix, sampleID='Merged.norm')




if __name__ == "__main__":
    # working dir
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')
    
    parser.add_argument('-c','--config', 
                        metavar='config',
                        required=False,
                        default=os.path.expanduser('~/git/mapseq-processing/etc/mapseq.conf'),
                        type=str, 
                        help='config file.')    
    
    parser.add_argument('-o','--outprefix', 
                    metavar='outprefix',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outfile prefix, e.g. M229, stdout if not given.')  

    parser.add_argument('-O','--outdir', 
                    metavar='outdir',
                    required=False,
                    default=None, 
                    type=str, 
                    help='outdir. input file base dir if not given.')     

    parser.add_argument('-s','--sampleinfo', 
                        metavar='sampleinfo',
                        required=True,
                        default=None,
                        type=str, 
                        help='XLS sampleinfo file. ')

    parser.add_argument('infiles',
                        metavar='infiles',
                        nargs ="+",
                        type=str,
                        help='"all" TSV from process_ssifasta. columns=(sequence, counts, type, label')
       
    args= parser.parse_args()

    


