import scanpy as sc
import pandas as pd
import muon as mu
from muon import atac as ac
import numpy as np


def read_data(path):

	mdata = mu.read_10x_h5(path)
	return mdata


def atac_qc_metrics(atac):

	sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)

	atac.obs.rename(
		columns={
		'n_genes_by_counts': 'n_features_per_cell',
		'total_counts': 'total_fragment_counts',
		},
		inplace=True
	)


	atac.obs['log_total_fragment_counts'] = np.log10(atac.obs['total_fragment_counts'])


