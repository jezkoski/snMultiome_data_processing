import scanpy as sc
import pandas as pd
import muon as mu
from muon import atac as ac
import numpy as np


def read_data(path):
	"""Read data to Muon object.

	"""

	mdata = mu.read_10x_h5(path)
	return mdata


def atac_qc_metrics(atac):
	"""Calcualte QC metrics for ATAC data.

	Rename columns 'n_genes_by_counts' -> 'n_features_per_cell' 
	and 'total_counts' -> 'total_fragment_counts'

	"""

	sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)

	atac.obs.rename(
		columns={
		'n_genes_by_counts': 'n_features_per_cell',
		'total_counts': 'total_fragment_counts',
		},
		inplace=True
	)

	atac.obs['log_total_fragment_counts'] = np.log10(atac.obs['total_fragment_counts'])




def calculate_nucleosome_signal(atac, fragment_file_path):
	"""The ratio of mono-nucleosome cut fragments to nucleosome-free fragments is called nucleosome signal.

	Nucleosome binding patterns show usually enrichment around values corresponding to
	fragments bound to a single nucleosome as well as nucleosome-free fragments.

	Args:
		atac: Muon ATAC object
		fragment_file_path: path to the fragments.tsv file, required for this analysis
	"""

	mu.atac.tl.locate_fragments(atac, fragment_file_path)

	ac.pl.fragment_histogram(atac, region='chr1:1-2000000', save='fragment_length_distribution.png')

	ac.tl.nucleosome_signal(atac, n=10e3 * atac.n_obs)

	mu.pl.histogram(atac, 'nucleosome_signal', linewidth=0, save='nucleosome_signal.png')

	
def label_pass_and_fail_NS(atac, threshold):
	"""Label cells that have good quality nucleosome signal (pass or fail the threshold).
	"""

	atac.obs['nuc_signal_filter'] = [
		'NS_FAIL' if ns > threshold else 'NS_PASS'
		for ns in atac.obs['nucleosome_signal']
	]

	atac.obs['nuc_signal_filter'].value_counts()

def plot_PASS_and_FAIL_NS_fragment_size_dist(atac):

	p1 = ac.pl.fragment_histogram(
		atac[atac.obs['nuc_signal_filter'] == 'NS_PASS'], region='chr1:1-2000000'
	)
	p2 = ac.pl.fragment_histogram(
		atac[atac.obs['nuc_signal_filter'] =='NS_FAIL'], region='ch1:1-2000000'
	)


def plot_TSS_enrichment_score_distribution(atac):

	fig, axs = plt.subplots(1,2, figsize=(7,3.5))

	p1 = sns.histplot(atac.obs, x='tss_score', ax=axs[0])
	p1.set_title('Full range')

	p2 = sns.histplot(
		atac.obs,
		x='tss_score',
		binrange(0, atac.obs['tss_score'].quantile(0.995)),
		ax=axs[1],
	)
	p2.set_title('Up to 99.5% percentile')

	plt.suptitle('Distribution of the TSS score')
	plt.tight_layout()

	fig.savefig('TSS_score_distribution.png', dpi=300)



def calculate_TSS_enrichment(atac, tss_regions_path):
	"""Calculate the enrichment of fragments around transcription start sites.

	This assess the signal-to-noise ratio in each cell. We expect chromatin accessibility 
	around TSSs compared to accessibility of flanking regions. 

	Args:
		atac: MUON ATAC seq object
		tss_regions_path: path to file containing the genomic locations of gene TSSs
	"""

	regions = pd.read_csv(tss_regions_path, sep='\t')
	regions.rename(columns={'hg38.refGene.chrom': 'Chromosome', 'hg38.refGene.txStart':'Start', 'hg38.refGene.txEnd':'End'})

	tss = ac.tl.tss_enrichment(mdata, n_tss=3000, features=regions)

	label_pass_and_fail_NS(atac, threshold=2)

	ac.pl.tss_enrichment(tss)

