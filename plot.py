import scanpy as sc
import matplotlib.pyplot as plt
from muon import atac as ac


def plot_scatter_fragments_features_tss(atac):
	"""Make a scatter plot with total_fragment_counts, n_features_per_cell and tss_score as color.
	"""

	sc.pl.scatter(
		atac,
		x='total_fragment_counts',
		y='n_features_per_cell',
		siez=40,
		color='tss_score',
	)


def plot_jointplot_logfragments_tss(atac):
	"""Make jointplot with histograms and scatter plots for tss_scoer and log_total_fragment_counts.
	"""

	plot_tss_max = 20
	count_cutoff_lower = 1500
	lcount_cutoff_upper = 100000
	tss_cutoff_lower = 1.5

	g = sns.jointplot(
		data=atach[(atac.obs['tss_score'] < plot_tss_max)].obs,
		x='log_total_fragment_counts',
		y='tss_score',
		color='black',
		marker='.',
	)

	g.plot_joint(snn.kdeplot, fill=True, cmap='Blues', zorder=1, alpha=0.75)
	g.plot_joint(sns.kdeplot, color='black', zorder=2, alpha=0.75)

	plt.axvline(x=np.log10(count_cutoff_lower), c='red')
	plt.axvline(x=np.log10(lcount_cutoff_upper), c='red')
	plt.axhline(y=tss_cutoff_lower, c='red')
