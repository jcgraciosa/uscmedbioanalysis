#----------------------------------------------------------------
# Class is for analyzing the MSD of the detected particles 
# 
#
#Revision History:
#
#Date			Author		Version		Description
#----------------------------------------------------------------
# 08.22.2018	peaBrane	0.1			- initial implementation
#										 
#                            
#----------------------------------------------------------------

import os, sys
import numpy as np
import pandas as pd
from time import time
import matplotlib.pyplot as plt
import trackpy as tp
#from scipy.stats import linregress

class TrackAnalyzer(object):

	def __init__(self, in_dir, out_dir):
		"""
		Initialization:
		in_dir = directory containing the MSD files
		out_dir = directory where output data are saved
		"""
		self.in_dir = in_dir
		self.out_dir = out_dir


	def read_msd_files(self, min_t = 0, max_t = 2, save_collated = False, show_msd = False):
		"""
		reads the MSD files in in_dir.
		max_t = maximum time interval
		save_collated = set to True if you want a csv file containing the MSD of all particles
		"""
		cnt = 0
		for filename in os.listdir(self.in_dir):
			if filename.endswith('.txt'): # FIXME: finalize the value of this
				data = pd.read_csv(os.path.join(self.in_dir, filename), sep="\t", header = None) # assume sep is tab and header is present (not None)
				if cnt == 0:
					self.msd_data = data.iloc[:,5] #change this if name is consistent
					time_df = data.iloc[:, 0] # or change if name is consistent
				else: # only append the msd data 
				# FIXME: add possibility of changing the time column in
					if data.shape[0] > self.msd_data.shape[0]:
						time_df = data.iloc[:, 0]
						#print('changed', time_df.shape[0], self.msd_data.shape[0])

					self.msd_data = pd.concat([self.msd_data, data.iloc[:, 5]], axis = 1)

				cnt = cnt + 1

		self.msd_data.insert(0, 'time', time_df)

		# remove data for t outside the desired range
		# less than desired range
		time_idx = time_df.iloc[(time_df - min_t).abs().argsort()[:2]]
		idx = time_idx.index.tolist()[0]
		self.msd_data = self.msd_data.truncate(before = idx)

		# after desired range
		time_idx = time_df.iloc[(time_df - max_t).abs().argsort()[:2]]
		idx = time_idx.index.tolist()[0]
		self.msd_data = self.msd_data.truncate(after = idx)
		self.msd_data.set_index('time')
		self.ensemble_msd = self.msd_data.mean(axis = 1)

		if save_collated:
			self.msd_data.to_csv(os.path.join(self.out_dir, 'msd_data.csv'), sep=',')

		if show_msd:
			#print(time_df)
			print(self.msd_data)

		return

	def plot_single_msd(self, save_plot = True, plot_suffix = '1'):

		fig, ax = plt.subplots(dpi = 300)
		ax.plot(self.msd_data.index, self.msd_data, 'k-', alpha = 0.1)  # black lines, semitransparent
		ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
		       xlabel='lag time $t$')
		ax.set_xscale('log')
		ax.set_yscale('log')

		if save_plot:
			fig.savefig(os.path.join(self.out_dir, 'single_msd_' + plot_suffix + '.png'), dpi = 300)

		return

	def plot_ensemble_msd(self, save_plot = True, plot_suffix = '1'):

		fig, ax = plt.subplots(dpi = 300)
		ax.plot(self.ensemble_msd.index, self.ensemble_msd, 'o')
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
		       xlabel='lag time $t$')

		if save_plot:
			fig.savefig(os.path.join(self.out_dir, 'single_msd_' + plot_suffix + '.png'), dpi = 300)

		return

	def get_diffusivity_ensemble(self):

		lag_t = np.asarray(self.ensemble_msd.index)
		msd = np.asarray(self.ensemble_msd)

		self.ensemble_diff, self.ensemble_b = np.polyfit(lag_t, msd, 1)
		self.ensemble_diff = np.round(self.ensemble_diff, 6)
		self.ensemble_b = np.round(self.ensemble_b, 6)

		print('Ensemble diffusivity: ', self.ensemble_diff)
		print('Ensemble intercept: ', self.ensemble_b)

		return

	def get_diffusivity_single(self, min_diff = 0, max_diff = 2, save_param = False, filename = 'diffusivity'):
		"""
		fits MSD and t -> input t and y are already in log
		save_collated = set to True if you want a csv file containing the MSD of all particles
		"""
		num_tracks = self.msd_data.shape[1] - 1
		params = np.zeros((num_tracks, 2))

		for i in range(num_tracks):
			df = self.msd_data.iloc[:, [0, i + 1]]
			df = df.dropna()
			t = np.asarray(df.iloc[:, 0])
			y = np.asarray(df.iloc[:, 1])
			#out_param = linregress(t, y)
			#params[i, 0] = out_param[0]
			#params[i, 1] = out_param[1]
			params[i, 0], params[i, 1] = np.polyfit(t, y, 1)
		
		params = np.round(params, 6)
		self.params = pd.DataFrame(data = params, columns = ['diffusivity', 'b']) 
		self.params = self.params[self.params.diffusivity < max_diff]
		self.params = self.params[self.params.diffusivity > min_diff]
		
		self.diff_mean = np.mean(self.params[:, 0]) 
		self.diff_std = np.std(self.params[:, 0])

		self.b_mean = np.mean(self.params[:, 1]) 
		self.b_std = np.std(self.params[:, 1])

		print('Average diffusivity: ', self.diff_mean)
		print('Standard deviation of diffusivity: ', self.diff_std)
		if save_param:
			self.params.to_csv(os.path.join(self.out_dir, filename + '.csv'), sep = '\t')	
		
		return

	def plot_diffusivity_dist(self, save_plot = True, plot_suffix = '1'):

		plt.rc('font', family = 'serif')
		fig = plt.figure(1, dpi = 300)

		diffusivity = np.asarray(self.params.iloc[:, 0])
		constant = np.asarray(self.params.iloc[:, 1])

		# plot histogram of diffusivity values
		self.diff_hist, self.diff_edges = np.histogram(diffusivity, density = False, bins = 'fd')
		#self.diff_hist, self.diff_edges = np.histogram(diffusivity, density = False)
		self.diff_edges = self.diff_edges[:-1]
		self.diff_edges = self.diff_edges[self.diff_hist > 0]
		self.diff_hist = self.diff_hist[self.diff_hist > 0]
		width = (self.diff_edges[1] - self.diff_edges[0])
		plt.bar(self.diff_edges, self.diff_hist, width = width)
		#plt.plot(self.diff_edges, self.diff_hist, 'o-', linewidth = 1, markersize = 2, label = 'Diffusivity')
		plt.legend()
		plt.xlabel('Diffusivity')
		plt.ylabel('Counts')
		plt.xlim(left = self.diff_edges.min(), right = self.diff_edges.max())
		plt.ylim(bottom = 0, top = self.diff_hist.max())
		fig.savefig(os.path.join(self.out_dir, 'diffusivity_' + plot_suffix + '.png'), dpi = 300)


		fig = plt.figure(2, dpi = 300)
		# plot histogram of constant
		#self.const_hist, self.const_edges = np.histogram(constant, density = False, bins = 'fd')
		self.const_hist, self.const_edges = np.histogram(constant, density = False)
		self.const_edges = self.const_edges[:-1]
		self.const_edges = self.const_edges[self.const_hist > 0]
		self.const_hist = self.const_hist[self.const_hist > 0]
		width = (self.const_edges[1] - self.const_edges[0])
		plt.bar(self.const_edges, self.const_hist)
		#plt.plot(self.const_edges, self.const_hist, 'o-', linewidth = 1, markersize = 2, label = 'Constant')
		plt.legend()
		plt.xlabel('Constant')
		plt.ylabel('Counts')
		plt.xlim(left = self.const_edges.min(), right = self.const_edges.max())
		plt.ylim(bottom = 0, top = self.const_hist.max())
		if save_plot:
			fig.savefig(os.path.join(self.out_dir, 'constant_' + plot_suffix + '.png'), dpi = 300)
		
		return 
