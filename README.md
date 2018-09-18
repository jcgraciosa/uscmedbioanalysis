# USC Medical Biophysics Tools for Analysis
## Description
This package contains a set of tools for analyzing Medical-Biophysics experimental data.

## Libraries/ Requirements
- Python 3.6
- Windows 10
- numpy
- trackpy

## Installation Guide
1. `pip3 install -r requirements.txt`
2. `cd clone_dir`
3. `pip3 install ./uscmedbioanalysis`

## Main Functions

1. read_msd_files
```
Reads the particle MSD data present in the input directory.
min_t - minimum time to include in reading (default = 0)
max_t - maximum time to include in the reading (default = 2)
save_collated - set to True if you want to save the collated MSD into a csv file (default = False)
show_msd - set to True if you want to print the collated MSD after reading (default = False)
```

2. plot_single_msd
```
For visualizing the individual MSD (use this plot to decide the time interval of interest)
save_plot - set to True if you want to save the generated plot (default = False)
plot_suffix - file name suffix and is only needed when you want to save the plot (default = '1')
```

3. plot_ensemble_msd
```
For visualizing the ensemble MSD (use this plot to decide the time interval of interest)
save_plot - set to True if you want to save the generated plot (default = False)
plot_suffix - file name suffix and is only needed when you want to save the plot (default = '1')
```

4. get_diffusivity_single
```
perform fitting to obtain the diffusivity for each particle
min_diff - minimum diffusion value to include in the final parameter list to remove spurious tracks (default = 0)
max_diff - maximum diffusion value to include in the final parameter list to remove spurious tracks (default = 2)
save_param - set to True if you want to save the fitting parameters (diffusion and intercept) to a csv file (default = False)
filename - filename of csv file containing parameters (default - 'diffusivity')
```

5. get_diffusivity_ensemble
```
 perform fitting to obtain the diffusivity for the ensemble

```

6. plot_diffusivity_dist
```
save_plot - set to True if you want to save the generated plot (default = False)
plot_suffix - file name suffix and is only needed when you want to save the plot (default = '1')
```
						
## Scripts and Notebooks

1. uscmedbioanalysis_demo 
```
Demonstration of the package's methods.
```

## Sample Usage

1. Reading track files then analyzing the individual properties
```python
import uscmedbioanalysis as medbio

indir = 'input/directory/here'
outdir = 'output/directory/here'

analyzer = medbio.TrackAnalyzer(indir, outdir)
analyzer.read_msd_files(save_collated = True, show_msd = False)
analyzer.plot_single_msd(save_plot = True, plot_suffix = '1')
analyzer.get_diffusivity_single(save_param = True)
analyzer.plot_diffusivity_dist()
```
2. Reading track files then analyzing the ensemble properties
```python
import uscmedbioanalysis as medbio

indir = 'input/directory/here'
outdir = 'output/directory/here'

analyzer = medbio.TrackAnalyzer(indir, outdir)
analyzer.read_msd_files(save_collated = True, show_msd = False)
analyzer.plot_ensemble_msd(save_plot = True, plot_suffix = '1')
analyzer.get_diffusivity_ensemble()
```

## Test data
Sample data in `./test/test_track`. To test, please follow entries in sample usage.


