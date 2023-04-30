# Space-Missions

Dante Basile and Dalia Bekdache

# Environment Specification

## General (History) Specification

A general specification of the environment obtained from running:

```
conda env export --from-history > environment_history.yml
```

can be found in `environment_history.yml`.

## Full Specification for win64 System

A full specification of the environment for a win64 system obtained from running:

```
conda env export > environment.yml
```

can be found in `environment_win64.yml`.

## Building Environment

Running:

```
conda env create -f <environment.yml>
```

with the appropriate file for the system platform will build the environment.

# `clipper3.py` Instructions

## Description

Improved version of the original script for visualizing the Europa Clipper mission along with the relevant bodies of the solar system

## Summary of Additions

### New Data

The latest mission trajectory file has been provided at `./naif/21F31_MEGA_L241010_A300411_LP01_V5_pad_scpse.bsp`. This is used as a replacement for the previous `./naif/19F23_VEEGA_L230511_A290930_LP01_V2_scpse.bsp` Hardcoded parameters have been adjusted for:

1. Loading the new file instead of the old one
2. Empirically determined new launch time
3. Empirically determined new arrival time

### Event Viewer

A new GUI element was added that allows the user to select key mission events for viewing. After a selection is made, all relevant simulation parameters are loaded to allow the user to easily view the event.

### Metrics and GUI Elements

Several relevant metrics (Clipper-Europa distance, Clipper velocity, and Clipper acceleration) are now computed and displayed on the HUD (heads up display). The HUD also shows a status indicating whether Clipper has launched.

### Plotting and Saving Metrics

The metrics can be saved as PNG plots and/or a CSV file after the simulation terminates.

## Usage

This script can be run as before:

```
python clipper3.py
```

Some additional command line flags have been added. `--plt_metrics` will cause the simulation to save plots of the metrics after the simulation terminates. The files are called `clipper_europa_dist.png`, `clipper_vel.png`, and `clipper_acc.png`. These plots are saved in the same directory as `clipper3.py` unless an optional directory path is supplied after the flag. In this case, they will be saved in the chosen directory. `--data_metrics` will cause the simulation to save the metrics in CSV format after the simulation terminates. For example, to save both the plots and CSV in the default location, the script would be run as:

```
python clipper3.py --plt_metrics --data_metrics
```

# `plot_metrics_custom.py` Instructions

## Description

This script allows custom plots to be made based on the metrics collected from a `clipper3.py` run. Specifically, it allows the user to select a range of the total timeframe of a run to plot. This was primarily used to obtain clear plots of Clipper's orbit of Jupiter, as the timescale of such an orbit is much briefer than the simulation as a whole. The metrics are accessed from a CSV file that can be produced by `clipper3.py`. The plots produced are otherwise similar to those that can be produced by running:

```
python clipper3.py --plt_metrics
```

## Usage

The plots taken as input can be produced from `clipper3.py` by running:

```
python clipper3.py --data_metrics
```

The user then edits the file directly to specify the hardcoded parameters:

1. `metrics_path`: a path to the CSV file storing the metrics to plot

2. `plt_metrics_path`: a path to the folder where the plots will be saved to. If given an empty string, the plots will be saved in the same directory as `plot_metrics_custom.py`.

For example, the parameters used to obtain the plots of Clipper's orbit of Jupiter are:

``` python
metrics_path = "clipper_metrics_jup.csv"
plt_metrics_path = ""
metrics_time_st = .994e9
metrics_time_ed = .999e9
```

The script can then be run with:

```
python plot_metrics_custom.py
```

# Generated Data in `outputs`

The folder `outputs` has been added to the repository. This folder contains some plots and CSV files that were produced from `clipper3.py` and `plot_metrics_custom.py`.
