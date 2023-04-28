import matplotlib.pyplot as plt
import os.path
import pandas as pd

metrics_path = "clipper_metrics_jup.csv"
plt_metrics_path = ""
metrics_time_st = .994e9
metrics_time_ed = .999e9

metrics_df = pd.read_csv(metrics_path)

if metrics_time_st:
    metrics_df = metrics_df.loc[metrics_df['timestep'] >= metrics_time_st]
if metrics_time_ed:
    metrics_df = metrics_df.loc[metrics_df['timestep'] <= metrics_time_ed]

plt.figure()
plt.plot(metrics_df['timestep'], metrics_df['Clipper-Europa distance'], color=(1.0, 1.0, 0))
plt.xlabel("clock time (s)")
plt.ylabel("Clipper-Europa distance (km)")
plt.title("Clipper-Europa Distance as a Function of Time")
plt.savefig(os.path.join(plt_metrics_path, "clipper_europa_dist_custom.png"))

plt.figure()
plt.plot(metrics_df['timestep'], metrics_df['velocity'], color=(0.016, 0.85, 1.0))
plt.xlabel("clock time (s)")
plt.ylabel("Clipper velocity (km/s)")
plt.title("Clipper Velocity as a Function of Time")
plt.savefig(os.path.join(plt_metrics_path, "clipper_vel_custom.png"))

plt.figure()
plt.plot(metrics_df['timestep'], metrics_df['acceleration'], color=(1, 0.063, 0.96))
plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 2])
plt.xlabel("clock time (s)")
plt.ylabel("Clipper acceleration (km/s$^2$)")
plt.title("Clipper Acceleration as a Function of Time")
plt.savefig(os.path.join(plt_metrics_path, "clipper_acc_custom.png"))

fig, ax1 = plt.subplots()
l1, = ax1.plot(metrics_df['timestep'], metrics_df['acceleration'], color=(1, 0.063, 0.96))
ax1.ticklabel_format(axis='y', style='sci', scilimits=[-2, 2])
ax1.set_xlabel("clock time (s)")
ax1.set_ylabel("unitless for comparison (scaled to acceleration)")
ax1.set_title("Combined Metric Plots as a Function of Time")
ax2 = ax1.twinx()
l2, = ax2.plot(metrics_df['timestep'], metrics_df['Clipper-Europa distance'], color=(1.0, 1.0, 0))
ax2.set_yticks([])
ax3 = ax1.twinx()
l3, = ax3.plot(metrics_df['timestep'], metrics_df['velocity'], color=(0.016, 0.85, 1.0))
ax3.set_yticks([])
box1 = ax1.get_position()
ax1.set_position([box1.x0, box1.y0, box1.width * 0.9, box1.height])
ax1.legend([l2, l3, l1], ["dist", "vel", "acc"], loc='center left', bbox_to_anchor=(1, 0.5))
fig.savefig(os.path.join(plt_metrics_path, "clipper_combined_custom.png"))
