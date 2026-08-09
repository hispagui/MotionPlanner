"""
pose_interpolation.py -- practical demo of the so3 / se3 methods
    A rigid body (a coordinate frame) travels from start pose T0 to end pose T1.
    2 different paths, screw_interpolate and decoupled_interpolate
    Produces pose_interpolation.gi
"""

import os
import sys
import numpy as np
import so3
import se3
import matplotlib
matplotlib.use("Agg") # saves to file
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D      # noqa: F401  (registers 3d proj)


# definition of poses, measurments and precomputations

T0 = se3.make(np.eye(3), [0.0, 0.0, 0.0]) # first pose (identity)
T1 = se3.make(so3.exp(np.array([0.3, 1.0, 0.5]) / np.linalg.norm([0.3, 1.0, 0.5]) * 2.4), [3.0, 1.5, 2.0]) # kind of a random pose
# check geod.distance between bodies and position of body at half time 0.5
print("geodesic distance T0 -> T1 (length_scale = 1.0) :", round(se3.geodesic_distance(T0, T1, 1.0), 4))
print("pose at t=0.5 (screw), translation part:", np.round(se3.split(se3.screw_interpolate(T0, T1, 0.5))[1], 4))

N = 60 # number of frames in animation
ts = np.linspace(0.0, 1.0, N)
# precompute and store every pose along both paths
screw = [se3.screw_interpolate(T0, T1, t) for t in ts]
decpl = [se3.decoupled_interpolate(T0, T1, t) for t in ts]
screw_pos = np.array([se3.split(T)[1] for T in screw])  
decpl_pos = np.array([se3.split(T)[1] for T in decpl])    



# below is all the plotting
# set up the 3D space
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection="3d")

allpts = np.vstack([screw_pos, decpl_pos])
lo, hi = allpts.min(0) - 1.2, allpts.max(0) + 1.2
ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
ax.set_box_aspect((hi - lo))
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
ax.view_init(elev=22, azim=-60)

AXIS_LEN = 0.9
AXIS_COLORS = ("#D14520", "#3B8BD4", "#1D9E75")   # body x, y, z

def draw_triad(style, alpha):
    # draws 3 Line3D for a body frame's (x,y,z axes)
    return [ax.plot([], [], [], style, color=c, lw=2.2, alpha=alpha)[0] for c in AXIS_COLORS]

screw_triad = draw_triad("-", 1.0)     # solid  = screw motion
decpl_triad = draw_triad("--", 0.65)   # dashed = decoupled motion

# faint full reference paths (so the arc vs chord shapes read even at rest)
ax.plot(screw_pos[:, 0], screw_pos[:, 1], screw_pos[:, 2],
        color="#D85A30", lw=1, alpha=0.25)
ax.plot(decpl_pos[:, 0], decpl_pos[:, 1], decpl_pos[:, 2],
        color="#1D9E75", lw=1, alpha=0.25, ls="--")

# trails from origin
screw_trail, = ax.plot([], [], [], color="#D85A30", lw=2.5, label = "screw (geodesic)")
decpl_trail, = ax.plot([], [], [], color="#1D9E75", lw=2.5, ls="--", label = "decoupled (slerp)")

# endpoint markers
ax.scatter(*screw_pos[0], color="k", s=40)
ax.scatter(*screw_pos[-1], color="k", s=40)
ax.text(*screw_pos[0], "  T0", fontsize=11)
ax.text(*screw_pos[-1], "  T1", fontsize=11)
ax.legend(loc="upper left")
ax.set_title("Moving a rigid body from T0 to T1: screw vs decoupled")

def set_triad(handles, T):
    R, p = se3.split(T)
    for i, h in enumerate(handles):
        tip = p + AXIS_LEN * R[:, i]
        h.set_data([p[0], tip[0]], [p[1], tip[1]])
        h.set_3d_properties([p[2], tip[2]])

def update(frame):
    set_triad(screw_triad, screw[frame])
    set_triad(decpl_triad, decpl[frame])
    screw_trail.set_data(screw_pos[:frame + 1, 0], screw_pos[:frame + 1, 1])
    screw_trail.set_3d_properties(screw_pos[:frame + 1, 2])
    decpl_trail.set_data(decpl_pos[:frame + 1, 0], decpl_pos[:frame + 1, 1])
    decpl_trail.set_3d_properties(decpl_pos[:frame + 1, 2])
    ax.view_init(elev=22, azim=-60 + 0.35 * frame)   # slow orbit for depth
    return (*screw_triad, *decpl_triad, screw_trail, decpl_trail)

# animate
anim = FuncAnimation(fig, update, frames=N, interval=55, blit=False)
out = os.path.join(os.path.dirname(__file__), "pose_interpolation.gif")
anim.save(out, writer=PillowWriter(fps=20))
print("wrote in", out)