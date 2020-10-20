""" 
filename: plot_hall_cube.py

The following functions are used for plotting data extracted from the Hall Sensor Cube. Mainly for 3D plots of the vector field.

Author: Jona Buehler 2020

Documentation and Updates by Nicholas Meinhardt, Maxwell Guerne-Kieferndorf (QZabre)
                             nmeinhar@student.ethz.ch, gmaxwell@student.ethz.ch
        
Date: 09.10.2020
"""

########## Standard library imports ##########
import numpy as np
from numpy.linalg import norm
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import to_hex
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os
import sys



def plot_set(SensorData, fig=None, ax=None, Pos=np.array([0, 0, 0]), Vecs=True,
             Mag_label=False, Cont=False, Scat_Mag=False, Show=True,
             title_on=True, cmin=None, cmax=None, single_comp=None, omit_64=False):
    """
    Generate 3d plot of measurement data, where a single (mean) magnetic field vector is provided for each sensor. 

    Args:
    - SensorData (ndarray of size (#sensors, 3)): measured (mean) values for all sensors.
    Note that the 64th senser used to be faulty, hence #sensor can be 63 or 64.
    - Pos (1d-array of length 3): contains the offset position of cube (which is position of sensor 49)
    - fig, ax: Figure and Axis objects of matplotlib. 
    - cmin, cmax (float): min and max of colormap for displaying magnetic field strength. 
    If both any of them is None, the scale is chosen automatically based on the data
    - Vecs (bool): flag to switch on/off plotting of vectors
    - Mag_label (bool): flag to switch on/off dislaying the magnitudes of vectors next 
    as labels next to the data points/arrows
    - Scat_Mag (bool): flag to switch on/off scatter plot (only points) of magnetic field
    - Con (boolt: flag to switch on/off plotting contours for z=const, i.e. flat surfaces with interpolation 
    between data points and the surfaces are stacked on top of each other.
    - Show (bool): flag to switch on/off plt.show() command 
    - title_on (bool): flag to switch title on or off
    If fig is not provided, a new figure with new axes are generated. 
    Make sure to also provide an axes object if a figure is passed.
    - single_comp (string): flag to choose how data should be plotted, possible values are None,'x','y','z','xy':
    None: euclidian norm/magnitude of magnetic field vector. 'x', 'y', 'z': plot the according vector component.
    'xy': euclidian norm of in-plane component of field (i.e. x,y components)
    - omit_64 (bool): flag to include (False) or exclude (True) the 64th sensor

    Return: fig, ax, (cf, cmin, cmax), u, v, w
    - fig, ax: Figure and Axis objects of the created plot. 
    - cf, cmin, cmax: colormap object and the (float) min and max values of the colorbar
    - u, v, w (ndarrays): vector components along x,y,z-axis in a suitable format for 3d-plotting
    """
    # initialize figure if not existing already
    if fig == None:
        # 3D Vector plot of magnetic flux
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_zlabel("z [mm]")
        ax.set_xlim(0, 15)
        ax.set_ylim(0, 15)
        ax.set_zlim(0, 15)

    # initialize positions of sensors
    x, y, z = np.meshgrid(np.arange(Pos[0], 20+Pos[0], 5),
                          np.arange(Pos[1], 20+Pos[1], 5),
                          np.arange(Pos[2], 20+Pos[2], 5))

    # reshuffle the vector components of magnetic field, such that they are in the correct format
    # u v, w: vector components along x,y,z-axis
    u = np.zeros(np.shape(x))
    v = np.zeros(np.shape(y))
    w = np.zeros(np.shape(z))
    for i in range(np.shape(x)[0]):
        for j in range(np.shape(x)[1]):
            for k in range(np.shape(x)[2]):
                # if omit_64 and int(48+j+4*i-16*k)==63: #account for read failure of sensor 64
                # account for read failure of sensor 64
                if omit_64 and int(60+j-4*i-16*k) == 63:
                    continue
                else:
                    # u[i,j,k] = SensorData[int(48+j+4*i-16*k),0]
                    # v[i,j,k] = SensorData[int(48+j+4*i-16*k),1]
                    # w[i,j,k] = SensorData[int(48+j+4*i-16*k),2]
                    u[i, j, k] = SensorData[int(60+j-4*i-16*k), 0]
                    v[i, j, k] = SensorData[int(60+j-4*i-16*k), 1]
                    w[i, j, k] = SensorData[int(60+j-4*i-16*k), 2]

    # choose how the data should be plotted
    if single_comp == None:
        mag = np.sqrt(u*u + v*v + w*w)
    elif single_comp == 'x':
        mag = u
    elif single_comp == 'y':
        mag = v
    elif single_comp == 'z':
        mag = w
    elif single_comp == 'xy':
        mag = np.sqrt(u*u + v*v)

    # set limits for colormap if they are not provided yet
    if cmin == None or cmax == None:
        cmin = np.amin(mag) - 1
        cmax = np.amax(mag) + 1

    # create the actual plot
    title = "Magnetic"
    if Vecs:  # add vectors into plot, magnetic field strength used for color code
        # get the correct color
        cmap = cm.get_cmap('viridis')
        normalized_mag = (mag - cmin) / (cmax - cmin)
        print((cmin, cmax))
        color_hex = np.zeros(np.shape(mag), dtype='U7')
        for i in range(np.shape(mag)[0]):
            for j in range(np.shape(mag)[1]):
                for k in range(np.shape(mag)[2]):
                    color_hex[i, j, k] = to_hex(
                        cmap(normalized_mag)[i, j, k, :])
        cf = ax.quiver(x, y, z, u, v, w, length=0.2,
                       colors=color_hex.flatten(), arrow_length_ratio=0.3)
        fig.colorbar(cf, ax=ax, label="|B| in mT")

        title += " field /"
        #ax.scatter(x,y,z, s=1)
        title += " field /"
    if Mag_label:  # add magnitude of magnetic field strength into plot
        for i in range(np.shape(x)[0]):
            for j in range(np.shape(x)[1]):
                for k in range(np.shape(x)[2]):
                    # if omit_64 and int(48+j+4*i-16*k)==63: #account for read failure of sensor 64
                    # account for read failure of sensor 64
                    if omit_64 and int(60+j-4*i-16*k) == 63:
                        continue
                    else:
                        M = "%.2f" % (
                            np.sqrt(u[i, j, k]**2 + v[i, j, k]**2 + w[i, j, k]**2))
                        label = '#{}: {} mT'.format(60+j-4*i-16*k + 1, M)
                        # label = '#{}: {} mT'.format(48+j+4*i-16*k + 1, M)
                        ax.text(x[i, j, k], y[i, j, k], z[i, j, k], label)
        title += " magnitude labels /"

    cf = None
    if Cont:  # add contours into plot
        cf = ax.contourf(x[0], z[0], np.transpose(
            mag[:, :, 0]), zdir='z', offset=0, vmin=cmin, vmax=cmax)  # , levels=20)
        ax.contourf(x[1], z[1], np.transpose(mag[:, :, 1]),
                    zdir='z', offset=5, vmin=cmin, vmax=cmax)  # , levels=20)
        ax.contourf(x[2], z[2], np.transpose(mag[:, :, 2]), zdir='z',
                    offset=10, vmin=cmin, vmax=cmax)  # , levels=20)
        ax.contourf(x[3], z[3], np.transpose(mag[:, :, 3]), zdir='z',
                    offset=15, vmin=cmin, vmax=cmax)  # , levels=20)

        fig.colorbar(cf, ax=ax, boundaries=(cmin, cmax), label="|B| in mT")
        title += " flux magnitudes"
    if Scat_Mag:  # scatter plot
        mag = mag.flatten()
        cf = ax.scatter(x, y, z, s=1, c=mag, vmin=cmin, vmax=cmax)
        fig.colorbar(cf, ax=ax, label="|B| in mT")
        title += " flux magnitudes"

    # if title ends with /, remove it
    if title[-1] == '/':
        title = title[:-1]

    if title_on:
        plt.title(title)
    plt.tight_layout()
    if Show:
        plt.show()
    return fig, ax, (cf, cmin, cmax), u, v, w


def plot_many_sets(directory, filename="means_grid_points.csv", Vecs=True, Cont=False, Scat_Mag=False, save=False, omit_64=False):
    """
    Generate 3d plot of many datasets, for example originating from running grid-measurements.

    Remark: if Scat_Mag = True or Cont = True, there are lots of colorbars, such that plot isn't 
    visible if there are too many data. Needs to be fixed later!

    Args:
    - directory (string): valid path of folder where the data file(s) that need(s) to be read is (are) located
    - filename (string): name of the data file that needs to be read
    - Vecs (bool): flag to switch on/off plotting of vectors
    - Scat_Mag (bool): flag to switch on/off scatter plot (only points) of magnetic field
    - Cont (bool): flag to switch on/off plotting contours for z=const, i.e. flat surfaces with interpolation 
    between data points and the surfaces are stacked on top of each other. 
    (This option is not reasonable if many z-layers are present)
    - save (bool): flag to save or not save the plot at the end
    - omit_64 (bool): flag to include (False) or exclude (True) the 64th sensor
    """
    # 3D Vector plot of magnetic flux
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    ax.set_zlim(0, 15)

    cwd = os.getcwd()
    access_rights = 0o755
    filedir = cwd + directory
    os.chmod(filedir, access_rights)    # won't have an effect on Windows

    # read in the measurement data from csv files
    mpoints = pd.read_csv(filedir + filename)
    if sys.version_info[0] == 3:
        mpoints = mpoints.to_numpy()  # use mpoints.values() if you use python 2
    else:
        mpoints = mpoints.values
    if omit_64:
        only_values = np.zeros((63, 3))  # ((np.shape(means))[0], 3))
    else:
        only_values = np.zeros((64, 3))

    cmin = None
    cmax = None
    cbar_par = None
    all_axes = []
    for i in range(np.shape(mpoints)[0]):
        meanname = "means_" + str(i+1) + ".csv"
        means = pd.read_csv(filedir + meanname)
        if sys.version_info[0] == 3:
            means = means.to_numpy()  # use means.values() if you use python 2
        else:
            means = means.values
        only_values[:, 0] = means[:, 0]
        only_values[:, 1] = means[:, 2]
        only_values[:, 2] = means[:, 4]
        fig, ax, cbar_par, _, _, _ = plot_set(only_values, fig=fig, ax=ax, Pos=mpoints[i, 1:], Vecs=Vecs,
                                              Cont=Cont, Scat_Mag=Scat_Mag, Show=False, title_on=False,
                                              cmin=cmin, cmax=cmax, omit_64=omit_64)
        cmin = cbar_par[1]
        cmax = cbar_par[2]
        all_axes.append(cbar_par[0])
    if Scat_Mag:
        fig.colorbar(all_axes[0], ax=ax, label="|B| in mT")

    plt.title("Magnetic field")
    if save:
        if Scat_Mag:
            plt.savefig(filedir + "plot_mags.png")
        elif Vecs:
            plt.savefig(filedir + "plot_vecs.png")
    plt.show()


def plot_stage_positions(directory, filename="means_grid_points.csv", save=False):
    """
    Plot the positions 

    Args:
    - directory (string): valid path of the folder that contains the csv file 
    - filename (string): name of the csv file 
    - save (bool): flag to save or not save the plot at the end
    """
    # extract stage positions from csv file
    cwd = os.getcwd()
    access_rights = 0o755
    filedir = cwd+directory
    os.chmod(cwd+directory, access_rights)  # won't have an effect on Windows
    mpoints = pd.read_csv(filedir+filename)
    if sys.version_info[0] == 3:
        mpoints = mpoints.to_numpy()  # use mpoints.values() if you use python 2
    else:
        mpoints = mpoints.values

    # create plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    ax.scatter(mpoints[:, 1], mpoints[:, 2], mpoints[:, 3])
    for i in range(1, np.shape(mpoints)[0]):
        vec = mpoints[i, 1:] - mpoints[i-1, 1:]
        ax.quiver(mpoints[i-1, 1], mpoints[i-1, 2],
                  mpoints[i-1, 3], vec[0], vec[1], vec[2])
    plt.title("Stage Positions: Movement of each sensor w.r.t its initial position")
    if save:
        plt.savefig(filedir+"plot_"+filename[:-4]+".png")
    plt.show()


def plot_sensor_positions(Pos=np.array([0, 0, 0])):
    """
    Generate 3d plot of the positions of the Hall sensors 

    Args:
    - Pos (1d-array of length 3): containx the offset position of cube (which is position of sensor 49)

    Return: fig, ax: Figure and Axis objects of the created plot. 
    """
    # initialize figure if not existing already
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 15)
    ax.set_zlim(0, 15)

    # initialize positions of sensors
    x, y, z = np.meshgrid(np.arange(Pos[0], 20+Pos[0], 5),
                          np.arange(Pos[1], 20+Pos[1], 5),
                          np.arange(Pos[2], 20+Pos[2], 5))

    z_offsets = np.arange(Pos[2], 20+Pos[2], 5)
    colors = np.array(['r', 'g', 'b', 'k'])

    # reshuffle the vector components of magnetic field, such that they are in the correct format
    # u v, w: vector components along x,y,z-axis
    u = np.zeros(np.shape(x), dtype=int)
    for i in range(np.shape(x)[0]):
        for j in range(np.shape(x)[1]):
            for k in range(np.shape(x)[2]):
                u[i, j, k] = int(48+j+4*i-16*k)

    # choose how the data should be plotted

    for i in range(np.shape(x)[0]):
        for j in range(np.shape(x)[1]):
            for k in range(np.shape(x)[2]):
                label = '#{}'.format(60+j-4*i-16*k + 1)
                ax.text(x[i, j, k], y[i, j, k], z[i, j, k], label,
                        color=colors[z_offsets == z[i, j, k]][0])

    plt.title('sensor positions', pad=20)
    plt.tight_layout()
    # plt.show()

    return fig, ax


def plot_angle(vec):
    """
    Generate 3d plot of the normalized input vector and show the angle with respect to z axis

    Arg: 
    - vec (ndarray of length 3): nonzero vector with x, y and z components
    """
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # normalize vector
    mag = norm(vec)
    vecn = vec/mag
    ang = np.arccos(vec[2]/mag)
    print(vecn)

    ax.quiver(0.0, 0.0, 0.0, vecn[0], vecn[1], vecn[2], arrow_length_ratio=0.3)
    ax.text(0.2, 0.2, 0.5, 'Angle: {:.1f} °'.format(np.degrees(ang)))

    # add red triangle between vector and z-axis to illustrate angle
    x = [0, 0, 0.8*vecn[0]]
    y = [0, 0, 0.8*vecn[1]]
    z = [0, 0.8, 0.8*vecn[2]]
    vtx = [list(zip(x, y, z))]
    tri = Poly3DCollection(vtx)
    tri.set_color('red')
    ax.add_collection3d(tri)

    # axis settings
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.5, 1.5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    plt.show()


def plot_angle_spherical(vec):
    """
    Generate 3d spherical plot of the normalized input vector and show the corresponding angles theta and phi.

    Arg: 
    - vec (ndarray of length 3): nonzero vector with x, y and z components
    """
    # normalize vector and estimate angles
    vecn = vec / norm(vec)
    theta = np.arccos(vec[2]/norm(vec))
    phi = np.arctan2(vec[1], vec[0])

    fig = plt.figure(figsize=plt.figaspect(1.))
    ax = fig.add_subplot(111, projection='3d')

    # plot arrows for x, y, z axis
    length_axes = 2.4
    ax.quiver(length_axes/2, 0, 0, length_axes, 0, 0, color='k',
              arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
    ax.quiver(0, length_axes/2, 0, 0, length_axes, 0, color='k',
              arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
    ax.quiver(0, 0, length_axes/2, 0, 0, length_axes, color='k',
              arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
    ax.text(1.45, 0, 0, 'x')
    ax.text(0, 1.35, 0, 'y')
    ax.text(0, 0, 1.3, 'z')

    # create a sphere
    u, v = np.mgrid[0:2*np.pi:16j, 0:np.pi:40j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_surface(x, y, z, color='k', rstride=1, cstride=1,
                    alpha=0.05, antialiased=False, vmax=2)  # cmap=cm.gray,

    # plot equator
    u, v = np.mgrid[0:2*np.pi:40j, np.pi/2:np.pi/2:1j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color='k', linewidth=0.5)

    # plot actual vector
    ax.quiver(0.0, 0.0, 0.0, vecn[0], vecn[1], vecn[2],
              arrow_length_ratio=0.2, color='b', linewidth=3)

    # add red triangle between vector and z-axis to illustrate theta-angle
    scaling_factor = 1.0
    x = [0, 0, scaling_factor*vecn[0]]
    y = [0, 0, scaling_factor*vecn[1]]
    z = [0, scaling_factor, scaling_factor*vecn[2]]
    vtx = [list(zip(x, y, z))]
    tri = Poly3DCollection(vtx)
    tri.set_color('red')
    tri.set_alpha(0.3)
    tri.set_linewidth(0.2)
    ax.add_collection3d(tri)

    # add red triangle between vector and x-axis in xy-plane to illustrate phi-angle
    x = [0, scaling_factor, scaling_factor*vecn[0]]
    y = [0, 0, scaling_factor*vecn[1]]
    z = [0, 0, 0]
    vtx = [list(zip(x, y, z))]
    tri = Poly3DCollection(vtx)
    tri.set_color('green')
    tri.set_alpha(0.3)
    tri.set_linewidth(0.2)
    ax.add_collection3d(tri)

    # print angles as title
    ax.set_title('$\\theta$ = {:.1f} °\n $\\phi$ =  {:.1f} °'.format(
        np.degrees(theta), np.degrees(phi)))

    # switch off axes and planes in background, rotate to nice position.
    ax.set_axis_off()
    ax.view_init(30, 45)
    plt.show()


if __name__ == "__main__":
    directory = '\\data_sets\\set_4_centered\\'
    plot_stage_positions(directory)
    plot_many_sets(directory, Vecs=True, Cont=False, Scat_Mag=False)

   
