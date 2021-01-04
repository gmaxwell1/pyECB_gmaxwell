""" 
filename: generate_configurations.py

This script can be used to generate configuration files containing the ratios 
between the currents in all three coils, which can be passed to the main menu. 

The first part of this file contains some function definitions
In the second part, one can generate configurations with current ratios ranging 
from -1 to 1, where equidistant ratios are chosen. 
In the third part, one can generate configurations that should cover a number of roughly 
equidistant magnetic field values on the upper half sphere. Here, the desired
field values are chosen in the beginning, then they are transformed to current values
for the three coils and the current ratios are estimated at the end.


Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch
        
Date: 27.10.2020
"""

#%%
# standard library imports
import numpy as np
import os
import pandas as pd
from itertools import product, permutations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# local imports
try:
    import transformations as tr
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
    import transformations as tr
finally:
    from modules.general_functions import ensure_dir_exists
    from modules.analysis_tools import get_phi, get_theta
    from modules.interpolation_tools import delaunay_triangulation_spherical_surface, add_triangles_to_3dplot
<<<<<<< HEAD
=======

>>>>>>> a81a7624081bbf5abd530e7e1ca1c46f4d45ae23

#%%
# Part 1 --------------------------------------------------------
# function definitions
def generate_unique_combinations(num_vals, remove_negative=False):
    """
    Generate all possible cominations of current ratios (r1, r2, r3) between the three coils.
    The set of combinations fulfills fulfills following properties:
    - at least one value is 1, all other values are <= 1
    - each triplet of ratios is unique, i.e. the set contains each triplet only once
    - no two triplets are the equivalent up to their sign, i.e. multiplying a triplet by -1 
    does not produce a different element of the set
    - all possible combinations with num_vals values between -1 and 1 are contained in the set

    Arg:
    - num_vals: there are num_vals possible values of r1,r2,r3 that are equally 
    spaced between (incluively) -1 and 1. 
    num_vals must be odd to contain 0.

    Return: 
    -unique_combinations (ndarray): array containing all possible combinations of shape (number of combis, 3)

    Note:
    This may not be the most efficient implementation to generate the set, 
    better do not use it for large values of num_vals.
    """
    # generate all configurations in raw numbers
    generator_combinations = product(np.linspace(-1,1, num=num_vals), 
                                    np.linspace(-1,1, num=num_vals), 
                                    np.linspace(-1,1, num=num_vals))
    raw_ratios = np.array(list(generator_combinations))

    combinations_ratios = np.zeros((len(raw_ratios)*6,3))
    for i in range(len(raw_ratios)):
        combinations_ratios[i*6:(i+1)*6] = np.array(list(permutations(raw_ratios[i])))

    # remove duplicates
    unique_combinations = np.unique(combinations_ratios, axis=0)

    # remove combinations that are equal up to the sign
    indices_equivalent_pairs = []
    for i in range(len(unique_combinations)//2):
        for j in range(len(unique_combinations)):
            if np.all(unique_combinations[i] == -1*unique_combinations[j]):
                # print('{},{}: {}, {}'.format(i,j, unique_combinations[i], unique_combinations[j]))
                indices_equivalent_pairs.append(i)
    if remove_negative:
        unique_combinations = np.delete(unique_combinations, indices_equivalent_pairs, axis=0)

    # reverse order for convenience to have 111 at the beginning 
    unique_combinations = np.flip(unique_combinations, axis=0)

    return unique_combinations

def generate_configs_half_sphere(n_sectors, windings = 508, resistance = 0.47, 
                                elevation_factor_equator=None, magnitude=1,
                                upper = True, include_equator=True):
    """
    Generate current configurations that correspond to magnetic field vectors 
    in the upper or lower half sphere. The associated current values are estimated using the 
    linear model with an actuation matrix. Afterwards, the configurations are estimated,
    which are the ratios between the currents in the coils. 

    Note: Since the ratios are used to sweep along the defined field direction from 
    negative to positive, the sweeps along the equator may be problematic. 
    Thus, the elevation_factor_equator parameter allows to elevate the equator by a 
    given amount towards the north pole.

    Args:
    - n_sectors (int >= 4): number of points on the equator. The number of different
    latitude levels is n_sectors // 4. While there are n_sectors points on the equator 
    and a single point on the north pole, the number of points on a fixed latitude 
    (or elevation from the equator) is adapted automatically, such that all points are
    approximately uniformly distributed on the sphere. 
    Note: To ensure that points along x and y axis are considered, too, it is recommended
    to set n_sectors to multiples of 4.
    - windings (int), resistance (float): properties of coils, used for the 
    transformation of field vectors to current values. 
    - elevation_factor_equator (None or float): If a value is provided, the vectors along 
    the equator are elevated by elevation_equator * (next higher latitude). Hence, 
    elevation_factor_equator should be in [0,1].
    - magnitude (float): magnitude of vectors that are generated
    - upper (bool): If True, the upper hemisphere is considered, else the lower 
    - include_equator (bool): If False, all vectors that are on the equator will be skipped, 
    which includes the points that are lifted from the equator by passing elevation_factor_equator.
    If True (default), the vectors within the equator plane are considered. 

    Return: 
    - ratios (ndarray of shape(N, 3)): Contain the ratios of the three current values.
    For normalization, at least one of the ratios is set to an absolute value of 1. 
    The sign remains the same, such that a single set of ratios contains a 1 or a -1.
    - vectors (ndarray of shape(N,3)): Array containing the Carthesian coordinates
    of all considered vectors, posing a simple way for plotting of the latter
    - thetas (ndarray of length 3): Contain latitude angles theta of all vectors in radians
    - phis (ndarray of length 3): Contain longitude angles phi of all vectors in radians
    """
    # estimate the distance between two points on equator along the circle
    distance_equator = 2 * np.pi / n_sectors

    # devide latitudes in n_sectors//4 sections, but ensure that latitudes contain at least 
    # two values if the equator should be included. 
    # Note that different values are required for upper and lower hemisphere
    if include_equator:
        n_lat = n_sectors//4 + 1
    else:
        n_lat = n_sectors//4 

    if upper:
        latitudes = np.linspace(0, np.pi/2, n_lat, endpoint=include_equator)
    else:
        latitudes = np.linspace(np.pi, np.pi/2, n_lat, endpoint=include_equator)

    # if desired, lift latitude of vectors on equator
    if elevation_factor_equator is not None and include_equator:
        latitudes[-1] -= elevation_factor_equator * (latitudes[-1] - latitudes[-2])

    # prepare lists to collect field vectors and current ratios
    ratios = []
    vectors = []
    thetas = []
    phis = []
    for theta in latitudes:
        # estimate radius of circle at given latitude, when the equator has radius 1
        radius = np.sqrt(1 - np.cos(theta)**2)
        if (2*np.pi*radius) < distance_equator:
            n_longitudinal = 1
        else:
            n_longitudinal = int(np.ceil( 2*np.pi*radius / distance_equator))
        longitudes = np.linspace(0, 2 * np.pi, n_longitudinal, endpoint=False)

        for phi in longitudes:
            # estimate vector in Carthesian coordinates from angles
            B_vector = magnitude * np.array([np.sin(theta) * np.cos(phi),
                                            np.sin(theta) * np.sin(phi),
                                            np.cos(theta)])

            # estimate the corresponding currents using linear model
            I_coils = tr.computeCoilCurrents(B_vector, windings, resistance)

            # collect the ratios of the three currents, where at least one value has absolute value 1
            # and the current directions remain the same 
            i_max = np.argmax(np.abs(I_coils))
            ratios.append(I_coils / I_coils[i_max] * np.sign(I_coils[i_max]))
            vectors.append(B_vector)

            # collect angular configuration
            thetas.append(theta)
            phis.append(phi)

    return np.array(ratios), np.array(vectors), np.array(thetas), np.array(phis)

def generate_test_points_whole_sphere(n_sectors, magnitude):
    """
    Generate a test set of points that are approximately equally spaced on a sphere with radius
    magnitude and return a set of vectors.

    Args:
    - n_sectors (int >= 4): number of points on the equator. The number of different
    latitude levels is n_sectors // 4. While there are n_sectors points on the equator 
    and a single point on the north pole, the number of points on a fixed latitude 
    (or elevation from the equator) is adapted automatically, such that all points are
    approximately uniformly distributed on the sphere. 
    Note: To ensure that points along x and y axis are considered, too, it is recommended
    to set n_sectors to multiples of 4.
    - magnitude (float): magnitude of vectors that are generated
    """
    ratios1, vectors_upper, _, _ = generate_configs_half_sphere(n_sectors, magnitude=magnitude, 
                                        upper=True, include_equator=True)
    ratios2, vectors_lower, _, _ = generate_configs_half_sphere(n_sectors, magnitude=magnitude, 
                                        upper=False, include_equator=False)

    # combine both hemispheres and return vectors
    return np.append(ratios1, ratios2, axis=0), np.append(vectors_upper, vectors_lower, axis=0)



def rng_test_points_whole_sphere(N=10, magnitude_range=[0,50], theta_range=[0,np.pi], phi_range=[0,2*np.pi], seed=None):
    """
    Generate magnetic field vectors in in random directions and with random magnitudes. The angles
    theta and phi follow a uniform distribution between 0 and pi or 0 and 2pi respectively (ranges can
    be customized). The magnitudes follow a power distribution on a customizable interval.
    The seed can be set to generate new but reproducible combinations of vectors.

    Args:
    - N (int >= 1): number of directions generated
    - magnitude_range (list): magnitude of vectors that are generated are between the two specified numbers.
    they are randomly distributed according to the 2nd order power distribution. P(x) = 3*x^2, for 0<x<=1
    default: [0,50]
    - theta_range (list): polar angle of vectors that are generated are between the two specified numbers.
    they are randomly distributed according to the uniform distribution.
    default: [0,np.pi]
    - phi_range (list): azimuthal angle of vectors that are generated are between the two specified numbers.
    they are randomly distributed according to the uniform distribution.
    default: [0,2*np.pi]
    - seed (int): A seed to initialize the `BitGenerator`. If None, then fresh, unpredictable entropy will 
    be pulled from the OS. If an `int` or `array_like[ints]` is passed, then it will be passed to `SeedSequence`
    to derive the initial `BitGenerator` state.
    
    Return: 
    - vectors (ndarray of shape(N,3)): Array containing the Cartesian coordinates
    of all considered vectors, posing a simple way for plotting of the latter
    - B_magnitudes (ndarray of length N): Array containing the magnitudes of all considered
    vectors, posing a simple way for plotting the latter
    - thetas (ndarray of length N): Contains latitude angles theta of all vectors in radians
    - phis (ndarray of length N): Contains longitude angles phi of all vectors in radians
    """
    #initialize random generator
    rng = np.random.default_rng(seed)
    
    # prepare list to collect field vectors
    vectors = []
    # randomly generated vectors
    B_magnitudes = magnitude_range[0] + rng.power(2, N) * (magnitude_range[1] - magnitude_range[0])
    thetas = theta_range[0] + rng.random((N,)) * (theta_range[1] - theta_range[0])
    phis = phi_range[0] + rng.random((N,)) * (phi_range[1] - phi_range[0])
    
    for i in range(N):
        # estimate vector in Cartesian coordinates from angles
        B_vector = B_magnitudes[i] * np.array([np.sin(thetas[i]) * np.cos(phis[i]),
                                        np.sin(thetas[i]) * np.sin(phis[i]),
                                        np.cos(thetas[i])])
        
        vectors.append(B_vector)

    return np.array(vectors), np.array(B_magnitudes), np.array(thetas), np.array(phis)


def plot_vectors(vectors, magnitude = 1, phis= None, thetas=None, add_tiangulation = False):
    """
    Generate and show 3d plot of sphere and the provided vectors.

    Args:
    - vectors (ndarray of shape(N, 3)): Contains normal vectors that should be plotted.
    - magnitude (float): magnitude of sphere that is plotted for a better visualization.
    If magnitude = 1, it corresponds to the unit sphere
    - phis, thetas (ndarrays of length N): Spherical angles in degrees of passed vectors. Kind of redundant,
    but it is easier to reuse the already defined angles than computing them from the vectors again.
    - add_tiangulation (bool): If True, triangulation between points is plotted as well.
    """
    # plot the generated vectors
    # generate figure with 3d-axis
    fig = plt.figure(figsize = 1.5*plt.figaspect(1.))
    ax = fig.add_subplot(111, projection='3d')

    # plot arrows for x, y, z axis
    length_axes = 2.8 * magnitude 
    ax.quiver(length_axes/2, 0, 0, length_axes, 0, 0, color='k',
            arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
    ax.quiver(0, length_axes/2, 0, 0, length_axes, 0, color='k',
            arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
    ax.quiver(0, 0, length_axes/2, 0, 0, length_axes, color='k',
            arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
    ax.text(1.6*magnitude, 0, 0, 'x')
    ax.text(0, 1.65*magnitude, 0, 'y')
    ax.text(0, 0, 1.6*magnitude, 'z')

    # create a sphere
    u, v = np.mgrid[0:2*np.pi:16j, 0:np.pi:40j]
    x = magnitude * np.cos(u)*np.sin(v)
    y = magnitude * np.sin(u)*np.sin(v)
    z = magnitude * np.cos(v)
    ax.plot_surface(x, y, z, color='k', rstride=1, cstride=1,
                        alpha=0.05, antialiased=False, vmax=2)  
    
    # plot all vectors as red dots
    ax.scatter(vectors[:,0], vectors[:,1], vectors[:,2], color='r')
    
    if add_tiangulation:
        if phis is None or thetas is None:
            phis = get_phi(vectors)
            thetas = get_theta(vectors)
        # apply Delaunay triangulation on the surface of the sphere, which corresponds to finding convex Hull
        _, inidces_simplices, points = delaunay_triangulation_spherical_surface(phis, 
                                                                                thetas, 
                                                                                radius=magnitude)

        # add triangles as lines
        add_triangles_to_3dplot(ax, points, inidces_simplices, spherical = True, colored_triangles = False,
                                    color='r')

    ax.set_xlabel('$B_x$ [mT]')
    ax.set_ylabel('$B_y$ [mT]')
    ax.set_zlabel('$B_z$ [mT]')
    # ax.set_axis_off()

    ax.view_init(30, 45)

    if add_tiangulation:
        return fig, ax, inidces_simplices, points
    else:
        return fig, ax


def plot_vectors_simple(vectors, magnitudes = 1):
    """
    Generate and show 3d plot of sphere and the provided vectors.

    Args:
    - vectors (ndarray of shape(N, 3)): Contains normal vectors that should be plotted.
    - magnitude (float/ndarray of shape(N,)): magnitude of each vector that is plotted.
    """
    # plot the generated vectors
    # generate figure with 3d-axis
    fig = plt.figure(figsize = 1.5*plt.figaspect(1.))
    ax = fig.add_subplot(111, projection='3d')

    # plot arrows for x, y, z axis
    length_axes = 50
    ax.quiver(length_axes/2, 0, 0, length_axes, 0, 0, color='k',
            arrow_length_ratio=0.1, pivot='tip', linewidth=1.)
    ax.quiver(0, length_axes/2, 0, 0, length_axes, 0, color='k',
            arrow_length_ratio=0.1, pivot='tip', linewidth=1.5)
    ax.quiver(0, 0, length_axes/2, 0, 0, length_axes, color='k',
            arrow_length_ratio=0.1, pivot='tip', linewidth=1.5)
    ax.text(40, 0, 0, 'x')
    ax.text(0, 40, 0, 'y')
    ax.text(0, 0, 40, 'z')
    
    # plot all vectors as arrows
    for index in range(len(vectors)):
        ax.quiver(0, 0, 0, vectors[index, 0], vectors[index, 1], vectors[index, 2], length=magnitudes[index], color='C2',
            arrow_length_ratio=0.1, pivot='tail', linewidth=0.9)
    
    
 
    ax.set_xlabel('$B_x$ [mT]')
    ax.set_ylabel('$B_y$ [mT]')
    ax.set_zlabel('$B_z$ [mT]')
    # ax.set_axis_off()
    
    ax.view_init(30, 45)

    return fig, ax



#%%
# Part 2 -----------------------------------------------------------
# generate configurations based on equidistant current ratios 

# if __name__ == '__main__':
#     # set the number of values between (incl) -1 and 1 that should be considered 
#     num_vals = 3

#     # for num_vals in range(2,10):
#         # generate the set
#     unique_combinations = generate_unique_combinations(num_vals, remove_negative=False)

#         # save the combinations to csv file
#     directory = r'.\config_files'

#     df = pd.DataFrame({ 'ratio coil 1': unique_combinations[:,0], 
#                         'ratio coil 2': unique_combinations[:,1], 
#                         'ratio coil 3': unique_combinations[:,2]})

#     output_file_name = 'configs_numvals{}_length{}.csv'.format(num_vals, len(unique_combinations))
#     data_filepath = os.path.join(directory, output_file_name)
#     df.to_csv(data_filepath, index=False, header=True)

# # %%
# Part 3 -----------------------------------------------------------
# generate configurations based on (approximately) equidistant 
# magnetic fields in upper half plane

if __name__ == '__main__':
    # generate configurations
    n_vectors = 5000
    magnitude_range = [0,70]
    seed = 1414
    vectors,magnitudes,thetas,phis = rng_test_points_whole_sphere(n_vectors, magnitude_range=magnitude_range, seed=seed)

    # plot all considered vectors on a sphere 
    plot_vectors(vectors,50)
    
    plt.show()

    thetas_deg = thetas * 180/np.pi
    phis_deg = phis * 180/np.pi
    # save the combinations to csv file
    directory = '../config_files/RNG_test_vectors'
        
    df = pd.DataFrame({ 'B_x': vectors[:,0], 
                        'B_y': vectors[:,1], 
                        'B_z': vectors[:,2],
                        'B_mag': magnitudes,
                        'theta (deg)': thetas_deg,
                        'phi (deg)': phis_deg})

    output_file_name = f'vectors_rng{seed}_{magnitude_range[0]}-{magnitude_range[1]}mT_size{len(vectors)}.csv'
    data_filepath = os.path.join(directory, output_file_name)
    df.to_csv(data_filepath, index=False, header=True)



#%%
def generate_grid(max_value, points_per_dim, threshold_magnitude = np.inf):
    """
    Generate 3d rectangular grid with points_per_dim vertices per dimension, which 
    are equally spaced in [-max_value, +max_value]. The order of the returned array
    is such that neighboring vertices in the array are also neighbors in real space.

    Notes:
    - The origin is included if points_per_dim is odd and excluded if 
    points_per_dim is even.
    - max_value specifies the maximum value for each dimension. This means that
    the points at the corners have a magnitude of sqrt(3)*max_value! 
    - Depending on the valueof threshold_magnitude, the returned grid points
    may represent a cube or the intersection of a cube with a ball. 

    Args:
    - max_value (float): maximum value along each dimension
    - points_per_dim (int): Number of vertices along each dimension
    - threshold_magnitude (float): Maximum magnitude that is allowed. Default is np.inf, 
    meaning that the entire rectangular grid is returned. If a value is provided, 
    only vectors with magnitudes below the threshold are added to the final array.

    Returns: 
    - grid_pts (ndarray of shape (points_per_dim**3, 3))
    """
    x = np.linspace(-max_value, max_value, points_per_dim)

    xx, yy, zz = np.meshgrid(x,x,x)

    below_thresh = np.sqrt(xx**2 + yy**2 + zz**2) <= threshold_magnitude

    # initialize grid movement
    grid_pts = []

    for k in range(points_per_dim): # z-axis
        if k % 2 == 0:
            y_range = np.arange(points_per_dim)
        else: 
            y_range = np.flip(np.arange(points_per_dim))

        for j in y_range: # y-axis

            if (-1)**(k+j) == 1:
                x_range = np.arange(points_per_dim)
            else: 
                x_range = np.flip(np.arange(points_per_dim))

            for i in x_range: # x-axis

                # only add a grid point if its magnitude is below the threshold
                if below_thresh[i,j,k]:
                    grid_pts.append([xx[i,j,k], yy[i,j,k], zz[i,j,k]])

    grid_pts = np.array(grid_pts)
    
    return grid_pts

# %%
# if __name__ == '__main__':
    
#     # generate grid points
#     # max_value = 2
#     # points_per_dim = 51
#     # grid_pts = generate_grid(max_value, points_per_dim, threshold_magnitude=np.inf)

#     # directory = './config_files/grid/'
#     # output_file_name = f'walk_on_grid_max{max_value}_PointsPerDim{points_per_dim}.csv'     
#     # ensure_dir_exists(directory)
#     # data_filepath = os.path.join(directory, output_file_name)
    
#     # df = pd.DataFrame({ 'x': grid_pts[:, 0], 
#     #                     'y': grid_pts[:, 1], 
#     #                     'z': grid_pts[:, 2]})
        
#     # df.to_csv(data_filepath, index=False, header=True)
    
#     ratios,vectors = generate_test_points_whole_sphere(24, 50)
    
#     plot_vectors(vectors)
#     plt.show()
#     # save the combinations to csv files
#     directory = r'.\config_files'
    


#     df = pd.DataFrame({ 'ratio coil 1': ratios[:,0], 
#                         'ratio coil 2': ratios[:,1], 
#                         'ratio coil 3': ratios[:,2]})

#     output_file_name = 'configs_wholeSphere_length{}.csv'.format(len(ratios))
#     data_filepath = os.path.join(directory, output_file_name)
#     df.to_csv(data_filepath, index=False, header=True)

