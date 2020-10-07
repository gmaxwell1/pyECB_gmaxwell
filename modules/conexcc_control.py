# Author Jona Buehler 2020
# Documentation and Updates by Nicholas Meinhardt

#%%
import numpy as np
from time import sleep
import os

from conexcc.conexcc_class import *


def all_ready(CC1:ConexCC, CC2=None, CC3=None, timeout = 30):
    """
    Tests whether all controllers are ready, if they aren't, tries to achieve the READY state.

    Args:
    - CC1, CC2, CC3 are instances of conexcc_class
    - timeout [s] is maximum waiting time for trials. 

    Returns True if this is the case and False else. 
    
    Note: 
    - repeads checking state of controller every 0.2 seconds
    - if current state is DISABLED, calls to exit DISABLED. 
        Else, if current state is NOT REFERENCED, start homing.
    """
    print(CC1.axis, " axis:")
    ready1 = CC1.wait_for_ready(timeout = timeout)
    ready2 = True
    ready3 = True
    if CC2 != None:
        print(CC2.axis, " axis:")
        ready2 = CC2.wait_for_ready(timeout = timeout)
    if CC3 != None:
        print(CC3.axis, " axis:")
        ready3 = CC3.wait_for_ready(timeout = timeout)

    if ready1 and ready2 and ready3:
        print("All good: Ready for action!")
    elif not ready1:
        print("There was a problem with the ", CC1.axis, " axis!")
    elif not ready2:
        print("There was a problem with the ", CC2.axis, " axis!")
    elif not ready3:
        print("There was a problem with the ", CC3.axis, " axis!")

    return ready1 and ready2 and ready3

def get_coords(CC1, CC2=None, CC3=None):
    """
    Returns current position [mm] of actuators, sorted in the usual order [x,y,z].

    Args:
    - CC1, CC2, CC3 are instances of conexcc_class
    
    Returns:
    - pos: 1d-ndarray containing current positions [mm]
    - 1d-array containing the labels of axes in sorted order (same order as used for pos)
    """
    axes=[]
    vals=[]
    
    axes.append(ord(CC1.axis))
    vals.append(CC1.read_cur_pos())
    if CC2 != None:
        vals.append(CC2.read_cur_pos())
        axes.append(ord(CC2.axis))
    if CC3 != None:
        vals.append(CC3.read_cur_pos())
        axes.append(ord(CC3.axis))

    axes = np.asarray(axes)
    vals = np.asarray(vals)
    ind=np.unravel_index(np.argsort(axes), axes.shape)
    pos=vals[ind]
    ret_ax = []
    axes = axes[ind]
    for i in range(len(axes)):
        ret_ax.append(chr(axes[i]))
    return pos, np.asarray(ret_ax)

def reset_to(position, CC1, CC2=None, CC3=None):
    """
    Moves the actuators to a desired, absolute position. 

    Args:
    - position is a list or array of the desired new positions of actuators [mm]
    - CC1, CC2, CC3 are instances of conexcc_class

    Note:
    - The movement can only be successful if the controllers are in READY state! 
    Thus, keep in mind to ensure this beforehand
    """
    CC1.move_absolute(position[0])
    if CC2 != None:
        CC2.move_absolute(position[1])
    if CC3 != None:
        CC3.move_absolute(position[2])
    return 0

def is_at_goal(goal, CC1, CC2=None, CC3=None,  eps=1e-3):
    """
    Checks whether all current coordinates match the target coordinates (goal) up to a tolerance eps.

    Args: 
    - goal: target position as array 
    - CC1, CC2, CC3 are instances of conexcc_class
    - eps is acceptable tolerance

    Returns:
    - ok: True if all coordinates are within eps-range of goal, False else
    - dist: computed distance between current and goal coordinates (componentwise) 
    """
    real = get_coords(CC1, CC2=CC2, CC3=CC3)[0]
    try:
        diff= np.subtract(goal, real)
        ok = (abs(diff) < eps).all()
    except:
        diff = np.zeros(np.shape(goal))
        ok = False
        print("Error in computing difference: goal and real might not have the same shapes.")
    return ok, diff

def correct_reset(goal, CC1, CC2=None, CC3=None,  eps=1e-3, end_after=3):
    """
    Recursively compares the current position with the desired goal and updates position accordingly 
    if both positions deviate by more than eps in any of the coordinates. 

    Args:
    - goal: target position as array 
    - CC1, CC2, CC3 are instances of conexcc_class
    - eps is acceptable tolerance
    - end_after: maximum number the function calls itself recursively

    Returns 0 if the final coordinates are within eps-range of goal, else returns 1
    """
    ok, diff = is_at_goal(goal, CC1, CC2=CC2, CC3=CC3,  eps=eps)
    run=end_after
    if ok:
        print("Verification and correction of position", goal, " successful!")
        return 0
    elif run != 0:
        print("New correction necessary! Difference: ", diff)
        CC1.move_relative(diff[0])
        if CC2 != None:
            CC2.move_relative(diff[1])
        if CC3 != None:
            CC3.move_relative(diff[2])
        all_ready(CC1, CC2=CC2, CC3=CC3)
        correct_reset(goal, CC1, CC2=CC2, CC3=CC3, end_after=(run-1))
    else:
        print("Correction was unsuccessfull! Please adjust epsilon!")
        return 1

def check_no_motion(CC1, CC2=None, CC3=None, eps=1e-3, wait=0.5, end_after=3):
    """
    Recursively check whether the (absolute value of) current velocity of motors is 
    below an acceptable level eps. If this is the case, 0 is returned. 
    Else the system is set to sleep for a given waiting time [s] before checking again on velocity
    
    Args: 
    - CC1, CC2, CC3 are instances of conexcc_class
    - eps is acceptable tolerance for zero velocity
    - wait is the waiting time [s] between recursions 
    - end_after: maximum number the function calls itself recursively

    Returns True if all velocities are below acceptable limit and False else.
    """
    v1 = CC1.read_cur_vel()
    ok= bool(abs(v1)<eps)
    run=end_after
    if CC2 != None:
        v2 = CC2.read_cur_vel()
        ok= bool(ok and abs(v2)<eps)
    if CC3 != None:
        v3 = CC3.read_cur_vel()
        ok=bool(ok and abs(v3)<eps)
    if  ok:
        return True
    elif run != 0:
        sleep(wait)
        check_no_motion(CC1, CC2=CC2, CC3=CC2, eps=eps, wait=wait, end_after=(run-1))
    else:
        print("No-Motion conditions could not be met. Please adjust epsilon.")
        return False

def setup(reset_position, COM_ports = ['COM7', 'COM6', 'COM5']):
    """
    Set up connection to the three controllers of the x,y,z stages, activate them and 
    move actuators to reset_position. 
    
    For the last step, the actuators are moved to
    roughly the correct position first and fine-tuned afterwards. Eventually, it ensures 
    that all actuators are not in motion (up to an acceptably small, nonzero velocity).

    Args:
    - reset_position: array or list of length 3, containing the desired initial positions [mm]
    - COM_ports: array or list of length 3, containing the COM-PORT numbers of controllers,
    sorted as [x, y, z]

    Returns CC_X, CC_Y, CC_Z: instances of conexcc_class

    Note:
    - Check which controller is assigned to which COM_Port before running this function,
    else establishing of connection might fail or axes can be interchanged!
    - If controllers fail to reach READY state, the function still works without crashing, 
    p.e. when the power supply of controllers is disconnected. 
    """
    # establish connection
    print("==========================")
    print("Establishing connection...")
    CC_X = ConexCC(com_port = COM_ports[0], velocity=0.4, set_axis='x')
    print()
    CC_Y = ConexCC(com_port = COM_ports[1], velocity=0.4, set_axis='y')
    print()
    CC_Z = ConexCC(com_port = COM_ports[2], velocity=0.4, set_axis='z')
    print()

    # try to achieve READY state for all controllers
    all_ready(CC_X, CC2=CC_Y, CC3=CC_Z)

    # reset positions
    print("\nResetting...")
    reset_to(reset_position, CC_X, CC2=CC_Y, CC3=CC_Z)
    print("resetting to:", reset_position)

    # reactive controllers after MOVING state
    all_ready(CC_X, CC2=CC_Y, CC3=CC_Z)

    # verify that reset was successfull

    print("\n===========================")
    print("Fine Tuning of starting position...")
    print("Updated current Position:", get_coords(CC_X, CC2=CC_Y, CC3=CC_Z))
    correct_reset(reset_position, CC_X, CC2=CC_Y, CC3=CC_Z)
    all_ready(CC_X, CC2=CC_Y, CC3=CC_Z)

    #checking for sufficient motion stability
    check_no_motion(CC_X, CC2=CC_Y, CC3=CC_Z)

    print("\n============================")
    print("Setup successfull! Now starting: Measurement sequence...\n")

    return CC_X, CC_Y, CC_Z

def save_in_dir(means, directory, label, stds=None, coords=False):
    """
    Write the provided means and (optional) standard deviations to the file 'means_'+label+'.csv' 
    in provided directory.

    The coords flag distinguishes between magnetic field strength (False) and spatial coordinates (True) as source of the data.

    Args:
    - means: array or list of measured mean values of B-field or coordinates 
    - directory: path of the directory in which the file should be stored
    - label: a number or string to label the csv file  
    - stds: array or list of measured standard deviations of B-field or coordinates. 
    Should have at least the same size as means
    - coords (bool): Flag to switch between B-field (False) and spatial coordinates (True)
    """
    # Under Linux, user rights can be set with the scheme below, 
    # where 755 means read+write+execute for owner and read+execute for group and others. 
    # However, note that you can only set the file’s read-only flag with it under Windows.
    access_rights = 0o755
    os.chmod(directory, access_rights)

    output_file_name = 'means_{}.csv'.format(label)
    with open(os.path.join(directory, output_file_name), 'w') as f:
        if stds is not None and not coords:
            # no idea why we need carriage return (\r) is required here
            f.write('Mean Bx (mT), +- std dev x (mT), Mean By (mT), +- std dev y (mT), Mean Bz (mT), +- std dev z (mT)\r\n') #header line
            for i in range(np.shape(means)[0]):
                f.write('{},{},{},{},{},{}'.format(means[i,0], stds[i,0], means[i,1], stds[i,1], means[i,2], stds[i,2]) + '\r\n')
        elif not coords:
            f.write('Mean Bx (mT), Mean By (mT), Mean Bz (mT)\r\n') #header line
            for i in range(np.shape(means)[0]):
                f.write('{},{},{}'.format(means[i,0], means[i,1], means[i,2]) + '\r\n')
        elif coords and stds is None:
            f.write('Index, x [mm], y [mm], z [mm]\r\n') #header line
            for i in range(np.shape(means)[0]):
                f.write('{},{},{},{}'.format((i+1), means[i,0], means[i,1], means[i,2]) + '\r\n')
        elif coords and stds is not None:
            f.write('Index, x [mm], +- std x [mm], y [mm], +- std y [mm], z [mm], +- std z [mm]\r\n') #header line
            for i in range(np.shape(means)[0]):
                f.write('{},{},{},{},{},{},{}'.format((i+1), means[i,0], stds[i,0], means[i,1], stds[i,1], means[i,2], stds[i,2]) + '\r\n')
        else:
            f.write('Failed to save!\r\n')
            print("Failed to save to" , directory, "!")

def grid(CC_X:ConexCC, CC_Y:ConexCC, CC_Z:ConexCC, step_size=1, sweep_range=2, 
            start=np.array([1,1,1]), measurement_function=None, N=None, filename=None, cube=None):
    """
    Move the actuator on grid and wait for 1 s on each lattice point.
    
    Note:
    - The grid's origin is at the starting point, from where the grid extends towards the 
    positive direction along all three axes. In this way, a cube of length = srange is formed.
    All positions at the lattice points are saved in a separate file at the end.

    - If measurement_function is provided, measurements are performed at each lattice point.
    The measurement outcomes are saved in separate files, which are consecutively labeled, starting from 0. 

    - Currently, the start position should have some buffer along x and y compared to (0,0,0), 
    since relative motion might fail sometimes due to inprecision (moving relatively by -1 
    from 0.999 is not allowed, since -0.001 is out of the allowed range). -> Fix this in a later step

    Args:
    - CC1, CC2, CC3 are instances of conexcc_class
    - step_size (float): step size (lattice constant) of the grid [mm]
    - srange (float): maximum range of grid [mm], i.e. the length of grid along one dimension.
    - start: array or list of length 3, starting point of the grid. 
    - measurement_function is a function that takes following arguments: 
    measurement_function(N, filename=filename, cube=cube, no_enter=True, on_stage=True) and returns 
    mean_data, std_data, _, directory
    - N: number of times all 64 (or 63) sensor are read out in series,
    required if measurement_function is provided
    - filename: optional argument passed to measurement_function
    - cube: optional argument passed to measurement_function

    Returns: nd array of shape (number lattice points, 3), containing all lattice points of the grid
    """
    # make sure that starting position has already been reached, if not reset to this position
    ok = is_at_goal(start, CC_X, CC2=CC_Y, CC3=CC_Z)[0]
    if not ok:
        print("\nResetting...")
        reset_to(start, CC_X, CC2=CC_Y, CC3=CC_Z)

        print("resetting to:", start)
        all_ready(CC_X, CC2=CC_Y, CC3=CC_Z)

        # verify reset successful
        print("\n===========================")
        print("Fine Tuning of starting position...")
        print("Updated current Position:", get_coords(CC_X, CC2=CC_Y, CC3=CC_Z))
        correct_reset(start, CC_X, CC2=CC_Y, CC3=CC_Z)
        all_ready(CC_X, CC2=CC_Y, CC3=CC_Z)

        # checking for sufficient motion stability
        check_no_motion(CC_X, CC2=CC_Y, CC3=CC_Z)

    # initialize grid movement
    mpoints = []
    print("GRID RUN INITIALIZED \n\n")
    l=1 # counter
    norm_range = int(sweep_range/step_size)
    
    # first point of grid
    print("================POS ", l, "==================")
    check_no_motion(CC_X, CC2=CC_Y, CC3=CC_Z) # redundant if ok==False
    now = get_coords(CC_X, CC2=CC_Y, CC3=CC_Z)[0] # if everything works fine also redundant, but still good to consider cases where things went wrong
    print("Started at Position:", now)
    mpoints.append(now)
    # perform measurement
    if measurement_function is not None:
        print("Waiting for measurement...")
        mean_data, std_data, _, directory = measurement_function(N, filename=filename, cube=cube, 
                                                                    no_enter=True, on_stage=True)
        save_in_dir(mean_data, directory, l, stds=std_data)
    # wait for one second to finish measurement
    sleep(1)

    # subsequent points of grid
    l+=1
    for i in range(norm_range+1): # z-axis
        for j in range(norm_range+1): # x,y-plane
            for _ in range(norm_range):
                # move one step along X axis, alternatingly in opposite directions -> (-1)**j
                print("\n================POS ", l, "==================")
                #  ensure that the direction of motion in X-direction alternates:
                if norm_range % 2 ==0: # even number of steps 
                    direction_x = (-1)**(i+j)
                else:   # odd number of steps 
                    direction_x = (-1)**j
                CC_X.move_relative(direction_x * step_size)
                all_ready(CC_X)
                check_no_motion(CC_X, CC2=CC_Y, CC3=CC_Z)
                now=get_coords(CC_X, CC2=CC_Y, CC3=CC_Z)[0]
                mpoints.append(now)
                print("Moved to Position:", now)
                if measurement_function is not None:
                    print("Waiting for measurement...")
                    mean_data, std_data, _, directory = measurement_function(N, filename=filename, cube=cube, 
                                                                                no_enter=True, on_stage=True)
                    save_in_dir(mean_data, directory, l, stds=std_data)
                sleep(1)
                l+=1
            # after finishing all steps along x-axis, make one step along y-axis 
            # alternatingly in opposite directions -> (-1)**i
            if j != norm_range:
                print("\n================POS ", l, "==================")
                CC_Y.move_relative((-1)**i * step_size)
                all_ready(CC_Y)
                check_no_motion(CC_X, CC2=CC_Y, CC3=CC_Z)
                now=get_coords(CC_X, CC2=CC_Y, CC3=CC_Z)[0]
                mpoints.append(now)
                print("Moved to Position:", now)
                if measurement_function is not None:
                    print("Waiting for measurement...")
                    mean_data, std_data, _, directory = measurement_function(N, filename=filename, cube=cube, 
                                                                            no_enter=True, on_stage=True)
                    save_in_dir(mean_data, directory, l, stds=std_data)
                sleep(1)
                l+=1
            else:
                continue
        # after finishing sweeping in xy-plane, make one step along z direction
        if i != norm_range:
            print("\n================POS ", l, "==================")
            CC_Z.move_relative(step_size) 
            all_ready(CC_Z)
            check_no_motion(CC_X, CC2=CC_Y, CC3=CC_Z)
            now=get_coords(CC_X, CC2=CC_Y, CC3=CC_Z)[0]
            mpoints.append(now)
            print("Moved to Position:", now)
            if measurement_function is not None:
                print("Waiting for measurement...")
                mean_data, std_data, _, directory = measurement_function(N, filename=filename, cube=cube, no_enter=True, on_stage=True)
                save_in_dir(mean_data, directory, l, stds=std_data)
            sleep(1)
            l+=1
        else:
            continue

    # note that directory has not been defined if measurement_function = None
    if measurement_function is None:
        directory = os.getcwd()

    # save the locations of lattice points
    mpoints = np.asarray(mpoints)
    save_in_dir(mpoints, directory, 'grid_points', coords=True)
    [print(np.round(mpoints[i])) for i in range(len(mpoints))]

    print("\n\nFINISHED GRID RUN")
    print("Finished: Measurement sequence!")
    return mpoints

def close_connection(CC_X, CC_Y=None, CC_Z=None):
    """
    Closes the communication with controllers. 

    Args: CC_X, CC_Y, CC_Z: instances of conexcc_class

    Note: Closing the connection does not stop ongoing motion of the motor!
    """
    print("Now closing connection:")
    #end connection to actuators
    CC_X.close()
    if CC_Y is not None:
        CC_Y.close()
    if CC_Z is not None:
        CC_Z.close()
    print("Goodbye!")

#%%

if __name__ == '__main__':
    # initial parameters
    reset= np.array([8.264, 4.248, 1.0])#np.array([0,0,0])
    COM_ports = ['COM7', 'COM6', 'COM5']

    # set things up
    CC_X, CC_Y, CC_Z = setup(reset, COM_ports = COM_ports)

    #grid run
    #meas_points = grid(CC_X, CC_Y, CC_Z, steps=1, srange=2, start=reset)
    #close connection
    # close_connection(CC_X, CC_Y=CC_Y, CC_Z=CC_Z)



    # ------------------------Testing area---------------------------------------------------
    #%%
    # initial parameters
    reset= np.array([3, 3, 3])
    COM_ports = ['COM7', 'COM6', 'COM5']

    # set things up
    CC_X, CC_Y, CC_Z = setup(reset, COM_ports = COM_ports)
    reset_to(np.array([3, 3, 10]), CC_X, CC2=CC_Y, CC3=CC_Z)

    #%%
    lattice_points = grid(CC_X, CC_Y, CC_Z, step_size=3, sweep_range=9, start=[3,3,3])


    # %%
    close_connection(CC_X, CC_Y=CC_Y, CC_Z=CC_Z)

