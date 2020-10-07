# NewPort Conex-CC class
# Author Itay Shahak 2019
# Free code for the community!
# GitHub: https://gist.github.com/ishahak/bd025295dd8f6976dc962c6a02dec86b

# dependant on 'clr' which is PythonNet package
import clr
from time import sleep, time

# We assume Newport.CONEXCC.CommandInterface.dll is copied to our folder
if __name__ == '__main__':
    clr.AddReference("Newport.CONEXCC.CommandInterface")
else:
    clr.AddReference("./conexcc/Newport.CONEXCC.CommandInterface")

import CommandInterfaceConexCC

DEV = 1                # hardcoded here to the first device
MAX_VELOCITY = 0.4     # mm/s, by spec of NewPort TRA25CC DC Servo Motor


class ConexCC:
    """
    Class for Conex Controller and Actuator
    """

    def __init__(self, com_port, velocity, set_axis):
        self.min_limit = -1
        self.max_limit = -1
        self.cur_pos = -1
        self.controller_state = ''
        self.positioner_error = ''
        self.axis=set_axis
        self.driver = CommandInterfaceConexCC.ConexCC()
        ret = self.driver.OpenInstrument(com_port)
        if ret != 0:
            print('Oops: error opening port %s' % com_port)
            self.positioner_error = 'init failed'
        else:
            print('ConexCC: Successfully connected to %s' % com_port)
            self.read_velocity()
            self.set_velocity(velocity)
            self.set_homing_velocity(velocity)
            self.read_limits()
            self.read_cur_pos(direct_print=True)

    def wait_for_ready(self, timeout=60):
        """
        Repeatedly asks controller whether it is ready until it is either ready or the timeout limit has been reached.
        If controller is in:
        - NOT REFERENCED state, self.is_ready tries to activate it via homing first. 
        - DISABLED state, self.is_ready tries to enable controller again

        Returns True if controller is in READY state, and False else.
        """
        print('waiting for ready state...')#, end='')
        count = 0
        sleep_interval = 0.2
        last_count = (1 / sleep_interval) * timeout
        while not self.is_ready():
            count += 1
            if count % 30 == 0:
                print('<%s>' % self.controller_state)
            else:
                print('<%s>' % self.controller_state, end='', flush=True)
            sleep(sleep_interval)
            if count >= last_count:
                print('\nfailed to get ready. exiting for timeout = %d seconds.' % timeout)
                return False
        print('ok')
        return True

    def is_ready(self):
        """
        Checks current state of controller and tries achieve READY state. 

        Returns True if controller is in READY state and False else.
        """
        self.read_controller_state(silent=True)

        if self.controller_state in ('3D', '3C'):  # in DISABLE state
            self.exit_disable_state()
            sleep(0.2)
            self.read_controller_state()
        elif self.controller_state.startswith('0'):  # not referenced state
            self.init_positioner()
            sleep(0.4)

        # ('32','33','34') means in READY state
        ready = self.positioner_error == '' and self.controller_state in ('32', '33', '34')
        return ready

    @classmethod
    def dump_possible_states(cls):
        # https://www.newport.com/mam/celum/celum_assets/resources/CONEX-CC_-_Controller_Documentation.pdf#page=54

        help_text = 'see: dump_possible_states function in conexcc_class.py for further information'

        for s in help_text.split('\n'):
            print(s.strip(' '))

    def read_limits(self):
        err_str = ''
        resp = 0
        res, resp, err_str = self.driver.SL_Get(DEV, resp, err_str)
        if res != 0 or err_str != '':
            print('Oops: Negative SW Limit: result=%d,response=%.2f,errString=\'%s\'' % (res, resp, err_str))
        else:
            print('Negative SW Limit = %.1f' % resp)
            self.min_limit = resp

        res, resp, err_str = self.driver.SR_Get(DEV, resp, err_str)
        if res != 0 or err_str != '':
            print('Oops: Positive SW Limit: result=%d,response=%.2f,errString=\'%s\'' % (res, resp, err_str))
        else:
            print('Positive SW Limit = %.1f' % resp)
            self.max_limit = resp

    def read_cur_pos(self, direct_print=False):
        """
        Returns current position in mm and updates self.cur_pos accoringly. 
        """
        err_str = ''
        resp = 0
        res, resp, err_str = self.driver.TP(DEV, resp, err_str)
        if res != 0 or err_str != '':
            print('Oops: Current Position: result=%d,response=%.2f,errString=\'%s\'' % (res, resp, err_str))
        else:
            self.cur_pos = resp
            if direct_print:
                print('Current Position = %.3f' % resp)
        return resp


    def read_velocity(self):
        err_str = ''
        resp = 0
        res, resp, err_str = self.driver.VA_Get(DEV, resp, err_str)
        if res != 0 or err_str != '':
            print('Oops: Current Velocity: result=%d,response=%.2f,errString=\'%s\'' % (res, resp, err_str))
        else:
            print('Current Velocity = %.3f' % resp)

    def read_cur_vel(self, direct_print=False):
        """
        Estimates average velocity [mm/s] based on two conecutive queries of current position and 
        the elaped time between both queries. 
        """
        err_str1=''
        resp1=0
        time1=time()
        res1, resp1, err_str1 = self.driver.TP(DEV, resp1, err_str1)
        time2=time()
        if res1 != 0 or err_str1 != '':
            print('Oops: Current Real Velocity: result=%d,response=%.2f,errString=\'%s\'' % (res1, resp1, err_str1))
        else:
            err_str2=''
            resp2=0
            time3=time()
            res2, resp2, err_str2 = self.driver.TP(DEV, resp2, err_str2)
            time4=time()
            if res2 != 0 or err_str2 != '':
                print('Oops: Current Real Velocity: result=%d,response=%.2f,errString=\'%s\'' % (res2, resp2, err_str2))
            else:
                elapsed= ((time4-time1) + (time3-time2))/2
                av_velocity=(resp2-resp1)/elapsed
                if direct_print:
                    print('Average Velocity = %.3f during %.6f seconds' % (av_velocity, elapsed))
                return av_velocity


    def read_controller_state(self, silent=False):
        err_str = ''
        resp = ''
        resp2 = ''
        res, resp, resp2, errString = self.driver.TS(DEV, resp, resp2, err_str)
        if res != 0 or err_str != '':
            print('Oops: Read controller Err/State: result=%d,response=Err=\'%s\'/State=\'%s\',err_str=\'%s\'' % (
                res, resp, resp2, err_str))
        else:
            if not silent:
                print('Controller State = \'%s\', Error = \'%s\'' % (resp2, resp))
            self.positioner_error = resp
            self.controller_state = resp2

    def exit_disable_state(self):
        err_str = ''
        state = 1  # enable
        res, err_str = self.driver.MM_Set(DEV, state, err_str)
        if res != 0 or err_str != '':
            print('Oops: Leave Disable: result=%d,errString=\'%s\'' % (res, err_str))
        else:
            print('Exiting DISABLE state')

    def init_positioner(self):
        err_str = ''
        res, err_str = self.driver.OR(DEV, err_str)
        if res != 0 or err_str != '':
            print('Oops: Find Home: result=%d,errString=\'%s\'' % (res, err_str))
        else:
            print('Finding Home')

    def set_homing_velocity(self, velocity):
        if velocity > MAX_VELOCITY:
            velocity = MAX_VELOCITY
        err_str = ''
        res, err_str = self.driver.OH_Set(DEV, velocity, err_str)
        if res != 0 or err_str != '':
            print('Oops: Homing velocity: result=%d,errString=\'%s\'' % (res, err_str))
        else:
            print('Homing velocity set to %.1f mm/s' % velocity)

    def set_velocity(self, velocity):
        if velocity > MAX_VELOCITY:
            velocity = MAX_VELOCITY
        err_str = ''
        res, err_str = self.driver.VA_Set(DEV, velocity, err_str)
        if res != 0 or err_str != '':
            print('Oops: Set velocity: result=%d,errString=\'%s\'' % (res, err_str))
        else:
            print('velocity Set to %.1f mm/s' % velocity)

    def move_relative(self, distance):
        """
        Move the actuator by distance [mm]
        """
        self.read_cur_pos()
        if (self.cur_pos + distance > self.max_limit) or (self.cur_pos + distance < self.min_limit):
            print("Would move out of bounds: Movement interrupted!")
            return -1

        elif self.is_ready():
            err_str = ''
            res, err_str = self.driver.PR_Set(DEV, distance, err_str)
            if res != 0 or err_str != '':
                print('Oops: Move Relative: result=%d,errString=\'%s\'' % (res, err_str))
            else:
                print('Moving Relative %.3f mm' % distance)


    def move_absolute(self, new_pos):
        """
        Moves actuator to new_pos [mm]. 

        Note:
        - new_pos is the absolute position of desired position.
        - Error message will be printed if new_pos is invalid
        - Moving only works if controller is in READY state or 
            if READY can be reached in a single step from DISABLED or NOT REFERENCED
        """
        if (new_pos > self.max_limit) or (new_pos < self.min_limit):
            print("Would move out of bounds: Movement interrupted!")
            return -1

        elif self.is_ready():
            err_str = ''
            res, err_str = self.driver.PA_Set(DEV, new_pos, err_str)
            if res != 0 or err_str != '':
                print('Oops: Move Absolute: result=%d,errString=\'%s\'' % (res, err_str))
            else:
                print('Moving to position %.3f mm' % new_pos)

    def close(self):
        # note that closing the communication will NOT stop the motor!
        self.driver.CloseInstrument()


if __name__ == '__main__':
    ConexCC.dump_possible_states()
    conex_cc = ConexCC(com_port='com5', velocity=0.4,  set_axis='x')
    ready = conex_cc.wait_for_ready(timeout=60)
    if ready:
        conex_cc.move_absolute(conex_cc.max_limit / 5)
        ready = conex_cc.wait_for_ready(timeout=60)
        if ready:
            conex_cc.move_relative(-0.3)
            ready = conex_cc.wait_for_ready(timeout=60)
            if ready:
                print('ok!')
            else:
                print('not ok 2!')
        else:
            print('not ok 1!')
        conex_cc.close()
    else:
        print('something went wrong')