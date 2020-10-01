// corrections added by gmaxwell on 28-09-2020
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <ecb_api.h>

static PyObject* py_initECBApi(PyObject* self, PyObject* args) {
    char * server_name;
    char * port;
    ecb_api_error_t result;
    PyObject * ret;

    // parse arguments
    if (!PyArg_ParseTuple(args, "ss", &server_name, &port)) {
        return NULL;
    }

    // run the actual function
    result = initECBApi(server_name, port);

    // build the resulting string into a Python object.
    ret = Py_BuildValue("l",result);

    return ret;
}

static PyObject* py_exitECBApi(PyObject* self){
    ecb_api_error_t result;
    PyObject * ret;

    // run the actual function
    result = exitECBApi();

    // build the resulting string into a Python object.
    ret=Py_BuildValue("l", result);

    return ret;
}

static PyObject* py_enableECBCurrents(PyObject* self) {
    ecb_api_error_t result;
    PyObject * ret;

    // run the actual function
    result = enableECBCurrents();

    // build the resulting string into a Python object.
    ret=Py_BuildValue("l", result);

    return ret;
}

static PyObject* py_disableECBCurrents(PyObject* self) {
    ecb_api_error_t result;
    PyObject * ret;

    // run the actual function
    result = disableECBCurrents();

    // build the resulting string into a Python object.
    ret=Py_BuildValue("l", result);

    return ret;
}

static PyObject* py_setDesCurrents(PyObject* self, PyObject* args) {

    /* the O! parses for a Python object (listObj) checked
     *        to be of type PyList_Type */
    PyObject* list_obj;
    unsigned char direct;

    if (! PyArg_ParseTuple( args, "O!c", &PyList_Type, &list_obj, &direct ))
        return NULL;

    // for now only support 8 coil configuration
    int curr_arr[8*40];

    // check if the arg is a block of current vectors rather than a single vector
    if (PyList_Size(list_obj) == 0) {
        printf("currents must a nonempty list\n");
        return NULL;
    }

    PyObject* first_list_obj = PyList_GetItem(list_obj, 0);

	//if the first object is a list object
    if ( PyList_Check(first_list_obj)) {
        size_t num_blocks = PyList_Size(list_obj);
        if( num_blocks != 40) {
            printf("number of blocks must be 40\n");
            return NULL;
        }
        int b = 0;
        for(; b < num_blocks; b++) {
            PyObject* block_obj;
            block_obj = PyList_GetItem(list_obj, b);

            size_t num_currents = PyList_Size(block_obj);
            if (num_currents != 8)
            {
                printf("number of currents must be 8\n");
                return NULL;
            }

            int c = 0;
            for(; c < num_currents; c++) {
                PyObject* cur_obj = PyList_GetItem(block_obj, c);
                curr_arr[c + b*8] = PyLong_AsLong(cur_obj);
            }
        }

    } else {
        // single vector, we then pad block
        size_t num_currents = PyList_Size(list_obj);

        // for now only support 8 coil configuration
        if (num_currents != 8) {
            printf("number of currents must be 8\n");
            return NULL;
        }

        int i=0;
        for (; i < num_currents; i++) {
            PyObject* cur_obj;
            cur_obj = PyList_GetItem(list_obj, i);

            // if (!PyLong_Check(cur_obj)) {
            //     printf("not long\n");
            //     return NULL;
            // }

            int curr = PyLong_AsLong(cur_obj);
            int j=0;
            for(; j < 40; j++) {
                curr_arr[i + j*8] = curr;
            }
        }
    }

    ecb_api_error_t result = setDesCurrents(curr_arr, direct);

    return Py_BuildValue("l",result);
}

static PyObject* py_getActCurrents( PyObject* self ) {
    PyObject* curr_list = PyList_New(8);
    int act_curr[8];

    ecb_api_error_t result = getActCurrents(act_curr);

    PyObject* ret_code  = Py_BuildValue("l", result);

    if (result == ERROR_NOERROR) {
        int c=0;
        for(; c<8; c++) {
            // correction added by gmaxwell on 28-09-2020
            PyList_SetItem(curr_list, c, Py_BuildValue("l", act_curr[c]));
        }
    }

    return PyTuple_Pack(2, ret_code, curr_list);

} 

static PyObject* py_setHeartbeatTimeout(PyObject* self, PyObject* args) {

    unsigned int timeout;
    if (! PyArg_ParseTuple(args, "i", &timeout)) 
        return NULL;

    ecb_api_error_t ret_val = setHeartbeatTimeout(&timeout);

    return Py_BuildValue("l", (ret_val));
}

static PyObject* py_getHeartbeatTimeout(PyObject* self) {

    unsigned int timeout;
    ecb_api_error_t ret_val = getHeartbeatTimeout(&timeout);

    PyObject* py_ret = Py_BuildValue("l", ret_val);
    PyObject* py_timeout = PyLong_FromUnsignedLong(timeout);
    return PyTuple_Pack(2, py_ret, py_timeout);
}

static PyObject* py_getECBRevision(PyObject* self) {

    unsigned int rev;
    ecb_api_error_t ret_val = getECBRevision(&rev);

    return PyTuple_Pack(2, PyLong_FromUnsignedLong(ret_val), PyLong_FromUnsignedLong(rev));
} 

static PyObject* py_setMaxCurrent(PyObject* self, PyObject* args) {

    unsigned int max_curr;
    if (! PyArg_ParseTuple(args, "i", &max_curr))
        return NULL;

    ecb_api_error_t ret = setMaxCurrent(&max_curr);

    return Py_BuildValue("l", ret);
}

static PyObject* py_getMaxCurrent(PyObject* self) {

    unsigned int max_curr[8];

    ecb_api_error_t ret = getMaxCurrent(max_curr);

    PyObject* py_ret = Py_BuildValue("l", ret);
    PyObject* py_max_curr = PyList_New(8);

    int c;
    for (c=0; c < 8; c++) {
        PyList_SetItem(py_max_curr, c, Py_BuildValue("l", max_curr[c]));
    }

    return PyTuple_Pack(2, py_ret, py_max_curr);
}

static PyObject* py_getCoilStatus(PyObject* self) {

    unsigned int coil_status[8];

    ecb_api_error_t ret = getCoilStatus(coil_status);

    PyObject* py_ret = Py_BuildValue("l", ret);
    PyObject* py_coil_status_list = PyList_New(8);

    if (ret == ERROR_NOERROR) {
        int c;
        for (c = 0; c < 8; c++) {
            PyList_SetItem(py_coil_status_list, c, PyLong_FromUnsignedLong(coil_status[c]));
        }
    }

    return PyTuple_Pack(2, py_ret, py_coil_status_list);
}

static PyObject* py_getCoilValues(PyObject* self) {

    unsigned int coil_temp[8];
    int hall[8];
    int currents[8];
    unsigned int coil_status[8];

    ecb_api_error_t ret = getCoilValues(coil_temp, hall, currents, coil_status);

    PyObject* py_ret = Py_BuildValue("l", ret);
    PyObject* py_temp_list = PyList_New(8);
    PyObject* py_hall_list = PyList_New(8);
    PyObject* py_currents_list = PyList_New(8);
    PyObject* py_coil_status = PyList_New(8);

    if (ret == ERROR_NOERROR) {
        int c;
        for (c = 0; c < 8; c++) {
            PyList_SetItem(py_temp_list, c, PyLong_FromUnsignedLong(coil_temp[c]));
            PyList_SetItem(py_hall_list, c, Py_BuildValue("l", hall[c]));
            PyList_SetItem(py_currents_list, c, Py_BuildValue("l", currents[c]));
            PyList_SetItem(py_coil_status, c, PyLong_FromUnsignedLong(coil_status[c]));
        }
    }

    return PyTuple_Pack(5, py_ret, py_temp_list, py_hall_list, py_currents_list, py_coil_status);

}

static PyObject* py_getECBStatus(PyObject* self) {
    unsigned int ecb_status;
    ecb_api_error_t ret = getECBStatus(&ecb_status);

    return PyTuple_Pack(2, Py_BuildValue("l", ret), PyLong_FromUnsignedLong(ecb_status));
}

static PyObject* py_checkECBComm(PyObject* self) {

    ecb_api_error_t ret = checkECBComm();

    return Py_BuildValue("l", ret);
}

static PyObject* py_errorToString(PyObject* self, PyObject* args) {

    ecb_api_error_t e;

    if (! PyArg_ParseTuple(args, "i", &e)) 
        return NULL;

    const char* err_str = errorToString(e);

    return Py_BuildValue("s", err_str);
}

static PyMethodDef ECBMethods[] = {
    { "initECBapi", py_initECBApi, METH_VARARGS, "Initializes the ECB API and the communication to the ECB" },
    { "exitECBapi", py_exitECBApi, METH_NOARGS, "Exits the ECB API and closes the ECB communication"},
    { "enableECBCurrents", py_enableECBCurrents, METH_NOARGS, "Enables the ECB currents on all coils"},
    { "disableECBCurrents", py_disableECBCurrents, METH_NOARGS, "Disables the ECB currents on all coils"},
    { "setDesCurrents", py_setDesCurrents, METH_VARARGS, "Sets the desired currents of the ECB"},
    { "getActCurrents", py_getActCurrents, METH_NOARGS, "Fetch measured current values" },
    { "setHeartbeatTimeout", py_setHeartbeatTimeout, METH_VARARGS, ""},
    { "getHeartbeatTimeout", py_getHeartbeatTimeout, METH_NOARGS, ""},
    { "getECBRevision", py_getECBRevision, METH_NOARGS, "" },
    { "setMaxCurrent", py_setMaxCurrent, METH_VARARGS, "" },
    { "getMaxCurrent", py_getMaxCurrent, METH_NOARGS, "" },
    { "getCoilStatus", py_getCoilStatus, METH_NOARGS, ""},
    { "getCoilValues", py_getCoilValues, METH_NOARGS, ""},
    { "getECBStatus", py_getECBStatus, METH_NOARGS, ""},
    { "checkECBComm", py_checkECBComm, METH_NOARGS, ""},
    { "errorToString", py_errorToString, METH_VARARGS, ""},
    { NULL, NULL, 0, NULL }
};

static struct PyModuleDef ecb_definition = {
    PyModuleDef_HEAD_INIT,
    "ecb",
    "A python wrapper module for controlling the Pantec ECB.",
    -1,
    ECBMethods
};

/*
void initECB(void)
{
    PyObject* m = Py_InitModule("ECB", ECBMethods);
    PyModule_AddIntConstant(m, "BLOCK_SIZE", 40);
    PyModule_AddIntConstant(m, "NUM_COILS", 8);
    PyModule_AddIntConstant(m, "INTERVAL_MICROSEC", 250);
}
*/

PyMODINIT_FUNC PyInit_ECB(void) {
    Py_Initialize();
    return PyModule_Create(&ecb_definition);
}
