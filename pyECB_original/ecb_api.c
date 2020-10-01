/**
 * \file		ecb_api.c
 * \author		Pantec Engineering AG
 * \author		biscchr
 *
 * Main source file for accessing functions for ECB communication
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ecb_func.h"
#include "ecb_api.h"
#ifdef _USRDLL
#include <windows.h>
#endif

#ifdef _USRDLL
BOOL WINAPI DllMain(HINSTANCE hinDLL, DWORD dwFunction, LPVOID lpNot)
{
	switch(dwFunction)
	{
	case DLL_PROCESS_ATTACH:
		break;
	case DLL_THREAD_ATTACH:
		break;
	case DLL_PROCESS_DETACH:
		break;
	case DLL_THREAD_DETACH:
		break;
	default:
		break;
	}
	return TRUE;
}
#endif

/**
 * \brief	Initializes the ECB API and the communication to the ECB
 * \note	This function needs to be called at the very beginning
 * \param	servername	ip address of the ECB
 * \param	port		port for the ECB communication
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t initECBApi(char* servername, char* port)
{
	unsigned int coilNr;
	int errorNo = 0;
	int desCurrent [NUMCOILS*NUMCURPNTS];


	// initialize ecbState structure
	ecbState.actTempSensorMask = TEMPSENSORALL;
	ecbState.maxTemp = 50;
	ecbState.ecbStatus = ECBSTATUS_UNINIT;
	ecbState.heartbeatTimeout = 0;
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
	{
		ecbState.coilState[coilNr].coilStatus = COIL_OK;
		ecbState.coilState[coilNr].coilTemp = 0;
		ecbState.coilState[coilNr].actCurrent = 0;
		ecbState.coilState[coilNr].hall = 0;
		ecbState.coilState[coilNr].drvTemp = 0;
		ecbState.coilState[coilNr].maxCurrent = MAXCUR_ALLOWED;
	}
	memset(&desCurrent[0], 0, sizeof(int)*NUMCOILS*NUMCURPNTS);
	ecbState.desCurrents = &desCurrent[0];

	errorNo = connectECB(servername, port);
	if(errorNo != 0)
		return ERROR_ECB_SOCKET;

	ecbState.ecbStatus = ECBSTATUS_INIT;

	// first get the revision from the ECB and check the value
	errorNo = getRevision();
	if(errorNo != 0)
		return ERROR_ECB_GETREVISION;
	if(ECB_REVISION_MIN > ecbState.revision)
		return ERROR_ECB_WRONGREVISION;

	// set initial values on ECB
	errorNo = setECBMaxTemp();
	if(errorNo != 0)
		return ERROR_ECB_SETMAXTEMP;
	errorNo = setECBTempMask();
	if(errorNo != 0)
		return ERROR_ECB_SETTEMPMSK;
	errorNo = setECBCurrents(0);
	if(errorNo != 0)
		return ERROR_ECB_SETCURRENT;
	errorNo = setECBMaxCurrents();
	if(errorNo != 0)
		return ERROR_ECB_SETMAXCUR;
	errorNo = setECBTimeout();
	if(errorNo != 0)
		return ERROR_ECB_SETTIMEOUT;

	// clear all open events by disabling the currents
	if(ERROR_NOERROR != disableECBCurrents())
		return ERROR_ECB_CUROFF;

	// get actual drive temperature
	errorNo = getECBDrvTemp();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;

	// get actual coil status
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;

	updateECBStatus();

	return ERROR_NOERROR;
}

/**
 * \brief	Exits the ECB API and closes the ECB communication
 * \note	This function needs to be called at the very end
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t exitECBApi()
{
	int errorNo;
	if(checkECBSocket() < 0)
		return ERROR_ECB_SOCKET;
	if(ERROR_NOERROR != disableECBCurrents())
		return ERROR_ECB_CUROFF;
	errorNo = disconnectECB();
	if(0 != errorNo)
		return ERROR_ECB_SOCKET;
	return ERROR_NOERROR;
}

/**
 * \brief	Sets the desired currents of the ECB
 * \param	currents	pointer to the desired currents buffer for all coils
 * \param	direct		flag for applying the values directly without appending to existing buffer
 * \note	The desired current is given in mA
 * \note	The parameter current must point to a memory array
			containing 8*40 integer current samples in mA.
			current(i,j) has the array index (i+j*8)
			where i is the coil index [0,7] and 
			j is the desired value index [0,39]
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t setDesCurrents(int* currents, unsigned char direct)
{
	int errorNo;
	int coilNr;
	int pnts;
	ecb_api_error_t apiErrorNo = ERROR_NOERROR;
	unsigned int ecbStatus = 0;

	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;

	// check for slew rate or max current exceeding
	for(pnts = 0; pnts < NUMCURPNTS; pnts++)
		for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
			if((unsigned int)abs(currents[coilNr+pnts*NUMCOILS]) > ecbState.coilState[coilNr].maxCurrent)
				return ERROR_MAXCUR_EXCEEDED;

	// values are OK, so copy the pointer into the ecb struct
	ecbState.desCurrents = currents;

	// set the currents on the ECB when buffer is free
	do
	{
		errorNo = setECBCurrents(direct);
		// buffer is full, check if ecbStatus is OK
		if(DESBUF_FULL == errorNo)
		{
			apiErrorNo = getECBStatus(&ecbStatus);
			if(ERROR_NOERROR != apiErrorNo)
				return apiErrorNo;
			if(ecbStatus != 0x3)
				return ERROR_ECB_NOTENABLED;
		}
	} while(DESBUF_FULL == errorNo);


	if(0 != errorNo)
		return ERROR_ECB_SETCURRENT;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the actual currents of the ECB
 * \param	currents	pointer where the actual currents of all coils get stored
 * \note	The actual current is given in mA
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getActCurrents(int* currents)
{
	int coilNr;
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	updateECBStatus();
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
		currents[coilNr] = ecbState.coilState[coilNr].actCurrent;
	return ERROR_NOERROR;
}

/**
 * \brief	Sets the timeout of the ECB heartbeat feature
 * \param	timeout		pointer to the timeout value
 * \note	If the timeout is set to 0, the heartbeat feature gets disabled,
 * 			otherwise the ECB turns off the current controllers if no communication is
 * 			performed after timeout ms
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t setHeartbeatTimeout(unsigned int* timeout)
{
	int errorNo;
	ecbState.heartbeatTimeout = *timeout;
	errorNo = setECBTimeout();
	if(0 != errorNo)
		return ERROR_ECB_SETTIMEOUT;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the timeout value of the ECB heartbeat feature
 * \param	timeout		pointer where the timeout value gets stored
 * \note	If the timeout is set to 0, the heartbeat feature gets disabled,
 * 			otherwise the ECB turns off the current controllers if no communication is
 * 			performed after timeout ms
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getHeartbeatTimeout(unsigned int* timeout)
{
	int errorNo;
	errorNo = getECBTimeout();
	if(0 != errorNo)
		return ERROR_ECB_GETTIMEOUT;
	*timeout = ecbState.heartbeatTimeout;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the revision of the connected ECB
 * \param	revision		pointer where the revision value gets stored
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getECBRevision(unsigned int* revision)
{
	int errorNo;
	errorNo = getRevision();
	if(0 != errorNo)
		return ERROR_ECB_GETREVISION;
	*revision = ecbState.revision;
	if(ECB_REVISION_MIN > ecbState.revision)
		return ERROR_ECB_WRONGREVISION;
	return ERROR_NOERROR;
}

/**
 * \brief	Sets the maximum currents of all ECB coils
 * \param	maxCurrent	pointer to the maximum currents for all coils
 * \note	The maximum current is given in mA and is limited by 20A due to ECB hardware.
 * 			A higher current than this value will lead to a over current event and disables the
 * 			ECB current controllers.
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t setMaxCurrent(unsigned int* maxCurrent)
{
	int errorNo;
	int coilNr;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
		if(maxCurrent[coilNr] > MAXCUR_ALLOWED)
			return ERROR_MAXCUR_NOTALLOWED;
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
		ecbState.coilState[coilNr].maxCurrent = maxCurrent[coilNr];
	errorNo = setECBMaxCurrents();
	if(0 != errorNo)
		return ERROR_ECB_SETMAXCUR;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the maximum currents of all ECB coils
 * \param	maxCurrent	pointer where the maximum currents of all coils get stored
 * \note	The maximum current is given in mA
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getMaxCurrent(unsigned int* maxCurrent)
{
	int coilNr;
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
		maxCurrent[coilNr] = ecbState.coilState[coilNr].maxCurrent;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the temperatures of all ECB coils
 * \param	coilTemp	pointer where the temperatures values of all coils get stored
 * \note	The coil temperatures are given in degree Celsius
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getCoilTemp(unsigned int* coilTemp)
{
	int coilNr;
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	updateECBStatus();
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
		coilTemp[coilNr] = ecbState.coilState[coilNr].coilTemp;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the temperatures of all ECB drives
 * \param	drvTemp	pointer where the temperatures values of all drives get stored
 * \note	The drive temperatures are given in degree Celsius
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getDriveTemp(unsigned int* drvTemp)
{
	int coilNr;
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	errorNo = getECBDrvTemp();
	if(0 != errorNo)
		return ERROR_ECB_GETDRVTEMP;
	updateECBStatus();
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
		drvTemp[coilNr] = ecbState.coilState[coilNr].drvTemp;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the Hall values of all ECB coils
 * \param	hall	pointer where the Hall values of all coils get stored
 * \note	The Hall values are given in uT
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getHall(int* hall)
{
	int coilNr;
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	updateECBStatus();
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
		hall[coilNr] = ecbState.coilState[coilNr].hall;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the status of all ECB coils
 * \param	coilStatus	pointer where the status values of all coils get stored
 * \note	The coil status values are given in \ref coil_status
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getCoilStatus(unsigned int* coilStatus)
{
	int coilNr;
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	updateECBStatus();
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
		coilStatus[coilNr] = ecbState.coilState[coilNr].coilStatus;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the ECB coils values temperature, Hall, actual current and the status of the ECB
 * \param	coilTemp	pointer where the temperatures of all coils get stored
 * \param	hall		pointer where the Hall values of all coils get stored
 * \param	currents	pointer where the actual ciol current values of all coils get stored
 * \param	ecbStatus	pointer where the ECB status gets stored
 * \note	The coil temperatures are given in degree Celsius
 * \note	The Hall values are given in uT
 * \note	The coil currents are given in mA
 * \note	The ECB status values are given in \ref ecb_status
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t getCoilValues(unsigned int* coilTemp, int* hall, int* currents, unsigned int* ecbStatus)
{
	int coilNr;
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	updateECBStatus();
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
	{
		coilTemp[coilNr] = ecbState.coilState[coilNr].coilTemp;
		hall[coilNr] = ecbState.coilState[coilNr].hall;
		currents[coilNr] = ecbState.coilState[coilNr].actCurrent;
	}
	*ecbStatus = ecbState.ecbStatus;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the maximal temperature of the ECB coils
 * \return	Maximal temperature
 */
ECBDLL unsigned int getMaxTemp()
{
	return ecbState.maxTemp;
}

/**
 * \brief	Sets the maximal temperature of the ECB coils
 * \note	The temperature is given in degree Celsius
 * 			A higher coil temperature than this value will lead to an over temperature
 * 			event and disables the ECB current controllers.
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t setMaxTemp(unsigned int maxTemp)
{
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	ecbState.maxTemp = maxTemp;
	errorNo = setECBMaxTemp();
	if(0 != errorNo)
		return ERROR_ECB_SETMAXTEMP;
	return ERROR_NOERROR;
}

/**
 * \brief	Returns the mask for enabling the coil temperature sensors
 * \note	The temperature sensor mask are defined by \ref temp_sensor_mask.
 * 			The LSB is the first coil, whereas the MSB is the last coil
 * 			1 means that the sensor is enabled, 0 disabled the sensor.
 * \return	Actual temperature sensor mask
 */
ECBDLL unsigned char getTempMask()
{
	return ecbState.actTempSensorMask;
}

/**
 * \brief	Sets the mask for enabling the coil temperature sensors
 * \note	The temperature sensor mask are defined by \ref temp_sensor_mask.
 * 			The LSB is the first coil, whereas the MSB is the last coil
 * 			1 means that the sensor is enabled, 0 disabled the sensor.
 * \param	tempMask	Temperature sensor mask
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t setTempMask(char tempMask)
{
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	ecbState.actTempSensorMask = tempMask;
	errorNo = setECBTempMask();
	if(0 != errorNo)
		return ERROR_ECB_SETTEMPMSK;
	return ERROR_NOERROR;
}

/**
 * \brief	Evaluates the actual status of the ECB
 * \param	ecbStatus	pointer where the ECB status gets stored
 * \return	\ref error_codes
 * \note	The ECB status values are given in \ref ecb_status
 */
ECBDLL ecb_api_error_t getECBStatus(unsigned int* ecbStatus)
{
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	updateECBStatus();
	*ecbStatus = ecbState.ecbStatus;
	return ERROR_NOERROR;
}

/**
 * \brief	Enables the ECB currents on all coils
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t enableECBCurrents()
{
	int errorNo;
	int desCurrent [NUMCOILS*NUMCURPNTS];

	memset(&desCurrent[0], 0, sizeof(int)*NUMCOILS*NUMCURPNTS);
	ecbState.desCurrents = &desCurrent[0];

	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;

	// update the ECB status
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	updateECBStatus();
	// check if ECB is OK
	if(ecbState.ecbStatus != ECBSTATUS_INIT)
		return ERROR_ECB_COILEVENT;

	// set all desired currents to zero, otherwise current-step possible
	errorNo = setECBCurrents(0);
	if(0 != errorNo)
		return ERROR_ECB_SETCURRENT;

	// enable the ECB current controllers
	errorNo = setECBCurrentOn();
	if(0 != errorNo)
		return ERROR_ECB_CURON;

	// update the status and set current on flag
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	ecbState.ecbStatus |= ECBSTATUS_CURON;
	updateECBStatus();

	return ERROR_NOERROR;
}

/**
 * \brief	Disables the ECB currents on all coils and clears all active errors
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t disableECBCurrents()
{
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;

	errorNo = setECBCurrentOff();
	if(0 != errorNo)
		return ERROR_ECB_CUROFF;
	errorNo = getECBCoilStatus();
	if(0 != errorNo)
		return ERROR_ECB_GETCOILSTATUS;
	ecbState.ecbStatus &= ~ECBSTATUS_CURON;
	updateECBStatus();
	return ERROR_NOERROR;
}

/**
 * \brief	Factory reset of the ECB, sets all ECB parameters to factory default
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t factoryResetECB()
{
	int errorNo;
	if(ERROR_NOERROR != checkECBComm())
		return ERROR_ECB_UNINITIALIZED;

	errorNo = resetECB();
	if(0 != errorNo)
		return ERROR_ECB_RESET;
	return ERROR_NOERROR;
}

/**
 * \brief	Updates the status of the ECB depending on the coil status
 */
ECBDLL void updateECBStatus()
{
	int coilNr;
	int errorFound;
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
	{

		// check if event is only coil over temperature and error is masked out
		if(((ecbState.coilState[coilNr].coilStatus & COIL_ERR_MSK) == COIL_OTCOIL) && ((ecbState.actTempSensorMask & (1<<coilNr)) == 0))
			errorFound = 0;
		// check if other error event on coil is present
		else if((ecbState.coilState[coilNr].coilStatus & COIL_ERR_MSK) != COIL_OK)
			errorFound = 1;
		else
			errorFound = 0;

		if(errorFound != 0)
		{
			ecbState.ecbStatus |= (ECBSTATUS_ERR_COIL0 << coilNr);
			ecbState.ecbStatus &= ~ECBSTATUS_CURON;
		}
		else
			ecbState.ecbStatus &= ~(ECBSTATUS_ERR_COIL0 << coilNr);
	}
}

/**
 * \brief	Checks the status of the ECB communication
 * \return	\ref error_codes
 */
ECBDLL ecb_api_error_t checkECBComm()
{
	if(!(ecbState.ecbStatus & ECBSTATUS_INIT))
		return ERROR_ECB_UNINITIALIZED;
	return ERROR_NOERROR;
}

/**
 * \brief	Debugging function for printing all currently stored values
 */
ECBDLL void printECBStruct()
{
	int coilNr;
	printf("Max Temp:\t\t%d\n", ecbState.maxTemp);
	printf("Temp Sensor Mask:\t0x%02x\n", ecbState.actTempSensorMask);
	printf("Heartbeat timeout:\t%04x\n", ecbState.heartbeatTimeout);
	printf("ECB Status:\t\t0x%02x\n", ecbState.ecbStatus);
	for(coilNr = 0; coilNr < NUMCOILS; coilNr++)
	{
		printf("\nCoil Number %d:\n", coilNr);
		printf("\tAct Current:\t%d\n", ecbState.coilState[coilNr].actCurrent);
		printf("\tMax Current:\t%d\n", ecbState.coilState[coilNr].maxCurrent);
		printf("\tCoil Temp:\t%d\n", ecbState.coilState[coilNr].coilTemp);
		printf("\tHall:\t\t%d\n", ecbState.coilState[coilNr].hall);
		printf("\tCoil Status:\t%d\n", ecbState.coilState[coilNr].coilStatus);
	}
}

ECBDLL const char* errorToString(ecb_api_error_t e)
{
    switch (e)
    {
    case ERROR_NOERROR:             return "0, Successful";
    case ERROR_ECB_UNINITIALIZED:   return "1, ECB uninitialized";
    case ERROR_ECB_SOCKET:          return "2, No valid ECB socket";
    case ERROR_ECB_CURON:           return "3, ECB current on communication error";
    case ERROR_ECB_CUROFF:          return "4, ECB current off communication error";
    case ERROR_ECB_GETREVISION:     return "5, ECB revision communication error";
    case ERROR_ECB_WRONGREVISION:   return "6, ECB revision value error";
    case ERROR_ECB_SETCURRENT:      return "20, ECB set current communication error";
    case ERROR_ECB_GETCURRENT:      return "21, ECB get current communication error";
    case ERROR_MAXCUR_EXCEEDED:     return "22, ECB set desired current, maximum current exceeded";
    case ERROR_MAXCUR_NOTALLOWED:   return "23, ECB set maximum current, value too high";
    case ERROR_ECB_SETMAXCUR:       return "24, ECB set maximal current communication error";
    case ERROR_ECB_GETDRVTEMP:      return "26, ECB get drive temperature communication error";
    case ERROR_ECB_NOTENABLED:      return "27, ECB cannot take new current values, since not enabled";
    case ERROR_ECB_GETCOILSTATUS:   return "30, ECB get coil status communication error";
    case ERROR_ECB_COILEVENT:       return "31, ECB coil has event";
    case ERROR_ECB_SETMAXTEMP:      return "40, ECB set maximal temperature communication error";
    case ERROR_ECB_SETTEMPMSK:      return "50, ECB set temperature sensor mask communication error";
    case ERROR_ECB_RESET:           return "60, ECB factory reset communication error";
    case ERROR_ECB_SETTIMEOUT:      return "70, ECB set hearbeat timeout communication error";
    case ERROR_ECB_GETTIMEOUT:      return "71, ECB get hearbeat timeout communication error";
    default:                        return "[Unknown Error]";
    }
}
