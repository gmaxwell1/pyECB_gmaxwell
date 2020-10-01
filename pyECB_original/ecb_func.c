/**
 * \file		ecb_func.c
 * \author		Pantec Engineering AG
 * \author		biscchr
 *
 * Source file for internal ECB functions
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if ECBLINUX == 1
#include "ecb_socket.h"
#else
#include "ecb_socket_win.h"
#endif

#include "ecb_msg.h"
#include "ecb_func.h"


ecb_state_t ecbState;

/// ECB Message buffer
static unsigned char msgBuf[MSGBUFSIZE];

/**
 * \brief	Sets the desired currents in the ECB
 * \return	0 on success, -1 on error
 */
int setECBCurrents(unsigned char direct)
{
	unsigned int cnt;
	unsigned int dataLength;
	unsigned int packetLength;

	dataLength = NUMCURPNTS*NUMCOILS*sizeof(int);
	packetLength = dataLength + 8;

	// construct message
	msgBuf[0] = MSG_SPEC_SET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	// use direct message when requested
	if(direct > 0)
		msgBuf[3] = MSG_CUR_DIRECT_SUBIND;
	else
		msgBuf[3] = MSG_CUR_SUBIND;

	// set the data length
	memcpy(&msgBuf[4], &dataLength, sizeof(int));

	for(cnt = 8; cnt < packetLength; cnt++)
		msgBuf[cnt] = 0x00;

	// load desired values into message buffer
	memcpy(&(msgBuf[8]), ecbState.desCurrents, sizeof(int)*NUMCOILS*NUMCURPNTS);

	// communicate with ECB
	if(transmitECBCommand(packetLength, MSG_SETCUR_RESPSIZE) < 0)
		return -1;

	return (unsigned int) msgBuf[4];
}

/**
 * \brief	Sets the maximal currents in the ECB
 * \return	0 on success, -1 on error
 */
int setECBMaxCurrents()
{
	int cnt;

	// construct message
	msgBuf[0] = MSG_SPEC_SET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_MAXCUR_SUBIND;
	msgBuf[4] = MSG_MAXCUR_DATABYTES;
	for(cnt = 5; cnt < MSG_MAXCUR_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// load desired values into message buffer
	for(cnt = 0; cnt < NUMCOILS; cnt++)
		memcpy(&(msgBuf[8+cnt*4]), &(ecbState.coilState[cnt].maxCurrent), sizeof(int));

	// communicate with ECB
	if(transmitECBCommand(MSG_MAXCUR_REQSIZE, MSG_MAXCUR_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	return 0;
}

/**
 * \brief	Sets the timeout of the hearbeat feature
 * \return	0 on success, -1 on error
 */
int setECBTimeout()
{
	int cnt;

	// construct message
	msgBuf[0] = MSG_DEF_SET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_HEARTBEAT_SUBIND;
	for(cnt = 4; cnt < MSG_HEARTBEAT_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// load timeout value into message buffer
	memcpy(&(msgBuf[4]), &(ecbState.heartbeatTimeout), sizeof(unsigned int));

	// communicate with ECB
	if(transmitECBCommand(MSG_HEARTBEAT_REQSIZE, MSG_HEARTBEAT_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	return 0;
}

/**
 * \brief	Gets the timeout of the hearbeat feature
 * \return	0 on success, -1 on error
 */
int getECBTimeout()
{
	int cnt;

	// construct message
	msgBuf[0] = MSG_DEF_GET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_HEARTBEAT_SUBIND;
	for(cnt = 4; cnt < MSG_HEARTBEAT_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// communicate with ECB
	if(transmitECBCommand(MSG_HEARTBEAT_REQSIZE, MSG_HEARTBEAT_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	// load timeout value from message buffer
	memcpy(&(ecbState.heartbeatTimeout), &(msgBuf[4]), sizeof(unsigned int));

	return 0;
}

/**
 * \brief	Sets the maximal temperature in the ECB
 * \return	0 on success, -1 on error
 */
int setECBMaxTemp()
{
	int cnt;

	// construct message
	msgBuf[0] = MSG_DEF_SET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_MAXTEMP_SUBIND;
	for(cnt = 4; cnt < MSG_MAXTEMP_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// load desired value into message buffer
	memcpy(&(msgBuf[4]), &(ecbState.maxTemp), sizeof(int));

	// communicate with ECB
	if(transmitECBCommand(MSG_MAXTEMP_REQSIZE, MSG_MAXTEMP_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	return 0;
}

/**
 * \brief	Sets the temperature sensor mask in the ECB
 * \return	0 on success, -1 on error
 */
int setECBTempMask()
{
	int cnt;

	// construct message
	msgBuf[0] = MSG_DEF_SET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_TEMPMSK_SUBIND;
	for(cnt = 4; cnt < MSG_TEMPMSK_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// load desired value into message buffer
	memcpy(&(msgBuf[4]), &(ecbState.actTempSensorMask), sizeof(char));

	// communicate with ECB
	if(transmitECBCommand(MSG_TEMPMSK_REQSIZE, MSG_TEMPMSK_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	return 0;
}

/**
 * \brief	Reads the actual drive temperatures from the ECB and stores
 * 			them in the ECB State struct
 * \return	0 on success, -1 on error
 */
int getECBDrvTemp()
{
	int cnt;

	// construct message
	msgBuf[0] = MSG_SPEC_GET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_DRVTEMP_SUBIND;
	for(cnt = 4; cnt < MSG_GETDRVTEMP_REQSIZE; cnt++)
		msgBuf[cnt] = 0x0;

	// communicate with ECB
	if(transmitECBCommand(MSG_GETDRVTEMP_REQSIZE, MSG_GETDRVTEMP_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	// save received values
	for(cnt = 0; cnt < NUMCOILS; cnt++)
		memcpy(&(ecbState.coilState[cnt].drvTemp), &(msgBuf[8+cnt*4]), sizeof(int));
	return 0;
}

/**
 * \brief	Reads the revision from the ECB and stores it in the ECB State struct
 * \return	0 on success, -1 on error
 */
int getRevision()
{
	int cnt;

	// construct message
	msgBuf[0] = MSG_DEF_GET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_REV_SUBIND;
	for(cnt = 4; cnt < MSG_REV_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// communicate with ECB
	if(transmitECBCommand(MSG_REV_REQSIZE, MSG_REV_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	// load timeout value from message buffer
	memcpy(&(ecbState.revision), &(msgBuf[4]), sizeof(unsigned int));

	return 0;
}

/**
 * \brief	Reads the actual coil status from the ECB and stores
 * 			them in the ECB State struct
 * \return	0 on success, -1 on error
 */
int getECBCoilStatus()
{
	int cnt;

	// construct message
	msgBuf[0] = MSG_SPEC_GET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_COILSTATUS_SUBIND;
	for(cnt = 4; cnt < MSG_COILSTATUS_REQSIZE; cnt++)
		msgBuf[cnt] = 0x0;

	// communicate with ECB
	if(transmitECBCommand(MSG_COILSTATUS_REQSIZE, MSG_COILSTATUS_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	// save received values in the struct
	for(cnt = 0; cnt < NUMCOILS; cnt++)
	{
		memcpy(&(ecbState.coilState[cnt].hall), &(msgBuf[8+cnt*4]), sizeof(int));
		memcpy(&(ecbState.coilState[cnt].coilTemp), &(msgBuf[8+32+cnt*4]), sizeof(int));
		memcpy(&(ecbState.coilState[cnt].actCurrent), &(msgBuf[8+64+cnt*4]), sizeof(int));
		memcpy(&(ecbState.coilState[cnt].coilStatus), &(msgBuf[8+96+cnt*4]), sizeof(int));
	}

	return 0;
}

/**
 * \brief	Enables the ECB currents, enables all coil controllers
 * \return	0 on success, -1 on error
 */
int setECBCurrentOn()
{
	int cnt;
	// construct message
	msgBuf[0] = MSG_DEF_SET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_CURON_SUBIND;
	for(cnt = 4; cnt < MSG_CURON_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// communicate with ECB
	if(transmitECBCommand(MSG_CURON_REQSIZE, MSG_CURON_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	return 0;
}

/**
 * \brief	Disables the ECB currents, disables all coil controllers and clear all error events
 * \return	0 on success, -1 on error
 */
int setECBCurrentOff()
{
	int cnt;
	// construct message
	msgBuf[0] = MSG_DEF_SET_HEAD;
	msgBuf[1] = MSG_INDLSB;
	msgBuf[2] = MSG_INDMSB;
	msgBuf[3] = MSG_CUROFF_SUBIND;
	for(cnt = 4; cnt < MSG_CUROFF_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// communicate with ECB
	if(transmitECBCommand(MSG_CUROFF_REQSIZE, MSG_CUROFF_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	return 0;
}

/**
 * \brief	Resets the ECB to factory default parameters
 * \return	0 on success, -1 on error
 */
int resetECB()
{
	int cnt;
	// construct message
	msgBuf[0] = MSG_RESET_HEAD;
	msgBuf[1] = MSG_RESET_INDLSB;
	msgBuf[2] = MSG_RESET_INDMSB;
	msgBuf[3] = MSG_RESET_LD_SUBIND;
	for(cnt = 4; cnt < MSG_RESET_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// communicate with ECB
	if(transmitECBCommand(MSG_RESET_REQSIZE, MSG_RESET_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	// save defaults in flash
	msgBuf[0] = MSG_RESET_HEAD;
	msgBuf[1] = MSG_RESET_INDLSB;
	msgBuf[2] = MSG_RESET_INDMSB;
	msgBuf[3] = MSG_RESET_SV_SUBIND;
	for(cnt = 4; cnt < MSG_RESET_REQSIZE; cnt++)
		msgBuf[cnt] = 0x00;

	// communicate with ECB
	if(transmitECBCommand(MSG_RESET_REQSIZE, MSG_RESET_RESPSIZE) < 0)
		return -1;

	// check for abort transfer response
	if(MSG_ABORT_HEAD == msgBuf[0])
		return -1;

	return 0;
}


/**
 * \brief	Transmitting the message buffer to the ECB
 * \param	txlength	Length of the transmit message in bytes
 * \param	rxlength	Length of the received message in bytes
 * \return	0 on success, -1 on error
 */
int transmitECBCommand(unsigned int txlength, unsigned int rxlength)
{
	// communicate with ECB
	if(sendECBMessage((char*)&msgBuf[0], txlength) < 0)
		return -1;
	if(readECBMessage((char*)&msgBuf[0], rxlength) < 0)
		return -1;
	return 0;
}

/**
 * \brief	Connect with the ECB
 * \param	servername	ip address of the ECB
 * \param	port		port for the ECB communication
 * \return	0 on success, -1 on error
 */
int connectECB(char* servername, char* port)
{
	return openECBConnection(servername, port);
}

/**
 * \brief	Disonnect with the ECB
 * \return	0 on success, -1 on error
 */
int disconnectECB()
{
	return closeECBConnection();
}


/**
 * \brief	Check the status of the TCP/IP socket
 * \return	0 on success, -1 on error
 */
int checkECBSocket()
{
	return checkSocket();
}
