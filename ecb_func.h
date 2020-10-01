/**
 * \file		ecb_func.h
 * \author		Pantec Engineering AG
 * \author		biscchr
 *
 * Header file for internal ECB functions
 */

#ifndef ECB_FUNC_H_
#define ECB_FUNC_H_

#ifndef ECBLINUX
#include <winsock2.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/// Number of coils controlled by the ECB
#define NUMCOILS	8

/// Message buffer size in bytes
#define MSGBUFSIZE	1500

/// Number of current values per coil in one current message
#define NUMCURPNTS	40

/// Full desired value buffer response
#define DESBUF_FULL	0x63


/// Structure holding the current state of a coil
typedef struct ecb_coil_state_s
{
	int				actCurrent;					///< Actual coil current in mA
	unsigned int	maxCurrent;					///< Maximum current in mA
	unsigned int	coilTemp;					///< Actual coil temperature in degree Celsius
	int				hall;						///< Actual Hall value in uT
	unsigned int	drvTemp;					///< Actual drive temperatur in degree Celsius
	unsigned int	coilStatus;					///< Actual coil state

} ecb_coil_state_t;

/// Structure holding the current state of the ECB
typedef struct ecb_state_s
{
	unsigned int		maxTemp;				///< Maximal temperature
	unsigned char		actTempSensorMask;		///< Actual temperature sensor mask
	unsigned int		heartbeatTimeout;		///< Timeout of the heartbeat feature in ms
	unsigned int		ecbStatus;				///< Current ECB status
	unsigned int		revision;				///< Revision of the ECB
	int*				desCurrents;			///< Pointer to the desired current array
	ecb_coil_state_t	coilState[NUMCOILS];	///< Array of coil state structs

} ecb_state_t;

/// Structure which holds the current state of the ECB
extern ecb_state_t ecbState;

// internal API function declarations
int checkECBSocket();
int transmitECBCommand(unsigned int txlength, unsigned int rxlength);
int setECBCurrents(unsigned char direct);
int setECBMaxCurrents();
int setECBTimeout();
int getECBTimeout();
int getRevision();
int getECBCoilStatus();
int setECBMaxTemp();
int setECBTempMask();
int getECBDrvTemp();
int setECBCurrentOn();
int setECBCurrentOff();
int resetECB();
int connectECB(char* servername, char* port);
int disconnectECB();

#ifdef __cplusplus
}
#endif

#endif /* ECB_FUNC_H_ */
