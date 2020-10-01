/**
 * \file		ecb_msg.h
 * \author		Pantec Engineering AG
 * \author		biscchr
 *
 * Header file for ECB message definitions
 */

#ifndef ECB_MSG_H_
#define ECB_MSG_H_

#ifdef __cplusplus
extern "C" {
#endif

// Abort transfer header
# define MSG_ABORT_HEAD				0x80

// Common message parameters
#define MSG_INDLSB					0x03	// special ECB object 0x2103
#define MSG_INDMSB					0x21	// special ECB object 0x2103

// message header ccs & 000 & e & s
#define MSG_SPEC_GET_HEAD			0x42	// ccs=10, e=1, s=0
#define MSG_SPEC_SET_HEAD			0x22	// ccs=01, e=1, s=0
#define MSG_DEF_GET_HEAD			0x43	// ccs=10, e=1, s=1
#define MSG_DEF_SET_HEAD			0x23	// ccs=01, e=1, s=1


// Current set/get message
#define MSG_CUR_SUBIND				0x01	// Current set/get
#define MSG_CUR_DIRECT_SUBIND		0x0C	// Current set direct
#define MSG_CUR_DATABYTES			32		// 32 Bytes data
#define MSG_SETCUR_REQSIZE			40
#define MSG_SETCUR_RESPSIZE			8


// MaxCurrent set/get message
#define MSG_MAXCUR_SUBIND			0x02	// MaxCurrent set/get
#define MSG_MAXCUR_DATABYTES		32		// 32 Bytes data
#define MSG_MAXCUR_REQSIZE			40
#define MSG_MAXCUR_RESPSIZE			8

// Slewrate set/get message
#define MSG_HEARTBEAT_SUBIND		0x09	// Heartbeat set/get
#define MSG_HEARTBEAT_REQSIZE		8
#define MSG_HEARTBEAT_RESPSIZE		8

// MaxTemp set/get message
#define MSG_MAXTEMP_SUBIND			0x03	// MaxTemp set/get
#define MSG_MAXTEMP_REQSIZE			8
#define MSG_MAXTEMP_RESPSIZE		8

// TempSensorMask set/get message
#define MSG_TEMPMSK_SUBIND			0x04	// TempSensorMask set/get
#define MSG_TEMPMSK_REQSIZE			8
#define MSG_TEMPMSK_RESPSIZE		8

// DriveTemp get message
#define MSG_DRVTEMP_SUBIND			0x0A	// DriveTemp get
#define MSG_DRVTEMP_DATABYTES		32		// 32 Bytes data
#define MSG_GETDRVTEMP_REQSIZE		8
#define MSG_GETDRVTEMP_RESPSIZE		40

// Revision get message
#define MSG_REV_SUBIND				0x0B	// Revision get
#define MSG_REV_REQSIZE				8
#define MSG_REV_RESPSIZE			8

// Current on message
#define MSG_CURON_SUBIND			0x05	// Enable ECB
#define MSG_CURON_REQSIZE			8
#define MSG_CURON_RESPSIZE			8

// Current off message
#define MSG_CUROFF_SUBIND			0x06	// Disable ECB
#define MSG_CUROFF_REQSIZE			8
#define MSG_CUROFF_RESPSIZE			8

// Coil status message
#define MSG_COILSTATUS_SUBIND		0x07	// Coil values
#define MSG_COILSTATUS_REQSIZE		8
#define MSG_COILSTATUS_RESPSIZE		136

// Factory reset message
#define MSG_RESET_INDLSB			0x01	// SYSTEM_PARAMS 0x2001
#define MSG_RESET_INDMSB			0x20	// SYSTEM_PARAMS 0x2001
#define MSG_RESET_LD_SUBIND			0x01	// LOAD_DEFAULT
#define MSG_RESET_REQSIZE			8
#define MSG_RESET_RESPSIZE			8

#define MSG_RESET_SV_SUBIND			0x03	// SAVE_FLASH

// message header factory reset ccs & 000 & e & s
#define MSG_RESET_HEAD				0x23	// ccs=01, e=1, s=1


#ifdef __cplusplus
}
#endif

#endif /* ECB_MSG_H_ */
