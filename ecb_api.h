/**
 * \file		ecb_api.h
 * \author		Pantec Engineering AG
 * \author		biscchr
 *
 * Main header file for accessing functions for ECB communication
 */

/**
 * \mainpage ECB Host API
 *

 * \section sec_desc	ECB Host API Description
 *
 * \n The ECB Host API is a software library that is used to communicate with the Electronic Control Box (ECB) by Pantec Engineering AG.\n
 *
 * \n This documentation of the ECB Host API describes all necessary functions and definitions for interfacing\n
 * the ECB Host API and therefore communicate with the Electronic Control Box.
 * \n\n
 * The list and description of the ECB Host API functions can be found in page \ref ecb_api.\n
 * The return values of the ECB Host API functions can be found in page \ref error_codes.\n
 * The masks of the ECB status can be found in page \ref ecb_status.\n
 * The masks of the coil status can be found in page \ref coil_status.\n\n
 *
 *
 * \section sec_hist Revision History
 *
 * Version	|	Date		|	Author	|	Comment
 * ---------|---------------|-----------|----------------------------------------------------------
 * 1.0.00	|	20.06.2014	|	biscchr	|	first release
 * 1.1.00	|	02.07.2014	|	biscchr	|	added slewrate trajectory functionality
 * 1.1.01	|	07.07.2014	|	biscchr	|	ported to Windows systems and supporting DLL generation
 * 1.2.00	|	18.07.2014	|	biscchr	|	heartbeat value can now be set by the ECB Host API
 * 1.2.01	|	29.07.2014	|	biscchr	|	temperature sensor mask bug fixed
 * 3.0.00	|	15.12.2014	|	biscchr	|	desired value buffer communication; drive temperature message; revision message
 * 3.0.01	|	20.01.2015	|	biscchr	|	default heartbeat timeout set to 0 in initECBApi
 * 3.1.00	|	30.01.2015	|	biscchr	|	implemented direct-set of desired current values
 * \n\n\n
 *
 */



#ifndef ECB_API_H_
#define ECB_API_H_

#ifdef _USRDLL
#ifdef ECB_API_DLL_EXPORTS
#define ECBDLL __declspec(dllexport)
#else
#define ECBDLL __declspec(dllimport)
#endif
#else
#define ECBDLL
#endif

#ifdef __cplusplus
extern "C" {
#endif


/**
 * \defgroup ecb_status ECB status
 * @{
 */
#define ECBSTATUS_UNINIT		0x0000	///< ECB uninitialized
#define ECBSTATUS_INIT			0x0001	///< ECB initialized
#define ECBSTATUS_CURON			0x0002	///< ECB coil currents on

#define ECBSTATUS_ERR_COIL0		0x0010	///< ECB error on coil 0
#define ECBSTATUS_ERR_COIL1		0x0020	///< ECB error on coil 1
#define ECBSTATUS_ERR_COIL2		0x0040	///< ECB error on coil 2
#define ECBSTATUS_ERR_COIL3		0x0080	///< ECB error on coil 3
#define ECBSTATUS_ERR_COIL4		0x0100	///< ECB error on coil 4
#define ECBSTATUS_ERR_COIL5		0x0200	///< ECB error on coil 5
#define ECBSTATUS_ERR_COIL6		0x0400	///< ECB error on coil 6
#define ECBSTATUS_ERR_COIL7		0x0800	///< ECB error on coil 7

/** @}*/


/**
 * \defgroup coil_status Coil status
 * @{
 */
#define COIL_OK			0x00000000		///< Coil status OK
#define COIL_OCUR		0x00000001		///< Coil status over current
#define COIL_OVOLT		0x00000002		///< Coil status over voltage
#define COIL_UVOLT		0x00000004		///< Coil status under voltage
#define COIL_OTDRV		0x00000080		///< Coil status over temperature drive
#define COIL_OTCOIL		0x00000100		///< Coil status over temperature coil
#define COIL_EMERG1		0x00000200		///< Coil status emergency stop 1
#define COIL_EMERG2		0x00000400		///< Coil status emergency stop 2
#define COIL_HEARTB		0x00000800		///< Coil status timeout heartbeat
#define COIL_OTHER		0x00200000		///< Coil status error due to error on other coil
/** @}*/
#define COIL_ERR_MSK	0x001FFFFF		// mask for errors on actual coil

/**
 * \defgroup temp_sensor_mask Temperature sensor mask
 * @{
 */
#define TEMPSENSOR0		0x01	///< Mask for sensor 0
#define TEMPSENSOR1		0x02	///< Mask for sensor 1
#define TEMPSENSOR2		0x04	///< Mask for sensor 2
#define TEMPSENSOR3		0x08	///< Mask for sensor 3
#define TEMPSENSOR4		0x10	///< Mask for sensor 4
#define TEMPSENSOR5		0x20	///< Mask for sensor 5
#define TEMPSENSOR6		0x40	///< Mask for sensor 6
#define TEMPSENSOR7		0x80	///< Mask for sensor 7
#define TEMPSENSORALL	0xFF	///< Mask for all sensors
/** @}*/

#define MAXCUR_ALLOWED		19800	///< Maximum allowed current for the ECB in mA

#define ECB_REVISION_MIN	346		///< Minimum revision value of the ECB

/**
 * \defgroup error_codes ECP Host API error codes
 * @{
 */
/// Error return codes of ECB Host API
typedef enum ecb_api_error_e
{
    ERROR_NOERROR					= 0,		///< 0, Successful
    ERROR_ECB_UNINITIALIZED			= 1,		///< 1, ECB uninitialized
    ERROR_ECB_SOCKET				= 2,		///< 2, No valid ECB socket
    ERROR_ECB_CURON					= 3,		///< 3, ECB current on communication error
    ERROR_ECB_CUROFF				= 4,		///< 4, ECB current off communication error
    ERROR_ECB_GETREVISION			= 5,		///< 5, ECB revision communication error
    ERROR_ECB_WRONGREVISION			= 6,		///< 6, ECB revision value error

    ERROR_ECB_SETCURRENT			= 20,		///< 20, ECB set current communication error
    ERROR_ECB_GETCURRENT			= 21,		///< 21, ECB get current communication error
    ERROR_MAXCUR_EXCEEDED			= 22,		///< 22, ECB set desired current, maximum current exceeded
    ERROR_MAXCUR_NOTALLOWED			= 23,		///< 23, ECB set maximum current, value too high
    ERROR_ECB_SETMAXCUR				= 24,		///< 24, ECB set maximal current communication error

    ERROR_ECB_GETDRVTEMP			= 26,		///< 26, ECB get drive temperature communication error
    ERROR_ECB_NOTENABLED			= 27,		///< 27, ECB cannot take new current values, since not enabled

    ERROR_ECB_GETCOILSTATUS			= 30,		///< 30, ECB get coil status communication error
    ERROR_ECB_COILEVENT				= 31,		///< 31, ECB coil has event

    ERROR_ECB_SETMAXTEMP			= 40,		///< 40, ECB set maximal temperature communication error

    ERROR_ECB_SETTEMPMSK			= 50,		///< 50, ECB set temperature sensor mask communication error

    ERROR_ECB_RESET					= 60,		///< 60, ECB factory reset communication error

    ERROR_ECB_SETTIMEOUT			= 70,		///< 70, ECB set hearbeat timeout communication error
    ERROR_ECB_GETTIMEOUT			= 71,		///< 71, ECB get hearbeat timeout communication error

} ecb_api_error_t;
/** @}*/



/**
 * \defgroup ecb_api ECB API functions
 * @{
 */
ECBDLL ecb_api_error_t initECBApi(char* servername, char* port);
ECBDLL ecb_api_error_t exitECBApi();
ECBDLL ecb_api_error_t enableECBCurrents();
ECBDLL ecb_api_error_t disableECBCurrents();
ECBDLL ecb_api_error_t setDesCurrents(int* currents, unsigned char direct);
ECBDLL ecb_api_error_t getActCurrents(int* currents);
ECBDLL ecb_api_error_t setHeartbeatTimeout(unsigned int* timeout);
ECBDLL ecb_api_error_t getHeartbeatTimeout(unsigned int* timeout);
ECBDLL ecb_api_error_t getECBRevision(unsigned int* revision);
ECBDLL ecb_api_error_t setMaxCurrent(unsigned int* maxCurrent);
ECBDLL ecb_api_error_t getMaxCurrent(unsigned int* maxCurrent);
ECBDLL ecb_api_error_t getCoilTemp(unsigned int* coilTemp);
ECBDLL ecb_api_error_t getDriveTemp(unsigned int* drvTemp);
ECBDLL ecb_api_error_t getHall(int* hall);
ECBDLL ecb_api_error_t getCoilStatus(unsigned int* coilStatus);
ECBDLL ecb_api_error_t getCoilValues(unsigned int* coilTemp, int* hall, int* currents, unsigned int* coilStatus);

ECBDLL int getSocketFD();
ECBDLL unsigned int getMaxTemp();
ECBDLL ecb_api_error_t setMaxTemp(unsigned int maxTemp);
ECBDLL unsigned char getTempMask();
ECBDLL ecb_api_error_t setTempMask(char tempMask);
ECBDLL ecb_api_error_t factoryResetECB();
ECBDLL ecb_api_error_t getECBStatus(unsigned int* ecbStatus);
ECBDLL void updateECBStatus();
ECBDLL ecb_api_error_t checkECBComm();



ECBDLL void printECBStruct();

ECBDLL const char* errorToString(ecb_api_error_t e);
/** @}*/

#ifdef _USRDLL
BOOL WINAPI DllMain(HINSTANCE hinDLL, DWORD dwFunction, LPVOID lpNot);
#endif

#ifdef __cplusplus
}
#endif

#endif /* ECB_API_H_ */
