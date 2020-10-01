/**
 * \file		ecb_socket_win.h
 * \author		Pantec Engineering AG
 * \author		biscchr
 *
 * Header file for the socket handling (for Windows systems)
 * as well as for sending and receiving messages
 */

#ifndef ECBLINUX
#ifndef ECB_SOCKET_WIN_H_
#define ECB_SOCKET_WIN_H_

#ifdef __cplusplus
extern "C" {
#endif


int openECBConnection(char* servername, char* port);
int closeECBConnection();
int sendECBMessage(char* msg, unsigned int length);
int readECBMessage(char* msg, unsigned int length);
int checkSocket();

#ifdef __cplusplus
}
#endif

#endif /* ECB_SOCKET_WIN_H_ */
#endif /* ECBLINUX */
