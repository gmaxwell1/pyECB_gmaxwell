/**
 * \file		ecb_socket_win.c
 * \author		Pantec Engineering AG
 * \author		biscchr
 *
 * Source file for the socket handling (for Windows systems)
 * as well as for sending and receiving messages
 */

#ifndef ECBLINUX
#include <stdio.h>
#include <winsock2.h>
#include "ecb_socket_win.h"


SOCKET sock;


/**
 * \brief	Creates an TCP/IP socket
 * \param	servername	ip address of the ECB
 * \param	port		port for the ECB communication
 * \return				0 on success, otherwise error code
 */
int openECBConnection(char* servername, char* port)
{
	int portno;
	int errorNo;
	WSADATA wsaData;
	struct sockaddr_in serv_addr;
	struct hostent *server;

	sock = INVALID_SOCKET;
    portno = atoi(port);
	if(WSAStartup(MAKEWORD(2,2), &wsaData) != NO_ERROR)
		return -1;
	sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
	if (sock == INVALID_SOCKET)
    	return -1;
    server = gethostbyname(servername);
    if (server == NULL)
        return -1;
    memset((void *) &serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
	serv_addr.sin_port = htons(7070);
	memcpy((void *)&serv_addr.sin_addr.s_addr, (void *)server->h_addr, server->h_length);
    serv_addr.sin_port = htons(portno);
	//function call hangs here
	errorNo = connect(sock, (SOCKADDR*) &serv_addr,sizeof(serv_addr));
    if (SOCKET_ERROR == errorNo){
		printf("Error creating socket: %d\n",WSAGetLastError());
        return -1;
	}
    return 0;
}

/**
 * \brief	Closes an TCP/IP socket
 * \return	0 on success, otherwise error code
 */
int closeECBConnection()
{
	return closesocket(sock);
}

/**
 * \brief	Sends a message over a TCP/IP socket
 * \param	msg		pointer to the message buffer
 * \param	length	number of bytes to send
 * \return				0 on success, otherwise error code
 */
int sendECBMessage(char* msg, unsigned int length)
{
	return send(sock, msg, length, 0);
}

/**
 * \brief	Reads a message from a TCP/IP socket
 * \param	msg		pointer to the message buffer
 * \param	length	number of bytes to read
 * \return				0 on success, otherwise error code
 */
int readECBMessage(char* msg, unsigned int length)
{
	return recv(sock, msg, length, 0);
}


/**
 * \brief	Check if TCP/IP socket is valid
 * \return	0 is socket is valid, otherwise -1
 */
int checkSocket()
{
	if(sock == INVALID_SOCKET)
		return -1;
	else
		return 0;
}

#endif /* ECBLINUX */
