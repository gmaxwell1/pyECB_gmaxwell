/**
 * \file		ecb_socket.c
 * \author		Pantec Engineering AG
 * \author		biscchr
 *
 * Source file for the socket handling (for Linux systems)
 * as well as for sending and receiving messages
 */
#ifdef ECBLINUX
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <netdb.h>
#include <sys/socket.h>
#include <netinet/in.h>

#include "ecb_socket.h"

int sockfd;

/**
 * \brief	Creates an TCP/IP socket
 * \param	servername	ip address of the ECB
 * \param	port		port for the ECB communication
 * \return				0 on success, otherwise error code
 */
int openECBConnection(char* servername, char* port)
{
	int portno;
	struct hostent *server;
	struct sockaddr_in serv_addr;

    portno = atoi(port);
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
    	return -1;
    server = gethostbyname(servername);
    if (server == NULL)
        return -1;
    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *)server->h_addr,
         (char *)&serv_addr.sin_addr.s_addr,
         server->h_length);
    serv_addr.sin_port = htons(portno);
    if (connect(sockfd,(struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0)
        return -1;
    return 0;
}

/**
 * \brief	Closes an TCP/IP socket
 * \return				0 on success, otherwise error code
 */
int closeECBConnection()
{
	return close(sockfd);
}

/**
 * \brief	Sends a message over a TCP/IP socket
 * \param	msg		pointer to the message buffer
 * \param	length	number of bytes to send
 * \return				0 on success, otherwise error code
 */
int sendECBMessage(char* msg, unsigned int length)
{
	return write(sockfd,msg,length);
}

/**
 * \brief	Reads a message from a TCP/IP socket
 * \param	msg		pointer to the message buffer
 * \param	length	number of bytes to read
 * \return				0 on success, otherwise error code
 */
int readECBMessage(char* msg, unsigned int length)
{
	return read(sockfd,msg,length);
}

/**
 * \brief	Checks the TCP/IP socket
 * \return	0 on success, otherwise -1
 */
int checkSocket()
{
    if (sockfd < 0)
    	return -1;
    else
    	return 0;
}
#endif /* ECBLINUX */
