# Installation instructions
## Installing on Windows vs Linux
In order to choose whether to compile for Windows or Linux, choose the appropriate setup code in `setup.py`
## Installation Commands
To install 
    `python setup.py install`

Compile locally 
    `python setup.py build`
To build the pyECB extension module.
Then navigate to the build folder where the module was compiled and run python from there

# Using pyECB
## Network Settings
### Linux
The following network settings have been verified on Ubuntu 18.04:
* IP of Python-PC : `192.168.237.1`
* Subnet : `192.168.237.0/24`
* Gateway : `<EMPTY>`
**HINT** : In the Ubunut network manager, the gateway field should be left blank, otherwise there is an issue building up the connection.
### Windows
The following network settings have been verified on Windows 10:
* IP of Python-PC : `192.168.237.1`
* Subnet : `192.168.237.0/24`
* Gateway : `192.168.237.255`
In the ethernet settings of Windows, a gateway had to be specified even though none is used.

## Sample Code
The bindings are pretty much one to one with the C API. See the following example.

     from ECB import *
     
     initECBapi("192.168.237.47", "7070")
	 enableECBCurrents()
     setDesCurrents([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000],b'1')
     
     # Do experiment
     
     disableECBCurrents()
     
     exitECBapi()
