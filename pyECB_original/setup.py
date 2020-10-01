#from distutils.core import setup, Extension

from distutils.core import setup, Extension

# chose this for compilation on Windows
extension_mod = Extension("ECB", 
                          sources=["pyECB.c",
                                   "ecb_api.c",
                                   "ecb_func.c",
                                   "ecb_socket_win.c",
                                   ],
                          include_dirs=["."],
                          libraries=["ws2_32"]
                          )
                          
'''    
# chose this for compilation on Linux              
extension_mod = Extension("ECB", 
                          sources=["pyECB.c", 
                              "ecb_api.c",
                              "ecb_func.c",
                              "ecb_socket.c",
                              ], 
                          include_dirs=["."],
			define_macros = [('ECBLINUX', '1')],
                          )
''' 

setup(name="ECB",
      ext_modules=[extension_mod],
      version="0.0.1")
