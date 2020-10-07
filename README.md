# pyECB_gmaxwell
A Python wrapper based on the pyECB extention (by Denis von Arx, @dvarx)

## Instructions on how to use this package:
After cloning this repository follow these steps (this is for Windows, but it shouldn't differ too much on Linux/Unix).
* **Create a virtual environment:**
1. open your current working directory in the command line/terminal (I would recommend using _this_ folder, i.e. `.\pyECB_gmaxwell\`)
2. type `virtualenv --version` and ensure that this package is installed
  (2a. if not, install it by typing `pip install virtualenv`)
3. use the command `python -m venv [env]`, replace env with any name you like
4. activate the virtual environment by typing `[env]\Scripts\activate`
5. make sure you have the latest vesion of pip: `python -m pip install --upgrade pip`
6. install required packets (non-exhaustive): `numpy, pyserial, pythonnet, matplotlib, ...`

* **install pyECB extension:**
1. follow the instructions in `.\pyECB_gmaxwell\pyECB_original\README.md` (within your virtual environment!!)

## Notes:
If you are using an IDE, make sure to choose the `python.exe` file in your virtual environment as an interpreter!
