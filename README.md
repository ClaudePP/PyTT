#  -----------------     PyTT     -----------------    # 

 Finite difference code to predict the thermal evolution of different classes
 of thin target detectors. For applications in particle accelerators.
 
 Developed by: Araceli Navarro & Mariusz Sapinski
 
 Contact: n.f.araceli@gmail.com or mariusz.sapinski@psi.ch
 
![plot](./HelpFolder/PyTTScreanshot.png)
![plot](./HelpFolder/PyTTresults.png)


see also: https://sapinski.web.cern.ch/sapinski/soft/pyTT/index.html



# -------- Quick start ----------------------------  #

Running from command line:

> python3 MAIN.py Simulations/PSI_RRL123MeV.txt

Running GUI (currently GUI is not working):

> python3 MAIN.py


from IPython console (eg. spyder 5)

cd the_top_directory

> run MAIN.py 

Running tests from command line (temporary solution):

> python3 Test.py


# -------- PyTT virtual environment ----------------  #
2025.01.29
in: /Users/sapinski/Phys/ThinTarget/PyTT
> python3 -m venv pytt.venv - create
> source pytt.venv/bin/activate
> deactivate (at the end)

installing additional modules via brew withion activated environment:
brew install python-matplotlib
- no, use pip3 install matplotlib


---------------------------------------

detector types:

WIRESCAN - wire scanner, single wire moving throught the beam
SEM - SEM-grid or harf, a set of fixed, parallel wires


