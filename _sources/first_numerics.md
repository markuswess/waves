(first_numerics)=
# Basic numerical methods for time-domain waves

In this section we work our way through the steps leading to a first finite element method for the wave-type equations from {numref}`modelling` using the high-order finite element library [Netgen/NGSolve](https://ngsolve.org).

Note that this is not the time and place to present a full-sized introduction the the finite element method in all its glory. What follows is a shortened and adapted version of the first chapters of Joachim Schöberl's [iFEM](https://jschoeberl.github.io/iFEM/) introduction, which is also a highly recommended read for the interested student.

## Running NGSolve

By now there exist various ways to run NGSolve

### Local installation
Install a recent Python. Then it should be easy to install NGSolve using

    pip install jupyter numpy scipy matplotlib
    pip install --pre ngsolve
    pip install webgui_jupyter_widgets


To check the installation of NGSolve run in a terminal:

    python3 -c "import ngsolve; print(ngsolve.__version__)"

Then, open jupyter-notebook (or jupyter-lab or VS Code), create a new notebook, create and execute a cell with

    from ngsolve import *
    from ngsolve.webgui import Draw
    Draw (unit_cube.shape);


Known issues are
- Use pip3 instead of pip if there is no pip
- If you get an error like `externally-managed-environment`, then either use
virtual environments, or add the flag `--break-system-packages` to the pip command, see [explanation](https://veronneau.org/python-311-pip-and-breaking-system-packages.html)

- If you have conflicts with other packages, you may install NGSolve in a [virtual environment](https://docs.python.org/3/library/venv.html#creating-virtual-environments). For example 

      python3 -m venv /Users/XXXXX/numpde
      source /Users/XXXXX/numpde/bin/activate

- If NGSolve compuatations are working, but you don't get the rendering: For jupyter notebook version < 7.0.0 you have to run additionally

      jupyter nbextension install --user --py webgui_jupyter_widgets
      jupyter nbextension enable --user --py webgui_jupyter_widgets


### Run NGSolve remote or via jupyter-lite

If local installation does not work, there are alternatives:

- login to a jupyter server from your browser:

  [jupyterhub.cerbsim.com](https://jupyterhub.cerbsim.com) <br>
  user: **ngshub_xx** <br>
  pwd:  **solve!xx** <br>
  with xx number from 01 to 31

  

- run NGSolve online within jupyter-lite:

  [https://markuswess.github.io/waves-lite](https://markuswess.github.io/waves_lite)

  The first time it might take a few minutes to start, and then again to import ngsolve.
