from setuptools import setup

setup(
   name='python_scripts',
   version='1.0.1',
   description='Personal Python Scripts for Research',
   author='Dustin Tracy',
   author_email='dtracy.uf@gmail.com',
   packages=['python_scripts'],  #same as name
   install_requires=['numpy', 'matplotlib'], #external packages as dependencies
   scripts=[
            'scripts/dipole.py',
            'scripts/solvent_dipole_rotation.py',
            'scripts/dihedral.py',
            'scripts/bla.py',
            'scripts/omega1.py',
            'scripts/simpleplot.py',
            'scripts/myaverage.py',
            'scripts/amber_compare_total_energy.py',
           ]
)
