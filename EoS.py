import numpy as np

from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.calculators.emt import EMT

a = 3.615  # approximate lattice constant
b = a / 2
cu = Atoms('Cu',
           cell=[(0, b, b), (b, 0, b), (b, b, 0)],
           pbc=1,
           calculator=EMT())  # use EMT potential
cell = cu.get_cell()
traj = Trajectory('Cu.traj', 'w')
for x in np.linspace(0.95, 1.05, 5):  # creates a loop that iterates over five evenly spaced values between 0.95 and 1.05 (inclusive)
    cu.set_cell(cell * x, scale_atoms=True)
    cu.get_potential_energy()
    traj.write(cu)

from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
configs = read('Cu.traj@0:5')  # read 5 configurations
# Extract volumes and energies:
volumes = [cu.get_volume() for cu in configs]
energies = [cu.get_potential_energy() for cu in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
eos.plot('Cu-eos.png')
