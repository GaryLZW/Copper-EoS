def EOS(name='Cu', crystalstructure='fcc', a=3.615):
    import numpy as np

    from ase import Atoms
    from ase.build import bulk
    from ase.io.trajectory import Trajectory
    from ase.calculators.emt import EMT

    a = a  # approximate lattice constant
    b = a / 2
    elem = bulk(name, crystalstructure, a)
    elem.calc = EMT()
    cell = elem.get_cell()
    traj = Trajectory('elem.traj', 'w')
    for x in np.linspace(0.95, 1.05, 5):
        elem.set_cell(cell * x, scale_atoms=True)
        elem.get_potential_energy()
        traj.write(elem)

    from ase.io import read
    from ase.units import kJ
    from ase.eos import EquationOfState
    configs = read('elem.traj@0:5')  # read 5 configurations
    # Extract volumes and energies:
    volumes = [elem.get_volume() for elem in configs]
    energies = [elem.get_potential_energy() for elem in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    print(B / kJ * 1.0e24, 'GPa')
    eos.plot('elem-eos.png')
    vol_cubic = (v0 * 4)
    a = vol_cubic ** (1 / 3)
    print(a)
    return a


def generate_slabs(name, lattice, a, layers=3, vacuum=10, indicesList=[(1, 0, 0), (0, 1, 0), (0, 0, 1),
                                                                       (1, 1, 1), (2, 1, 1), (1, 2, 1),
                                                                       (1, 1, 2)]):
    from ase.build import bulk
    a1 = bulk(name=name, crystalstructure=lattice, a=a, cubic=True)

    from ase.build import surface
    slabList = []
    for indices in indicesList:
        # create slabs
        slab = surface(lattice=a1, indices=indices, layers=layers, vacuum=vacuum)
        slabList.append((indices, slab))

    # All slabs created
    return slabList


def workflow_slab(name, crystalstructure, a, layers, vacuum, indices):
    a_eos = EOS(name, crystalstructure, a)

    slabs = generate_slabs(name, crystalstructure, a_eos, layers, vacuum, indices)

    from ase.io.trajectory import Trajectory
    traj = Trajectory('slabs.traj', 'w')
    for slab in slabs:
        traj.write(slab[1])

    return slabs


slabs = workflow_slab(name='Cu', crystalstructure='fcc', a=3.615, layers=3, vacuum=10.0,
                      indices=[(1, 0, 0), (1, 1, 0), (1, 1, 1)])

from ase.visualize import view
from ase.io import read
slabs0 = read('slabs.traj@:')
view(slabs0)
print(slabs[-1][0])
