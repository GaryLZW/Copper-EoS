def workflow_slab():
   #You can define this EOS function outside the workflow function.
   #You did not set an parameters for EOS(). What if the user want to change the element and lattice parameter then?
   #There is no return value. The optimised lattice parameter is only printed.
   def EOS():      
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
       for x in np.linspace(0.95, 1.05, 5):
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
       vol_cubic = (v0 * 4)
       a = vol_cubic ** (1 / 3)
       print(a)

   #bulk
   #You can define this function outside the workflow function too. This would make your code looks tidy and easier to maintain
   #Why do you want this bulk parameter? It seems redundant to me.
   #Parameter fcc and cubic not used in the function. They are supposed 
   #to be used when you call the bulk function, right?

   #Not really sure why you delete the slab part of this function. We should build
   #the bulk structure inside this function and use it to construct slabs.
   def generate_slabs(bulk, Cu,  fcc, a, cubic=True):
       from ase.build import bulk
       from ase.visualize import view

       if bulk:
        a1 = bulk('Cu', 'fcc', a=a)
        view(a1)
 
   EOS() #This line only prints out the lattice parameter you get. Add a return statement to the function. 
         #Then a = EOS() stores the value of your optimised lattice parameter
   generate_slabs(bulk=True, Cu='Cu', fcc='fcc', a=a, cubic=True)

workflow_slab()
