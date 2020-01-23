from H_5j64_WaterTranslationRotation import WaterTranslationRotationMove
from blues.moves import MoveEngine
from blues.simulation import *
import json
from blues.settings import Settings

opt = Settings('H_5j64_mup_md.yaml').asDict()
structure = opt['Structure']
restart = parmed.amber.Rst7('/oasis/tscc/scratch/bergazin/water/H_5j64_equil_NPT.rst7')#change to dir where rst is

#structure.positions = restart.positions
#structure.box = restart.box

print(json.dumps(opt, sort_keys=True, indent=2, skipkeys=True, default=str))

#Select move type
ligand = WaterTranslationRotationMove(structure, water_name='WAT')

#Iniitialize object that selects movestep
ligand_mover = MoveEngine(ligand)

#Generate the openmm.Systems outside SimulationFactory to allow modifications
systems = SystemFactory(structure, ligand.atom_indices, opt['system'])

#Freeze atoms in the alchemical system
#systems.md = systems.freeze_atoms(structure, systems.md, **opt['freeze'])
#systems.alch = systems.freeze_atoms(structure, systems.alch, **opt['freeze'])

# Restrain atoms in the MD and alchemical system
#systems.md = systems.restrain_positions(structure, systems.md, **opt['restraints'])
#systems.alch = systems.restrain_positions(structure, systems.alch, **opt['restraints'])

#Generate the OpenMM Simulations
simulations = SimulationFactory(systems, ligand_mover, opt['simulation'], opt['md_reporters'], opt['ncmc_reporters'])

#Energy minimize system
simulations.md.minimizeEnergy(maxIterations=0)
#simulations.md.step(50000)
#state = simulations.md.context.getState(getPositions=True, getEnergy=True)
#print('Minimized energy = {}'.format(state.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole)))

# MD simulation
#simulations.md.step(opt['simulation']['nstepsMD'])

#BLUES Simulation
blues = BLUESSimulation(simulations)
blues.run(**opt['simulation'])

# MC Simulation
#mc = MonteCarloSimulation(simulations, opt['simulation'])
#mc.run()
