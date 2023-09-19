from orgdisord.utils import molecule_collide
from ase.build import molecule
from ase.visualize import view
from ase import Atom

atoms = molecule("CH3CH2OCH3", cell=[10, 10, 10], pbc=True)
initial_fractional_positions = atoms.get_scaled_positions()
cell = atoms.cell
false_test_positions = [
    initial_fractional_positions[0] + 0.21,
    [0.25, 0.25, 0.25],
    [0.5, 0.5, 0.5],
]
true_test_positions = [
    initial_fractional_positions[0] + 0.01,
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
]

assert all(
    molecule_collide(Atom("H", position=cell.T.dot(pos)), atoms, tolerance=0.2) == False
    for pos in false_test_positions
)
assert all(
    molecule_collide(Atom("H", position=cell.T.dot(pos)), atoms, tolerance=0.2) == True
    for pos in true_test_positions
)

# for pos in false_test_positions:
#     print(molecule_collide(pos, atoms, tolerance=0.2))
# for i, pos in enumerate(test_positions):
#     atoms.append(Atom('X', cell.T.dot(pos), magmom=5))
# atoms.set_tags(range(len(atoms)))
# atoms.translate([5,5,5])
# atoms.wrap()
# view(atoms)
