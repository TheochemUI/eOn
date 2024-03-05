import click
import numpy as np
import ase.io


def find_climbing_image(dat_file):
    # Read the .dat file to find the image with the maximum energy
    data = np.loadtxt(dat_file, skiprows=1)  # Skip the header
    energies = data[:, 2]  # Potential energy column
    ci_index = np.argmax(energies)
    return ci_index


def generate_dimer_inputs(
    neb_path_file, dat_file, displacement_magnitude=0.01, box_size=10.0
):
    ci_index = find_climbing_image(dat_file)
    neb_path = ase.io.read(neb_path_file, index=":")

    ci_atoms = neb_path[ci_index]
    # Set a large box around the atoms to prevent unwanted interactions
    ci_atoms.set_cell([box_size, box_size, box_size])
    ci_atoms.center()
    ase.io.write("pos.con", ci_atoms)

    # Calculate displacement vector based on difference with next or previous image
    if ci_index > 0:
        prev_image = neb_path[ci_index - 1]
        displacement_vector = ci_atoms.positions - prev_image.positions
    else:
        next_image = neb_path[ci_index + 1]
        displacement_vector = next_image.positions - ci_atoms.positions

    # Normalize and save the displacement vector as direction.dat
    norm_displacement_vector = (
        displacement_vector / np.linalg.norm(displacement_vector, axis=1)[:, None]
    )

    # Scale to fixed magnitude
    scaled_displacement_vector = norm_displacement_vector * displacement_magnitude

    # Apply displacement to CI atoms for displacement.con
    displaced_atoms = ci_atoms.copy()
    displaced_atoms.set_cell([box_size, box_size, box_size])
    displaced_atoms.center()
    displaced_atoms.positions += scaled_displacement_vector
    ase.io.write("displacement.con", displaced_atoms)

    # Save the scaled displacement vector
    np.savetxt("direction.dat", norm_displacement_vector, fmt="%f")


@click.command()
@click.argument("neb_path_file", type=click.Path(exists=True))
@click.argument("dat_file", type=click.Path(exists=True))
def main(neb_path_file, dat_file):
    """Generate Dimer input files from an EON NEB calculation."""
    generate_dimer_inputs(neb_path_file, dat_file)
    click.echo(
        f"Generated Dimer inputs using CI from {dat_file} and path from {neb_path_file}"
    )


if __name__ == "__main__":
    main()
