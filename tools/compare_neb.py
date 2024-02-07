import click
import re

@click.command()
@click.argument('file1', type=click.File('r'))
@click.argument('file2', type=click.File('r'))
def compare_extrema(file1, file2):
    """
    This script takes two files as arguments and finds the lines about the extrema,
    then compares them and prints out the differences.
    """

    def parse_extrema(file_content):
        # Regex to find lines with extrema data
        pattern = re.compile(r"extrema #(\d+) at image position ([\d\.]+) with energy ([\d\.\-]+) and curvature ([\d\.\-]+)")
        extrema = {}

        # Find all matches and store them in a dictionary
        for line in file_content:
            match = pattern.search(line)
            if match:
                extrema[int(match.group(1))] = {
                    'position': float(match.group(2)),
                    'energy': float(match.group(3)),
                    'curvature': float(match.group(4))
                }
        return extrema

    # Parse both files
    extrema1 = parse_extrema(file1)
    extrema2 = parse_extrema(file2)

    # Compare extrema and print the differences
    for extrema_number in extrema1:
        if extrema_number in extrema2:
            diff_position = abs(extrema1[extrema_number]['position'] - extrema2[extrema_number]['position'])
            diff_energy = abs(extrema1[extrema_number]['energy'] - extrema2[extrema_number]['energy'])
            diff_curvature = abs(extrema1[extrema_number]['curvature'] - extrema2[extrema_number]['curvature'])

            click.echo(f"Extrema #{extrema_number} differences:")
            click.echo(f"Position: {diff_position}")
            click.echo(f"Energy: {diff_energy}")
            click.echo(f"Curvature: {diff_curvature}")
        else:
            click.echo(f"Extrema #{extrema_number} found in file1 but not in file2")

if __name__ == '__main__':
    compare_extrema()
