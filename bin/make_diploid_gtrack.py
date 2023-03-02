import sys
import re
import argparse
import pathlib

def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "input",
        type=pathlib.Path,
        help="Path to the input Gtrack file.",
    )

    return cli
if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    gtrack_fname = str(args[1])

    with open(gtrack_fname) as gf:
        gtrack_file = gf.readlines()

    for line in gtrack_file:
        line=line.rstrip()
        if line.startswith("#"):
            print(line)
        else:
            line1=re.sub(r'(chr\w+)',r'\1_A',line)
            line2=re.sub(r'(chr\w+)',r'\1_B',line)
            print(line1)
            print(line2)
