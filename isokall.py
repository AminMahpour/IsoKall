#! python
import argparse
import os.path
from argparse import RawTextHelpFormatter

__version__ = "0.0.1"
__prog_name__ = ""


def file_exist(filename):
    if os.path.isfile(filename):
        return True
    else:
        return False


def run_args():
    usage = """
---------------------------------------------------------
             .-'''-.                                      
            '   _    \                         .---..---. 
.--.      /   /` '.   \     .                  |   ||   | 
|__|     .   |     \  '   .'|                  |   ||   | 
.--.     |   '      |  '.'  |                  |   ||   | 
|  |     \    \     / /<    |            __    |   ||   | 
|  |     _`.   ` ..' /  |   | ____    .:--.'.  |   ||   | 
|  |   .' |  '-...-'`   |   | \ .'   / |   \ | |   ||   | 
|  |  .   | /           |   |/  .    `" __ | | |   ||   | 
|__|.'.'| |//           |    /\  \    .'.''| | |   ||   | 
  .'.'.-'  /            |   |  \  \  / /   | |_'---''---' 
  .'   \_.'             '    \  \  \ \ \._,\ '/           
                       '------'  '---'`--'  `"            
    Ver: {0}
    Support: https://gitlab.com/mahpour/Isokall
----------------------------------------------------------
    """.format(__version__)

    # Instantiate the parser
    parser = argparse.ArgumentParser(prog=__prog_name__, description=usage, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='{version}'.format(version=__version__),
                        help='Get version number of this software.')

    # Required positional argument
    parser.add_argument('input_bed', type=str,
                        help='An input Bed12 formatted file from the ISO-seq pipeline')

    # Required positional argument
    parser.add_argument('input_tsv', type=str,
                        help='The Kallisto tsv-formatted output file')

    parser.add_argument('output_bed', type=str,
                        help='The name for output bed12-formatted file')

    args = parser.parse_args()

    return args


def run_conversion(input_bed, input_abund, output_bed):
    l = dict()
    b = True
    with open(input_abund, mode="r") as x:
        for line in x.readlines():
            if b is False:
                line = line.strip("\n").split("\t")
                l.update({line[0]: float(line[4])})
            b = False
    out = open(output_bed, mode="w")
    chr_list = ["chr{0}".format(x) for x in list(range(0, 23))] + ["chrX", "chrY"]

    with open(input_bed, mode="r") as y:
        for line in y.readlines():
            line = line.strip("\n").split("\t")
            if line[0] not in chr_list: continue
            if line[0] == "chrMT": line[0] = "chrM"

            if line[3] in l.keys():
                line[4] = l[line[3]]
                o = str.join("\t", [str(i) for i in line])

                out.write(f"{o}\n")

    out.close()


def main():
    args = run_args()
    run_conversion(args.input_bed, args.input_tsv, args.output_bed)


if __name__ == "__main__":
    main()
