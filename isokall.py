#! python
import argparse
import os.path
from argparse import RawTextHelpFormatter


class MyArgumentParser(argparse.ArgumentParser):

    def __init__(self, *args, **kwargs):
        super(MyArgumentParser, self).__init__(*args, **kwargs)

        self.error_message = ''

    def error(self, message):
        self.error_message = message

    def parse_args(self, *args, **kwargs):
        # catch SystemExit exception to prevent closing the application
        result = None
        try:
            result = super().parse_args(*args, **kwargs)
        except SystemExit:
            pass
        return result


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
    Author: Amin Mahpour
    Support: https://gitlab.com/AminMahpour/Isokall
----------------------------------------------------------
    """.format(__version__)

    # Instantiate the parser
    parser = MyArgumentParser(prog=__prog_name__, description=usage, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='{version}'.format(version=__version__),
                        help='Get version number of this software.')

    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')

    parser_quant = subparsers.add_parser('quant', help='Add Kallisto quant data to BED12 files')
    # parser_quant.add_argument('bar', type=int, help='bar help')
    # Required positional argument
    parser_quant.add_argument('input_bed', type=str,
                              help='An input Bed12 formatted file from the ISO-seq pipeline')

    # Required positional argument
    parser_quant.add_argument('input_tsv', type=str,
                              help='The Kallisto tsv-formatted output file')

    parser_quant.add_argument('output_bed', type=str,
                              help='The name for output bed12-formatted file')

    parser_filter = subparsers.add_parser('filter', help='Filter by TSS distance.')
    parser_filter.add_argument('-dist_TSS', '-g', type=int, help='Distance from TSS (e.g. 100 nucleotide)')
    parser_filter.add_argument('-classification_file', '-c', type=str, help='Input classification file from Qanti2')
    parser_filter.add_argument('-fasta_file', '-f', type=str, help='Input FASTA file from Qanti2')
    parser_filter.add_argument('-input_genepred_file', '-gp', type=str, help='Input GenePred file from Qanti2')
    parser_filter.add_argument('-output_name', '-o', type=str, help='Prefix for output files')

    args = parser.parse_args()

    return args, parser


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
            if line[0] == "chrMT":
                line[0] = "chrM"

            if line[3] in l.keys():
                line[4] = l[line[3]]
                o = str.join("\t", [str(i) for i in line])

                out.write(f"{o}\n")

    out.close()


def main():
    args, parser = run_args()

    if args.subparser_name == "quant":
        if args.input_bed == None or args.input_tsv == None or args.output_bed:
            parser.print_help()
            print("Please provide all arguments.")
        else:
            run_conversion(args.input_bed, args.input_tsv, args.output_bed)

    elif args.subparser_name == "filter":
        if args.dist_TSS == None or args.classification_file == None or args.fasta_file == None or args.input_genepred_file == None or args.output_name:
            parser.print_help()
            print("Please provide all arguments.")
        else:
            filter_TSS(args.dist_TSS, args.classification_file, args.fasta_file, args.input_genepred_file,
                       args.output_name)
    else:
        parser.print_help()

def filter_TSS(dist, classification_txt, input_fasta, input_genepred, output_file):
    from Bio import SeqIO

    TSS_50 = []

    skip = False
    with open(classification_txt, "r") as data:

        for line in data:
            if skip and line.split("\t")[12] != "NA":
                # print(line.split("\t")[10])
                diff = abs(int(line.split("\t")[12]))
                # print(diff)
                if diff < dist:
                    TSS_50.append(line.split("\t")[0])
            skip = True

    number_of_records = len(TSS_50)

    print(f"{number_of_records} records were matched!")

    list_fasta_record = []
    with open(input_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = record.id.split("|")[0]

            if id in TSS_50:
                list_fasta_record.append(record)

    gp_record = []
    with open(input_genepred, "r") as handle:
        for line in handle:
            if line.split("\t")[0] in TSS_50:
                gp_record.append(line)

    print(gp_record)

    SeqIO.write(list_fasta_record, "{}.fasta".format(output_file), "fasta")

    with open("{}.gp".format(output_file), "w") as handle:
        handle.writelines([i for i in gp_record])


if __name__ == "__main__":
    main()
