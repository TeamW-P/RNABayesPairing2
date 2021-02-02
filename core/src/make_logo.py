#!/var/www/html/cgi-bin/carnaval/Src/VEnv/bin/python
import sys
from weblogo import *


def help():
    return """Required args:
            -i <input file> list of sequences, one per line
            -o <output file name> name of output, will add .eps
        Optionals:
            -t <string> title
            -h          this help displayed
           """


def make_logo(in_file, out_file):
    fin = open(in_file)
    seqs = read_seq_data(fin)
    data = LogoData.from_seqs(seqs)
    options = LogoOptions()

    options = LogoOptions()
    options.color_scheme = colorscheme.nucleotide
    options.unit_name = 'probability'
    options.fineprint = ''
    options.creator_text = ''

    format = LogoFormat(data, options)
    out = png_formatter(data, format)
    with open('%s.png' % out_file, 'wb')  as f:
        f.write(out)


if __name__ == '__main__':
    in_file = None
    out_file = None
    title = None
    opts = sys.argv

    for i, o in enumerate(opts):
        if o == '-i':
            in_file = opts[i + 1]
        elif o == '-o':
            out_file = opts[i + 1]
        elif o in ('-h', '-help'):
            print(help())

    if not in_file or not out_file:
        # print "need args -i and -o"
        sys.exit(1)

    make_logo(in_file, out_file)

