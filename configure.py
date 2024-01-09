import logging, os, _collections, click, gzip, shutil

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)


root_dir = os.path.dirname(os.path.abspath(__file__))
executables = dict(
    pigz = 'gzip',
    gzip = 'gzip',
    getorf = shutil.which('getorf'),
    transeq = shutil.which('transeq'),
    hmmsearch = shutil.which('hmmsearch'), #'/usr/bin/hmmsearch',
    bindash = shutil.which('bindash'), #'/titan/softwares/bin/bindash', 
    minimap2 = shutil.which('minimap2'), #'/titan/softwares/bin/minimap2',
    samtools = shutil.which('samtools'), #'/titan/softwares/bin/samtools',
    iqtree = shutil.which('iqtree'), #'/titan/softwares/bin/iqtree'
    EnFlt    = os.path.join(root_dir, '_EnFlt.py'),
)

complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', '-':'-',
              'a':'t', 't':'a', 'g':'c', 'c':'g', '-':'-'}
def rc(seq) :
    return ''.join([ complement.get(s, 'N') for s in seq[::-1] ])


def readFasta(fname) :
    seqs = _collections.OrderedDict()
    if fname.lower().endswith('.gz') :
        with gzip.open(fname, 'rt') as fin :
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    seqs[name] = []
                else :
                    seqs[name].extend(line.strip().split())
    else :
        with open(fname, 'rt') as fin :
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    seqs[name] = []
                else :
                    seqs[name].extend(line.strip().split())

    for n, s in seqs.items() :
        seqs[n] = ''.join(s)
    return seqs    



if __name__ == '__main__' :
    main()
