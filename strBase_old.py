import os, gzip, click, glob, re, igraph
import subprocess, tempfile, multiprocessing
import numpy as np, pandas as pd


try :
    from configure import executables, logging, rc, readFasta
except :
    from .configure import executables, logging, rc, readFasta



def detranseq(d) :
    d[0], f = d[0].rsplit('_', 1)
    f = int(f)
    if f <= 3 :
        d[3], d[4] = (d[3]-1)*3+f, d[4]*3-1+f
        d[7] = 1
    else :
        qlen = d[9] - (1 if f>4 else 0)
        d[3], d[4] = qlen - d[3] + 1, qlen - d[4] + 1
        d[3], d[4] = -(d[3]*3-1+f-3), -((d[4]-1)*3+f-3)
        if d[4] > -1 :
            d[4] -= 3
        d[7] = -1
    return d


def parseHits(data, c) :
    data.sort(key=lambda x:[-x[8], x[0], x[3], x[4]])
    outputs = []
    for d in data :
        ingroup = False
        for oo in outputs :
            if d[0] == oo[0][0] and (d[3] > 0) == (oo[0][3] > 0) and \
                    min(abs(oo[0][3] - d[4]), abs(oo[0][4] - d[3])) < 10000 :
                ingroup = True
                minus = 0
                for o in oo :
                    if (d[3] - o[3]) * (d[5] - o[5]) < 0 or (d[4] - o[4]) * (d[6] - o[6]) < 0 :
                        ingroup = False
                        break
                    o1 = min(o[4], d[4]) - max(o[3], d[3]) + 1
                    o2 = min(o[6], d[6]) - max(o[5], d[5]) + 1
                    if o1 >= 0.9 * min(o[4] - o[3]+1, d[4] - d[3]+1) or o2 >= 0.9 * min(o[6] - o[5]+1, d[6] - d[5]+1) :
                        ingroup = False
                        break
                    if o2 > 0 :
                        minus += o2/(d[6] - d[5]+1)*d[8]
                if ingroup :
                    d[8] = d[8] - minus if d[8] > minus else 0.
                    oo.append(d)
                    break
        if not ingroup :
            outputs.append([d])
    o2 = []
    for dat in outputs :
        score = 0
        cov = np.zeros(dat[0][10])
        for d in dat :
            cov[d[5]-1:d[6]] = 1
            score += d[8]
        n_cov = np.sum(cov)
        if score >= c[0] and n_cov >= c[2]*0.4 and n_cov >= dat[0][10]*0.5 :
            o2.append(sorted(dat, key=lambda d:[d[3], d[4]]))
    return o2

def get_uscgs(query, dirname, acc) :
    subprocess.Popen('{transeq} -frame 6 -table 4 -sequence {0} -nomethionine -outseq {1}'.format(
        query, os.path.join(dirname,'{0}.aa'.format(acc)), **executables).split(), stderr=subprocess.PIPE).communicate()

    dbname = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'uscg_profile')
    dataset = []
    cutoffs, toMerge = {}, {}

    with open(os.path.join(dbname, 'merged_ortho')) as fin :
        for line in fin :
            p = line.strip().split()
            toMerge[p[1]] = p[0]
    with open(os.path.join(dbname, 'scores_cutoff')) as fin :
        for line in fin :
            p = line.strip().split()
            cutoffs[p[0]] = [float(p[1]), 0, 0]
    with open(os.path.join(dbname, 'lengths_cutoff')) as fin :
        for line in fin :
            p = line.strip().split()
            cutoffs[p[0]][1:] = [float(p[2]), float(p[3])]

    hmms = os.path.join(dbname, 'hmms', '*.hmm')
    for hmm in sorted(glob.glob(hmms)) :
        data = []
        key = os.path.basename(hmm)[:-4]
        subprocess.Popen('{hmmsearch} --notextw --noali -T {3} --cpu {4} --domT {3} --domtblout {2} {0} {1} '.format(
            hmm, os.path.join(dirname,'{0}.aa'.format(acc)), os.path.join(dirname,'{0}.hmm'.format(acc)),
            cutoffs[key][0]*0.4, 2, **executables
        ).split(), stdout=subprocess.PIPE).communicate()
        with open(os.path.join(dirname, '{0}.hmm'.format(acc))) as fin :
            for line in fin :
                if line.startswith('#') :
                    continue
                p = line.strip().split()
                d = [p[0], toMerge.get(p[3], p[3]), float(p[21]), int(p[17]), int(p[18]), \
                     int(p[15]), int(p[16]), float(p[12]), float(p[13]), int(p[2]), int(p[5])]
                data.append(detranseq(d))
        if len(data) :
            dataset.extend(parseHits(data, cutoffs[key]))

    for d in dataset :
        if d[0][3] < 0 :
            for dd in d :
                dd[3:5] = -dd[4], -dd[3]
            d.sort(key=lambda dd:dd[3])
    dataset.sort(key=lambda d:[d[0][0], d[0][3]] )


    for i, d1 in enumerate(dataset) :
        if d1[0][0] == '' :
            continue
        for d2 in dataset[i+1:] :
            if d2[0][0] == '' or d1[0][0] != d2[0][0] :
                break
            overlap = min(d1[-1][4], d2[-1][4]) - max(d1[0][3], d2[0][3]) + 1
            if overlap >= 0.6 * (d1[-1][4] - d1[0][3]+1) or overlap >= 0.6 * (d2[-1][4] - d2[0][3]+1) :
                s1 = sum([ d[8] for d in d1 ])
                s2 = sum([ d[8] for d in d2 ])
                if s1 >= s2 :
                    d2[0][0] = ''
                else :
                    d1[0][0] = ''
                    break
    dataset = sorted([ d for d in dataset if d[0][0] != '' ], key=lambda x:[x[0][1], x[0][0], x[0][3]])

    sequences = readFasta(query)
    ids = {}
    outfile = os.path.join(dirname, '{0}.USCGs.ffn.gz'.format(acc))
    with gzip.open(outfile, 'wt') as ffn_out :
        for data in dataset :
            s = sequences[data[0][0]][data[0][3]-1:data[-1][4]]
            if data[0][7] < 0 :
                s = rc(s)
            ids[data[0][1]] = ids.get(data[0][1], 0) + 1
            n = '{0}__{1}__{2}'.format(data[0][1], acc, ids[data[0][1]])
            ffn_out.write('>{0} {2} {3} {4} {5}\n{1}\n'.format(
                n, s, data[0][0], data[0][3], data[-1][4], data[0][7]))
    return outfile


def prepare_uscg(data) :
    in_fna, acc, tmpdir = data

    if in_fna.lower().endswith('.gz') :
        subprocess.Popen('{gzip} -cd {0} > {1}'.format(
            in_fna,
            os.path.join(tmpdir, '{0}.fna'.format(acc)),
            **executables), shell=True).communicate()
        in_fna = os.path.join(tmpdir, '{0}.fna'.format(acc))

    uscg_ffn = get_uscgs(in_fna, tmpdir, acc)
    return acc, uscg_ffn


def get_representative(references, threads, tmpdir, min_dist) :
    ref_ids = {}
    with open(os.path.join(tmpdir, 'reference.list'), 'wt') as fout :
        for rid, (fname, acc, species) in enumerate(references) :
            fout.write('{0}\n'.format(fname))
            ref_ids[fname] = rid
    
    subprocess.Popen('{bindash} sketch --nthreads={0} --listfname=reference.list --sketchsize64=128 --outfname=reference'.format(threads, **executables).split(), cwd=tmpdir).communicate()
    subprocess.Popen('{bindash} dist --mthres={1} --nthreads={0} --outfname=reference.dist reference reference'.format(threads, min_dist, **executables).split(), cwd=tmpdir).communicate()
    
    data = pd.read_csv(os.path.join(tmpdir, 'reference.dist'), sep='\t', dtype=str)
    data = data.values[:, :3]
    data = data[data.T[0] != data.T[1]]
    edges = np.vectorize(ref_ids.get)(data[:, :2])
    g = igraph.Graph(len(ref_ids))
    g.add_edges(edges)
    m = np.array(g.community_leiden(resolution=1e-6).membership)
    
    reprs = []
    for n in np.unique(m) :
        members = references[m == n]
        s, x = np.unique(members.T[2], return_counts=True)
        reprs.append([members[0, 0], members[0, 1], s[np.argmax(x)]])
    return np.array(reprs)        





@click.command()
@click.option('-i', '--fna_list', help='comma-delimited list of the included fasta files in a format of "file_name,accession,species".', required=True)
@click.option('-d', '--min_dist', help='minimum distance between representatives. [default: 0.01]', default=0.01, type=float)
@click.option('-c', '--min_cov', help='minimum presences of genes and genomes. [default: 0.75]', default=0.75, type=float)
@click.option('-p', '--prefix', help='prefix for the output [required].', required=True)
@click.option('-t', '--threads', help='number of threads to use [Default: 20]', type=int, default=20)
def main(fna_list, min_dist, min_cov, prefix, threads) :
    strBase(fna_list, min_dist, min_cov, prefix, threads)

def strBase(fna_list, min_dist, min_cov, prefix, threads) :
    pool = multiprocessing.Pool(int(threads/2))
    meta_tab = pd.read_csv(fna_list, header=None).values
    
    references = []
    for fname, acc, species in meta_tab :
        if os.path.isfile(fname) :
            references.append([fname, acc, species])
    references = np.array(references)

    uscgs = {}
    files = {}
    with tempfile.TemporaryDirectory(dir='.') as tmpdir :
        representatives = get_representative(references, threads, tmpdir, min_dist)
        for id, (acc, uscg_ffn) in enumerate(pool.imap_unordered(prepare_uscg, [[r[0], r[1], tmpdir] for r in representatives])) :
            if acc :
                files[acc] = uscg_ffn
                seq = readFasta(uscg_ffn)
                if acc not in uscgs  :
                    uscgs[acc] = {}
                for n in seq :
                    ortho = re.findall('^(.+)__.+__\d+$', n)[0]
                    uscgs[acc][ortho] = uscgs[acc].get(ortho, 0) + 1
            if id % 1000 == 0 :
                logging.info('Processing {0} genomes.'.format(id))
        genes = {}
        for acc, orthos in uscgs.items() :
            for ortho, cnt in orthos.items() :
                if cnt == 1 :
                    genes[ortho] = genes.get(ortho, 0) + 1
        genes = { o:c for o, c in genes.items() if c >= min_cov * len(uscgs) }
        uscgs = { acc: {o:1 for o, c in orthos.items() if o in genes and c == 1 } for acc, orthos in uscgs.items() }
        uscgs = { acc: orthos for acc, orthos in uscgs.items() if len(orthos) >= min_cov * len(genes) }
        
        with gzip.open(prefix+'.ffn.gz', 'wt') as fout :
            for acc, orthos in sorted(uscgs.items()) :
                fname = files[acc]
                seqs = readFasta(fname)
                for n, s in sorted(seqs.items()) :
                    ortho = re.findall('^(.+)__.+__\d+$', n)[0]
                    if ortho in genes :
                        fout.write('>{0}\n{1}\n'.format(n, s))
        with open(prefix + '.list', 'wt') as fout :
            for fname, acc, species in references :
                if acc in uscgs :
                    fout.write('{0},{1}\n'.format(acc, species))

if __name__ == '__main__' :
    strBase()
