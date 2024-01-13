import os, click, gzip, json, numpy as np, ete3, pandas as pd
import subprocess, tempfile, re
from multiprocessing import Pool
import ete3_extensions

try :
    from configure import executables, logging, rc, readFasta
except :
    from .configure import executables, logging, rc, readFasta




cc = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'-',
      'a':'t', 't':'a', 'g':'c', 'c':'g', 'n':'n'}
def rc(s) :
    return [cc.get(x, 'N') for x in s[::-1]]


def each_minimap(data) :
    tmpdir, target, min_iden, min_presence, idx = data
    r = os.path.basename(target).rsplit('.', 1)[0]
    if target.lower().endswith('fq.gz') or target.lower().endswith('fastq.gz') :
        with gzip.open(target, 'rt') as fin, open(os.path.join(tmpdir, '{0}.fas'.format(idx)), 'wt') as fout :
            for lid, line in enumerate(fin) :
                if lid % 4 == 0 :
                    fout.write('>' + line[1:])
                elif lid % 4 == 1 :
                    fout.write(line)
        target = '{0}.fas'.format(idx)
    elif target.lower().endswith('.gz') :
        with gzip.open(target, 'rt') as fin, open(os.path.join(tmpdir, '{0}.fas'.format(idx)), 'wt') as fout :
            for line in fin :
                fout.write(line)
        target = '{0}.fas'.format(idx)
    t_seq = readFasta(os.path.join(tmpdir, target))

    map_cmd = '{0} -k13 -w5 -c -t1 --frag=yes --rmq -A1 -B4 -O8,16 -E2,1 -r20k,40k -g10k -P -N5000 -f1000,5000 -n2 -m20 -s30 -z200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(
        executables['minimap2'], 'query', target
    )

    p = subprocess.Popen(map_cmd.split(), cwd=tmpdir, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    res = []
    for line in p.stdout :
        if line.startswith('[') :
            continue
        p = line.strip().split('\t')
        p[9:11] = int(p[9]), 100. * float(p[9])/float(p[10])
        if p[10] <= min_iden :
            continue

        p[1:4] = [int(p[1]), int(p[2])+1, int(p[3])]
        p[6:9] = [int(p[6]), int(p[7])+1, int(p[8])]

        if p[4] == '+' :
            qi, ri, cigar, d = p[2], p[7], p[-1][5:], 1
        else :
            qi, ri, cigar, d = p[3], p[7], p[-1][5:], -1
        p[11] = cigar
        p[12] = int(p[14][5:])
        res.append(p[:13])
    res.sort(key=lambda x:(x[5], -x[12]))
    
    res_seqs = {}
    for p in res :
        if p[0].find('__') >= 0 :
            gid = re.split('__', p[0])[1]
        else :
            gid = r
        key = (gid, p[5])
        if key in res_seqs :
            s = res_seqs[key]
            if np.sum(np.array(s[p[7]-1:p[8]]) != '-') >= 0.1 * (p[8] - p[7] + 1) :
                continue
        else :
            s = ['-'] * p[6]
        if p[4] == "+" :
            i0, i1, d = p[2]-1, p[7]-1, 1
        else :
            i0, i1, d = p[3], p[7]-1, -1
        for n, t in re.findall('(\d+)([MDI])', p[11]) :
            n = int(n)
            if t == 'M' :
                if p[4] == '+' :
                    ss = t_seq[p[0]][i0:(i0+n)]
                else :
                    ss = rc(t_seq[p[0]][(i0-n):i0])
                s[i1:(i1+n)] = ss
                i0, i1 = i0 + n*d, i1 + n
            elif t == 'D' :
                i1 += n
            elif t == 'I' :
                i0 += n*d
        res_seqs[key] = s
    res = {}
    for (strain, gene), s in res_seqs.items() :
        s = re.sub('[^ACGT]', '-', ''.join(s))
        if len(s.replace('-', '')) >= len(s) * min_presence :
            if strain not in res :
                res[strain] = {}
            res[strain][gene] = s
    return idx, res



def minimap_align(tmpdir, ref_acc, ref_seqs, db_seqs, items, genomes, min_identity, min_presence, pool) :
    logging.info('Building phylogeny.')
    with open(os.path.join(tmpdir, 'query'), 'wt') as fout :
        for n, s in ref_seqs.items() :
            fout.write('>{0}\n{1}\n'.format(n, s))

    with open(os.path.join(tmpdir, 'from_db'), 'wt') as fout :
        for n, s in db_seqs.items() :
            fout.write('>{0}\n{1}\n'.format(n, s))

    with open(os.path.join(tmpdir, 'target'), 'wt') as fout:
        for key, seqx in items.items() :
            for n, s in seqx.items() :
                ns = re.split('__', n)
                try :
                    fout.write('>{0}__{1}|{2}__{3}\n{4}\n'.format(ns[0], key[0], key[2], ns[-1], s))
                except :
                    pass

    genes = sorted(ref_seqs.keys())
    qryseqs = {}
    for idx, s in pool.imap_unordered(each_minimap, [ [tmpdir, genome, min_identity, min_presence if idx==0 else 0.5, idx]
                                      for idx, genome in enumerate(['target', 'from_db']+genomes) ]) :
        for tag, seq in s.items() :
            if idx == 0 and len(seq) >= min_presence * len(genes) :
                qryseqs[tag] = seq
            elif idx > 0 and len(seq) >= 0.6 * len(genes) :
                qryseqs[tag] = seq
        if len(qryseqs) % 100 == 0 :
            logging.info('Extracted {0} USCGs from both the db and the samples.'.format(len(qryseqs)))
    logging.info('Extracted {0} USCGs from both the db and the samples.'.format(len(qryseqs)))
    concatenated_seqs = {ref_acc: ''.join([ref_seqs[g].replace('n', '-').replace('N', '-') for g in genes])}

    for qry, qryseq in qryseqs.items() :
        concatenated_seqs[qry] = ''.join([ qryseq.get(g, '-'*len(ref_seqs[g])) for g in genes ])
    logging.info('Identified {0} samples with good sequences.'.format(len(concatenated_seqs)))
    with open(os.path.join(tmpdir, 'USCG_align.fas'), 'w') as fout :
        for n, s in concatenated_seqs.items() :
            fout.write('>{0}\n{1}\n'.format(n, s))

    subprocess.Popen('{iqtree} -nt 8 -fast -redo -s USCG_align.fas -m GTR -pre USCG_align --runs 5 --polytomy'.format(
        **executables).split(), stdout=subprocess.PIPE, cwd=tmpdir).communicate()
    tre = ete3.Tree(os.path.join(tmpdir, 'USCG_align.treefile'), format=0)
    tre.set_outgroup(tre.get_midpoint_outgroup())
    for i, n in enumerate(tre.traverse('postorder')) :
        if not n.name :
            n.name = 'N_{0}'.format(i)
    for leaf in tre.get_leaves() :
        leaf.annotations = {'gene_presence' : len(qryseqs[leaf.name]), 'nucl_presence' : len(concatenated_seqs[leaf.name].replace('-', ''))}
    return tre


def readJson(param) :
    uscgs, risky = param
    items = {}
    for uscg in uscgs :
        tag = os.path.basename(uscg).rsplit('.', 1)[0]
        data = json.load(open(uscg))['OTU']
        for d in data :
            if len(d) < 7 or  not d[6] :
                continue
            for n, s in d[6].items() :
                if risky :
                    s = s.upper()
                g = re.split('__', n)[1]
                key = (tag, d[3], g)
                if key not in items :
                    items[key] = {}
                items[key][n] = s
    return items



@click.command()
@click.argument('uscgs', nargs=-1)
@click.option('-d', '--dbname', help='absolute path of the database. [required]', required=True)
@click.option('-r', '--ref_acc', help='accession in the database to be used as reference. [required]', required=True)
@click.option('--min_identity', help='minumum identity of a sequence comparing to ref_acc. default: 0.94 [0. - 1.]', default=0.94, type=float)
@click.option('--min_presence', help='minumum coverages of genes for a strain to be evaluated. [default: 0.2]', default=0.2, type=float)
@click.option('--risky', help='use low depth sites for tree. default: False', default=False, is_flag=True)
@click.option('-N', '--no_db', help='do not include genomes in the database', default=False, is_flag=True)
@click.option('-g', '--genome_list', help='additional genomes to be included in the analysis', default='')
@click.option('-j', '--json_list', help='list of uscgs as a file [default: None]', default='')
@click.option('-o', '--output', help='prefix of the outputs.', required=True)
@click.option('-n', '--n_proc', help='number of processes [default: 8]', default=8, type=int)
def main(uscgs, dbname, ref_acc, min_identity, min_presence, no_db, genome_list, json_list, output, n_proc, risky) :
    strPhylo(uscgs, dbname, ref_acc, min_identity, min_presence, no_db, genome_list, json_list, output, n_proc, risky)

def strPhylo(uscgs, dbname, ref_acc, min_identity, min_presence, no_db, genome_list, json_list, output, n_proc, risky) :
    pool = Pool(n_proc)
    if min_identity <= 1 :
        min_identity *= 100.

    seqs = readFasta('{0}.ffn.gz'.format(dbname))
    ref_seqs =  {n:s for n, s in seqs.items() if re.split('__', n)[1] == ref_acc}
    if no_db :
        seqs = {}
    
    ref_cnt = len(ref_seqs)
    if json_list :
        uscgs = list(uscgs)
        with open(json_list, 'rt') as fin :
            for line in fin :
                uscgs.append(line.strip().split()[0])
    items = {}
    for items0 in pool.imap_unordered(readJson, [ [uscgs[i::n_proc], risky] for i in np.arange(n_proc) ]) :
        for k, s in items0.items() :
            if len(s) >= min_presence * ref_cnt :
                items[k] = items0[k]

    logging.info('Kept {0} sets of USCGs after initial filtering.'.format(len(items)))

    genomes = []
    if genome_list :
        with open(genome_list, 'rt') as fin :
            for line in fin :
                genomes.append( os.path.abspath(line.strip().split()[0]) )

    with tempfile.TemporaryDirectory(prefix='qry_', dir='.') as tmpdir :
        tre = minimap_align(tmpdir, ref_acc, ref_seqs, seqs, items, genomes, min_identity, min_presence, pool)
        with open(output + '.aln', 'wt') as aln_out, open(os.path.join(tmpdir, 'USCG_align.fas'), 'rt') as fin :
            aln_out.write(fin.read())
    if os.path.isfile('{0}.list'.format(dbname)):
        metadata = {genome:species for genome, species in pd.read_csv('{0}.list'.format(dbname), sep=',', header=None).values}
        for leaf in tre.get_leaves() :
            leaf.annotations['species'] = metadata.get(leaf.name, '')
    clusters = get_tre_cluster(tre, ref_acc)
    for leaf in tre.get_leaves() :
        leaf.annotations['cluster_id'] = clusters.get(leaf.name, 0)
    with open(output + '.nwk', 'wt') as nwk_out, open(output + '.nex', 'wt') as nex_out :
        nwk_out.write(tre.write(format=1)+'\n')
        nex_out.write(ete3_extensions.write_nexus([tre]))


def get_tre_dist(tre) :
    dists = {}
    for n in tre.traverse('postorder') :
        if n.is_leaf() :
            n.d = [[n.name, n.dist]]
        else :
            for i1, c1 in enumerate(n.children) :
                for c2 in n.children[:i1] :
                    for n1, d1 in c1.d :
                        for n2, d2 in c2.d :
                            dists[tuple(sorted([n1, n2]))] = d1 + d2
            n.d = [[x, d+n.dist] for c in n.children for x, d in c.d ]
    return dists


def dist2leiden(dists, cutoff=0.01, nodes={}) :
    import igraph as ig
    edges = []
    if nodes :
        nodes = {n:i for i, n in enumerate(nodes) }
    for (d1, d2), dist in dists.items() :
        if d1 not in nodes :
            nodes[d1] = len(nodes)
        if d2 not in nodes :
            nodes[d2] = len(nodes)
        if dist <= cutoff :
            edges.append([nodes[d1], nodes[d2]])
    # nodes = [[0] for n in sorted(nodes.items(), key=lambda n:n[1])]
    
    if len(edges) :
        g = ig.Graph(len(nodes))
        g.add_edges(edges)
        
        best_module = [-10, None]
        for _ in range(100) :
            m = g.community_leiden(resolution=0.005, n_iterations=-1)
            if m.modularity > best_module[0] :
                best_module = [m.modularity, m]
        return nodes, best_module[1].membership
    else :
        return nodes, np.arange(len(nodes)).tolist()

def get_tre_cluster(tre, ref_acc) :
    from scipy.stats import fisher_exact
    dists = get_tre_dist(tre)
    
    nodes = []
    for (k1, k2), d in dists.items() :
        if k1 == ref_acc :
            nodes.append([d, k2])
        elif k2 == ref_acc :
            nodes.append([d, k1])
    nodes = [ref_acc] + [n[1] for n in sorted(nodes)]
    nodes, clusters = dist2leiden(dists, 0.01, nodes)
    return {n:clusters[i]+1 for n, i in nodes.items()}


if __name__ == '__main__' :
    strPhylo()
