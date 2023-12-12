import os, ete3, pandas as pd, numpy as np, subprocess, click, tempfile, pickle, json, re, time
import ete3_extensions, configure, copy
from scipy.stats import norm

executables = configure.executables
rc = configure.rc

def get_sites(aln_fas, nwk) :
    names = list(aln_fas.keys())
    fas = np.array([list(aln_fas[n]) for n in names]).T
    sites, vseqs = [], []
    for i, s in enumerate(fas) :
        stype = np.unique(s)
        if stype.size > 2 or (stype.size == 2 and stype[0] != '-') :
            sites.append(i+1)
            vseqs.append(s)
    vseqs = [''.join(s) for s in np.array(vseqs).T]
    with tempfile.TemporaryDirectory(prefix='se_', dir='.') as tmpdir :
        with open(os.path.join(tmpdir, 'aln'), 'wt') as fout :
            for n, s in zip(names, vseqs) :
                fout.write('>{0}\n{1}\n'.format(n, s))
        subprocess.Popen('{0} --ancestral -s {1} -te {2} -m GTR+G4 -redo --prefix {3}'.format(configure.executables['iqtree'], os.path.join(tmpdir, 'aln'), nwk, os.path.join(tmpdir, 'aln')).split(),
                         stdout=subprocess.PIPE).communicate()
        tre = ete3.Tree(os.path.join(tmpdir, 'aln.treefile'), format=1)
        nodes = {n.name:[] for n in tre.traverse() if not n.is_leaf() }
        with open(os.path.join(tmpdir, 'aln.state'), 'rt') as fin :
            for line in fin :
                if line.startswith('#') or line.startswith('Node\t') :
                    continue
                p = line.strip().split()
                nodes[p[0]].append([ float(v) for v in p[3:] ])
        for n, s in zip(names, vseqs) :
            nodes[n] = [{'A':[1., 0., 0., 0.], 'C':[0., 1., 0., 0.], 'G':[0., 0., 1., 0.], 'T':[0., 0., 0., 1.]}.get(b, [0., 0., 0., 0.]) for b in s ]
        for n, s in nodes.items() :
            nodes[n] = np.array(s)

        for n in tre.iter_descendants('postorder') :
            if n.is_leaf() :
                nodes[n.name][np.sum(nodes[n.name], 1) == 0] = nodes[n.up.name][np.sum(nodes[n.name], 1) == 0]
    return nodes, [[s, []] for s in sites], tre

def map_qry(aln_fas, otu, sites) :
    with tempfile.TemporaryDirectory(prefix='se_', dir='.') as tmpdir :
        with open(os.path.join(tmpdir, 'ref'), 'wt') as fout :
            n, s = list(aln_fas.items())[0]
            s = s.upper()
            fout.write('>ref\n{1}\n'.format(n, s))
            ref_seq = {'ref':s}
        qry_seq = {}
        with open(os.path.join(tmpdir, 'qry'), 'wt') as fout :
            for n, s in otu[6].items() :
                s = s.upper()
                qry_seq[n] = list(s)
                fout.write('>{0}\n{1}\n'.format(n, s))

        map_cmd = '{0} -k13 -w5 -c -t1 --frag=yes --rmq -A1 -B4 -O8,16 -E2,1 -r20k,40k -g10k -P -N5000 -f1000,5000 -n2 -m20 -s30 -z200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(
            executables['minimap2'], 'ref', 'qry'
        )

        p = subprocess.Popen(map_cmd.split(), cwd=tmpdir, stdout=subprocess.PIPE, universal_newlines=True) #, stderr=subprocess.PIPE, )
        for line in p.stdout :
            if line.startswith('[') :
                continue
            p = line.strip().split('\t')
            p[9:11] = int(p[9]), 100. * float(p[9])/float(p[10])
            p[1:4] = [int(p[1]), int(p[2])+1, int(p[3])]
            p[6:9] = [int(p[6]), int(p[7])+1, int(p[8])]

            if p[4] == '+' :
                qi, ri, cigar, d = p[2], p[7], p[-1][5:], 1 
            else :
                qi, ri, cigar, d = p[3], p[7], p[-1][5:], -1
            xi = 0
            while xi < len(sites) and sites[xi][0] < ri :
                xi += 1
            for s, t in re.findall('(\d+)([MDI])', cigar) :
                s = int(s)
                if t != 'I' :
                    rj = ri + s
                if t != 'D' :
                    qj = qi + s*d
                while xi < len(sites) and sites[xi][0] >= ri and sites[xi][0] < rj :
                    rd = sites[xi][0] - ri
                    rx = ri + (rd if t != 'I' else 0) - 1
                    qx = qi + (rd*d if t != 'D' else 0) - 1
                    rseq = ref_seq[p[5]][rx]
                    qseq = qry_seq[p[0]][qx] if d > 0 else configure.rc(qry_seq[p[0]][qx])
                    if qseq != '-' :
                        rr = sites[xi][1]
                        if len(rr) == 0 or rr[5] < p[10] :
                            sites[xi][1] = [p[0], qx+d, rseq, qseq, p[9], p[10], p[4]]
                    xi += 1
                ri, qi = rj, qj
            qry_seq[p[0]][p[2]-1:p[3]] = ['-']* (p[3] - p[2] + 1)
    return sites

def parse_bam(bam, ref_acc, sites) :
    base_comp = {}
    if bam :
        qry_sites = {}
        for site, var in sites :
            if len(var) :
                qry_sites[(var[0], var[1])] = ('ref', site)

        p = subprocess.Popen('{0} mpileup -AB -q 0 -Q 0 {1}'.format(executables['samtools'], bam).split(),
                            universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        contigs = {}

        for line in p.stdout :
            p = line.strip().split('\t')
            if re.split('__', p[0])[1] != ref_acc :
                continue
            if p[0] not in contigs :
                contigs[p[0]] = int(p[1]) - 1

            s = re.sub(r'\^.', '', p[4]).upper().replace('$', '')
            s = list(''.join([b[int(n):] for n, b in re.findall('[+-](\d+)(.+)', '+0' + s)]))
            base, cdp = min(zip(*np.unique(s, return_counts=True)), key=lambda x:-x[1])
            if base in ('*', 'N') :
                contigs[p[0]] -= 1
                continue
            
            site = int(p[1]) - contigs[p[0]]
            key = (p[0], site)
            if key not in qry_sites :
                continue
            
            bases = dict(zip(*np.unique(s, return_counts=True)))
            base_comp[key] = [int(bases.get(b, 0.)) for b in ('A', 'C', 'G', 'T')]

    for i, (site, v) in enumerate(sites) :
        if len(v) :
            x = base_comp.get((v[0], v[1]), [0, 0, 0, 0])
            sites[i].append(x if v[6] == '+' else list(x[::-1]))
        else :
            sites[i].append([0, 0, 0, 0])
    return sites

def updateGenotypes(genotypes, genolocs, seqs, branches, neighbors, nChoice=1, big_move=0.1) :
    genotypes = np.copy(genotypes)
    genolocs = np.copy(genolocs)
    seqs = np.copy(seqs)

    ids = np.random.choice(genotypes.shape[0], replace=False, size=nChoice)
    for i in ids :
        if np.random.rand() >= big_move :
            sets = set(neighbors[genotypes[i]]) 
        else : 
            sets = set(np.arange(len(neighbors)))
        g = np.random.choice(list(sets - set([n for x in np.arange(genotypes.shape[0]) if x != i for n in neighbors[genotypes[x]]])))
        genotypes[i] = g
        if g >= len(neighbors)-1 :
            genolocs[i] = 1.
            ss = branches[g]
        else :
            genolocs[i] = np.random.rand()
            ss = branches[g]*genolocs[i] + branches[neighbors[g][0]]*(1-genolocs[i])
            seqs[i] = np.array([np.random.choice(np.arange(4), p = s/np.sum(s)) for s in ss])
    return genotypes, genolocs, seqs


def updateOnlyGenotypes(genotypes, genolocs, seqs, branches, neighbors) :
    seqs = np.copy(seqs)

    for i in np.random.choice(genotypes.shape[0], size=1) :
        g = genotypes[i]
        if g >= len(neighbors)-1 :
            ss = branches[g]
        else :
            x = genolocs[i]
            ss = branches[g]*x + branches[neighbors[g][0]]*(1-x)
            seqs[i] = np.array([np.random.choice(np.arange(4), p = s/np.sum(s)) for s in ss])
    return seqs



def getLikelihood(seqs, covs, sigma, sample, w) :
    seq = np.zeros([seqs.shape[1], 4], dtype=float)
    for s, c in zip(seqs, covs/np.sum(covs)) :
        seq[np.arange(seq.shape[0]), s] += c
    return np.sum(norm.logpdf(np.abs(sample - seq), 0, sigma)*w)


def estimate(nodes, sites, tre, max_nGenotype=4) :
    sample = nodes['__query__'].astype(float)
    cov = np.sum(sample, 1)
    max_cov = np.max([np.quantile(cov, 0.75) + 3.*(np.quantile(cov, 0.75) - np.quantile(cov, 0.25)), 6.])
    
    sample_idx = (cov > 0) & (cov <= max_cov)
    cov = cov[sample_idx].reshape([-1, 1])
    sample = sample[sample_idx]/cov

    neighbors = [set([i]) for i, n in enumerate(tre.traverse())]
    branches = np.array([nodes[n.name][sample_idx] for n in tre.traverse('postorder')])
    for i, n in enumerate(tre.traverse('postorder')) :
        n.i = i
        if not n.is_leaf() :
            group = {c.i for c in n.children}
            for c in n.children :
                neighbors[c.i] = [i] + sorted(neighbors[c.i] | group)
            neighbors[n.i] |= group
    neighbors[-1] = sorted(neighbors[-1])
    nBranch, nSite, nType = branches.shape
    
    n_map = {n.i:n.name for n in tre.traverse('postorder')}

    bestModel = {'aic':999999999999.}

    for nGenotype in np.arange(max_nGenotype, 0, -1) :
        if nGenotype == max_nGenotype :
            globalParams = dict(
                sigma = 1., 
                pGenotype = np.random.rand(nGenotype), 
                genotypes = np.zeros(nGenotype, dtype=int), 
                genolocs = np.random.rand(nGenotype), 
                seqs = np.zeros([nGenotype, nSite], dtype=int), 
                lk = -999999999999., 
                aic = 999999999999.,
            )
            globalParams['pGenotype'] = sorted(globalParams['pGenotype'])/np.sum(globalParams['pGenotype'])
            if min(globalParams['pGenotype']) < 0.02 :
                globalParams['pGenotype'] += (0.02 - min(globalParams['pGenotype']))/(1-0.02*globalParams['pGenotype'])
                globalParams['pGenotype'] = globalParams['pGenotype']/np.sum(globalParams['pGenotype'])
                
            for i in np.arange(1, nGenotype) :
                globalParams['genotypes'][i] = np.random.choice(list(set(range(nBranch)) - {x for g in globalParams['genotypes'][:i] for x in neighbors[g]}))
            
            globalParams['genotypes'], globalParams['genolocs'], globalParams['seqs'] = updateGenotypes(globalParams['genotypes'], globalParams['genolocs'], globalParams['seqs'], branches, neighbors, nChoice=nGenotype)
        else :
            idx = np.ones(nGenotype+1, dtype=bool)
            idx[np.argmin(globalParams['pGenotype'])] = False
            globalParams['pGenotype'] = globalParams['pGenotype'][idx]
            globalParams['genotypes'] = globalParams['genotypes'][idx]
            globalParams['genolocs'] = globalParams['genolocs'][idx]
            globalParams['seqs'] = globalParams['seqs'][idx]
            globalParams['pGenotype'] = (globalParams['pGenotype'])/np.sum(globalParams['pGenotype'])
            
        globalParams['lk'] = getLikelihood(globalParams['seqs'], globalParams['pGenotype'], globalParams['sigma'], sample, cov)
        globalParams['aic'] = 1 + 2*(nGenotype-1)*2 - 2*globalParams['lk']
        if globalParams['aic'] < bestModel['aic'] :
            bestModel = copy.deepcopy(globalParams)
        
        nChain = 500 * nGenotype
        for ite in range(nChain) :
            # update initState
            genotypes, genolocs, seqs = updateGenotypes(globalParams['genotypes'], globalParams['genolocs'], globalParams['seqs'], branches, neighbors, nChoice=1 if nGenotype == 1 else np.random.choice(int(nGenotype/2.))+1)
            newLK = getLikelihood(seqs, globalParams['pGenotype'], globalParams['sigma'], sample, cov)
            if newLK >= globalParams['lk'] or np.random.rand() < np.exp(newLK-globalParams['lk']) :
                globalParams['genotypes'][:] = genotypes
                globalParams['genolocs'][:] = genolocs
                globalParams['seqs'][:] = seqs
                globalParams['lk'] = newLK
                globalParams['aic'] = 1 + 2*(nGenotype-1)*2 - 2*globalParams['lk']
                if globalParams['aic'] < bestModel['aic'] :
                    bestModel = copy.deepcopy(globalParams)

            seqs = updateOnlyGenotypes(globalParams['genotypes'], globalParams['genolocs'], globalParams['seqs'], branches, neighbors)
            newLK = getLikelihood(seqs, globalParams['pGenotype'], globalParams['sigma'], sample, cov)
            if newLK >= globalParams['lk'] or np.random.rand() < np.exp(newLK-globalParams['lk']) :
                globalParams['seqs'][:] = seqs
                globalParams['lk'] = newLK
                globalParams['aic'] = 1 + 2*(nGenotype-1)*2 - 2*globalParams['lk']
                if globalParams['aic'] < bestModel['aic'] :
                    bestModel = copy.deepcopy(globalParams)

            if globalParams['pGenotype'].size > 1 :
                pGenotype = np.copy(globalParams['pGenotype'])
                id = np.random.choice(nGenotype)
                pGenotype[id] = (np.random.rand()-0.5) + pGenotype[id]
                pGenotype = pGenotype/np.sum(pGenotype)
                if min(pGenotype) < 0.02 :
                    pGenotype += (0.02 - min(pGenotype))/(1 - 0.02*nGenotype)
                    pGenotype = pGenotype/np.sum(pGenotype)
                
                newLK = getLikelihood(globalParams['seqs'], pGenotype, globalParams['sigma'], sample, cov)
                if newLK >= globalParams['lk'] or np.random.rand() < np.exp(newLK-globalParams['lk']) :
                    globalParams['pGenotype'][:] = pGenotype
                    globalParams['lk'] = newLK
                    globalParams['aic'] = 1 + 2*(nGenotype-1)*2 - 2*globalParams['lk']
                    if globalParams['aic'] < bestModel['aic'] :
                        bestModel = copy.deepcopy(globalParams)
            if ite % 50 == 0 :
                print(ite, globalParams['aic'], bestModel['aic'], {'genolocs':bestModel['genolocs'], 'pGenotype':bestModel['pGenotype'], 'sigma':bestModel['sigma'], 'nodes':[n_map[i] for i in bestModel['genotypes']]})
    seqs = np.zeros([bestModel['pGenotype'].size, sample_idx.size], dtype=int)
    seqs[:] = -1
    seqs[:, sample_idx] = bestModel['seqs']
    bestModel['seqs'] = seqs
    return bestModel


def seperate_strains(bam, ref_acc, sites, best_model) :
    best_model['seqs'] = best_model['seqs'].astype(int)
    
    nGenotype = best_model['pGenotype'].size
    
    qry_sites = {}
    for idx, (site, var, bases) in enumerate(sites) :
        if len(var) :
            qry_sites[(var[0], var[1])] = best_model['seqs'][:, idx] if var[6] == '+' else 3-best_model['seqs'][:, idx]

    p = subprocess.Popen('{0} mpileup -AB -q 0 -Q 0 {1}'.format(executables['samtools'], bam).split(),
                        universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    contigs = {}
    seqs = [{} for _ in np.arange(nGenotype)]
    for line in p.stdout :
        p = line.strip().split('\t')
        if re.split('__', p[0])[1] != ref_acc :
            continue

        if p[0] not in contigs :
            contigs[p[0]] = int(p[1]) - 1
            for s in seqs :
                s[p[0]] = []

        s = re.sub(r'\^.', '', p[4]).upper().replace('$', '')
        s = list(''.join([b[int(n):] for n, b in re.findall('[+-](\d+)(.+)', '+0' + s)]))
        cons, cdp = min(zip(*np.unique(s, return_counts=True)), key=lambda x:-x[1])
        if cons in ('*', 'N') :
            contigs[p[0]] -= 1
            continue

        bases = dict(zip(*np.unique(s, return_counts=True)))
        bases = np.array([int(bases.get(b, 0.)) for b in ('A', 'C', 'G', 'T')])
        n_bases = np.sum(bases)
        bi = set(np.where(bases>0)[0])

        site = int(p[1]) - contigs[p[0]]
        key = (p[0], site)
        if key in qry_sites :
            bi |= set(qry_sites[key].tolist())
        bi = np.array(list(bi))
        
        max_comb = None
        for i in np.arange(bi.size**nGenotype) :
            pp_idx = np.array([ bi[int(i/(bi.size**pi))%bi.size] for pi in np.arange(nGenotype) ])
            base2 = np.zeros(4)
            for pp, pt in zip(best_model['pGenotype'], pp_idx) :
                base2[pt] += pp
            base2 *= n_bases
            diff = np.sum(np.abs(bases-base2))*0.5
            same = n_bases - diff
            prob = np.log(0.01)*diff + np.log(0.97)*same
            if key in qry_sites :
                for gt, pt in zip(qry_sites[key], pp_idx) :
                    if gt == pt :
                        prob += np.log(9997.)

            if not max_comb or prob > max_comb[0] :
                max_comb = [prob, pp_idx, same/n_bases]
        for i, (b, c) in enumerate(zip(max_comb[1], best_model['pGenotype'] * n_bases)) :
            if max_comb[2] < 0.8 :
                c = 0
            elif bases[b] < c :
                c = bases[b]
            seqs[i][p[0]].append('ACGT'[b] if c>=3 else 'acgt'[b])
            
    return [ {n:''.join(s) for n, s in seq.items()} for seq in seqs ]



@click.command()
@click.option('-n', '--nwk', help='nwk file')
@click.option('-a', '--aln', help='aln file')
@click.option('-u', '--uscg', help='uscg file')
@click.option('-p', '--prefix', help='prefix for the output')
@click.option('-b', '--bam', help='bam file')
def explore(nwk, aln, uscg, prefix, bam) :
    aln_fas = configure.readFasta(aln)
    nodes, sites, tre = get_sites(aln_fas, nwk)
    uscgs = json.load(open(uscg))
    t = ete3.Tree(nwk, format=1)
    
    res = {'OTU':[], 'profile':[]}
    for g in tre.get_leaf_names() :
        if g.find('|') <= 0 :
            continue
        ref_acc = g.split('|')[-1]
        
        otu = []
        for o in uscgs['OTU'] :
            if ref_acc in o[4] :
                otu.extend(o[:6])
                otu.append({})
                for n, s in o[6].items() :
                    if n.find(ref_acc) > 0 :
                        otu[6][n] = s

        sites = map_qry(aln_fas, otu, sites)
        sites = parse_bam(bam, ref_acc, sites)
        nodes['__query__'] = np.array([s[2] for s in sites])
        best_model = estimate(nodes, sites, tre)
        # pickle.dump([nodes, sites, tre, best_model], open('tmp.dump', 'wb'))
        # nodes, sites, tre, best_model = pickle.load(open('tmp.dump', 'rb'))
        if best_model['pGenotype'].size == 1 :
            res['OTU'].append(otu)
        else :
            seqs = seperate_strains(bam, ref_acc, sites, best_model)
            otus = []
            for sid, (p, seq) in enumerate(zip(best_model['pGenotype'], seqs)) :
                new_acc = ref_acc+'.{0}'.format(sid+1)
                otus.append([otu[0]*p, otu[1], otu[2], otu[3], [new_acc], otu[5], {}])
                for n, s in seq.items() :
                    nn = re.split('__', n)
                    n = '__'.join([nn[0], new_acc, nn[2]])
                    otus[-1][6][n] = s
            res['OTU'].extend(otus)
    json.dump(res, open(prefix + '.resolved.json', 'wt'))
      
if __name__ == '__main__' :
    explore()
