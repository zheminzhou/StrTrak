import os, click, numpy as np, pandas as pd, gzip, json
import subprocess, re
from multiprocessing import Pool
import configure

executables = configure.executables
logging = configure.logging
readFasta = configure.readFasta



def get_matches(metadata, uscg_info, genome_ids, read_maps, min_gene_match=3, block_size=300, single_sequence=False) :
    genome_info, gene_info = {}, {}
    for gene, (s_id, g_id, size, o_id) in uscg_info.items() :
        if s_id not in genome_info :
            genome_info[s_id] = []
        genome_info[s_id].append([g_id, size, o_id])
        gene_info[g_id] = o_id
    
    '''genome_id, gene_id, read_id, start, size, mutation*6, diff*6, cur_diff*6
        0            1       2       3      4        5        6         7   '''
    genomes = dict(zip(*np.unique(read_maps[read_maps.T[6] < 15, 0], return_counts=True)))
    genomes = { g:c for g, c in genomes.items() if c >= min_gene_match and g in genome_info }
    read_maps = read_maps[[g in genomes for g in read_maps.T[0]]] # read_maps[np.in1d(read_maps.T[0], list(genomes.keys()))]
    read_maps.T[3] -= 1
    if read_maps.shape[0] <= 0 :
        return [], []
    
    max_block = int(np.ceil(np.max([size for g in genomes for gene, size, ortho in genome_info[g]])/block_size))
    
    genome_info2 = {}
    regions = {}
    for g in genomes :
        genome_info2[g] = []
        for gene, size, ortho in genome_info[g] :
            for i in np.arange(np.ceil(size / block_size).astype(int)):
                regions[gene*max_block + i] = len(regions)
                region_size = (block_size if i > 0 else block_size+100) if (i+1)*block_size <= size else (i+1)*block_size-size + 100
                region_id = len(regions)-1
                if single_sequence :
                    genome_info2[g].append([region_id, region_size, region_id, ortho])
                else :
                    genome_info2[g].append([region_id, region_size, gene, ortho])

    genome_info2 = {g:np.array(info) for g, info in genome_info2.items()}
    if single_sequence :
        genome_info = genome_info2
    
    coords = np.vstack([(read_maps.T[3]/block_size).astype(int), ((read_maps.T[3]+read_maps.T[4]-1)/block_size).astype(int)]).T
    coords = np.random.randint(coords.T[0], coords.T[1]+1)
    read_maps.T[3, :] = read_maps.T[1]
    read_maps.T[1, :] = np.vectorize(regions.get)(read_maps.T[1]*max_block + coords)
    read_maps = np.hstack([read_maps, np.zeros([read_maps.shape[0], 1], dtype=np.int32)])

    match_results = []
    summed_reads = np.empty([np.max(read_maps.T[2])+1, 4], dtype=np.int32)
    summed_reads[:, 1:].fill(99999)
    summed_reads[:, 0].fill(-1)
    coverages = [ [-1, -1, genome, []] for genome in np.unique(read_maps.T[0]) ]
    while len(coverages) > 0 :
        gene_cov = np.bincount(read_maps.T[1], weights=np.power(0.368, read_maps.T[6].astype(np.float64) / 6.))
        max_i = -1
        for i, (depth, n_gene, genome, good_genes) in enumerate(coverages) :
            if depth == -1 or max_i < 0 or depth >= coverages[max_i][0] :
                genes = genome_info2[genome].astype(int)
            else :
                break
            covs = np.array([ [((gene_cov[g]+.5)/(s+.5)), gene_cov[g],s, og] \
                                  if g < gene_cov.size else [0., 0, s, 0] for g, s, og, oo in genes ])
            n_gene = np.unique(covs[covs.T[1] > 0.5, 3]).size
            if n_gene < min_gene_match and n_gene < covs.shape[0]*0.5 :
                coverages[i] = [0, n_gene, genome, np.zeros(0, dtype=int)]
                continue
            cv2 = covs[covs.T[2]>=0.5*block_size] if np.sum(covs.T[2]>=0.5*block_size) > 1 else covs
            q1, q3 = np.quantile(cv2.T[0], 0.25), np.quantile(cv2.T[0], 0.75)
            delta_q = max(q3 - q1, 1./block_size)
            idx = (covs.T[1]/covs.T[2] <= q3 + 3*delta_q)
            covs = covs[ idx ]
            if len(covs) < 1 :
                coverages[i] = [0, n_gene, genome, np.zeros(0, dtype=int)]
                continue
            cov = np.sum(covs.T[1]) / np.sum(covs.T[2])
            n_gene = np.unique(covs[covs.T[1] > 0.5, 3]).size
            coverages[i] = [cov, n_gene, genome, genes[idx, 0]]
            if (max_i < 0 or cov > coverages[max_i][0]) and n_gene > min_gene_match :
                max_i = i
        (depth, n_gene, match, good_genes) = coverages[max_i]
        coverages.sort(reverse=True)
        if n_gene < min_gene_match and n_gene < 0.5 * len(genome_info[match]) :
            break

        new_matches, to_del = [], []
        matched_reads = np.empty([summed_reads.shape[0], 3], dtype=int)
        matched_reads.fill(99999)

        x = read_maps[read_maps.T[0] == match]
        matched_reads[x.T[2], :] = x[:, 4:7]
        matched_ortho = set(np.unique(np.vectorize(gene_info.get)(np.unique(x.T[3]))))
        dist = 0.0001 + np.sum(matched_reads[matched_reads.T[2] <= 12, 1])/np.sum(matched_reads[matched_reads.T[2] <= 12, 0])/6. \
                    if np.sum(matched_reads[matched_reads.T[2] <= 12, 0]) > 0 else 0.01
        dist2 = 0.0001 + np.sum(matched_reads[(matched_reads.T[2] > 0) & (matched_reads.T[2] < 1000), 1])/np.sum(matched_reads[(matched_reads.T[2] > 0) & (matched_reads.T[2] < 1000), 0])/6. \
                    if np.sum(matched_reads[(matched_reads.T[2] > 0) & (matched_reads.T[2] < 1000), 0]) > 0 else 0.05

        read_maps.T[7] = matched_reads[read_maps.T[2], 2] - read_maps.T[6]

        total_reads = np.bincount(read_maps.T[0])
        m_reads = np.bincount(read_maps[read_maps.T[7] < 1000, 0], minlength=total_reads.size)

        similar_genomes = list(set(np.where(m_reads > total_reads * 0.75)[0]) - set([match]) - set(to_del))
        if len(similar_genomes) :
            all_reads = read_maps[np.in1d(read_maps.T[0], similar_genomes)]
            all_reads = all_reads[[ o in matched_ortho for o in np.vectorize(gene_info.get)(all_reads.T[3])]]
            total_reads = np.bincount(all_reads.T[0], minlength=max(similar_genomes) + 1)
            strict_pos2 = all_reads[(all_reads.T[7] > 12)]
            strict_pos = all_reads[(all_reads.T[7] > 12) & (all_reads.T[7] < 1000)]
            pos_reads = all_reads[(all_reads.T[7] >= 0) & (all_reads.T[7] < 1000)]
            unq_reads = all_reads[all_reads.T[6] <= 12]
            neg_reads = all_reads[all_reads.T[6] > 0]
            shared_dists = dict(zip(similar_genomes, np.vstack([ np.bincount(pos_reads.T[0], weights=pos_reads.T[4], minlength=max(similar_genomes) + 1),
                                                                 np.bincount(pos_reads.T[0], weights=pos_reads.T[5], minlength=max(similar_genomes) + 1) / 6.,
                                                                 np.bincount(neg_reads.T[0], weights=neg_reads.T[4], minlength=max(similar_genomes) + 1),
                                                                 np.bincount(neg_reads.T[0], weights=neg_reads.T[5], minlength=max(similar_genomes) + 1) / 6.,
                                                                 np.bincount(unq_reads.T[0], weights=unq_reads.T[4], minlength=max(similar_genomes) + 1),
                                                                 np.bincount(unq_reads.T[0], weights=unq_reads.T[5], minlength=max(similar_genomes) + 1) / 6.,
                                                                 np.bincount(strict_pos.T[0], minlength=max(similar_genomes) + 1),
                                                                 np.bincount(strict_pos2.T[0], minlength=max(similar_genomes) + 1),
                                                             ]).T[list(similar_genomes)]))

            for g, n in sorted(shared_dists.items(), key=lambda n:(n[1][5]+0.01)/(n[1][4]+0.01)) :
                pd = 0.0001+n[1]/n[0] if n[0] else 1.
                nd = 0.0001+n[3]/n[2] if n[2] else 1.
                ud = 0.0001+n[5]/n[4] if n[4] else 1.
                if (m_reads[g] >= np.sum(matched_reads.T[2]<1000) * 0.75) and ((pd/dist <= dist2/pd) or (ud/dist <= dist2/ud)) and np.unique(strict_pos2[strict_pos2.T[0] == g, 1]).size >= min_gene_match :
                    if pd < dist or ud < dist :
                        new_matches.append(match)
                        match, g, dist, dist2 = g, match, min(ud, pd), nd
                        continue
                    elif (n[6] >= m_reads[g] * 0.05) or (n[7] >= total_reads[g] * 0.05) :
                        new_matches.append(g)
                        continue
                if (n[6] < m_reads[g] * 0.05) or (n[7] < total_reads[g] * 0.05) :
                    to_del.append(g)

        if len(new_matches) :
            x = read_maps[np.in1d(read_maps.T[0], new_matches)]
            x = x[np.argsort(x.T[6])]
            x = x[np.unique(x.T[2], return_index=True)[1]]
            x.T[7] = matched_reads[x.T[2], 2]
            x = x[x.T[7] > x.T[6]]
            matched_reads[x.T[2], :] = x[:, 4:7]
        matches = [match] + new_matches

        if np.sum(matched_reads.T[2] <= 3) >= min_gene_match :
            gene_cov = np.bincount(read_maps[read_maps.T[7] < 21, 1]) #np.bincount(read_maps.T[1])
            good_genes = []
            for gid, genome in enumerate(matches) :
                genes = genome_info2[genome].astype(int)
                covs = np.array([ [((gene_cov[g] + .5)/(s + .5)), gene_cov[g],s, og] \
                                  if g < gene_cov.size else [0., 0, s, 0] for g, s, og, oo in genes ])

                cv2 = covs[covs.T[2]>=0.5*block_size] if np.sum(covs.T[2]>=0.5*block_size) > 1 else covs
                q1, q3 = np.quantile(cv2.T[0], 0.25), np.quantile(cv2.T[0], 0.75)
                delta_q = max(q3 - q1, 1./block_size)
                idx = (covs.T[1] > 0) & (covs.T[1]/covs.T[2] <= q3 + 3*delta_q)
                if gid == 0 :
                    remained = np.unique(genes[idx, 2]).size
                good_genes.extend(genes[idx, 0].tolist())
            bad_reads = list(set(np.where(matched_reads.T[2] < 1000)[0]) - set(read_maps[np.in1d(read_maps.T[1], good_genes), 2]))
            matched_reads[bad_reads, :] = 99999
            if remained > 0 :
                match_results.append([matches, depth])
                idx = matched_reads.T[2] < summed_reads.T[3]
                summed_reads[idx, 0] = len(match_results)-1
                summed_reads[idx, 1:] = matched_reads[idx, :]
        read_maps = read_maps[(read_maps.T[7] >= 3) & np.in1d(read_maps.T[0], to_del, invert=True)]
        mdel = set(matches) | set(to_del)
        coverages = [ c for c in coverages if c[2] not in mdel and c[0] != 0 ]


    for m in np.unique(summed_reads[summed_reads.T[0] >= 0, 0]):
        p = summed_reads[summed_reads.T[0] == m, 2] / summed_reads[summed_reads.T[0] == m, 1]
        q1, q3 = np.quantile(p, 0.25), np.quantile(p, 0.75)
        delta_q = max(q3-q1, 0.1)
        scope = q3 + 3*delta_q
        summed_reads[(summed_reads.T[0] == m) & (summed_reads.T[2] / summed_reads.T[1] > scope), :] = -1

    grp_cnt = np.bincount(summed_reads[summed_reads.T[0] >= 0, 0])
    match_results = np.array([tuple(match_results[i][0]) for i in np.where(grp_cnt >= min_gene_match)[0]], dtype=object)
    i1 = 0
    for i0, t in enumerate(grp_cnt >= min_gene_match) :
        if t == False :
            summed_reads[summed_reads.T[0] == i0, :] = [-1, 99999, 99999, 99999]
        else :
            if i0 != i1 :
                summed_reads[summed_reads.T[0] == i0, 0] = i1
            i1 += 1

    return match_results, summed_reads



def generate_outputs(sam_files, metadata, uscg_info, genome_ids, file_reads, matches, reads, tmpdir, output) :
    primary_refs = { match[0]:i for i, match in enumerate(matches) }
    uscg_map = {k:{ v[0]:v[2] } for k, v in uscg_info.items() if  v[0] in primary_refs}
    read_matches = { i:r[0] for i, r in enumerate(reads) }
    results = [[] for match in matches]
    species_res = {}
    for i, match in enumerate(matches) :
        match_reads = reads[reads.T[0] == i]
        n_bases = np.sum(match_reads.T[1]) + 1e-6
        n_diffs = np.sum(match_reads.T[2])/6.
        references = [genome_ids[m][0] for m in match]
        species = [ metadata[r] for r in references ]
        species = [species[0]] + list(set(species[1:]) - set([species[0]]))

        results[i] = [float(match_reads.shape[0]),
                      float(n_bases),
                      100. - float(n_diffs/n_bases*100),
                      species[0] + ('' if len(species) == 1 else ' ({0})'.format(','.join(species[1:]))),
                      references, '']
        if species[0] not in species_res :
            species_res[species[0]] = results[i][:3] + [species[0], [references[0]], '']
        else :
            species_res[species[0]][0] += results[i][0]
            species_res[species[0]][2] = (species_res[species[0]][2]*species_res[species[0]][1] + results[i][2]*results[i][1])/(species_res[species[0]][1]+results[i][1])
            species_res[species[0]][1] += results[i][1]
            species_res[species[0]][4].append(references[0])
    results = [ r for r in results if r[0] > 0 ]
    species_res = { s:r for s, r in species_res.items() if r[0] > 0 }
    for r in results :
        r[1] = float(r[1])/r[0] if r[0] > 0 else 0.
    for s, r in species_res.items() :
        r[1] = float(r[1])/r[0] if r[0] > 0 else 0.
    outputs = {'profile':sorted(species_res.values(), reverse=True), 'OTU':sorted(results, reverse=True)}
    json.dump(outputs, open(output+'.json', 'wt'))
    with gzip.open('{0}/primary.sam.gz'.format(tmpdir), 'wt') as pout, gzip.open('{0}.group_reads.gz'.format(output), 'wt') as sout :
        read_name = ''
        lxx = {}
        for fn_id, (fname, read_id) in enumerate(zip(sam_files, file_reads)) :
            p = subprocess.Popen('{pigz} -cd {0}'.format(fname, **executables).split(), cwd=tmpdir, stdout=subprocess.PIPE,
                                 universal_newlines=True)
            for i, line in enumerate(p.stdout) :
                if line.startswith('@') :
                    if fn_id > 0 :
                        continue
                    if line.startswith('@SQ') :
                        p = line.split('\t')
                        refname = p[1][3:]
                        ref_ids = uscg_map.get(refname, {})
                        for ref_id, g in ref_ids.items() :
                            lx = '\t'.join(p)
                            if lx not in lxx :
                                pout.write(lx)
                                lxx[lx] = 1
                    else :
                        pout.write(line)
                else :
                    p = line.split('\t')
                    if int(p[1]) >= 256 :
                        p[1] = str(int(p[1]) - 256)
                    if p[2] not in uscg_info :
                        continue
                    if p[0] != read_name :
                        read_name, read_id = p[0], read_id+1
                        if read_id >= len(read_matches) :
                            break
                        read_grp = read_matches[read_id]
                        if read_grp >= 0 :
                            sout.write('@{0}_{1}\n{2}\n+\n{3}\n'.format(read_grp, p[0], p[9], p[10]))
                    if p[2] not in uscg_map :
                        continue
                    read_grp = read_matches[read_id]
                    if read_grp >= 0 :
                        p[0] = '{0}_{1}'.format(read_grp, p[0])
                        for ref_id, g in uscg_map[p[2]].items() :
                            if primary_refs[ref_id] == read_grp :
                                pout.write('\t'.join(p))
            # try :
            #     os.unlink(os.path.join(tmpdir, fname))
            # except :
            #     pass
    subprocess.Popen('{pigz} -cd {0}/primary.sam.gz | {samtools} sort -m 4G -@ 8 -O bam -l 0 -T {0}/tmp - > {1}.primary.bam'.format(
        tmpdir, output, **executables
    ), shell=True).communicate()

    return outputs, '{0}.primary.bam'.format(output)


def get_uscg_info(sam_file, metadata) :
    p = subprocess.Popen('{samtools} view -H {0}'.format(sam_file, **executables).split(), stdout=subprocess.PIPE, universal_newlines=True)
    uscg_info, s_id, o_id = {}, {}, {}
    for line in p.stdout :
        if not line.startswith('@') :
            break
        if line.startswith('@SQ') :
            p = line.strip().split('\t')
            gene, size = p[1][3:], int(p[2][3:])
            ortho, strain, tag = re.split('__', gene)
            if strain not in metadata :
                continue
            if strain not in s_id :
                s_id[strain] = len(s_id)
            if ortho not in o_id :
                o_id[ortho] = len(o_id)
            uscg_info[gene] = [s_id[strain], len(uscg_info), size, o_id[ortho]]
    return uscg_info, s_id


def map_to_uscgs(sam_files, uscg_info, tmpdir, max_dist, pool) :
    read_maps, read_id = [], 0
    file_ids = []
    for rmap_file in pool.map(parse_sam, [ [sfile, tmpdir, uscg_info, max_dist] for sfile in sam_files ] ) :
        x = np.load(rmap_file)['reads']
        if x.shape[0] > 0 :
            x.T[2] += read_id
            file_ids.append(read_id)
            read_maps.append(x)
            read_id = x[-1, 2]
        try :
            os.unlink(rmap_file)
        except :
            pass
    if len(read_maps) :
        read_maps = np.vstack(read_maps)
    return read_maps, file_ids



def parse_sam(data) :
    outfile, tmpdir, uscg_info, max_dist = data

    prev, read_id = ['', 0], 0
    read_maps, rmaps = [], []

    p = subprocess.Popen('{pigz} -cd {0}'.format(outfile, **executables).split(), stdout=subprocess.PIPE, universal_newlines=True)
    for lid, line in enumerate(p.stdout):
        if line.startswith('@'):
            continue
        p = line.strip().split('\t')
        if p[2] not in uscg_info:
            continue
        score = int(p[12][5:])
        matches = [ [int(mat[0]), mat[1]] for mat in re.findall(r'(\d+)([SMD])', p[5])]
        mlen = sum([ mat[0] for mat in matches if mat[1] != 'S' ])
        if matches[-1][1] == 'S' :
            mlen += min(int(matches[-1][0]), uscg_info[p[2]][2] - (int(p[3]) + mlen - 1))
        if matches[0][1] == 'S' :
            mlen_s = min(int(matches[0][0]), int(p[3])-1)
            mlen += mlen_s
            p[3] = int(p[3]) - mlen_s
        if p[0] == prev[0]:
            diff = prev[1] - score
        else:
            if len(rmaps):
                rmaps = np.array(rmaps, dtype=int)
                rmaps = rmaps[np.lexsort((np.random.rand(rmaps.shape[0]), rmaps.T[6], rmaps.T[0]))]
                read_maps.extend([m for m in rmaps[np.concatenate([[True], rmaps.T[0][1:] != rmaps.T[0][:-1]])]])
                rmaps = []
            prev, diff = [p[0], score], 0
            read_id += 1
        if 2*mlen-score <= max_dist * mlen * 6 :
            matches = uscg_info[p[2]]
            rmaps.append( matches[:2] + [read_id, int(p[3]), mlen, 2*mlen - score, diff] )
    if len(rmaps):
        rmaps = np.array(rmaps, dtype=int)
        rmaps = rmaps[np.lexsort((np.random.rand(rmaps.shape[0]), rmaps.T[6], rmaps.T[0]))]
        read_maps.extend([m for m in rmaps[np.concatenate([[True], rmaps.T[0][1:] != rmaps.T[0][:-1]])]])
    np.savez_compressed( os.path.join(tmpdir, '{0}.npz'.format(outfile)), reads=np.array(read_maps, dtype=np.int32))
    return os.path.join(tmpdir, '{0}.npz'.format(outfile))




def map_reads(query, uscg_db, tmpdir, max_dist) :
    outputs = []
    for qid, qry in enumerate(query):
        out_sam = 'uscg.{0}.sam.gz'.format(qid)
        if qid == 0 :
            subprocess.Popen(
                '{minimap2} -t16 -ax sr --sr --frag=yes -A2 -B4 -O8,16 -E2,1 -r50 -p0.1 -N20000 -f2000,10000 -Y -n2 -m20 -s40 -g200 --end-bonus 12 -2K50m --heap-sort=yes --secondary=yes --sam-hit-only {0} {1}|{EnFlt} {3}|{pigz} -c > {2}'.format(
                    uscg_db, qry, out_sam, max_dist, **executables,
                ), cwd=tmpdir, shell=True).communicate()
        else :
            subprocess.Popen(
                '{minimap2} -t16 -ax sr --sr --frag=yes -A2 -B4 -O8,16 -E2,1 -r50 -p0.1 -N20000 -f2000,10000 -Y -n2 -m20 -s40 -g200 --end-bonus 12 -2K50m --heap-sort=yes --secondary=yes --sam-hit-only {0} {1}|{EnFlt} {3}|grep -v "^@"|{pigz} -c > {2}'.format(
                    uscg_db, qry, out_sam, max_dist, **executables,
                ), cwd=tmpdir, shell=True).communicate()

        outputs.append(os.path.abspath(os.path.join(tmpdir, out_sam)))
    return outputs



def query_uscg(query, uscg_db, output, metadata, max_dist, pool, min_depth, min_consensus, debug=[False, False, False, False]) :
    if not os.path.isdir(output) :
        os.makedirs(output)

    if not debug[0] :
        sam_files = map_reads(query, uscg_db, output, max_dist)
    else :
        sam_files = [os.path.abspath(os.path.join(output, 'uscg.{0}.sam.gz'.format(id))) for id, _ in enumerate(query)]
    
    uscg_info, genome_ids = get_uscg_info(sam_files[0], metadata)

    if not debug[1] :
        read_maps, file_reads = map_to_uscgs(sam_files, uscg_info, output, max_dist, pool)
        if len(read_maps) == 0 :
            return {'profile':[], 'OTU':[]}
        np.savez_compressed(os.path.join(output, 'uscg.npz'), reads=read_maps, file_reads=file_reads)
    else :
        data = np.load(os.path.join(output, 'uscg.npz'))
        file_reads = data['file_reads']
        if not debug[2] :
            read_maps = data['reads']
    
    if not debug[2] :
        matches, reads = get_matches(metadata, uscg_info, sorted(genome_ids.items(), key=lambda g:g[1]), read_maps)
        if len(matches) == 0 :
            return {'profile':[], 'OTU':[]}
        np.savez_compressed(os.path.join(output, 'out1.npz'), matches=matches, reads=reads)
    else :
        data = np.load(os.path.join(output, 'out1.npz'), allow_pickle=True)
        matches, reads = data['matches'], data['reads']
    if not debug[3] :
        outputs, bam = generate_outputs(sam_files, metadata, uscg_info, sorted(genome_ids.items(), key=lambda g:g[1]), file_reads, matches, reads, output, output)
        json.dump(dict(outputs=outputs, bam=bam), open(os.path.join(output, 'out2.json'), 'wt'))
    else :
        data = json.load(open(os.path.join(output, 'out2.json'), 'rt'))
        outputs, bam = data['outputs'], data['bam']
    write_seq(output, outputs, bam, min_depth, min_consensus)
    return outputs


def write_seq(output, outputs, bam, min_depth=3, min_consensus=0.8) :
    sequences = {}
    prev = 0
    p = subprocess.Popen('{samtools} mpileup -AB -q 0 -Q 0 {0}'.format(bam, **executables).split(), universal_newlines=True, stdout=subprocess.PIPE)
    for line in p.stdout :
        p = line.strip().split('\t')
        p[1] = int(p[1])
        if p[0] not in sequences :
            sequences[p[0]] = []
            prev = p[1]
            dist = 0
        else :
            dist = p[1] - prev - 1
            prev = p[1]
        s = re.sub(r'\^.', '', p[4]).upper().replace('$', '')
        s = list(''.join([b[int(n):] for n, b in re.findall('[+-](\d+)(.+)', '+0' + s)]))
        base, cdp = min(zip(*np.unique(s, return_counts=True)), key=lambda x:-x[1])
        depth = len(s)
        if base in ('*', 'N') :
            base = ''
        elif cdp < min_depth or float(cdp) <= min_consensus * float(depth) :
            base = base.lower()
        sequences[p[0]].append('n'*dist + base)
    seq_out = {}
    for n, s in sequences.items() :
        t = re.split(r'__', n)[1]
        if t not in seq_out :
            seq_out[t] = {}
        seq_out[t][n] = ''.join(s).replace('*', '').replace('+', '').replace('-', '')
    for otu in outputs.get('OTU', []) :
        refs = otu[4]
        otu.append(seq_out.get(refs[0], ''))
    for profile in outputs.get('profile', []) :
        refs = profile[4]
        s = {}
        for ref in refs :
            s.update(seq_out.get(ref, ''))
        profile.append(s)
    json.dump(outputs, open(output + '.json', 'wt'))
    return outputs



@click.command()
@click.option('-q', '--query', help='fastq file [specify multiple times for multiple fastq files]', required=True, multiple=True)
@click.option('-d', '--dbname', help='name of the USCG database [required]', required=True)
@click.option('-o', '--output', help='prefix of the output. [required]', required=True)
@click.option('-t', '--num_threads', help='number of threads [Default: 16]', default=16, type=int)
@click.option('-m', '--max_dist', help='maximum distance of alignment [Default: 0.08]', default=0.08, type=float)
@click.option('--min_depth', help='minimum number of coverage to call a base. [Default: 3]', default=3, type=int)
@click.option('--min_consensus', help='minimum proportion of consensus base to call [Default: 0.8]', default=0.8, type=float)
def main(query, dbname, output, max_dist, num_threads, min_depth, min_consensus) :
    strProfile(query, dbname, output, max_dist, num_threads, min_depth, min_consensus)

def strProfile(query, dbname, output, max_dist, num_threads, min_depth, min_consensus) :
    pool = Pool(max(int(num_threads/4), 1))
    np.random.seed(42)
    query = [os.path.abspath(qry) for qry in query]
    dbname = os.path.abspath(dbname)
    uscg_db = dbname + '.ffn.gz'
    metadata = dbname + '.list'
    if not (os.path.isfile(uscg_db) and os.path.isfile(metadata) ) :
        raise Exception('Cannot find all essential database files.')
        return 1

    if not os.path.isdir(output) :
        os.makedirs(output)

    metadata = {genome:species for genome, species in pd.read_csv(metadata, sep=',', header=None).values}
    outputs = query_uscg(query, uscg_db, output, metadata, max_dist, pool, min_depth, min_consensus)
    json.dump(outputs, open(output + '.json', 'wt'))


if __name__ == '__main__' :
    main()