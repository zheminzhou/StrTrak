import ete3, re, numpy as np, numba as nb, pandas as pd
import _collections

@nb.njit(fastmath=True)
def __pd(c1, c2, res) :
    for i1 in range(c1.shape[0]) :
        n1 = c1[i1]
        for i2 in range(c2.shape[0]) :
            n2 = c2[i2]
            d = n1[1] + n2[1]
            res[int(n1[0]), int(n2[0])] = d
            res[int(n2[0]), int(n1[0])] = d
    return res
    

def pairwise_distance(tre) :
    names = { leaf:id for id, leaf in enumerate(tre.get_leaf_names()) }
    distance = np.zeros([len(names), len(names)], dtype=float)
    
    for node in tre.traverse('postorder') :
        if node.is_leaf() :
            node.__d = np.array([[names[node.name], node.dist]], dtype=float)
        else :
            desc = [c.__d for c in node.children]
            for i2, c2 in enumerate(desc) :
                for c1 in desc[:i2] :
                    __pd(c1, c2, distance)
            
            node.__d = np.vstack(desc)
            node.__d.T[1] += node.dist
    name_list = [n for n, i in sorted(names.items(), key=lambda n:n[1])]
    return pd.DataFrame(distance, columns=name_list)
    

def prune(tre, kept_tips) :
    for node in tre.iter_descendants('postorder') :
        if node.is_leaf() :
          if node.name not in kept_tips :
                up = node.up
                up.remove_child(node)
                node.up = None
        else :
            if len(node.children) < 2 :
                up = node.up
                dist = node.dist
                for c in node.children :
                    c.dist += dist
                    up.add_child(c)
                    c.up = up
                up.remove_child(node)
                node.up = None
    while len(tre.children) < 2 :
        tre = tre.children[0]
    return tre


def write_nexus(trees) :
    names = _collections.OrderedDict([[n, i + 1] for i, n in enumerate(sorted(trees[0].get_leaf_names()))])

    res = '''#NEXUS

Begin taxa;
        Dimensions ntax={0};
        Taxlabels
                {1}
                ;
End;

Begin trees;
        Translate
                {2}
                ;
'''.format(len(names),
           '\n                '.join(names.keys()),
           '\n                '.join(['{0} {1},'.format(v, k) for k, v in names.items()]), )

    for tid, tre in enumerate(trees):
        try :
            tre = tre.copy()
        except :
            pass
        for n in tre.get_leaves():
            n.name = str(names[n.name])

        features = {}
        for nid, node in enumerate(tre.traverse('preorder')):
            ann = []
            if 'annotations' in node.__dict__ :
                for k, v in sorted(node.annotations.items()):
                    if isinstance(v, (list,tuple,dict,set)):
                        if len(v) > 0 and isinstance(v[0], str):
                            v = '{{{0}}}'.format(','.join(['"{0}"'.format(vv) for vv in v]))
                        else:
                            v = '{{{0}}}'.format(','.join(['{0}'.format(vv) for vv in v]))
                    elif isinstance(v, str):
                        v = '"{0}"'.format(v)
                    else:
                        v = str(v)
                    ann.append('{0}={1}'.format(k, v))
            node.name += '___{0}___'.format(nid)
            features['___{0}___'.format(nid)] = '[&{0}]'.format(','.join(ann))

        tree_str = 'tree TREE{0} = [&R] '.format(tid + 1) + re.sub(r'(___\d+___)', lambda g: features[g.group(1)],
                                                                   tre.write(format=1)[:-1]+tre.name+';') + '\n'
        res += tree_str
    res += 'End;\n'
    return res


def read_nexus(tre_file) :
    tree_list = []

    inTranslate = False
    translate_table = {}

    annotation_tags = {}
    def get_tag(g) :
        ann = g.group(0)
        key = '___{0}___'.format(len(annotation_tags))
        annotation_tags[key] = ann
        return key

    def repos(g) :
        return '{0}:{1}'.format(g.group(2), g.group(1))

    with open(tre_file) as fin:
        for line in fin :
            line = line.strip()
            if line.lower() == 'translate' :
                inTranslate = True
            elif inTranslate :
                p = line.split()
                if p[0] == ';' :
                    inTranslate = False
                else :
                    translate_table[p[0]] = p[1].strip(',')
            elif line.lower().strip().startswith('tree') :
                prefix, tree_line = line.split('(', 1)
                try :
                    tre_name = re.findall('tree (\S+) *=', prefix.strip())[0]
                except :
                    tre_name = ''
                tree_line = '(' + tree_line
                annotation_tags = {}
                tree_line = re.sub('\[[^\]]+\]',get_tag, tree_line)
                tree_line = re.sub(':([\d\+\.eE-]+)(___\d+___)', repos, tree_line)
                tre = ete3.Tree(tree_line, format=1)
                for node in tre.traverse() :
                    names = re.findall('^(.*)(___\d+___)$', node.name)
                    node.annotations = {}
                    if len(names) :
                        name, ann_key = names[0]
                        node.name = name
                        annotation = annotation_tags[ann_key]
                        x = annotation[1:-1].split('=')
                        annotations = [an.strip('&') for ann in x[:-1] for an in ann.rsplit(',', 1)] + [x[-1]]
                        for k, v in zip(annotations[::2], annotations[1::2]) :
                            if v == '{}' :
                                v = []
                            elif v.startswith('{') :
                                v = [vv.strip('"') if vv.startswith('"') else float(vv) for vv in v[1:-1].split(',')]
                            else :
                                v = v.strip('"') if v.startswith('"') else float(v)
                            node.annotations[k] = v
                    if node.is_leaf() and node.name in translate_table :
                        node.name = translate_table[node.name]
                tre.tree_name = tre_name
                tree_list.append(tre)
    return tree_list


def read_trees(tre_file, format=0) :
    trees = []
    with open(tre_file) as fin :
        for line in fin :
            if line.lower().startswith('#nexus') :
                break
            tre = ete3.Tree(line.strip(), format=format)
            trees.append(tre)
    if len(trees) == 0 :
        trees = iter_nexus(tre_file)
    return trees


def iter_trees(tre_file, format=0) :
    with open(tre_file) as fin :
        for line in fin :
            if line.lower().startswith('#nexus') :
                break
            tre = ete3.Tree(line.strip(), format=format)
            yield tre
    yield iter_nexus(tre_file)


def iter_nexus(tre_file) :
    inTranslate = False
    translate_table = {}

    annotation_tags = {}
    def get_tag(g) :
        ann = g.group(0)
        key = '___{0}___'.format(len(annotation_tags))
        annotation_tags[key] = ann
        return key

    def repos(g) :
        return '{0}:{1}'.format(g.group(2), g.group(1))

    with open(tre_file) as fin:
        for line in fin :
            line = line.strip()
            if line.lower() == 'translate' :
                inTranslate = True
            elif inTranslate :
                p = line.split()
                if p[0] == ';' :
                    inTranslate = False
                else :
                    translate_table[p[0]] = p[1].strip(',')
            elif line.lower().strip().startswith('tree') :
                prefix, tree_line = line.split('(', 1)
                try :
                    tre_name = re.findall('tree (\S+) *=', prefix.strip())[0]
                except :
                    tre_name = ''
                tree_line = '(' + tree_line
                annotation_tags = {}
                tree_line = re.sub('\[[^\]]+\]',get_tag, tree_line)
                tree_line = re.sub(':([\d\+\.eE-]+)(___\d+___)', repos, tree_line)
                tre = ete3.Tree(tree_line, format=1)
                for node in tre.traverse() :
                    names = re.findall('^(.*)(___\d+___)$', node.name)
                    node.annotations = {}
                    if len(names) :
                        name, ann_key = names[0]
                        node.name = name
                        annotation = annotation_tags[ann_key]
                        annotations = [ an.strip('&') for ann in annotation[1:-1].split('=') for an in ann.rsplit(',', 1)]
                        for k, v in zip(annotations[::2], annotations[1::2]) :
                            if v.startswith('{') :
                                v = v[1:-1].split(',')
                            node.annotations[k] = v
                    if node.is_leaf() and node.name in translate_table :
                        node.name = translate_table[node.name]
                tre.tree_name = tre_name
                yield tre
