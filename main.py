import random
from graph import Graph
import re

K_MER_SIZE: int = 15
FILE = 'example.fastq'

class Correction_ratio:
    def __init__(self, fn=0,fp=0,tp=0):
        self.fn = fn
        self.fp = fp
        self.tp = tp

    def add_false_negative(self):
        return Correction_ratio(fn=self.fn+1,fp=self.fp,tp=self.tp)
    def add_false_positive(self):
        return Correction_ratio(fn=self.fn,fp=self.fp+1,tp=self.tp)
    def add_true_positive(self):
        return Correction_ratio(fn=self.fn,fp=self.fp,tp=self.tp+1)

class Code:
    def __init__(self, left, corrected, right):
        self.left = left
        self.corrected = corrected
        self.right = right
    def correct_right(self, correction: str):
        return Code(self.left, self.corrected+correction,self.right[1:])
    def correct_left(self, correction: str):
        return Code(self.left[:-1],correction+self.corrected,self.right)
    def right_is_empty(self):
        return not self.right
    def left_is_emtpy(self):
        return not self.left
    def is_corrected(self):
        return self.left_is_emtpy() and self.right_is_empty()


def get_random_base():
    return random.choice('AGCT')



def get_genetic_string(repeats: int):
    result = ''
    for x in range(repeats):
        result = result + get_random_base()
    return result

def introduce_errors(read: str, chance: float):
    result = ''
    for s in read:
        if random.random() <= chance:
            result = result + get_random_base()
        else:
            result = result + s
    return result

def get_k_mers(read: str):
    k_mers = list()
    for i in range(len(read) - K_MER_SIZE + 1):
        k_mers.append(read[i:i + K_MER_SIZE])
    return k_mers

def get_graph_from_k_mers(k_mers: list, g: Graph):
    for k_mer in k_mers:
        g.add_node(k_mer)
        prefix = k_mer[:-1]
        predecessors = ['A'+prefix, 'G'+prefix, 'C'+prefix, 'T'+prefix]
        suffix = k_mer[1:]
        successors = [suffix+'A',suffix+'G',suffix+'C',suffix+'T']
        for pred in predecessors:
            if g.has_node(pred):
                g.add_edge((pred, k_mer))
        for succ in successors:
            if g.has_node(succ):
                g.add_edge((k_mer,succ))
    return g

def count_k_mers(k_mers, d: dict):
    for k_mer in k_mers:
        if k_mer in d:
            d[k_mer] = d[k_mer] + 1
        else:
            d[k_mer] = 1
    return d

def k_mers_count_above_solidity_threshold(dictionary: dict, threshold: int):
    result = dict()
    for d in dictionary.items():
        if d[1] > threshold:
            result[d[0]] = d[1]
    return result

def k_mers_above_solidity_threshold(dictionary: dict, threshold: int):
    result = list()
    for d in dictionary.keys():
        if dictionary[d] > threshold:
            result.append(d)
    return result

def compact_graph(g: Graph):
    paths = g.get_simple_paths()
    l = len(paths)
    c = 0
    for path in paths:
        unitig = ''
        for node in path:
            if node == path[-1]:
                unitig = unitig + node
            else:
                unitig = unitig + node[0]
        first, last = path[0], path[-1]
        preds, succs = g.get_neighbours(first,incoming=True), g.get_neighbours(last, incoming=False)
        g.add_node(unitig)
        for p in preds:
            g.add_edge((p,unitig))
        for s in succs:
            g.add_edge((unitig, s))
        for node in path:
            g.delete_node(node)

def remove_dead_ends(g: Graph, limit: int):
    def is_dead_end(n):
        zero_in = g.count_edges(n, incoming=True) == 0
        zero_out = g.count_edges(n, incoming=False) == 0
        one_in = g.count_edges(n, incoming=True) == 1
        one_out = g.count_edges(n, incoming=False) == 1
        return (zero_in and one_out) or (zero_out and one_in)
    dead_ends = [n for n in g.nodes if is_dead_end(n)]
    short_dead_ends = [n for n in dead_ends if len(n) < limit]
    for n in short_dead_ends:
        g.delete_node(n)
    return g

def remove_all_below_unitig_threshold(g: Graph, threshold: float, k_mer_count: dict):
    for n in g.nodes:
        k_mers: list = get_k_mers(n)
        unitig_abundance = sum([k_mer_count[k_mer] for k_mer in k_mers])/len(k_mers)
        if unitig_abundance < threshold:
            g.delete_node(n)
    return g

def clean_graph(g: Graph, iterations: int, limit: int, threshold: float, k_mer_count: dict):
    for i in range(iterations):
        remove_dead_ends(g, limit)
        remove_all_below_unitig_threshold(g, threshold, k_mer_count)

def read_mapping(g: Graph, read: str, hamming_distance: int):
    result = list()
    def extend_right(node: str, unitig_left, errors_left, corrected_code: Code):
        no_errors_left = (errors_left == 0)
        read_is_finished = corrected_code.right_is_empty()
        need_next_unitig = (not unitig_left)
        read_left = corrected_code.right
        if no_errors_left:
            result.append(False)
        elif read_is_finished:
            result.append((errors_left, corrected_code))
        elif need_next_unitig:
            neighbours = g.get_neighbours(node, incoming=False)
            for n in neighbours:
                extend_right(n, n, errors_left, corrected_code)
        else:
            first_letter_matches = (unitig_left[0] == read_left[0])
            if first_letter_matches:
                code = corrected_code.correct_right(unitig_left[0]) #change nothing + move right
                extend_right(node, unitig_left[1:], errors_left, code)
            else:
                correction = unitig_left[0]
                code = corrected_code.correct_right(correction)
                extend_right(node, unitig_left[1:], errors_left-1, code)

    def extend_left(seed: str, node: str, unitig_left, errors_left, corrected_code: Code):
        no_errors_left = (errors_left == 0)
        read_is_finished = corrected_code.left_is_emtpy()
        need_next_unitig = (not unitig_left)
        read_left = corrected_code.left
        if no_errors_left:
            result.append(False)
        elif read_is_finished:
            seed_node = [n for n in g.nodes if (seed in n)][0]
            extend_right(seed_node, seed_node.split(seed)[1],errors_left, corrected_code)
        elif need_next_unitig:
            neighbours = g.get_neighbours(node, incoming=True)
            for n in neighbours:
                extend_left(seed, n, n, errors_left, corrected_code)
        else:
            last_letter_matches = (unitig_left[-1] == read_left[-1])
            if last_letter_matches:
                code = corrected_code.correct_left(unitig_left[-1])  # change nothing + move left
                extend_left(seed, node, unitig_left[:-1], errors_left, code)
            else:
                correction = unitig_left[-1]
                code = corrected_code.correct_left(correction)
                extend_left(seed, node, unitig_left[:-1], errors_left - 1, code)

    for i in range(len(read)-K_MER_SIZE+1):
        seed = read[i:i+K_MER_SIZE]
        left = read[:i]
        right = read[i + K_MER_SIZE:]
        code = Code(left, seed, right)
        nodes = [n for n in g.nodes if (seed in n)]
        if nodes:
            node = nodes[0]
            extend_left(seed,node,node.split(seed)[0],hamming_distance,code)
        else:
            result.append(False)
    return result

def read_file():
    file = open(FILE,'r')
    counts = dict()
    while True:
        line = file.readline()
        if not line:
            break
        else:
            x = re.match("[ACGT]+",line)
            if x:
                k_mers = get_k_mers(line)
                count_k_mers(k_mers, counts)
    return counts

def divide_into_reads(code: str):
    l = len(code)
    reads = set()
    for i in range(200):
        read_len = random.randrange(20,100)
        start = random.randrange(0,l-read_len)
        real_read = code[start:start + read_len]
        simulated_read = introduce_errors(real_read, 0.05)
        print(simulated_read)
        reads.add(simulated_read)
    print(reads)
    return reads


real_code = get_genetic_string(4000)
print('real sequence:  ', real_code)
print('READS:')
reads = divide_into_reads(real_code)
print('K-MERS:')
count = dict()
for read in reads:
    k_mers = get_k_mers(read)
    count_k_mers(k_mers,count)
    print(k_mers)
print(count)
filtered = k_mers_above_solidity_threshold(count, 0)
print(filtered)
g = get_graph_from_k_mers(filtered, Graph())
compact_graph(g)
clean_graph(g, 5, (K_MER_SIZE-1),4,count)
print(g.nodes)
print(g.edges)
for r in reads:
    res = read_mapping(g, r, 5)
    corrected = [x for x in res if x]
    print(corrected)
    if corrected:
        for c in corrected:
            print(str(c[0]) + " remaining: ", c[1].corrected)


# genetic_read = get_genetic_string(16000)
# errors_read = introduce_errors(genetic_read, 0.1)
# k_mers = get_k_mers(errors_read)
# k_mer_count = count_k_mers(k_mers, dict())
# filtered_k_mer_count = k_mers_count_above_solidity_threshold(k_mer_count, 2)
# filtered_k_mers = k_mers_above_solidity_threshold(k_mer_count, 2)
# g = get_graph_from_k_mers(filtered_k_mers, Graph())
# compact_graph(g)
# clean_graph(g,10, 2*(K_MER_SIZE-1),4,k_mer_count)
# res = read_mapping(g, errors_read, 10)
# print('real sequence:  ', genetic_read)
# print('erroneous read: ',errors_read)
# print('k-mers:         ',k_mers)
# print('counts:         ',k_mer_count)
# print("after threshold:",filtered_k_mers)
# print('nodes:          ',g.nodes)
# print('edges:          ',g.edges)
# print(res)
#
#
# count = read_file()
# print(count)
# k_mers = k_mers_above_solidity_threshold(count, 15)
# k_mer_count = k_mers_count_above_solidity_threshold(count, 15)
# print(k_mer_count)
# print(len(k_mers))
# g = get_graph_from_k_mers(k_mers, Graph())
# compact_graph(g)
# clean_graph(g, 5, 2*(K_MER_SIZE-1), 40, count)
# file = open(FILE,'r')
# while True:
#     line = file.readline()
#     if not line:
#         break
#     else:
#         x = re.match("[ACGT]+", line)
#         if x:
#             x = read_mapping(g,line,10)
#             print(x)