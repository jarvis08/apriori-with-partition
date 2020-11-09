import sys
from collections import Counter
from itertools import combinations


def load_xacts(input_filename):
    xacts = []
    with open(input_filename, "r") as f:
        read = f.readlines()
        for line in read:
            splitted = line.split("\t")
            splitted[-1] = splitted[-1].replace("\n", '')
            xacts.append(list(map(int, splitted)))
    return xacts


def make_candidates_of_candidates(prev_fq, len_c):
    if len_c == 2:
        return [set(x) for x in combinations(prev_fq, len_c)]
    else:
        cur_items = []
        for item_set in prev_fq:
            for item in item_set:
                cur_items.append(item)
        items = set(cur_items)
        return [set(x) for x in combinations(items, len_c)]


def make_candidates(prev_fq, cc, len_c):
    candidates = {}
    if len_c == 2:
        for item_set in cc:
            for item in [x for x in combinations(item_set, 1)]:
                if item[0] not in prev_fq:
                    break
            else:
                candidates[tuple(item_set)] = 0
        return candidates

    else:
        prev = [set(x) for x in prev_fq]
        for item_set in cc:
            for subset in [set(x) for x in combinations(item_set, len_c - 1)]:
                if subset not in prev:
                    break
            else:
                # tuple can be key for dict
                candidates[tuple(item_set)] = 0
        return candidates


def check_sup(candidates, min_cnt):
    frequent_set = {}
    for item_set, cnt in candidates.items():
        if cnt >= min_cnt:
            frequent_set[item_set] = cnt
    return frequent_set


def save_association_rules(tdb, len_tdb, frequent_set, len_fq, output_filename):
    with open(output_filename, 'w') as f:
        for item_set, cnt in frequent_set.items():
            cur_len = len_fq
            if not len_fq:
                cur_len = len(item_set)
            while cur_len > 1:
                for item in list(combinations(item_set, cur_len-1)):
                    item = set(item)
                    diff = set(item_set).difference(item)
                    support = cnt / len_tdb * 100

                    condition = 0
                    for xact in tdb:
                        if set(xact) >= item:
                            condition += 1
                    confidence = cnt / condition * 100

                    rule = str(item) + '\t' + str(diff) + '\t' + str('%.2f' % round(support, 2)) + '\t' + str('%.2f' % round(confidence, 2)) + '\n'
                    f.write(rule)
                cur_len = cur_len - 1
    

def partitioned_apriori(min_sup, input, num_partition, output):
    tdb = load_xacts(input)
    n_xacts = len(tdb)
    min_cnt = n_xacts * min_sup
    if num_partition < 2:
        num_partition = 1

    global_frequent_patterns = {}

    n_partition = num_partition
    n_local_xacts = int(n_xacts / n_partition)
    local_min_sup = min_sup / n_partition
    for i in range(n_partition):
        partitioned = []
        if i == n_partition - 1:
            partitioned = tdb[i * n_local_xacts:]
            n_local_xacts = len(partitioned)
        else:
            partitioned = tdb[i * n_local_xacts : (i+1) * n_local_xacts]

        local_min_cnt = local_min_sup * n_local_xacts
        frequent_patterns = []

        # make frequent 1-itemset
        all_items = {}
        for xact in partitioned:
            for item in xact:
                if item not in all_items.keys():
                    all_items[item] = 0
                all_items[item] += 1

        frequent_1 = check_sup(all_items, local_min_cnt)
        frequent_patterns.append(frequent_1)

        while True:
            l = len(frequent_patterns) + 1
            candidates_of_candidates = make_candidates_of_candidates(frequent_patterns[-1], l)
            candidates = make_candidates(frequent_patterns[-1], candidates_of_candidates, l)

            for xact in partitioned:
                for c in candidates.keys():
                    if set(c) <= set(xact): 
                       candidates[c] += 1 
            frequent_l = check_sup(candidates, local_min_cnt)

            if not len(frequent_l):
                break
            else:
                frequent_patterns.append(frequent_l)
                tuple_key = list(global_frequent_patterns.keys())
                set_key = [set(k) for k in tuple_key]
                for pattern, cnt in frequent_l.items():
                    tmp = set(pattern)
                    if tmp in set_key:
                        key_idx = set_key.index(tmp)
                        global_frequent_patterns[tuple_key[key_idx]] += cnt
                    else:
                        global_frequent_patterns[tuple(pattern)] = cnt
    sup_checked = check_sup(global_frequent_patterns, min_cnt)
    save_association_rules(tdb, n_xacts, sup_checked, 0, output)
    

if __name__ == '__main__':
    argv = sys.argv
    MIN_SUPPORT = 0.05
    INPUT_FILENAME = "input.txt"
    OUTPUT_FILENAME = "output.txt"
    N_PARTITION = 4
    if len(argv) < 4:
        print("Not enough arguments given.")
        print("Run apriori algorithm with default paramters")
        print("Minimum support : 5 %")
        print("Input filename : input.txt")
        print("Output filename : output.txt")
    else:
        MIN_SUPPORT = float(argv[1]) / 100
        INPUT_FILENAME = argv[2]
        OUTPUT_FILENAME = argv[3]
        print("Run apriori algorithm with given paramters")
        print("Minimum support : {}".format(MIN_SUPPORT))
        print("Input filename : {}".format(INPUT_FILENAME))
        print("Output filename : {}".format(OUTPUT_FILENAME))
    partitioned_apriori(MIN_SUPPORT, INPUT_FILENAME, N_PARTITION, OUTPUT_FILENAME)
