import sys
import argparse
import time
import zipfile


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need

"""

def SuffixArray(text):
    suffix_to_pos = dict()
    for i in range(len(text)):
        suffix_to_pos[text[i:]] = i
    return suffix_to_pos

def BurrowsWheelerTransform(text):
    matrix = []
    for i in range(len(text)):
        temp = text[len(text) - i:] + text[:len(text) - i]
        matrix.append(temp)

    matrix.sort()
    res = ''
    for e in matrix:
        res += e[len(text) - 1]
    return res

def Count(lastcol):
    res = dict()
    for i in range(len(lastcol)):
        if lastcol[i] not in res:
            res[lastcol[i]] = [0]
    for i in range(len(lastcol)):
        for k in res.keys():
            if k == lastcol[i]:
                res[k].append(res[k][-1] + 1)
            else:
                res[k].append(res[k][-1])
    return res

def FirstOccurrence(lastcol):
    firstcol = [c for c in lastcol]
    firstcol.sort()
    res = dict()
    for i in range(len(firstcol)):
        if firstcol[i] not in res:
            res[firstcol[i]] = i
    return res

# returns number of matches of pattern in text, given only BW transform
def BetterBWMatching(lastcol, pattern):
    count = Count(lastcol)
    fo = FirstOccurrence(lastcol)
    lastcol = [c for c in lastcol]
    top = 0
    bottom = len(lastcol) - 1
    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if(symbol in lastcol[top:bottom+1]):
                top = fo[symbol] + count[symbol][top]
                bottom = fo[symbol] + count[symbol][bottom + 1] - 1
            else:
                return 0
        else:
            return bottom - top + 1

def HammingDistance(s1, s2):
    if len(s1) != len(s2):
        return math.inf
    diff = []
    num_diff = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            num_diff += 1
    return num_diff

def divide(pattern, d):
    length = len(pattern)
    l = length // (d+1)
    k = length % (d+1)
    result = []
    i = 0
    while i < length:
        if k > 0:
            result.append((pattern[i:i+l+1], i))
            k -= 1
            i += l+1
        else:
            result.append((pattern[i:i+l], i))
            i += l
    return result

def APM(text, p, d, lastcol, sa):
    res = []
    l = len(p)
    splits = divide(p, d)
    perfectMatch = False
    for s, i in splits:
        if perfectMatch:
            break
        rest = p[:i] + p[i+len(s):]

        # there is at least one match for the pattern
        if BetterBWMatching(lastcol, s) > 0:

            # retrieve the positions of each match in the reference
            matches = [val for key, val in sa.items() if key.startswith(s)]

            # check if there are mutations
            for m in matches:
                target = text[m-i:m] + text[m+len(s):m-i+l]
                hd = HammingDistance(rest, target)
                if hd == 0:
                    perfectMatch = True
                    if (m-i) not in res:
                        res.append(m-i)
                # should this be else???
                elif hd <= d:
                    if (m-i) not in res:
                        res.append(m-i)
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')

    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    """
        TODO: Call functions to do the actual read alignment here
        
    """
    # split paired-end reads into single reads
    single_reads = []
    for read in input_reads:
        h1 = read[0]
        # must reverse second half of paired-end read
        h2 = read[1][::-1]
        single_reads.append(h1)
        single_reads.append(h2)

    print('Finished splitting paired-end reads')

    # construct Burrows-Wheeler transform for reference
    bw_text = reference + '$'
    bwt = BurrowsWheelerTransform(bw_text)

    print('Finished constructing BWT')

    # construct suffix array for reference to obtain positions of matches
    suffix_arr = SuffixArray(reference)

    print('Finished constructing suffix array')

    # perform BWT pattern matching with at most one mismatch in read
    matches = []
    for r in single_reads:
        positions = APM(reference, r, 1, bwt, suffix_arr)
        matches.append([r, positions])

    print('Finished BWT pattern matching')

    # find mismatches in reads matched with reference positions
    mismatches = dict()
    for m in matches:
        read = m[0]
        positions = m[1]
        l = len(read)
        for p in positions:
            ref_segment = reference[p:p+l]
            for i in range(l):
                if read[i] != ref_segment[i]:
                    if read in mismatches:
                        mismatches[read].append((ref_segment[i], read[i], p+i))
                    else:
                        mismatches[read] = [(ref_segment[i], read[i], p+i)]

    print('Finished finding mismatch positions')

    # only keep mutations that appear multiple times (to avoid sequencing errors)
    mutations = dict()
    seen_mismatches = []
    for read, muts in mismatches.items():
        for mut in muts:
            if mut in seen_mismatches:
                if mut[2] in mutations:
                    mutations[mut[2]].append((mut[0], mut[1]))
                else:
                    mutations[mut[2]] = [(mut[0], mut[1])]
            else:
                seen_mismatches.append(mut)

    snps = []

    for pos, muts in mutations.items():
        if len(muts) < 2:
            continue
        else:
            # use consensus algorithm to find correct mutation
            diffs = 0
            same = True
            for m in muts:
                if m != muts[0]:
                    same = False
                    diffs += 1

            if same:
                snps.append([muts[0][0], muts[0][1], pos])
            else:
                # tie between two different mutations, don't include as SNP
                if len(muts) - diffs <= 1:
                    continue
                # one mutation has majority vote, final majority mutation and include as SNP
                else:
                    occurrences = dict()
                    for m in muts:
                        if m in occurrences:
                            occurrences[m] += 1
                        else:
                            occurrences[m] = 1
                    max = 0
                    final_mut = None
                    for m in list(occurrences.keys()):
                        if occurrences[m] > max:
                            final_mut = m
                            max = occurrences[m]
                    snps.append([final_mut[0], final_mut[1], pos])


    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)

