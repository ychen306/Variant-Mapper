import re

VARIATIONS = [
            r'(?P<sub>[ATCG]>(?P<sub_base>[ATCG])$)',
            r'(?P<del>del(?:[ATCG]*)$)',
            r'(?P<dup>(?:dup)?(?P<dup_base>[ATCG]+)(?:\[(?P<dup_count>\d+)\])?$)',
            r'(?P<ins>ins(?P<ins_base>[ATCG]+)$)',
            r'(?P<indel>del(?:[ATCG]*)ins(?P<indel_base>[ATCG]+)$)'
        ]
variation_re = re.compile('|'.join(VARIATIONS))
cigar_re = re.compile(r'(\d+)([MIDNSHPX=])')


def get_complements(bases):
    '''
    Get complements of a DNA sequence
    '''
    comp_bases = ''
    complements = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
    for base in bases:
        comp_bases += complements[base]
    return comp_bases


def cigar_ops(cigar):
    '''
    Yield cigar operations. E.g. yield M, M, M, I, I, D for '3M2I1D'
    '''
    cigar_groups = ((int(cigar_group.group(1)), cigar_group.group(2))
                    for cigar_group in cigar_re.finditer(cigar))

    for count, op in cigar_groups:
        for _ in xrange(count):
            yield op


def is_reverse(read):
    '''
    check whether a read is from the reverse (-) strand
    '''
    return (read['flags']>>4) % 2 == 1


def get_ref_length(cigar):
    '''
    Get length of portion of reference sequence that maps to the read.
    '''
    length = 0
    for cigar_group in cigar_re.finditer(cigar):
        count, op = int(cigar_group.group(1)), cigar_group.group(2)
        if op in ('M', 'D', 'N', 'S', 'H', '=', 'X'):
            length += count
    return length

def get_bases_from_read(read, start, end):
    '''
    Get aligned bases from reads; this is different from read.alignedBases
    in that this takes insertion into account.
    start and end indicate the region of interest.
    '''
    #operations that increment the reference sequence
    inc_ops = ('M', 'D', 'N', 'S', 'H', '=', 'X')
    #operations that "occupy" the SEQ
    oc_ops = ('M', 'I', 'S', '=', 'X')
    position = read['position']
    bases = iter(read['originalBases'])
    cigar = read['cigar']
    #base of reference sequence
    ref_base = 0
    end_base = end - position
    aligned_bases = ''
    operations = cigar_ops(cigar)
    while ref_base <= end_base:
        current_op = operations.next()
        if current_op in oc_ops:
            #current base of SEQ
            current_base = bases.next()
            if ref_base >= start-position:
                aligned_bases += current_base
        if current_op in inc_ops:
            ref_base += 1
    return aligned_bases




def get_bases_from_hgvs(hgvs):
    '''
    Translate HGVS into bases
    '''
    variation = variation_re.match(hgvs)
    if not variation:
        raise ValueError('Unable to interpret HGVS notation: %s' % hgvs)
    var_type = variation.lastgroup
    if var_type in ('sub', 'ins', 'indel'):
        bases = variation.group('%s_base' % var_type)
    elif var_type == 'dup':
        dup_count = variation.group('dup_count')
        if not dup_count:
            dup_count = 2
        else:
            dup_count = int(dup_count)
        bases = variation.group('dup_base') * dup_count
    else:
        #means it's deletion
        bases = ''
    return bases

def match(report, read):
    '''
    Check if a variation in a report mathces a read.
    '''
    try: 
        read_bases = get_bases_from_read(read, report['seqStart'], report['seqEnd'])
        hgvs_bases = get_bases_from_hgvs(report['variation'])
    except (StopIteration, ValueError):
        return False
    reverse_report = (report['strand'] == '-')
    reverse_read = is_reverse(read)
    if reverse_report != reverse_read:
        return read_bases == get_complements(hgvs_bases)
    else:
        return read_bases == hgvs_bases