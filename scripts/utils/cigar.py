import re


cigar_re = re.compile(r'([0-9]*)([M=XID])')
FINDALL_CMP=re.compile(r'([0-9]*)([MID=])')

def expandCigar( cigar):
    cigstr = ""
    for m in cigar_re.finditer(cigar):
        if m.group(1):
            cigstr +=  m.group(2) * int(m.group(1))
        else:
            cigstr +=  m.group(2)
    return cigstr


def compressCigar(cigar):
    """ 
    We build the cigar as a string of I,M, and D.  We need to do a RLE on this 
    to compact it. 
    """

    cigar = cigar.rstrip("D")
    c = None
    cnt = 0
    ciglst = []
    for b in cigar:
        if b != c:
            if c!= None:
                ciglst.append("%s%s"%(cnt, c))
            cnt = 1
            c = b
        else:
            cnt += 1
    if c!= None:
        ciglst.append("%s%s"%(cnt, c))
    return "".join(ciglst)



# I='D' or I='N'        
CIGREMAP = dict(M='M', D='I', I='D')
def interpretAln(parts):
    index = 0
    rm = 0
    pos = 0
    for idx, f in enumerate(parts):
        index = idx
        if f[1] == 'I':
            pos += int(f[0])
        elif f[1] == 'D':
            rm += int(f[0])
        else:
            break
    return index, rm, pos


def cleanupCigar(pos, cig, seqlen):
    parts = FINDALL_CMP.findall(cig)
    # in all its wisdom.. usearch can return = for 
    # the alignment, if we have a perfect match.
    if len(parts) == 1:
        if parts[0][1] == '=':
            return "%sM"%(seqlen)
        else:
            return cig, pos, 0,0 
    # fix the RLE so that everything has both elements
    parts = [ (p[0] if p[0] else '1', p[1],) for p in parts]
    f_idx, rmfront, shiftright = interpretAln(parts)
    # trim and reverse the parts so that we pull off the tail end for processing and dont overlap the forward
    r_idx, rmback, noshift = interpretAln(parts[f_idx:][::-1]) 
    # r_idx == 0 would make the final result empty
    if r_idx:
        parts = parts[f_idx: -r_idx]
    else:
        parts = parts[f_idx:]
    newcig = ''.join( ( "".join([p0, CIGREMAP[p1]]) for p0,p1 in parts) )
    return newcig, pos + shiftright, rmfront, rmback
