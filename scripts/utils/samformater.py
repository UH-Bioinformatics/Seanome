


def makeSAMHdr(ref, seqlen, rgroups = None):
    if rgroups:
        return dict(HD = dict(VN = '1.0'), SQ = [ {'SN': ref, 'LN': seqlen }], RG = rgroups )
    else:
        return dict(HD = dict(VN = '1.0'), SQ = [ {'SN': ref, 'LN': seqlen }])


def generateReadGroups(groups):
    RG_template = { 'ID': '',
                    'LB': '',
                    'SM': '',
                    'PL': 'ILLUMINA'} # Platform/technology used to produce the reads.  Spoof it
    readgroups = []
    idtomap = {}
    for idx, g in enumerate(groups):
        rgrp = RG_template.copy()
        rgrp['ID'] = "GROUP-%s"%(idx)
        rgrp['LB'] = g
        rgrp['SM'] = g
        readgroups.append(rgrp)
        idtomap[g] = "GROUP-%s"%(idx)
    return readgroups, idtomap
