import numpy
from ssqueezepy import cwt
from Bio.Seq import Seq
from ssqueezepy import Wavelet

wavlet = Wavelet(('gmw', {'gamma': 3, 'beta': 60, 'norm': 'bandpass', 'centered_scale': True, 'dtype': 'float32'}))


def get_numseq(seq, is_dna):
    seq = seq.strip().upper()
    seq.strip()
    if is_dna:
        code = {'A': 0.1260, 'C': 0.1340, 'G': 0.0806, 'T': 0.1335}
    else:
        code = {'D': 2.77, 'A': 6.00, 'R': 10.76, 'N': 5.41, 'C': 5.05, 'Q': 5.65, 'E': 3.22, 'G': 5.97,
                'H': 7.59, 'I': 6.02, 'L': 5.98, 'K': 9.74, 'M': 5.74, 'F': 5.48, 'P': 6.30, 'S': 5.68, 'T': 6.16,
                'W': 5.89, 'Y': 5.66, 'V': 5.96, 'B': 0}
        s_len = len(seq)
        if s_len % 3 != 0:
            seq = seq[0:(s_len - s_len % 3)]
        seq = Seq(seq).translate(stop_symbol='B')

    seq_list = numpy.array(seq)
    res = numpy.zeros(len(seq_list))
    for c_i, c_v in code.items():
        res[seq_list == c_i] = c_v
    return res


def findnzeros(num):
    res = []
    p_num = num.copy()
    p_num[p_num > 0] = 1
    num_h = numpy.append(p_num, 0)
    num_t = numpy.append(0, p_num)
    num_index = num_h - num_t

    seg_count = sum(abs(num_index)) / 2
    seg_start = numpy.where(num_index == 1)[0]
    seg_end = numpy.where(num_index == -1)[0]
    seg_len = seg_end - seg_start
    seg_mid = 0.5 * (seg_start + seg_end) / seg_end
    orseg_len = sorted(seg_len)
    if seg_count != 1:
        split_p = numpy.mean(orseg_len[0:round(seg_count / 3)])
        seg_3u_lenav = numpy.mean(seg_len[orseg_len >= split_p])
        segs = zip(seg_start[seg_len >= split_p], seg_end[seg_len >= split_p])
        sum_3u_pw = 0
        segs_num = sum(orseg_len >= split_p)
        for s_i, e_i in segs:
            sum_3u_pw = sum_3u_pw + sum(num[s_i:e_i])
        seg_3u_pav = sum_3u_pw / sum(seg_len >= split_p)
    else:
        seg_3u_lenav = seg_len[0]
        seg_3u_pav = sum(num[seg_start[0]:seg_end[0]])
        segs_num = 1

    res.append(seg_count)
    res.append(numpy.mean(seg_len))
    res.append(numpy.std(seg_len))
    res.append(numpy.mean(seg_mid))
    res.append(numpy.std(seg_mid))
    res.append(seg_3u_lenav)
    res.append(seg_3u_pav)
    res.append(segs_num)
    res.append(sum(seg_len >= 30))
    seg30 = seg_mid[seg_len >= 30]
    if sum(seg30) != 0:
        res.append(numpy.mean(seg_mid[seg_len >= 30]))
        res.append(numpy.std(seg_mid[seg_len >= 30]))
    else:
        res.append(0)
        res.append(0)
    return res


def get_Ewtf_back(numseq, wavlet):
    res = []
    ns = numseq
    for i in [0, 1]:
        w, s = cwt(ns, wavlet, cache_wavelet=True)
        csum_w = sum(abs(w) ** 2)
        ns = csum_w.copy()
        pw_av = numpy.mean(csum_w)
        csum_w[csum_w < pw_av] = 0
        t = findnzeros(csum_w)
        res = res + t
    return res



def get_Ewtf(numseq, wavlet):
    res = []
    w, s = cwt(numseq, wavlet, cache_wavelet=True)
    csum_w = sum(abs(w) ** 2)
    pw_av = numpy.mean(csum_w)
    csum_w[csum_w < pw_av] = 0
    res = findnzeros(csum_w)
    return res


if __name__ == '__main__':
    seq = Seq(
        'TACTTTTCTAATATCACGAGGACTTACATGGCCTCAAGTCACCTGTGGTGTTGTGCAAGAAGGAGAAGCAAAGTCTGTCTATGTATTATGAGATAGCTACTTCTATGGCTAGGATATATGTTGTACAAGACCGGCTTTTCTTCTACTTCTTGCACAACCTGAGGTTATTGAGGCTATACAAGTCTTCTTCTATAATGTTATTTATTA')
    out = get_Ewtf(get_numseq(seq, 1), wavlet)
    for o in out:
        print("%0.8f" % o, end="\t")
