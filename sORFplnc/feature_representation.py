import numpy as np
import math
import re
from Bio import SeqIO
import random

# feature extraction

def countnum(seq, nuacid):
    return len([1 for i in range(len(seq)) if seq.startswith(nuacid, i)])

def construct_kmer():
    ntarr = ("A", "C", "G", "T")

    kmerArray = []

    for n in range(4):
        kmerArray.append(ntarr[n])

    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            kmerArray.append(str2)
    #############################################
    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            for x in range(4):
                str3 = str2 + ntarr[x]
                kmerArray.append(str3)
    #############################################
    # change this part for 3mer or 4mer
    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            for x in range(4):
                str3 = str2 + ntarr[x]
                for y in range(4):
                    str4 = str3 + ntarr[y]
                    kmerArray.append(str4)
    ############################################
    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            for x in range(4):
                str3 = str2 + ntarr[x]
                for y in range(4):
                    str4 = str3 + ntarr[y]
                    for z in range(4):
                        str5 = str4 + ntarr[z]
                        kmerArray.append(str5)
    # ####################### 6-mer ##############
    # for n in range(4):
    #     str1 = ntarr[n]
    #     for m in range(4):
    #         str2 = str1 + ntarr[m]
    #         for x in range(4):
    #             str3 = str2 + ntarr[x]
    #             for y in range(4):
    #                 str4 = str3 + ntarr[y]
    #                 for z in range(4):
    #                     str5 = str4 + ntarr[z]
    #                     for t in range(4):
    #                         str6 = str5 + ntarr[t]
    #                         kmerArray.append(str6)
    return kmerArray

def get_kmer(seq,kmerArray):
    rst = []
    total = 0.0
    for n in range(len(kmerArray)):
        item = countnum(seq, kmerArray[n])
        total = total + item
        rst.append(item)
    for n in range(len(rst)):
        if total!=0:
            rst[n] = rst[n]/total

    return rst

def kmer_encode(seq, kmerArray):
    seq_data = []
    for i in range(len(seq)):
        seq_feature = get_kmer(seq[i], kmerArray)
        seq_data.append(seq_feature)
        if i == 5000 or i == 10000 or i == 20000:
            print(i)
    seq_data = np.array(seq_data)
    return seq_data

def get_Value(arrpc, arrnc, i, j):
    m = 0.0
    if arrnc[i, j] == 0.0:
        m == 1.0
    else:
        m = arrpc[i, j] / arrnc[i, j]
    if m == 0.0:
        value = 0.0
    else:
        value = math.log(m)
    return value

def nucleic_biasin(seq, arrnc, arrpc):
    if len(seq) == 6:
        bia = np.zeros((6))
        for i in range(6):
            if seq[i] == 'A':
                bia[i] = get_Value(arrpc, arrnc, i, j=0)
            elif seq[i] == 'C':
                bia[i] = get_Value(arrpc, arrnc, i, j=1)
            elif seq[i] == 'G':
                bia[i] = get_Value(arrpc, arrnc, i, j=2)
            elif seq[i] == 'T':
                bia[i] = get_Value(arrpc, arrnc, i, j=3)
        bia0 = sum(bia) / 6
    else:
        bia0 = 0
    return bia0

def nucleic_bia(seq, arrnc, arrpc):
    bia = np.zeros((len(seq)))
    for i in range(len(seq)):
        bia[i] = nucleic_biasin(seq[i], arrnc, arrpc)
    bia0 = np.expand_dims(bia, axis=1)
    return bia0

def orf_single(seq):
    startss = []
    stopss = []
    starts_ = []
    stop_s = []
    start = 0
    newseq = seq
    newseq_6 = []
    l = len(seq)
    for i in range(len(seq)):
        if (seq[i:i + 3] == "ATG"):
            start = 1  # has start codon
            for j in range(int((len(seq) - (i + 3)) / 3)):
                if (seq[i + 3 + 3 * j:i + 3 + 3 * j + 3] == "TAA") or (
                        seq[i + 3 + 3 * j:i + 3 + 3 * j + 3] == "TAG") or (
                        seq[i + 3 + 3 * j:i + 3 + 3 * j + 3] == "TGA"):
                    startss.append(i)
                    stopss.append(i + 3 + 3 * j + 3)
                    break
            if len(startss) == 0:
                starts_.append(i)

    if start == 0:
        for k in range(len(seq)):
            if (seq[k:k + 3] == "TAA") or (seq[k:k + 3] == "TAG") or (seq[k:k + 3] == "TGA"):
                stop_s.append(k + 3)

    if len(startss) != 0:
        startss = np.array(startss)
        stopss = np.array(stopss)
        coding_len = stopss - startss
        max_len_position = np.argmax(coding_len)
        newseq = seq[(startss[max_len_position]):(stopss[max_len_position])]
        if (startss[max_len_position] - 3) >= 0 and (startss[max_len_position] + 5) < l:
            newseq_6 = seq[(startss[max_len_position] - 3): (startss[max_len_position])] + seq[(startss[max_len_position] + 3):(startss[max_len_position] + 6)]

    elif len(starts_) != 0:
        starts_ = np.array(starts_)
        newseq = seq[(starts_[0]):len(seq)]
        if (starts_[0] - 3) >= 0 and (starts_[0] + 5) < l:
            newseq_6 = seq[(starts_[0] - 3):(starts_[0])] + seq[(starts_[0] + 3):(starts_[0] + 6)]

    elif len(stop_s) != 0:
        stop_s = np.array(stop_s)
        newseq = seq[0:(stop_s[-1])]

    return newseq, newseq_6


# single nucleic ggap
def g_gap_single(seq, ggaparray, g):
    # seq length is fix =23

    rst = np.zeros((16))
    for i in range(len(seq) - 1 - g):
        str1 = seq[i]
        str2 = seq[i + 1 + g]
        idx = ggaparray.index(str1 + str2)
        rst[idx] += 1

    for j in range(len(ggaparray)):
        rst[j] = rst[j] / (len(seq) - 1 - g)  # l-1-g

    return rst

def ggap_encode(seq, ggaparray, g):
    result = []
    for x in seq:
        temp = g_gap_single(x[21:], ggaparray, g)
        result.append(temp)
    result = np.array(result)
    return result

# binucleic ggap
# kmerarray[64:340]
def big_gap_single(seq, ggaparray, g):
    # seq length is fix =23

    rst = np.zeros((256))
    for i in range(len(seq) - 1 - g):
        str1 = seq[i] + seq[i + 1]
        str2 = seq[i + g] + seq[i + 1 + g]
        idx = ggaparray.index(str1 + str2)
        rst[idx] += 1

    for j in range(len(ggaparray)):
        rst[j] = rst[j] / (len(seq) - 1 - g)  # l-1-g

    return rst

def biggap_encode(seq, ggaparray, g):
    result = []
    for x in seq:
        temp = big_gap_single(x[:200], ggaparray, g)
        result.append(temp)
    result = np.array(result)
    return result

def nucl_bias_arrays(RNA_pos, RNA_neg):
    nc_arr, pc_arr = list(), list()
    RNA_pos = [i for i in RNA_pos]
    RNA_neg = [i for i in RNA_neg]
    RNA_pos = np.array(RNA_pos).T.tolist()
    RNA_neg = np.array(RNA_neg).T.tolist()
    for i in RNA_pos:
        nucl = [i.count('A'), i.count('T'), i.count('G'), i.count('C')]
        prob = [i/sum(nucl) for i in nucl]
        pc_arr.append(prob)
    for i in RNA_neg:
        nucl = [i.count('A'), i.count('T'), i.count('G'), i.count('C')]
        prob = [i/sum(nucl) for i in nucl]
        nc_arr.append(prob)
    return nc_arr, pc_arr

def sorf_length(RNA_seq):
    length = [[len(i)] for i in RNA_seq]
    return length

def Hexamer_score(fastas, c_m, nc_m, k=6):
    hexamer_score_encoding = []
    # header = ['#', 'hexamer_score']
    # hexamer_score_encoding.append(header)
    ntarr = 'ACGT'
    hexnuc = [nt1 + nt2 + nt3 + nt4 + nt5 + nt6 for nt1 in ntarr for nt2 in ntarr for nt3 in ntarr for nt4 in ntarr for nt5 in ntarr for nt6 in ntarr]
    for seq in fastas:
        # name, seq = i[0], re.sub('-', '', i[1])
        # code = [name]
        code = []
        if len(seq) > 12:
            l = len(seq) - k + 1
            log_r = np.zeros((l-6))
            for j in range(3, l-3):
                tempseq = seq[j: j + k]
                idx = hexnuc.index(tempseq)
                Fc = float(c_m[idx])
                Fnc = float(nc_m[idx])
                if Fc == 0.0 and Fnc == 0.0:
                    log_r[j-3] = 0
                elif Fc == 0.0 and Fnc != 0.0:
                    log_r[j-3] = -1
                elif Fc != 0.0 and Fnc == 0.0:
                    log_r[j-3] = 1
                else:
                    m = Fc / Fnc
                    log_r[j-3] = math.log(m)
            miu = sum(log_r) / (l-6)
            code.append(miu)
        else:
            code.append(0)
        hexamer_score_encoding.append(code)
    return hexamer_score_encoding

def nucl_bias_arrays(RNA_pos, RNA_neg):
    lnc_arr, pc_arr = list(), list()
    pos, neg = [], []
    pos_seq, pos_ids = [i[0] for i in RNA_pos], [i[1] for i in RNA_pos]
    neg_seq, neg_ids = [i[0] for i in RNA_neg], [i[1] for i in RNA_neg]
    for i, j in zip(pos_ids, pos_seq):
        stat, _, _ = get_Position(i)
        if stat >= 174:
            result = '>' + i + '\n' + j[171:174] +j[177:180] + '\n'
        else:
            result = '>' + i + '\n' + j[stat-3:stat] + j[stat+3:stat+6] + '\n'
        pos.append(result)
    for i, j in zip(neg_ids, neg_seq):
        stat, _, _ = get_Position(i)
        if stat >= 174:
            result = '>' + i + '\n' + j[171:174] +j[177:180] + '\n'
        else:
            result = '>' + i + '\n' + j[stat-3:stat] + j[stat+3:stat+6] + '\n'
        neg.append(result)
    RNA_pos = np.array(pos).T.tolist()
    RNA_neg = np.array(neg).T.tolist()
    for i in RNA_pos:
        nucl = [i.count('A'), i.count('T'), i.count('G'), i.count('C')]
        prob = [i/sum(nucl) for i in nucl]
        pc_arr.append(prob)
    for i in RNA_neg:
        nucl = [i.count('A'), i.count('T'), i.count('G'), i.count('C')]
        prob = [i/sum(nucl) for i in nucl]
        lnc_arr.append(prob)
    return lnc_arr, pc_arr

def get_Position(i):
    stat = int(i.strip().split("_")[2])
    end = int(i.strip().split("_")[3])
    leg = end - stat
    return stat, end, leg

def RNA_devide(seq, id):
    upstream, nuc6, sorfs, downstream = [], [], [], []
    for i, j in zip(id, seq):
        stat, end, leg = get_Position(i)
        if stat >= 174:
            upstream.append(j[:174])
            nuc6.append(j[171:174] + j[177:180])
            sorfs.append(j[174:174+leg])
            downstream.append(j[174+leg:])
        else:
            upstream.append(j[:stat])
            nuc6.append(j[stat-3:stat] + j[stat+3:stat+6])
            sorfs.append(j[stat:end])
            downstream.append(j[end:])
    return upstream, nuc6, sorfs, downstream

def get_features(nc_arr, pc_arr, RNA_seq, c_m, nc_m):
    nc_arr, pc_arr = np.array(nc_arr), np.array(pc_arr)
    kmerArray = construct_kmer()
    seq, ids = [i[0] for i in RNA_seq], [i[1] for i in RNA_seq]
    upstream, nuc6, sorfs, downstream = RNA_devide(seq, ids)

    cpf = cppred(sorfs)
    print(cpf[:2])

    # upstream
    up_RNA_biggap3 = biggap_encode(upstream, kmerArray[84:340], 3)
    nucbia = nucleic_bia(nuc6, nc_arr, pc_arr)
    # sORFs
    RNA_kmer1 = kmer_encode(sorfs, kmerArray[0:4])
    RNA_kmer2 = kmer_encode(sorfs, kmerArray[4:20])
    RNA_kmer3 = kmer_encode(sorfs, kmerArray[20:84])
    sORFs_length = sorf_length(sorfs)
    Hscore = Hexamer_score(sorfs, c_m, nc_m, k=6)
    # downstream
    down_biggap2 = biggap_encode(downstream, kmerArray[84:340], 2)

    RNA_data = np.concatenate((up_RNA_biggap3, nucbia, RNA_kmer1, RNA_kmer2, RNA_kmer3, sORFs_length, Hscore, down_biggap2), axis=1)
    # RNA_data = np.concatenate(
    #     (sORFs_length, Hscore, cpf), axis=1)
    print(RNA_data.shape)
    return RNA_data

def read_Hexamer_files(flag):
    if flag == 'ath':
        rf = "..\data\ath\\feature\\ath_Hexamer.tsv"
    elif flag == 'gma':
        rf = "..\data\\gma\\feature\\gma_Hexamer.tsv"
    elif flag == 'osa':
        rf = "..\data\\osa\\feature\\osa_Hexamer.tsv"
    elif flag == 'zma':
        rf = "..\data\\zma\\feature\\zma_Hexamer.tsv"

    c_m, nc_m = [], []
    with open(rf, 'r', encoding='utf-16') as files:
        if True:
            next(files)
        for file in files:
            s = file.strip().split("\t")
            c_m.append(s[1])
            nc_m.append(s[2])

    return c_m, nc_m

def read_files_ath():
    pos_sorfs = list()
    pos_sorfs_IDs = list()
    for seq_record in SeqIO.parse("..\data\\ath\\ath_pos.fa", "fasta"):
        pos_sorfs_IDs.append(seq_record.id)
        pos_sorfs.append(str(seq_record.seq))

    neg_sorfs = list()
    neg_sorfs_IDs = list()
    for seq_record in SeqIO.parse("..\data\\ath\\ath_neg.fa", "fasta"):
        neg_sorfs_IDs.append(seq_record.id)
        neg_sorfs.append(str(seq_record.seq))

    pos_sorfs = np.column_stack([pos_sorfs, pos_sorfs_IDs]).tolist()
    neg_sorfs = np.column_stack([neg_sorfs, neg_sorfs_IDs]).tolist()

    # training-validation
    pos_sorfs_training = random.sample(pos_sorfs[:40], 40)
    neg_sorfs_training = random.sample(neg_sorfs[:100], 100)
    train_data = pos_sorfs_training + neg_sorfs_training
    train_label = [1] * len(pos_sorfs_training) + [0] * len(neg_sorfs_training)

    # testing
    test_pos_sorfs = random.sample(pos_sorfs[40:80], 40)
    test_neg_sorfs = random.sample(neg_sorfs[100:200], 100)
    test_data = test_pos_sorfs + test_neg_sorfs
    test_label = [1] * len(test_pos_sorfs) + [0] * len(test_neg_sorfs)

    return pos_sorfs_training, neg_sorfs_training, train_data, train_label, test_data, test_label

def read_files_gma():
    pos_sorfs = list()
    pos_sorfs_IDs = list()
    for seq_record in SeqIO.parse("..\data\\gma\\gma_pos.fa", "fasta"):
        pos_sorfs_IDs.append(seq_record.id)
        pos_sorfs.append(str(seq_record.seq))

    neg_sorfs = list()
    neg_sorfs_IDs = list()
    for seq_record in SeqIO.parse("..\data\\gma\\gma_neg.fa", "fasta"):
        neg_sorfs_IDs.append(seq_record.id)
        neg_sorfs.append(str(seq_record.seq))

    pos_sorfs = np.column_stack([pos_sorfs, pos_sorfs_IDs]).tolist()
    neg_sorfs = np.column_stack([neg_sorfs, neg_sorfs_IDs]).tolist()

    # training-validation
    pos_sorfs_training = random.sample(pos_sorfs[:30], 30)
    neg_sorfs_training = random.sample(neg_sorfs[:40], 40)
    train_data = pos_sorfs_training + neg_sorfs_training
    train_label = [1] * len(pos_sorfs_training) + [0] * len(neg_sorfs_training)

    # testing
    test_pos_sorfs = random.sample(pos_sorfs[30:54], 24)
    test_neg_sorfs = random.sample(neg_sorfs[40:64], 24)
    test_data = test_pos_sorfs + test_neg_sorfs
    test_label = [1] * len(test_pos_sorfs) + [0] * len(test_neg_sorfs)

    return pos_sorfs_training, neg_sorfs_training, train_data, train_label, test_data, test_label

def read_files_zma():
    pos_sorfs = list()
    pos_sorfs_IDs = list()
    for seq_record in SeqIO.parse("..\data\\zma\\zma_pos.fa", "fasta"):
        pos_sorfs_IDs.append(seq_record.id)
        pos_sorfs.append(str(seq_record.seq))

    neg_sorfs = list()
    neg_sorfs_IDs = list()
    for seq_record in SeqIO.parse("..\data\\zma\\zma_neg.fa", "fasta"):
        neg_sorfs_IDs.append(seq_record.id)
        neg_sorfs.append(str(seq_record.seq))

    pos_sorfs = np.column_stack([pos_sorfs, pos_sorfs_IDs]).tolist()
    neg_sorfs = np.column_stack([neg_sorfs, neg_sorfs_IDs]).tolist()

    # training-validation
    pos_sorfs_training = random.sample(pos_sorfs[:100], 100)
    neg_sorfs_training = random.sample(neg_sorfs[:1000], 1000)
    train_data = pos_sorfs_training + neg_sorfs_training
    train_label = [1] * len(pos_sorfs_training) + [0] * len(neg_sorfs_training)

    # testing
    test_pos_sorfs = random.sample(pos_sorfs[100:200], 100)
    test_neg_sorfs = random.sample(neg_sorfs[1000:2000], 1000)
    test_data = test_pos_sorfs + test_neg_sorfs
    test_label = [1] * len(test_pos_sorfs) + [0] * len(test_neg_sorfs)

    return pos_sorfs_training, neg_sorfs_training, train_data, train_label, test_data, test_label

def read_files_osa():
    pos_sorfs = list()
    pos_sorfs_IDs = list()
    for seq_record in SeqIO.parse("..\data\\osa\\osa_pos.fa", "fasta"):
        pos_sorfs_IDs.append(seq_record.id)
        pos_sorfs.append(str(seq_record.seq))

    neg_sorfs = list()
    neg_sorfs_IDs = list()
    for seq_record in SeqIO.parse("..\data\\osa\\osa_neg.fa", "fasta"):
        neg_sorfs_IDs.append(seq_record.id)
        neg_sorfs.append(str(seq_record.seq))

    pos_sorfs = np.column_stack([pos_sorfs, pos_sorfs_IDs]).tolist()
    neg_sorfs = np.column_stack([neg_sorfs, neg_sorfs_IDs]).tolist()

    # training-validation
    pos_sorfs_training = random.sample(pos_sorfs[:30], 30)
    neg_sorfs_training = random.sample(neg_sorfs[:70], 70)
    train_data = pos_sorfs_training + neg_sorfs_training
    train_label = [1] * len(pos_sorfs_training) + [0] * len(neg_sorfs_training)

    # testing
    test_pos_sorfs = random.sample(pos_sorfs[30:60], 30)
    test_neg_sorfs = random.sample(neg_sorfs[70:109], 39)
    test_data = test_pos_sorfs + test_neg_sorfs
    test_label = [1] * len(test_pos_sorfs) + [0] * len(test_neg_sorfs)

    return pos_sorfs_training, neg_sorfs_training, train_data, train_label, test_data, test_label


def run_feature_representation(read_files, hex_file):
    if read_files == 'ath':
        RNA_pos, RNA_neg, RNA_data, RNA_label_train, RNA_seq_test, RNA_label_test = read_files_ath()
    elif read_files == 'gma':
        RNA_pos, RNA_neg, RNA_data, RNA_label_train, RNA_seq_test, RNA_label_test = read_files_gma()
    elif read_files == 'osa':
        RNA_pos, RNA_neg, RNA_data, RNA_label_train, RNA_seq_test, RNA_label_test = read_files_osa()
    elif read_files == 'zma':
        RNA_pos, RNA_neg, RNA_data, RNA_label_train, RNA_seq_test, RNA_label_test = read_files_zma()
    c_m, nc_m = read_Hexamer_files(flag=hex_file)
    nc_arr, pc_arr = nucl_bias_arrays(RNA_pos, RNA_neg)
    train_feature = get_features(nc_arr, pc_arr, RNA_data, c_m, nc_m)
    test_feature = get_features(nc_arr, pc_arr, RNA_seq_test, c_m, nc_m)

    return train_feature, RNA_label_train, test_feature, RNA_label_test