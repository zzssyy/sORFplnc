from Bio import SeqIO
import numpy as np
import random

def read_files():
    pos_sorfs = list()
    pos_sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\ath\\ath_pos.fa", "fasta"):
        pos_sorfs_IDs.append(seq_record.id)
        pos_sorfs.append(str(seq_record.seq))

    neg_sorfs = list()
    neg_sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\ath\\ath_neg.fa", "fasta"):
        neg_sorfs_IDs.append(seq_record.id)
        neg_sorfs.append(str(seq_record.seq))

    pos_sorfs = np.column_stack([pos_sorfs, pos_sorfs_IDs]).tolist()
    neg_sorfs = np.column_stack([neg_sorfs, neg_sorfs_IDs]).tolist()

    # training-validation
    pos_sorfs_training = random.sample(pos_sorfs[:40], 40)
    neg_sorfs_training = random.sample(neg_sorfs[:100], 100)

    pos_sorfs = list()
    pos_sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\gma\\gma_pos.fa", "fasta"):
        pos_sorfs_IDs.append(seq_record.id)
        pos_sorfs.append(str(seq_record.seq))

    neg_sorfs = list()
    neg_sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\gma\\gma_neg.fa", "fasta"):
        neg_sorfs_IDs.append(seq_record.id)
        neg_sorfs.append(str(seq_record.seq))

    pos_sorfs = np.column_stack([pos_sorfs, pos_sorfs_IDs]).tolist()
    neg_sorfs = np.column_stack([neg_sorfs, neg_sorfs_IDs]).tolist()

    # training-validation
    pos_sorfs_training = pos_sorfs_training + random.sample(pos_sorfs[:30], 30)
    neg_sorfs_training = neg_sorfs_training + random.sample(neg_sorfs[:40], 40)

    pos_sorfs = list()
    pos_sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\zma_pos.fa", "fasta"):
        pos_sorfs_IDs.append(seq_record.id)
        pos_sorfs.append(str(seq_record.seq))

    neg_sorfs = list()
    neg_sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\zma_neg.fa", "fasta"):
        neg_sorfs_IDs.append(seq_record.id)
        neg_sorfs.append(str(seq_record.seq))

    pos_sorfs = np.column_stack([pos_sorfs, pos_sorfs_IDs]).tolist()
    neg_sorfs = np.column_stack([neg_sorfs, neg_sorfs_IDs]).tolist()

    # training-validation
    pos_sorfs_training = pos_sorfs_training + random.sample(pos_sorfs[:100], 100)
    neg_sorfs_training = neg_sorfs_training + random.sample(neg_sorfs[:1000], 1000)

    pos_sorfs = list()
    pos_sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\osa\\osa_pos.fa", "fasta"):
        pos_sorfs_IDs.append(seq_record.id)
        pos_sorfs.append(str(seq_record.seq))

    neg_sorfs = list()
    neg_sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\osa\\osa_neg.fa", "fasta"):
        neg_sorfs_IDs.append(seq_record.id)
        neg_sorfs.append(str(seq_record.seq))

    pos_sorfs = np.column_stack([pos_sorfs, pos_sorfs_IDs]).tolist()
    neg_sorfs = np.column_stack([neg_sorfs, neg_sorfs_IDs]).tolist()

    # training-validation
    pos_sorfs_training = pos_sorfs_training + random.sample(pos_sorfs[:30], 30)
    neg_sorfs_training = neg_sorfs_training + random.sample(neg_sorfs[:70], 70)

    train_data = pos_sorfs_training + neg_sorfs_training
    train_label = [1] * len(pos_sorfs_training) + [0] * len(neg_sorfs_training)
    print(len(pos_sorfs_training))
    print(len(neg_sorfs_training))

    return train_data, train_label


def write_files(data, label):
    results_pos = []
    results_neg = []
    for i, j in zip(data, label):
        if j == 1:
            results_pos.append(">" + i[1] + "\n" + i[0] + '\n')
        else:
            results_neg.append(">" + i[1] + "\n" + i[0] + '\n')

    print(len(results_pos))
    print(len(results_neg))

    f = open('G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\all\\all_pos.fa', 'w')
    f.writelines(results_pos)
    f.close()

    f = open('G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\all\\all_neg.fa', 'w')
    f.writelines(results_neg)
    f.close()

data, label = read_files()
write_files(data, label)