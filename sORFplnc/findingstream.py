from Bio import SeqIO

def get_pos_sORFs():
    RNA = []
    IDs = []
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\positive_cdhit.fasta",   #ppa gma osa zma
                                  "fasta"):
        IDs.append(seq_record.id)
        RNA.append(str(seq_record.seq))

    results = []
    for i, j in zip(IDs, RNA):
        pos = int(i.strip().split("_")[2])
        end = int(i.strip().split("_")[3])
        length = end-pos
        if pos >= 174:
            result = '>' + i + '\n' + j[174:174+length] + '\n'
        else:
            result = '>' + i + '\n' + j[pos:end] + '\n'
        results.append(result)

    print(results[:5])
    print(len(results))

    f = open('G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\pos_sorfs_candidate.fa',
             'w')
    f.writelines(results)
    f.close()

def get_IRES():
    RNA = list()
    IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\gma\\negative_cdhit.fasta", "fasta"): # G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\negative_cdhit.fa    G:\赵思远\miRNA-encoded peptides\FslHID\数据集\positive_cdhit.fa
        IDs.append(seq_record.id)
        RNA.append(str(seq_record.seq))

    print(len(IDs))

    results = []
    for i, j in zip(IDs, RNA):
        pos = int(i.strip().split("_")[2])
        if pos >= 174:
            result = '>' + i + '\n' + j[:174] + '\n'
        else:
            result = '>' + i + '\n' + j[:pos] + '\n'
        results.append(result)

    print(len(results))

    f = open('G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\gma\\stream\\upstream_neg.fa', 'w') # G:\赵思远\miRNA-encoded peptides\FslHID\数据集\stream\\upstream_neg.fa    G:\赵思远\miRNA-encoded peptides\FslHID\数据集\stream\\upstream_pos.fa
    f.writelines(results)
    f.close()

def filter_sORFs_from_ires():
    sorfs = list()
    sorfs_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\negative_cdhit.fasta", "fasta"): # G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\negative_cdhit.fa    G:\赵思远\miRNA-encoded peptides\FslHID\数据集\positive_cdhit.fa
        sorfs_IDs.append(seq_record.id)
        sorfs.append(str(seq_record.seq))

    ires = list()
    ires_IDs = list()
    with open("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\stream\\upstream_mode_1_neg.result", "r") as files: # G:\赵思远\miRNA-encoded peptides\FslHID\数据集\stream\\upstream_mode_1_neg.result    G:\赵思远\miRNA-encoded peptides\FslHID\数据集\stream\\upstream_mode_1.result
        if True:
            next(files)
        for file in files:
            s = file.strip().split('\t')
            if file[0] == '>':
                ires_IDs.append(s[0])
            else:
                ires.append(s[1])

    # results = ['>' + i + '\n' + j + '\n' for i, j, k in zip(sorfs_IDs, sorfs, ires) if k != 'non-IRES']
    results = ['>' + i + '\n' + j + '\n' for i, j, k in zip(sorfs_IDs, sorfs, ires) if k == 'non-IRES']
    print(len(results))

    f = open('G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\zma_neg.fa', 'w')

    f.writelines(results)
    f.close()

def filter_pos_sORFs_from_ires_and_blast():
    sorfs_IDs = list()
    with open("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\blastx_pos_report.tsv", "r") as files:  #ppa gma osa zma
        for file in files:
            if file[0] == '#':
                continue
            else:
                ids = file.strip().split("\t")[0]
                sorfs_IDs.append(ids)

    sorfs_IDs = list(set(sorfs_IDs))
    print(len(sorfs_IDs))

    alls = list()
    alls_IDs = list()
    for seq_record in SeqIO.parse("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\positive_cdhit.fasta", "fasta"):  #ppa gma osa zma
        alls_IDs.append(seq_record.id)
        alls.append(str(seq_record.seq))

    ires = list()
    ires_IDs = list()
    with open("G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\stream\\upstream_mode_1_pos.result", "r") as files:  #ppa gma osa zma
        if True:
            next(files)
        for file in files:
            s = file.strip().split('\t')
            if file[0] == '>':
                ires_IDs.append(s[0])
            else:
                ires.append(s[1])

    results = ['>' + i + '\n' + j + '\n' for i, j, k in zip(alls_IDs, alls, ires) if i in sorfs_IDs and k != 'non-IRES']
    print(len(results))

    f = open('G:\赵思远\miRNA-encoded peptides\FslHID\数据集\\zma\\zma_pos.fa', 'w')   #ppa gma osa zma

    f.writelines(results)
    f.close()

# get_IRES()
# filter_sORFs_from_ires()
# get_pos_sORFs()
# filter_pos_sORFs_from_ires_and_blast()