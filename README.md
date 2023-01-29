# sORFplnc
Identification of small open reading frames in plant lncRNA using class-imbalance learning

An alternative SMOTE based on weighted cosine distance (WCDSMOTE) which enables interaction with feature selection is put forward to synthesize minority class samples. Weighted edited nearest neighbor (WENN) is applied to clean up majority class samples, thus, hybrid sampling WCDSMOTE-ENN is proposed. Therefore, we present a novel computational method which is based on class-imbalance learning to identify the sORFs with coding potential in plant lncRNA (sORFplnc).

# OS
win10

# Dependencies
Language dependency: Python 3 (Please do not use Python 2 to run the code.)

# Library dependency

python==3.6.15

numpy==1.19.2

scikit_learn==0.24.2

# Usage

python sORFplnc.py

Also，you can retrain the model using your data.

# Graphical abstract

![image](https://github.com/zzssyy/sORFplnc/blob/main/Graphical_Abstract.png)
