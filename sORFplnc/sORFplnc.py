from sklearn.svm import SVC
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sORFplnc.feature_representation import *
from sklearn.linear_model import LogisticRegression as LR
from sklearn.externals import joblib
from CDSSMOTEENN.HybridSampling import *
from sklearn.preprocessing import minmax_scale
from sklearn.ensemble import StackingClassifier

def G_mean(matrics):
    import math
    TP, TN, FP, FN = matrics[1][1], matrics[0][0], matrics[0][1], matrics[1][0]
    value = 100*math.sqrt(TP*TN/((TP+FN)*(TN+FP)))
    return value

# single classifier
def svc(traindata, testdata, trainlabel, testlabel):
    print("Start training SVM...")
    svc = SVC(kernel="linear", C=1.0, cache_size='3000')
    svc.fit(traindata, trainlabel)

    pred_testlabel = svc.predict(testdata)
    g = G_mean(metrics.confusion_matrix(testlabel, pred_testlabel))
    # print('accuracy', metrics.accuracy_score(testlabel, pred_testlabel))
    print(metrics.confusion_matrix(testlabel, pred_testlabel))
    print('recall', 100*metrics.recall_score(testlabel, pred_testlabel, average='weighted'))
    print('precision', 100*metrics.precision_score(testlabel, pred_testlabel, average='weighted'))
    print('f1-score:', 100*metrics.f1_score(testlabel, pred_testlabel, average='weighted'))
    print("G_mean:", g)

def lgb(traindata, testdata, trainlabel, testlabel):
    print("Start training lightGBM...")
    lg = LGBMClassifier(n_estimators=100)
    lg.fit(traindata, trainlabel)

    pred_testlabel = lg.predict(testdata)
    g = G_mean(metrics.confusion_matrix(testlabel, pred_testlabel))
    print(metrics.confusion_matrix(testlabel, pred_testlabel))
    print('recall', 100*metrics.recall_score(testlabel, pred_testlabel, average='weighted'))
    print('precision', 100*metrics.precision_score(testlabel, pred_testlabel, average='weighted'))
    print('f1-score:', 100*metrics.f1_score(testlabel, pred_testlabel, average='weighted'))
    print("G_mean:", g)

# ensemble learning
def ensemble_pred(x_train, x_test, y_train, y_test, flag = 0, species='ath'):
    print('starting ensemble learning...')
    print("Heterogeneous ...")
    estimators = [('svr', SVC(kernel="linear", cache_size=3000, probability=True)),
                  ('lg', LGBMClassifier(n_estimators=100))]

    # ensemble stacking method
    print('starting stacking ensemble...')
    sclf = StackingClassifier(estimators=estimators)
    sclf.fit(np.array(x_train), y_train)
    pred_testlabel = sclf.predict(np.array(x_test))
    g = G_mean(metrics.confusion_matrix(y_test, pred_testlabel))
    print(metrics.confusion_matrix(y_test, pred_testlabel))
    print('recall', 100*metrics.recall_score(y_test, pred_testlabel, average='weighted'))
    print('precision', 100*metrics.precision_score(y_test, pred_testlabel, average='weighted'))
    print('f1-score:', 100*metrics.f1_score(y_test, pred_testlabel, average='weighted'))
    print("G_mean:", g)
    fileb = species + '_FslHID_stacking_balanced.m'
    fileib = species + '_FslHID_stacking_imbalanced.m'
    if flag == 0:
        joblib.dump(sclf, fileb)
    else:
        joblib.dump(sclf, fileib)


# on the training set
def training(read_files, hex_file):
    sp = read_files
    train_feature, train_label, test_feature, test_label = run_feature_representation(read_files, hex_file)
    train_feature, train_label, test_feature, test_label = minmax_scale(train_feature), train_label, minmax_scale(
        test_feature), test_label
    x_train, x_test, y_train, y_test = train_test_split(train_feature, train_label, test_size=0.2)


    x_train_hs, y_train_hs, u = HybridSampling(x_train, y_train)
    tf = test_feature.T
    xt = x_test.T
    x_test = np.array([xt[i] for i in range(len(xt)) if i in u]).T
    test_feature_1 = np.array([tf[i] for i in range(len(tf)) if i in u]).T
    ensemble_pred(x_train_hs, x_test, y_train_hs, y_test, flag=1, species=sp)

    return test_feature, test_label, test_feature_1, test_label


# on the ath test set
def independent_test(test_feature, test_label, flag=0, species='ath'):
    fileb = species + '_FslHID_stacking_balanced.m'
    fileib = species + '_FslHID_stacking_imbalanced.m'
    if flag == 0:
        model = joblib.load(fileb)
    else:
        model = joblib.load(fileib)
    y_pred = model.predict(test_feature)
    g = G_mean(metrics.confusion_matrix(test_label, y_pred))
    print(metrics.confusion_matrix(test_label, y_pred))
    print("recall:", 100*metrics.recall_score(test_label, y_pred))
    print("MCC:", 100 * metrics.matthews_corrcoef(test_label, y_pred))
    print("precision:", 100*metrics.precision_score(test_label, y_pred))
    print("f1_score:", 100*metrics.f1_score(test_label, y_pred))
    print("G_mean:", g)
    print("AUC:", 100*metrics.roc_auc_score(test_label, y_pred))


if __name__ == '__main__':
    read_files, hex_file = 'zma', 'zma'
    test_feature, test_label, test_feature_im, test_label = training(read_files, hex_file)
    independent_test(test_feature_im, test_label, flag=1, species=read_files)

