import numpy as np
import pandas as pd
from math import log,exp,inf,pi
from sklearn.metrics import mean_squared_error
from statistics import mean,variance
import shap
from sklearn.model_selection import cross_validate, train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

def revstrand(inseq):
    outseq=inseq[::-1]
    return outseq


def complement(inseq):
    '''
    complement for both complete and incomplete nucleotides.
    '''
    inseq.upper()
    complement_dict = {'A':'T','T':'A','C':'G','G':'C',
                       'B':'V','D':'H','H':'D','K':'M','M':'K',
                       'S':'S','V':'B','W':'W','N':'N','R':'Y','Y':'R' }
    clist=[]
    for ntide in inseq:
        clist.append( complement_dict[ ntide ] )
        
    complement=''.join(clist)
    return complement

#calculate hamming distance
def hmd(seq1,seq2):
    hamd = 0
    for i in range(0,len(seq1)):
        if seq1[i] != seq2[i]:
            hamd += 1
    return hamd

def score_motif(s, M):
    '''
    score nucleotide sequence: https://www.bioinformatics.org/sms/iupac.html
    '''
    if len(s) != len(M):
        return -inf
        
    out = 0
    for i in range( len(s)):
        if s[i] == 'N':
            out = out + (M[i]['A'] + M[i]['C'] + M[i]['G']+M[i]['T'])/4
        elif s[i] == 'R':
            out = out + (M[i]['A'] + M[i]['G'] )/2
        elif s[i] == 'Y':
            out = out + (M[i]['C'] + M[i]['T'] )/2
        elif s[i] == 'S':
            out = out + (M[i]['C'] + M[i]['G'] )/2
        elif s[i] == 'W':
            out = out + (M[i]['A'] + M[i]['T'] )/2
        elif s[i] == 'K':
            out = out + (M[i]['T'] + M[i]['G'] )/2
        elif s[i] == 'M':
            out = out + (M[i]['C'] + M[i]['A'] )/2
        elif s[i] == 'B':
            out = out + (M[i]['C'] + M[i]['G'] +  M[i]['T'] )/3
        elif s[i] == 'D':
            out = out + (M[i]['A'] + M[i]['G'] +  M[i]['T'] )/3
        elif s[i] == 'H':
            out = out + (M[i]['A'] + M[i]['C'] +  M[i]['T'] )/3
        elif s[i] == 'V':
            out = out + (M[i]['C'] + M[i]['G'] +  M[i]['A'] )/3
        else:
            out = out + M[i][s[i]]
    return out

def rpwm(fname):
    logodds = {}
    ffile = open(fname, "rt")
    lines = ffile.readlines()
    ffile.close()
    
    for line in lines:
        sline = line.strip()
        if sline[0] == '>':
            name = sline.split()[0][1:]
            size = int(sline.split()[1])
            logodds[name] = [ {'A':0,'C':0,'G':0,'T':0} for i in range(size)]
        else:
            slist = sline.split()
            letter = slist[0]
            vlist = slist[1:]
            for i in range(size):
                logodds[name][i][letter] = float(vlist[i])
            
    return logodds

def rlogodds(fname):
    logodds = {}
    ffile = open(fname, "rt")
    lines = ffile.readlines()
    ffile.close()
    
    for line in lines:
        sline = line.strip()
        if sline[0] == '>':
            name = sline[2:].strip()
            logodds[name] = []
        elif sline[0:3] == 'log':
            continue
        else:
            temp = {}
            lineList = sline.split()
            temp['A']= int(lineList[0])
            temp['C']=int(lineList[1])
            temp['G']= int(lineList[2])
            temp['T']=int(lineList[3])
            logodds[name].append( temp )   
    
    return logodds


def get_shape( start, end, shape_table ):
    '''
    shape table should imported.
    '''
    start = start -3
    end = end -2
    result = shape_table.iloc[start:end,:]
    return result


def multi_gaussian(x,mean,sigma):
    '''
    compute multi variate gaussian likelihood,
    to avoid non-invertible covariance matrix, pseudo-inverse is used.
    '''
    if len(x) != len(mean):
        print("error");
        return -inf
    
    k = len(x)
    xc = np.subtract(x,mean)
    sig_inv = np.linalg.pinv(sigma);
    score = (2*pi)**(-k/2)*(np.linalg.det(sigma))**(-1/2)*exp( (-1/2) * np.dot( np.dot(xc.T,sig_inv), xc  ) )
    return score
    
def mse_population( subbox, group ):
    '''
    average normalized mse over the whole group of shape vectors.
    '''
    result = []
    for vlist in group:
        range_vlist = max(vlist) - min(vlist)
        result.append(  (1/range_vlist)*(mean_squared_error( subbox, vlist ))**0.5    )
    mean_result = sum(result)/len(result)
    
    return mean_result
    
def nmse( subbox, consensus ):
    range_con = max( consensus ) - min( consensus )
    result = (  (1/range_con)*(mean_squared_error( subbox, consensus ))**0.5    )
    
    return result  
    
def nmse_border( subbox, consensus ):
    '''
    mse to compare shape features on binding box's left and right border.
    '''
    sub_border = subbox[0:4] + subbox[len(subbox)-4:len(subbox)]
    con_border = consensus[0:4] + consensus[len(consensus)-4:len(consensus)]
    range_con = max( con_border ) - min( con_border )
    result = (  (1/range_con)*(mean_squared_error( sub_border, con_border  ))**0.5    )
    
    return result   

def PSSMgenerator( count:list ):
    M = []
    width = len(count);
    nu = ['A','T','C','G'];
    for i in range(width):
        temp_dict = {};
        temp_sum = 4 + count[i]['A'] + count[i]['T'] + count[i]['C'] + count[i]['G']
        for n in nu:
            temp_dict[n] = log( (count[i][n] + 1)/temp_sum )
        M.append( temp_dict)
    
    return M        
            
    
#Sizhe Qiu
def motif_match_noshape( tss:int , strand:int,N_up:int, N_down:int,bitome_obj,pwm:dict, name:str, isIM: bool )->dict:
    #M is the log odds matrix for the motif.
    """
    Return log odds score of matched motif, and whether it locates upstream to TSS, range is -100 to +50.
    :param dict M: PWM matrix of motifs
    :param TrancriptionUnit tu: the transcription unit to match the motif
    :param str name: name of the motif, either TF or im motif
    :return dict result: log odds score of matched motif, and whether it locates upstream to TSS
    """

    M = pwm[name]
    motif_len = len(M)


    if int(strand) == -1:
        flag = False
    else:
        flag = True

    TSS = int( tss )   
    if flag:
        left = TSS - N_up
        right = TSS+ N_down
        s = str(bitome_obj.sequence)[left:right]
    else:
        left = TSS-N_down
        right = TSS+ N_up
        s = complement(revstrand( str( bitome_obj.sequence)[left:right] ))


    scores = []
    for i in range(0, len(s)+1 - motif_len ):
        scores.append( score_motif(s[i:i+motif_len], M) )

    smax = max(scores)

    for i in range ( 0, len(s)+1 - motif_len ):
        if scores[i] == smax:
            s_motif = s[i:i+motif_len]
            index = i

    if flag:
        start = index + left
        end = start + motif_len - 1 #inclusive
    else:
        end = right - index
        start = end - motif_len + 1 #inclusive

    if flag:
        if (end+start)/2 <= int( tss ):
            up = 1
        else:
            up = 0
    else:
        if  (end+start)/2 >= int( tss ):
            up = 1
        else:
            up = 0

    
    if isIM:
        name = name + ' im'


    result = { name+' Matched Motif': s_motif ,name+' score': smax, 
              name +' upstreamTSS': up, name+' distToTSS': abs( TSS-(end+start)/2 ) }
    

    return result   

def ld( X, labels ):
    '''
    labels is binary
    '''
    X1,X0 = [],[]
    for i in range( len(labels) ):
        if labels[i] == 1:
            X1.append( X[i] )
        else:
            X0.append( X[i] )
    X1 = np.array(X1).T
    X0 = np.array(X0).T
    Sw = ( np.cov(X1) + np.cov(X0)  )
    mu1 = np.mean(X1, axis = 1)
    mu0 = np.mean( X0, axis = 1)
    w =  np.matmul (  np.linalg.pinv( Sw ), np.subtract(mu0, mu1).T    ) 
              
    return w
    
def get_seq_shape( seq_whole, left:int, right:int, data_dir ):
    '''
    left and right are cardinal position, including left and right bps
    '''
    if left < 4 or right > len(seq_whole) - 4 :
        return None
    cardinal_pos = np.arange(left,right+1)
    
    shape5mer_table = {}
    for shape_name in ['HelT','MGW','Roll','ProT']:
        shape5mer_table[shape_name] = pd.read_csv(data_dir + '/shapeTables/'+shape_name+'.5mer.csv',header=None)
    L = left -1
    R = right -1
    shape_vectors = { 'HelT':[],'MGW':[],'Roll':[],'ProT':[] }
    
    inter_tuple = { 'HelT':[],'Roll':[] }
    for i in range(L-1, right+1 ):
        kmer = seq_whole[ i-2:i+3 ]
        for shape_name in ['HelT','MGW','Roll','ProT']:
            temp = shape5mer_table[shape_name][shape5mer_table[shape_name][0] == str(kmer)]
            if shape_name in ['MGW','ProT']:
                shape_vectors[shape_name].append( list(temp[1])[0] )
            else:
                inter_tuple[shape_name].append( tuple( (list(temp[1])[0], list(temp[3])[0]) ) )
    shape_vectors['MGW'] = shape_vectors['MGW'][1:-1]
    shape_vectors['ProT'] = shape_vectors['ProT'][1:-1]
    
    for i in range( len(inter_tuple['Roll'])-2 ):
        for  shape_name in ['HelT','Roll']:
            shape_vectors[shape_name].append(  0.25*(inter_tuple[shape_name][i][1] + inter_tuple[shape_name][i+1][0] +
                                               inter_tuple[shape_name][i+1][1] + inter_tuple[shape_name][i+2][0] )  )
    table = pd.DataFrame( shape_vectors )
    table['cardinal_pos'] =  cardinal_pos
    table = table.set_index( 'cardinal_pos' )
    
    return table


def motif_match_withshape( tss:int , strand:int,N_up:int, N_down:int,bitome_obj,pwm:dict, name:str, isIM: bool, data_dir )->dict:
    #M is the log odds matrix for the motif.
    """
    Return log odds score of matched motif, and whether it locates upstream to TSS, range is -100 to +50.
    :param dict M: PWM matrix of motifs
    :param TrancriptionUnit tu: the transcription unit to match the motif
    :param str name: name of the motif, either TF or im motif
    :return dict result: log odds score of matched motif, and whether it locates upstream to TSS
    """

    M = pwm[name]
    motif_len = len(M)


    if int(strand) == -1:
        flag = False
    else:
        flag = True

    TSS = int( tss )   
    if flag:
        left = TSS - N_up
        right = TSS+ N_down
        s = str(bitome_obj.sequence)[left:right]
    else:
        left = TSS-N_down
        right = TSS+ N_up
        s = complement(revstrand( str( bitome_obj.sequence)[left:right] ))


    scores = []
    for i in range(0, len(s)+1 - motif_len ):
        scores.append( score_motif(s[i:i+motif_len], M) )

    smax = max(scores)

    for i in range ( 0, len(s)+1 - motif_len ):
        if scores[i] == smax:
            s_motif = s[i:i+motif_len]
            index = i

    if flag:
        start = index + left
        end = start + motif_len - 1 #inclusive
    else:
        end = right - index
        start = end - motif_len + 1 #inclusive

    if flag:
        if (end+start)/2 <= int( tss ):
            up = 1
        else:
            up = 0
    else:
        if  (end+start)/2 >= int( tss ):
            up = 1
        else:
            up = 0
            
    table = get_seq_shape( str(bitome_obj.sequence),start,end, data_dir)
    helt =  list(table['HelT'])
    roll =  list(table['Roll'])
    mgw = list(table['MGW'])
    prot=  list(table['ProT'])
    
    if isIM:
        name = name + ' im'


    result = { name+' Matched Motif': s_motif ,name+' score': smax, 
              name +' upstreamTSS': up, name+' distToTSS': abs( TSS-(end+start)/2 ),
              name+' HelT_avg': mean(helt),name+' HelT_range': max(helt)-min(helt),
                  name+' MGW_avg': mean(mgw),name+' MGW_range': max(mgw) - min(mgw),
                  name+' ProT_avg': mean(prot), name+' ProT_range': max(prot) - min(prot),
                  name+' Roll_avg': mean(roll), name+' Roll_range': max(roll) - min(roll),
                 name+' HelT_max': max(helt),name+' HelT_min': min(helt),
                  name+' MGW_max': max(mgw),name+' MGW_min': min(mgw),
                  name+' ProT_max': max(prot), name+' ProT_min': min(prot),
                  name+' Roll_max': max(roll), name+' Roll_min': min(roll)
             }
    

    return result   
        
def feature_importance(x: pd.DataFrame, y: pd.Series, model, model_type: str = 'tree'):
    """
    Given a specific xy and a model, use Shapley values to visualize the feature importances. This
    function will automatically split the provided x and y into training and validation sets, a
    step required for using the Shapley values package.

    :param pd.DataFrame x: the X matrix
    :param pd.Series y: the targets to use for cross-validated training
    :param model: any sklearn model for which to determine feature importances
    :param str model_type: the type of model in use; options are 'tree', 'other'; this
        will influence the choice of explainer from the shap package that is used
    """

    x_train, x_val, y_train, y_val = train_test_split(x, y, test_size=0.2)
    fit_model = model.fit(x_train, y_train)

    if model_type == 'tree':
        explainer = shap.TreeExplainer(fit_model)
    else:
        explainer = shap.Explainer(fit_model)
    shap_values = explainer.shap_values(x_val)
    shap.summary_plot(shap_values, x_val)  
    
    
    
    
    