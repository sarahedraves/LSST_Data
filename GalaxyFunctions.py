import numpy as np

def galaxy_match(galaxy_calc, galaxy_truth):
    print(galaxy_calc)
    print(galaxy_truth)
    condlist = [np.logical_and(galaxy_calc==1, galaxy_truth==1),
            np.logical_and(galaxy_calc==1, galaxy_truth!=1),
            np.logical_and(galaxy_calc==0, galaxy_truth==1),
            np.logical_and(galaxy_calc==0, galaxy_truth!=1)]
    choicelist = ['true positive','false positive','false negative','true negative']
    matched_array=np.select(condlist, choicelist, np.NaN)
    print(matched_array)
    return matched_array

def confusion_matrix(matched_array, style='percents'):
    total=matched_array.size
    tp=np.count_nonzero(matched_array=='true positive')
    fp=np.count_nonzero(matched_array=='false positive')
    fn=np.count_nonzero(matched_array=='false negative')
    tn=np.count_nonzero(matched_array=='true negative')
    sums_matrix={'total':total, 'true positive':tp, 'false positive':fp, 'false negative':fn, 'true negative':tn}
    percents_matrix={'true positive %':("{:.1%}".format(tp/total)), 'false positive %':("{:.1%}".format(fp/total)), 
                     'false negative %':("{:.1%}".format(fn/total)), 'true negative %':("{:.1%}".format(tn/total))}
    if style=='percents':
        return percents_matrix
    elif style=='sums':
        return sums_matrix
    else:
        pass