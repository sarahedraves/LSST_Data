"""
Function from Irene. Description:
Here’s the code I use to get stuff from GCRCatalogs. The first bit (above #ADD grabbing of truth info) gets stuff from the DC2 catalog you’re using. The next bit does the truth matching, and pulls out the truth quantities you want for comparison. The end just does some various cleaning of non-detections that are probably not of interest to you.
"""

def make_catalog_for_tract(gc,tract,verbose=False):
    """
    function to make a pandas dataframe with some basic info for a tract/patch for all six bands 
    inputs:
    gc: catalog reader
    tract: int; tract number
    """
    bands = ['u','g','r','i','z','y']
    columns = ['ra','dec','extendedness','blendedness','tract','patch','objectId']
    for band in bands:
        columns.append('mag_{}_cModel'.format(band))
        columns.append('magerr_{}_cModel'.format(band))
        columns.append('snr_{}_cModel'.format(band))
        columns.append('cModelFlux_{}'.format(band))
        columns.append('cModelFluxErr_{}'.format(band))
    tractfilt = f"tract=={tract}"
    #data = gc.get_quantities(columns, filters = simple_cuts, native_filters=tractfilt)
    data = gc.get_quantities(columns, native_filters=tractfilt)
    df = pd.DataFrame(data)
    #df.info()
    id = df['objectId']
    print(df['dec'][:5])
    ebv_vec = compute_ebv(df['ra'],df['dec'])
    df['ebv']=ebv_vec
    
    #ADD grabbing of truth info!
    truth_dir = "/global/cfs/cdirs/lsst/gsharing/lsstdesc-public/dc2/run2.2i-dr6-v2/truth_match"
    truth_file = f"truth_tract{tract}.parquet"
    truthpath = os.path.join(truth_dir, truth_file)
    truth_data = pd.read_parquet(truthpath)
    truthsz = np.array(truth_data['redshift'])
    matchid = np.array(truth_data['match_objectId'])
    isgoodmatch = np.array(truth_data['is_good_match'])
    isuniquematch = np.array(truth_data['is_unique_truth_entry'])
    true_r = np.array(truth_data['mag_r'])
    szdict = dict(zip(matchid, truthsz))
    magdict = dict(zip(matchid, true_r))
    isgooddict = dict(zip(matchid, isgoodmatch))
    isuniquedict = dict(zip(matchid, isuniquematch))
    
    goodtsz = [szdict[idx] for idx in id]
    goodtr = [magdict[idx] for idx in id]
    goodtgood = [isgooddict[idx] for idx in id]
    goodtunique = [isuniquedict[idx] for idx in id]
    df['true_redshift'] = goodtsz
    df['true_mag_r'] = goodtr
    df['is_truth_match_good'] = goodtgood
    df['is_truth_match_unique'] = goodtunique
    
    band_meanlam = [3671., 4827.,6223.,7546.,8691.,9710.] #just the mean for now
    band_a_ebv = np.array([4.81,3.74,2.70,2.06,1.58,1.31]) #A/E(B-V) quick calculation from CCM model!!!
    #grab the list of all available patches for the tract
    patches = list(set(df['patch']))
    
    for ii,band in enumerate(bands):
        #add dereddened magnitudes and re-calculate log version of errors    
        deredden_mag = ebv_vec*band_a_ebv[ii]
        cmod_dered =df[f"mag_{band}_cModel"] - deredden_mag
        df[f"cModel_{band}_dered"]=cmod_dered
        invsn = 1./df[f"snr_{band}_cModel"]
        logmagerr = 2.5*np.log10(1.+invsn)
        df[f"magerrlog_{band}_dered"] = logmagerr
        
        #now, replace the non-detections in each band with 99.0 for mag and
        #the 1 sigma limit determined from 1 sigma objects in the same band
        #do this patch by patch to partially account for variable depth/conditions!
        siglow = 0.73
        sighigh = 0.77
        defaultlimmag = 25.8
        for patch in patches:
            goodpatch = True
            sigselect = ((df[f'magerrlog_{band}_dered']>siglow) & (df[f'magerrlog_{band}_dered']<sighigh) \
                         & (df['patch']==patch)\
                         & (np.isfinite(df[f'cModel_{band}_dered'])))
            if np.sum(sigselect)==0:
                siglow = 0.71
                sighigh = 0.79
                sigselect = ((df[f'magerrlog_{band}_dered']>siglow) & (df[f'magerrlog_{band}_dered']<sighigh) \
                             & (df['patch']==patch)\
                             & (np.isfinite(df[f'cModel_{band}_dered'])))
            if np.sum(sigselect)==0:
                print(f"bad data in patch {patch} for band {band}: no 1 sig objects, put in hard coded 25.8 limit")
                goodpatch = False
            if verbose: print(f"{np.sum(sigselect)} total in cut for patch {patch}")
            if goodpatch:
                sigmag = df[f'cModel_{band}_dered'][sigselect]
                limmag = np.median(sigmag)
                defaultlimmag = limmag
            else:
                limmag = 25.8 #hardcoded temp solution
            if verbose: print(f"1 sigma mag for patch {patch} in band {band} is {limmag}")
            #find all NaN and Inf and replace with mag = 99 and magerr = 1 sigma limit
            nondet = ((~np.isfinite(df[f'cModel_{band}_dered']) | (~np.isfinite(df[f'magerrlog_{band}_dered']))) \
                      & (df['patch']==patch))
            df[f'cModel_{band}_dered'][nondet] = 99.0
            df[f'magerrlog_{band}_dered'][nondet] = limmag
            if verbose: print(f"replacing inf and nan for {np.sum(nondet)} {band} band detects for patch {patch}")
    return df