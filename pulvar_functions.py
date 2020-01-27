import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.spatial as sp
import math
import os
import scipy.signal as ss
import sys
from tqdm import tqdm_notebook

def plot_profile(profile_data):
    fig, ax = plt.subplots(figsize=(15, 3))
    ax.grid()
    ax.plot(profile_data,'k')
    ax.set_xlim([0,profile_data.shape[0]-1])
    ax.set_ylim([-np.max(profile_data)*0.05,np.max(profile_data)*1.05])
    ax.set_xlabel('Phase Bins')
    ax.set_ylabel('Flux Density (mJy)')
    return ax

def removebaseline(data, outliers):
    # chop profile into 8 parts and check the part with the lowest rms.
    # Remove the mean of that from everything. Remove outliers based on rms.
    nbins = data.shape[0]
    nprofiles = data.shape[1]
    # initialize output array
    baselineremoved = data
    smallestrms = np.zeros(nprofiles)
    smallestmean = np.zeros(nprofiles)
    peak = np.zeros(nprofiles)
    for i in range(nprofiles):
        rms = np.zeros(8)
        mean = np.zeros(8)
        section = int(nbins/8)
        for j in range(8):
            rms[j] = np.std(data[j*section:(j+1)*section,i])
            mean[j] = np.mean(data[j*section:(j+1)*section,i])
        smallestrms[i] = np.min(rms)
        peak[i] = np.max(data[:,i]) # remove low snr not just the std of noise
        baseindex = np.argmin(rms)
        baseline = mean[baseindex]
        smallestmean[i] = baseline
        baselineremoved[:,i] = data[:,i] - baseline
    medianrms = np.median(smallestrms)
    medianpeak = np.median(peak)
    outlierindex = []
    inlierindex = []
    for i in range(nprofiles):
	    if smallestrms[i]/np.max(data[:,i]) > outliers * medianrms/medianpeak:
#        if smallestrms[i] > outliers * medianrms or smallestrms[i] * outliers < medianrms or np.min(baselineremoved[:,i]) < - 5 * smallestrms[i]:
		    outlierindex.append(i)
	    else:
		    inlierindex.append(i)

    ou = np.array(outlierindex)
    inl = np.array(inlierindex)
    

    
    removedprofiles = np.delete(baselineremoved,inl,1)
    baselineoutlierremoved = np.delete(baselineremoved, ou, 1)
    rmsremoved = np.delete(smallestrms, ou)
    print( 'Number of noisy profiles removed: ', ou.shape[0])
    return baselineoutlierremoved, removedprofiles, rmsremoved, ou, inl


def findbrightestprofile(data,rmsdata):
    snr = np.zeros(rmsdata.shape[0])
    for i in range(data.shape[1]):
        snr[i] = np.max(data[:,i])/rmsdata[i]
    brightestindex = np.argmax(snr)
    return brightestindex


def smart_align(template,observations):

    aligned_obs = np.zeros((observations.shape[1],observations.shape[0]))
    aligned_obs_norm = np.zeros((observations.shape[1],observations.shape[0]))

    for n in tqdm_notebook(range(observations.shape[1])):
        #print('Aligning observation',n+1,'of',observations.shape[1])

        obs = observations[:,n]/np.max(observations[:,n])
        bins = template.shape[0]
        first_try = 0
        no_scale_incr = 100
        template_noise_list = []
        obs_noise_list = []
        bins_with_signal_test = []
        list_of_means = []
        list_of_stds = []
        list_list_means = []
        list_list_stds = []
        list_list_no_points = []
        min_arg_list = []
        min_val_list = []
        std_min_val_list = []
        mean_times_points = []


# Correlate to find rough alignment and then start with a fractional offset before rolling the observations past each other                                                               

        # make sure observations don't span the edge
        # rotate template to put peak at 1/4                                                      
        peak_bin = np.argmax(obs)
        initial_shift = int(bins/4)-peak_bin
        obs = np.roll(obs, initial_shift)
                                          
        xcorr = np.correlate(template,obs,"full")
        lag = np.argmax(xcorr)
        obs = np.roll(obs,lag)
        # obs = np.roll(obs,-int(bins/7.0))

        # Break the template into 8 parts and find the rms of each. Then find the smallest. Do with the observation too.                                                                                                            
        for z in range(8):
            template_noise_list.append(np.std(template[z*int(bins/8.0):(z+1)*int(bins/8.0)]))
            obs_noise_list.append(np.std(obs[z*int(bins/8.0):(z+1)*int(bins/8.0)]))

# Find the approximate peaks of template and observation so give an idea of the range over which to scale the observations.                                                               

        temp_peak = np.mean(np.sort(template)[-10:])
        obs_peak = np.mean(np.sort(obs)[-10:])
        rough_scale = temp_peak / obs_peak
        rough_scale_upp = rough_scale * 1.1
        rough_scale_low = rough_scale * 0.9
        scale_incr = (rough_scale_upp - rough_scale_low)/no_scale_incr

        # Keep a copy of the observation in its original state.                                                   
        obs_original = obs[:]

        # Phase shift over all bins.                                                                              
        for roll in range(int(bins/3.5)):
#            if (roll+1)%100 == 0:
#                print( 'Bin',roll+1,'out of',int(bins/3.5)
            closest_to_one =1e10
            bins_with_signal = []
            bins_with_signal_test = []
            list_mean_each_scale_shift = []
            list_std_each_scale_shift = []
            list_points_each_scale_shift = []
        # No shift needed for first try.                                                                          
            if roll != 0:
                obs = np.roll(obs,1)
# If the level is too low in either template or observation, don't include the bin in further analysis.    
            for r in range(bins):
                #print( r,obs[r],obs_peak,np.min(obs_noise_list),template[r],temp_peak,np.min(template_noise_list) 
                if obs[r] > obs_peak/3. and template[r] > temp_peak/3.:
                    bins_with_signal.append(r)
                    bins_with_signal_test.append(1)
                else:
                    bins_with_signal_test.append(0)

        # For each roll, only proceed if there are more than 20 bins that have signal in them.                    
            if len(bins_with_signal) >= 10.0:
        # Loop over each of the 100 scale attempts to find which is the best fit.                                 
                for s in range(no_scale_incr):
                    if s == 0:
                        first_scale_val = rough_scale_low+s*scale_incr
                    if s == no_scale_incr-1:
                        last_scale_val = rough_scale_low+s*scale_incr
                    diff = []
                    escape = 0
                    scaled_obs=obs*(rough_scale_low+s*scale_incr)
         # Loop over all the bins with signal and find the absolute difference between template and observation.  
                    for each in bins_with_signal:
                        diff.append(abs(scaled_obs[each] - template[each]))
        # Save this difference list before outliers are removed.                                                  
                    orig_diff = diff[:]
        # Remove outliers (over 2 sigma) and re-evaluate the mean. If mean doesn't change much, exit the loop. Record the last set of data that had outliers removed.                                                               
                    while escape == 0:
                        diff_mean = np.mean(diff)
                        diff_std = np.std(diff)
                        outlier_list = []
                        for y in range(len(diff)):
                            if abs(diff[y]-diff_mean) > 2*diff_std:
                                outlier_list.append(y)
                        latest_diff_list = diff[:]
                        for index in sorted(outlier_list, reverse=True):
                            del diff[index]

                        if np.mean(diff) == 0:
                            escape = 1
                            diff = latest_diff_list[:]
                        else:
                            if np.mean(diff)/diff_mean < 1.001 and np.mean(diff)/diff_mean > 0.999 and first_try == 1:
                                escape = 1
                                diff = latest_diff_list[:]
                        first_try = 1
        # In lists - For any phase, record the mean and std and number of data points after all outliers removed at each scale attempt.                                                                                             
                    list_mean_each_scale_shift.append(abs(np.mean(diff)))
                    list_std_each_scale_shift.append(np.std(diff))
                    list_points_each_scale_shift.append(len(diff))

        # Make a list containing the above lists. 1 per phase shift.                                              
                list_list_means.append(list_mean_each_scale_shift)
                list_list_stds.append(list_std_each_scale_shift)
                list_list_no_points.append(list_points_each_scale_shift)

            else:
                # If the number of bins with signal is not high enough, just put 1s into the list of lists. We will find minimum later, and 1 is >>.                                                                                
                list_list_means.append([1]*no_scale_incr)
                list_list_stds.append([1]*no_scale_incr)
                list_list_no_points.append([1]*no_scale_incr)

        # Calculate the mean / number of points. This should be minimised to find the best fit.                   
        for h in range(len(list_list_means)):
            for y in range(len(list_list_means[0])):
                mean_times_points.append(list_list_means[h][y]/list_list_no_points[h][y])

        min_arg_final = np.argmin(mean_times_points)
        the_scale = min_arg_final%no_scale_incr
        the_roll = min_arg_final/no_scale_incr
        min_val_final = np.min(mean_times_points)
        std_min_val_final = list_list_stds[int(min_arg_final/no_scale_incr)][int(min_arg_final%no_scale_incr)]

# Return the aligned and scaled observations.                                                                     

        aligned_obs_norm[n,:] = np.roll(obs_original*(rough_scale_low+the_scale*scale_incr),int(the_roll))
#        obs = observations[:,n]/np.max(observations[:,n])                                                        


        aligned_obs[n,:] = np.roll(obs_original*np.max(observations[:,n]),int(the_roll))                         
    aligned_obs = np.transpose(aligned_obs)
    aligned_obs_norm = np.transpose(aligned_obs_norm)
    
    return aligned_obs_norm, aligned_obs

def save_profiles(pulsar, data, mjd, allbins,directory,ylabel='Flux Density (mJy)',template=None,no_y_lim=False):
    nbins = data.shape[0]
    nprofiles = data.shape[1]
    xaxis = np.linspace(0,float(nbins)/allbins,nbins)
    xlocs = np.linspace(0,nbins-1,11)
    xticklabels = []
    for i in xlocs:
            xticklabels.append(np.round(xaxis[int(i)],2))
    if not (os.path.exists('./{0}/{1}'.format(pulsar,directory))):
        os.mkdir('./{0}/{1}'.format(pulsar,directory))
    for i in tqdm_notebook(range(nprofiles)):
        fig, ax = plt.subplots(figsize=(15, 4))
        ax.grid()
        ax.plot(data[:,i],'k',alpha=0.7,linewidth=2)
        if template is None:
            pass
        else:
            ax.plot(np.roll(template,0),'r--',linewidth=2)
        #plt.yticks([])
        ax.set_xlim(0,data.shape[0]-1)
        if no_y_lim is False:
            ax.set_ylim([-np.max(data)*0.05,np.max(data)*1.05])
        ax.set_ylabel(r'Flux Density (mJy)',fontsize=14)
        ax.set_xlabel(r'Phase Bins',fontsize=14)
        ax.set_xticks(xlocs,xticklabels)
        fig.tight_layout()
        plt.savefig('./{0}/{3}/{1}_{2}.png' .format(pulsar,int(math.floor(mjd[i])),i,directory))
        plt.clf()
        plt.close(fig)

def jd_to_date(jd):
    """
    Convert Julian Day to date.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)

    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
        
    return year, month, day

def variability_map_plot(obs_mjd,gp_data,flux_data,template,window_begin,window_end,total_bins):
    days_from_start = obs_mjd-obs_mjd[0]
    max_gp_data = np.max(gp_data)
    min_gp_data = np.min(gp_data)
    limitdifference = np.max((max_gp_data, np.abs(min_gp_data)))
    my_max = limitdifference
    my_min = -limitdifference

    fig = plt.figure(figsize=(20,8))
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    ax1 = plt.subplot2grid((3,10),(0,0),colspan = 9)
    ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False, top=True, labeltop=True)
    ax1.plot(obs_mjd,flux_data,'kx',markersize=10)
    ax1.set_xlim(obs_mjd[0],gp_data.shape[1]+obs_mjd[0])
    ax1.set_ylim(-0.15*np.max(flux_data),1.15*np.max(flux_data))
    plt.grid()
    fig.text(0.05, 0.765, 'Peak Flux\nDensity (mJy)', ha='center', va='center', rotation='vertical', size=18)
    fig.text(0.47, 0.98, 'Modified Julian Date', ha='center', va='center', rotation='horizontal', size=18)
    
    ax2 = plt.subplot2grid((3,10),(1,0), rowspan=2, colspan = 9)
#plt.imshow(gp_data, origin='upper', interpolation='nearest', aspect='auto', cmap=plt.cm.Greys);
    plt.imshow(gp_data, aspect='auto', cmap = "RdBu_r", vmin = my_min, vmax = my_max)
    cbaxes = fig.add_axes([0.87, 0.649, 0.03, 0.232])
    cb = plt.colorbar(cax = cbaxes,orientation="vertical")
    cb.update_ticks()
    cbaxes.tick_params(axis='y', which='both', left=True, labelleft=True, right=False, labelright=False)
    cbaxes.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
#cb.tick_params(labelsize=12)
    for obs in range(days_from_start.shape[0]):
        ax2.vlines(days_from_start[obs],0,window_end-window_begin,linestyles='solid',alpha=0.1)
    ax2.set_xlim(0,gp_data.shape[1]);
    ax2.set_ylim(window_end-window_begin-1,0)
    date_labels = np.linspace(0,gp_data.shape[1],5)+obs_mjd[0]
    date_label_str = []
    for d in range(date_labels.shape[0]):
        date = (jd_to_date(date_labels[d]+2400000.5))
        date_label_str.append(str(int(date[0]))+'/'+str(int(date[1]))+'/'+str(int(date[2])))
    ax2.set_xticklabels(date_label_str)
    ax2.set_xticks(np.linspace(0,gp_data.shape[1],5))
    y_fractions = np.linspace(0,float(gp_data.shape[0])/total_bins,gp_data.shape[0])
    mylocs = np.linspace(0,gp_data.shape[0]-1,5)
    myticklabels = []
    for y in mylocs:
        myticklabels.append(np.round(y_fractions[int(y)],2))
    ax2.set_yticks(mylocs)
    ax2.set_yticklabels(myticklabels)
    fig.text(0.05, 0.380, 'Fraction of\nPulse Period', ha='center', va='center', rotation='vertical', size=18)
    fig.text(0.475, 0.025, 'Calendar Date', ha='center', va='center', rotation='horizontal', size=18)
    
    ax3 = plt.subplot2grid((3,10),(1,9),rowspan = 2, colspan = 1)
    ax3.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax3.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax3.plot(template[window_begin:window_end],np.arange(window_end-window_begin-1,-1,-1),'k')
    ax3.set_ylim(0,window_end-window_begin)
    ax3.set_xlim(-np.max(template[window_begin:window_end])*0.05,np.max(template[window_begin:window_end])*1.05)
    
    fig.subplots_adjust(hspace=0.0)
    fig.subplots_adjust(wspace=0.0)
    

def aligndata(baselineremoved, brightest, pulsar):
        nbins = baselineremoved.shape[0]
        nprofiles = baselineremoved.shape[1]
        template = baselineremoved[:,brightest]
    # rotate template to put peak at 1/4
        peakbin = np.argmax(template)
        fixedlag = int(nbins/4)-peakbin
    #fixedlag = 0
        aligned = np.zeros((nbins,nprofiles))
        aligned_temp = np.zeros((nbins,nprofiles))
        aligned_temp2 = np.zeros((nbins,nprofiles))	
        newtemplate = np.roll(template, fixedlag)
        template = newtemplate
        #plt.plot(newtemplate)
        #plt.savefig('./{0}/{0}_brightest.png' .format(pulsar))
        #plt.clf()
        for i in range(nprofiles):
            xcorr = np.correlate(template,baselineremoved[:,i],"full")
            lag = np.argmax(xcorr)
            aligned[:,i] = np.roll(baselineremoved[:,i],lag)
        template = np.median(aligned,1)
    # repeat with better template now and shift peak to 1/4 of the profile
        peakbin = np.argmax(template)
        fixedlag = int(nbins/4)-peakbin
    #fixedlag = 0
        double = np.zeros(2*nbins)
        for i in range(nprofiles):
            double[0:nbins] = baselineremoved[:,i]
            double[nbins:2*nbins] = baselineremoved[:,i]
        #xcorr = np.correlate(template,baselineremoved[:,i],"full")
            xcorr = np.correlate(template,double,'full')
            lag = np.argmax(xcorr) + fixedlag
            lag_temp = np.argmax(xcorr) + fixedlag - 1
            aligned[:,i] = np.roll(baselineremoved[:,i],lag)
            aligned_temp[:,i] = np.roll(baselineremoved[:,i],lag_temp)
            newtemplate = np.median(aligned,1)
            run_tot = 0
        
        return np.array(aligned), np.array(newtemplate)
