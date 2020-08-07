
import numpy as np
from scipy.interpolate import interp1d

def series_to_N_points(series,N):
    #find the total length of the series
    totalDistance=0
    for s in range(len(series[:,0])-1):
        totalDistance+=((series[s,0]-series[s+1,0])**2+(series[s,1]-series[s+1,1])**2)**0.5
    intervalDistance=totalDistance/(N-1)

    #make the list of points
    newSeries=series[0,:]
    currentS = 0
    currentPoint1=series[currentS,:]
    currentPoint2=series[currentS+1,:]
    for p in range(N-2):
        distanceAccrued = 0
        while distanceAccrued<intervalDistance:
            currentLineDistance=((currentPoint1[0]-currentPoint2[0])**2+(currentPoint1[1]-currentPoint2[1])**2)**0.5
            if currentLineDistance<intervalDistance-distanceAccrued:
                distanceAccrued+=currentLineDistance
                currentS+=1
                currentPoint1 = series[currentS, :]
                currentPoint2 = series[currentS + 1, :]
            else:
                distance=intervalDistance-distanceAccrued
                newX=currentPoint1[0]+(distance/currentLineDistance)*(currentPoint2[0]-currentPoint1[0])
                newY = currentPoint1[1] + (distance / currentLineDistance) * (currentPoint2[1] - currentPoint1[1])
                distanceAccrued=intervalDistance+1
                newSeries=np.vstack([newSeries,np.array([newX,newY])])
                currentPoint1=np.array([newX,newY])
    newSeries = np.vstack([newSeries, series[-1,:]])
    return(newSeries)

def sliding_average(series,N):
    if N%2==0:
        N+=1
    averagedSeries=np.zeros_like(series)
    averagedSeries[:,0]=series[:,0]
    for i in range(len(series)):
        minIndex=int(np.max([0,i-(N-1)/2]))
        maxIndex=int(np.min([len(series),i+(N-1)/2]))
        averagedSeries[i,1]=np.mean(series[minIndex:maxIndex,1])
    return(averagedSeries)

def smooth_timeseries_annually(timeseries):
    step=1/365.25
    time = np.arange(np.min(timeseries[:,0]),np.max(timeseries[:,0]),step)

    set_int = interp1d(timeseries[:,0],timeseries[:,1])
    dense_timeseries = np.hstack([np.reshape(time,(len(time),1)),
                                           np.reshape(set_int(time),(len(time),1))])
    smooth_timeseries = sliding_average(dense_timeseries,365)

    return(smooth_timeseries)

def smooth_timeseries_seasonally(timeseries):
    step=1/365.25
    time = np.arange(np.min(timeseries[:,0]),np.max(timeseries[:,0]),step)

    set_int = interp1d(timeseries[:,0],timeseries[:,1])
    dense_timeseries = np.hstack([np.reshape(time,(len(time),1)),
                                           np.reshape(set_int(time),(len(time),1))])
    smooth_timeseries = sliding_average(dense_timeseries,90)

    return(smooth_timeseries)