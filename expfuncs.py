from TimeTagger import *
from tqdm import tqdm, tqdm_notebook, tnrange
import numpy as np
import TimeTagger as TT
import pandas as pd
import scipy.io as sio
import sys
import qcodes as qc
from swabian import SwabianLaser
from pulse_streamer_grpc import PulseStreamer
from pulse_streamer_grpc import Start, Serial, Mode
from Sequence import Sequence
from Sequence import random_ana, random_digi
import numpy as np
from TimeTagger import *
from plotly import graph_objs as go
from qcodes.instrument_drivers.rohde_schwarz import SGS100A
from time import sleep
from qcodes.dataset.measurements import Measurement
from qcodes import ParamSpec
import qcodes.dataset.experiment_container as exc
import datetime
import plotly.graph_objs as go


def laserCW():
    """
    sets laser to on permanently via pulsestreamer
    :return:
    None
    """
    pulser.constant((1, [2], 0, 0))
    return None


def laserTrig():
    """
    sets laser to pulsed mode by setting pulsestreamer to off
    :return:
    """
    pulser.constant((1, [], 0, 0))
    return None


def fluorescentMeasurement(pulser, tagger, measTime, xPoints, yPoints, dataSetName='Fluorescence', origin=[0, 0],
                           scanRange=10e-6):
    """
    Performs a fluorescence measurement over a specified grid
    :param pulser: Pulsestreamer being used
    :param tagger: Timetagger being used
    :param measTime: How long to record at each pixel
    :param xPoints: Number of pixels on x axis
    :param yPoints: Number of pixels on y axis
    :param dataSetName: Dataset name for qcodes
    :param origin: x,y origin for saving exact scan position
    :param scanRange: scan range for correct x,y index values
    :return:
    """
    data_set = exc.new_data_set(name=dataSetName)
    data_set.add_parameter(ParamSpec(name='x', unit='um', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='y', unit='um', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Counts', unit='', paramtype='array'))
    xData = np.linspace(origin[0], origin[0] + scanRange, xPoints)
    yData = np.linspace(origin[1], origin[1] + scanRange, yPoints)
    countsData = np.zeros([yPoints, xPoints])
    data_set.add_result({'x': xData, 'y': yData, 'Counts': countsData})
    gateStart = 0
    #     gateFinish = 4
    #     APD_gate = 6
    Laser = 2
    Atto = 7
    #     recSeq.insert(0, (int(100+measTime/2)))
    #     recSeq[0] = (int(100+measTime/2), 1)
    # Define pulse sequence and send it to pulsestreamer
    pulser.reset()
    print('pulser_reset')
    seq = Sequence()
    #     seq.setDigitalChannel(laser, [(int(10), 0), (int(2*measTime), 1), (int(10), 0)])
    #     seq.setDigitalChannel(laser, recSeq)
    seq.setDigitalChannel(gateStart, [(int(0 + measTime / 2), 0), (int(measTime), 1), (int(measTime / 2), 0)])
    seq.setDigitalChannel(Laser, [(int(0 + measTime / 2), 0), (int(measTime), 1), (int(measTime / 2), 0)])
    #     seq.setDigitalChannel(APD_gate, [(int(0+measTime/2), 0), (int(measTime), 1), (int(100), 0)])
    #     seq.setDigitalChannel(gateFinish, [(int(0+3*measTime/2), 0), (int(10e6), 1), (int(measTime/2), 0)])
    seq.setDigitalChannel(Atto, [(int(2 * measTime), 0), (int(1e4), 1), (int(10e3), 0), (int(1e4), 1)])
    finalSeq = seq.getSequence()
    print(finalSeq)

    pulser.setTrigger(start=Start.HARDWARE_RISING_AND_FALLING, mode=Mode.NORMAL)
    pulser.stream(finalSeq, n_runs=1, final='CONSTANT_ZERO')

    # Set up timetagger for gated measurement
    countbetweenmarkers = CountBetweenMarkers(tagger, 1, 7, -7, xPoints * yPoints)
    while np.count_nonzero(countbetweenmarkers.getBinWidths()) < xPoints * yPoints:
        sleep(0.2)
        data_set.modify_result(0, {'Counts': countbetweenmarkers.getData().reshape(yPoints, xPoints)})
    return countbetweenmarkers


def g2Measurement(pulser, tagger, dataSetName='g2', time=100, binwidth=100, nbins=200):
    """
    Performs a g2 measurement and saves it to a qcodes database
    :param pulser: Pulsestreamer being used
    :param tagger: Timetagger being used
    :param dataSetName: name for the qcodes dataset, default 'g2'
    :param time: How long to perform the g2 measurement for in seconds
    :param binwidth: binwidth in picoseconds
    :param nbins: Number of bins to use
    :return: None
    """
    data_set = exc.new_data_set(name=dataSetName)
    data_set.add_parameter(ParamSpec(name='Time', unit='ns', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Counts', unit='ns', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Norm_Counts', unit='ns', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Count_Rate', unit='ns', paramtype='array'))

    g2 = Correlation(tagger, 1, 2, binwidth=binwidth, n_bins=nbins)
    count_rate = Countrate(tagger, [1, 2])
    count_rate_val = count_rate.getData()
    g2.stop()
    g2.startFor(int(time * 1e12))
    data_set.add_result({'Time': np.squeeze(g2.getIndex()), 'Counts': np.squeeze(g2.getData()),
                         'Norm_Counts': np.squeeze(g2.getData()) / (np.max(g2.getData())),
                         'Count_Rate': count_rate_val})
    i = 0
    while g2.isRunning() == True:
        i += 1
        count_rate_val = (count_rate_val * (i - 1) + count_rate.getData()) / i
        sleep(2)
        data = np.squeeze(g2.getData())
        norm = np.mean(data[:10])
        if norm == 0.:
            norm = 1
        dataNorm = data / norm
        data_set.modify_result(0, {'Counts': np.squeeze(data), 'Norm_Counts': np.squeeze(dataNorm),
                                   'Count_Rate': count_rate_val})


def counting(readtime=50e6,samples=300,name='counting'):
    data_set = exc.new_data_set(name = name)
    data_set.add_parameter(ParamSpec(name = 'Time', unit = '', paramtype='array'))
    data_set.add_parameter(ParamSpec(name = 'Rebased_Counts', unit = '', paramtype='array'))
    seq = Sequence()
    gateStart = 0
    laserNum = 2
    waitLen = 1e3
    lasPulse = int(readtime)
    countLen = int(readtime)
    seq.setDigitalChannel(laserNum, [(int(waitLen), 0), (int(lasPulse), 1), (int(waitLen), 0)])
    seq.setDigitalChannel(gateStart, [(int(waitLen), 0), (int(countLen), 1), (int(waitLen), 0)])
    finalSequence = seq.getSequence()
    pulser.reset()
    pulser.setTrigger(start= Start.SOFTWARE, mode=Mode.NORMAL)
    pulser.stream(finalSequence, n_runs = 1, final = 'CONSTANT_ZERO')
    counter = CountBetweenMarkers(tagger, 1, 7, -7, n_values=1)
    totalData = np.zeros([samples])
    time = np.linspace(readtime,readtime*samples,samples)
    sleep(1)
    data_set.add_result({'Time': time, 'Rebased_Counts': totalData})
    for i in range(samples):
        counter.clear()
        pulser.startNow()
        while np.count_nonzero(counter.getBinWidths()) < len(counter.getBinWidths()):
            pass
        totalData[i] = 1e6*counter.getData()/countLen
        counter.clear()
        data_set.modify_result(0, {'Rebased_Counts': totalData})
    print('Done')
    return None


def constantOut(digital = [], analog = [], totalTime = 10, onTime = 1e6, offTime = 1e6):
    seqChan = Sequence()
    if np.size(digital) != 0:
        for i in range(np.size(digital)):
            chanD = digital[i]
            seqChan.setDigitalChannel(chanD, [(int(onTime), 1),(int(offTime), 0)])
    if np.size(analog) != 0:
        for i in range(np.size(analog)):
            chanA = analog[i]
            seqChan.setAnalogChannel(chanA, [(int(onTime), 1),(int(offTime), 0)])
    chanSeq = seqChan.getSequence()
    pulser.reset()
    pulser.setTrigger(start=Start.SOFTWARE)
    pulser.stream(chanSeq, n_runs = int(totalTime*1e9/(int(onTime) + int(offTime))), final = 'CONSTANT_ZERO')
    pulser.startNow()
    sleep(int(totalTime+2))
    print('Done')
    return None


def cwODMR(start, end, steps, samples, averages, mwPulse=10e6, lasPulse=10e6, waitLen=2e6, name='cwODMR', cwMW=0):
    data_set = exc.new_data_set(name=name)
    data_set.add_parameter(ParamSpec(name='Frequency', unit='Hz', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Counts', unit='', paramtype='array'))
    seq = Sequence()
    gateStart = 0
    sourceNum = 1
    Ichan = 0
    laserNum = 2
    countLen = lasPulse - 5e6
    frequencies = np.linspace(start, end, steps)
    #    seq.setDigitalChannel(sourceNum, [(int(waitLen), 0), (int(lasPulse), 1), (int(waitLen), 0)])
    seq.setAnalogChannel(Ichan, [(int(waitLen), 0), (int(lasPulse), 1), (int(waitLen), 0)])
    seq.setDigitalChannel(laserNum, [(int(waitLen), 0), (int(lasPulse), 1), (int(waitLen), 0)])
    seq.setDigitalChannel(gateStart, [(int(waitLen + 5e6), 0), (int(countLen), 1), (int(waitLen), 0)])
    #     seq.setDigitalChannel(gateFinish, [(int(waitLen+5e6+countLen), 0), (int(10e3), 1), (int(waitLen-10e3),0)])
    finalSequence = seq.getSequence()
    #    print(finalSequence)
    print('Applying cw-MW for ' + str(cwMW) + ' s')
    constantOut(analog=[0], totalTime=cwMW)
    print('Starting the sequence')
    pulser.reset()
    pulser.setTrigger(start=Start.SOFTWARE, mode=Mode.NORMAL)
    totalData = np.zeros([len(frequencies), samples])
    binWidthData = np.zeros([len(frequencies), samples])
    pulser.stream(finalSequence, n_runs=int(samples), final='CONSTANT_ZERO')
    counter = CountBetweenMarkers(tagger, 1, 7, -7, n_values=samples)
    sleep(1)
    counts_data = np.zeros([len(frequencies), 1])
    data_set.add_result({'Frequency': frequencies, 'Counts': counts_data})
    data_set.add_metadata('sourcePower', source.power())
    #     with meas.run() as datasaver:
    for av in tnrange(averages):
        #         print(av)
        countsList = np.zeros([len(frequencies), 1])
        #         print(counts_data)
        #         print(countsList)
        for i in range(len(frequencies)):
            #         print(frequencies[i])
            source.frequency(frequencies[i])
            #             print(source.frequency())
            counter.clear()
            pulser.startNow()
            #         sleep(0.2)
            while np.count_nonzero(counter.getBinWidths()) < len(counter.getBinWidths()):
                pass
            #         print(counter.getBinWidths())
            totalData[i] = counter.getData()
            binWidthData[i] = counter.getBinWidths()
            countsList[i] = np.mean(counter.getData())
            #             print(countsList)
            if av == 0:
                counts_data[i] = countsList[i]
            #                 print(counts_data)
            else:
                counts_data[i] = (av * counts_data[i] + countsList[i]) / (av + 1)
            data_set.modify_result(0, {'Counts': counts_data})

            #                 datasaver.add_result(('Frequency', source.frequency()), ('Counts', np.mean(counter.getData())))
            counter.clear()
    #             print('counts_data is: \n', counts_data)
    return counts_data, frequencies


def countSweep(start, end, steps, samples, lasInit = 5e3, readTime = 1e2, tWait = 20e3):
#     gap = [(int(1e6),0)]
    lasChan = 2
    testChan = 6
    testChan2 = 5
    mwChan = 1
    gateChan = 0
#     stopChan = 2
    seqLas = []
    seqMW = []
#     seqStart = []
#     seqStop = []
    seqRead = []
    readTime = int(readTime)
    lengths = np.linspace(start, end, steps)
    tWait = int(tWait)
    for length in lengths:
        for i in range(samples):
            seqLas += [(tWait, 0), (int(lasInit), 1), (tWait, 0)]
            seqRead += [(int(tWait+length), 0), (readTime, 1), (int(lasInit-length-readTime+tWait), 0)]
#             seqStop += [(int(tWait + length + readTime), 0), (500, 1), (int(lasInit - length +tWait - 500), 0)]
#             totLas = 0
#             for x in seqRead:
#                 totLas += int(x[0])
#             print(totLas)
    seq = Sequence()
    seq.setDigitalChannel(lasChan, seqLas)
#     seq.setDigitalChannel(mwChan, seqMW)
#     seq.setDigitalChannel(testChan, seqLas)
#     seq.setDigitalChannel(testChan2, seqRead)
    seq.setDigitalChannel(gateChan, seqRead)
#     seq.setDigitalChannel(stopChan, seqStop)
    return seq.getSequence()


def sweepRead(start, end, steps, samples, lasInit = 5e3, lasRead = 5e3, lasRef = 5e3, piPulse = 1.4e2):
#     gap = [(int(1e6),0)]
    lasChan = 2
    mwChan = 1
    actChan = 0
    refChan = 4
    seqLas = []
    seqMW = []
    seqAct = []
    seqRef = []
#     readTime = int(readTime)
    lengths = np.linspace(start, end, steps)
    tWait = int(5e2)
    delay = int(2.2e2)
    for length in lengths:
        for i in range(samples):
            seqLas += [(int(tWait), 0), (int(lasInit), 1), (int(tWait + piPulse + tWait), 0), (int(lasRead), 1), (int(tWait), 0), (int(lasRef), 1), (int(tWait), 0)]
            seqMW += [(int(tWait+lasInit+tWait), 0), (int(piPulse), 1), (int(tWait+lasRead+tWait+lasRef+tWait), 0)]
            seqAct += [(int(tWait+lasInit+tWait+piPulse+tWait+delay), 0), (int(length), 1), (int(lasRead-length+tWait), 0), (int(lasRef+tWait-delay), 0)]
            seqRef += [(int(tWait+lasInit+tWait+piPulse+tWait+delay+lasRead+tWait), 0), (int(length), 1), (int(lasRef-length - delay + tWait), 0)]
#             totLas = 0
#             for x in seqStop:
#                 totLas += int(x[0])
#             print(totLas)
    seq = Sequence()
    seq.setDigitalChannel(lasChan, seqLas)
    seq.setDigitalChannel(mwChan, seqMW)
    seq.setDigitalChannel(actChan, seqAct)
    seq.setDigitalChannel(refChan, seqRef)
    return seq.getSequence()


def sweepDelay(start, end, steps, samples, lasInit = 5e3, lasRead = 5e3, lasRef = 5e3, piPulse = 1.4e2, readTime = 200):
#     gap = [(int(1e6),0)]
    lasChan = 2
    mwChan = 1
    actChan = 0
    refChan = 4
    seqLas = []
    seqMW = []
    seqAct = []
    seqRef = []
    readTime = int(readTime)
    lengths = np.linspace(start, end, steps)
    tWait = int(5e2)
    delay = int(2e2)
    for length in lengths:
        for i in range(samples):
            seqLas += [(int(tWait), 0), (int(lasInit), 1), (int(tWait + piPulse + tWait), 0), (int(lasRead), 1), (int(tWait), 0), (int(lasRef), 1), (int(tWait), 0)]
            seqMW += [(int(tWait+lasInit+tWait), 0), (int(piPulse), 1), (int(tWait+lasRead+tWait+lasRef+tWait), 0)]
            seqAct += [(int(tWait+lasInit+tWait+piPulse+tWait+length), 0), (readTime, 1), (int(lasRead-readTime+tWait), 0), (int(lasRef+tWait-length), 0)]
            seqRef += [(int(tWait+lasInit+tWait+piPulse+tWait+length+lasRead+tWait), 0), (readTime, 1), (int(lasRef-readTime - length + tWait), 0)]
#             totLas = 0
#             for x in seqStop:
#                 totLas += int(x[0])
#             print(totLas)
    seq = Sequence()
    seq.setDigitalChannel(lasChan, seqLas)
    seq.setDigitalChannel(mwChan, seqMW)
    seq.setDigitalChannel(actChan, seqAct)
    seq.setDigitalChannel(refChan, seqRef)
    return seq.getSequence()


def seqInit(lasChan=2, lasInit=3e3):
    seqLas = []
    tWait = int(5e2)
    seqLas = [(int(tWait), 0), (int(lasInit), 1), (int(tWait), 0)]
    seq_init = Sequence()
    seq_init.setDigitalChannel(lasChan, seqLas)
    return seq_init.getSequence()


def seqRead(lasChan=2, countReadChan=0, lasRead=3e3, readTime=3e2):
    seqLas = []
    seqAct = []
    tWait = int(5e2)
    delay = int(2.2e2)
    seqLas = [(int(tWait), 0), (int(lasRead), 1), (int(tWait), 0)]
    seqAct = [(int(tWait+delay), 0), (int(readTime), 1), (int(lasRead-readTime-delay), 0), (int(tWait), 0)]
    seq_read = Sequence()
    seq_read.setDigitalChannel(lasChan, seqLas)
    seq_read.setDigitalChannel(countReadChan, seqAct)
    return seq_read.getSequence()


def seqRabi(start, stop, steps, samples=1, lasChan=2, countReadChan=0, countRefChan=4, lasInit=3e3, lasRead=3e3, readTime=3e2):
    seq1 = seqInit(lasChan=lasChan, lasInit=lasInit)
    seq2 = seqRead(lasChan=lasChan, countReadChan=countReadChan, lasRead=lasRead, readTime=readTime)
    seq3 = seqRead(lasChan=lasChan, countReadChan=countRefChan, lasRead=lasRead, readTime=readTime)
    Ichan = 0
    Qchan = 1
    seqI = []
    seq_operate = Sequence()
    seq_final = []
    lengths = np.linspace(start, stop, steps)
    for length in lengths:
        seqI = [(int(1e1), 0), (int(length), 1), (int(1e1), 0)]
        seq_operate.setAnalogChannel(Ichan, seqI)
        seq_test = seq_operate.getSequence()
        seq_final += (seq1 + seq_test + seq2 + seq3)*samples
    return seq_final


def rabiAnalog(start, end, steps, samples=1, lasInit = 5e3, lasRead = 5e3, lasRef = 5e3, readTime = 3e2):
    lasChan = 2
    mwChan = 1
    Ichan = 0
    actChan = 0
    refChan = 4
    seqLas = []
    seqMW = []
    seqAct = []
    seqRef = []
    readTime = int(readTime)
    lengths = np.linspace(start, end, steps)
    tWait = int(5e2)
    delay = int(2.2e2)
    for length in lengths:
        for i in range(samples):
            seqLas += [(int(tWait), 0), (int(lasInit), 1), (int(tWait + length + tWait), 0), (int(lasRead), 1), (int(tWait), 0), (int(lasRef), 1), (int(tWait), 0)]
            seqMW += [(int(tWait+lasInit+tWait), 0), (int(length), 1), (int(tWait+lasRead+tWait+lasRef+tWait), 0)]
            seqAct += [(int(tWait+lasInit+tWait+length+tWait+delay), 0), (int(readTime), 1), (int(lasRead-readTime+tWait), 0), (int(lasRef+tWait-delay), 0)]
            seqRef += [(int(tWait+lasInit+tWait+length+tWait+delay+lasRead+tWait), 0), (int(readTime), 1), (int(lasRef-readTime - delay + tWait), 0)]

    seq = Sequence()
    seq.setDigitalChannel(lasChan, seqLas)
    seq.setAnalogChannel(Ichan, seqMW)
    seq.setDigitalChannel(actChan, seqAct)
    seq.setDigitalChannel(refChan, seqRef)
    return seq.getSequence()


def seqOdmrPulse(piPulse, samples=1, lasChan=2, countReadChan=0, countRefChan=4, lasInit=5e3, lasRead=5e3, readTime=5e2):
    seq1 = seqInit(lasChan=lasChan, lasInit=lasInit)
    seq2 = seqRead(lasChan=lasChan, countReadChan=countReadChan, lasRead=lasRead, readTime=readTime)
    seq3 = seqRead(lasChan=lasChan, countReadChan=countRefChan, lasRead=lasRead, readTime=readTime)
    Ichan = 0
    Qchan = 1
    seqI = []
    seq_operate = Sequence()
    seq_final = []
    seqI = [(int(1e1), 0), (int(piPulse), 1), (int(1e1), 0)]
    seq_operate.setAnalogChannel(Ichan, seqI)
    seq_test = seq_operate.getSequence()
    seq_final += (seq1 + seq_test + seq2 + seq3)*samples
    return seq_final


def seqRamsey(start, stop, steps, piPulse, samples=1, lasChan=2, countReadChan=0, countRefChan=4, lasInit=5e3,
              lasRead=5e3, readTime=5e2, readState=0):
    seq1 = seqInit(lasChan=lasChan, lasInit=lasInit)
    seq2 = seqRead(lasChan=lasChan, countReadChan=countReadChan, lasRead=lasRead, readTime=readTime)
    seq3 = seqRead(lasChan=lasChan, countReadChan=countRefChan, lasRead=lasRead, readTime=readTime)
    Ichan = 0
    Qchan = 1
    seqI0 = []
    seqI1 = []
    seq_final = []
    tWait = 0
    halfpiPulse = int(piPulse / 2)
    seq_operate0 = Sequence()
    seq_operate1 = Sequence()
    lengths = np.linspace(start, stop, steps)
    if lengths[0] == 0:
        lengths[0] = 100
    for length in lengths:
        seqI0 = [(int(1e1), 0), (halfpiPulse, 1), (int(tWait + length), 0), (halfpiPulse, 1), (int(1e1), 0)]
        seqI1 = [(int(1e1), 0), (halfpiPulse, 1), (int(tWait + length), 0), (halfpiPulse, -1), (int(1e1), 0)]

        seq_operate0.setAnalogChannel(Ichan, seqI0)
        seq_test0 = seq_operate0.getSequence()

        seq_operate1.setAnalogChannel(Ichan, seqI1)
        seq_test1 = seq_operate1.getSequence()

        if readState == 0:
            seq_final += (seq1 + seq_test0 + seq2 + seq3) * samples
        elif readState == 1:
            seq_final += (seq1 + seq_test1 + seq2 + seq3) * samples
        else:
            seq_final += (seq1 + seq_test0 + seq2 + seq1 + seq_test1 + seq2) * samples
    return seq_final


def seqSpinecho(start, stop, steps, piPulse, samples=1, lasChan=2, countReadChan=0, countRefChan=4, lasInit=5e3,
                lasRead=5e3, readTime=3e2, readState=0):
    "In this version one can select the readout state of NV to be 0, 1, or both"
    seq1 = seqInit(lasChan=lasChan, lasInit=lasInit)
    seq2 = seqRead(lasChan=lasChan, countReadChan=countReadChan, lasRead=lasRead, readTime=readTime)
    seq3 = seqRead(lasChan=lasChan, countReadChan=countRefChan, lasRead=lasRead, readTime=readTime)
    Ichan = 0
    Qchan = 1
    seqI0 = []
    seqI1 = []
    seq_final = []
    tWait = 0
    halfpiPulse = int(piPulse / 2)
    seq_operate0 = Sequence()
    seq_operate1 = Sequence()
    lengths = np.linspace(start, stop, steps)
    if lengths[0] == 0:
        lengths[0] = 100
    for length in lengths:
        seqI0 = [(int(1e1), 0), (halfpiPulse, 1), (int(tWait + length), 0), (int(piPulse), 1), (int(tWait + length), 0),
                 (halfpiPulse, 1), (int(1e1), 0)]
        seqI1 = [(int(1e1), 0), (halfpiPulse, 1), (int(tWait + length), 0), (int(piPulse), 1), (int(tWait + length), 0),
                 (halfpiPulse, -1), (int(1e1), 0)]

        seq_operate0.setAnalogChannel(Ichan, seqI0)
        seq_test0 = seq_operate0.getSequence()

        seq_operate1.setAnalogChannel(Ichan, seqI1)
        seq_test1 = seq_operate1.getSequence()

        if readState == 0:
            seq_final += (seq1 + seq_test0 + seq2 + seq3) * samples
        elif readState == 1:
            seq_final += (seq1 + seq_test1 + seq2 + seq3) * samples
        else:
            seq_final += (seq1 + seq_test0 + seq2 + seq1 + seq_test1 + seq3) * samples
    return seq_final


def spinechoAnalog(start, end, steps, piPulse, samples=1, lasInit=5e3, lasRead=5e3, lasRef=5e3, readState=0,
                   readTime=300, sequence=None):
    "In this version one can select the readout state of NV to be 0, 1, or both"
    lasChan = 2
    mwChan = 1
    Ichan = 0
    actChan = 0
    refChan = 4
    seqLas = []
    seqMW = []
    seqAct = []
    seqRef = []
    readTime = int(readTime)
    lengths = np.linspace(start, end, steps)
    if sequence is not None:
        lengths = sequence
    halfpiPulse = int(piPulse / 2)
    tWait = int(2e2)
    tWait2 = int(50)
    delay = int(2.2e2)
    for length in lengths:
        for i in range(samples):
            if readState == 0:
                seqLas += [(int(tWait), 0), (int(lasInit), 1), (
                int(tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait), 0),
                           (int(lasRead), 1), (int(tWait), 0), (int(lasRef), 1), (int(tWait), 0)]
                seqMW += [(int(tWait + lasInit + tWait), 0), (int(halfpiPulse), 1), (int(tWait2 + length), 0),
                          (int(piPulse), 1), (int(tWait2 + length), 0), (int(halfpiPulse), 1),
                          (int(tWait + lasRead + tWait + lasRef + tWait), 0)]
                seqAct += [(int(
                    tWait + lasInit + tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait + delay),
                            0), (readTime, 1), (int(lasRead - readTime - delay + tWait + lasRef + tWait), 0)]
                seqRef += [(int(
                    tWait + lasInit + tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait + lasRead + tWait + delay),
                            0), (readTime, 1), (int(lasRef - readTime - delay + tWait), 0)]
            elif readState == 1:
                seqLas += [(int(tWait), 0), (int(lasInit), 1), (
                int(tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait), 0),
                           (int(lasRead), 1), (int(tWait), 0), (int(lasRef), 1), (int(tWait), 0)]
                seqMW += [(int(tWait + lasInit + tWait), 0), (int(halfpiPulse), 1), (int(tWait2 + length), 0),
                          (int(piPulse), 1), (int(tWait2 + length), 0), (int(halfpiPulse), -1),
                          (int(tWait + lasRead + tWait + lasRef + tWait), 0)]
                seqAct += [(int(
                    tWait + lasInit + tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait + delay),
                            0), (readTime, 1), (int(lasRead - readTime - delay + tWait + lasRef + tWait), 0)]
                seqRef += [(int(
                    tWait + lasInit + tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait + lasRead + tWait + delay),
                            0), (readTime, 1), (int(lasRef - readTime - delay + tWait), 0)]
            else:
                seqLas += [(int(tWait), 0), (int(lasInit), 1), (
                int(tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait), 0),
                           (int(lasRead), 1), (int(tWait), 0), \
                           (int(lasInit), 1), (
                           int(tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait),
                           0), (int(lasRead), 1), (int(tWait), 0)]
                seqMW += [(int(tWait + lasInit + tWait), 0), (int(halfpiPulse), 1), (int(tWait2 + length), 0),
                          (int(piPulse), 1), (int(tWait2 + length), 0), (int(halfpiPulse), 1),
                          (int(tWait + lasRead + tWait), 0), \
                          (int(lasInit + tWait), 0), (int(halfpiPulse), 1), (int(tWait2 + length), 0),
                          (int(piPulse), 1), (int(tWait2 + length), 0), (int(halfpiPulse), -1),
                          (int(tWait + lasRead + tWait), 0)]
                seqAct += [(int(
                    tWait + lasInit + tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait + delay),
                            0), (readTime, 1), (int(lasRead - readTime - delay + tWait), 0), \
                           (int(
                               lasInit + tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait + lasRead + tWait),
                            0)]
                seqRef += [(int(
                    tWait + lasInit + tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait + lasRead + tWait),
                            0), \
                           (int(
                               lasInit + tWait + halfpiPulse + tWait2 + length + piPulse + length + tWait2 + halfpiPulse + tWait + delay),
                            0), (readTime, 1), (int(lasRead - readTime - delay + tWait), 0)]

    seq = Sequence()
    seq.setDigitalChannel(lasChan, seqLas)
    seq.setAnalogChannel(Ichan, seqMW)
    seq.setDigitalChannel(actChan, seqAct)
    seq.setDigitalChannel(refChan, seqRef)
    return seq.getSequence()


def seqCpmg(start, stop, steps, piPulse, samples=1, nPulses=2, lasChan=2, countReadChan=0, countRefChan=4, lasInit=5e3,
            lasRead=5e3, readTime=3e2, readState=0):
    seq1 = seqInit(lasChan=lasChan, lasInit=lasInit)
    seq2 = seqRead(lasChan=lasChan, countReadChan=countReadChan, lasRead=lasRead, readTime=readTime)
    seq3 = seqRead(lasChan=lasChan, countReadChan=countRefChan, lasRead=lasRead, readTime=readTime)
    Ichan = 0
    Qchan = 1
    seqI0 = []
    seqI1 = []
    seq_final = []
    tWait = 0
    halfpiPulse = int(piPulse / 2)
    seq_operate0 = Sequence()
    seq_operate1 = Sequence()
    lengths = np.linspace(start, stop, steps)
    if lengths[0] == 0:
        lengths[0] = 100
    for length in lengths:
        seqI0 = [(int(1e1), 0), (halfpiPulse, 1), (int((tWait + length) / 2), 0)]
        seqI1 = [(int(1e1), 0), (halfpiPulse, 1), (int((tWait + length) / 2), 0)]
        for i in range(nPulses):
            if i == nPulses - 1:
                seqI0 += [(int(piPulse), 1)]
                seqI1 += [(int(piPulse), 1)]
            else:
                seqI0 += [(int(piPulse), 1), (int(tWait + length), 0)]
                seqI1 += [(int(piPulse), 1), (int(tWait + length), 0)]

        seqI0 += [(int((tWait + length) / 2), 0), (halfpiPulse, 1), (int(1e1), 0)]
        seqI1 += [(int((tWait + length) / 2), 0), (halfpiPulse, -1), (int(1e1), 0)]

        seq_operate0.setAnalogChannel(Ichan, seqI0)
        seq_test0 = seq_operate0.getSequence()

        seq_operate1.setAnalogChannel(Ichan, seqI1)
        seq_test1 = seq_operate1.getSequence()

        if readState == 0:
            seq_final += (seq1 + seq_test1 + seq2 + seq3) * samples
        elif readState == 1:
            seq_final += (seq1 + seq_test0 + seq2 + seq3) * samples
        else:
            seq_final += (seq1 + seq_test1 + seq2 + seq1 + seq_test0 + seq3) * samples
    return seq_final


def seqXY(start, stop, steps, piPulse, samples=1, nPulses=8, lasChan=2, countReadChan=0, countRefChan=4, lasInit=5e3,
          lasRead=5e3, readTime=3e2, readState=0):
    seq1 = seqInit(lasChan=lasChan, lasInit=lasInit)
    seq2 = seqRead(lasChan=lasChan, countReadChan=countReadChan, lasRead=lasRead, readTime=readTime)
    seq3 = seqRead(lasChan=lasChan, countReadChan=countRefChan, lasRead=lasRead, readTime=readTime)
    Ichan = 0
    Qchan = 1
    seqI0 = []
    seqI1 = []
    seqQ = []
    seq_final = []
    tWait = 0
    halfpiPulse = int(piPulse / 2)
    seq_operate0 = Sequence()
    seq_operate1 = Sequence()
    lengths = np.linspace(start, stop, steps)
    if lengths[0] == 0:
        lengths[0] = 100
    if nPulses % 2 != 0:
        raise Exception('number of pulses must be even')
    for length in lengths:
        seqI0 = [(int(1e1), 0), (halfpiPulse, 1), (int((tWait + length) / 2), 0)]
        seqI1 = [(int(1e1), 0), (halfpiPulse, 1), (int((tWait + length) / 2), 0)]
        seqQ = [(int(1e1), 0), (halfpiPulse, 0), (int((tWait + length) / 2), 0)]
        for i in range(int(nPulses / 2)):
            #            print(i)
            #             if i == nPulses - 1:
            # #                 break
            #                 seqQ += [(int(piPulse), 1)]
            #                 seqI0 += [(int(piPulse), 0)]
            #                 seqI1 += [(int(piPulse), 0)]
            if i % 2 == 0:
                #                print('x')
                seqI0 += [(int(piPulse), 1), (int(tWait + length), 0)]
                seqI1 += [(int(piPulse), 1), (int(tWait + length), 0)]
                seqQ += [(int(piPulse), 0), (int(tWait + length), 0)]
            elif i % 2 != 0:
                #                print('y')
                seqQ += [(int(piPulse), 1), (int(tWait + length), 0)]
                seqI0 += [(int(piPulse), 0), (int(tWait + length), 0)]
                seqI1 += [(int(piPulse), 0), (int(tWait + length), 0)]
        for i in range(int(nPulses / 2)):
            if i == int(nPulses / 2) - 1:
                #                 break
                seqQ += [(int(piPulse), 0)]
                seqI0 += [(int(piPulse), 1)]
                seqI1 += [(int(piPulse), 1)]
            elif i % 2 != 0:
                #                print('x')
                seqI0 += [(int(piPulse), 1), (int(tWait + length), 0)]
                seqI1 += [(int(piPulse), 1), (int(tWait + length), 0)]
                seqQ += [(int(piPulse), 0), (int(tWait + length), 0)]
            elif i % 2 == 0:
                #                print('y')
                seqQ += [(int(piPulse), 1), (int(tWait + length), 0)]
                seqI0 += [(int(piPulse), 0), (int(tWait + length), 0)]
                seqI1 += [(int(piPulse), 0), (int(tWait + length), 0)]
        seqI0 += [(int((tWait + length) / 2), 0), (halfpiPulse, 1), (int(1e1), 0)]
        seqI1 += [(int((tWait + length) / 2), 0), (halfpiPulse, -1), (int(1e1), 0)]
        seqQ += [(int((tWait + length) / 2), 0), (halfpiPulse, 0), (int(1e1), 0)]

        seq_operate0.setAnalogChannel(Ichan, seqI0)
        seq_operate0.setAnalogChannel(Qchan, seqQ)
        seq_test0 = seq_operate0.getSequence()

        seq_operate1.setAnalogChannel(Ichan, seqI1)
        seq_operate1.setAnalogChannel(Qchan, seqQ)
        seq_test1 = seq_operate1.getSequence()

        if readState == 0:
            seq_final += (seq1 + seq_test1 + seq2 + seq3) * samples
        elif readState == 1:
            seq_final += (seq1 + seq_test0 + seq2 + seq3) * samples
        else:
            seq_final += (seq1 + seq_test1 + seq2 + seq1 + seq_test0 + seq3) * samples
    return seq_final


def plotSequence(seq, divisor = 10):
    dataArray = np.zeros([10,1])
    for x in seq:
        chanString = x[1]
        chanString = '0'*(8-len(bin(x[1])[2:]))+bin(x[1])[2:]
#         print((chanString))
        dataHolder = np.zeros([8, int(x[0]/divisor)])
        an_0 = np.zeros([1, int(x[0]/divisor)])
        an_1 = np.zeros([1, int(x[0]/divisor)])
        for i in range(len(chanString)):
            if chanString[7-i] == '1':
#                 print(chanString[7-i])
                dataHolder[i]+=1.
        an_0 += x[2]/32767
        an_1 += x[3]/32767
        dataHolder = np.concatenate((dataHolder, an_0, an_1), axis = 0)
#         print(f'dataHolder shape is: {np.shape(dataHolder)}')
        dataArray = np.concatenate((dataArray, dataHolder), axis = 1)
        seqPlot = go.FigureWidget()
        for i in range(len(dataArray[:,0])):
            if i == 8:
                seqPlot.add_scatter(y = dataArray[i],name='Channel '+ 'I')
            elif i == 9:
                seqPlot.add_scatter(y = dataArray[i],name='Channel '+ 'Q')
            else:
                seqPlot.add_scatter(y = dataArray[i],name='Channel '+ str(i))
    return seqPlot


def runExperiment(sequence, averages, data_set_name, steps, samples, start, stop):
    # Create a new dataset and add parameters - these are default parameters, they should ideally be changed
    data_set = exc.new_data_set(data_set_name)
    data_set.add_parameter(ParamSpec(name='Time', unit='ns', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Act_Counts', unit='', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Ref_Counts', unit='', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Rebased_Counts', unit='', paramtype='array'))
    data_set.add_metadata('Laser_Power', laser.power())
    data_set.add_metadata('Source_Power', source.power())
    data_set.add_metadata('Source_Frequency', source.frequency())
    times = np.linspace(start, stop, steps)
    pulser.reset()
    actCount = CountBetweenMarkers(tagger, 1, 7, -7, steps * samples)
    refCount = CountBetweenMarkers(tagger, 1, 8, -8, steps * samples)
    pulser.setTrigger(start=Start.SOFTWARE)
    pulser.stream(sequence, n_runs=int(samples), final='CONSTANT_ZERO')
    actFinal = np.zeros([steps, 1])
    refFinal = np.zeros([steps, 1])
    rebased = np.zeros([steps, 1])
    data_set.add_result({'Time': times, 'Act_Counts': actFinal, 'Ref_Counts': refFinal, 'Rebased_Counts': rebased})
    # data_set.add_result({'Act_Counts': actFinal})
    # data_set.add_result({'Ref_Counts': refFinal})
    for i in tnrange(averages):
        #         print(i)
        sleep(0.5)
        actCount.clear()
        sleep(0.5)
        refCount.clear()
        sleep(1)
        pulser.startNow()
        #         print(f'pulser started on average {i}')
        while np.count_nonzero(actCount.getBinWidths()) < steps * samples:
            sleep(1)
            #             print(np.count_nonzero(actCount.getBinWidths()))
            pass
        actualData = actCount.getData()
        refData = refCount.getData()
        actualDataArr = actualData.reshape(samples, steps)
        refDataArr = refData.reshape(samples, steps)
        actDataMean = np.mean(actualDataArr, axis=0).reshape(steps, 1)
        refDataMean = np.mean(refDataArr, axis=0).reshape(steps, 1)
        if i == 0:
            actFinal += actDataMean
            refFinal += refDataMean
            #         rebased = actFinal/refFinal
            print(np.count_nonzero(actFinal))
        else:
            actFinal = np.concatenate((actFinal, actDataMean), axis=1)
            refFinal = np.concatenate((refFinal, refDataMean), axis=1)
            rebased = actFinal / refFinal
        #     data_set.modify_results(0, [{'Act_Counts': np.mean(actFinal, axis = 1), 'Ref_Counts': np.mean(refFinal, axis=1)}])
        data_set.modify_result(0, {'Act_Counts': np.mean(actFinal, axis=1), 'Ref_Counts': np.mean(refFinal, axis=1),
                                   'Rebased_Counts': np.mean(rebased, axis=1)})
    #     data_set.modify_result(0, {'Ref_Counts': np.mean(refFinal, axis = 1)})


def runExpFreqSweep(sequence, averages, data_set_name, steps, samples, start, stop):
    data_set = exc.new_data_set(data_set_name)
    data_set.add_parameter(ParamSpec(name='Frequency', unit='ns', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Act_Counts', unit='', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Ref_Counts', unit='', paramtype='array'))
    data_set.add_parameter(ParamSpec(name='Rebased_Counts', unit='', paramtype='array'))
    data_set.add_metadata('Laser_Power', laser.power())
    data_set.add_metadata('Source_Power', source.power())
    data_set.add_metadata('Source_Frequency', source.frequency())
    pulser.reset()
    pulser.setTrigger(start=Start.SOFTWARE)
    pulser.stream(sequence, n_runs=int(samples), final='CONSTANT_ZERO')
    # Create frequency list and counters
    frequencies = np.linspace(start, stop, steps)
    actCount = CountBetweenMarkers(tagger, 1, 7, -7, samples)
    refCount = CountBetweenMarkers(tagger, 1, 8, -8, samples)
    # Create counter containers
    actCountsData = np.zeros([len(frequencies), 1])
    refCountsData = np.zeros([len(frequencies), 1])
    rebCountsData = np.zeros([len(frequencies), 1])
    data_set.add_result({'Frequency': frequencies, 'Act_Counts': actCountsData, 'Ref_Counts': refCountsData,
                         'Rebased_Counts': rebCountsData})
    for av in tnrange(averages):
        #         print(av)
        countsListAct = np.zeros([len(frequencies), 1])
        countsListRef = np.zeros([len(frequencies), 1])
        countsListReb = np.zeros([len(frequencies), 1])
        #         print(counts_data)
        #         print(countsList)
        for i in range(len(frequencies)):
            source.frequency(frequencies[i])
            actCount.clear()
            refCount.clear()
            sleep(0.1)
            pulser.startNow()
            #         sleep(0.2)

            while np.count_nonzero(actCount.getBinWidths()) < len(actCount.getBinWidths()):
                pass
            countsListAct[i] = np.mean(actCount.getData())
            #             print('Counts Act is: ', countsListAct[i])
            countsListRef[i] = np.mean(refCount.getData())
            #             print('Counts Ref is: ', countsListRef[i])
            countsListReb[i] = countsListAct[i] / countsListRef[i]
            #             print('Counts Reb is: ', countsListReb[i])
            #             print(countsList)
            if av == 0:
                actCountsData[i] = countsListAct[i]
                refCountsData[i] = countsListRef[i]
                rebCountsData[i] = countsListReb[i]
            #                 print(counts_data)
            else:
                actCountsData[i] = (av * actCountsData[i] + countsListAct[i]) / (av + 1)
                refCountsData[i] = (av * refCountsData[i] + countsListRef[i]) / (av + 1)
                rebCountsData[i] = (av * rebCountsData[i] + countsListReb[i]) / (av + 1)
            data_set.modify_result(0, {'Act_Counts': actCountsData, 'Ref_Counts': refCountsData,
                                       'Rebased_Counts': rebCountsData})

            #                 datasaver.add_result(('Frequency', source.frequency()), ('Counts', np.mean(counter.getData())))
            actCount.clear()
            refCount.clear()
    #             print('counts_data is: \n', counts_data)
    return frequencies, rebCountsData


class Barney(object):

    def __init__(self,
                 pulser,
                 tagger,
                 exp_name):

        self.pulser = pulser
        self.tagger = tagger
        self.exp_name = exp_name
        self.__sequence = None
        self.__grid = None
        self.__dims = None

    def assign_sequence(self, sequence: Sequence):
        self.__sequence = sequence

    def get_sequence(self):
        return self.__sequence

    def assign_grid(self, origin: np.array, range: np.array, points: np.array):
        """
        Creates grid of spacial points for fluorescence measurements
        :param origin: length 2 array - x,y coords for start of grid
        :param range: length 2 array - grid dimensions
        :param points: length 2 array - number of pixels per axis
        :return: None
        """
        self.__dims = np.array([np.linspace(origin[0], origin[0]+range[0], points[0]),
                                np.linspace(origin[1], origin[1]+range[1], points[1])])
        self.__grid = np.zeros([])

    def run_experiment(self):
        if self.__sequence == None:
            raise ""




