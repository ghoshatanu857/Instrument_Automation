# !pip install pulsestreamer
# !pip install nidaqmx
import numpy as np
import time
from tqdm import trange
import pulsestreamer
import pulsestreamer as psl
import nidaqmx
import nidaqmx.stream_readers
from pulsestreamer import PulseStreamer,findPulseStreamers,OutputState,TriggerStart,Sequence,TriggerRearm


# Delay Measurement Sequence
def seqDelay(pulser,laserNum=1,gateStart=5,source=7,rising_delay=2,gatelen = 6, laserontime = 31,delay_pad = 2,delay_shift = 2,gatesourcedelay=2):
    
    seq = pulser.createSequence()
   
    laserNum = 1
    gateStart = 5
    source=7
    
    totaltime= 2*delay_pad + laserontime +2*rising_delay
    # steps=int((totaltime-gatelen-2*rising_delay)/delay_shift)
    steps=int((totaltime-2*gatelen-2*rising_delay)/delay_shift)
    
        
    # i=0
    for i in range(steps):
    # while i<steps:
        seq.setDigital(
           laserNum,
           [
               (int(delay_pad+rising_delay), 0),
               (int(laserontime), 1),
               (int(rising_delay+delay_pad),0),

           ],
        )
        totaltime= 2*delay_pad + laserontime +2*rising_delay
        gatingofftime=totaltime - gatelen - i*delay_shift-rising_delay
        
        seq.setDigital(
           gateStart,
           [
               (int(rising_delay),0),
               (int(gatelen),1),
               (int(i*delay_shift+rising_delay), 0),
               (int(gatelen), 1),
               (int(totaltime-2*rising_delay-2*gatelen-i*delay_shift), 0),
           ],
        )
        time = int(rising_delay+gatelen+rising_delay+i*delay_shift)
        seq.setDigital(
           source,
           [
               (int(rising_delay),0),
               (int(gatelen-gatesourcedelay),1),
               (int(i*delay_shift+rising_delay+gatesourcedelay), 0),
               (int(gatelen-gatesourcedelay), 1),
               (int(totaltime-2*rising_delay-2*gatelen-i*delay_shift-gatesourcedelay), 0),
           ],
        )
        yield seq,time,steps
        # i=i+1


# SNR Measurement Sequence
def seqSNR(pulser,laserNum=1,gateStart=5,source=7,rising_delay = 50,gatelen = 50, laserontime = 3e3,delay_pad = 50,delay_shift = 0.1e3,gatesourcedelay = 5,evolution_time = 5e6):  
    
    seq = pulser.createSequence()
   
    laserNum = 1
    gateStart = 5
    source=7
    
    steps=int((laserontime-gatelen)/delay_shift)
    # print(f'Number of Steps : {steps}')
    
    
    for i in range(steps):
        
        seq.setDigital(
           laserNum,
           [
               (int(delay_pad+rising_delay), 0),
               (int(laserontime), 1),
               (int(rising_delay), 0),
               (int(laserontime), 1),
               (int(rising_delay+evolution_time),0),
               (int(laserontime), 1),
               (int(delay_pad+rising_delay), 0),

           ],
        )
        
        seq.setDigital(
           gateStart,
            [
               (int(delay_pad+rising_delay), 0),
               (int(laserontime), 0),
               (int(rising_delay), 0),
               (int(gatelen+i*delay_shift), 1),
               (int(laserontime-gatelen-i*delay_shift+rising_delay+evolution_time),0),
               (int(gatelen+i*delay_shift), 1),
               (int(delay_pad+rising_delay+laserontime-gatelen+i*delay_shift), 0),

           ],
        )
        
        time = int(gatelen+i*delay_shift)
        
        seq.setDigital(
           source,
            [
               (int(delay_pad+rising_delay), 0),
               (int(laserontime), 0),
               (int(rising_delay), 0),
               (int(gatelen+i*delay_shift-gatesourcedelay), 1),
               (int(laserontime-gatelen-i*delay_shift+gatesourcedelay+rising_delay+evolution_time),0),
               (int(gatelen+i*delay_shift-gatesourcedelay), 1),
               (int(delay_pad+rising_delay+laserontime-gatelen-i*delay_shift+gatesourcedelay), 0),

           ],
        )
        yield seq,time,steps

# T1 Measurement Sequence
def seqT1(pulser,laserNum=1,gateStart=5,source=7,rising_delay = 100,gatelen = 2e3, laserontime = 20e3,delay_pad = 100,delay_shift = 100e3,gatesourcedelay = 5,evolution_time = 5e6):  
    
    seq = pulser.createSequence()
   
    laserNum = 1
    gateStart = 5
    source=7
    
    total_time= delay_pad+rising_delay+laserontime+rising_delay+laserontime+rising_delay+evolution_time+laserontime+rising_delay+delay_pad
    steps=int(evolution_time/delay_shift)
    
    
    for i in range(steps):
        laser_offtime = total_time - delay_pad -3*rising_delay-3*laserontime-i*delay_shift
        seq.setDigital(
           laserNum,
           [
               (int(delay_pad+rising_delay), 0),
               (int(laserontime), 1),
               (int(rising_delay), 0),
               (int(laserontime), 1),
               (int(rising_delay+i*delay_shift),0),
               (int(laserontime), 1),
               (int(delay_pad+rising_delay), 0),

           ],
        )
        
        gate_offtime = total_time - delay_pad -3*rising_delay-2*laserontime-gatelen-i*delay_shift
        seq.setDigital(
           gateStart,
            [
               (int(delay_pad+rising_delay), 0),
               (int(laserontime), 0),
               (int(rising_delay), 0),
               (int(gatelen), 1),
               (int(laserontime-gatelen+rising_delay+i*delay_shift),0),
               (int(gatelen), 1),
               (int(delay_pad+rising_delay+laserontime-gatelen), 0),

           ],
        )
        
        time = int(rising_delay+i*delay_shift)
        
        seq.setDigital(
           source,
            [
               (int(delay_pad+rising_delay), 0),
               (int(laserontime), 0),
               (int(rising_delay), 0),
               (int(gatelen-gatesourcedelay), 1),
               (int(laserontime-gatelen+gatesourcedelay+rising_delay+i*delay_shift),0),
               (int(gatelen-gatesourcedelay), 1),
               (int(delay_pad+rising_delay+laserontime-gatelen+gatesourcedelay), 0),

           ],
        )
        yield seq,time,steps

# getting time axis 
def get_time(pulser,exp_name,specifications): 
    
    delay_time = []; steps=0
    
    if exp_name.lower()=='t1':
        sequence_time=seqT1(pulser,**specifications)
    if exp_name.lower()=='snr':
        sequence_time=seqSNR(pulser,**specifications)
    if exp_name.lower()=='delay':
        sequence_time=seqDelay(pulser,**specifications)
    if exp_name.lower()=='lifetime':
        sequence_time=seqT1(pulser,**specifications)  #change this sequence
        
    for t in sequence_time:
        delay_time.append(t[1])
    delay_time = np.array(delay_time)

    for st in sequence_time:
        steps=st[2]
        break
        
    return delay_time,steps

# measuremet function 
def measure(pulser,device_name,specifications,samples=1000,averages=5):
    
    time_axis,steps=get_time(pulser,exp_name,**specifications)
    print(f'number of steps : {steps}')
    
    numberofpoints=samples*2 
    
    pixel=numberofpoints*steps 
    print(f'Pixel : {pixel}')
    DAQ_device.reset_device()
    pulser.reset()
    print("creating sequence")
   
    # Counter
    CountWidth = nidaqmx.Task()
    ciChannel = CountWidth.ci_channels.add_ci_count_edges_chan('/Dev1/ctr1',edge=nidaqmx.constants.Edge.RISING, initial_count=0,
                                                               count_direction=nidaqmx.constants.CountDirection.COUNT_UP) # which specification are we measuring here?

    CountWidth.triggers.pause_trigger.dig_lvl_src='/Dev1/PFI4'
    CountWidth.triggers.pause_trigger.trig_type=nidaqmx.constants.TriggerType.DIGITAL_LEVEL
    CountWidth.triggers.pause_trigger.dig_lvl_when=nidaqmx.constants.Level.LOW


    #CountWidth.timing.cfg_implicit_timing(sample_mode = nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan=(pixel)*averages)#samps per channel defines the buffer size for the memory
    CountWidth.timing.cfg_samp_clk_timing(rate=1e8,source='/Dev1/PFI5',active_edge=nidaqmx.constants.Edge.FALLING,
                                          sample_mode = nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan=(pixel)*averages )
    cps = []
    callback=[]  
   
    #Pulse streamer gating
    # Digital output
    DigChannel = 'Dev1/port0/line7' #connect this to PFI 4 #this is ctr 1 gate
    DigTask = nidaqmx.Task()
    DigTask.do_channels.add_do_chan(lines = DigChannel)
    DigChannel = 'Dev1/port0/line7' #Defining the port for taking the output
   
   
    def readBuffer(task_handle, every_n_samples_event_type, number_of_samples, callback_data):
        CountWidth.in_stream.read_all_avail_samp = True
        readPixels=readerWidth.read_many_sample_uint32(highCount, number_of_samples_per_channel=- 1, timeout=10.0)
        cps.extend(highCount)
        callback.extend([1])
        return 0

    buffersamplecount=numberofpoints
    
    # Counter read task
    readerWidth = nidaqmx.stream_readers.CounterReader(CountWidth.in_stream)    
    highCount = np.zeros(buffersamplecount, dtype = np.uint32)
    lowCount =  np.zeros(buffersamplecount,dtype = np.uint32)


    # Read after filling the buffer with given number of samples
    CountWidth.register_every_n_samples_acquired_into_buffer_event(buffersamplecount,readBuffer) #after every pixel it will trigger the callback

   
    # Start tasks (digital output will be triggered by analog output)
    print("starting DAQ")
    CountWidth.start()
    
    #Adding infinite loop
    t=0
    run=0
    data=[]
    finaldata=[]
    print("Preparing Ni Daq for the experiment")
    print("callback number in beginning:",len(callback))

    i=0
    for run in trange(averages):

        if exp_name.lower()=='t1':
            sequence_time=seqT1(pulser,**specifications)
        if exp_name.lower()=='snr':
            sequence_time=seqSNR(pulser,**specifications)
        if exp_name.lower()=='delay':
            sequence_time=seqDelay(pulser,**specifications)
        if exp_name.lower()=='lifetime':
            sequence_time=seqT1(pulser,**specifications)  #change this sequence

        pulser.setTrigger(start=psl.TriggerStart.HARDWARE_RISING,rearm=psl.TriggerRearm.AUTO)
        seq_num=0
       
        for s in sequence_time:
            t1=len(callback)
           
            seq_num=seq_num+1
            print(seq_num)
            i+=1
            pulser.stream(s[0], n_runs=samples)         

            DigTask.write(True)
            while len(callback)==t1:
                time.sleep(0.001)         
            DigTask.write(False)
         
        run=run+1
        print(f"callback number after {run}-th average end: {len(callback)}")
        
    print(f'Total Run : {i}')
    data=signal_counts(cps,pixel,numberofpoints,steps)   
    
    CountWidth.close()
    DigTask.close()    

    print('returning averaged counts and time_axis')
    return data,time_axis

# Function to Modify the Data
def signal_counts(pulser,device_name,all_counts,counts_in_one_average,numberofpoints,steps,*args):
    all_counts=np.array(all_counts)
    print(f'Total Counts & Counts in one average : {len(all_counts), counts_in_one_average}')
    no_of_averages=int(len(all_counts)/counts_in_one_average)
    print("Crosscheck number of averges=",no_of_averages)

    # Changing the cumulative counts to actual counts
    cumulative_counts = np.reshape(all_counts,(no_of_averages,counts_in_one_average))
    modified_matrix = np.delete(cumulative_counts, -1, 1)
    zero_array = np.zeros(no_of_averages, dtype=int)
    new_matrix = np.hstack((zero_array[:, np.newaxis], modified_matrix))
    actual_counts = np.subtract(cumulative_counts,new_matrix)
    averaged_actual_counts = np.mean(actual_counts,axis=0)

    
    return averaged_actual_counts 