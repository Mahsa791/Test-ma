from potentiostat import Potentiostat
import matplotlib.pyplot as plt

port = '/dev/ttyACM2'    # Serial port for potentiostat device
datafile = 'data_Constant.txt'    # Name of output data file

test_name = 'constant'   # Name of test to run - constant volate voltammetery
curr_range = '100uA'     # Name of current range for test [-10uA, +10uA]
sample_rate = 100.0      # Rate (samples/sec) at which to collect samples
num_cycles = 3             # 
conc=[]
test_param = { 
        'quietValue' : -0.3,        # Output voltage during quiet peroid
        'quietTime'  : 500,       # Duration of quiet period (ms)
        'value'      : 0.7,        # Output volatage (V) durring constant voltage test
        'duration'   : 2000,       # Duration of constant voltage test (ms)
        'numCycles'  : num_cycles,
        'conc'       : []
        }

# Create Device object and set sample rate, current range and test parameters
dev = Potentiostat(port)             
dev.set_sample_rate(sample_rate)   
dev.set_curr_range(curr_range)     
dev.set_param(test_name,test_param)

# Run cyclic voltammetry test

t,volt,curr = dev.run_test(test_name,display='pbar',filename=datafile)

for k in range(len(curr)):
    
     conc=(curr[k]*4.802)+0.739
     
     print conc
               
 
# plot results using matplotlib
plt.subplot(211)
plt.title('Voltage and current vs time')
plt.plot(t,volt)
plt.ylabel('potential (V)')
plt.ylim(0,test_param['value']*1.1)
#plt.xlim(0,12)
plt.grid('on')

plt.subplot(212)
plt.plot(t,curr)
plt.ylabel('current (uA)')
plt.xlabel('time (sec)')
#plt.xlim(0,12)
plt.grid('on')

plt.show()






