import os.path
import RPi.GPIO as GPIO
import time
from datetime import datetime
from decimal import *
##import Adafruit_ADS1x15
import math
import random
import smbus
import serial
import serial.tools.list_ports
import struct
import time
import shelve
import traceback
from potentiostat import Potentiostat
import threading
import csv
import numpy
import pandas as pd
import pandas.stats.moments
from time import gmtime, strftime

try:
    
    
                
    useLoadedVars = 0       # Set to 1 unless debugging code. A 1 will tell the RPi to use previous parameters.


    ### Vessel Parameters ###
    # Setup run parameters for each reaction module. Each module should have parameters at the same index in each list.
    experimentName = ['RAA_1', 'RAA_2-27_07_18_2']
    experimentEndDates = ["11/11/2018 10:00", "11/11/2019 10:00"]   # "DD/MM/YYYY hh:mm"
    experimentRunParameters = [[15, 67], [15, 87]]                  # TARGET [H2O2], TEMP(C)
    experimentCurrToConcFunctions = [lambda x: 3.8351*x -1.4538, lambda x: 4.6705*x -0.3215]
    pumpPins = [36,37]                                              # GPIO pins: [RAA 1 Pump, RAA 2 Pump]
    defaultPath = '/home/pi/Documents/RAA/RAA_Run_Data'                 # Parent directory for run data and saved parameters
    sampleRate = float(10**-1)                                          # Hz
    dataWriteRate = float(10**-1)                                       # Hz

    raa_1_Set_temps=experimentRunParameters[0][1]
    raa_1_Set_H2O2=experimentRunParameters[0][0]
    raa_1_electrode=['BASi_New100um']
    raa_1_equation=['5.557x-2.017']
    
    raa_2_Set_temps=experimentRunParameters[1][1]
    raa_2_Set_H2O2=experimentRunParameters[1][0]
    raa_2_electrode=['BASi_Old100um']
    raa_2_equation=['4.8729x-0.0461']



 
##    pump_on_resolution = 15             # seconds



    bus = smbus.SMBus(1)    # Setup GPIO object
##    adc = Adafruit_ADS1x15.ADS1115()    # On board clock object
    numberOfRAAs = len(experimentName);     # Number of reaction modules based on length of list created above
    allTemps = [[] for k in range(numberOfRAAs)]        # Initiate run temp data 
    allH2O2Conc = [[] for k in range(numberOfRAAs)]     # Initiate run concentration data
    all_current = [float('nan') for k in range(numberOfRAAs)]
    new_current_available = [False for k in range(numberOfRAAs)]
    all_threads = [[] for k in range(numberOfRAAs)]
    all_thread_init = [False for k in range(numberOfRAAs)]
    _H2O2Conc_ave_f = [[] for k in range(numberOfRAAs)]
    
    # Check RAA Dates. This checks that dates are in the correct format ("DD/MM/YYYY hh:mm")
    for k in range(numberOfRAAs):
        try:
           x = datetime.strptime(experimentEndDates[k], "%d/%m/%Y %H:%M")
        except:
            print 'Check date format for RAA #%i'%(k+1)
            quit()
    errorNum = 0
    dataLogFileName = list()
    # Create log file if it does not exist
    for k in range(numberOfRAAs):
        dataLogFileName.append(defaultPath + '/' + experimentName[k] + '.csv')

    with open(defaultPath + '/' + experimentName[0] + '.csv', 'wb') as outcsv:
         writer = csv.DictWriter(outcsv,delimiter='\n', fieldnames = ['Experiment Name = %s' %(str(experimentName[0])) ,'raa_1_Set_temps = %s' %(str(raa_1_Set_temps)) ,'raa_1_Set_H2O2 = %s' %(str(raa_1_Set_H2O2)),'raa_1_electrode = %s\r' %(str(raa_1_electrode)),'Equation=%s'%(str(raa_1_equation))])
         writer.writeheader()
         with open(defaultPath + '/' + experimentName[0] + '.csv', 'rb') as incsv:
             reader = csv.reader(incsv)
             writer.writerows(row + [0.0] for row in reader)

         with open(defaultPath + '/' + experimentName[0] + '.csv', 'rb') as incsv:
              reader = csv.reader(incsv)
              writer.writerows(row[:1] + [0.0] + row[1:] for row in reader)

                                                                                                                      
    with open(defaultPath + '/' + experimentName[1] + '.csv', 'wb') as outcsv:
         writer = csv.DictWriter(outcsv, delimiter='\n',fieldnames = ["Experiment Name = %s" %(str(experimentName[1])) ,"raa_2_Set_temps =%s " %(str(raa_2_Set_temps)) ,"raa_2_Set_H2O2 = %s" %(str(raa_2_Set_H2O2)),"raa_2_electrode = %s" %(str(raa_2_electrode)),"Equation=%s" %(str(raa_2_equation))])
         writer.writeheader()
         with open(defaultPath + '/' + experimentName[1] + '.csv', 'rb') as incsv:
             reader = csv.reader(incsv)
             writer.writerows(row + [0.0] for row in reader)

         with open(defaultPath + '/' + experimentName[1] + '.csv', 'rb') as incsv:
              reader = csv.reader(incsv)
              writer.writerows(row[:1] + [0.0] + row[1:] for row in reader)

    with open(defaultPath + '/' + experimentName[0] + '.csv', 'a') as outcsv:
         writer = csv.DictWriter(outcsv, delimiter=',', fieldnames = ['Time' ,'temp', 'raa1_conc_raw','raa2_conc_ave','current','raa1_pump'])
         writer.writeheader()
         
    with open(defaultPath + '/' + experimentName[1] + '.csv','a') as outcsv:
         writer = csv.DictWriter(outcsv, delimiter=',', fieldnames = ['Time' ,'temp', 'raa2_conc_raw','raa2_conc_ave','current','raa2_pump'])
         writer.writeheader()     ,
        
    GPIO.setmode(GPIO.BOARD)    # Setup GPIO map (BOARD = physical pin layout, BCM = GPIO pin numbering)
    
    for k in range(numberOfRAAs):
        print 'Experiment: %s'%experimentName[k]
        print 'File Dir: %s'%dataLogFileName[k]
        print 'End Date: %s'%experimentEndDates[k]
        t_end = datetime.strptime(experimentEndDates[k], "%d/%m/%Y %H:%M")
        t_end = (t_end - datetime.now()).total_seconds()
        print 'Seconds Until %s Ends: %s\n\n'%(experimentName[k],t_end)
    for k in range(len(pumpPins)):
        GPIO.setup(pumpPins[k], GPIO.OUT, initial=GPIO.LOW)

    ### Command List ###
    # The following two commands are for assigning a USB serial port to a reaction module temp or echem
    # controller. This is necessary when all 4 USB ports are being used by controllers. The run can be
    # initialized using a connected mouse and keyboard followed by disconnecting those peripherals and
    # connecting the 4 temp or echem controllers.
    def find_PID_serial_ports(num_RAAs):
        ser_ports = list(serial.tools.list_ports.comports())
        serialPID = [[] for x in range(num_RAAs)]
        for k in range(len(ser_ports)):
            for k1 in range(num_RAAs):
                ser = serial.Serial(ser_ports[k][0], 9600, timeout = 0)
                ser.flushInput()
                ser.flushOutput()
                cT = ser.readlines()
                string = '*00'+str(k1+1)+'G110 \r\r'
                ser.flushInput()
                ser.flushOutput()
                ser.write(string)
                ser.flushInput()
                ser.flushOutput()
                time.sleep(.1)
                cT = ser.readlines()
                ser.flushInput()
                ser.flushOutput()
                ser.close()
                if cT:
                    serialPID[k1] = ser_ports[k][0]
        return serialPID

    # Find and link Arduino to a module number
    def find_arduino_serial_ports(num_RAAs):
        _Ard_Address = [[] for x in range(num_RAAs)]
        ports = list(serial.tools.list_ports.comports())
        for p in ports:
            try:
                dev = Potentiostat(p[0])
                cT = dev.get_device_id()
                print cT
                if float(cT) == 1:
                    print 'Connecting RAA 1 to Arduino via port %s.'%(p[0])
                    _Ard_Address[0] = dev
                if float(cT) == 2:
                    print 'Connecting RAA 2 to Arduino via port %s.'%(p[0])
                    _Ard_Address[1] = dev
            except:
                print 'Could not connect to arduino'
        return _Ard_Address

    # Talk to the temperature controller and get the current temp value.
    def get_current_temperature(_serialPID):
        ser = serial.Serial(_serialPID, 9600, timeout = 0)
        ser.write("*G110 \r")
        ser.flushInput()
        ser.flushOutput()
        time.sleep(.1)
        cT = ser.readlines()
        ser.close()
        return float(cT[0])


    # Log the current experimental values
    def Write_To_RAA_File(log_file_name, new_line):
        current_time = datetime.now()
        current_time=current_time.strftime("%Y-%m-%d-%H-%M-%S")
##        current_time = current_time.timetuple()
##        current_time = str(current_time[0:6])
##        current_time = current_time[1:(len(current_time)-1)]
        current_time = current_time+'-'+str(int(math.floor(math.modf(time.time())[0]*1000)))
        
        
        for k in range(len(new_line)):
            new_line[k] = str(new_line[k])
        new_line = [current_time] + new_line 
        new_line = ','.join(new_line) 
        new_line = new_line + '\n'
        with open(log_file_name, 'a') as RAAFile:
            #print new_line
            RAAFile.write(new_line)            



    # Manage a list of values for temperature and concentration. This is so a running average can be used if
    # necessary.
    def Add_Value_To_Data_List(old_vals, new_val, lenLimit):
        if len(old_vals) >= lenLimit:
            old_vals = old_vals[1:(len(old_vals))]
        old_vals = old_vals + [new_val]
        return old_vals


##    class peristaltic_thread (threading.Thread):
##        def __init__(self, on_fraction, gpio_pin):
##            threading.Thread.__init__(self)
##            self.on_fraction = on_fraction
##            self.gpio_pin = gpio_pin
##        def run(self):
##            try:
##                operate_pump(self.on_fraction, self.gpio_pin)
##            except:
##                pass
##
##    def operate_pump(on_fraction, pin_number):
##        GPIO.output(pin_number, GPIO.HIGH)
##        time.sleep(pump_on_resolution*on_fraction)
##        GPIO.output(pin_number, GPIO.LOW)
        
    # Create thread to run each potentiostat independantly without interrupting main script
    class potentiostat_thread (threading.Thread):
        def __init__(self, threadID, device_address):
            threading.Thread.__init__(self)
            self.threadID = threadID
            self.portNum = device_address
        def run(self):
##            print "Starting " + str(self.threadID)
            try:
                curr = run_rodeo(self.threadID, self.portNum)
             ##  print "Finishing " + str(self.threadID)
            except:
                curr = float('nan')
           ## print 'RAA %i Current: %f'%(self.threadID + 1, curr)
            global all_current
            all_current[self.threadID] = curr
          

    # Potentiostat method used by threads
    def run_rodeo(threadID, device_ID):
        curr_range = '10uA'
        volt_range = '1V'
        low_volt = -0.3             # Volts 
        high_volt = 0.7             # Volts
        low_volt_time = .5          # Seconds
        high_volt_time = 2          # Seconds
        sample_period = 0.1

##        print 'Attempting to talk to Potentiostat %i'%(thread_ID + 1)
        device_ID.set_curr_range(curr_range)
        device_ID.set_volt_range(volt_range)
##        print 'Successfully talked to Potentiostat %i'%(thread_ID + 1)

        device_ID.set_volt(low_volt)
        time.sleep(low_volt_time)
        device_ID.set_volt(high_volt)
        time.sleep(high_volt_time)
        current = device_ID.get_curr()
        device_ID.set_volt(0)
        return current

    # Initialize communication with peripherals
    RAA_Temp_Comm = find_PID_serial_ports(numberOfRAAs)
    Ard_Address = find_arduino_serial_ports(numberOfRAAs)

    
    currentRAA = 0;
    lastSampleTime = time.time()
    lastWriteTime = lastSampleTime
    allRAANotComplete = True
    temperatureErrorCount = 0;
    ardErrorCount = 0;
    
    while allRAANotComplete:
        currentRAA = currentRAA + 1
        if currentRAA > numberOfRAAs:
            currentRAA = 1
        # Check if run is finished
        if datetime.now() > datetime.strptime(experimentEndDates[currentRAA-1], "%d/%m/%Y %H:%M"):
            for key in dir():
                if globals()[key].__class__ is list:
                    if len(globals()[key]) == numberOfRAAs:
                        try:
                            del globals()[key][currentRAA-1]
                        except:
                            pass
            numberOfRAAs = numberOfRAAs - 1
            if numberOfRAAs == 0:
                allRAANotComplete = False
                print "Run over"
            continue
        
        if useLoadedVars == 1 and os.path.isfile(defaultPath + '/wsvars.out'):
            print 'Loading previous session...'
            recovFile = shelve.open(defaultPath + '/wsvars.out', 'r')
            for key in recovFile:
                globals()[key] = recovFile[key]
            recovFile.close()
            # Reset start timer using the elapsed time from previous session
            startTimer = time.time() - (currentTimer - startTimer)
            currentTimer = time.time()
            useLoadedVars = 0

            
        pumpPinVals = [0 for k in range(numberOfRAAs)]
        pumpStat = [0 for k in range(numberOfRAAs)]
        
        # Determine State of Pins
        for k in range(numberOfRAAs):
            try:
##                _H2O2Temp = []
##                k1 = len(allH2O2Conc[k]) - 1
##                while _H2O2Temp == []:
##                    if not math.isnan(allH2O2Conc[k][k1]):
##                        _H2O2Temp = allH2O2Conc[k][k1]
##                    elif k1 == 0:
##                        _H2O2Temp = float('nan')
##                    k1 -= 1
                _H2O2Temp = allH2O2Conc[k][-1]
                if _H2O2Temp < experimentRunParameters[k][0]:
                    pumpPinVals[k] = 1;
            except:
                pass
        GPIO.output(pumpPins, pumpPinVals)
        
        
        
        
        # Get data at frequency as defined in parameters
        if (time.time()-lastSampleTime) > (1/sampleRate):
            for k in range(numberOfRAAs):
                try:
                    # Check if potentiostat is connected
                    if not Ard_Address[k] == []:
                        # If potentiostat is connected, initialize thread as necessary
                        if not all_thread_init[k]:
                            all_threads[k] = potentiostat_thread(k, Ard_Address[k])
                            all_threads[k].start()
                            all_thread_init[k] = True
                        else:
                            # If initialized, check if currently running
                            if not all_threads[k].isAlive():
                                all_threads[k] = potentiostat_thread(k, Ard_Address[k])
                                all_threads[k].start()
                    H2O2Val = all_current[k]
                    H2O2Val = experimentCurrToConcFunctions[k](H2O2Val)
                    
                    
##                    concTemp = [x for x in allH2O2Conc[k] + [H2O2Val] if not math.isnan(x)]
                    concTemp = allH2O2Conc[k] + [H2O2Val]
                    H2O2ValLast = []
                    k1 = len(concTemp) - 1
                    while H2O2ValLast == []:
                        if not math.isnan(concTemp[k1]):
                            H2O2ValLast = concTemp[k1]
                        elif k1 == 0:
                            H2O2ValLast = float('nan')
                        k1 -= 1
                        
                   ##                    if concTemp == []:
##                        concTemp = float('nan')
##                    else:
##                        concTemp = float(sum(concTemp)/max(len(concTemp), 1))
##                    print allH2O2Conc[k]
##                    print 'RAA %i Conc: %s'%(k+1, str(H2O2Val))
                    if math.isnan(H2O2ValLast):
                        ardErrorCount += 1
                    if ardErrorCount > 5:
                        print 'Attempting to reconnect to Arduinos...\n'
                        Ard_Address = find_arduino_serial_ports(numberOfRAAs)
                        all_thread_init = [False for k in range(numberOfRAAs)]
                        ardErrorCount = 0
                except:
                    H2O2Val = float('nan')
                    ardErrorCount = ardErrorCount + 1
                    if ardErrorCount > 5:
                        print 'Attempting to reconnect to Arduinos...\n'
                        Ard_Address = find_arduino_serial_ports(numberOfRAAs)
                        all_thread_init = [False for k in range(numberOfRAAs)]
                        ardErrorCount = 0
                    print 'Failed to read concentration for %s.'%(experimentName[k])
                allH2O2Conc[k] = Add_Value_To_Data_List(allH2O2Conc[k], H2O2Val, 20)

##                _H2O2Conc_ave_f[k]= pd.DataFrame(data= allH2O2Conc[k])
##                _H2O2_size=len(_H2O2Conc_ave_f[k])
##               
##                if _H2O2_size >1:   
##               
##                        _H2O2Conc_ave=pd.rolling_mean( _H2O2Conc_ave_f[k],10, center=True)
##                        
##                        if len(_H2O2Conc_ave)== 14:
##                             for i in range(9,14):
##                               mov_ave=_H2O2Conc_ave.iloc[i]
##                             print 'Moving average: %s' %(mov_ave)
##                        elif len(_H2O2Conc_ave) == 20:
##                               mov_ave= _H2O2Conc_ave.iloc[15]
##                               print 'Moving average: %s' %(mov_ave)

                        

##                elif _H2O2_size ==0:
##                  
##                        _H2O2Conc_ave=0
##                    

                try:
                    currentTemp = get_current_temperature(RAA_Temp_Comm[k])
                except:
                    currentTemp = float('nan')
                    temperatureErrorCount = temperatureErrorCount + 1
                    if temperatureErrorCount > 5:
                        print 'Attempting to reconnect to PID controllers...\n'
                        RAA_Temp_Comm = find_PID_serial_ports(numberOfRAAs)
                        temperatureErrorCount = 0
                    print 'Failed to read temp for %s.'%(experimentName[k])
                allTemps[k] = Add_Value_To_Data_List(allTemps[k], currentTemp, 20)
            lastSampleTime = time.time()

            ## Calculating the moving average for the concentrations

       

        # Write data to file at frequency as defined in parameters
        if (time.time()-lastWriteTime) > (1/dataWriteRate):
            print 'Time(HH/MM/SS): %s\n' %strftime("%H %M %S")
            for k in range(numberOfRAAs):
                _H2O2Conc = [x for x in allH2O2Conc[k] if not math.isnan(x)]
                _H2O2Conc = float(sum(_H2O2Conc)/max(len(_H2O2Conc),1))
                _Temperature = [x for x in allTemps[k] if not math.isnan(x)]
                _Temperature = float(sum(_Temperature)/max(len(_Temperature),1))
                _Pump=pumpPinVals[k]
                _current= all_current[k]
               
                print 'RAA%i Conc_raw: %.2f Conc_ave: %.2f Current: %.2f uA   Temp: %.2f    Pump(ON/OFF):%.2f\n'%(k+1,H2O2ValLast,_H2O2Conc, _current, _Temperature,_Pump)

                _H2O2Conc = float("{0:.3f}".format(_H2O2Conc))
                H2O2ValLast = float("{0:.3f}".format(H2O2ValLast))
                _current= float("{0:.1f}".format(_current))
                _Temperature = float("{0:.1f}".format(_Temperature))                                       
                
                # Write variables to file in case of failure
                if os.path.isfile(defaultPath + '/wsvars.out'):
                    os.remove(defaultPath + '/wsvars.out')
                recovFile = shelve.open(defaultPath + '/wsvars.out', 'n')
                for key in dir():
                    try:
                        recovFile[key] = globals()[key]
                    except:
##                        print('Error shelving: {0}'.format(key))
                        pass
                recovFile.close()
                useLoadedVars = 0

                # Write RAA vessel data to file
                Write_To_RAA_File(dataLogFileName[k], [_Temperature , \
                                                       
                                                       H2O2ValLast , \
                                                       
                                                       _H2O2Conc ,\
                                                       
                                                       _current , \
                                                       
                                                       _Pump])
                                                       
                                                       
            lastWriteTime = time.time()
            





except:
    traceback.print_exc()
    pass
finally:
    GPIO.cleanup()
    ports = list(serial.tools.list_ports.comports())
    for p in ports:
        try:
            dev = Potentiostat(p[0])
            dev.stop_test()
        except:
            pass

