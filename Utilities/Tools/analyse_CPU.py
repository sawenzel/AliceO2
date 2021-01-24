#!/usr/bin/env python3
# A python script to calculate time evolution
# of CPU usage for processes (that were supervised by O2_ROOT/share/scripts/jobutils.sh)

# algorithm: we go trough the file line by line
# each line represents: timestamp pid cpucount_total user_count_total system_count_total
import sys

f=open(sys.argv[1])
corenumber=int(sys.argv[2])

data={}
cpuusage={}
names={}
totalcpu_timesorted={}
totalcpucounter={}


for line in f:
    fields=line.split()
    if len(fields)!=6:
       continue

    timestamp=int(fields[0])
    pid=int(fields[1])
    cputotal=int(fields[2])
    utime=int(fields[3])
    stime=int(fields[4])
    name=fields[5]
    # for each pid create a dictionary that contains an increasing list of entries
    if data.get(pid) == None:
        data[pid]={} # create empty hashmap
        cpuusage[pid]={}
        names[pid]=name

    if totalcpu_timesorted.get(timestamp) == None:
        totalcpu_timesorted[timestamp]=0
        totalcpucounter[timestamp]=cputotal
        
    # fill data
    data[pid][timestamp]=[utime, stime, name]
    if timestamp > 1:
            current=utime
            previous=0
            if data[pid].get(timestamp-1) != None:
                previous=data[pid][timestamp-1][0]
            tmp=(current - previous)/(1.*(totalcpucounter[timestamp] - totalcpucounter[timestamp-1]))*100*corenumber
            cpuusage[pid][timestamp]=tmp
            totalcpu_timesorted[timestamp]+=tmp
        
# for pid in data:
#     print ("ANALYSIS PID " + str(pid) + " " + names[pid])
#     meanusage=0
#     maxusage=0
#     minusage=corenumber*10
#     counter=0
#     for time, value in cpuusage[pid].items():
#         # print (str(time) + " : " + str(cpuusage[pid][time]))
#         meanusage+=value
#         maxusage=max(value,maxusage)
#         minusage=min(value,minusage)
#         counter=counter+1
#     if counter>0:
#         print("STATISTICS " + str(pid) + " " + names[pid] + " mean " + str(meanusage/counter) + " max " + str(maxusage) + " min " + str(minusage) + " counts " + str(counter))

# get a global time cpu util average
globalmean=0
globalmax=0
globalmin=corenumber*10
counter=0
for time, value in totalcpu_timesorted.items():
    if time > 1:
        globalmean+=value
        globalmax=max(value,globalmax)
        globalmin=min(value,globalmin)
        counter=counter+1

if counter>0:
    globalmean/=counter
    print("CPU_MEAN/MAX/MIN:{:.2f}/{:.2f}/{:.2f}".format(globalmean,globalmax,globalmin))
else:
    print("CPU_MEAN/MAX/MIN:NaN/NaN/NaN")


