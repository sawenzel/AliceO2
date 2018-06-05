# starts a simulation device setup
NSIMWORKERS=$1

killall -9 xterm
killall -9 O2PrimaryServerDeviceRunner
killall -9 O2SimDeviceRunner
killall -9 O2HitMergerRunner

# we have one primary distributor 
xterm -geometry 80x25+0+0 -e "O2PrimaryServerDeviceRunner --control static --id primary-server --mq-config ~/alisw_new/O2/run/primary-server.json -n 10 -m TOF TRD TPC -e TGeant3 | tee serverlog;bash" &

for i in `seq 1 ${NSIMWORKERS}`; do
  xterm -geometry 80x25+500+0 -e "O2SimDeviceRunner --control static --id worker${i} --config-key worker --mq-config ~/alisw_new/O2/run/primary-server.json --verbosity INFO  | tee simlog${i};bash" &
done


# the its digitizer
#xterm -geometry 80x25+1000+0 -hold -e "O2ITSDigitizerDeviceRunner --control static --id itsdigitizer --mq-config  ~/alisw_new/O2/run/primary-server.json" &


# one hit merger -> the time measures the walltime of the complete session
time O2HitMergerRunner --id hitmerger --control static --mq-config ~/alisw_new/O2/run/primary-server.json | tee mergelog 

