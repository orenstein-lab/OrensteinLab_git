#!/usr/bin/env python3
#
#  qd_socket_client.py
#  

import time, socket, sys
from parse_inputs import inputs

# function to send a TEMP request
def doTempSend(count):
    while count > 0:
        sock.send(b'temp?\r')
        count -= 1

# function to send a FIELD request
def doFieldSend(count):
    while count > 0:
        sock.send(b'field?\r')
        count -= 1

# function to send a CHAMB request
def doChambSend(count):
    while count > 0:
        sock.send(b'chamb?\r')
        count -= 1
        
# set up socket port
timeWas = time.time()
instrumentInfo = inputs(instrumentRequired=False)
instrument, simulateMode, host = instrumentInfo.parseInput(sys.argv[1:])
port = 5000
addr = (host, port)
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.connect(addr)
print(sock.recv(128).decode('utf-8'))

# send a series of requests
doTempSend(334)
doFieldSend(333)
doChambSend(333)

# close connection
sock.send(b'close\r')

# process the reply
accumData = ''
while True:
    data = sock.recv(1).decode('utf-8')
    # look for newline characters
    if data == '\r': data = '\n'
    if data == '\n':
        # process completed reply string
        if len(accumData):
            print(accumData)
        # start a new reply string
        accumData = ''
    else:
        # append to reply string
        accumData += data
    if 'Closing' in accumData:
        break

# show summary of activity
print ('\n"{0}",{1}'.format(accumData, time.time() - timeWas))
time.sleep(0.01)
sock.close()
