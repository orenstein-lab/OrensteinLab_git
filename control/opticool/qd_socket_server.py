#!/usr/bin/env python3
#
# qd_socket_server.py
#

import socket, select, sys
if sys.platform == 'win32':
    import msvcrt
from parse_inputs import inputs, instrumentType
from qdcommandparser import QdCommandParser

PORT = 5000
LINE_TERM = '\r\n'

def main():
    
    # set up socket connection
    instrumentInfo = inputs()
    instrument, simulateMode, host = instrumentInfo.parseInput(sys.argv[1:])
    
    qdCommand = QdCommandParser(instrument, line_term = LINE_TERM, simulateMode=simulateMode)
    
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server_socket.bind((host, PORT))
    server_socket.listen(10)
    
    # Dictionary to keep track of sockets and addresses.
    # Keys are sockets and values are addresses.
    # Add server socket to the dictionary first.
    socket_dict = {server_socket: (host, PORT)}
    
    print('Server started at {0} on port {1}.'.format(host, PORT))
    if sys.platform == 'win32':
        quitKeys = "ESC"
    else:
        quitKeys = "ctrl-c"
    print(f'Press {quitKeys} to exit.')
    
    keep_going = True
    cmd_buffer = ''
    while keep_going:
        # Get the list of sockets which are ready to be read through select
        read_sockets = select.select(socket_dict.keys(), [], [], 1)[0]

        # Windows looks for the ESC key to quit.
        if sys.platform == 'win32':
            if msvcrt.kbhit() and msvcrt.getch().decode() == chr(27):
                print ('Server exiting')
                break
    
        for sock in read_sockets:
            # New connection
            if sock == server_socket:
                sock_fd, address = server_socket.accept()
                socket_dict[sock_fd] = address
                print('Client ({0}, {1}) connected.'.format(*address))
                sock_fd.send(b'Connected to QDInstrument Socket Server.\r\n')
    
            # Incoming message from existing connection
            else:
                char = sock.recv(1)
                if char:
                    if (char == b'\r') or (char == 'b\n'):
                        command = cmd_buffer.upper().strip(' ')
                        if command == 'EXIT':
                            sock.send(bytes('Server exiting.' + LINE_TERM, 'utf-8'))
                            print('Server exiting.')
                            keep_going = False
                        elif command == 'CLOSE':
                            sock.send(b'Closing connection.\r\n')
                            print('Client ({0}, {1}) disconnected.'.format(*socket_dict[sock]))
                            socket_dict.pop(sock, None)
                            sock.close()
                        elif len(command):
                            sock.send(bytes(qdCommand.parse_cmd(command) + LINE_TERM, 'utf-8'))
                        command = cmd_buffer = ''
                    elif (char == chr(0xff)):
                        char = sock.recv(2) # discard escape characters from Putty
                    elif (ord(char) > 25) and (ord(char) < 125):
                        cmd_buffer += (char.decode('utf-8')) # keep printable characters
    
    server_socket.close()

if __name__ == '__main__':
    main()