'''
classes and methods for multithreading.
'''

from threading import Thread, Lock, Event
import threading
import multiprocessing as mp
import queue
import keyboard
import time

###
### Concurrency Classes
###

def queue_write(q, message):
    '''
    helper function that takes in a Queue object, clears the contents of the Queue and then rewrites to it in a thread and process safe manner.
    '''
    while not q.empty():
        q.get()
    q.put(message)

def queue_read(q):
    '''
    helper function that takes in a Queue object, reads the contents, and then puts the message back so that future calls can still read the value. Thread and process safe.
    '''
    val = q.get()
    q.put(val)
    return val

class LockedVar:
    '''
    Minimal class to implement a locking variable. Contains two private attributes, a value and a lock, and a few methods for safetly reading writing value via the lock.
    '''
    def __init__(self, val):
        self._value = val
        self._lock = Lock()

    def locked_read(self):
        with self._lock:
            return self._value

    def locked_update(self, val):
        with self._lock:
            self._value = val

class StoppableThread(Thread):
    '''
    Thread class with a stop() method. The thread itself has to check
    regularly for the stopped() condition.
    '''
    def __init__(self,  *args, **kwargs):
        super(StoppableThread, self).__init__(*args, **kwargs)
        self._stop_event = Event()

    def stop(self):
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()

class LockedDict:
    '''
    Minimal class to implement a locking dictionary. Contains two private attributes, a value and a lock, and a few methods for safetly reading writing value via the lock.
    '''
    def __init__(self, dict):
        self._dict = dict
        self._lock = Lock()

    def locked_read(self, key):
        with self._lock:
            return self._dict[key]

    def locked_update(self, key, val):
        with self._lock:
            self._dict[key] = val

    def locked_get_dict(self):
        with self._lock:
            return self._dict

###
### Multithreading Methods
###

def press_any_key_to_stop(run_var):
    print('Press any key to stop:')
    while run_var.locked_read():
        if keyboard.is_pressed():  # Check if the any key is pressed
            run_var.locked_update(False)
        time.sleep(0.1)

def press_esc_to_stop(run_var):
    print('Press \"Esc\" to stop:')
    while run_var.locked_read():
        if keyboard.is_pressed('esc'):  # Check if the 'Esc' key is pressed
            run_var.locked_update(False)
        time.sleep(0.1)
