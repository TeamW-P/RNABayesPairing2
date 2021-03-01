import threading
import os

class Counter(object):

    CURRENT_DIRECTORY = os.path.dirname(__file__)

    def __init__(self):
        self._lock = threading.Lock()
        self.state_file = open(os.path.join(self.CURRENT_DIRECTORY, "counter.txt"),"r+") 
        self.value = int(self.state_file.read()[0])
        
    def increment_and_get(self):
        '''
        Acquires a lock, updates the counter and returns the current value. 

        :returns: the current counter value
        '''
        with self._lock:
            self.value += 1
            self.state_file.seek(0)
            self.state_file.write(str(self.value))
            return self.value

