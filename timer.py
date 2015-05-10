#!/usr/bin/python
import datetime

class Timer:
    is_running = False
    time_start = None
    time_end = None
    delta = None
    duration = 0.0
    def start(self, restart = False):
        if self.is_running and not restart:
            raise RuntimeError("Timer is already running.")
        self.time_start = datetime.datetime.now()
        self.is_running = True

    def stop(self):
        self.time_end = datetime.datetime.now()
        self.is_running = False
        self.delta = self.time_end - self.time_start
        self.duration = self.delta.seconds + self.delta.microseconds/1e6

    def seconds(self):
        if self.delta is None:
            return 0.0
        return self.delta.seconds + self.delta.microseconds * 1.0e-6

def timed_execute(function, *args, **kwargs):
    t = Timer()
    t.start()
    if callable(function):
        result = function(*args, **kwargs)
    else:
        result = eval(function)
    t.stop()
    return t.delta, result

if __name__ == '__main__':
    duration, result = timed_execute(range,12,20)
    print(duration.microseconds)
    print(result)
