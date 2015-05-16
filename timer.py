# encoding: utf-8
##############################################################################
#                                                                            #
#                  Ezeplot - Dynamical systems visualisation                 #
#                                                                            #
#   Copyright (C) 2015 Raj Kiran Grandhi <rajkiran@aero.iitkgp.ernet.in>     #
#                                                                            #
#   This program is free software: you can redistribute it and/or modify     #
#   it under the terms of the GNU General Public License as published by     #
#   the Free Software Foundation, either version 3 of the License, or        #
#   (at your option) any later version.                                      #
#                                                                            #
#   This program is distributed in the hope that it will be useful,          #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#   GNU General Public License for more details.                             #
#                                                                            #
#   You should have received a copy of the GNU General Public License        #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                            #
##############################################################################
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
