import datetime

__start = datetime.datetime.now()

def uptime():
    now = datetime.datetime.now()
    dt = now - __start
    return dt.total_seconds()
