from datetime import datetime as dt
from datetime import timedelta
import math
import time


def toYearFraction(datestring):
    difference = 0
    if int(datestring.split('-')[0]) < 1900:
        real_date = datestring.split('-')[0]
        fake_date = '19' + datestring.split('-')[0][2:]
        difference = int(fake_date) - int(real_date)
        datestring = datestring.replace(real_date, fake_date)
    date = dt.strptime(datestring, '%Y-%m-%d')
    def sinceEpoch(date):
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction - difference


def toDatestring(yearFraction):
    year = math.floor(yearFraction)
    remaining_days = (yearFraction-year)*365
    datestring = str(dt(year,1,1) + timedelta(remaining_days))[:10]
    return datestring
