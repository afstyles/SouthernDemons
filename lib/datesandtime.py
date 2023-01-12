"""
datesandtime.py

This file contains all functions related to the management of time variables
"""

import numpy as np

def sec_to_year_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
    """
    Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the year using the Met Office 360 day calendar

    year in form: 2001
    """

    secs_in_day = 86400.
    secs_in_month = 2592000.
    secs_in_year = 31104000.

    #Calculate the number of seconds since the beginning of year0 for the starting date
    t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

    deltayear = int(np.floor((t+t0)/secs_in_year))
    year = year0 + deltayear

    return int(year)  


def sec_to_month_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
    """
    Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the month using the Met Office's 360 day calendar
    month is an integer ranging from 1 (Jan) to 12 (Dec)
    """

    secs_in_day = 86400.
    secs_in_month = 2592000.
    secs_in_year = 31104000.

    #Calculate the number of seconds since the beginning of year0 for the starting date
    t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

    dayofyear = int(np.floor(((t + t0) % (360 * 24 * 60 * 60))/(24*60*60) + 1))
    month = int(np.floor((dayofyear-1)/30)+1)   
    return month  

def sec_to_day_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
    """
    Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the day of the month using the Met Office 360 day calendar.

    day is an integer from 1 to 30
    """

    secs_in_day = 86400.
    secs_in_month = 2592000.
    secs_in_year = 31104000.

    #Calculate the number of seconds since the beginning of year0 for the starting date
    t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

    dayofyear = int(np.floor(((t + t0) % (360 * 24 * 60 * 60))/(24*60*60) + 1))
    month = int(np.floor((dayofyear-1)/30)+1)   
    day = dayofyear - (month - 1) * 30

    return int(day)   

def sec_to_dayofyear_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
    """
    Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the day of the year
    in the Met Office's 360 day calendar (12 x 30 day months)

    dayofyear is an integer ranging from 1 to 360

    """
    secs_in_day = 86400.
    secs_in_month = 2592000.
    secs_in_year = 31104000.

    #Calculate the number of seconds since the beginning of year0 for the starting date
    t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

    dayofyear = int(np.floor(((t + t0) % (360 * 24 * 60 * 60))/(24*60*60) + 1))

    return dayofyear  


    #  def sec_to_month_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
    # """
    # Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the date 
    # (year, month, day, sec) in the Met Office's 360 day calendar (12 x 30 day months)

    # year in form: 2001
    # month is an integer ranging from 1 (Jan) to 12 (Dec)
    # day is an integer from 1 to 30
    # dayofyear is an integer ranging from 1 to 360

    # """

    # secs_in_hour = 60 * 60
    # secs_in_day = 24 * secs_in_hour
    # secs_in_month = 30 * secs_in_day
    # secs_in_year = 12 * secs_in_month

    # #Calculate the number of seconds since the beginning of year0 for the starting date
    # t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

    # deltayear = int(np.floor((t+t0)/secs_in_year))
    # dayofyear = int(np.floor(((t + t0) % (360 * 24 * 60 * 60))/(24*60*60) + 1))

    # month = int(np.floor((dayofyear-1)/30)+1)   
    # day = dayofyear - (month - 1) * 30
    # year = year0 + deltayear

    # return (year, month, day, dayofyear)  