"""
datesandtime.py

This file contains all functions related to the management of time variables
"""

import numpy as np

def sec_to_datetime_365day(t, year0=2000, month0=1, day0=1):
    import cftime
    import datetime

    """
    Returns a DatetimeNoLeap object for a given elapsed time, t (seconds)
    """

    t0 = cftime.datetime(year0, month0, day0, calendar='noleap')
    dt = datetime.timedelta(seconds=t)
    return t0 + dt

def sec_to_year(t, year0=2000, month0=1, day0=1):
    t = sec_to_datetime_365day(t, year0=year0, month0=month0, day0=day0)
    return t.year

def sec_to_month(t, year0=2000, month0=1, day0=1):
    t = sec_to_datetime_365day(t, year0=year0, month0=month0, day0=day0)
    return t.month

def sec_to_day(t, year0=2000, month0=1, day0=1):
    t = sec_to_datetime_365day(t, year0=year0, month0=month0, day0=day0)
    return t.day

def sec_to_dayofyear(t, year0=2000, month0=1, day0=1):
    t = sec_to_datetime_365day(t, year0=year0, month0=month0, day0=day0)
    return t.timetuple().tm_yday


# def sec_to_year_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
#     """
#     Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the year using the Met Office 360 day calendar

#     year in form: 2001
#     """

#     secs_in_day = 86400.
#     secs_in_month = 2592000.
#     secs_in_year = 31104000.

#     #Calculate the number of seconds since the beginning of year0 for the starting date
#     t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

#     deltayear = int(np.floor((t+t0)/secs_in_year))
#     year = year0 + deltayear

#     return int(year)  


# def sec_to_month_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
#     """
#     Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the month using the Met Office's 360 day calendar
#     month is an integer ranging from 1 (Jan) to 12 (Dec)
#     """

#     secs_in_day = 86400.
#     secs_in_month = 2592000.
#     secs_in_year = 31104000.

#     #Calculate the number of seconds since the beginning of year0 for the starting date
#     t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

#     dayofyear = int(np.floor(((t + t0) % (360 * 24 * 60 * 60))/(24*60*60) + 1))
#     month = int(np.floor((dayofyear-1)/30)+1)   
#     return month  

# def sec_to_day_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
#     """
#     Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the day of the month using the Met Office 360 day calendar.

#     day is an integer from 1 to 30
#     """

#     secs_in_day = 86400.
#     secs_in_month = 2592000.
#     secs_in_year = 31104000.

#     #Calculate the number of seconds since the beginning of year0 for the starting date
#     t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

#     dayofyear = int(np.floor(((t + t0) % (360 * 24 * 60 * 60))/(24*60*60) + 1))
#     month = int(np.floor((dayofyear-1)/30)+1)   
#     day = dayofyear - (month - 1) * 30

#     return int(day)   

# def sec_to_dayofyear_30day(t, year0=2000, month0=1, day0=1, leapcheck=False):
#     """
#     Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into the day of the year
#     in the Met Office's 360 day calendar (12 x 30 day months)

#     dayofyear is an integer ranging from 1 to 360

#     """
#     secs_in_day = 86400.
#     secs_in_month = 2592000.
#     secs_in_year = 31104000.

#     #Calculate the number of seconds since the beginning of year0 for the starting date
#     t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

#     dayofyear = int(np.floor(((t + t0) % (360 * 24 * 60 * 60))/(24*60*60) + 1))

#     return dayofyear  

# def sec_to_date365(t, year0=2000, month0=1, day0=1, secs_in_day=86400, secs_in_year=365*86400):
#     """
#     Converts the number of seconds elapsed since a certain date (year0,month0,day0,sec0) into a date assuming a 365 day calendar (no leap year)

#     date in form: 2001, 11, 30
#     """
#     import numpy as np

    
#     dayofyear0 = monthday_to_dayofyear((month0, day0))

#     #Calculate the number of whole years that have elapsed since date0
#     nyears = int(np.sign(t)*np.floor(np.abs(t)/secs_in_year))
#     tremain = t - nyears * secs_in_year

#     dayofyear = int(dayofyear + np.sign(tremain)*np.floor(  np.abs(tremain)/secs_in_day ))

#     if dayofyear < 1:
#         nyear = nyear + int(np.floor((dayofyear-1)/365))
#         dayofyear = 
        



#     #Calculate the number of seconds since the beginning of year0 for the starting date
#     t0 = (day0 - 1)*secs_in_day + (month0-1)*secs_in_month

#     deltayear = int(np.floor((t+t0)/secs_in_year))
#     year = year0 + deltayear

#     return int(year)  

# def dayofyear_to_monthday(dayofyear):
#     """
#     Calculates the month and date from the day of year

#     Only accepts a day of year between 1 and 365
#     """

#     if (dayofyear > 365) or (dayofyear < 1):
#         print(f"Cannot have a {dayofyear}th day in a 365 day calendar")
#         return None

#     daysinmonths = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30 ,31, 30 ,31],dtype=int)
#     enddayofmonth = np.cumsum(daysinmonths)

#     endday = enddayofmonth[0]
#     month = 1

#     while (int(dayofyear) > endday):
#         month = month + 1
#         endday = enddayofmonth[month-1]

#     day = int(daysinmonths[month-1] - (enddayofmonth[month-1]-dayofyear))

#     return month, day
    

# def monthday_to_dayofyear(monthday):
#     """
#     Calculates the day of year (1-365) from month, day date  (MM, DD)

#     Only accepts a day of year between 1 and 365
#     """
#     month, day = monthday


#     daysinmonths = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30 ,31, 30 ,31],dtype=int)
#     enddayofmonth = np.cumsum(daysinmonths, dtype=int)

#     if day > daysinmonths[month-1]:
#         print(f"{day}th of the {month}th month does not exist in 365 day calendar")
#         return None

#     else:
#         dayofyear = enddayofmonth[month-1] - ( daysinmonths[month-1] - day )
#         return dayofyear

    



    





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