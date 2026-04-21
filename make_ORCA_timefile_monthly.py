################################################################
# make_timefile_monthly.py
# --------------------------------------------------------------
# Description: Script to create timeFile seeding input file for
# TRACMASS simulations. Use to seed water parcels on the
# first available day of each month using the TRACMASS seeding
# timefile format.
# --------------------------------------------------------------
# Created On: 2022-07-19
# Created By: Ollie Tooth
# Email: oliver.tooth@seh.ox.ac.uk
###############################################################
# Import packages.
import numpy as np
import pandas as pd

# --------------------------------------------------------------
# Setting the start and end dates of our TRACMASS simulation -
# MODIFY: dates as required for Forward-in-Time.
start = '1992-12-01'
end = '2012-12-16'
target_months = [12] # List of calendar months to include, e.g. [10,11,12]. Set to None to include all calendar months

# Choosing Forward-in-Time or Backward-in-Time for seeding Time-Levels.
# MODIFY: Set at True for evaluating trajectories backward-in-time.
BinT = True

# Setting the interval between the output velocity field (e.g 5 days) -
# MODIFY: interval of tracer fields as required.
interval = '5D'

# Creating dates DatetimeIndex object to store the dates of our
# model velocity fields/output times of TRACMASS simulation.
dates = pd.DataFrame(pd.date_range(start,
                                   end,
                                   freq=interval), columns=['Date'])

# Splitting the date in YYYY-MM-DD format into Year, Month and Day variables.
dates['Year'] = dates['Date'].dt.year
dates['Month'] = dates['Date'].dt.month
dates['Day'] = dates['Date'].dt.day
dates['Time-Level'] = dates.index

# Finding the first available date for each month for all years.
df = dates.groupby(["Year", "Month"]).agg({'Date': np.min})

# If target months specified, only include dates within the specified calendar months
if target_months is not None:
      df = df[df.index.get_level_values('Month').isin(target_months)]

print(df)

# For BinT experiments reverse the order of first available dates.
if BinT == True:
	# Re-order dates from end to start.
	dates = dates.sort_index(axis=0, ascending=False)
	# Re-order Time-Level variable.
	dates['Time-Level'] = np.max(dates['Time-Level'].values) - dates['Time-Level']
	# Re-order first available dates from end to start.
	df = df.sort_index(axis=0 ,ascending=False)

# --------------------------------------------------------------

# Initialising date_index variable with zeros.
date_index = np.zeros(len(df))

# Initialising time_level variable as empty list.
time_level = []

# Find the time-level (row no.) of the first day of available
# output each month once each decade and appending to time_level list.

# Initialising count variable for indexing.
count = 0

# Outer for loop ensures monthly time levels are drawn from df.
for i in range(0, len(df.Date.values)):

    # Add 1 to location given by np.where for Fotran 1-based indexing.
    date_index[count] = np.where(dates.Date.values == df.Date.values[i])[0] + 1

    # Populating the time_level empty list with date_index as a
    # string padded with zeros to the left (i12 format).
    time_level.append(str(int(date_index[count])).zfill(12))

    count += 1

# --------------------------------------------------------------
# Using a DataFrame to export time_level variable to a .csv file with
# commas terminating each line.
# MODIFY: output filename as required.
filename = '/gws/nopw/j04/bas_pog/astyles/DRAKKAR_SET/seeding/seeding_subsampled_3x398/ORCA_timefile_20ydecember.csv'
pd.DataFrame(time_level).to_csv(filename,
                                header=None,
                                index=None,
                                lineterminator=',\n'
                                )
