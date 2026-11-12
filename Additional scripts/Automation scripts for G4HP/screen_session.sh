#!/bin/bash

# Get the current date in YYYYMMDD format
DATE=$(date +"%Y%m%d")

# Define a counter for session numbers
SESSION_NUM=1

# Loop to ensure a unique session name
while screen -list | grep -q "${DATE}_$SESSION_NUM"; do
    ((SESSION_NUM++))  # Increment session number if name exists
done

# Create the screen session with a unique name
SESSION_NAME="${DATE}_$SESSION_NUM"
screen -S "$SESSION_NAME"

# Display the session name
echo "Screen session created: $SESSION_NAME"
