import re
from datetime import datetime
while True:
    lines = []
    while True:
        line = input("Enter a line of data (press Enter for a new line, or leave it empty to finish):\n")
        if not line:
            break
        lines.append(line)

    for line in lines:
        match = re.search(r'Fix:\s+([\d.]+)\s+([\d.]+)\s+(\d{2}/\d{2}/\d{4})\s+(\d{2})(\d{2})(\d{2})', line)
        if match:
            latitude = float(match.group(2))
            longitude = float(match.group(1))
            date_str = match.group(3)
            hour = int(match.group(4))
            minute = int(match.group(5))
            second = int(match.group(6))

            timestamp = datetime.strptime(date_str, "%m/%d/%Y").replace(hour=hour, minute=minute, second=second)
            formatted_timestamp = timestamp.strftime("%d-%b-%Y %H:%M")

            print(f"Latitude:{latitude} Longitude: {longitude} time: {formatted_timestamp} UTC")

