def parse_general_info(data):
    date = None
    time = None

    lines = data.split('\n')
    for line in lines:
        if "Float time : Hour" in line:
            hour = int(line.split(';')[1])
        elif "Float time : Minute" in line:
            minute = int(line.split(';')[1])
        elif "Float time : Second" in line:
            second = int(line.split(';')[1])
        elif "Float time : Day" in line:
            day = int(line.split(';')[1])
        elif "Float time : Month" in line:
            month = int(line.split(';')[1])
        elif "Float time : Year" in line:
            year = int(line.split(';')[1])

    date = f"{year}-{month}-{day}"
    time = f"{hour}:{minute}:{second}"

    return date, time


def parse_coordinates(data):
    latitude_degrees = None
    latitude_minutes = None
    latitude_fractions = None
    latitude_orientation = None
    longitude_degrees = None
    longitude_minutes = None
    longitude_fractions = None
    longitude_orientation = None

    lines = data.split('\n')
    for line in lines:
        if "GPS latitude (°)" in line:
            latitude_degrees = int(line.split(';')[1])
        elif "GPS latitude (minutes)" in line:
            latitude_minutes = int(line.split(';')[1])
        elif "GPS latitude (minutes fractions (4th))" in line:
            latitude_fractions = int(line.split(';')[1])
        elif "GPS latitude orientation" in line:
            latitude_orientation = int(line.split(';')[1])
        elif "GPS longitude (°)" in line:
            longitude_degrees = int(line.split(';')[1])
        elif "GPS longitude (minutes)" in line:
            longitude_minutes = int(line.split(';')[1])
        elif "GPS longitude (minutes fractions (4th))" in line:
            longitude_fractions = int(line.split(';')[1])
        elif "GPS longitude orientation" in line:
            longitude_orientation = int(line.split(';')[1])

    latitude_decimal = latitude_degrees + (latitude_minutes + latitude_fractions / 10000) / 60
    if latitude_orientation == 1:
        latitude_decimal = -latitude_decimal

    longitude_decimal = longitude_degrees + (longitude_minutes + longitude_fractions / 10000) / 60
    if longitude_orientation == 1:
        longitude_decimal = -longitude_decimal

    latitude_degrees_minutes_seconds = convert_to_degrees_minutes_seconds(latitude_degrees, latitude_minutes, latitude_fractions)
    longitude_degrees_minutes_seconds = convert_to_degrees_minutes_seconds(longitude_degrees, longitude_minutes, longitude_fractions)

    return latitude_decimal, latitude_degrees_minutes_seconds, longitude_decimal, longitude_degrees_minutes_seconds


def convert_to_degrees_minutes_seconds(degrees, minutes, fractions):
    total_seconds = (minutes + fractions / 10000) * 60
    seconds = int(total_seconds)
    minutes = int(total_seconds / 60)
    remaining_seconds = total_seconds % 60
    return degrees, minutes, remaining_seconds


while True:
    print("Copy and paste the GPS data (enter 'q' to quit):")
    data = ""

    # Read the input until an empty line is encountered
    while True:
        line = input()
        if line.strip() == "":
            break
        data += line + "\n"

    if data.strip() == "q":
        break

    # Find the index of the first relevant line
    start_index = data.find("-----GPS DATA-----")
    if start_index == -1:
        print("Invalid input. Please make sure the GPS data is included.")
        continue

    # Find the index of the last relevant line
    end_index = data.find("-----END OF LIFE-----")
    if end_index == -1:
        print("Invalid input. Please make sure the GPS data is properly terminated.")
        continue

    # Extract the relevant portion of the data
    relevant_data = data[start_index:end_index]

    latitude_decimal, latitude_degrees_minutes_seconds, longitude_decimal, longitude_degrees_minutes_seconds = parse_coordinates(relevant_data)

    print("Decimal Coordinates:")
    print("Latitude:", latitude_decimal)
    print("Longitude:", longitude_decimal)

    print("Coordinates in Degrees, Minutes, and Seconds:")
    print("Latitude:", latitude_degrees_minutes_seconds)
    print("Longitude:", longitude_degrees_minutes_seconds)

    # Check if the "GENERAL INFORMATION" section exists
    general_info_start_index = data.find("-----GENERAL INFORMATIONS-----")
    general_info_end_index = data.find("-----GPS DATA-----")
    if general_info_start_index != -1 and general_info_end_index != -1:
        general_info_data = data[general_info_start_index:general_info_end_index]
        date, time = parse_general_info(general_info_data)

        print("Date:", date)
        print("Time:", time)

    print("--------------------")

