import sys
import os

# Get the username from the command-line argument
if len(sys.argv) < 2:
    print("Usage: python script.py username")
    sys.exit(1)
username = sys.argv[1]

# Define the source and destination directories
source_dir = '/home/{}/'.format(username)
dest_dir = '{}@apex.webbresearch.com:/home/{}/'.format(username,username)

# Define the rsync command to copy files from source to destination
rsync_cmd_1 = 'rsync -avz --update --include="*.msg" --include="*.log" --exclude="*" {} {}'.format(source_dir, dest_dir)


# Execute the first rsync command to copy files from source to destination
print("Copying to webb:")
print(rsync_cmd_1)
os.system(rsync_cmd_1)

# Execute the second rsync command to copy files from destination to source if there are new files in the destination directory
# Define the rsync command to copy files from destination to source
rsync_cmd_2 = 'rsync -avz --update --include="*.msg" --include="*.log" --exclude="*" {} {}'.format(dest_dir, source_dir)
print(rsync_cmd_2)
print("Copying from webb:")
os.system(rsync_cmd_2)

