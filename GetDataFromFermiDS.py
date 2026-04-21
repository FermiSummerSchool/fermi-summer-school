import subprocess
import re
import time
import requests
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
import sys

# Read the XML file and extract the needed values for the Data Server Query
def read_config(file_path):
    # Parse the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Define a helper function to safely get text or None if tag is missing
    def get_val(tag):
        element = root.find(tag)
        return element.text if element is not None else None

    # Extract your specific fields
    data = {
        "coordfield": get_val("CoordField"),
        "coordsystem": get_val("CoordSystem"),
        "shapefield": get_val("ShapeField"),
        "timefield": get_val("TimeField"),
        "timetype": get_val("TimeType"),
        "energyfield": get_val("EnergyField"),
        "photonOrExtendedOrNone": get_val("photonOrExtendedOrNone").capitalize(),
        "destination": get_val("Destination"),
        "spacecraft": get_val("Spacecraft"),
        "end_url": get_val("EndURL"),
        "output_dir": get_val("OutputDir")
    }

    return data

# Usage
config_data = read_config(sys.argv[1])

# Data Server query with the XML values included.
# We only need to do this once.
command = ['curl', '-X', 'POST','-d','coordfield='+config_data["coordfield"],'-d','coordsystem='+config_data["coordsystem"],'-d','shapefield='+config_data["shapefield"],'-d','timefield='+config_data["timefield"],'-d','timetype='+config_data["timetype"],'-d','energyfield='+config_data["energyfield"],'-d','photonOrExtendedOrNone='+config_data["photonOrExtendedOrNone"],'-d','destination=query','-d','spacecraft='+config_data["spacecraft"], config_data["end_url"] ]

print("cmd ", command)

try:
    # Run the CURL command and capture output
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    
    # Print the standard output (the response)
    print("Response Body:")
    # print(result.stdout)
    # links = re.findall(r'href=["\']?(https?://[^"\'>\s]+)', result.stdout)
    # Regex pattern to match href containing "fermi"
    # 1. href=["\'] : Matches href=" or href='
    # 2. (.*?) : Non-greedy capture of the URL
    # 3. fermi : Matches the required text
    # 4. .*? : Capture any characters after fermi until quote
    # 5. ["\'] : Matches closing " or '
    # pattern = r'href=["https://fermi.gsfc.nasa\'](.*?\bQueryResults\b.*?)["\']'
    pattern = r'href=["\'](https://fermi[^"\'>\s]+)["\']'
    # Find all matches
    links = re.findall(pattern, result.stdout)


    print(links) # Output: ['https://www.example.com']
    # Normally the queries
    time.sleep(30)
    # There's only one result url, this one should be it.
    url = links[0]

    # Step 1: Send an HTTP GET request for the results URL.
    response = requests.get(url)

    # Step 2: Parse the HTML content using BeautifulSoup
    # We use 'html.parser' which is built into Python
    soup = BeautifulSoup(response.text, 'html.parser')

    # Step 3: Find all anchor tags (<a>) and extract the 'href' attribute with FTP in them.
    # Download the data in individual subprocesses.
    for link in soup.find_all('a'):
        href = link.get('href')
        if "FTP" in href:
           subprocess.run(["wget", "-O", "my_data.bin/", href])
           # Using subprocess for security and to catch errors
           subprocess.run(["wget", "-P", config_data["output_dir"], href], check=True) 
           # respone2 = requests.get(href)
           print (href)

except subprocess.CalledProcessError as e:
    print(f"Error occurred: {e}")
