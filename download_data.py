import os
from google_drive_downloader import GoogleDriveDownloader

dl = GoogleDriveDownloader()

URL = '1uw_amOFLsv-RT1fc3Wzt4WntkfGhQ_GW'
DEST = os.sep.join(['.', 'resources', 'data.zip'])

print(DEST)

dl.download_file_from_google_drive(URL, DEST, unzip=True)
