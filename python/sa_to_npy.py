#!/usr/bin/env python3
"""
Convert raw data from a scan archive file to a numpy array.

- Assumes resolution of each readout is the same
"""
import argparse
import logging
from pathlib import Path

import numpy as np
from GERecon import Archive

logging.basicConfig(level=logging.INFO)

def get_useful_packets(archive_filename: str) -> int:
    """
    Read a ScanArchive file and count useful opcodes.

    Args:
        archive_filename (str): Path to an existing ScanArchive*.h5 file.

    Returns:
        sum where opcode = 1

    """
    archive = Archive(archive_filename)
    metadata = archive.Metadata()
    num_control = metadata["controlCount"]
    opcode = []

    # Loop over all the control packets
    for i in range(num_control):
        # Retrieve the next control packet
        control = archive.NextControl()
        opcode.append(control["opcode"])

    return opcode.count(1)

def read_archive(archive_filename: str) -> np.complex64:
    """
    Read a ScanArchive file and convert to a numpy array.

    Args:
        archive_filename (str): Path to an existing ScanArchive*.h5 file.

    Returns:
        data_out: numpy array

    """
    # First figure out how many useful control packets we have got in order to determine the data size
    num_packets = get_useful_packets(archive_filename)
    print(f'Found {num_packets} useful packets')

    # Re-open the archive
    archive = Archive(archive_filename)
    metadata = archive.Metadata()
    print('Metadata: ',metadata)

    xres = metadata["acquiredXRes"]
    yres = metadata["acquiredYRes"]
    num_control = metadata["controlCount"]
    num_channels = metadata["numChannels"]

    # initialise array to store acquired data
    rawdata = np.zeros((xres, num_packets, num_channels), dtype=np.complex64)

    # Loop over all the control packets
    for view in range(num_control):

        # Retrieve the next control packet
        control = archive.NextControl()
        # this is a raw control packet, but still get the next frame so that the control and frames are in sync
        if control["opcode"] == 16:
            next_frame = np.squeeze(archive.NextFrame())

        # this is a programmable control packet, so use next frame to fill a line of rawdata
        elif control["opcode"] == 1:
            # Assume operation == 0
            if control["operation"] == 0:       # view = frame
                rawdata[:, view, :] = archive.NextFrame()
            else:
                raise ValueError('Unknown control operation')

    return rawdata

if __name__=='__main__':
    # parse arguments
    parser = argparse.ArgumentParser(
         description='Extract raw data from a scan archive file and export as numpy array.')

    parser.add_argument('scan_arc', type=Path, help='Path to Scan Archive')
    args = parser.parse_args()

    scan_arc = args.scan_arc

    # Read ScanArchive using Orchestra Python SDK",
    rawdata = read_archive(str(scan_arc))
    print(f'Saving k-space with dimensions: {rawdata.shape}')

    np.save(scan_arc.parent / 'rawdata',rawdata)

