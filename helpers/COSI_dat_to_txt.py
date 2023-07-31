#!/usr/bin/env python

import numpy as np
import struct, getopt, sys

def read(fname):
    #Alex Lowell, 2020
    float_t = np.double
    offset = 0
    meta = {}

    #set dtype based on filename extension
    if fname.endswith('wv'):
        reclen = int(fname.split('.')[-1].replace('wv',''))
        if reclen:
            dtype = [('time',np.uint64),('amplitude',float_t),('waveform',np.int8,(reclen,))]
        else:
            dtype = [('time',np.uint64),('amplitude',float_t)]
    elif fname.endswith('cbdump'):
        dtype = float_t
    elif fname.endswith('.dat'):
        #read first four bytes to determine file format identifier
        with open(fname,'rb') as f:
            bb = f.read(4)
            fid = struct.unpack('I',bb)[0]
            print('format identifier: %s' % hex(fid))
        if (fid == 0x7c2b9a9f) or (fid == 0xd7d92ac6):
            nb = 37 + 32 #29 bytes of format, 32 bytes of filter info
            with open(fname,'rb') as f:
                bb = f.read(nb)
                fields = struct.unpack('=IIBIddddddd',bb)
                meta.update(zip(['format_id','raw_npoints','bit_depth','shaper_npoints','timestep','vertical','offset','rcpole','crzero','hpfpole','gain'],fields))
                print(meta)
            if fid == 0x7c2b9a9f:
                dtype = [('time',np.uint64),('amplitude',float_t)]
            else:
                dtype = [('time',np.uint64),('amplitude',float_t),('tot',np.uint16)]
            if meta['raw_npoints']:
                npint = np.int16
                dtype.append(('waveform',npint,(meta['raw_npoints'],)))
            if meta['shaper_npoints']:
                dtype.append(('shaper',np.double,((2*meta['shaper_npoints']) + 1,)))
            offset = nb
        else:
            raise IOError('unknown format id')

    else:
        raise ValueError('bad filename extension')

    return np.fromfile(fname,dtype=dtype,offset=offset),meta

def main(argv):
    inputfile = ''
    outputfile = ''
    opts, args = getopt.getopt(argv,"hi:o:",["infile=","outfile="])
    for opt, arg in opts:
        if opt == '-h':
            print ('COSI_dat_to_txt.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            inputfile = arg
        elif opt in ("-o", "--outfile"):
            outputfile = arg
    if not inputfile or not outputfile:
        print ('Expected usage: \"COSI_dat_to_txt.py -i <inputfile> -o <outputfile>\"')
        sys.exit()

    data, meta = read(inputfile)
    bincenters = np.arange(int(np.max(data['amplitude']) + 2))
    hist, binedges = np.histogram(data['amplitude'], bincenters-0.5)
    np.savetxt(outputfile, np.array([bincenters[:-1].astype(int), hist.astype(int)]).T, header ='Channel\t Counts')

if __name__ == "__main__":
   main(sys.argv[1:])

