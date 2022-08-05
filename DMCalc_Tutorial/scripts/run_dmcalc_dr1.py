#!$HOME/anaconda3/envs/pulsar/bin/

import numpy as np
import psrchive
import sys
import argparse
import os
import glob
import shutil

def check_file(filename):
    if os.path.isfile(filename):
        return filename
    else:
        raise ValueError("File {} does not exist.".format(filename))

def check_folder(dirname):
    if os.path.isdir(dirname):
        return dirname
    else:
        raise ValueError("Folder {} does not exist.".format(dirname))

def process_arch_dirs(b3dir, b5dir):
    b3files = dict()
    b5files = dict()
    mjds_set = set()
    
    for fname in glob.glob(f"{b3dir}/*"):
        #print(fname)
        
        if fname.split('.')[-1].lower() in ['pdf','txt','dat','ps','gpt','fil']:
            continue
        
        try:
            arch = psrchive.Archive_load(fname)
        except:
            continue
        freq = arch.get_centre_frequency()
        mjd = arch.get_Integration(0).get_start_time().intday()
        
        #print(mjd, fname, freq)
        
        if freq>=300 and freq<=500:
            b3files[mjd] = fname
            mjds_set.add(mjd)
            
    for fname in glob.glob(f"{b5dir}/*"):
        #print(fname)
        
        if fname.split('.')[-1].lower() in ['pdf','txt','dat','ps','gpt','fil']:
            continue
        
        try:
            arch = psrchive.Archive_load(fname)
        except:
            continue
        freq = arch.get_centre_frequency()
        mjd = arch.get_Integration(0).get_start_time().intday()
        
        #print(mjd, fname, freq)
        
        if freq>=1260 and freq<=1460:
            b5files[mjd] = fname
            mjds_set.add(mjd)
    
    mjds_set = sorted(mjds_set)
    b35_files = dict()
   
    for mjd in mjds_set:
        b35_files[mjd] = (b3files.get(mjd), b5files.get(mjd))
    
    return b35_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ephem", help="Pulsar ephemeris (par file)")
    parser.add_argument("b3templ", help="Band 3 Template")
    parser.add_argument("b5templ", help="Band 5 Template")
    parser.add_argument("b3n", help="No. of subbands for band 3", type=int)
    parser.add_argument("b5n", help="No. of subbands for band 5", type=int)
    parser.add_argument("b3_arch_dir", help="Folder containing band 3 PSRFITS archives")
    parser.add_argument("b5_arch_dir", help="Folder containing band 5 PSRFITS archives")
    args = parser.parse_args()
    
    parfile = check_file(args.ephem)
    
    b3templfile = check_file(args.b3templ)
    b5templfile = check_file(args.b5templ)
    
    b3n_sub = args.b3n
    b5n_sub = args.b5n
    
    b3dir = check_folder(args.b3_arch_dir)
    b5dir = check_folder(args.b5_arch_dir)
    
    print("Par file:", parfile)
    print("Template files:", b3templfile, ';', b5templfile)
    print("Frequency subbands:", b3n_sub, b5n_sub)
    print("Archive Directories:", b3dir, ';', b5dir)
    
    # Process archive directories
    b35_archive_files = process_arch_dirs(b3dir, b5dir)
    
    print('{} epochs to process'.format(len(b35_archive_files.keys())))
    print('')

    # Checking band 3 template
    b3templ_archive = psrchive.Archive_load(b3templfile)
    b3templ_bw = b3templ_archive.get_bandwidth()
    b3templ_nchan = b3templ_archive.get_nchan()
    b3templ_dm = b3templ_archive.get_dispersion_measure()
    b3templ_epoch = b3templ_archive.get_Integration(0).get_start_time().intday()
    
    # if b3templ_bw==200 and b3templ_nchan==128:
    if b3templ_bw==200:
        print("Template epoch is {}.".format(b3templ_epoch))
        print("Template has 200 MHz bandwidth.")
        print("Fiducial DM present in the template is {}.".format(b3templ_dm))
    else:
        raise ValueError("Template should have 200 MHz bandwidth. Check this file.")    
    
    # Checking band 5 template
    b5templ_archive = psrchive.Archive_load(b5templfile)
    b5templ_bw = b5templ_archive.get_bandwidth()
    b5templ_nchan = b5templ_archive.get_nchan()
    b5templ_dm = b5templ_archive.get_dispersion_measure()
    b5templ_epoch = b5templ_archive.get_Integration(0).get_start_time().intday()
    
    # if b5templ_bw==200 and b5templ_nchan==128:
    if b5templ_bw==200:
        print("Template epoch is {}.".format(b5templ_epoch))
        print("Template has 200 MHz bandwidth.")
        print("Fiducial DM present in the template is {}.".format(b5templ_dm))
    else:
        raise ValueError("Template should have 200 MHz bandwidth. Check this file.")
    
    if b5templ_epoch != b3templ_epoch:
        raise ValueError("The band 3 and band 5 templates should be from the same epoch.")
    
    # Splitting the band 3 template into two for 100 MHz observations.
    split_b3nchan = b3templ_nchan/2
    os.system("psrsplit -c {} {}".format(split_b3nchan, b3templfile))
    # Pick the correct file from the two split files.
    b3templfile_100 = '.'.join(b3templfile.split('.')[:-1] + ["0001_0000"] + b3templfile.split('.')[-1:])
    b3templfile_100_other = '.'.join(b3templfile.split('.')[:-1] + ["0000_0000"] + b3templfile.split('.')[-1:])
    os.unlink(b3templfile_100_other)
    
    b3templfiles = { 100 : b3templfile_100,
                     200 : b3templfile      }
    
    # Splitting the band 5 template into two for 100 MHz observations.
    split_b5nchan = b5templ_nchan/2
    os.system("psrsplit -c {} {}".format(split_b5nchan, b5templfile))
    # Pick the correct file from the two split files.
    b5templfile_100 = '.'.join(b5templfile.split('.')[:-1] + ["0001_0000"] + b5templfile.split('.')[-1:])
    b5templfile_100_other = '.'.join(b5templfile.split('.')[:-1] + ["0000_0000"] + b5templfile.split('.')[-1:])
    os.unlink(b5templfile_100_other)
    
    b5templfiles = { 100 : b5templfile_100,
                     200 : b5templfile      }
    
    for mjd, (b3_arch_file, b5_arch_file) in b35_archive_files.items():
        print('')
        print("Processing", mjd, str(b3_arch_file), str(b5_arch_file))
        
        b3_arch_file_tmp = b5_arch_file_tmp = None
        
        if b3_arch_file is not None:
            b3_archive = psrchive.Archive_load(check_file(b3_arch_file))
            
            # Check 1
            b3_nsubints = b3_archive.get_nsubint()
            if b3_nsubints != 1:
                raise ValueError("The band 3 archive should be time scrunched. Use `$ pam -T` to do this.")
            
            # Check 2
            b3_bw = b3_archive.get_bandwidth()
            if b3_bw < 0:
                raise ValueError("The band 3 archive should have positive bandwidth. Use `$ pam --reverse_freq` to do this.")
            
            # Check 3
            b3_ra = b3_archive.get_coordinates().ra().getRadians()
            b3_dec = b3_archive.get_coordinates().dec().getRadians()
            if b3_ra==0 and b3_dec==0:
                raise ValueError("The archive coordinates should be corrected. Use `$ update_coords.sh` to do this.")
            
            # Check 4
            b3_freq = b3_archive.get_centre_frequency()
            if not (b3_freq>=300 and b3_freq<=500):
                raise ValueError("The archives are not of the correct Band(s). Check these files.")
            
            # Template selection
            b3_templfile_select = b3templfiles[int(b3_bw)]
            
            b3_arch_file_tmp = b3_arch_file + ".tmp"
            shutil.copyfile(b3_arch_file, b3_arch_file_tmp)
            
            # DM correction
            # b3_pdmp_dm = os.popen(f"pdmp -g /NULL {b3_arch_file_tmp} | grep \"Best DM\" | tr -s ' '", 'r', 1).read().split()[3]
            # os.system(f"pam -d {b3_pdmp_dm} -m {b3_arch_file_tmp}")
            
            # Frequency scrunching
            b3_nchan_new = int(b3_bw*split_b3nchan/100)
            os.system(f"pam --setnchn {b3_nchan_new} -m {b3_arch_file_tmp}")
            
        if b5_arch_file is not None:
            
            b5_archive = psrchive.Archive_load(check_file(b5_arch_file))
        
            # Check 1
            b5_nsubints = b5_archive.get_nsubint()
            if b5_nsubints != 1:
                raise ValueError("The archive should be time scrunched. Use `$ pam -T` to do this.")
            
            # Check 2
            b5_bw = b5_archive.get_bandwidth()
            if b5_bw < 0:
                raise ValueError("The archive should have positive bandwidth. Use `$ pam --reverse_freq` to do this.")
        
            # Check 4
            b5_ra = b5_archive.get_coordinates().ra().getRadians()
            b5_dec = b5_archive.get_coordinates().dec().getRadians()
            if b5_ra==0 and b5_dec==0:
                raise ValueError("The archive coordinates should be corrected. Use `$ update_coords.sh` to do this.")
        
            # Check 5
            b5_freq = b5_archive.get_centre_frequency()
            if not (b5_freq>=1260 and b5_freq<=1460):
                raise ValueError("The archives are not of the correct Band(s). Check these files.")
            
            # Template selection
            b5_templfile_select = b5templfiles[int(b5_bw)]
        
            b5_arch_file_tmp = b5_arch_file + ".tmp"
            shutil.copyfile(b5_arch_file, b5_arch_file_tmp)
        
            # DM correction
            # if b3_arch_file is not None:
            #     os.system(f"pam -d {b3_pdmp_dm} -m {b5_arch_file_tmp}")
            # else:
            #     b5_pdmp_dm = os.popen(f"pdmp -g /NULL {b5_arch_file_tmp} | grep \"Best DM\" | tr -s ' '", 'r', 1).read().split()[3]
            #     os.system(f"pam -d {b5_pdmp_dm} -m {b5_arch_file_tmp}")   
        
            # Frequency scrunching
            b5_nchan_new = int(b5_bw*split_b5nchan/100)
            os.system(f"pam --setnchn {b5_nchan_new} -m {b5_arch_file_tmp}")
           
        # Now run dmcalc.
        if b3_arch_file_tmp is not None and b5_arch_file_tmp is not None:
          dmcalc_cmd = f"python3.9 ../scripts/dmcalc_dr1.py {b3_arch_file_tmp} {b5_arch_file_tmp} -E {parfile} -M {b3_templfile_select} {b5_templfile_select} -b3n {b3n_sub} -b5n {b5n_sub} "
        elif b5_arch_file_tmp is None:
          dmcalc_cmd = f"python3.9 ../scripts/dmcalc_dr1.py {b3_arch_file_tmp} -E {parfile} -M {b3_templfile_select} -nch {b3n_sub} "
        elif b3_arch_file_tmp is None:
          dmcalc_cmd = f"python3.9 ../scripts/dmcalc_dr1.py {b5_arch_file_tmp} -E {parfile} -M {b5_templfile_select} -nch {b5n_sub} "
            
        print("[CMD]",dmcalc_cmd)
        os.system(dmcalc_cmd)
        
        if b3_arch_file_tmp is not None:
            os.unlink(b3_arch_file_tmp)
        if b5_arch_file_tmp is not None:
            os.unlink(b5_arch_file_tmp)
    
    os.unlink(b3templfile_100)
        
    
