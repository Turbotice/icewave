import numpy as np
import xarray as xr
import pandas as pd
import os
from datetime import datetime, time


# TODO:
# - create a class for anemometer ?
# - get the location of the anemometers
# - get the height of each anemometers
# - get the notes associated with the measurement


def parse_log_to_xarray(filename, 
                        reference_date=None, 
                        verbose=False,
                        log_errors=True,
                        check_gaps=True,
                        sampling_rate_hz=5):
    """
    Parse a log file with timestamped variable measurements into an xarray Dataset.
    
    Parameters:
    -----------
    filename : str
        Path to the log file
    reference_date : str
        date of the log, format 'DD/MM/YYYY'
    verbose: bool
    log_errors: bool
    
    Returns:
    --------
    xr.Dataset
        Dataset with variables indexed by timestamp
    
    Example line format:
    10:30:33.493 S  04.41 D  321 U  02.69 V -03.30 W -01.18 T -00.60 H  36.69 P  1011.82 PI  000.1 RO -001.7 MD  106
    """
    
    timestamps = []
    data_dict = {}
    errors = []  # List to store error information
    
    get_date(filename)

    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            # Line Validation
            # ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
            is_valid, error_msg, parts = validate_line(line, line_num)

            if not is_valid:
                # Log not valid line
                # ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
                errors.append({
                    'line_num': line_num,
                    'error': error_msg,
                    'raw_line': line[:100]  # First 100 chars
                })
                if verbose:
                    print(f"Warning: Line {line_num} - {error_msg}: {line[:80]}")
                # No continue needed - we just skip to next iteration
            else:
                # Process valid line
                # ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
                # First part is the timestamp
                timestamp_str = parts[0]
                timestamps.append(timestamp_str)
                
                # Parse the remaining parts as variable-value pairs
                i = 1
                while i < len(parts):
                    var_name = parts[i]
                    if i + 1 < len(parts):
                        try:
                            value = float(parts[i + 1])
                            
                            # Initialize list for this variable if not exists
                            if var_name not in data_dict:
                                data_dict[var_name] = []
                            
                            data_dict[var_name].append(value)
                        except ValueError:
                            # If we can't convert to float, skip this pair
                            pass
                    
                    i += 2
    
    # Convert timestamps to pandas datetime
    # ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    # Assuming timestamps are time-only (HH:MM:SS.fff)
    # We'll use today's date as a reference
    if reference_date:
        ref_date = pd.to_datetime(reference_date, format='%d/%m/%Y').date()
    else:
        ref_date = pd.Timestamp.now().date()
    time_objects = pd.to_datetime(timestamps, format='%H:%M:%S.%f').time
    datetime_index = pd.to_datetime([
        pd.Timestamp.combine(ref_date, t) for t in time_objects])
        

    # Create xarray Dataset
    # ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    data_vars = {}
    for var_name, values in data_dict.items():
        data_vars[var_name] = (['time'], values)
    
    ds = xr.Dataset(
        data_vars=data_vars,
        coords={'time': datetime_index},
        attrs=trisonica_infos(),
    )
    
    # Error logging
    # ‾‾‾‾‾‾‾‾‾‾‾‾‾
    # Check for gaps if requested
    gap_info = None
    if check_gaps and len(ds['time']) > 0:
        gap_info = check_missing_data(ds, sampling_rate_hz=sampling_rate_hz)

    # Write error log if there were errors
    if log_errors and (errors or (gap_info and gap_info['has_gaps'])):
        log_filename = os.path.splitext(filename)[0] + '_errors.log'
        print(log_filename)
        with open(log_filename, 'w') as log_file:
            log_file.write(f"Parsing log for file: {filename}\n")
            log_file.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_file.write("="*80 + "\n\n")
            

            if errors:
                log_file.write(f"PARSING ERRORS\n")
                log_file.write(f"Total errors: {len(errors)}\n")
                log_file.write("-"*80 + "\n\n")
                
                for error in errors:
                    log_file.write(f"Line {error['line_num']}: {error['error']}\n")
                    log_file.write(f"  Content: {error['raw_line']}\n\n")
            
            # Write gap section
            if gap_info:
                log_file.write("\n" + "="*80 + "\n")
                log_file.write("MISSING DATA ANALYSIS\n")
                log_file.write("-"*80 + "\n")
                log_file.write(f"Expected sampling interval: {gap_info['expected_interval_ms']:.1f} ms\n")
                log_file.write(f"Total samples found: {gap_info['total_samples']}\n")
                log_file.write(f"Expected samples: {gap_info['expected_samples']}\n")
                log_file.write(f"Missing samples: {gap_info['total_missing_samples']}\n")
                log_file.write(f"Number of gaps: {gap_info['num_gaps']}\n\n")
                
                if gap_info['has_gaps']:
                    log_file.write("Gap details:\n")
                    log_file.write("-"*80 + "\n")
                    for i, gap in enumerate(gap_info['gaps'], 1):
                        log_file.write(f"\nGap {i}:\n")
                        log_file.write(f"  Time before: {gap['time_before']}\n")
                        log_file.write(f"  Time after:  {gap['time_after']}\n")
                        log_file.write(f"  Duration: {gap['gap_duration_ms']:.1f} ms\n")
                        log_file.write(f"  Missing samples: {gap['missing_samples']}\n")
                else:
                    log_file.write("✓ No gaps detected!\n")
        
        if verbose:
            print(f"\nLog file saved to: {log_filename}")
            if errors:
                print(f"Total parsing errors: {len(errors)}")
            if gap_info and gap_info['has_gaps']:
                print(f"Total gaps found: {gap_info['num_gaps']}")
                print(f"Total missing samples: {gap_info['total_missing_samples']}")


    return ds










def get_date(filename):
    """
    Get date from filename, output is string of format DD/MM/YYYY

    Raise an exception if the format isnt recognized.
    """    
    fichier = filename.split('/')[-1].split('_')[1]
    day, month, year = fichier[6:], fichier[4:6], fichier[0:4]
    reference_date = day+'/'+month+'/'+year 
    default_date = '01/01/2000'
    stop = False
    # sanity check
    if not fichier.isdigit():
        print(f'Error in filename: not a date ({fichier})')
        print(' expected format YYYYMMDD')
        reference_date = default_date
    else:
        if year != '2026':
            stop=True
            print(f'year is not 2026 ({year}), check the file ...')

        if month != '02':
            stop=True
            print(f'month is not Februrary ({month}), check the file ...')
        if int(day) < 1 or int(day) > 31:
            stop=True
            print(f'day is not inside [1,31] ({day}), check the file ...')

    if stop:
        raise Exception(f'Filename format is not recognized,\n filename is:\n {filename}')
    else:
        return reference_date


def trisonica_infos():
    dict = {'description':'Anemometer trisonica mini',
            'frequency':'5hz',
            'location':'',
            'height':'',
            'notes':''}

def flag_trisonica(ds):
    """
    Flag the uncorrect data from ds

    Uncorrect is:
    - unphysical value
    - 
    """
    return ds

def check_missing_data(ds, sampling_rate_hz=5, tolerance_ms=50):
    """
    Check for gaps in the time series data.
    
    Parameters:
    -----------
    ds : xr.Dataset
        Dataset with time coordinate
    sampling_rate_hz : float
        Expected sampling rate in Hz (default: 5)
    tolerance_ms : float
        Tolerance in milliseconds for considering a gap (default: 50)
        Gaps larger than (1/sampling_rate + tolerance) are flagged
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'has_gaps': bool, whether gaps were found
        - 'expected_interval_ms': expected interval between samples
        - 'num_gaps': number of gaps found
        - 'gaps': list of dicts with gap details
        - 'total_missing_samples': estimated total number of missing samples
    """
    
    # Calculate expected interval
    expected_interval_ms = 1000.0 / sampling_rate_hz
    max_allowed_interval_ms = expected_interval_ms + tolerance_ms
    
    # Get time differences
    time_diffs = ds['time'].diff(dim='time')
    time_diffs_ms = time_diffs.values / pd.Timedelta(milliseconds=1)
    
    # Find gaps
    gap_indices = time_diffs_ms > max_allowed_interval_ms
    
    gaps = []
    total_missing = 0
    
    for i, is_gap in enumerate(gap_indices):
        if is_gap:
            time_before = ds['time'].values[i]
            time_after = ds['time'].values[i + 1]
            gap_duration = time_diffs_ms[i]
            
            # Estimate number of missing samples
            missing_samples = int(round(gap_duration / expected_interval_ms)) - 1
            total_missing += missing_samples
            
            gaps.append({
                'time_before': pd.Timestamp(time_before),
                'time_after': pd.Timestamp(time_after),
                'gap_duration_ms': gap_duration,
                'missing_samples': missing_samples
            })
    
    result = {
        'has_gaps': len(gaps) > 0,
        'expected_interval_ms': expected_interval_ms,
        'num_gaps': len(gaps),
        'gaps': gaps,
        'total_missing_samples': total_missing,
        'total_samples': len(ds['time']),
        'expected_samples': len(ds['time']) + total_missing
    }
    
    return result

def print_gap_report(gap_info):
    """
    Print a formatted report of the gap analysis.
    
    Parameters:
    -----------
    gap_info : dict
        Output from check_missing_data function
    """
    print(f"=== Missing Data Report ===")
    print(f"Expected sampling interval: {gap_info['expected_interval_ms']:.1f} ms")
    print(f"Total samples found: {gap_info['total_samples']}")
    print(f"Expected samples: {gap_info['expected_samples']}")
    print(f"Missing samples: {gap_info['total_missing_samples']}")
    print(f"Number of gaps: {gap_info['num_gaps']}")
    
    if gap_info['has_gaps']:
        print(f"\nGap details:")
        for i, gap in enumerate(gap_info['gaps'], 1):
            print(f"\n  Gap {i}:")
            print(f"    Before: {gap['time_before']}")
            print(f"    After:  {gap['time_after']}")
            print(f"    Duration: {gap['gap_duration_ms']:.1f} ms")
            print(f"    Missing samples: {gap['missing_samples']}")
    else:
        print("\n✓ No gaps detected!")


def validate_line(line, line_num):
    """
    Validate a single line from the log file.
    
    Parameters:
    -----------
    line : str
        The line to validate
    line_num : int
        Line number in the file
    
    Returns:
    --------
    tuple : (is_valid, error_message, parts)
        is_valid: bool, whether the line is valid
        error_message: str or None, description of the error if invalid
        parts: list or None, split line parts if valid
    """

    N_total = 23 # the file should have 23 parts for each line

    if not line:
        return False, "Empty line", None
    
    parts = line.split()
    
    # Check minimum number of elements
    if len(parts) < N_total:
        return False, f"Too few elements ({len(parts)} parts)", None
    
    if len(parts) > N_total:
        return False, f"Too much elements ({len(parts)} parts)", None

    # Check if number of parts is odd (timestamp + pairs of var-value)
    if len(parts) % 2 != 1:
        return False, f"Odd number of elements ({len(parts)} parts)", None
    
    

    # Validate timestamp format (HH:MM:SS.fff)
    timestamp_str = parts[0]
    try:
        pd.to_datetime(timestamp_str, format='%H:%M:%S.%f')
    except ValueError:
        return False, f"Invalid timestamp format '{timestamp_str}'", None
    
    return True, None, parts
