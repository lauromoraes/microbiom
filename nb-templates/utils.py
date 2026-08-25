#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import os
import tempfile

from qiime2.plugins.feature_table.methods import filter_samples
from qiime2.plugins.feature_table.methods import filter_seqs

def filter_samples(metadata, tabs, reps)
    # Filter FeatureTable[Frequency | RelativeFrequency | PresenceAbsence | Composition] based on Metadata sample ID values
    tabs = filter_samples(
        table=tabs,
        metadata=metadata,
    ).filtered_table
    # Filter SampleData[SequencesWithQuality | PairedEndSequencesWithQuality | JoinedSequencesWithQuality] based on Metadata sample ID values; returns FeatureData[Sequence | AlignedSequence]
    reps = filter_seqs(
        data=reps,
        table=tabs,
    ).filtered_data
    return tabs, reps



def get_direction_info(manifest_file_path):
    """Process and gets information about fastq directions using manifest file
    
    Example:
        d_type, v_type, direction = get_direction_info(manifest_file)
        print('\n'.join([d_type, v_type, direction]))

    Parameters:
    manifest_file_path (str): manifest file with all sample input paths

    Returns:
    str: QIIME2 input data type
    str: QIIME2 Manifest file type
    str: Direction type

    """
    d_type, v_type, direction = None, None, None
    manifest_df = pd.read_csv(manifest_file)
    n_directions = len(manifest_df['direction'].unique())
    if n_directions == 1:
        d_type = 'SampleData[SequencesWithQuality]'
        v_type = 'SingleEndFastqManifestPhred33'
        direction = 'single'
    elif n_directions == 2:
        d_type = 'SampleData[PairedEndSequencesWithQuality]'
        v_type = 'PairedEndFastqManifestPhred33'
        direction = 'paired'
    else:
        print(f'ERROR: invalid number of directions {n_directions}')
    return d_type, v_type, direction

def setup_temp_environment(
    base_dir: str | Path,
    qiime_subdir: str = "qiime2",
    joblib_subdir: str = "joblib",
    cache_subdir: str = "qiime2-cache",
    verbose: bool = True,
) -> dict[str, Path]:

    base_dir = Path(base_dir).expanduser().resolve()

    tmp_dir = base_dir / qiime_subdir
    joblib_dir = base_dir / joblib_subdir
    cache_dir = base_dir / cache_subdir

    for directory in (tmp_dir, joblib_dir, cache_dir):
        directory.mkdir(parents=True, exist_ok=True)

    os.environ["TMPDIR"] = str(tmp_dir)
    os.environ["TMP"] = str(tmp_dir)
    os.environ["TEMP"] = str(tmp_dir)
    os.environ["JOBLIB_TEMP_FOLDER"] = str(joblib_dir)

    tempfile.tempdir = str(tmp_dir)

    paths = {
        "tmp": tmp_dir,
        "joblib": joblib_dir,
        "qiime_cache": cache_dir,
    }

    if verbose:
        print("Temporary environment configured:")
        print(f"  TMPDIR             = {tmp_dir}")
        print(f"  JOBLIB_TEMP_FOLDER = {joblib_dir}")
        print(f"  QIIME cache        = {cache_dir}")

    return paths

