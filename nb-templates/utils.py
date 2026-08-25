#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import os
import tempfile
import pandas as pd
import numpy as np

from qiime2.plugins.feature_table.methods import filter_samples
from qiime2.plugins.feature_table.methods import filter_seqs

def filter_samples_seqs(metadata, tabs, reps):
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
    """Process and get information about FASTQ directions using the manifest file
    
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

def get_deepest_taxonomic_level(taxonomy) -> int:
    """
    Return the deepest taxonomic level found in QIIME 2
    FeatureData[Taxonomy] artifact.

    Taxonomic levels:
        1 = Domain
        2 = Phylum
        3 = Class
        4 = Order
        5 = Family
        6 = Genus
        7 = Species

    Parameters
    ----------
    taxonomy: QIIME 2.Artifact
        QIIME 2 FeatureData[Taxonomy] artifact.

    Returns
    -------
    int
        Deepest taxonomic level found.

    Raises
    ------
    TypeError
        If the artifact is not FeatureData[Taxonomy].
    ValueError
        If no valid taxonomic classification is found.
    """

    if str(taxonomy.type) != "FeatureData[Taxonomy]":
        raise TypeError(
            f"Expected FeatureData[Taxonomy], got {taxonomy.type}"
        )

    taxonomy_df = taxonomy.view(pd.DataFrame)

    if "Taxon" not in taxonomy_df.columns:
        raise ValueError("Taxonomy table does not contain a 'Taxon' column.")

    prefixes = {
        "d__": 1,
        "k__": 1,  # compatibility with older classifiers
        "p__": 2,
        "c__": 3,
        "o__": 4,
        "f__": 5,
        "g__": 6,
        "s__": 7,
    }

    deepest_level = 0

    for taxon in taxonomy_df["Taxon"].dropna():
        for rank in str(taxon).split(";"):
            rank = rank.strip()

            for prefix, level in prefixes.items():
                if rank.startswith(prefix):
                    # Ignore empty classifications such as g__ or s__
                    value = rank[len(prefix):].strip()

                    if value:
                        deepest_level = max(deepest_level, level)

    if deepest_level == 0:
        raise ValueError("No valid taxonomic levels were found.")

    return deepest_level

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

