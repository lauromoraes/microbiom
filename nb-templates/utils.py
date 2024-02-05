#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
