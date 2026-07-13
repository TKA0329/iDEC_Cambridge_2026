#!/usr/bin/env python3
"""
Test script to identify performance bottlenecks.
"""

import time
import sys

print("Testing import speed...")
start = time.time()
from protein_pipeline import analyze_sequence, CIDER_CACHE
elapsed = time.time() - start
print(f"  protein_pipeline import: {elapsed:.2f}s\n")

print("Testing first sequence analysis (warmup)...")
test_seq = "MKSFQNWREQLAKRFQGGKDDDDEEEEKKKKRRRR"
start = time.time()
result1 = analyze_sequence(test_seq)
elapsed = time.time() - start
print(f"  First sequence: {elapsed:.3f}s")
print(f"  Result keys: {list(result1.keys())}\n")

print("Testing 10 sequences (should be faster after warmup)...")
start = time.time()
for i in range(10):
    result = analyze_sequence(test_seq)
elapsed = time.time() - start
avg = elapsed / 10
print(f"  10 sequences: {elapsed:.2f}s total")
print(f"  Average per sequence: {avg*1000:.1f}ms")
print(f"  Cache size: {len(CIDER_CACHE)} entries\n")

print("Testing 10 MORE sequences (should benefit from caching)...")
start = time.time()
for i in range(10):
    result = analyze_sequence(test_seq)
elapsed2 = time.time() - start
avg2 = elapsed2 / 10
print(f"  10 cached sequences: {elapsed2:.2f}s total")
print(f"  Average per cached sequence: {avg2*1000:.1f}ms")
print(f"  Speedup factor: {avg/avg2:.1f}x\n")

print("Estimated time for 3.7M sequences:")
estimated_seconds = 3_747_989 * avg
estimated_hours = estimated_seconds / 3600
print(f"  Single-threaded: {estimated_hours:.1f} hours")
print(f"  With 16 threads: {estimated_hours/16:.1f} hours")
