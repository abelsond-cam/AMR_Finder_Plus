#!/usr/bin/env python
"""
Test script to verify multiprocessing is working correctly
"""

import os
import time
from multiprocessing import Pool, cpu_count, current_process
from datetime import datetime

def worker_function(x):
    """Simple worker that sleeps and returns info"""
    worker_name = current_process().name
    pid = os.getpid()
    
    # Simulate work
    time.sleep(2)
    
    return {
        'input': x,
        'worker': worker_name,
        'pid': pid,
        'timestamp': datetime.now().strftime('%H:%M:%S')
    }

def main():
    n_cpus = cpu_count()
    n_tasks = 10
    
    print(f"Testing multiprocessing with {n_cpus} CPUs")
    print(f"Main process PID: {os.getpid()}")
    print(f"Creating {n_tasks} tasks...")
    print(f"Expected: {min(n_cpus, n_tasks)} processes running in parallel\n")
    
    start = datetime.now()
    
    with Pool(processes=n_cpus) as pool:
        results = pool.map(worker_function, range(n_tasks))
    
    elapsed = (datetime.now() - start).total_seconds()
    
    print(f"\nResults:")
    for r in results:
        print(f"  Task {r['input']}: PID {r['pid']}, Worker {r['worker']}, Time {r['timestamp']}")
    
    # Count unique PIDs
    unique_pids = len(set(r['pid'] for r in results))
    
    print(f"\nSummary:")
    print(f"  Total tasks: {n_tasks}")
    print(f"  Unique worker PIDs: {unique_pids}")
    print(f"  Total time: {elapsed:.1f}s")
    print(f"  Expected time (parallel): ~{2 * (n_tasks / n_cpus):.1f}s")
    print(f"  Expected time (sequential): ~{2 * n_tasks:.1f}s")
    
    if unique_pids > 1:
        print(f"\n✓ Multiprocessing is WORKING! ({unique_pids} workers used)")
    else:
        print(f"\n✗ Multiprocessing NOT working (only 1 process)")
    
    if elapsed < (2 * n_tasks * 0.8):
        print(f"✓ Parallel execution confirmed (fast)")
    else:
        print(f"✗ Sequential execution (slow)")

if __name__ == "__main__":
    main()

